module CUTEstWorkerTests

using CUTEst
using Pkg
using Test

include("cutest_benchmark_utils.jl")

# Keep the solver budget unchanged while testing deliberate hangs.
hard_timeout_seconds(::Float32) = 25.0

function child_running(pid)
    status = read(ignorestatus(`ps -o stat= -p $pid`), String)
    return !isempty(strip(status)) && !startswith(strip(status), "Z")
end

try
    @testset "Worker startup from a document module" begin
        worker = ensure_benchmark_worker()
        @test worker > 1
        @test remotecall_fetch(getpid, worker) == BENCHMARK_WORKER_PGID[]
        @test remotecall_fetch(Core.eval, worker, Main, :(ccall(:getpgrp, Cint, ()))) ==
            BENCHMARK_WORKER_PGID[]
        @test remotecall_fetch(Core.eval, worker, Main, :(Threads.nthreads())) == Threads.nthreads()
        @test run_single_solve_guarded("ROSENBR", "LBFGS").status == "OK"
    end

    @testset "Warm-up survives a mid-list hang" begin
        worker = ensure_benchmark_worker()
        warmup_solvers(["LBFGS"]; hard_timeout = true)
        @test solver_warmed_here("LBFGS")
        remotecall_fetch(
            Core.eval, worker, Main, quote
                function run_single_solve(problem::String, solver::String; maxtime)
                    @assert solver == "ConjugateGradient"
                    sleep(120)
                end
                nothing
            end
        )
        rows = warmup_solvers(
            ["LBFGS", "ConjugateGradient", "NelderMead"];
            hard_timeout = true, maxtime = Float32(90),
        )
        @test rows.retcode[1] == "HARD_TIMEOUT"
        @test rows.status[2] == "OK"
        @test BENCHMARK_WORKER[] != worker
        @test WARMED_SOLVERS == Set(["NelderMead"])
        @test WARMUP_FAILED_SOLVERS == Set(["ConjugateGradient"])
        ensure_solver_warmed!("LBFGS")
        @test solver_warmed_here("LBFGS")
        @test nrow(warmup_solvers(["ConjugateGradient"]; hard_timeout = true)) == 0
    end

    @testset "Non-yielding solve timeout kills descendants" begin
        worker = ensure_benchmark_worker()
        mktemp() do pidfile, _
            remotecall_fetch(
                Core.eval, worker, Main, quote
                    function run_single_solve(problem::String, solver::String; maxtime)
                        child = run(`sleep 120`; wait = false)
                        write($pidfile, string(getpid(child)))
                        ccall(:sleep, Cuint, (Cuint,), 120)
                    end
                    nothing
                end
            )
            rows = run_benchmarks(
                "worker test", fill("ROSENBR", 5), ["LBFGS"];
                warmup = false, maxtime = Float32(90),
            )
            @test rows.retcode[1] == "HARD_TIMEOUT"
            @test rows.status[2:end] == fill("OK", 4)
            @test BENCHMARK_WORKER[] != worker
            @test isempty(WARMED_SOLVERS)
            pid = parse(Int, read(pidfile, String))
            @test timedwait(() -> !child_running(pid), 5.0) == :ok
        end
    end

    @testset "Failed initialization removes the worker" begin
        kill_benchmark_worker(BENCHMARK_WORKER[])
        previous_project = Base.active_project()
        mktempdir() do empty_project
            try
                Pkg.activate(empty_project)
                withenv("JULIA_LOAD_PATH" => "@:@stdlib") do
                    @test_throws TaskFailedException ensure_benchmark_worker()
                end
                @test nprocs() == 1
                @test BENCHMARK_WORKER[] == 0
                @test BENCHMARK_WORKER_PGID[] == 0
            finally
                Pkg.activate(dirname(previous_project))
            end
        end
    end
finally
    for worker in filter(!=(1), workers())
        if worker == BENCHMARK_WORKER[]
            kill_benchmark_worker(worker)
        else
            rmprocs(worker; waitfor = 60)
        end
    end
end

end
