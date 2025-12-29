using ExplicitImports
using SciMLBenchmarks
using Test

@testset "ExplicitImports" begin
    @test check_no_implicit_imports(SciMLBenchmarks) === nothing
    @test check_no_stale_explicit_imports(SciMLBenchmarks) === nothing
end
