#/////////////////////////////////////////////////////////////////////////////////
# INTERFACE TO RUN MUPLTIPLE EXAMPLES WITH DIFFERENT STRATEGIES / SETTINGS
#/////////////////////////////////////////////////////////////////////////////////
using Plots
# Import all the examples
#include("./nernst_planck_3D.jl")
#include("./level_set.jl")
#include("./allen_cahn.jl")
include("./hamilton_jacobi.jl")


# Settings:
maxIters   = 5  #number of iterations


strategies = [NeuralPDE.QuadratureTraining(quadrature_alg = CubaCuhre(), reltol = 1e-4, abstol = 1e-3, maxiters = 10, batch = 10),
              NeuralPDE.QuadratureTraining(quadrature_alg = HCubatureJL(), reltol = 1e-4, abstol = 1e-3, maxiters = 10, batch = 0),
              NeuralPDE.QuadratureTraining(quadrature_alg = CubatureJLh(), reltol = 1e-4, abstol = 1e-3, maxiters = 10, batch = 10),
              NeuralPDE.QuadratureTraining(quadrature_alg = CubatureJLp(), reltol = 1e-4, abstol = 1e-3, maxiters = 10, batch = 10),
              NeuralPDE.GridTraining(0.1),
              NeuralPDE.StochasticTraining(100),
              NeuralPDE.QuasiRandomTraining(100; sampling_alg = UniformSample(), minibatch = 100)]

strategies_short_name = ["CubaCuhre",
                        "HCubatureJL",
                        "CubatureJLh",
                        "CubatureJLp",
                        "GridTraining",
                        "StochasticTraining",
                        "QuasiRandomTraining"]

minimizers = [GalacticOptim.ADAM(0.01)]
             #[GalacticOptim.BFGS()]
              #GalacticOptim.LBFGS()]

minimizers_short_name = ["ADAM"]
                        # "BFGS",
                        # "LBFGS"]

# Run models
numeric_res = Dict()
prediction_res =  Dict()
benchmark_res = Dict()
error_res =  Dict()
domains = Dict()
params_res = Dict()  #to use same params for the next run


print("\n Starting first run for time benchmark")
## Time Benchmark
for strat=1:length(strategies) # strategy
      for min =1:length(minimizers) # minimizer
            res = hamilton_jacobi_time(strategies[strat], minimizers[min], maxIters)
            push!(benchmark_res, string(strat,min)  => res[1])
            push!(params_res, string(strat,min)  => res[2])
            print(string("Training time (", strategies_short_name[strat], " ", minimizers_short_name[min], ") = ",(res[1])))
      end
end


benchmark_res_name = Dict()
for strat=1:length(strategies) # strategy
      for min =1:length(minimizers)
            push!(benchmark_res_name, string(strategies_short_name[strat]) => benchmark_res[string(strat,min)])#, " + " , minimizers_short_name[min]) => benchmark_res[string(strat,min)])
      end
end

print("\n Plotting time benchmark")
bar = Plots.bar(collect(keys(benchmark_res_name)), collect(values(benchmark_res_name)),
      title = string("HJ training time (ADAM(0.01))"), titlefontsize = 12, yrotation = 30, orientation= :horizontal,
      legend = false, xlabel = "Seconds", xguidefontsize=9, size = (600,400))

#Plots.savefig("HJ_training_time.pdf")



print("\n Starting second run for convergence")
## Convergence

for strat=1:length(strategies) # strategy
      for min =1:length(minimizers) # minimizer
            println(string(strategies_short_name[strat], "  ", minimizers_short_name[min]))
            res = hamilton_jacobi(strategies[strat], minimizers[min], maxIters, params_res[string(strat,min)])
            push!(error_res, string(strat,min)     => res[1])
            push!(prediction_res, string(strat,min) => res[2])
            push!(numeric_res, string(strat,min)    => res[3])
            push!(domains, string(strat,min)        => res[4])
      end
end

print("\n Plotting error vs iters")
#Plotting the first strategy with the first minimizer out from the loop to initialize the canvas
current_label = string(strategies_short_name[1], " + " , minimizers_short_name[1])
error = Plots.plot(1:(maxIters + 1), error_res["11"], yaxis=:log10, title = string("HJ"), ylabel = "log(error)", label = current_label)#legend = true)#, size=(1200,700))
#=for strat=2:length(strategies) # strategy
      for min =1:length(minimizers) # minimizer
            # Learning curve plots with different strategies, minimizer
            current_label = string(strategies_short_name[strat], " + " , minimizers_short_name[min])
            plot!(error, 1:(maxIters + 1), error_res[string(strat,min)], yaxis=:log10, title = string("Allen_Cahn"), ylabel = "log(error)", label = current_label)#, size=(2000,900))
      end
end
=#

plot!(error, 1:(maxIters + 1), error_res[string(2,1)], yaxis=:log10, title = string("HJ"), ylabel = "log(error)", label = string(strategies_short_name[2], " + " , minimizers_short_name[1]))
plot!(error, 1:(maxIters + 1), error_res[string(3,1)], yaxis=:log10, title = string("HJ"), ylabel = "log(error)", label = string(strategies_short_name[3], " + " , minimizers_short_name[1]))
plot!(error, 1:(maxIters + 1), error_res[string(4,1)], yaxis=:log10, title = string("HJ"), ylabel = "log(error)", label = string(strategies_short_name[4], " + " , minimizers_short_name[1]))
plot!(error, 1:(maxIters + 1), error_res[string(5,1)], yaxis=:log10, title = string("HJ"), ylabel = "log(error)", label = string(strategies_short_name[5], " + " , minimizers_short_name[1]))
plot!(error, 1:(maxIters + 1), error_res[string(6,1)], yaxis=:log10, title = string("HJ"), ylabel = "log(error)", label = string(strategies_short_name[6], " + " , minimizers_short_name[1]))
plot!(error, 1:(maxIters + 1), error_res[string(7,1)], yaxis=:log10, title = string("HJ convergence (ADAM (0.01))"), ylabel = "log(error)", label = string(strategies_short_name[7], " + " , minimizers_short_name[1]))

#Plots.plot!(1:(maxIters + 1), error_res["11"], yaxis=:log10, title = string("Allen_Cahn"), ylabel = "log(loss)")# size=(600,350))

#Plots.savefig("HJ_error.pdf")

Plots.plot(error, bar, layout = Plots.grid(1, 2, widths=[0.6 ,0.4]), size = (1500,500))

Plots.savefig("HJ_comparison.pdf")

#=
## Comparison Predicted solution vs Numerical solution
integrals = Dict()

to_print = [] #["11", "32"]

for strat=1:5 #length(strategies) # strategy
      for min =1:length(minimizers)
            u_predict = prediction_res[string(strat,min)][1]
            u_real = numeric_res[string(strat,min)][1]

            # dimensions
            ts = collect(domains[string(strat,min)][1])
            xs = collect(domains[string(strat,min)][2])
            ys = collect(domains[string(strat,min)][3])
            zs = collect(domains[string(strat,min)][4])
            XS = [xs,ys,zs]

            diff_u = abs.(u_predict .- u_real)

            if string(strat,min) ∈ to_print
                  """
                  p1 = Plots.plot(XS, ts, u_real, linetype=:contourf,title = "analytic");
                  p2 = Plots.plot(u_predict[1,:,:,1], linetype=:contourf, title = "predict");
                  p3 = Plots.plot(XS, ts, diff_u,linetype=:contourf,title = "error");
                  savefig(p1,"p1.pdf")
                  savefig(p2,"p2.pdf")
                  savefig(p3,"p3.pdf")
                  """
            end

            integral = sum(diff_u)[1]
            push!(integrals, string(strategies_short_name[strat], " + " , minimizers_short_name[min]) => integral)

      end
end

Plots.bar(collect(keys(integrals)), collect(values(integrals)), title = string("Nernst-Planck Error Integrals"), xrotation = 90)
savefig("Nernst-Planck_strategies_MAEs.pdf")

=#
