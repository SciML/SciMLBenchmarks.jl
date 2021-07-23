# benchmarks for various grids
using DiffEqJump, BenchmarkTools, LightGraphs
using ProgressLogging

function ABC_setup(spatial_system, linear_size, starting_site, u0, diffusivity, end_time)
    domain_size = 1.0 #μ-meter
    mesh_size = domain_size/linear_size
    # ABC model A + B <--> C
    reactstoch = [[1 => 1, 2 => 1],[3 => 1]]
    netstoch = [[1 => -1, 2 => -1, 3 => 1],[1 => 1, 2 => 1, 3 => -1]]
    rates = [0.1/mesh_size, 1.]
    majumps = MassActionJump(rates, reactstoch, netstoch)

    # spatial system setup
    hopping_rate = diffusivity / mesh_size^2
    num_nodes = DiffEqJump.num_sites(spatial_system)

    # Starting state setup
    starting_state = zeros(Int, length(u0), num_nodes)
    starting_state[:,starting_site] = copy(u0)
    prob = DiscreteProblem(starting_state,(0.0,end_time), rates)

    hopping_constants = [hopping_rate for i in starting_state]

    alg = NSM()
    return JumpProblem(prob, alg, majumps, hopping_constants=hopping_constants, spatial_system = spatial_system, save_positions=(false,false))
end

# benchmark solving pure diffusion
diffusivity = 0.1
linear_size = 64
dim = 3
dims = Tuple([linear_size for i in 1:dim])
num_species = 1

majumps = MassActionJump([], [], [])
tf = 0.5
tspan = (0.0, tf)
u0 = [10^4]
domain_size = 1.0 #μ-meter
num_nodes = prod(dims)

starting_state = zeros(Int, length(u0), num_nodes)
center_node = trunc(Int,(num_nodes+1)/2)
starting_state[:,center_node] = copy(u0)
prob = DiscreteProblem(starting_state,(0.0,tf), [])

hopping_rate = diffusivity * (linear_size/domain_size)^2
hopping_constants = [hopping_rate for i in starting_state]

algs = [NSM(), RDirect(), RSSA(), RSSACR()]

grids = [DiffEqJump.CartesianGrid1(dims), DiffEqJump.CartesianGrid2(dims), DiffEqJump.CartesianGrid3(dims), LightGraphs.grid(dims)]
jump_problems = JumpProblem[JumpProblem(prob, NSM(), majumps, hopping_constants=hopping_constants, spatial_system = grid, save_positions=(false,false)) for grid in grids]
# setup flattenned jump prob
graph = LightGraphs.grid(dims)
hopping_constants = Vector{Matrix{Float64}}(undef, num_nodes)
for site in 1:num_nodes
    hopping_constants[site] = hopping_rate*ones(num_species, DiffEqJump.num_neighbors(graph, site))
end
for alg in algs
    push!(jump_problems, DiffEqJump.flatten([], [], [], graph, starting_state, tspan, alg, hopping_constants, save_positions=(false,false)))
end
solve_diffusion_benchmarks = Vector{BenchmarkTools.Trial}(undef, length(jump_problems))

@progress "diffusion benchmarking on a $dims grid" for (i,spatial_jump_prob) in enumerate(jump_problems)
    b = @benchmarkable solve($spatial_jump_prob, SSAStepper())
    solve_diffusion_benchmarks[i] = run(b, samples = 15, seconds = 600)
end

# Benchmarking sampling
dims = (256,256,256)
sites = 1:prod(dims)
grids = [DiffEqJump.CartesianGrid1(dims), DiffEqJump.CartesianGrid2(dims), DiffEqJump.CartesianGrid3(dims), LightGraphs.grid(dims)]
sample_benchmarks = Vector{BenchmarkTools.Trial}(undef, length(grids))
for (i,grid) in enumerate(grids)
    sample_benchmarks[i] = @benchmark DiffEqJump.rand_nbr($grid, site) setup = (site = rand(sites))
end

# Benchmarking solving A+B <-> C
lin_size = 4
dims = Tuple((lin_size,lin_size,lin_size))
starting_site = trunc(Int,(prod(dims) + 1)/2)
u0 = [500,500,0]
end_time = 10.0
diffusivity = 1.0

grids = [DiffEqJump.CartesianGrid1(dims), DiffEqJump.CartesianGrid2(dims), DiffEqJump.CartesianGrid3(dims), LightGraphs.grid(dims)]
solve_ABC_benchmarks = Vector{BenchmarkTools.Trial}(undef, length(grids))

for (i,grid) in enumerate(grids)
    spatial_jump_prob = ABC_setup(grid, lin_size, starting_site, u0, diffusivity, end_time)
    solve(spatial_jump_prob, SSAStepper())
    solve_ABC_benchmarks[i] = @benchmark solve($spatial_jump_prob, SSAStepper())
end
