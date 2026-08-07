# The list of Turing.jl models benchmarked in TuringADTests.jmd, grouped by
# category. Each name must match both `models/<name>.jl` and the name of the
# `@model function` defined inside that file.
#
# The model files in `models/` and the list below are taken from
# https://github.com/TuringLang/ADTests, whose results are published at
# https://turinglang.org/ADTests/.

const MODEL_CATEGORIES = [
    "Base Julia features" => [
        "control_flow",
        "threaded_assume",
        "threaded_observe",
    ],
    "Core Turing syntax" => [
        "assume_submodel",
        "broadcast_macro",
        "dot_assume",
        "dot_observe",
        "dynamic_constraint",
        "multiple_constraints_same_var",
        "observe_index",
        "observe_literal",
        "observe_multivariate",
        "observe_submodel",
    ],
    "Distributions" => [
        "assume_beta",
        "assume_dirichlet",
        "assume_lkjcholu",
        "assume_mvnormal",
        "assume_normal",
        "assume_wishart",
        "observe_bernoulli",
        "observe_categorical",
        "observe_von_mises",
    ],
    "DynamicPPL arXiv paper" => [
        "dppl_gauss_unknown",
        "dppl_hier_poisson",
        "dppl_high_dim_gauss",
        "dppl_hmm_semisup",
        "dppl_lda",
        "dppl_logistic_regression",
        "dppl_naive_bayes",
        "dppl_sto_volatility",
    ],
    "DynamicPPL demo models" => [
        "demo_assume_dot_observe",
        "demo_assume_dot_observe_literal",
        "demo_assume_index_observe",
        "demo_assume_matrix_observe_matrix_index",
        "demo_assume_multivariate_observe",
        "demo_assume_multivariate_observe_literal",
        "demo_assume_observe_literal",
        "demo_assume_submodel_observe_index_literal",
        "demo_dot_assume_observe",
        "demo_dot_assume_observe_index",
        "demo_dot_assume_observe_index_literal",
        "demo_dot_assume_observe_matrix_index",
        "demo_dot_assume_observe_submodel",
    ],
    "Effect of model size" => [
        "n010",
        "n050",
        "n100",
        "n500",
    ],
    "PosteriorDB" => [
        "pdb_arma11",
        "pdb_earnings",
        "pdb_earnings_male",
        "pdb_eightsch_centered",
        "pdb_eightsch_noncentered",
        "pdb_garch11",
        "pdb_kidiq",
        "pdb_radon",
        "pdb_rats",
        "pdb_sblrc",
        "pdb_sblri",
    ],
    "External libraries" => [
        "abstractgps",
        "delaydiffeq",
        "lux_nn",
        "ordinarydiffeq",
    ],
]

const MODEL_NAMES = reduce(vcat, last.(MODEL_CATEGORIES))

const CATEGORY_OF = Dict(
    name => category for (category, names) in MODEL_CATEGORIES for name in names
)

# AD backends, in the order they are reported. FiniteDifferences doubles as the
# reference against which every other backend's gradient is checked.
const REFERENCE_BACKEND = "FiniteDifferences"

const BACKEND_NAMES = [
    "FiniteDifferences",
    "ForwardDiff",
    "ReverseDiff",
    "ReverseDiffCompiled",
    "MooncakeRvs",
    "MooncakeFwd",
    "EnzymeFwd",
    "EnzymeRvs",
]
