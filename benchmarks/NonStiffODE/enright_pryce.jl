using ModelingToolkit
using Symbolics
using Symbolics: unwrap
using DiffEqBase, StaticArrays, LinearAlgebra

@independent_variables t
@variables y(t)[1:10]
y = collect(y)
D = Differential(t)

# Stiff
sa1eqs = [
    D(y[1]) ~ -0.5 * y[1], D(y[2]) ~ -y[2], D(y[3]) ~ -100 * y[3], D(y[4]) ~ -90 * y[4]]
@mtkcompile sa1sys = ODESystem(sa1eqs, t)

sa1prob = ODEProblem{false}(sa1sys, y[1:4] .=> 1.0, (0, 20.0), dt = 1e-2; u0_constructor = x -> SVector(x...))

sa2eqs = [D(y[1]) ~ -1800 * y[1] + 900 * y[2]]
for i in 2:8
    push!(sa2eqs, D(y[i]) ~ y[i - 1] - 2 * y[i] + y[i + 1])
end
push!(sa2eqs, D(y[9]) ~ 1000 * y[8] - 2000 * y[9] + 1000)
@mtkcompile sa2sys = ODESystem(sa2eqs, t)

sa2prob = ODEProblem{false}(sa2sys, y[1:9] .=> 0.0, (0, 120.0), dt = 5e-4; u0_constructor = x -> SVector(x...))

sa3eqs = [
    D(y[1]) ~ -1e4 * y[1] + 100 * y[2] - 10 * y[3] + y[4],
    D(y[2]) ~ -1e3 * y[2] + 10 * y[3] - 10 * y[4],
    D(y[3]) ~ -y[3] + 10 * y[4],
    D(y[4]) ~ -0.1 * y[4]
]
@mtkcompile sa3sys = ODESystem(sa3eqs, t)

sa3prob = ODEProblem{false}(sa3sys, y[1:4] .=> 1.0, (0, 20.0), dt = 1e-5; u0_constructor = x -> SVector(x...))

sa4eqs = [D(y[i]) ~ -i^5 * y[i] for i in 1:10]
@named sa4sys_raw = ODESystem(sa4eqs, t)
sa4sys = structural_simplify(sa4sys_raw)

sa4prob = ODEProblem{false}(sa4sys, y[1:10] .=> 1.0, (0, 1.0), dt = 1e-5; u0_constructor = x -> SVector(x...))

const SA_PROBLEMS = [sa1prob, sa2prob, sa3prob, sa4prob]

sb1eqs = [D(y[1]) ~ -y[1] + y[2],
    D(y[2]) ~ -100y[1] - y[2],
    D(y[3]) ~ -100y[3] + y[4],
    D(y[4]) ~ -10_000y[3] - 100y[4]]
@named sb1sys_raw = ODESystem(sb1eqs, t)
sb1sys = structural_simplify(sb1sys_raw)

sb1prob = ODEProblem{false}(
    sb1sys, [y[[1, 3]] .=> 1.0; y[[2, 4]] .=> 0.0], (0, 20.0), dt = 7e-3; u0_constructor = x -> SVector(x...))

@parameters α
sb2eqs = [D(y[1]) ~ -10y[1] + α * y[2],
    D(y[2]) ~ -α * y[1] - 10 * y[2],
    D(y[3]) ~ -4y[3],
    D(y[4]) ~ -y[4],
    D(y[5]) ~ -0.5y[5],
    D(y[6]) ~ -0.1y[6]]
@named sb2sys_raw = ODESystem(sb2eqs, t)
sb2sys = structural_simplify(sb2sys_raw)

sb2prob = ODEProblem{false}(sb2sys, y .=> 1.0, (0, 20.0), [α => 3], dt = 1e-2)
sb3prob = ODEProblem{false}(sb2sys, y .=> 1.0, (0, 20.0), [α => 8], dt = 1e-2)
sb4prob = ODEProblem{false}(sb2sys, y .=> 1.0, (0, 20.0), [α => 25], dt = 1e-2)
sb5prob = ODEProblem{false}(sb2sys, y .=> 1.0, (0, 20.0), [α => 100], dt = 1e-2)

sc1eqs = [
    D(y[1]) ~ -y[1] + y[2]^2 + y[3]^2 + y[4]^2,
    D(y[2]) ~ -10y[2] + 10 * (y[3]^2 + y[4]^2),
    D(y[3]) ~ -40y[3] + 40 * y[4]^2,
    D(y[4]) ~ -100y[4] + 2]
@named sc1sys_raw = ODESystem(sc1eqs, t)
sc1sys = structural_simplify(sc1sys_raw)

sc1prob = ODEProblem{false}(sc1sys, y .=> 1.0, (0, 20.0), dt = 1e-2; u0_constructor = x -> SVector(x...))

@parameters β
sc2eqs = [D(y[1]) ~ -y[1] + 2,
    D(y[2]) ~ -10y[2] + β * y[1]^2,
    D(y[3]) ~ -40y[3] + 4β * (y[1]^2 + y[2]^2),
    D(y[4]) ~ -100y[4] + 10β * (y[1]^2 + y[2]^2 + y[3]^2)]
@mtkcompile sc2sys = ODESystem(sc2eqs, t)

sc2prob = ODEProblem{false}(sc2sys, y .=> 1.0, (0, 20.0), [β => 0.1], dt = 1e-2; u0_constructor = x -> SVector(x...))
sc3prob = ODEProblem{false}(sc2sys, y .=> 1.0, (0, 20.0), [β => 1.0], dt = 1e-2; u0_constructor = x -> SVector(x...))
sc4prob = ODEProblem{false}(sc2sys, y .=> 1.0, (0, 20.0), [β => 10.0], dt = 1e-2; u0_constructor = x -> SVector(x...))
sc5prob = ODEProblem{false}(sc2sys, y .=> 1.0, (0, 20.0), [β => 20.0], dt = 1e-2; u0_constructor = x -> SVector(x...))

sd1sys = let
    sd1eqs = [D(y[1]) ~ 0.2 * (y[2] - y[1]),
        D(y[2]) ~ 10y[1] - (60 - 0.125y[3]) * y[2] + 0.125y[3],
        D(y[3]) ~ 1]

    structural_simplify(ODESystem(sd1eqs, t; name=Symbol(gensym("sys"))))
end

sd1prob = ODEProblem{false}(sd1sys, y .=> 0.0, (0, 400.0), [β => 0.1], dt = 1.7e-2; u0_constructor = x -> SVector(x...))

sd2sys = let
    sd2eqs = [D(y[1]) ~ -0.04y[1] + 0.01 * (y[2] * y[3]),
        D(y[2]) ~ 400y[1] - 100 * (y[2] * y[3]) - 3000 * y[2]^2,
        D(y[3]) ~ 30y[2]^2]

    structural_simplify(ODESystem(sd2eqs, t; name=Symbol(gensym("sys"))))
end

sd2prob = ODEProblem{false}(
    sd2sys, [y[1] => 1.0; y[2:3] .=> 0.0], (0, 40.0), dt = 1e-5, cse = true; u0_constructor = x -> SVector(x...))

sd3sys = let
    sd3eqs = [D(y[1]) ~ y[3] - 100 * (y[1] * y[2]),
        D(y[2]) ~ y[3] + 2y[4] - 100 * (y[1] * y[2]) - 2e4 * y[2]^2,
        D(y[3]) ~ -y[3] + 100 * (y[1] * y[2]),
        D(y[4]) ~ -y[4] + 1e4 * y[2]^2]

    structural_simplify(ODESystem(sd3eqs, t; name=Symbol(gensym("sys"))))
end

sd3prob = ODEProblem{false}(
    sd3sys, [y[1:2] .=> 1; y[3:4] .=> 0.0], (0, 20.0), dt = 2.5e-5, cse = true; u0_constructor = x -> SVector(x...))

sd4sys = let
    sd4eqs = [D(y[1]) ~ -0.013y[1] - 1000 * (y[1] * y[3]),
        D(y[2]) ~ -2500 * (y[2] * y[3]),
        D(y[3]) ~ -0.013y[1] - 1000 * (y[1] * y[3]) - 2500 * (y[2] * y[3])]

    structural_simplify(ODESystem(sd4eqs, t; name=Symbol(gensym("sys"))))
end

sd4prob = ODEProblem{false}(
    sd4sys, [y[1:2] .=> 1; y[3] => 0.0], (0, 50.0), dt = 2.9e-4, cse = true; u0_constructor = x -> SVector(x...))

sd5sys = let
    sd5eqs = [D(y[1]) ~ 0.01 - (1 + (y[1] + 1000) * (y[1] + 1)) * (0.01 + y[1] + y[2]),
        D(y[2]) ~ 0.01 - (1 + y[2]^2) * (0.01 + y[1] + y[2])]

    structural_simplify(ODESystem(sd5eqs, t; name=Symbol(gensym("sys"))))
end

sd5prob = ODEProblem{false}(sd5sys, y[1:2] .=> 0.0, (0, 100.0), dt = 1e-4, cse = true; u0_constructor = x -> SVector(x...))

sd6sys = let
    sd6eqs = [D(y[1]) ~ -y[1] + 10^8 * y[3] * (1 - y[1]),
        D(y[2]) ~ -10y[2] + 3e7 * y[3] * (1 - y[2]),
        D(y[3]) ~ -(-y[1] + 10^8 * y[3] * (1 - y[1])) - (-10y[2] + 3e7 * y[3] * (1 - y[2]))
    ]

    structural_simplify(ODESystem(sd6eqs, t; name=Symbol(gensym("sys"))))
end

sd6prob = ODEProblem{false}(
    sd6sys, [y[1] => 1.0; y[2:3] .=> 0.0], (0, 1.0), dt = 3.3e-8, cse = true; u0_constructor = x -> SVector(x...))

se1sys = let
    Γ = 100
    se1eqs = [D(y[1]) ~ y[2],
        D(y[2]) ~ y[3],
        D(y[3]) ~ y[4],
        D(y[4]) ~ (y[1]^2 - sin(y[1]) - Γ^4) * y[1] +
                  (y[2] * y[3] / (y[1]^2 + 1) - 4 * Γ^3) * y[2] + (1 - 6 * Γ^2) * y[3] +
                  (10 * exp(-y[4]^2) - 4Γ) * y[4] + 1
    ]

    structural_simplify(ODESystem(se1eqs, t; name=Symbol(gensym("sys"))))
end

se1prob = ODEProblem{false}(se1sys, y .=> 0.0, (0, 1.0), dt = 6.8e-3, cse = true; u0_constructor = x -> SVector(x...))

se2sys = let
    se2eqs = [D(y[1]) ~ y[2],
        D(y[2]) ~ 5 * (1 - y[1]^2) * y[2] - y[1]
    ]

    structural_simplify(ODESystem(se2eqs, t; name=Symbol(gensym("sys"))))
end

se2prob = ODEProblem{false}(
    se2sys, [y[1] => 2.0, y[2] => 0.0], (0, 1.0), dt = 1e-3, cse = true; u0_constructor = x -> SVector(x...))

se3sys = let
    se3eqs = [D(y[1]) ~ -(55 + y[3]) * y[1] + 65 * y[2],
        D(y[2]) ~ 0.0785 * (y[1] - y[2]),
        D(y[3]) ~ 0.1 * y[1]
    ]

    structural_simplify(ODESystem(se3eqs, t; name=Symbol(gensym("sys"))))
end

se3prob = ODEProblem{false}(se3sys, [y[1:2] .=> 1.0; y[3] => 0.0], (0, 500.0), dt = 0.02; u0_constructor = x -> SVector(x...))

se4sys = let y = y
    y = y[1:4]
    U = ones(4, 4)
    U[diagind(U)] .= -1
    U ./= 2
    Z = U * y
    G = U' * [(Z[1]^2 - Z[2]^2) / 2, Z[1] * Z[2], Z[3]^2, Z[4]^2]
    A = [-10 -10 0 0; 10 -10 0 0; 0 0 1000 0; 0 0 0 0.01]
    se4eqs = D.(y) .~ -(U' * A * Z) + G

    structural_simplify(ODESystem(se4eqs, t; name=Symbol(gensym("sys"))))
end

se4prob = ODEProblem{false}(se4sys, [y[1] => 0.0; y[2] => -2.0; y[3:4] .=> -1.0],
    (0, 1000.0), dt = 1e-3, cse = true; u0_constructor = x -> SVector(x...))

se5sys = let
    se5eqs = [
        D(y[1]) ~ -7.89e-10 * y[1] - 1.1e7 * y[1] * y[3],
        D(y[2]) ~ 7.89e-10 * y[1] - 1.13e9 * y[2] * y[3],
        D(y[3]) ~ 7.89e-10 * y[1] - 1.1e7 * y[1] * y[3] + 1.13e3 * y[4] -
                  1.13e9 * y[2] * y[3],
        D(y[4]) ~ 1.1e7 * y[1] * y[3] - 1.13e3 * y[4]
    ]

    structural_simplify(ODESystem(se5eqs, t; name=Symbol(gensym("sys"))))
end

se5prob = ODEProblem{false}(
    se5sys, [y[1] => 1.76e-3; y[2:4] .=> 0.0], (0, 1000.0), dt = 1e-3, cse = true; u0_constructor = x -> SVector(x...))

sf1sys = let
    k = exp(20.7 - 1500 / y[1])
    sf1eqs = [
        D(y[1]) ~ 1.3 * (y[3] - y[1]) + 10400 * k * y[2],
        D(y[2]) ~ 1880 * [y[4] - y[2] * (1 + k)],
        D(y[3]) ~ 1752 - 269 * y[3] + 267 * y[1],
        D(y[4]) ~ 0.1 + 320 * y[2] - 321 * y[4]
    ]

    structural_simplify(ODESystem(sf1eqs, t; name=Symbol(gensym("sys"))))
end

sf1prob = ODEProblem{false}(
    sf1sys, [y[1] => 761.0; y[2] => 0.0; y[3] => 600.0; y[4] => 0.1],
    (0, 1000.0), dt = 1e-4, cse = true)

sf2sys = let
    sf2eqs = [
        D(y[1]) ~ -y[1] - y[1] * y[2] + 294 * y[2],
        D(y[2]) ~ y[1] * (1 - y[2]) / 98 - 3 * y[2]
    ]

    structural_simplify(ODESystem(sf2eqs, t; name=Symbol(gensym("sys"))))
end

sf2prob = ODEProblem{false}(
    sf2sys, [y[1] => 1.0; y[2] => 0.0], (0, 240.0), dt = 1e-2, cse = true; u0_constructor = x -> SVector(x...))

sf3sys = let
    sf3eqs = [
        D(y[1]) ~ -1e7 * y[2] * y[1] + 10 * y[3],
        D(y[2]) ~ -1e7 * y[2] * y[1] - 1e7 * y[2] * y[5] + 10 * y[3] + 10 * y[4],
        D(y[3]) ~ 1e7 * y[2] * y[1] - 1.001e4 * y[3] + 1e-3 * y[4],
        D(y[4]) ~ 1e4 * y[3] - 10.001 * y[4] + 1e7 * y[2] * y[5],
        D(y[5]) ~ 10 * y[4] - 1e7 * y[2] * y[5]
    ]

    structural_simplify(ODESystem(sf3eqs, t; name=Symbol(gensym("sys"))))
end

sf3prob = ODEProblem{false}(
    sf3sys, [y[1] => 4e-6; y[2] => 1e-6; y[3:5] .=> 0.0], (0, 100.0), dt = 1e-6, cse = true; u0_constructor = x -> SVector(x...))

sf4sys = let
    sf4eqs = [
        D(y[1]) ~ 77.27 * (y[2] - y[1] * y[2] + y[1] - 8.375e-6 * y[1]^2),
        D(y[2]) ~ (-y[2] - y[1] * y[2] + y[3]) / 77.27,
        D(y[3]) ~ 0.161 * (y[1] - y[3])
    ]

    structural_simplify(ODESystem(sf4eqs, t; name=Symbol(gensym("sys"))))
end

sf4prob = ODEProblem{false}(
    sf4sys, [y[1] => 4.0; y[2] => 1.1; y[3] => 4.0], (0, 300.0), dt = 1e-3, cse = true; u0_constructor = x -> SVector(x...))

sf5sys = let
    sf5eqs = [
        D(y[1]) ~ 1e11 * (-3 * y[1] * y[2] + 0.0012 * y[4] - 9 * y[1] * y[3]),
        D(y[2]) ~ -3e11 * y[1] * y[2] + 2e7 * y[4],
        D(y[3]) ~ 1e11 * (-9 * y[1] * y[3] + 0.001 * y[4]),
        D(y[4]) ~ 1e11 * (3 * y[1] * y[2] - 0.0012 * y[4] + 9 * y[1] * y[3])
    ]

    structural_simplify(ODESystem(sf5eqs, t; name=Symbol(gensym("sys"))))
end

sf5prob = ODEProblem{false}(
    sf5sys, [y[1] => 3.365e-7; y[2] => 8.261e-3; y[3] => 1.642e-3; y[4] => 9.38e-6],
    (0, 100.0), dt = 1e-7, cse = true; u0_constructor = x -> SVector(x...))

# Non-stiff
y1 = y[1]
na1eqs = D(y1) ~ -y1
@named na1sys_raw = ODESystem(na1eqs, t)
na1sys = structural_simplify(na1sys_raw)

na1prob = ODEProblem{false}(na1sys, [y[1] => 1], (0, 20.0), cse = true; u0_constructor = x -> SVector(x...))

y2 = y[1]
na2eqs = D(y2) ~ -y2^3 / 2
@named na2sys_raw = ODESystem(na2eqs, t)
na2sys = structural_simplify(na2sys_raw)

na2prob = ODEProblem{false}(na2sys, [y[1] => 1], (0, 20.0), cse = true; u0_constructor = x -> SVector(x...))

na3sys = let y = y[1]
    na3eqs = D(y) ~ y * cos(t)

    structural_simplify(ODESystem(na3eqs, t; name=Symbol(gensym("sys"))))
end

na3prob = ODEProblem{false}(na3sys, [y[1] => 1], (0, 20.0), cse = true; u0_constructor = x -> SVector(x...))

na4sys = let y = y[1]
    na4eqs = D(y) ~ y / 4 * (1 - y / 20)

    structural_simplify(ODESystem(na4eqs, t; name=Symbol(gensym("sys"))))
end

na4prob = ODEProblem{false}(na4sys, [y[1] => 1], (0, 20.0), cse = true; u0_constructor = x -> SVector(x...))

na5sys = let y = y[1]
    na5eqs = D(y) ~ (y - t) / (y + t)

    structural_simplify(ODESystem(na5eqs, t; name=Symbol(gensym("sys"))))
end

na5prob = ODEProblem{false}(na5sys, [y[1] => 4], (0, 20.0), cse = true; u0_constructor = x -> SVector(x...))

const NA_PROBLEMS = [na1prob, na2prob, na3prob, na4prob, na5prob]

nb1sys = let
    nb1eqs = [
        D(y[1]) ~ 2 * (y[1] - y[1] * y[2]),
        D(y[2]) ~ -(y[2] - y[1] * y[2])
    ]

    structural_simplify(ODESystem(nb1eqs, t; name=Symbol(gensym("sys"))))
end

nb1prob = ODEProblem{false}(nb1sys, [y[1] => 1.0, y[2] => 3], (0, 20.0), cse = true; u0_constructor = x -> SVector(x...))

nb2sys = let
    nb2eqs = [
        D(y[1]) ~ -y[1] + y[2],
        D(y[2]) ~ y[1] - 2y[2] + y[3],
        D(y[3]) ~ y[2] - y[3]
    ]

    structural_simplify(ODESystem(nb2eqs, t; name=Symbol(gensym("sys"))))
end

nb2prob = ODEProblem{false}(
    nb2sys, [y[1] => 2.0, y[2] => 0.0, y[3] => 1.0], (0, 20.0), cse = true; u0_constructor = x -> SVector(x...))

nb3sys = let
    nb3eqs = [
        D(y[1]) ~ -y[1],
        D(y[2]) ~ y[1] - y[2]^2,
        D(y[3]) ~ y[2]^2
    ]

    structural_simplify(ODESystem(nb3eqs, t; name=Symbol(gensym("sys"))))
end

nb3prob = ODEProblem{false}(nb3sys, [y[1] => 1.0; y[2:3] .=> 0.0], (0, 20.0), cse = true; u0_constructor = x -> SVector(x...))

nb4sys = let
    r = sqrt(y[1]^2 + y[2]^2)
    nb4eqs = [
        D(y[1]) ~ -y[2] - (y[1] * y[3]) / r,
        D(y[2]) ~ y[1] - (y[2] * y[3]) / r,
        D(y[3]) ~ y[1] / r
    ]

    structural_simplify(ODESystem(nb4eqs, t; name=Symbol(gensym("sys"))))
end

nb4prob = ODEProblem{false}(nb4sys, [y[1] => 3.0; y[2:3] .=> 0.0], (0, 20.0), cse = true; u0_constructor = x -> SVector(x...))

nb5sys = let
    nb5eqs = [
        D(y[1]) ~ y[2] * y[3],
        D(y[2]) ~ -y[1] * y[3],
        D(y[3]) ~ -0.51 * y[1] * y[2]
    ]

    structural_simplify(ODESystem(nb5eqs, t; name=Symbol(gensym("sys"))))
end

nb5prob = ODEProblem{false}(nb5sys, [y[1] => 0.0; y[2:3] .=> 1.0], (0, 20.0), cse = true; u0_constructor = x -> SVector(x...))

const NB_PROBLEMS = [nb1prob, nb2prob, nb3prob, nb4prob, nb5prob]

nc1sys = let y = y
    n = 10
    y = y[1:n]
    A = Bidiagonal(fill(-1, n), fill(1, n - 1), :L)
    nc1eqs = D.(y) .~ A * y
    structural_simplify(ODESystem(nc1eqs, t; name=Symbol(gensym("sys"))))
end

nc1prob = ODEProblem{false}(nc1sys, [y[1] => 1.0; y[2:10] .=> 0.0], (0, 20.0); u0_constructor = x -> SVector(x...))

nc2sys = let y = y
    n = 10
    y = y[1:n]
    A = Bidiagonal([-1:-1:(-n + 1); 0], collect(1:(n - 1)), :L)
    nc2eqs = D.(y) .~ A * y
    structural_simplify(ODESystem(nc2eqs, t; name=Symbol(gensym("sys"))))
end

nc2prob = ODEProblem{false}(nc2sys, [y[1] => 1.0; y[2:10] .=> 0.0], (0, 20.0), cse = true; u0_constructor = x -> SVector(x...))

nc3sys = let y = y
    n = 10
    y = y[1:n]
    A = SymTridiagonal(fill(-2, n), fill(1, n - 1))
    nc3eqs = D.(y) .~ A * y
    structural_simplify(ODESystem(nc3eqs, t; name=Symbol(gensym("sys"))))
end

nc3prob = ODEProblem{false}(nc3sys, [y[1] => 1.0; y[2:10] .=> 0.0], (0, 20.0), cse = true; u0_constructor = x -> SVector(x...))

@variables y(t)[1:51]
y = collect(y)
nc4sys = let y = y
    n = 51
    y = y[1:n]
    A = SymTridiagonal(fill(-2, n), fill(1, n - 1))
    nc4eqs = D.(y) .~ A * y
    structural_simplify(ODESystem(nc4eqs, t; name=Symbol(gensym("sys"))))
end

nc4prob = ODEProblem{false}(nc4sys, [y[1] => 1.0; y[2:51] .=> 0.0], (0, 20.0), cse = true; u0_constructor = x -> SVector(x...))

@variables y(t)[1:3, 1:5]
y = collect(y)
# TODO: nc5sys - complex N-body system needs special handling for MTK v10
# nc5sys = let
#     k_2 = 2.95912208286
#     m_0 = 1.00000597682
#     ms = [0.000954786104043
#           0.000285583733151
#           0.0000437273164546
#           0.0000517759138449
#           0.00000277777777778]
# 
#     r = [sqrt(sum(i -> y[i, j]^2, 1:3)) for j in 1:5]
#     d = [sqrt(sum(i -> (y[i, k] - y[i, j])^2, 1:3)) for k in 1:5, j in 1:5]
#     ssum(i, j) =
#         sum(1:5) do k
#             k == j && return 0
#             ms[k] * (y[i, k] - y[i, j]) / d[j, k]^3
#         end
#     nc5eqs = [D(D(y[i, j])) ~ k_2 * (-(m_0 + ms[j]) * y[i, j]) / r[j]^3 + ssum(i, j)
#               for i in 1:3, j in 1:5]
#     structural_simplify(ODESystem(nc5eqs, t))
# end

# Placeholder for nc5 system - needs MTK v10 compatible implementation  
@variables temp_y(t)
@named nc5_placeholder = ODESystem([D(temp_y) ~ 0], t)
nc5sys = structural_simplify(nc5_placeholder)

ys = [3.42947415189, 3.35386959711, 1.35494901715,
    6.64145542550, 5.97156957878, 2.18231499728,
    11.2630437207, 14.6952576794, 6.27960525067,
    -30.1552268759, 1.65699966404, 1.43785752721,
    -21.1238353380, 28.4465098142, 15.3882659679]
ys′ = [-0.557160570446, 0.505696783289, 0.230578543901,
    -0.415570776342, 0.365682722812, 0.169143213293,
    -0.325325669158, 0.189706021964, 0.087726532278,
    -0.0240476254170, -0.287659532608, -0.117219543175,
    -0.176860753121, -0.216393453025, -0.0148647893090]
y0 = y .=> reshape(ys, 3, 5)
y0′ = D.(y) .=> reshape(ys′, 3, 5)
# The orginal paper has t_f = 20, but 1000 looks way better
# nc5prob = ODEProblem{false}(nc5sys, [y0; y0′], (0, 20.0), cse = true; u0_constructor = x -> SVector(x...))

const NC_PROBLEMS = [nc1prob, nc2prob, nc3prob, nc4prob] # nc5prob temporarily disabled for MTK v10

@variables y(t)[1:4]
y = collect(y)
@parameters ε
# Use intermediate variable to avoid complex symbolic expansion
# r_cubed = (x^2 + y^2)^(3/2) for the gravitational force law
r_squared = y[1]^2 + y[2]^2
r_cubed = r_squared * sqrt(r_squared)
nd1eqs = [D(y[1]) ~ y[3],
    D(y[2]) ~ y[4],
    D(y[3]) ~ (-y[1]) / r_cubed,
    D(y[4]) ~ (-y[2]) / r_cubed]
@named nd1sys_raw = ODESystem(nd1eqs, t)
nd1sys = structural_simplify(nd1sys_raw)

function make_ds(nd1sys, e)
    y = collect(@nonamespace nd1sys.y)
    y0 = [y[1] => 1 - e; y[2:3] .=> 0.0; y[4] => sqrt((1 + e) / (1 - e))]
    ODEProblem{false}(nd1sys, y0, (0, 20.0), [ε => e])
end
nd1prob = make_ds(nd1sys, 0.1)
nd2prob = make_ds(nd1sys, 0.3)
nd3prob = make_ds(nd1sys, 0.5)
nd4prob = make_ds(nd1sys, 0.7)
nd5prob = make_ds(nd1sys, 0.9)

const ND_PROBLEMS = [nd1prob, nd2prob, nd3prob, nd4prob, nd5prob]

ne1sys = let
    ne1eqs = [D(y[1]) ~ y[2],
        D(y[2]) ~ -(y[2] / (t + 1) + (1 - 0.25 / (t + 1)^2) * y[1])
    ]
    structural_simplify(ODESystem(ne1eqs, t; name=Symbol(gensym("sys"))))
end

y0 = [y[1] => 0.6713967071418030; y[2] => 0.09540051444747446]
ne1prob = ODEProblem{false}(ne1sys, y0, (0, 20.0); u0_constructor = x -> SVector(x...))

ne2sys = let
    ne2eqs = [D(y[1]) ~ y[2],
        D(y[2]) ~ (1 - y[1]^2) * y[2] - y[1]
    ]
    structural_simplify(ODESystem(ne2eqs, t; name=Symbol(gensym("sys"))))
end

y0 = [y[1] => 2.0; y[2] => 0.0]
ne2prob = ODEProblem{false}(ne2sys, y0, (0, 20.0), cse = true)

ne3sys = let
    ne3eqs = [D(y[1]) ~ y[2],
        D(y[2]) ~ y[1]^3 / 6 - y[1] + 2 * sin(2.78535t)
    ]
    structural_simplify(ODESystem(ne3eqs, t; name=Symbol(gensym("sys"))))
end

ne3prob = ODEProblem{false}(ne3sys, y[1:2] .=> 0, (0, 20.0); u0_constructor = x -> SVector(x...))

ne4sys = let
    ne4eqs = [D(y[1]) ~ y[2],
        D(y[2]) ~ 0.032 - 0.4 * y[2]^2
    ]
    structural_simplify(ODESystem(ne4eqs, t; name=Symbol(gensym("sys"))))
end

ne4prob = ODEProblem{false}(ne4sys, [y[1] => 30.0, y[2] => 0.0], (0, 20.0); u0_constructor = x -> SVector(x...))

ne5sys = let
    ne5eqs = [D(y[1]) ~ y[2],
        D(y[2]) ~ sqrt(1 + y[2]^2) / (25 - t)]
    structural_simplify(ODESystem(ne5eqs, t; name=Symbol(gensym("sys"))))
end

ne5prob = ODEProblem{false}(ne5sys, y[1:2] .=> 0.0, (0, 20.0); u0_constructor = x -> SVector(x...))

const NE_PROBLEMS = [ne1prob, ne2prob, ne3prob, ne4prob, ne5prob]

@inline myifelse(x, y, z) = ifelse(x, y, z)
nf1sys = let
    a = 0.1
    cond = term(iseven, term(floor, Int, unwrap(t), type = Int), type = Bool)
    b = 2a * y[2] - (pi^2 + a^2) * y[1]
    nf1eqs = [D(y[1]) ~ y[2],
        D(y[2]) ~ b + term(myifelse, cond, 1, -1, type = Real)]
    structural_simplify(ODESystem(nf1eqs, t; name=Symbol(gensym("sys"))))
end

nf1prob = ODEProblem{false}(nf1sys, y[1:2] .=> 0.0, (0, 20.0);u0_constructor = x -> SVector(x...))

# nf2sys = complete(let
#     cond = term(iseven, term(floor, Int, unwrap(t), type = Int), type = Bool)
#     nf2eqs = [D(y[1]) ~ 55 - term(myifelse, cond, 2y[1] / 2, y[1] / 2, type = Real)]
#     ODESystem(nf2eqs, t, name = :nf2)
# end)

# nf2prob = ODEProblem{false}(nf2sys, [y[1] .=> 110.0], (0, 20.0); u0_constructor = x -> SVector(x...))

nf3sys = let
    nf3eqs = [D(y[1]) ~ y[2],
        D(y[2]) ~ 0.01 * y[2] * (1 - y[1]^2) - y[1] - abs(sin(pi * t))]
    structural_simplify(ODESystem(nf3eqs, t; name=Symbol(gensym("sys"))))
end

nf3prob = ODEProblem{false}(nf3sys, y[1:2] .=> 0.0, (0, 20.0); u0_constructor = x -> SVector(x...))

nf4sys = let
    nf4eqs = [D(y[1]) ~ term(
        myifelse, t <= 10, -2 / 21 - (120 * (t - 5)) / (1 + 4 * (t - 5)^2),
        -2y[1], type = Real)]
    structural_simplify(ODESystem(nf4eqs, t; name=Symbol(gensym("sys"))))
end

nf4prob = ODEProblem{false}(nf4sys, [y[1] => 1.0], (0, 20.0); u0_constructor = x -> SVector(x...))

nf5sys = let
    c = sum(i -> cbrt(i)^4, 1:19)
    p = sum(i -> cbrt(t - i)^4, 1:19)
    nf5eqs = [D(y[1]) ~ inv(c) * Symbolics.derivative(p, t) * y[1]]
    structural_simplify(ODESystem(nf5eqs, t; name=Symbol(gensym("sys"))))
end

nf5prob = ODEProblem{false}(nf5sys, [y[1] => 1.0], (0, 20.0); u0_constructor = x -> SVector(x...))


const NF_PROBLEMS = [nf1prob, nf3prob, nf4prob, nf5prob]
