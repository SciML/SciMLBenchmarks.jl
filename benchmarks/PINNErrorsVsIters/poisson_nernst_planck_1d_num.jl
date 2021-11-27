################################################################################
#   This code solves a 1D Poisson-Nernst-Planck system using finite difference
#
#   Equations:
#        d2φ/dx2 = -F/epsilon * ( z_H * H + z_OH * OH)
#        dH/dt =  D_H * d2H/dx2
#                 + c_H * ( dH/dx * dφ/dx + H * d2φ/dx2 )
#                 + k_wb * H2O - k_wf * H * OH
#        dOH/dt = D_OH * d2OH/dx2
#                 + c_OH * ( dOH/dx * dφ/dx + OH * d2φ/dx2 )
#                 + k_wb * H2O - k_wf * H * OH
#
#   Initial conditions:
#        H(0,x) = H_0
#        OH(0,x) = OH_0
#        φ(0,x) = 0
#
#   Boundary conditions:
#        H(t,0) = H_anode_rate * t + H_0
#        dH(t,n)/dx = 0
#        dOH(t,0)/dx = 0
#        OH(t,n) = OH_cathode_rate * t + OH_0
#        φ(t,0) = φ_0
#        φ(t,n) = 0
#
#   How to run:
#        julia poisson-nernst-planck-finite-difference.jl
#
################################################################################

function nernst_planck_1d_num()
    # Maximum time
    time_max = 50.0 # s

    # Delta t [s]
    dt = 5e-9 # s

    # Time iterations
    time_it = floor(Int,time_max/dt)

    # Domain longitud
    x_max = 0.01 # m

    # Delta x
    dx = 2e-4 # m

    # No. of spatial domain nodes
    n = floor(Int,x_max/dx)

    # No. of convergence iterations
    it_max = 1000

    # No. of iterations to print results
    it_print = floor(1.0/dt/20.0)

    # Electrochemical parameters
    H_anode_rate = 4.2 # mol m^-3 s^-1
    OH_cathode_rate = 5.0 # mol m^-3 s^-1
    H_0 = 1e-4 # mol m^-3
    OH_0 = 1e-4 # mol m^-3
    φ_0 = 4.0 # V
    z_H = 1.0
    D_H = 9.31e-9 # m^2 s^-1
    z_OH = -1.0
    D_OH = 5.26e-9 # m^2 s^−1
    F = 96485.3415 # A s mol^-1
    R = 8.31 # kg m^2 K^-1 mol^-1 s^-2
    T = 298.0 # K
    σ_0 = 0.2 # S m^-1
    k_wf = 1.5e8 # m^3 mol^-1 s^-1
    k_wb = 2.7e-5 # s^-1
    H2O = 55500.0 # mol m^-3
    c_H = z_H * D_H * F / (R * T)
    c_OH = z_OH * D_OH * F / (R * T)
    epsilon = 78.5 # K

    # Set initial conditions
    H = zeros(n)
    OH = zeros(n)
    φ = zeros(n)
    H_aux = zeros(n)
    OH_aux = zeros(n)

    # Add boundary conditions
    H .= H_0
    OH .= OH_0
    φ[1] = φ_0
    H_aux .= H_0
    OH_aux .= OH_0

    # Solve Poisson-Nernst-Planck equations
    it = 1
    while it < time_it

        # One each iteration two time steps are computed to avoid copying variables

        # First time step

        # Solve Nernst-Planck equation
        for i = 2:n-1
            H_aux[i] =  (D_H * (H[i+1]-2*H[i]+H[i-1])/(dx^2)
                        + c_H * (H[i+1]-H[i])/(dx) * (φ[i+1]-φ[i-1])/(2*dx)
                        + c_H * H[i] * (φ[i+1]-2*φ[i]+φ[i-1])/(dx^2)
                        + k_wb * H2O - k_wf * H[i] * OH[i]
                        ) * dt + H[i]
            OH_aux[i] =  (D_OH * (OH[i+1]-2*OH[i]+OH[i-1])/(dx^2)
                         + c_OH * (OH[i]-OH[i-1])/(dx) * (φ[i+1]-φ[i-1])/(2*dx)
                         + c_OH * OH[i] * (φ[i+1]-2*φ[i]+φ[i-1])/(dx^2)
                         + k_wb * H2O - k_wf * H[i] * OH[i]
                         ) * dt + OH[i]
        end
        H_aux[1] = H_anode_rate * it * dt + H_0
        OH_aux[n] = OH_cathode_rate * it * dt + OH_0
        H_aux[n] = H_aux[n-1]
        OH_aux[1] = OH_aux[2]

        # Solve Poisson equation
        for k = 1:it_max
            for i = 2:n-1
                φ[i]  = ( φ[i+1] + φ[i-1] ) / 2.0 + (F * dx^2) / (2.0 * epsilon) * ( z_H * H_aux[i] + z_OH * OH_aux[i])
            end
        end

        # Print results
        if it % it_print == 0
            println(stdout, "time:", it*dt, " s")

            println("φ:")
            println(φ)

            println("H+:")
            println(H)

            println("OH-:")
            println(OH)

            u_numeric = [φ, H, OH]
            return u_numeric

            flush(stdout)
        end

        global it += 1

        # Second time step

        # Solve Nernst-Planck equation
        for i = 2:n-1
            H[i] =  (D_H * (H_aux[i+1]-2*H_aux[i]+H_aux[i-1])/(dx^2)
                     + c_H * (H_aux[i+1]-H_aux[i])/(dx) * (φ[i+1]-φ[i-1])/(2*dx)
                     + c_H * H_aux[i] * (φ[i+1]-2*φ[i]+φ[i-1])/(dx^2)
                     + k_wb * H2O - k_wf * H_aux[i] * OH_aux[i]
                     ) * dt + H_aux[i]
            OH[i] =  (D_OH * (OH_aux[i+1]-2*OH_aux[i]+OH_aux[i-1])/(dx^2)
                      + c_OH * (OH_aux[i]-OH_aux[i-1])/(dx) * (φ[i+1]-φ[i-1])/(2*dx)
                      + c_OH * OH_aux[i] * (φ[i+1]-2*φ[i]+φ[i-1])/(dx^2)
                      + k_wb * H2O - k_wf * H_aux[i] * OH_aux[i]
                      ) * dt + OH_aux[i]
        end
        H[1] = H_anode_rate * it * dt + H_0
        OH[n] = OH_cathode_rate * it * dt + OH_0
        H[n] = H[n-1]
        OH[1] = OH[2]

        # Solve Poisson equation
        for k = 1:it_max
            for i = 2:n-1
                φ[i]  = ( φ[i+1] + φ[i-1] ) / 2.0 + (F * dx^2) / (2.0 * epsilon) * ( z_H * H[i] + z_OH * OH[i])
            end
        end

        # Print results
        if it % it_print == 0
            println(stdout, "time:", it*dt, " s")

            println("φ:")
            println(φ)

            println("H+:")
            println(H)

            println("OH-:")
            println(OH)

            u_numeric = [φ, H, OH]
            return u_numeric

            flush(stdout)
        end

        global it += 1

    end
end
