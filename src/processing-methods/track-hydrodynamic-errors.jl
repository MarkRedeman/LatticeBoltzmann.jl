struct TrackHydrodynamicErrors{T} <: ProcessingMethod
    problem::FluidFlowProblem
    should_process::Bool
    n_steps::Int64
    stop_criteria::StopCriteria
    df::T
end
TrackHydrodynamicErrors(problem, should_process, n_steps) = TrackHydrodynamicErrors(
    problem,
    should_process,
    n_steps,
    StopCriteria(problem),
    Vector{NamedTuple{
        (
            :timestep,
            :error_ρ,
            :error_u,
            :error_p,
            :error_σ_xx,
            :error_σ_xy,
            :error_σ_yy,
            :error_σ_yx,
        ),
        Tuple{
            Int64,
            Float64,
            Float64,
            Float64,

            Float64,
            Float64,
            Float64,
            Float64,
        }
    }}()
)


function next!(process_method::TrackHydrodynamicErrors, q, f_in, t::Int64)
    if (! process_method.should_process)
        if (t != process_method.n_steps)
            return false
        end
    end

    problem = process_method.problem
    nx, ny, nf = size(f_in)
    x_range, y_range = range(problem)

    f = Array{Float64}(undef, size(f_in, 3))
    u = zeros(dimension(q))
    expected_u = zeros(dimension(q))

    error_ρ = 0.0
    error_u = 0.0
    error_p = 0.0
    error_σ_xx = 0.0
    error_σ_xy = 0.0
    error_σ_yy = 0.0
    error_σ_yx = 0.0

    time = t * delta_t(problem)
    Δ = Float64(y_range.step) * Float64(x_range.step)

    τ = q.speed_of_sound_squared * lbm.lattice_viscosity(problem) + 0.5
    # @show q.speed_of_sound_squared * lbm.lattice_viscosity(problem) + 0.5, problem.ν, 1/problem.ν, viscosity(problem)

    # τ = q.speed_of_sound_squared * lbm.lattice_viscosity(problem)
    D = dimension(q)
    N = div(lbm.order(q), 2)
    N = 2
    Hs = [
        [
            hermite(Val{n}, q.abscissae[:, i], q)
            for i = 1:length(q.weights)
        ]
        for n = 1:N
    ]

    @show problem.ν * delta_x(problem)^2 / delta_t(problem)
    @inbounds for x_idx in 1:nx, y_idx in 1:ny
        @inbounds for f_idx = 1 : size(f_in, 3)
            f[f_idx] = f_in[x_idx, y_idx, f_idx]
        end

        ρ = density(q, f)
        velocity!(q, f, ρ, u)
        p = pressure(q, f, ρ, u)

        # Adding the forcing term moves the optimal tau for poiseuille flows
        # F = cm.force(x_idx, y_idx, 0.0)
        # u += cm.τ * F

        ρ = dimensionless_density(problem, ρ)
        u = dimensionless_velocity(problem, u)
        p = dimensionless_pressure(q, problem, p)

        # Analytical
        x = x_range[x_idx]
        y = y_range[y_idx]

        # Compute statistics
        expected_ρ = lbm.density(q, problem, x, y, time)
        expected_u = lbm.velocity(problem, x, y, time)

        expected_p = lbm.pressure(q, problem, x, y, time)
        expected_ϵ = (dimension(q) / 2) * expected_p / expected_ρ
        expected_T = expected_p / expected_ρ

        error_ρ += Δ * (ρ - expected_ρ)^2
        # @show p expected_p
        # error_p += Δ * (p - expected_p)^2
        error_u += Δ * ((u[1] - expected_u[1])^2 + (u[2] - expected_u[2])^2)


        # if (x_idx == div(nx, 2) && y_idx == 1)
            a_f = [
                sum([f[idx] * Hs[n][idx] for idx = 1:length(q.weights)])
                for n = 1:N
            ]

            P = a_f[2] - (a_f[1] * a_f[1]') / ρ - ρ * I(2)
            P = P * (1 - 1 / (2 * τ))
            # σ_lb = momentum_flux(q, f, ρ, u) - I(2) * p
            σ_lb = (P - I(2) * tr(P) / (D))

            # Rescale to dimensionless number (TODO check why problem.u_max)
            σ_lb = σ_lb / (problem.u_max)

            # P = a_f[2] - (a_f[1] * a_f[1]') / ρ - ρ * I(2)
            P = (1 - 1 / (2 * τ))a_f[2] - (a_f[1] * a_f[1]') / ρ - ρ * I(2)
            # P = P * (1 - 1 / (2 * τ))
            # @show P

        # @show tr(P) / D

        error_p += Δ * ((-tr(P)/D) - expected_p)^2

            σ_lb = (P - I(2) * tr(P) / (D)) / (problem.u_max)

            a_eq_2 = equilibrium_coefficient(Val{2}, q, ρ, a_f[1] ./ ρ, 1.0)
            σ_a = (a_f[2] - a_eq_2) / (1 - 1 / (2 * τ))

            a_2 = (a_f[2] + (1 / (2 * τ)) * a_eq_2) / (1 + 1 / (2 * τ))

        # Vorige poging
        # @show "SIGMAS"
            σ_lb = (a_2 - (a_f[1] * a_f[1]') / ρ) / problem.u_max
            # @show σ_lb

        # Huidige poging
            P = a_2 - (a_f[1] * a_f[1]') / ρ - ρ * I
            σ_lb = (P - (1/D) * tr(P) * I)
            # @show σ_lb
            # @show P -(P - (1/D) * tr(P) * I)

            # @show (a_f[1] * a_f[1]') / ρ
            # @show σ_a
            # @show (σ_a - (a_f[1] * a_f[1]') / ρ)
            # @show (σ_a - (a_f[1] * a_f[1]') / ρ) / problem.u_max
            # @show a_f[2] a_eq_2 σ_a / problem.u_max
            # @show tr(σ_a)

            factor = 0.9549296889427464
        # factor = q.speed_of_sound_squared / pi
            factor = (1.0 / delta_x(problem)) / 16
            factor = (1.0 / (delta_t(problem)))
        # factor *= 6/8
        # factor *= 6/8
        factor *= 3/4
            σ_lb *= (factor)
            cs = 1 / (q.speed_of_sound_squared)
            σ_exact = deviatoric_tensor(q, problem, x, y, time)
            # σ_exact /= (problem.domain_size[1] / delta_x(problem)) / 16

        # P = (
        #     a_f[2] - ρ * (u * u' - I)
        # )
        # @show P P /(1 + 1 / (2 * τ)) P / problem.u_max P / (problem.u_max * (1 + 1 / (2 * τ))) σ_lb σ_exact

            σ_err = (σ_exact .- σ_lb)
        # @show σ_err
            # @show x y
            # @show σ_lb σ_exact σ_err
            # @show σ_exact ./ σ_lb
            # @show σ_lb ./ σ_exact
            # @show σ_exact[1,2] σ_lb[1,2]
            # @show Δ *(σ_exact .- σ_lb).^2
        if (x_idx == div(nx, 2) && y_idx == 1)
            @show factor
            factor =  σ_exact[1,2] ./ σ_lb[1,2]
            @show factor
            factor =  σ_exact[2,2] ./ σ_lb[2,2]
            @show factor
            @show σ_err
            # @show τ q.speed_of_sound_squared
            # @show pi * τ, pi * q.speed_of_sound_squared, factor
            # @show (1.0 / delta_x(problem)) / 16
            # @show factor * q.speed_of_sound_squared
            # @show factor / q.speed_of_sound_squared
            # @show factor * pi
            # @show factor / pi
            # @show factor * τ
            # @show factor / τ
            # @show factor * problem.ν
            # @show factor / problem.ν

            # factor =  σ_exact[1,2] ./ σ_lb[1,2]
            factor =  σ_exact[2,2] ./ σ_lb[2,2]
            # @show factor σ_lb σ_exact (σ_exact .- σ_lb).^2
            # @show factor * q.speed_of_sound_squared
            # @show factor / q.speed_of_sound_squared
            # @show factor * pi
            # @show factor / pi
            # @show factor * τ
            # @show factor / τ
            # @show factor * problem.ν
            # @show factor / problem.ν
            end

            error_σ_xx += Δ * σ_err[1, 1]^2
            error_σ_xy += Δ * σ_err[1, 2]^2
            error_σ_yx += Δ * σ_err[2, 1]^2
            error_σ_yy += Δ * σ_err[2, 2]^2
        # end
    end

    @show error_σ_xx
    @show error_σ_xy
    @show sqrt(error_σ_xx)
    @show sqrt(error_σ_xy)

    # @show error_σ_xy
    push!(process_method.df, (
        timestep = t,
        error_ρ = sqrt(error_ρ),
        error_u = sqrt(error_u),
        error_p = sqrt(error_p),
        error_σ_xx = sqrt(error_σ_xx),
        error_σ_xy = sqrt(error_σ_xy),
        error_σ_yy = sqrt(error_σ_yy),
        error_σ_yx = sqrt(error_σ_yx),
    ))

    if mod(t, 100) == 0
        if (should_stop!(process_method.stop_criteria, q, f_in))
            @info "Stopping after $t steps out of $process_method.n_steps"
            return true
        end
    end

    if mod(t, 1) == 0
        should_visualize = false
        if (process_method.should_process)
            if t == process_method.n_steps
                should_visualize = true
            end

            if mod(t, max(10, round(Int, process_method.n_steps / 5))) == 0
                should_visualize = true
            end
        end


        if (should_visualize)
            Δt = delta_t(process_method.problem)
            visualize(
                process_method.problem,
                q,
                f_in,
                time,
                process_method.df
            )
        end
    end

    return false
end
