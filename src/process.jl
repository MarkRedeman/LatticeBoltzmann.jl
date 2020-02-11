abstract type ProcessingMethod end
ProcessingMethod(problem, should_process, n_steps, stop_criteria = StopCriteria(problem)) =
    CompareWithAnalyticalSolution(problem, should_process, n_steps, stop_criteria)
ProcessingMethod(problem::TaylorGreenVortex, should_process, n_steps, stop_criteria = StopCriteria(problem)) =
    TrackHydrodynamicErrors(problem, should_process, n_steps, stop_criteria)
ProcessingMethod(problem::DecayingShearFlow, should_process, n_steps, stop_criteria = StopCriteria(problem)) =
    TrackHydrodynamicErrors(problem, should_process, n_steps, stop_criteria)
ProcessingMethod(problem::TGV, should_process, n_steps, stop_criteria = StopCriteria(problem)) =
    TrackHydrodynamicErrors(problem, should_process, n_steps, stop_criteria)
# ProcessingMethod(problem::TGV, should_process, n_steps) =
#     TrackHydrodynamicErrors(problem, should_process, n_steps)

include("processing_methods/track_hydrodynamic_errors.jl")
include("processing_methods/visualize.jl")
include("processing_methods/take_snapshots.jl")
struct CompareWithAnalyticalSolution{T} <: ProcessingMethod
    problem::FluidFlowProblem
    should_process::Bool
    n_steps::Int64
    stop_criteria::StopCriteria
    df::T
end
CompareWithAnalyticalSolution(
    problem,
    should_process,
    n_steps,
    stop_criteria = StopCriteria(problem)
) = CompareWithAnalyticalSolution(
        problem,
        should_process,
        n_steps,
        stop_criteria,
        Vector{
            NamedTuple{
                (
                    :density,
                    :momentum,
                    :total_energy,
                    :kinetic_energy,
                    :internal_energy,
                    :density_a,
                    :momentum_a,
                    :total_energy_a,
                    :kinetic_energy_a,
                    :internal_energy_a,
                    :error_u,
                    :error_p,
                    :error_σ_xx,
                    :error_σ_xy,
                    :error_σ_yy,
                    :error_σ_yx,
                ),
                Tuple{
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    Float64,

                    Float64,
                    Float64,
                    Float64,
                    Float64,
                },
            },
        }(),
    )

function next!(process_method::CompareWithAnalyticalSolution, q, f_in, t::Int64)
    if mod(t, 100) == 0
        if (should_stop!(process_method.stop_criteria, q, f_in))
            @info "Stopping after $t steps "

            Δt = delta_t(process_method.problem)
            process!(
                process_method.problem,
                q,
                f_in,
                t * Δt,
                process_method.df;
                should_visualize = false,
            )
            return true
        end
    end

    if (!process_method.should_process)
        if (t != process_method.n_steps)
            return false
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


        Δt = delta_t(process_method.problem)
        process!(
            process_method.problem,
            q,
            f_in,
            t * Δt,
            process_method.df;
            should_visualize = should_visualize,
        )
    end

    return false
end

function process!(
    problem::FluidFlowProblem,
    q::Quadrature,
    f_in,
    time,
    stats;
    should_visualize = false,
)
    f = Array{Float64}(undef, size(f_in, 3))
    u = zeros(dimension(q))
    expected_u = zeros(dimension(q))

    Nx = size(f_in, 1)
    Ny = size(f_in, 2)

    x_range, y_range = range(problem)

    total_density = 0.0
    total_momentum = 0.0
    total_energy = 0.0
    total_kinetic_energy = 0.0
    total_internal_energy = 0.0

    expected_total_density = 0.0
    expected_total_momentum = 0.0
    expected_total_energy = 0.0
    expected_total_kinetic_energy = 0.0
    expected_total_internal_energy = 0.0

    rho_error_squared = 0.0
    ux_error_squared = 0.0
    uy_error_squared = 0.0
    error_u = 0.0
    error_p = 0.0

    @inbounds for x_idx = 1:Nx, y_idx = 1:Ny
        # Calculated
        @inbounds for f_idx = 1:size(f_in, 3)
            f[f_idx] = f_in[x_idx, y_idx, f_idx]
        end
        ρ = density(q, f)
        velocity!(q, f, ρ, u)

        # Adding the forcing term moves the optimal tau for poiseuille flows
        # F = cm.force(x_idx, y_idx, 0.0)
        # u += cm.τ * F

        T = temperature(q, f, ρ, u)
        p = pressure(q, f, ρ, u)

        ρ = dimensionless_density(problem, ρ)
        u = dimensionless_velocity(problem, u)
        T = dimensionless_temperature(q, problem, T)
        p = dimensionless_pressure(q, problem, p)

        opp = 1.0/(Nx * Ny)
        # CHECK!
        opp = 1.0
        total_density += opp * ρ
        total_momentum += opp * (u[1] + u[2]) * ρ
        kinetic_energy = opp * (u[1]^2 + u[2]^2) * ρ

        internal_energy = T

        total_kinetic_energy += opp * kinetic_energy
        total_internal_energy += opp * internal_energy
        total_energy += opp * (kinetic_energy + internal_energy)

        # Analytical
        x = x_range[x_idx]
        y = y_range[y_idx]

        # Compute statistics
        expected_ρ = LatticeBoltzmann.density(q, problem, x, y, time)
        expected_p = LatticeBoltzmann.pressure(q, problem, x, y, time)
        expected_ϵ = (dimension(q) / 2) * expected_p / expected_ρ
        expected_v = LatticeBoltzmann.velocity(problem, x, y, time)
        expected_T = expected_p / expected_ρ

        expected_kinetic_energy = (expected_v[1]^2 + expected_v[2]^2)
        expected_internal_energy = expected_T

        expected_total_density += expected_ρ
        expected_total_momentum += expected_ρ * (expected_v[1] + expected_v[2])
        expected_total_kinetic_energy += expected_kinetic_energy
        expected_total_internal_energy += expected_internal_energy
        expected_total_energy += expected_kinetic_energy + expected_internal_energy

        rho_error_squared += (ρ - expected_ρ)^2
        ux_error_squared += (u[1] - expected_v[1])^2
        uy_error_squared += (u[2] - expected_v[2])^2
        # error_u += sqrt((u[1] - expected_v[1])^2 + (u[2] - expected_v[2])^2)
        opp = Float64(y_range.step) * Float64(x_range.step)
        # opp = 1.0
        # opp = Float64(y_range.step)
        error_u += (opp * (((u[1] - expected_v[1]))^2 + ((u[2] - expected_v[2]))^2))
        error_p += opp * ((p - expected_p)^2)
    end

    # Compare with analytical results?
    push!(
        stats,
        (
            density = total_density,
            momentum = total_momentum,
            total_energy = total_energy,
            kinetic_energy = total_kinetic_energy,
            internal_energy = total_internal_energy,
            density_a = expected_total_density,
            momentum_a = expected_total_momentum,
            total_energy_a = expected_total_energy,
            kinetic_energy_a = expected_total_kinetic_energy,
            internal_energy_a = expected_total_internal_energy,
            error_u = sqrt(error_u),
            error_p = sqrt(error_p),

            error_σ_xx = 0.0,
            error_σ_xy = 0.0,
            error_σ_yy = 0.0,
            error_σ_yx = 0.0,
        ),
    )

    if should_visualize
        visualize(problem, q, f_in, time, stats)
    end

    return false
end

struct ShowVelocityError{V} <: ProcessingMethod
    problem::FluidFlowProblem
    plot_every::Int
    # TODO: change to vector of vector?
    l2_norms::Vector{V}
    stop_criteria::StopCriteria
end
ShowVelocityError(problem, plot_every = 1, l2_norms = Float64[]) = ShowVelocityError(problem, plot_every, l2_norms, NoStoppingCriteria())

function next!(process_method::ShowVelocityError, q, f, t::Int64)
    should_stop = false
    if mod(t, 100) == 0
        if (should_stop!(process_method.stop_criteria, q, f))
            @info "Stopping after $t iterations"
            should_stop = true
        end
    end
    # if mod(t, process_method.plot_every) != 1
    #     return false
    # end

    problem = process_method.problem
    analytical_velocity = (y) -> velocity(problem, 1.0, y)
    x_domain = (0.0, problem.domain_size[1])
    y_domain = (0.0, problem.domain_size[2])

    x_idx = 1

    v_y = zeros(problem.NY)
    v_e = zeros(problem.NY)
    v_a = zeros(problem.NY)
    u = zeros(dimension(q))
    x_range, y_range = range(problem)
    total_expected_momentum = 0.0
    for y_idx = 1 : problem.NY
        ρ = density(q, f[x_idx, y_idx, :])
        velocity!(q, f[x_idx, y_idx, :], ρ, u)
        u = dimensionless_velocity(problem, u)
        # v_e[y_idx] = norm(u - analytical_velocity(y_range[y_idx]))
        v_e[y_idx] = ((u[1] - analytical_velocity(y_range[y_idx])[1])^2 + (u[1] - analytical_velocity(y_range[y_idx])[1])^2)
        v_y[y_idx] = u[1]
        v_a[y_idx] = analytical_velocity(y_range[y_idx])[1]

        total_expected_momentum += analytical_velocity(y_range[y_idx])[1]^2 + analytical_velocity(y_range[y_idx])[2]^2
    end

    push!(
        process_method.l2_norms,
        sqrt(sum(v_e) / total_expected_momentum)
    )
    # push!(process_method.l2_norms, v_e)

    # return false
    if mod(t, process_method.plot_every) != 1 && ! should_stop
        return should_stop
    end

    y_idx = problem.NY
    ρ = density(q, f[x_idx, y_idx, :])
    velocity!(q, f[x_idx, y_idx, :], ρ, u)
    # @show v_y[problem.NY], v_a[problem.NY], v_e[problem.NY], u[1], ρ

    velocity_plot = plot(y_range, v_y, label = "Numerical solution", legend=:topleft, title = string(q))
    plot!(velocity_plot, y_range, v_a, label = "Exact solution")

    plot(
        scatter(y_range, v_e, title = "Abs. Error", legend=nothing),
        velocity_plot
    )
    gui()


    return should_stop
end
