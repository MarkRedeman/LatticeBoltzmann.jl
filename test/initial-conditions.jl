import Base: OneTo
import LatticeBoltzmann:
    deviatoric_tensor,
    velocity!,
    temperature,
    ConstantDensity,
    AnalyticalVelocityAndStress,
    AnalyticalEquilibrium,
    AnalyticalEquilibriumAndOffEquilibrium,
    IterativeInitializationMeiEtAl

struct InitializationTestProblem <: FluidFlowProblem
    u_max::Float64
    ν::Float64
    NX::Int64
    NY::Int64
    domain_size::Tuple{Float64,Float64}
end
InitializationTestProblem(N = 4, ν = 1.0 / 6.0) = InitializationTestProblem(
    1.0,
    ν,
    N,
    N,
    # (2.0 * pi, 2.0 * pi)
    (1.0, 1.0)
)

# const U = 0.01
const U = 0.001
# const U = 1.00
LatticeBoltzmann.density(q::Quadrature, problem::InitializationTestProblem, x::T, y::T, timestep::Real = 0.0) where { T <: Real }= pressure(q, problem, x, y, timestep)
function LatticeBoltzmann.pressure(q::Quadrature, problem::InitializationTestProblem, x::T, y::T, timestep::Real = 0.0) where { T <: Real }
    return 1.0 - (1/4) * q.speed_of_sound_squared * U^2 * (cos(2x) + cos(2y))
end
LatticeBoltzmann.velocity(problem::InitializationTestProblem, x, y, timestep= 0.0) where { T <: Real } = U * [
    cos(2pi * x) * sin(2pi * y)
    -sin(2pi * x) * cos(2pi * y)
]
function LatticeBoltzmann.velocity_gradient(problem::InitializationTestProblem, x::T, y::T, timestep::Real = 0.0) where { T <: Real }
    u_x = -sin(2pi * x) * sin(2pi * y)
    v_y = sin(2pi * x) * sin(2pi * y)
    u_y = cos(2pi * x) * cos(2pi * y)
    v_x = -cos(2pi * x) * cos(2pi * y)

    return U * 2pi .*  [u_x v_x; u_y v_y]
end

function LatticeBoltzmann.deviatoric_tensor(
    q::Quadrature,
    problem::InitializationTestProblem,
    x::T,
    y::T,
    time::Real = 0.0,
) where { T <: Real }
    a = LatticeBoltzmann.velocity_gradient(problem, x, y, time)
    ν = LatticeBoltzmann.viscosity(problem)

    σ = - ν * [
        2 * a[1, 1] a[1, 2] + a[2, 1]
        a[1, 2] + a[2, 1] 2 * a[2, 2]
    ]

    return σ*(problem.NX)
end

@testset "Initialization strategies" for q in [D2Q9(), D2Q13(), D2Q17(), D2Q21()]
    # q = D2Q9()
    # q = D2Q13()
    # q = D2Q17()
    τ = 1.00 / q.speed_of_sound_squared
    problem = InitializationTestProblem(1 * 8, τ)
# (1 + 1 / (2 * τ))
    τ = 1.0
    # τ = 4.0 / q.speed_of_sound_squared + 0.5

    # @show problem viscosity(problem)
    @show problem.ν * q.speed_of_sound_squared + 0.5
    # @show q.speed_of_sound_squared
    # problem = LatticeBoltzmann.TGV(q, 1.0, 1, 2, 2)
    range_x, range_y = range(problem)

    @testset "Velocity + constant density" begin
        strategy = ConstantDensity()
        # @show strategy
        f = initialize(strategy, q, problem)
        f_ = Array{eltype(q.weights)}(undef, size(f, 3))
        u = zeros(dimension(q))

        for x_idx in OneTo(problem.NX), y_idx in OneTo(problem.NY)
            x = range_x[x_idx]
            y = range_y[y_idx]
            f_ = f[x_idx, y_idx, :]

            ρ = density(q, f_)
            velocity!(q, f_, ρ, u)
            T = temperature(q, f_, ρ, u)
            p = pressure(q, f_, ρ, u)
            σ = deviatoric_tensor(q, τ, f_, ρ, u)
            expected_σ = deviatoric_tensor(q, problem, x, y)

            if (density(q, problem, x, y) ≉ 1.0)
            @test ρ ≉ density(q, problem, x, y)
            @test p ≉ pressure(q, problem, x, y)
            end

            @test u ≈ velocity(problem, x, y)
            @test σ ≉ deviatoric_tensor(q, problem, x, y)
        end
    end

    @testset "Velocity + pressure" begin
        strategy = AnalyticalEquilibrium()
        # @show strategy
        f = initialize(strategy, q, problem)
        f_ = Array{eltype(q.weights)}(undef, size(f, 3))
        u = zeros(dimension(q))

        for x_idx in OneTo(problem.NX), y_idx in OneTo(problem.NY)
            x = range_x[x_idx]
            y = range_y[y_idx]
            f_ = f[x_idx, y_idx, :]

            ρ = density(q, f_)
            velocity!(q, f_, ρ, u)
            T = temperature(q, f_, ρ, u)
            p = pressure(q, f_, ρ, u)
            σ = deviatoric_tensor(q, τ, f_, ρ, u)
            expected_σ = deviatoric_tensor(q, problem, x, y)

            @test ρ ≈ density(q, problem, x, y)
            @test p ≈ pressure(q, problem, x, y)

            @test u ≈ velocity(problem, x, y)
            @test σ ≉ deviatoric_tensor(q, problem, x, y)
        end
    end

    @testset "Velocity + stress initialization" begin
        strategy = AnalyticalVelocityAndStress()
        # @show strategy
        f = initialize(strategy, q, problem)

        f_ = Array{eltype(q.weights)}(undef, size(f, 3))
        u = zeros(dimension(q))

        for x_idx in OneTo(problem.NX), y_idx in OneTo(problem.NY)
            x = range_x[x_idx]
            y = range_y[y_idx]
            f_ = f[x_idx, y_idx, :]

            ρ = density(q, f_)
            velocity!(q, f_, ρ, u)
            T = temperature(q, f_, ρ, u)
            p = pressure(q, f_, ρ, u)
            σ = deviatoric_tensor(q, τ, f_, ρ, u)
            expected_σ = deviatoric_tensor(q, problem, x, y)

            @warn deviatoric_tensor(q, problem, x, y) ./ σ

            if (density(q, problem, x, y) ≉ 1.0)
            @test ρ ≉ density(q, problem, x, y)
            @test p ≉ pressure(q, problem, x, y)
            end

            @test u ≈ velocity(problem, x, y)
            @test σ ≈ deviatoric_tensor(q, problem, x, y)
            break;
        end
    end

    @testset "Velocity + pressure + stress" begin
        strategy = AnalyticalEquilibriumAndOffEquilibrium()
        # @show strategy
        f = initialize(strategy, q, problem)

        f_ = Array{eltype(q.weights)}(undef, size(f, 3))
        u = zeros(dimension(q))

        for x_idx in OneTo(problem.NX), y_idx in OneTo(problem.NY)
            x = range_x[x_idx]
            y = range_y[y_idx]
            f_ = f[x_idx, y_idx, :]

            ρ = density(q, f_)
            velocity!(q, f_, ρ, u)
            T = temperature(q, f_, ρ, u)
            p = pressure(q, f_, ρ, u)
            σ = deviatoric_tensor(q, τ, f_, ρ, u)
            expected_σ = deviatoric_tensor(q, problem, x, y)

            @warn deviatoric_tensor(q, problem, x, y) ./ σ

            @test ρ ≈ density(q, problem, x, y)
            @test p ≈ pressure(q, problem, x, y)

            @test u ≈ velocity(problem, x, y)
            @test σ ≈ deviatoric_tensor(q, problem, x, y)

            break;
        end
    end

    @testset "Mei et al" begin
        strategy = IterativeInitializationMeiEtAl(τ, 1E-14)
        # @show strategy
        f = initialize(strategy, q, problem)

        f_ = Array{eltype(q.weights)}(undef, size(f, 3))
        u = zeros(dimension(q))

        for x_idx in OneTo(problem.NX), y_idx in OneTo(problem.NY)
            x = range_x[x_idx]
            y = range_y[y_idx]
            f_ = f[x_idx, y_idx, :]

            ρ = density(q, f_)
            velocity!(q, f_, ρ, u)
            T = temperature(q, f_, ρ, u)
            p = pressure(q, f_, ρ, u)
            σ = deviatoric_tensor(q, τ, f_, ρ, u)
            expected_σ = deviatoric_tensor(q, problem, x, y)

            # @show deviatoric_tensor(q, problem, x, y) ./ σ

            @test_broken ρ ≈ density(q, problem, x, y)
            @test_broken p ≈ pressure(q, problem, x, y)

            @test_broken u ≈ velocity(problem, x, y)
            @test_broken σ ≈ deviatoric_tensor(q, problem, x, y)
        end
    end
end
