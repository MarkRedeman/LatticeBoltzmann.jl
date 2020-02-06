import LatticeBoltzmann: LatticeBoltzmannMethod,
    ProcessingMethod

@testset "TRT with τ_a = τ_s is identical to SRT" begin
    τ = 0.8

    q = D2Q9()
    N = 10
    f_in = zeros(N, N, length(q.weights))
    nx, ny, nf = size(f_in)
    @inbounds for x = 1:nx, y = 1:ny
        @inbounds for f_idx = 1:nf
            f_in[x, y, f_idx] = q.weights[f_idx]
        end
    end
    f_srt = similar(f_in)
    srt = SRT(τ)
    collide!(srt, q, f_in, f_srt)

    f_trt = similar(f_in)
    trt = LatticeBoltzmann.TRT(τ, τ)
    collide!(trt, q, f_in, f_trt)

    @test f_srt ≈ f_trt

    for τ = [0.51, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
        cs² = q.speed_of_sound_squared
        ν = (τ - 0.5) / cs²

        problem = PoiseuilleFlow(ν, 1, static = false)

        process_method = ProcessingMethod(problem, false, 10)
        force_field = LatticeBoltzmann.has_external_force(problem) ? (x_idx, y_idx, t) -> lattice_force(problem, x_idx, y_idx, t) : nothing

        srt_model = LatticeBoltzmannMethod(problem, q, collision_model = SRT(τ), process_method = process_method)
        trt_model = LatticeBoltzmannMethod(problem, q, collision_model = LatticeBoltzmann.TRT(τ, τ), process_method = process_method)

        srt_result = simulate(srt_model, 0:10)
        trt_result = simulate(trt_model, 0:10)

        nx, ny, nf = size(srt_result.f_stream)
        @inbounds for x = 1:nx, y = 1:ny, f_idx = 1:nf
            @test srt_result.f_stream[x, y, f_idx] .≈ trt_result.f_stream[x, y, f_idx]
        end
    end
end

@testset "MRT with τ_i = τ is identical to SRT" begin
end
