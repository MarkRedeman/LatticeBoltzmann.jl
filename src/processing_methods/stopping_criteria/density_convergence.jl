struct DensityConvergence{FT <: Real, MT <: AbstractMatrix{FT}} <: StopCriteria
    ϵ::FT
    ρ_old::MT
    ρ::MT
end
function should_stop!(stop_criteria::DensityConvergence, q, f_in)
    nx, ny = size(stop_criteria.ρ_old)
    for x_idx in nx, y_idx in ny
        stop_criteria.ρ_old[x_idx, y_idx] = stop_criteria.ρ[x_idx, y_idx]
        stop_criteria.ρ[x_idx, y_idx] = density(q, f_in[x_idx, y_idx, :])
    end

    δρ = norm(stop_criteria.ρ - stop_criteria.ρ_old)

    # Stop procesing if the density converged or the method seems to be diverging
    return δρ < stop_criteria.ϵ || δρ > 100.0
end
