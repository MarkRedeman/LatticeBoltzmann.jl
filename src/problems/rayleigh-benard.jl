export RayleighBenard

struct RayleighBenard <: lbm.InitialValueProblem
    rho_0::Float64
    u_max::Float64
    ν::Float64
    κ::Float64
    NX::Int64
    NY::Int64
    k::Float64
    domain_size::Tuple{Float64, Float64}
end
