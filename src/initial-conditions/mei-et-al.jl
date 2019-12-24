"""
4. Initialization scheme from Mei et al
"""
struct IterativeInitialization <: InitializationStrategy
    τ
end
IterativeInitialization() = IterativeInitialization(0.8)
