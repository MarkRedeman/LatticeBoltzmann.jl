# https://github.com/grasingerm/lbxflow/blob/master/inc/lattice.jl
# http://mafija.fmf.uni-lj.si/seminar/files/2012_2013/Lattice_Boltzmann_method.pdf

# const original_order = [1, 2 ,3, 4]

struct D2Q4 <: Quadrature
    abscissae
    weights
    speed_of_sound_squared

    D2Q4() = new(
        [
            1  -1   0   0
            0   0   1  -1
        ],
        [
            1/4, 1/4, 1/4, 1/4
        ],
        2.0
    )
end

order(q::D2Q4) = 3
Base.show(io::IO, q::D2Q4)= show(io, "D2Q4")
