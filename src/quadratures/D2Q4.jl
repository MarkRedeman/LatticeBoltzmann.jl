# https://github.com/grasingerm/lbxflow/blob/master/inc/lattice.jl
# http://mafija.fmf.uni-lj.si/seminar/files/2012_2013/Lattice_Boltzmann_method.pdf

# const original_order = [1, 2 ,3, 4]

struct D2Q4 <: Quadrature
    abscissae
    weights
    speed_of_sound

    D2Q4() = new(
        [
            1  -1   0   0
            0   0   1  -1
        ],
        [
            1/4, 1/4, 1/4, 1/4
        ],
        1 / √2
    )
end
