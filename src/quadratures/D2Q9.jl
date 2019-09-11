# weights(lattice::D2Q9) = [4/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9]'
# abscissaes(lattice::D2Q9) = [0  0; -1  1; -1  0; -1 -1; 0 -1; 1 -1; 1  0; 1  1; 0  1]'
# opposite(lattice::D2Q9) = [1 2 3 4 5 6 7 8 9]
# opposite(lattice::D2Q9, idx) = idx + (9 - 1) / 2
# dimension(lattice::D2Q9) = 2
# order(lattice::D2Q9) = [1, 7, 9, 3, 5, 8, 2, 4, 6]

# D2Q9
# const original_order = [1, 7, 9, 3, 5, 8, 2, 4, 6]

struct D2Q9 <: Quadrature
    abscissae
    weights
    speed_of_sound_squared

    D2Q9() = new(
        [
            0   -1    -1   -1     0    1     1    1     0
            0    1     0   -1    -1   -1     0    1     1
        ],
        [
            4/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9
        ],
        3.0
    )
end

order(q::D2Q9) = 5

Base.show(io::IO, q::D2Q9)= show(io, "D2Q9")
