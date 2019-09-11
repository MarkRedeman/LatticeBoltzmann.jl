immutable Lattice{D, Quadrature}
    fs::Array{Node, D}
end

abstract Node

immutable InternalNode <: Node
end


# Can do communication
abstract GhostNode <: Node

immutable BoundaryNode <: Node

end
