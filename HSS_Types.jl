## HSS_Types.jl

module HSS_Types

export Leaf, Node, Hss, Leaf_Spine, Node_Spine, Tree

mutable struct Leaf{T1 <: Integer, T2 <: Number}
    m::T1
    n::T1
    depth::T1
    d::T1 #diagonal index
    U::Array{T2,2}
    V::Array{T2,2}
    D::Array{T2,2}
end

#Define node
mutable struct Node{T1 <: Integer, T2 <: Number}
    treeL::Union{Node,Leaf} #Union(Node{T1,T2},Leaf{T1,T2}) doesn't work for some reason
    treeR::Union{Node,Leaf}
    m::T1
    n::T1
    depth::T1
    d::T1
    Bu::Array{T2,2}
    Bl::Array{T2,2}
    Rl::Array{T2,2}
    Rr::Array{T2,2}
    Wl::Array{T2,2}
    Wr::Array{T2,2}
end

#Define HSS structure
const Hss = Union{Leaf,Node}

# #Define partition(spine) tree
mutable struct Leaf_Spine{T <: Integer}
    m::T
    n::T
    depth::T
    d::T
end

mutable struct Node_Spine{T <: Integer}
    treeL::Union{Node_Spine{T},Leaf_Spine{T}}
    treeR::Union{Node_Spine{T},Leaf_Spine{T}}
    m::T
    n::T
    depth::T
    d::T
end

Tree = Union{Leaf_Spine,Node_Spine}

end
