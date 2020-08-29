##
# goal : quick look-up for positions of a float number in a list
# for example : list <= [1.0, 1.2, 1.2, 1.3, 1.5]
# look up 1.2 returns [2,3], 1.5 returns [5], 0.9 returns []

# potential applications : search the index of a point in a coordinate list of N balls  


##
# helper functions

@inline f1e(n) = parse(Float64,"1e-$n")


@inline uniquen(l,n) = unique(x->round(x,digits=n), l)

# middle position of a list l
@inline mpos(l) = length(l)÷2

# value at middle position
@inline midval(l) = (length(l)==2 ? 0.5*(l[1]+l[2]) : (length(l)==1 ? l[1] : 0.5*(l[mpos(l)]+l[mpos(l)+1])))

# median of a list of real numbers, with precision 10^(-n)
@inline median(lsorted::Vector{T}, n::Int) where {T<:Real} = midval(uniquen(lsorted,n))

# index of median (to the right) 
@inline mid_pos_r(lsorted::Vector{T}, n::Int) where {T<:Real} = findfirst(lsorted.>median(lsorted,n))


##

## call only once for each list `R0_total`

function make_tree(
    R0_total::Vector{T},        # (unsorted) list of numbers
    subset_id::Vector{Int},     # indices in R0_total
    interval::Tuple,            # (min, max)
    ΔR::Real,
    tree::Vector;
    digit = 8
    ) where {T<:Real}

    eps = f1e(digit)

    # note
    # `tree` is a list of tuples (lc, rc, intvL, intvR, ros, ids)
    # a tree has leaves and branching points
    # `lc` and `rc` are the indices of left and right child in `tree`
    # `lc == rc == -1` means leaf, else it is a branching point
    # `intvL` and `intvR` are intervals for the `lchild` and `rchild`
    # `ros` is T[] for a branching point
    # and subset of `R0_total` within the interval for a leaf
    # `ids` is Int[] for a branching point
    # and indices in `R0_total` for the subset within the interval
    
    @assert length(subset_id) > 0
    
    if length(subset_id) == 1 # a leaf
        # interval contains only one element of `R0_total`
        push!(tree, (-1, -1, interval, interval, R0_total[subset_id], subset_id))
        return length(tree) 
        # length after push! is the position of element just pushed into the list 
    elseif (interval[2]-interval[1]) < ΔR # a leaf
        # interval too short
        push!(tree, (-1, -1, interval, interval, R0_total[subset_id], subset_id))
        # length after push! is the position of element just pushed into the list 
        return length(tree)
    else # a branching point
        # break the subset into two parts
        bisect = median(    R0_total[subset_id], digit)
        mp     = mid_pos_r( R0_total[subset_id], digit)
        Lorder = subset_id[1:mp-1]
        Lmax   = min(maximum(R0_total[Lorder])+eps, bisect) # tightening the bounds 
        Rorder = subset_id[mp:end]
        Rmin   = max(minimum(R0_total[Rorder])-eps, bisect) # tightening the bounds
        # construct subtrees recursively
        lchild = make_tree(R0_total, Lorder, (interval[1],Lmax), ΔR, tree)
        rchild = make_tree(R0_total, Rorder, (Rmin,interval[2]), ΔR, tree)
        # third component is [bisection_point_value,] for a branching point
        push!(tree, (lchild, rchild, (interval[1],Lmax), (Rmin,interval[2]), T[], Int[]))
        # length after push! is the position of element just pushed into the list 
        return length(tree)
    end
end

##

# for the entire list of numbers `R0_total`
function make_tree(
    R0_total::Vector{T},
    ΔR::Real;
    digit = 8
    ) where {T<:Real}
    eps = f1e(digit)
    tree = []
    make_tree( R0_total, 
               sortperm(R0_total), 
               (minimum(R0_total)-eps, maximum(R0_total)+eps), 
               ΔR, 
               tree,
               digit=digit )
    return tree
end

##

# split the tree into three lists to speed up the look-up operation 
mutable struct LookupTable{T}
    Children::Vector{Tuple{Int,Int}}
    Intervals::Vector{Tuple{Tuple{T,T},Tuple{T,T}}}
    Data::Vector{Vector{T}}
    IDs::Vector{Vector{Int}}
end

tree2lkpt(tree::Vector) = LookupTable{eltype(tree[1][3])}(
                                map(x->(x[1],x[2]),tree), 
                                map(x->(x[3],x[4]),tree), 
                                map(x->x[5],tree), 
                                map(x->x[6],tree) )

##

function look_up_value_in_tree(val, tree::Vector; eps=1e-16)
    pos = length(tree)
    while true
        (lc,rc) = tree[pos][1:2]
        if (lc==-1 && rc==-1)
            break
        else
            intvL = tree[pos][3]
            intvR = tree[pos][4]
            if (val<intvL[1] || val>intvR[2] || (val<intvR[1] && val>intvL[2]))
                return Int[]
            else
                pos = val<intvL[2] ? lc : rc
            end
        end
    end
    k = findall(abs.(tree[pos][5].-val).<eps)
    return (k==nothing) ? Int[] : tree[pos][6][k]
end


function look_up_value_in_tree(val, lkpt::LookupTable; eps=1e-16)
    pos = length(lkpt.Children)
    while true
        (lc,rc) = lkpt.Children[pos]
        if lc==-1 && rc==-1
            break
        else
            intvL, intvR = lkpt.Intervals[pos]
            if (val<intvL[1] || val>intvR[2] || (val<intvR[1] && val>intvL[2]))
                return Int[]
            else
                pos = val<intvL[2] ? lc : rc
            end
        end
    end
    k = findall(abs.(lkpt.Data[pos].-val).<eps)
    return (k==nothing) ? Int[] : lkpt.IDs[pos][k]
end


##
# test

R0 = ones(Float64,100) ;
tree = make_tree(R0,0.01) ;
LT = tree2lkpt(tree) ;
look_up_value_in_tree(1.0,LT)
##
look_up_value_in_tree(2.0,LT)

##
# test

R0 = [ones(Float64,100); 2ones(Float64,100)] ;
tree = make_tree(R0,0.01) ;
LT = tree2lkpt(tree) ;
look_up_value_in_tree(2.0,LT)
##
look_up_value_in_tree(1.0,LT)

##


# test

R0 = rand(400000) ;
tree = make_tree(R0,0.01) ;
LT = tree2lkpt(tree) ;

#
##

@time all([i==look_up_value_in_tree(R0[i],tree)[1] for i=1:length(R0)])

##

# 5x speed up for the version with `LookupTable` 
@time all([i==look_up_value_in_tree(R0[i],LT)[1] for i=1:length(R0)])



## -------------------------------------------------------------------------------------
 
# application

function find_position_point(
    pt::Vector{V},
    LTS;
    eps=1e-16
    ) where {V<:Real}
    @assert length(pt) == length(LTS)
    intersect([look_up_value_in_tree(pt[k],LTS[k],eps=eps) for k=1:length(pt)]...)
end


EulerMatrix(α::Real,β::Real,γ::Real) = hcat(
                [cos(α)*cos(β)*cos(γ) - sin(α)*sin(γ),  cos(β)*cos(γ)*sin(α) + cos(α)*sin(γ),  -(cos(γ)*sin(β))],
                [-(cos(γ)*sin(α)) - cos(α)*cos(β)*sin(γ),  cos(α)*cos(γ) - cos(β)*sin(α)*sin(γ),  sin(β)*sin(γ)],
                [cos(α)*sin(β),  sin(α)*sin(β),  cos(β)] )


## test with a 3d example ------

R = EulerMatrix(2π*rand(),π*rand(),2π*rand()) ; # random rotation to avoid aliasing

grid = R*hcat([[a,b,c] for a=0:0.2:10.0 for b=0:0.1:5.0 for c=0:0.15:7.5]...) ;

lookuptables = [make_tree(vec(grid[1,:]),0.01), 
                make_tree(vec(grid[2,:]),0.01),
                make_tree(vec(grid[3,:]),0.01)] .|> tree2lkpt ;



## test speed and verify results
@time all([i==find_position_point(grid[:,i],lookuptables)[1] for i=1:size(grid,2)])


## test with a 2d example ------

R = EulerMatrix(2π*rand(),0,0)[1:2,1:2] ; # random rotation to avoid aliasing

grid = R*hcat([[a,b] for a=0:0.2:200.0 for b=0:0.12:120.0]...) ;

lookuptables = [make_tree(vec(grid[1,:]),0.01), 
                make_tree(vec(grid[2,:]),0.01)] .|> tree2lkpt ;



## test speed and verify results
@time all([i==find_position_point(grid[:,i],lookuptables)[1] for i=1:size(grid,2)])

##
size(grid)