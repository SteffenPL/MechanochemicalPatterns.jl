struct PeriodicBoundary{T}
    data::T
end

struct NeumannBoundary{T}
    data::T
end

struct DirichletBoundary{T}
    data::T
end


bcindex(l, ax, bnd::NeumannBoundary) = clamp(l, ax)
bcindex(l, ax, bnd::PeriodicBoundary) = mod(l, ax)
bcindex(l, ax, bnd::DirichletBoundary) = l in ax ? l : 0

bc(u::Float64, inds) = u
function bc(u::AbstractArray, inds)    
    I = map(ind -> bcindex(ind...), inds)
    i_bnd = findfirst(I .== 0)
    return isnothing(i_bnd) ? u[I...] : inds[i_bnd][3].data::Float64
end

function laplace!(du::AbstractArray{Float64,2}, u::AbstractArray{Float64,2}, D, dV; factor = 1.0, boundaries = (NeumannBoundary(nothing), NeumannBoundary(nothing)))
    bnds = boundaries
    
    get_D(i,j) = D * bc(factor, ((i,axes(factor,1), NeumannBoundary(nothing)), (j, axes(factor,2), NeumannBoundary(nothing))))
    get_u(i,j) = bc(u, ((i,axes(u,1),bnds[1]), (j,axes(u,2),bnds[2])))

    @inbounds for i in axes(u,1), j in axes(u,2)
        Dc = get_D(i,j)
        Dl = (0.5 * Dc + 0.5 * get_D(i-1,j))
        Dr = (0.5 * Dc + 0.5 * get_D(i+1,j))
        Db = (0.5 * Dc + 0.5 * get_D(i,j-1))
        Dt = (0.5 * Dc + 0.5 * get_D(i,j+1))

        uc = get_u(i,j)
        ul = get_u(i-1,j)
        ur = get_u(i+1,j)
        ub = get_u(i,j-1)
        ut = get_u(i,j+1)

        du[i,j] += (
            Dl * (ul - uc) / dV[1]^2 + 
            Dr * (ur - uc) / dV[1]^2 +
            Db * (ub - uc) / dV[2]^2 +
            Dt * (ut - uc) / dV[2]^2
            )
    end
    return nothing
end

function laplace!(du::AbstractArray{Float64,3}, u::AbstractArray{Float64,3}, D, dV; factor = 1.0, boundaries = (NeumannBoundary(nothing), NeumannBoundary(nothing), NeumannBoundary(nothing)))
    bnds = boundaries
    
    get_D(i,j,k) = D * bc(factor, ((i,axes(factor,1), NeumannBoundary(nothing)), (j, axes(factor,2), NeumannBoundary(nothing)), (k, axes(factor,3), NeumannBoundary(nothing))))
    get_u(i,j,k) = bc(u, ((i,axes(u,1),bnds[1]), (j,axes(u,2),bnds[2]), (k,axes(u,3),bnds[3])))

    @inbounds for i in axes(u,1), j in axes(u,2), k in axes(u,3)
        Dc = get_D(i,j,k)
        Dl = (0.5 * Dc + 0.5 * get_D(i-1,j,k))
        Dr = (0.5 * Dc + 0.5 * get_D(i+1,j,k))
        Db = (0.5 * Dc + 0.5 * get_D(i,j-1,k))
        Dt = (0.5 * Dc + 0.5 * get_D(i,j+1,k))
        Di = (0.5 * Dc + 0.5 * get_D(i,j,k-1))
        Do = (0.5 * Dc + 0.5 * get_D(i,j,k+1))

        uc = get_u(i,j,k)
        ul = get_u(i-1,j,k)
        ur = get_u(i+1,j,k)
        ub = get_u(i,j-1,k)
        ut = get_u(i,j+1,k)
        ui = get_u(i,j,k-1)
        uo = get_u(i,j,k+1)


        du[i,j,k] += (
            Dl * (ul - uc) / dV[1]^2 + 
            Dr * (ur - uc) / dV[1]^2 +
            Db * (ub - uc) / dV[2]^2 +
            Dt * (ut - uc) / dV[2]^2 + 
            Di * (ui - uc) / dV[3]^2 +
            Do * (uo - uc) / dV[3]^2
            )
    end
    return nothing
end