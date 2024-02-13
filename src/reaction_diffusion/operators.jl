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

function laplace!(du::AbstractArray{Float64,3}, u::AbstractArray{Float64,3}, D, dV)
    cx(l) = clamp(l, axes(u,1))
    cy(l) = clamp(l, axes(u,2))
    cz(l) = clamp(l, axes(u,3))

    @inbounds for i in axes(u,1), j in axes(u,2), k in axes(u,3)
        du[i,j,k] += D*(
            (u[cx(i-1),j,k] + u[cx(i+1),j,k] - 2*u[i,j,k])/dV[1]^2 +
            (u[i,cy(j-1),k] + u[i,cy(j+1),k] - 2*u[i,j,k])/dV[2]^2 +
            (u[i,j,cz(k-1)] + u[i,j,cz(k+1)] - 2*u[i,j,k])/dV[3]^2
        )
    end
    return nothing
end

function laplace_periodic!(du::AbstractArray{Float64,2}, u::AbstractArray{Float64,2}, D, dV)
    cx(l) = mod(l, axes(u,1))
    cy(l) = mod(l, axes(u,2))

    @inbounds for i in axes(u,1), j in axes(u,2)
        du[i,j] += D*(
            (u[cx(i-1),j] + u[cx(i+1),j] - 2*u[i,j])/dV[1]^2 +
            (u[i,cy(j-1)] + u[i,cy(j+1)] - 2*u[i,j])/dV[2]^2 
        )
    end
    return nothing
end

function laplace_periodic!(du::AbstractArray{Float64,3}, u::AbstractArray{Float64,3}, D, dV)
    cx(l) = mod(l, axes(u,1))
    cy(l) = mod(l, axes(u,2))
    cz(l) = mod(l, axes(u,3))

    @inbounds for i in axes(u,1), j in axes(u,2), k in axes(u,3)
        du[i,j,k] += D*(
            (u[cx(i-1),j,k] + u[cx(i+1),j,k] - 2*u[i,j,k])/dV[1]^2 +
            (u[i,cy(j-1),k] + u[i,cy(j+1),k] - 2*u[i,j,k])/dV[2]^2 +
            (u[i,j,cz(k-1)] + u[i,j,cz(k+1)] - 2*u[i,j,k])/dV[3]^2
        )
    end
    return nothing
end