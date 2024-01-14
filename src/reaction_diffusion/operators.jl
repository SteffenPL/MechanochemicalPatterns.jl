
function laplace!(du::AbstractArray{Float64,2}, u::AbstractArray{Float64,2}, D, dV)
    cx(l) = clamp(l, axes(u,1))
    cy(l) = clamp(l, axes(u,2))

    @inbounds for i in axes(u,1), j in axes(u,2)
        du[i,j] += D*(
            (u[cx(i-1),j] + u[cx(i+1),j] - 2*u[i,j])/dV[1]^2 +
            (u[i,cy(j-1)] + u[i,cy(j+1)] - 2*u[i,j])/dV[2]^2 
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