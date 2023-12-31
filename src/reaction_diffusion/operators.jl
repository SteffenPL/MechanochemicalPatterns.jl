
function laplace!(du::Array{Float64,2}, u::Array{Float64,2}, D, dV, dt)
    Ddt = D*dt
    cx(l) = clamp(l, axes(u,1))
    cy(l) = clamp(l, axes(u,2))

    @inbounds for i in axes(u,1), j in axes(u,2)
        du[i,j] += Ddt*(
            (u[cx(i-1),j] + u[cx(i+1),j] - 2*u[i,j])/dV[1] +
            (u[i,cy(j-1)] + u[i,cy(j+1)] - 2*u[i,j])/dV[2] 
        )
    end
    return nothing
end

function laplace!(du::Array{Float64,3}, u::Array{Float64,3}, D, dV, dt)
    Ddt = D*dt
    cx(l) = clamp(l, axes(u,1))
    cy(l) = clamp(l, axes(u,2))
    cz(l) = clamp(l, axes(u,3))

    @inbounds for i in axes(u,1), j in axes(u,2), k in axes(u,3)
        du[i,j,k] += Ddt*(
            (u[cx(i-1),j,k] + u[cx(i+1),j,k] - 2*u[i,j,k])/dV[1] +
            (u[i,cy(j-1),k] + u[i,cy(j+1),k] - 2*u[i,j,k])/dV[2] +
            (u[i,j,cz(k-1)] + u[i,j,cz(k+1)] - 2*u[i,j,k])/dV[3]
        )
    end
    return nothing
end