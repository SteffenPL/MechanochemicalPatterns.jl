
function laplace(du::Array{Float64,3}, u::Array{Float64,3}, D, dV, dt)
    Ddt = D*dt
    @inbounds for i in 2:size(u,1)-1, j in 2:size(u,2)-1, k in 2:size(u,3)-1
        du[i,j,k] += Ddt*(
            (u[i-1,j,k] + u[i+1,j,k] - 2*u[i,j,k])/dV[1] +
            (u[i,j-1,k] + u[i,j+1,k] - 2*u[i,j,k])/dV[2] +
            (u[i,j,k-1] + u[i,j,k+1] - 2*u[i,j,k])/dV[3]
        )
    end
    return nothing
end