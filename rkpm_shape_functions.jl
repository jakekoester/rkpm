function H_compute(x,xi,N)
    # Compute the vector of basis function for a given
    # location, x, nodal point, xi, and order, N

    H = zeros(N+1)
    H[1] = 1.

    for i = 2:N+1
        H[i] = (x - xi)^(i-1)
    end

    return H

end

function phi_a(x,xi,a)
    # Compute the kernel value using the cubic spline 
    # kernel function
    
    z = abs((x-xi)/a) # abs was a typo (excluded) from the '96 paper
    
    if abs(z) <= 0.5
        phi_z = 2./3. - 4.*z^2. + 4.*z^3.
    elseif abs(z) <= 1.
        phi_z = 4./3. - 4.*z + 4.z^2. - (4./3.)*z^3.
    else
        phi_z = 0.
    end
    
    return phi_z
    
end

function shape_functions(x,xi,N,a)
    
    M = zeros(N+1,N+1,length(x))
    psi = zeros(length(x),length(xi))
    H_0 = H_compute(1.,1.,N)

    for j=1:length(x)
        for i=1:length(xi)
            H = H_compute(x[j],xi[i],N)
            phi = phi_a(x[j],xi[i],a)
            M[:,:,j] = M[:,:,j] + H * H' * phi
        end
        for i=1:length(xi)
            H = H_compute(x[j],xi[i],N)
            phi = phi_a(x[j],xi[i],a)
            psi[j,i]  = (H_0' * M[:,:,j]^-1 * H * phi)[1]
        end
    end
    
    return psi

end

function reproduce_u(psi,u)
    # Reproduce u with RK shape functions
    
    u_a = zeros(size(psi,1))
    
    for i=1:size(psi,1)
        u_a[i] = sum(psi[i,:]*u)
    end
    
    return u_a
end


function refinement_study(N,a,u_order,h_list,e_per_h)
    # h Refinement study of reproducing u with RK 
    #
    # N 	= Order of basis functions
    # a 	= Support size
    # u_order 	= Order of function to reproduce
    #             {zero, linear, quadratic, cubic}
    # h_list 	= array of different h sizes
    # e_per_h	= number of evaluation points in h interval

    x = [0.01:0.01:1.]  # Evaluation points
    xi = [0.:0.1:1.]      # Nodal Points
    
    psi = shape_functions(x,xi,N,a)
    
    # Input function (at nodal points)
    u=zeros(length(xi))
    for i=1:length(xi)
        u[i]=(xi[i]/xi[end])^2.
    end
    
    # True values (at evaluation points)
    u_true=zeros(length(x))
    for i=1:length(x)
        u_true[i]=(x[i]/x[end])^2.
    end
    
    u_a = reproduce_u(psi,u)
    u_error = u_a - u_true
    
    using PyPlot
    figure()
    subplot(2,1,1)
    plot(xi,u,label="in")
    plot(x,u_a,"gx",label="out")
    ylabel("u(x)")
    legend(loc="upper left")
    
    subplot(2,1,2)
    plot(x,u_error,label="error")
    ylabel("error")
    
    savefig("reproduced.pdf")
end


N = 2   # Order of basis functions
a = 3.1 # Support size
