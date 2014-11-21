using PyPlot

function H_compute(x,xi,N)
    # Compute the vector of basis function for a given
    # location, x, nodal point, xi, and order, N

    H = zeros(N+1)

    for i = 1:N+1
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
    # reproduces u(x) with x = 0:1
    #
    # N 	= Order of basis functions
    # a 	a * h = Support size
    # u_order 	= Order of function to reproduce
    # h_list 	= array of different h sizes
    # e_per_h	= number of evaluation points in h interval

    println("*****************************************************")
    println(string("Starting Convergence Study: N=",N,
        	   " a=",a, " uorder=",u_order," h_list=", 
		   h_list," e=",e_per_h))
    println("*****************************************************")

    u_enorm = zeros(length(h_list))  # allocate for storing error norm

    for k=1:length(h_list)

        tstart = time()

	h = h_list[k]   # this h
        println(string("h=",h))
        eh = h / e_per_h # spacing of eval points

        x = [eh/2:eh:1.]  # Evaluation points
        xi = [0.:h:1.]    # Nodal Points
       
        # Create shape functions
        psi = shape_functions(x,xi,N,a*h)
       
        # Input function (at nodal points)
        u=zeros(length(xi))
        for i=1:length(xi)
            u[i]=(xi[i]/xi[end])^u_order
        end
        
        # True values (at evaluation points)
        u_true=zeros(length(x))
        for i=1:length(x)
            u_true[i]=(x[i]/xi[end])^u_order
        end
       
        # Reproduce function, check error, error norm
        u_a = reproduce_u(psi,u)
        u_error = u_a - u_true
        u_enorm[k] = norm(u_error)/length(x)
        
        figure()
        subplot(2,1,1)
        plot(xi,u,label="in")
        plot(x,u_a,"gx",label="out")
        ylabel("u(x)")
        legend(loc="upper left")
        
        subplot(2,1,2)
        plot(x,u_error,label="error")
        ylabel("error")
        
        savefig(string("pdfs/reproduced_N=",N,"_a=",a,
		"_uorder=",u_order,"_h=",h,"_e=",e_per_h,".pdf"))
        tend = time()
        println(string("Time to complete: ",tend-tstart))
    end
    figure()
    loglog(h_list,u_enorm,label="L2 error")
    ylabel("L2 error")
    
    savefig(string("pdfs/error_reproduced_N=",N,"_a=",a,
    	"_uorder=",u_order,"_e=",e_per_h,".pdf"))
    return u_enorm
end

h_list = [0.1,0.01,0.001]

e_1 = refinement_study(1,1.1,2,h_list,10)
e_2 = refinement_study(1,1.1,3,h_list,10)
e_3 = refinement_study(2,2.1,3,h_list,10)

figure()
loglog(h_list,e_1,label=string("N=1 a=1.1*h uorder=2, e=10*h"))
loglog(h_list,e_3,label=string("N=1 a=1.1*h uorder=3, e=10*h"))
loglog(h_list,e_2,label=string("N=2 a=2.1*h uorder=3, e=10*h"))
ylabel("L2 error")
xlabel("h size")
legend(loc="upper left")
grid()

savefig(string("pdfs/all_error.pdf"))
