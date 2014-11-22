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

function d_H_compute(x,xi,N)
    # Compute the derivative of the vector of basis function 
    # for a given location, x, nodal point, xi, and order, N

    d_H = zeros(N+1)

    for i = 2:N+1
        d_H[i] = (i-1)*(x - xi)^(i-2)
    end

    return d_H

end

function phi_a(x,xi,a)
    # Compute the kernel value using the cubic spline 
    # kernel function
    
    z = abs((x-xi)/a) # abs was a typo (excluded) from the 96 paper
    
    if abs(z) <= 0.5
        phi_z = 2./3. - 4.*z^2. + 4.*z^3.
    elseif abs(z) <= 1.
        phi_z = 4./3. - 4.*z + 4.z^2. - (4./3.)*z^3.
    else
        phi_z = 0.
    end
    
    return phi_z
    
end

function d_phi_a(x,xi,a)
    # Compute the kernel derivative using the cubic spline 
    # kernel function
    
    z = abs((x-xi)/a) 
    
    if abs(z) <= 0.5
        d_phi_z = (-8./a*z + 12./a*z^2.)*sign(x-xi)
    elseif abs(z) <= 1.
        d_phi_z = (-4./a + 8./a*z - 4./a*z^2.)*sign(x-xi)
    else
        d_phi_z = 0.
    end
    
    return d_phi_z
    
end

function shape_functions(x,xi,N,a)
    # Compute the shape functions
    
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

function d_shape_functions(x,xi,N,a)
    # Compute the derivative of the shape functions
    
    M = zeros(N+1,N+1,length(x))
    d_M = zeros(N+1,N+1,length(x))
    d_psi = zeros(length(x),length(xi))
    H_0 = H_compute(1.,1.,N)

    for j=1:length(x)
        for i=1:length(xi)
            H = H_compute(x[j],xi[i],N)
            d_H = d_H_compute(x[j],xi[i],N)
            phi = phi_a(x[j],xi[i],a)
            d_phi = d_phi_a(x[j],xi[i],a)
            M[:,:,j] = M[:,:,j] + H * H' * phi
            d_M[:,:,j] = d_M[:,:,j] + d_H * H' * phi +
			              H * d_H' * phi +
			              H * H' * d_phi
        end
        for i=1:length(xi)
            H = H_compute(x[j],xi[i],N)
            d_H = d_H_compute(x[j],xi[i],N)
            phi = phi_a(x[j],xi[i],a)
            d_phi = d_phi_a(x[j],xi[i],a)
	    M_inv = M[:,:,j]^-1
            d_psi[j,i]  = (H_0' * (-M_inv * d_M[:,:,j] * M_inv * H * phi + 
			           M_inv * d_H * phi +
				   M_inv * H * d_phi))[1]
        end
    end
    
    return d_psi

end

function reproduce_u(psi,u)
    # Reproduce u with RK shape functions
    
    u_a = zeros(size(psi,1))
    
    for i=1:size(psi,1)
        u_a[i] = sum(psi[i,:]*u)
    end
    
    return u_a
end

function test()
    # Check phi and derivative of phi
    x=[-10:10]
    xi=0
    phi_z = zeros(length(x))
    num_d_phi_z = zeros(length(x))
    d_phi_z = zeros(length(x))
    for i=1:length(x)
        phi_z[i] = phi_a(x[i],xi,10)
        if i > 1
        	num_d_phi_z[i] = (phi_z[i]-phi_z[i-1])/(x[i]-x[i-1])
        else
        	num_d_phi_z[i] = 0
        end
        d_phi_z[i] = d_phi_a(x[i],xi,10)
    end
    figure()
    plot(x,phi_z)
    plot(x,d_phi_z)
    plot(x,num_d_phi_z,"bx")
    grid()
    savefig(string("pdfs/phi_z_10.pdf"))
    
    x=[-1:0.1:1]
    xi=0
    phi_z = zeros(length(x))
    num_d_phi_z = zeros(length(x))
    d_phi_z = zeros(length(x))
    d_phi_z = zeros(length(x))
    for i=1:length(x)
        phi_z[i] = phi_a(x[i],xi,1)
        if i > 1
        	num_d_phi_z[i] = (phi_z[i]-phi_z[i-1])/(x[i]-x[i-1])
        else
        	num_d_phi_z[i] = 0
        end
        d_phi_z[i] = d_phi_a(x[i],xi,1)
    end
    figure()
    plot(x,phi_z)
    plot(x,d_phi_z)
    plot(x,num_d_phi_z,"bx")
    grid()
    savefig(string("pdfs/phi_z.pdf"))
end

function refinement_study(N,a,u_order,h_list,e_points)
    # h Refinement study of reproducing u with RK 
    # reproduces u(x) with x = 0:1
    #
    # N 	= Order of basis functions
    # a 	a * h = Support size
    # u_order 	= Order of function to reproduce
    # h_list 	= array of different h sizes
    # e_points	= number of evaluation points

    println("*****************************************************")
    println(string("Starting Convergence Study: N=",N,
        	   " a=",a, " uorder=",u_order," h_list=", 
		   h_list," e=",e_points))
    println("*****************************************************")

    u_enorm = zeros(length(h_list))  # allocate for storing error norm

    for k=1:length(h_list)

        tstart = time()

	h = h_list[k]   # this h
        println(string("h=",h))
        eh = 1. / e_points # spacing of eval points

        x = [0:eh:1.]  # Evaluation points
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
        title(string("N=",N," a=",a,
		" uorder=",u_order," h=",h," e=",e_points))
        plot(xi,u,label="in")
        plot(x,u_a,"g--",label="out")
        ylabel("u(x)")
        xlabel("x")
        legend(loc="upper left")
        grid()
        
        subplot(2,1,2)
        plot(x,u_error,label="error")
        ticklabel_format(style="sci", axis="y", scilimits=(0,0))
        ylabel("error")
        xlabel("x")
        grid()
        
        savefig(string("pdfs/reproduced_N=",N,"_a=",a,
		"_uorder=",u_order,"_h=",h,"_e=",e_points,".pdf"))
        tend = time()
        println(string("Time to complete: ",tend-tstart))
    end
    figure()
    loglog(h_list,u_enorm,label=
	    string("N=",N," a=",a,"h uorder=",u_order," e=",e_points))
    ylabel("L2 error")
    xlabel("h size")
    legend(loc="best")
    grid()
    
    savefig(string("pdfs/error_reproduced_N=",N,"_a=",a,
    	"_uorder=",u_order,"_e=",e_points,".pdf"))
    return u_enorm
end

function d_refinement_study(N,a,u_order,h_list,e_points)
    # h Refinement study of reproducing derivative of u with RK 
    # reproduces d_u(x) with x = 0:1
    #
    # N 	= Order of basis functions
    # a 	a * h = Support size
    # u_order 	= Order of function to reproduce
    # h_list 	= array of different h sizes
    # e_points	= number of evaluation points

    println("*****************************************************")
    println(string("Starting Derivative Convergence Study: N=",N,
        	   " a=",a, " uorder=",u_order," h_list=", 
		   h_list," e=",e_points))
    println("*****************************************************")

    d_u_enorm = zeros(length(h_list))  # allocate for storing error norm

    for k=1:length(h_list)

        tstart = time()

	h = h_list[k]   # this h
        println(string("h=",h))
        eh = 1. / e_points # spacing of eval points

        x = [0:eh:1.]  # Evaluation points
        xi = [0.:h:1.]    # Nodal Points
       
        # Create shape functions
        d_psi = d_shape_functions(x,xi,N,a*h)
       
        # Input function (at nodal points)
        u=zeros(length(xi))
        for i=1:length(xi)
            u[i]=(xi[i]/xi[end])^u_order
        end
        
        # True derivative values (at evaluation points)
        d_u_true=zeros(length(x))
        for i=1:length(x)
            d_u_true[i]=u_order*(x[i]/xi[end])^(u_order-1)
        end
       
        # Reproduce function, check error, error norm
        d_u_a = reproduce_u(d_psi,u)
        d_u_error = d_u_a - d_u_true
        d_u_enorm[k] = norm(d_u_error)/length(x)
        
        figure()
        subplot(2,1,1)
        title(string("N=",N," a=",a,
		" uorder=",u_order," h=",h," e=",e_points))
        plot(x,d_u_true,label="in")
        plot(x,d_u_a,"g--",label="out")
        ylabel("du/dx")
        xlabel("x")
        legend(loc="upper left")
        grid()
        
        subplot(2,1,2)
        plot(x,d_u_error,label="error")
        ticklabel_format(style="sci", axis="y", scilimits=(0,0))
        ylabel("error")
        xlabel("x")
        grid()
        
        savefig(string("pdfs/derivative_reproduced_N=",N,"_a=",a,
		"_uorder=",u_order,"_h=",h,"_e=",e_points,".pdf"))
        tend = time()
        println(string("Time to complete: ",tend-tstart))
    end
    figure()
    loglog(h_list,d_u_enorm,label=
	    string("N=",N," a=",a,"h uorder=",u_order," e=",e_points))
    ylabel("L2 error")
    xlabel("h size")
    legend(loc="best")
    grid()
    
    savefig(string("pdfs/derivative_error_reproduced_N=",N,"_a=",a,
    	"_uorder=",u_order,"_e=",e_points,".pdf"))
    return d_u_enorm
end

function run_u_refine()
    h_list = [0.1,0.01,0.001]
    
    e_1 = refinement_study(1,1.1,2,h_list,10000)
    e_2 = refinement_study(1,1.1,3,h_list,10000)
    e_3 = refinement_study(2,2.1,3,h_list,10000)
    
    figure()
    loglog(h_list,e_1,label=string("N=1 a=1.1*h uorder=2, e=1.E4"))
    loglog(h_list,e_2,label=string("N=1 a=1.1*h uorder=3, e=1.E4"))
    loglog(h_list,e_3,label=string("N=2 a=2.1*h uorder=3, e=1.E4"))
    ylabel("L2 error")
    xlabel("h size")
    legend(loc="lower right")
    grid()
    
    savefig(string("pdfs/all_error.pdf"))
end

function run_du_refine()
    h_list = [0.1,0.01,0.001]
    
    e_1 = d_refinement_study(1,1.1,2,h_list,10000)
    e_2 = d_refinement_study(1,1.1,3,h_list,10000)
    e_3 = d_refinement_study(2,2.1,3,h_list,10000)
    
    figure()
    loglog(h_list,e_1,label=string("N=1 a=1.1*h uorder=2, e=1.E4"))
    loglog(h_list,e_2,label=string("N=1 a=1.1*h uorder=3, e=1.E4"))
    loglog(h_list,e_3,label=string("N=2 a=2.1*h uorder=3, e=1.E4"))
    ylabel("L2 error")
    xlabel("h size")
    legend(loc="lower right")
    grid()
    
    savefig(string("pdfs/derivative_all_error.pdf"))
end

run_u_refine()
run_du_refine()
