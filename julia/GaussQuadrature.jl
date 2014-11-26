using PyPlot
using GaussTables

function gauss_pts(intervals, order)
   gp,gw = GaussTables.pts_and_wts(order)

   pts = zeros(length(intervals-1)*order)
   wts = zeros(length(intervals-1)*order)

   for i=1:length(intervals)-1
      a = intervals[i]
      b = intervals[i+1]
      pts[(i-1)*order+1:i*order] = gp * (b-a)/2. + (b+a)/2.
      wts[(i-1)*order+1:i*order] = gw * (b-a)/2.
   end

   return  pts, wts
end
    
function gauss(vals,wts)
    I = vals'*wts

    return I
end    

function test(u_order,h,g_order)
    # Test Gauss Quadrature 
    #
    # u_order 	= Order of function to integrate
    # h 	= h size
    # g_order	= Order of Gauss Quadrature

    println("*****************************************************")
    println(string("Testing Gauss Quadratrue: ",
        	   " uorder=",u_order," h=", h," gorder=",g_order))

    xi = [0.:h:1.]    # Nodal Points
    
    x, w = gauss_pts(xi, g_order) # Gauss pts and wts


    # Evaluate function
    u = zeros(length(x))
    for i=1:length(x)
        u[i]=x[i]^u_order
    end

    # True Integral
    I_true = xi[end]^(u_order+1) / (u_order+1)

    # Gauss Integral
    I = gauss(u,w)
    
    # Check error, error norm
    I_error = I - I_true

    println(string("True Integral:   ",I_true))
    println(string("Gauss Ingtegral: ",I))
    println(string("Error:           ",I_error))
    println("*****************************************************")

end

test(2,0.1,2)
test(3,0.1,2)
test(4,0.1,2)
test(5,0.1,2)
test(6,0.1,2)
test(2,0.1,3)
test(3,0.1,3)
test(4,0.1,3)
test(5,0.1,3)
test(6,0.1,3)
