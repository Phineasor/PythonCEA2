# Bisect Function for root finding zeros of a defined function in python, the function needs to be only dependent on x
def Bisect(f, x1, x2, tol, *args):
    Error = 1
    xm = x2
    i = 0
    fx1 = f(x1, *args)
    while (Error > tol) and (i < 50000):
        i += 1
        xmP = xm
        xm = (x1 + x2) / 2
        fxm = f(xm, *args)
        #print("xm: " + str(f(xm, *args)))
        #print("x1: " + str(f(x1, *args)))
        if fxm * fx1 < 0:
            x2 = xm
            Error = abs(xmP - xm) / xm
        else:
            x1 = xm
            Error = abs(xmP - xm) / xm
            fx1 = f(x1, *args)
    #can probably make this function faster by computing f less.
    return xm

#this function gets the derivative of some other arbitrary function with respect to the first variable int the arbitarry function
#uses the standard trivial center difference method
dx = 10**(-5)
def Derivative(x, f, *args):
    D = (f(x+dx, *args)-f(x-dx, *args))/(2*dx)
    return D