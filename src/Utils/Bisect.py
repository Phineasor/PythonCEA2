# Bisect Function for root finding zeros of a defined function in python, the function needs to be only dependent on x
def Bisect(f, x1, x2, tol, *args):
    Error = 1
    xm = x2
    i = 0
    while (Error > tol) and (i < 50000):
        i += 1
        xmP = xm
        xm = (x1 + x2) / 2
        #print("xm: " + str(f(xm, *args)))
        #print("x1: " + str(f(x1, *args)))
        if f(xm, *args) * f(x1, *args) < 0:
            x2 = xm
            Error = abs(xmP - xm) / xm
        else:
            x1 = xm
            Error = abs(xmP - xm) / xm
    return xm