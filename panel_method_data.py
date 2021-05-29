import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import path
import scipy.interpolate as interpolate
from numpy import *

def PanelMethod(nodes, alpha, velocity, density, chord):

    #--------------------------------------------------------------
    # Number of panels
    #--------------------------------------------------------------
    n = len(nodes) - 1

    #--------------------------------------------------------------
    # Body points coordinates
    #--------------------------------------------------------------
    xb, yb = [], []

    for i in range(n):
        xb.append(nodes[i][0])
        yb.append(nodes[i][1])

    #--------------------------------------------------------------
    # Control points coordinates
    #--------------------------------------------------------------
    x, y = [], []

    for i in range(n):
        x.append((nodes[i][0] + nodes[i+1][0])/2)
        y.append((nodes[i][1] + nodes[i+1][1])/2)

    #--------------------------------------------------------------
    # Free-stream velocity and Angle of Attack
    #--------------------------------------------------------------
    V_inf = velocity
    rho_inf = density
    c = chord
    a = alpha*pi/180

    #--------------------------------------------------------------
    # Angles of panels (from "x" to the inner side of panel)
    #--------------------------------------------------------------
    def phi_f(i):
        dx = nodes[i+1][0] - nodes[i][0]
        dy = nodes[i+1][1] - nodes[i][1]
        phi = arctan2(dy,dx)
        if phi < 0:
            return phi + 2*pi
        else:
            return phi

    phi = []

    for j in range(n):
        phi.append(phi_f(j))

    #--------------------------------------------------------------
    # Angles between free-stream velocity and outward normal
    #--------------------------------------------------------------
    def beta_f(i):
        b = phi[i] + (pi/2) - a
        if b > 2*pi:
            return b - (2*pi)
        else:
            return b

    beta = []

    for j in range(n):
        beta.append(beta_f(j))

    #--------------------------------------------------------------
    # Panel lengths
    #--------------------------------------------------------------
    def s_f(i):
        dx = nodes[i+1][0] - nodes[i][0]
        dy = nodes[i+1][1] - nodes[i][1]
        return ((dx)**2 + (dy)**2)**0.5

    s = []

    for j in range(n):
        s.append(s_f(j))

    #--------------------------------------------------------------
    # I, J, K and L matrices
    #--------------------------------------------------------------
    def I_K(i,j):
        if i == j:
            return 0, 0
        else:
            A = -(x[i]-nodes[j][0]) * cos(phi[j]) - (y[i]-nodes[j][1]) * sin(phi[j])
            B = (x[i]-nodes[j][0])**2 + (y[i]-nodes[j][1])**2
            C_I = sin(phi[i]-phi[j])
            C_K = -cos(phi[i]-phi[j])
            D_I = -(x[i]-nodes[j][0]) * sin(phi[i]) + (y[i]-nodes[j][1]) * cos(phi[i])
            D_K = (x[i]-nodes[j][0]) * cos(phi[i]) + (y[i]-nodes[j][1]) * sin(phi[i])
            if (B-A**2) <= 0:
                return 0, 0
            else:
                E = (B-A**2)**0.5
                log_term = log(((s[j])**2 + 2*A*s[j] + B)/B)
                arctan_term = arctan2((s[j]+A),E) - arctan2(A,E)
                return (C_I/2) * (log_term) + ((D_I-A*C_I)/E) * (arctan_term), (C_K/2) * (log_term) + ((D_K-A*C_K)/E) * (arctan_term)

    I, K = [[] for _  in range(n)], [[] for _  in range(n)]
    J, L = [[] for _  in range(n)], [[] for _  in range(n)]

    for i in range(n):
        for j in range(n):
            I_value = I_K(i,j)[0]
            K_value = I_K(i,j)[1]
            I[i].append(I_value)
            L[i].append(-I_value)
            K[i].append(K_value)

    J = K

    #--------------------------------------------------------------
    # Generating influence matrix
    #--------------------------------------------------------------
    A = [[] for _ in range(n+1)]

    for i in range(n):                        # i, j -> N-1, N-1
        for j in range(n):
            if i == j:
                A[i].append(pi)
            else:
                A[i].append(I[i][j])

    for i in range(n):                        # Last column (j=N)
        K_sum = 0
        for j in range(n):
            if i == j:
                pass
            else:
                K_sum += K[i][j]
        A[i].append(-K_sum)

    for j in range(n):                        # Last row to j = N-1
        if j == 0:
            A[n].append(J[n-1][0])
        elif j == (n-1):
            A[n].append(J[0][n-1])
        else:
            A[n].append(J[0][j]+J[n-1][j])

    L_sum = 0

    for j in range(1,n-1):                    # Last row and column
        L_sum += L[0][j]+L[n-1][j]

    L_sum += L[n-1][0] + L[0][n-1]

    A[n].append(-(L_sum)+2*pi)


    #--------------------------------------------------------------
    # RHS terms matrix
    #--------------------------------------------------------------
    def b_i(i):
        return -(V_inf) * 2*pi * cos(beta[i])

    b = []

    for i in range(n):
        b.append(b_i(i))

    b.append(-(V_inf)*2*pi*(sin(beta[0])+sin(beta[n-1])))

    #--------------------------------------------------------------
    # Solving Linear system A * x = b for x
    #--------------------------------------------------------------
    invA = linalg.inv(A)   # Inverting influence coefficients matrix
    sol = dot(invA, b)     # Solving using dot product

    #--------------------------------------------------------------
    # Tangential velocity on panel "i" function
    #--------------------------------------------------------------
    def V_t(i):
        term1 = V_inf*sin(beta[i])
        term2 = 0
        for j in range(n):
            term2 += sol[j]*J[i][j]/(2*pi)
        g = sol[n]
        term3 = g/2
        term4 = 0
        for j in range(n):
            term4 += -(g*L[i][j]/(2*pi))

        return term1 + term2 + term3 + term4

    #--------------------------------------------------------------
    # Calculate tangential velocity and Cp on each panel
    #--------------------------------------------------------------
    Vt = []
    Cp = []

    for i in range(n):
        v = V_t(i)
        Vt.append(v)
        Cp.append(1-(v**2/V_inf**2))

    #--------------------------------------------------------------
    # Matrices to calculate Vx and Vy
    #--------------------------------------------------------------
    def M_N(xp,yp,j):
        A = -(xp-nodes[j][0]) * cos(phi[j]) - (yp-nodes[j][1]) * sin(phi[j])
        B = (xp-nodes[j][0])**2 + (yp-nodes[j][1])**2
        C_Mx = -cos(phi[j])
        C_My = -sin(phi[j])
        C_Nx = -C_My
        C_Ny = C_Mx

        D_Mx = (xp-nodes[j][0])
        D_My = (yp-nodes[j][1])
        D_Nx = -D_My
        D_Ny = D_Mx

        if (B-A**2) <= 0:
            return 0, 0, 0, 0
        else:
            E = (B-A**2)**0.5
            log_term = log(((s[j])**2 + 2*A*s[j] + B)/B)
            arctan_term = arctan2((s[j]+A),E) - arctan2(A,E)
            return (C_Mx/2) * (log_term) + ((D_Mx-A*C_Mx)/E) * (arctan_term), (C_My/2) * (log_term) + ((D_My-A*C_My)/E) * (arctan_term), (C_Nx/2) * (log_term) + ((D_Nx-A*C_Nx)/E) * (arctan_term), (C_Ny/2) * (log_term) + ((D_Ny-A*C_Ny)/E) * (arctan_term)

    #--------------------------------------------------------------
    # Velocities functions Vx and Vy
    #--------------------------------------------------------------
    def Vx(xp,yp):
        term1 = V_inf * cos(a)
        term2 = 0
        term3 = 0
        for j in range(n):
            term2 += sol[j]/(2*pi) * M_N(xp,yp,j)[0]
            term3 -= sol[n]/(2*pi) * M_N(xp,yp,j)[2]
        return term1 + term2 + term3

    def Vy(xp,yp):
        term1 = V_inf * sin(a)
        term2 = 0
        term3 = 0
        for j in range(n):
            term2 += sol[j]/(2*pi) * M_N(xp,yp,j)[1]
            term3 -= sol[n]/(2*pi) * M_N(xp,yp,j)[3]
        return term1 + term2 + term3

    #--------------------------------------------------------------
    # Defining grid points
    #--------------------------------------------------------------
    grid_separation = 50

    x_max  = [min(xb)-0.5, max(xb)+0.5]
    y_max  = [min(yb)-0.3, max(yb)+0.3]

    xp_rows  = linspace(x_max[0], x_max[1], grid_separation)
    yp_rows  = linspace(y_max[0], y_max[1], grid_separation)

    xps, yps = meshgrid(xp_rows, yp_rows)

    all_xps = [j for i in xps for j in i]
    all_yps = [j for i in yps for j in i]

    airfoil_body = path.Path(nodes)

    #--------------------------------------------------------------
    # Calculating velocities and Cp at points (xp, yp)
    #--------------------------------------------------------------
    vx, vy = [], []
    Cp_xy = []

    for i in range(len(xps)):
        for j in range(len(xps)):
            if airfoil_body.contains_points([(xps[i][j],yps[i][j])]):
                vx_comp = 0
                vy_comp = 0
                vx.append(0)
                vy.append(0)
            else:
                vx_comp = Vx(xps[i][j],yps[i][j])
                vy_comp = Vy(xps[i][j],yps[i][j])
                vx.append(vx_comp)
                vy.append(vy_comp)
            Cp_xy.append(1-((vx_comp**2+vy_comp**2)**0.5/V_inf)**2)

    U = interpolate.griddata((all_xps, all_yps), vx, (xps, yps), method='nearest')
    V = interpolate.griddata((all_xps, all_yps), vy, (xps, yps), method='nearest')
    Cp_XY = interpolate.griddata((all_xps, all_yps), Cp_xy, (xps, yps), method='nearest')

    #--------------------------------------------------------------
    # Normal and Axial Force Coefficients
    #--------------------------------------------------------------

    CN, CA = [], []

    for j in range(n):
        CN.append(-Cp[j]*s[j]*sin(beta[j]))
        CA.append(-Cp[j]*s[j]*cos(beta[j]))

    #--------------------------------------------------------------
    # Lift and Moment Coefficients
    #--------------------------------------------------------------

    Cl, Cm = 0, 0

    for j in range(n):
        Cl += CN[j]*cos(a) - CA[j]*sin(a)
        Cm += Cp[j]*(x[j]-0.25)*s[j]*cos(phi[j])

    L = 0.5*rho_inf*V_inf**2*c*Cl
    M = 0.5*rho_inf*V_inf**2*c*c*Cl

    return xps, yps, U, V, Cp_XY, x, Cp, Cl, L, Cm, M