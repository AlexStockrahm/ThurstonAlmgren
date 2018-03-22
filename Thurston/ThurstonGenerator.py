# -*- coding: utf-8 -*-
"""
Created on Wed Jul 06 10:58:31 2016

@author: alex.stockrahm
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math as m
import numpy as np



#assume ellipse described by x^2+y^2+(z/4)^2 =< 1 


#for z on [-c, c], c=4, return radius in x,y plane for given z, 
#z is a scalar
def gen_r_z(z):
    return (1-(z/4.0)**2)**.5

#Generate the increment dz over a spiral
#take from -z_o to z_o, divide by how many substeps you would like 
#you will have to repeat it for each winding, as this dz is for one winding
def gen_dz(winding, substeps, z_o):
    return 2.0*z_o/winding/substeps

#substeps is taken as per winding, same as for gen_dz    
def gen_dtheta(substeps):
    return (m.pi*2.0/substeps)

#same rules as above, substeps per winding
#z_0 gives unit with 
def pickPoints(ls, substeps):
    bezierControlIndex = int(.08*substeps)
    lastControl= ls[0]
    last = ls[bezierControlIndex]
    first = ls[-bezierControlIndex]
    firstControl = [-1]
    
    return ls[bezierControlIndex:-bezierControlIndex], lastControl, last, first, firstControl



def CubicBez(x0, x1, x2, x3, t):
    #cubic bezier from explicit form, t on [0,1]
    return list((1-t)**3*x0+3*(1-t)**2*t*x1+3*(1-t)*t**2*x2+t**3*x3)
    #you can your np.array.tolist() if this doesn't work

def bezier(first, last, FirstControl, SecondControl, steps):
    #give 4 points as a list, quantity of steps in between
    #return a cubic bezier curve xs, ys, zs between first, last
    t = np.linspace(0, 1, num=steps)
    x = CubicBez(first[0], FirstControl[0], SecondControl[0], last[0], t)
    y = CubicBez(first[1], FirstControl[1], SecondControl[1], last[1], t)
    z = CubicBez(first[2], FirstControl[2], SecondControl[2], last[2], t)
    return x, y, z
        
    
  
def generate_spiral(z_0, winding_no, substeps):
    dz = gen_dz(winding_no, substeps, z_0)
    dtheta = gen_dtheta(substeps)
    thetaCycle = [k*dtheta for k in range(substeps)]

    xList = []
    yList = []
    zList = []
    #pick out first, first_control etc.

    
    #winding no, subs are intengers, add 1 to be inclusinve on [-z_0, z_0]
    #k is an element from the index of each subdivision
    for k in range(winding_no*substeps+1):
        thetaIndex = k % substeps
        theta = thetaCycle[thetaIndex]
        #calculate z
        z = -1*z_0 + k* dz
        #calculate r assuming x and y sweep out a circle
        r_z = gen_r_z(z)
        x = r_z * m.cos(theta)
        y = r_z * m.sin(theta)
        
        xList.append(x)
        yList.append(y)
        zList.append(z)
        
    xList, xlastControl, xlast, xfirst, xfirstControl = pickPoints(xList, substeps)
    yList, ylastControl, ylast, yfirst, yfirstControl = pickPoints(yList, substeps)
    zList, zlastControl, zlast, zfirst, zfirstControl = pickPoints(zList, substeps)
    
    return xList, yList, zList     

def generate_ellipse(z_o, major_ax, wire_thickness, substeps):
    #Generate the bounding elipse
    larger_ax = 2*wire_thickness+major_ax
    l_minor_ax = 2*wire_thickness+1
    
    """
    This is  such that along the major and minor
    axes, the wire barely touches the bounding ellipse, not clear if this
    definitively keeps the spiral from touching the ellipse
    
    """
    jump = 5
    dz = 2.0*(z_o-wire_thickness) / substeps
    
    #Part (i), from z_o to -z_o, Thurston Almgren p back towards q along -z direction
    z_ret = [z_o - k*dz for k in range(substeps) ]
    x_ret = [l_minor_ax*(1-(z/larger_ax)**2)**.5 for z in z_ret]
    y_ret = [0 for k in range(substeps)]
    
#save elements for return journey
    z_ret3 = list( reversed(z_ret))
    z_ret3 = z_ret3[:-3]
    x_ret3 = [-1*x for x in reversed(x_ret)]
    x_ret3 = x_ret3[:-3]
    y_ret3 = [0 for k in range(len(x_ret3))]

    
#####################################################################    
#    #Part(ii)
    dz = (larger_ax + 2*wire_thickness-z_o )/substeps
    tight_ell_z = [-1*z_o - k*dz + 2*wire_thickness for k in range(substeps)]
    #curve starting to bend to cross z axis
    #1st half of this bend for z

    tight_x = [l_minor_ax*(1-(z/larger_ax)**2)**.5 for z in tight_ell_z]
    #same 1st half bend for x
    #append the same list with x switched to neg in reversed order
    
        #append z values in reverse order 
    tight_ell_z.extend(reversed(tight_ell_z)) 
    # second half of bend for z
    tight_x.extend([x*-1 for x in reversed(tight_x)])
    #append the first swath switched to negative in reversed order
    #second half of bend for x
    
    

    #stay in the x-z plane, 
    
    
    
    tight_y = [ 0 for k in range(4*(substeps-1))]
    tight_y.extend(list(np.linspace(0,2*wire_thickness, num=substeps-1, endpoint= False )))
    
    
    tight_y.extend([ 2*wire_thickness for k in range(substeps)])
    
    
    #print(tight_y)


    

    x_ret.extend(tight_x) 
    z_ret.extend(tight_ell_z)    
    #first loop to reverse direction
    
    
    x_ret.extend(x_ret3) # return journey
    z_ret.extend(z_ret3) # return journey
    
   
    #second loop where you need to rise
    x_ret.extend(reversed(tight_x))
    x_ret.extend(x_ret[jump:substeps-1])
    
    y_ret.extend(tight_y[:-4])
    
    
    
    z_ret.extend([z*-1 for z in tight_ell_z])
    z_ret.extend(z_ret[jump:substeps-1]) 
    
    
    
    

    return [x_ret, y_ret, z_ret]
    
    ##############################
    
x, y, z = generate_spiral(3.3, 7, 200)
half_ellipse_x, half_ellipse_y, half_ellipse_z = generate_ellipse(3.3, 4, .1, 100)

print(len(half_ellipse_x), len(half_ellipse_y), len(half_ellipse_z))

x.extend(half_ellipse_x)
y.extend(half_ellipse_y)
z.extend(half_ellipse_z)


##Test script for the bezier
#first=[0, 0, 0]
#last = [10, 10, 10]
#firstControl = [ -2, 3, 3]
#secondControl = [8, 11, 7]

#x, y, z = bezier(first, last, firstControl, secondControl, 200)





fig = plt.figure()

ax = fig.add_subplot(111, projection = '3d')
ax.plot(x, y, z)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
