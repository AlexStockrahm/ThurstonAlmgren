 
from sys import argv
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import cross,eye,dot
from scipy.linalg import expm3,norm
from math import sin,cos

sys.path.append('..')

from transformations import transformations as tra
from simplex_class import simplex
from Thurston import ThurstonGenerator as tg

#Set up parameters of knot: How many points per cross section (points_on_tube slice), how many points to include on the curve parameterization (steps)
pi=math.pi
steps = 30
points_on_tube_slice = 3
rotation_angle = 2*pi/points_on_tube_slice

#p_param , q_param = int(argv[1]), int(argv[2])
phi =np.linspace(0, 2*pi, steps)
phi = phi[:-1]
scaling_factor = .1
#print(phi)
#print(list(range(1,points_on_tube_slice)))


#Parameterization of spiral over elipse:
 
x, y, z = tg.generate_spiral(3.3, 7, 50)
half_ellipse_x, half_ellipse_y, half_ellipse_z = tg.generate_ellipse(3.3, 4, .1, 30)

x.extend(half_ellipse_x)
y.extend(half_ellipse_y)
z.extend(half_ellipse_z)

#update this later to allow for arbitrary torus knot parameterizations


#Obtain Unit vectors for the tangent of points on the knot
b=np.column_stack([x,y,z])
tangents = np.diff(b , axis=0)
#add on difference between the two end points of the knot
tangents= np.concatenate(( tangents, [b[0,:] - b[-1,:]] ), axis=0)



lengths = np.apply_along_axis(norm , 1 , tangents)
lengths = np.column_stack([lengths, lengths, lengths])


unit_rotation_axes=tangents/lengths
#unit rotation axis serve as the tangents at each point along curve


cross_vectors = np.zeros_like(tangents)
cross_vectors[:, 2]= 1
angle = math.pi/90
rotatingMat = tra.rotation_matrix(angle, [1, 0, 0])[:3,:3]


current_cross = cross_vectors[0]

for i in range(len(cross_vectors)):
    inLineMeasurement = np.abs(np.dot(unit_rotation_axes[i], current_cross))
#    print(i)
#    print( inLineMeasurement)
#    print("\n")
    while inLineMeasurement > .999:
        print("{step} cross yields no value, rotating".format(step=i))
        current_cross = np.dot(rotatingMat, current_cross)
        print(current_cross)
        inLineMeasurement = np.abs(np.dot(unit_rotation_axes[i], current_cross))
    cross_vectors[i] = current_cross
        
        



# cross each tangent w/ 0,0,1, leaves normal in x,y plane
normals_to_tangents=cross( cross_vectors , unit_rotation_axes )



lengths = np.apply_along_axis(norm , 1 , normals_to_tangents)
lengths = np.column_stack([lengths, lengths, lengths])
unit_normal_axes=normals_to_tangents/lengths


#construct tubular neighborhood of knot (Scaling factor gives radius of the tubular neighborhood)


normals_to_tangents= scaling_factor*unit_normal_axes
shifted_knot = b + normals_to_tangents

#shifted knot created by bumping the points on the tangent line along normal


####################################################################

############/***Rotation about unit tangent vectors*****/###########

####################################################################



#Convert the shifted knot to a 4-D vector for use with quaternions
points_on_shifted_knot = len(shifted_knot)
shifted_knot = np.concatenate((shifted_knot, np.ones((points_on_shifted_knot,1))), axis=1)


#Now ready to start building the simplicial complex:
tube=simplex.simplex()

#This loop adds the first "row" of points and specifies the edge direction
#as upwards
for i in list(range(points_on_shifted_knot)):

########!!!!*Heres the important part*!!!!!###########
   point=b[i,:] #the point you'll rotate about
   rot_mat = tra.rotation_matrix(rotation_angle, unit_rotation_axes[i,:], point )
   s_knot_point = shifted_knot[i,:]
   tube.addpoint(s_knot_point.tolist(), rot_mat , 'first')



points_thus_far = tube.get_points_with_state('first')
#Create the first row of edges and generate first row of points
for p in list(range(len(points_thus_far))):
    tube.addedge(points_thus_far[p-1] , points_thus_far[p], 'first')
    new_point = dot(points_thus_far[p].rotation_matrix,np.array(points_thus_far[p].coordinates))
    tube.addpoint(new_point.tolist(), points_thus_far[p].rotation_matrix, 'generated')
    
#generate the edges for the faces to be constructed in this row
#set the generated points to the current row as you go (for next round)
generated_points = tube.get_points_with_state('generated')
for i in list(range(len(generated_points))):
    tube.addedge(generated_points[i-1],generated_points[i],'down')
    #make sure to set the 'down' state with negative id's 
    #since they're actually pointing up
    #change these to 'up' as you make faces
    tube.addedge(points_thus_far[i] , generated_points[i-1], 'dd')
    tube.addedge(generated_points[i] , points_thus_far[i] , 'rl')
    generated_points[i].state_point = 'current'

#prepare the faces
first_edges = tube.get_edges_with_state('first')
down_edges = tube.get_edges_with_state('down')
right_left_edges = tube.get_edges_with_state('rl')
diagonals = tube.get_edges_with_state('dd')
for i in list(range(len(first_edges))):
    #lower triangle
    tube.addface(first_edges[i].id_number , diagonals[i].id_number , right_left_edges[i-1].id_number)
    #upper triangle
    tube.addface(-(down_edges[i].id_number), -(diagonals[i].id_number) , -(right_left_edges[i].id_number))
    down_edges[i].state_edge = 'up'
    right_left_edges[i].state_edge = 'spent'
    diagonals[i].state_edge = 'spent'

####################################################
###///***Loop for filling in the middle***////######
####################################################
for i in list(range(1, points_on_tube_slice-1)):
    current_row=tube.get_points_with_state('current')
    for p in list(range(len(current_row))):
        new_point = dot(current_row[p].rotation_matrix,np.array(current_row[p].coordinates))
        tube.addpoint(new_point.tolist(), points_thus_far[p].rotation_matrix, 'generated')
    generated_points = tube.get_points_with_state('generated')
    for i in list(range(len(generated_points))):
        tube.addedge(generated_points[i-1],generated_points[i],'down')
        #make sure to set the 'down' state with negative id's 
        #since they're actually pointing up
        #change these to 'up' as you make faces
        tube.addedge(current_row[i] , generated_points[i-1], 'dd')
        tube.addedge(generated_points[i] , current_row[i] , 'rl')
        generated_points[i].state_point = 'current'
        current_row[i].state_point = 'spent'

        #prepare the faces
    up_edges = tube.get_edges_with_state('up')
    down_edges = tube.get_edges_with_state('down')
    right_left_edges = tube.get_edges_with_state('rl')
    diagonals = tube.get_edges_with_state('dd')
    for i in list(range(len(first_edges))):
        #lower triangle
        tube.addface(up_edges[i].id_number , diagonals[i].id_number , right_left_edges[i-1].id_number)
        #upper triangle
        tube.addface(-(down_edges[i].id_number), -(diagonals[i].id_number) , -(right_left_edges[i].id_number))
        down_edges[i].state_edge = 'up'
        right_left_edges[i].state_edge = 'spent'
        diagonals[i].state_edge = 'spent'
        up_edges[i].state_edge = 'spent'

#####################################################
#####Sew Up the last bit#############################
#####################################################

current_row=tube.get_points_with_state('current')
generated_points = tube.get_points_with_state('first')
for i in list(range(len(generated_points))):

    #youre using the points from the original row, so you only need to 
    #make diagonals and crosses
    tube.addedge(current_row[i] , generated_points[i-1], 'dd')
    tube.addedge(generated_points[i] , current_row[i] , 'rl')
    generated_points[i].state_point = 'sewn'
    current_row[i].state_point = 'spent'

    #prepare the faces
up_edges = tube.get_edges_with_state('up')
down_edges = tube.get_edges_with_state('first')
right_left_edges = tube.get_edges_with_state('rl')
diagonals = tube.get_edges_with_state('dd')
for i in list(range(len(first_edges))):
    #lower triangle
    tube.addface(up_edges[i].id_number , diagonals[i].id_number , right_left_edges[i-1].id_number)
    #upper triangle
    tube.addface(-(down_edges[i].id_number), -(diagonals[i].id_number) , -(right_left_edges[i].id_number))
    down_edges[i].state_edge = 'up'
    right_left_edges[i].state_edge = 'spent'
    diagonals[i].state_edge = 'spent'
    up_edges[i].state_edge = 'spent'


#####################################################
## Plotting ###
#####################################################
"""
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111 , projection='3d')

for i in list(range(len(tube.faces))):
    e_ids =tube.faces[i].edge_ids
    p_ids = [tube.edges[abs(j)-1].point_ids for j in e_ids]
    coordinates = [tube.points[j[i]].coordinates for j in p_ids for i in range(2)]
    #ax.scatter(coordinates[:][0],coordinates[:][1], coordinates[:][2])
    xcoords = [coordinates[i][0] for i in range(len(coordinates))]
    ycoords = [coordinates[i][1] for i in range(len(coordinates))]
    zcoords = [coordinates[i][2] for i in range(len(coordinates))]
    #print('What youre printing:' , xcoords, len(coordinates))
    ax.plot(xcoords, ycoords , zcoords )

plt.show()
"""

###### Write to Geometry File (gmsh) ###########
characteristic_length = .1
#!!!!make this an optional parameter!!!!####
f = open('ThurAlm.geo' , 'w')

## Use this section to write the parameterization of the knot
## Should include number of steps, number of rotations,
## P, Q and (When you set it up) bounding sphere radius
description_string = "//Gmsh geometry file for Thurston Almgren curve"
f.write(description_string)
f.write("\n\n\n")
f.write("//List of Points:\n")

for i in tube.points:
    point_string = "Point ( {pid} ) = {{ {x} , {y} , {z} , {char_length} }} ;".format(
        pid=i.id_number , x=i.coordinates[0], y=i.coordinates[1] , 
        z=i.coordinates[2], char_length=characteristic_length)
    f.write(point_string + '\n')
f.write("\n\n\n")
f.write("//Construct Edges:\n")
for edge in tube.edges:
    edge_string = "Line( {lid} ) = {{ {pid_1} , {pid_2} }} ;".format( 
        lid=edge.id_number , pid_1=edge.point_ids[0], pid_2=edge.point_ids[1])
    f.write(edge_string + '\n')
f.write("\n\n\n")
f.write("//Construct Faces:\n")
for face in tube.faces:
    face_string = "Line Loop( {id_num} ) = {{ {edge_1} , {edge_2} , {edge_3} }};"\
                  .format(id_num=face.id_number, edge_1=face.edge_ids[0],
                          edge_2=face.edge_ids[1], edge_3=face.edge_ids[2])
    f.write(face_string +'\n')
    surface_string = "Plane Surface( {id_num} ) = {{ {id_num} }} ;".format(
        id_num = face.id_number)
    f.write(surface_string + '\n')
f.write("\n\n\n")
f.write("//Construct Knot Volume:\n")
f.write("Surface Loop(1) = {")
for i in list(range(len(tube.faces))):
    if i % 10 == 0:
        f.write('\n')
    if i != len(tube.faces)-1:
        surface = " {} ,".format(i) 
        f.write(surface)
    else:
        surface = " {} }};\n".format(i) 
        f.write(surface)

#Update 2/12
f.write("Physical Surface(1) = {")
for i in list(range(len(tube.faces))):
    if i % 10 == 0:
        f.write('\n')
    if i != len(tube.faces)-1:
        surface = " {} ,".format(i) 
        f.write(surface)
    else:
        surface = " {} }};\n".format(i) 
        f.write(surface)

boundingSphere = True

## Bounding Sphere 
if boundingSphere:
    f.write("\n\n\n")
    f.write("//Construct Bounding Sphere:\n")
    
    radius = 5
    last_id = tube.points[-1].id_number 
    origin_id = last_id+1 
    x_pos_id = last_id+2
    x_neg_id = last_id+3
    y_pos_id = last_id+4
    y_neg_id = last_id+5
    z_pos_id = last_id+6
    z_neg_id = last_id+7
    
    point_string = "Point ( {pid} ) = {{ {x} , {y} , {z} , {char_length} }} ;".format(
            pid=origin_id , x=0, y=0 , 
            z=0, char_length=characteristic_length)
    f.write(point_string + '\n')
    point_string = "Point ( {pid} ) = {{ {x} , {y} , {z} , {char_length} }} ;".format(
            pid=x_pos_id , x=radius, y=0 , 
            z=0, char_length=characteristic_length)
    f.write(point_string + '\n')
    point_string = "Point ( {pid} ) = {{ {x} , {y} , {z} , {char_length} }} ;".format(
            pid=x_neg_id , x=-radius, y=0 , 
            z=0, char_length=characteristic_length)
    f.write(point_string + '\n')
    point_string = "Point ( {pid} ) = {{ {x} , {y} , {z} , {char_length} }} ;".format(
            pid=y_pos_id , x=0, y=radius , 
            z=0, char_length=characteristic_length)
    f.write(point_string + '\n')
    point_string = "Point ( {pid} ) = {{ {x} , {y} , {z} , {char_length} }} ;".format(
            pid=y_neg_id , x=0, y=-radius , 
            z=0, char_length=characteristic_length)
    f.write(point_string + '\n')
    point_string = "Point ( {pid} ) = {{ {x} , {y} , {z} , {char_length} }} ;".format(
            pid=z_pos_id , x=0, y=0 , 
            z=radius, char_length=characteristic_length)
    f.write(point_string + '\n')
    point_string = "Point ( {pid} ) = {{ {x} , {y} , {z} , {char_length} }} ;".format(
            pid=z_neg_id , x=0, y=0 , 
            z=-radius, char_length=characteristic_length)
    f.write(point_string + '\n')
    
    last_id = tube.edges[-1].id_number 
    arc_1 = last_id+1 
    arc_2 = last_id+2
    arc_3 = last_id+3
    arc_4 = last_id+4
    arc_5 = last_id+5
    arc_6 = last_id+6
    arc_7 = last_id+7
    arc_8 = last_id+8
    
    edge_string = "Circle( {lid} ) = {{ {pid_1} , {pid_2}, {pid_3} }} ;".format( 
            lid=arc_1 , pid_1=x_pos_id, pid_2=origin_id , pid_3=y_pos_id)
    f.write(edge_string + '\n')
    
    edge_string = "Circle( {lid} ) = {{ {pid_1} , {pid_2}, {pid_3} }} ;".format( 
            lid=arc_2 , pid_1=y_pos_id, pid_2=origin_id , pid_3=x_neg_id)
    f.write(edge_string + '\n')
    
    edge_string = "Circle( {lid} ) = {{ {pid_1} , {pid_2}, {pid_3} }} ;".format( 
            lid=arc_3 , pid_1=x_pos_id, pid_2=origin_id , pid_3=y_neg_id)
    f.write(edge_string + '\n')
    
    edge_string = "Circle( {lid} ) = {{ {pid_1} , {pid_2}, {pid_3} }} ;".format( 
            lid=arc_4 , pid_1=y_neg_id, pid_2=origin_id , pid_3=x_neg_id)
    f.write(edge_string + '\n')
    
    edge_string = "Circle( {lid} ) = {{ {pid_1} , {pid_2}, {pid_3} }} ;".format( 
            lid=arc_5 , pid_1=x_pos_id, pid_2=origin_id , pid_3=z_pos_id)
    f.write(edge_string + '\n')
    
    edge_string = "Circle( {lid} ) = {{ {pid_1} , {pid_2}, {pid_3} }} ;".format( 
            lid=arc_6 , pid_1=z_pos_id, pid_2=origin_id , pid_3=x_neg_id)
    f.write(edge_string + '\n')
    
    edge_string = "Circle( {lid} ) = {{ {pid_1} , {pid_2}, {pid_3} }} ;".format( 
            lid=arc_7 , pid_1=x_pos_id, pid_2=origin_id , pid_3=z_neg_id)
    f.write(edge_string + '\n')
    
    edge_string = "Circle( {lid} ) = {{ {pid_1} , {pid_2}, {pid_3} }} ;".format( 
            lid=arc_8 , pid_1=z_neg_id, pid_2=origin_id , pid_3=x_neg_id)
    f.write(edge_string + '\n')
    
    last_id = tube.faces[-1].id_number 
    loop_1_id = last_id+1
    loop_2_id = last_id+2
    loop_3_id = last_id+3
    loop_4_id = last_id+4
    
    face_string = ("Line Loop( {id_num} ) = {{ {edge_1} , {edge_2} , {edge_3} ,"
                   "{edge_4} }};")\
                      .format(id_num=loop_1_id, edge_1=-arc_1, 
                              edge_2=-arc_2, edge_3=arc_8 , edge_4 = arc_7)
    f.write(face_string +'\n')
    
    surface_string = "Ruled Surface( {id_num} ) = {{ {id_num} }} ;".format(
            id_num = loop_1_id)
    f.write(surface_string + '\n')
    
    face_string = ("Line Loop( {id_num} ) = {{ {edge_1} , {edge_2} , {edge_3} ,"
                   "{edge_4} }};")\
                      .format(id_num=loop_2_id, edge_1=-arc_4, 
                              edge_2=arc_6, edge_3=arc_5 , edge_4 = -arc_3)
    f.write(face_string +'\n')
    
    surface_string = "Ruled Surface( {id_num} ) = {{ {id_num} }} ;".format(
            id_num = loop_2_id)
    f.write(surface_string + '\n')
    
    face_string = ("Line Loop( {id_num} ) = {{ {edge_1} , {edge_2} , {edge_3} ,"
                   "{edge_4} }};")\
                      .format(id_num=loop_3_id, edge_1=arc_3, 
                              edge_2=arc_4, edge_3=-arc_8 , edge_4 =-arc_7)
    f.write(face_string +'\n')
    
    surface_string = "Ruled Surface( {id_num} ) = {{ {id_num} }} ;".format(
            id_num = loop_3_id)
    f.write(surface_string + '\n')
    
    face_string = ("Line Loop( {id_num} ) = {{ {edge_1} , {edge_2} , {edge_3} ,"
                   "{edge_4} }};")\
                      .format(id_num=loop_4_id, edge_1=-arc_5, 
                              edge_2=-arc_6, edge_3=arc_2 , edge_4 = arc_1)
    f.write(face_string +'\n')
    
    surface_string = "Ruled Surface( {id_num} ) = {{ {id_num} }} ;".format(
            id_num = loop_4_id)
    f.write(surface_string + '\n')
    
    volume_string = "Surface Loop( 2 ) = {{ {a} , {b} , {c} , {d} }};".format(
       a = loop_1_id, b = loop_2_id , c = loop_3_id , d = loop_4_id )
    f.write(volume_string + '\n')

#update 2/12
    volume_string = "Physical Surface( 2 ) = {{ {a} , {b} , {c} , {d} }};".format(
       a = loop_1_id, b = loop_2_id , c = loop_3_id , d = loop_4_id )
    f.write(volume_string + '\n')
    
    f.write("Volume(1) = { 2, 1};\n")
    #f.write("Volume(2) = { 1};\n")
 
    #Update, no knot volume
    f.write("Physical Volume(2001) = {1};\n")
    #f.write("Physical Volume(2002) = {2};\n")
    
    
    f.write("Coherence;")
    f.close()
    
