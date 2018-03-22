#################################################

##########Class Definition for simplex class#####

#################################################

class simplex:
    'Simplex has points, edges, and faces'

    def __init__(self):
        self.point_count = 0
        self.points = []
        self.edge_count=0
        self.edges = []
        self.face_count = 0
        self.faces = []
   #Create an empty list of point objects, each point has coordinates and a number

    def addpoint(self, coordinates, rot_mat, state):
        dummy_point = point(self.point_count, coordinates, rot_mat, state)
        self.points.append(dummy_point)
        self.point_count += 1
    
    def addedge(self, point_1, point_2, state):
        dummy_edge = edge(self.edge_count+1, point_1, point_2, state)
        self.edges.append(dummy_edge)
        self.edge_count += 1

    def addface(self, edge_1_id, edge_2_id, edge_3_id):
        dummy_face = face(self.face_count, edge_1_id, edge_2_id, edge_3_id)
        self.faces.append(dummy_face)
        self.face_count += 1
    
    def get_points_with_state(self, state):
        list_of_points=[]
        for p in self.points:
            if p.state_point == state:
                list_of_points.append(p)
        return list_of_points
   
    def get_edges_with_state(self, state):
        list_of_edges=[]
        for e in self.edges:
            if e.state_edge == state:
                list_of_edges.append(e)
        return list_of_edges

class point:
    'Point has three ordered  x,y,z coordinates stored in list "coordinates" and has a single id number according to the order it was generated'
    
    def __init__(self, id_num, coords, rot_mat, state_point):
        self.id_number=id_num
        self.coordinates = coords
        self.rotation_matrix = rot_mat
        self.state_point = state_point
#Here, state_point specifies the state of a point

class edge:
    'edge consists of an id number and two (directed) points, each one having an id number, a negative specifying to reverse orientation when constructing faces'
    
    def __init__(self, id_num, point_1, point_2, state_edge):
        self.id_number=id_num
        self.point_ids=[point_1.id_number, point_2.id_number] 
        self.state_edge=state_edge

class face:
    'face consists of an id number and three (directed) edges, each one having an id number, a negative specifying to reverse orientation when constructing faces'
    
    def __init__(self, id_num, edge_1_id, edge_2_id, edge_3_id):
        self.id_number=id_num
        self.edge_ids=[edge_1_id, edge_2_id , edge_3_id] 
