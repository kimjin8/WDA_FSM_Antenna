# 11/11/17
# Kim Jin
# Folded Spherical Meander Antenna
# 3 - Spherical Coordinates, Curl, Multiple arms & Top circle, Export

# %% Import libraries
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d
import scipy.constants
import PyNEC

# %% Equations
def ka():               # electrical size of antenna
    ka = k * a
    return ka


def theta_to_pi(theta):    # theta in pi
    theta = [i * (2*np.pi/360) for i in theta]
    return theta


def phi_pi(phi):        # phi in pi
    phi = [(90 - i) * (2*np.pi/360) for i in phi]
    return phi


def spherical_to_Cartesian(r, theta, phi):
    x = r * np.sin(phi) * np.cos(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(phi)
    return x, y, z


def find_distance(r0, theta0, phi0, r1, theta1, phi1):
    global wire_length
    [x0, y0, z0] = spherical_to_Cartesian(r0, theta0, phi0)
    radius = np.sqrt(a**2 - z0**2)
    wire_length = abs(radius * (theta1-theta0))
    #print(theta0, phi0, theta1, phi1)
    #print(wire_length)
    #distance = np.sqrt(r0**2 + r1**2 - 2*r0*r1*(np.sin(theta0)*np.sin(theta1) \
                #*np.cos(theta0-theta1)+np.cos(theta0)*np.cos(theta1)))
    return wire_length


# %% Functions
def vertices_generator(arm_number):     # arm number starting from 0
    # Get the list of vertices
    # r_vertex, theta_vertex, phi_vertex, in spherical coordinates
    r_vertex = N_vertex * [a]
    
    theta_stair = [gamma/2, gamma/2, gamma/2 + alpha, gamma/2 + alpha] 
    theta_stair_arm = [(i + arm_number*(360/N_arm)) for i in theta_stair]
    theta_vertex = theta_stair_arm * int(N_stair / 2)
    theta_vertex = theta_to_pi(theta_vertex)
    
    phi_vertex_deg = []
    phi_start = 0
    phi_vertex_deg.append(phi_start)
    for i in range(1,N_stair):
        phi_vertex_deg.append(beta*i)
        phi_vertex_deg.append(beta*i)
    phi_vertex_deg.append(beta*N_stair)   
    phi_vertex = phi_pi(phi_vertex_deg)
    
    return r_vertex, theta_vertex, phi_vertex


def init_wire():
    global r_wire, theta_wire, phi_wire
    r_wire = []
    theta_wire = []
    phi_wire = []
    
    
def wire_segmentation(r0, theta0, phi0, r1, theta1, phi1):
    wire_length = find_distance(r0, theta0, phi0, r1, theta1, phi1)
    # find number of segments
    N_segment = int(wire_length / segment_length)
    #print(wire_length)
    print(N_segment)
    # generate segment coordinates & wire segment objects
    r_wire = [a] * (N_segment+1)
    theta_wire = []
    phi_wire = [phi0] * (N_segment+1)
    
    theta_segment = (theta1 - theta0) / N_segment
    for i in range(0, N_segment+1):
        theta_wire.append(theta0 + theta_segment * i)
        
    return r_wire, theta_wire, phi_wire


def plot_wire(r, theta, phi):
    global N_wire
    [x, y, z]= spherical_to_Cartesian(r, theta, phi)
    plot = ax.plot(x, y, z,'b')
    nec_wire(x, y, z)
    

def nec_wire(x, y, z):
    global N_wire
    for i in range(1, len(x)):
        #   GW	1(tag)	 3(segment)	0	0	g	0	0	h	0.0005(thickness)	'v1
        f.write("GW" + "\t" + str(N_wire) + "\t" + "1" + "\t" +  \
                str(x[i-1]) + "\t" + str(y[i-1]) + "\t" + str(z[i-1]) + "\t" + \
                str(x[i]) + "\t" + str(y[i]) + "\t" + str(z[i]) + "\t" + \
                str(R_wire) + "\n")
        N_wire = N_wire + 1


# %% Variables
# Antenna properties
f = 915         # antenna frequency, in MHz
lamda = scipy.constants.c / (f*1e6)
k = (2*np.pi) / lamda   # propagation constant
segment_length = lamda/150

# Antenna space constraints
    
# Antenna geometry parameters
a = 21.15e-3    # radius of the imagnary hemisphere, in meters
alpha = 77.8    # angle of one arm in the x-y plane, in degrees
beta = 6        # angle of one meander segment to the x-y plane, in degrees
R_wire = 0.0005 # radius of the wire, in meters
N_arm = 4       # number of arms
N_stair = 14   # number of parallel segments
N_vertex = N_stair * 2  # number of vertices in an arm
N_wire = 0

gamma = 360/N_arm - alpha       # angle between arms in the x-y plane, in degrees
delta = 90 - N_stair * beta    # phi angle of the circle on the top
R_top = a * (2*np.pi/360) * delta   # radius of the circle on the top

# plot setup
fig = plt.figure()
ax = fig.gca(projection = '3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

# nec file setup
f = open('nec_spherical.txt','w')
f.write("CM\tFolded Spherical Meander Antenna_Numbers\n")   # comment title
f.write("CE\n")


# %% Main Function
# plot arms
for arm_number in range(0, N_arm):
    # generate list of vertices -> wires
    [r_vertex, theta_vertex, phi_vertex] = vertices_generator(arm_number)
   
    # generate wire segments
    for wire_number in range(1, N_vertex):
        # initialize wire coordinates
        init_wire()
        # Vertical element wire coordinates; odd number index
        if (wire_number % 2):
            r_wire =[a, a]      
            theta_wire.append(theta_vertex[wire_number-1])
            theta_wire.append(theta_vertex[wire_number])
            phi_wire.append(phi_vertex[wire_number-1])
            phi_wire.append(phi_vertex[wire_number])
            
            #Wire([a, theta_vertex[i-1], phi_vertex[i-1]], [a, theta_vertex[i], phi_vertex[i]])
            
        # Horizontal element wire coordinates in segments; even number index
        else: 
            # find wire segment coordinates
            [r_wire, theta_wire, phi_wire] = wire_segmentation \
                (a, theta_vertex[wire_number-1], phi_vertex[wire_number-1], 
                 a, theta_vertex[wire_number], phi_vertex[wire_number])
        # add to the plot 
        plot_wire(r_wire, theta_wire, phi_wire)


# plot top circle
# initialize wire coordinates
init_wire()
segment_length = lamda/200
[r_wire, theta_wire, phi_wire] = wire_segmentation \
    (a, 0, phi_vertex[-1], 
     a, 2*np.pi, phi_vertex[-1])
# add to the plot 
plot_wire(r_wire, theta_wire, phi_wire)


# show plot
plt.show()
# print(N_wire)
# finish nec file
f.write("GE\t0\n")
#f.write("GN\t")
f.write("EK\n")
f.write("EX\t0\t0\t1\0\t1\t0\n")
f.write("FR\t0\t1\t0\t0\t915\t0\n")
f.write("EN")
f.close()