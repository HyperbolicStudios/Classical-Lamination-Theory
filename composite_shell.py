import numpy as np

#analyzing the strains on a compressed gas tank, made from carbon/epoxy at [30, -30] layup
#Shell has n symmetrical layers - find strain as a function of n

#tank properties
diameter = 300 #mm
pressure = 25 #MPa

#laminate properties
ply_t = 0.127 #mm
E_1 = 138 #GPa
E_2 = 8.96 #GPa
G_12 = 7.1  #GPa
v_12 = 0.3 #GPa

theta = 30 #degrees - angle of the ply wrt the axis of the tank (+/-)

v_21 = v_12*E_2/E_1

#hoop stress calculations
#x is along the central axis of the tank, y is perpendicular to the axis
N_x = pressure*diameter/4
N_y = pressure*diameter/2

print("N_x: ", N_x, " [N]")
print("N_y: ", N_y, " [N]")

stress_vector = np.array([[N_x], [N_y], [0]])

#Q matrix - stiffness matrix of an individual ply
Q_11 = E_1/(1-v_12*v_21)
Q_12 = v_12*E_2/(1-v_12*v_21)
Q_22 = E_2/(1-v_12*v_21)
Q_66 = G_12

Q = np.array([[Q_11, Q_12, 0],
                [Q_12, Q_22, 0],
                [0, 0, Q_66]])

print("Q matrix: ", Q, "[GPa]")

def T(theta): #transformation matrix
    theta = theta*np.pi/180
    transform_matrix = np.array([[np.cos(theta)**2, np.sin(theta)**2, 2*np.sin(theta)*np.cos(theta)],
                                [np.sin(theta)**2, np.cos(theta)**2, -2*np.sin(theta)*np.cos(theta)],
                                [-np.sin(theta)*np.cos(theta), np.sin(theta)*np.cos(theta), np.cos(theta)**2 - np.sin(theta)**2]])
    
    return transform_matrix

#Q_bar = T**-1 @ Q @ T**-T
Q_bar_plus = np.linalg.inv(T(theta)) @ Q @ np.linalg.inv(T(theta)).T
Q_bar_minus = np.linalg.inv(T(-1*theta)) @ Q @ np.linalg.inv(T(-1*theta)).T

print("Q_bar_plus: ", Q_bar_plus, "[GPa]")
print("Q_bar_minus: ", Q_bar_minus, "[GPa]")

#Calculate the 'A' matrix
#calculate A for a single n (where N = 4n) - i.e. 1n is a set of 4 symmetrical layers

#for layers of equal thickness, A = n * 2 (Q_bar_plus * thickness + Q_bar_minus * thickness)

A = 2*(Q_bar_plus*ply_t + Q_bar_minus*ply_t)

print("A: ", A, "[GPa * mm]")

#Apparently B is zero, so we don't need it. And we don't need D either.

#so stress_vector = n* A * strain_vector
#and strain vector = 1/n * stress_vector * A**-1

A = A*1000 #(convert from GPA * mm to MPa * mm)

strain_vector = np.linalg.inv(A) @ stress_vector

print("Strain vector: 1/n * ", strain_vector)