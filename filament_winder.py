import numpy as np

diameter = 1000 #mm
thickness = 20 #mm

theta = 53.1 #degrees

E_1 = 40 #MPa
E_2 = 10 #MPa
G_12 = 3.5 #MPa
v_12 = 0.25 #unitless

v_21 = v_12*E_2/E_1

#Q matrix - stiffness matrix of an individual ply.
Q_11 = E_1/(1-v_12*v_21)
Q_12 = v_12*E_2/(1-v_12*v_21)
Q_22 = E_2/(1-v_12*v_21)
Q_66 = G_12

#Q matrix
Q = np.array([[Q_11, Q_12, 0],
                [Q_12, Q_22, 0],
                [0, 0, Q_66]])

print("Q matrix: ", Q, "[GPa]")

def Trans(theta): #transformation matrix
    theta = theta*np.pi/180
    transform_matrix = np.array([[np.cos(theta)**2, np.sin(theta)**2, 2*np.sin(theta)*np.cos(theta)],
                                [np.sin(theta)**2, np.cos(theta)**2, -2*np.sin(theta)*np.cos(theta)],
                                [-np.sin(theta)*np.cos(theta), np.sin(theta)*np.cos(theta), np.cos(theta)**2 - np.sin(theta)**2]])
    
    return transform_matrix

T = Trans(theta)

print("T matrix: ", T)

Q_bar = np.linalg.inv(T) @ Q @ np.linalg.inv(T).T

print("Q_bar matrix: ", Q_bar, "[GPa]")

A = Q_bar * thickness * 1000

print("A matrix: ", A, " *[MPa]")

print(T @ np.linalg.inv(A))

#P*1000/4 * (T @ np.linalg.inv(A))[0][0] + P*1000/2 * (T @ np.linalg.inv(A))[0][1] = 0.001
P = 0.001/(1000/4 * (T @ np.linalg.inv(A))[0][0] + 1000/2 * (T @ np.linalg.inv(A))[0][1])

print("Pressure: ", P, "[MPa]")