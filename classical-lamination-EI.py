import numpy as np

def Trans(theta): #transformation matrix
    theta = theta*np.pi/180
    transform_matrix = np.array([[np.cos(theta)**2, np.sin(theta)**2, 2*np.sin(theta)*np.cos(theta)],
                                [np.sin(theta)**2, np.cos(theta)**2, -2*np.sin(theta)*np.cos(theta)],
                                [-np.sin(theta)*np.cos(theta), np.sin(theta)*np.cos(theta), np.cos(theta)**2 - np.sin(theta)**2]])
    
    return transform_matrix

def summarize_layers(layers):
    for i in layers:
        print("{}: Name: {} Thickness: {}mm Angle: {}".format(i, layers[i]['name'], layers[i]['t'], layers[i]['angle']))
    return

#define layers and properties
layers = {
    0: {'name': 'concrete', 'angle': 0, 't': 100, 'E1': 25000, 'E2': 25000, 'v12': 0.15, 'G12': 25000/(2*(1+0.15))},
}

#define CLT layers
#odd layers have a higher quality of wood than even layers
#odd layers are longitudinal, even layers are transversal
for i in range(1,6):
    name = 'CLT layer ' + str(i)
    if i % 2 != 0: #odd number - longitudinal layer
        E1 = 11700 #MPa
        E2 = 11700/30 #MPa

        angle = 0 #degreees
    else:
        E1 = 9000 #MPa
        E2 = 9000/30 #MPa

        angle = 90 #degrees

    G12 = E1/16 #based on CSA 086 approximation

    v12 = 0.449 #USDA value for Douglas Fir for mu-LT (longitudinal to transversal)

    layers[i] = {'name': name, 'angle': angle, 't': 35, 'E1': E1, 'E2': E2, 'v12': v12, 'G12': G12}

#Split layer dict down the centre. Middle layer is split into two layers
def split_layer_dict(layers):
    total_thickness = sum([layer['t'] for layer in layers.values()])
    t = 0
    
    for i in range(0, len(layers)):
        
        if layers[i]['t'] + t > total_thickness/2:
            new_layer_1 = layers[i].copy()
            new_layer_2 = layers[i].copy()

            new_layer_1['t'] = total_thickness/2 - t
            new_layer_2['t'] = layers[i]['t'] + t - total_thickness/2

            new_dict = {}
            
            j = 0
            while j < i:
                new_dict[j] = layers[j]
                j += 1
            
            new_dict[j] = new_layer_1
            new_dict[j+1] = new_layer_2

            while j < len(layers)-1:
                new_dict[j+2] = layers[j+1]
                j += 1

            break

        t += layers[i]['t']
    
    #edge case: if there are an even number of layers, this creates a layer with 0 thickness. delete it.
    layers = new_dict
    if layers[i]['t'] == 0:
        del layers[i]
        new_dict = {}
        i = 0
        for key in layers.keys():
            new_dict[i] = layers[key]
            i += 1
    
        layers = new_dict

    return layers

layers = split_layer_dict(layers)

Q_matrices = []
Q_bar_matrices = []

for i in range(0, len(layers)):
    #Q matrix - stiffness matrix of an individual ply.
    E_1 = layers[i]['E1']
    E_2 = layers[i]['E2']
    v_12 = layers[i]['v12']
    G_12 = layers[i]['G12']
    v_21 = (v_12 * E_2) / E_1

    Q_11 = E_1 / (1 - v_12 * v_21)
    Q_12 = v_12 * E_2 / (1 - v_12 * v_21)
    Q_22 = E_2 / (1 - v_12 * v_21)
    Q_66 = G_12

    Q = np.array([[Q_11, Q_12, 0],
                    [Q_12, Q_22, 0],
                    [0, 0, Q_66]])
    
    Q_matrices.append(Q)
    
    T = Trans(layers[i]['angle'])
    
    Q_bar = np.linalg.inv(T) @ Q @ np.linalg.inv(T).T

    Q_bar_matrices.append(Q_bar)

#Calculate A, B, and D matrices.

#begin by calculating z_k for each layer

z_middle = sum([layer['t'] for layer in layers.values()])/2

for i in range(0, len(layers)):
    layers[i]['z'] = z_middle - sum([layer['t'] for layer in layers.values()][:i])

#z values to construct A, B, and D matrices
z = []

total_thickness = sum([layer['t'] for layer in layers.values()])

z.append(-total_thickness/2)

for i in range(1, len(layers)+1):
    z.append(z[i-1] + layers[i-1]['t'])

#define blank A, B, and D matrices
A = np.zeros((3,3))
B = np.zeros((3,3))
D = np.zeros((3,3))

#define A matrix
for i in range(0, 3):
    for j in range(0, 3):
        for k in range(1, len(layers)+1):
            A[i][j] += Q_bar_matrices[k-1][i][j] * (z[k] - z[k-1])

#define B matrix
for i in range(0, 3):
    for j in range(0, 3):
        for k in range(1, len(layers)+1):
            B[i][j] += Q_bar_matrices[k-1][i][j] * (z[k]**2 - z[k-1]**2)

B = B/2

#define D matrix
for i in range(0, 3):
    for j in range(0, 3):
        for k in range(1, len(layers)+1):
            D[i][j] += Q_bar_matrices[k-1][i][j] * (z[k]**3 - z[k-1]**3)

D = D/3

#call the stiffness matrix K
#concatenate A, B, and D matrices
K = np.concatenate((np.concatenate((A, B), axis=1), np.concatenate((B, D), axis=1)), axis=0)

#define stress vector. These are arbitrary given that we later divide the result vector by Mx
Mx = 20

#vertical array of 0, 0, 0, Mx, 0, 0
stress_vector = np.array([[0], [0], [0], [Mx], [0], [0]])

reaction_vector = np.linalg.inv(K) @ stress_vector

kx = reaction_vector[3][0]

EI = Mx/kx * 1000

##output:
print("Layer summary:")
summarize_layers(layers)

print("Q - concrete:")
print(Q_matrices[0])

print("Q - CLT - longitudinal:")
print(Q_matrices[1])

print("Q - CLT - transversal:")
print(Q_matrices[2])

#print Q bar matrices
print("Q bar - concrete:")
print(Q_bar_matrices[0])

print("Q bar - CLT - longitudinal:")
print(Q_bar_matrices[1])

print("Q bar - CLT - transversal:")
print(Q_bar_matrices[2])

print('Stiffness Matrix, K:')
print(K)

print('Reaction Vector:')
print(reaction_vector)

print(kx)
print("EI: {:.2e}".format(EI))

w = 2.4 #MPa/1000mm
L = 9000 #mm
delta = 5*w*L**4 / (384*EI)
print("Deflection: {}mm".format(delta))