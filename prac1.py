#Program to determine Huckel energies and degeneracies for an n length
#linear or cyclic polyene and the Platonic solids: cube, tetrahedron and dodecahedron 

#Assuming alpha = 0, beta = -1

import numpy as np

#Prints eigenvalues for a defined matrix
def get_evals(matrix):
    evals, evecs=np.linalg.eig(matrix)            
    print(sorted(evals))


#Produces an n by n Huckel matrix for a length n poly-ene
def lin(size):
    n=int(size)
    m = np.zeros((n,n))
    b = 0
    c = 1
    while b != n and c != n:
            m[b][c] = -1
            m[c][b] = -1
            b += 1
            c += 1
    return m
    

#Produces an n by n Huckel marix for a length n cyclic polyene
def cyc(size):
    n = int(size)
    m = np.zeros((n,n))
    b = 0
    c = 1
    while b != n and c != n:
            m[b][c] = -1
            m[c][b] = -1
            b += 1
            c += 1
    l = int(n-1)
    m[0][l] = m[l][0] = -1
    return m
    

#Working out degeneracy of each eigenvalue
#Degeneracies upto 2
def degeneracy(matrix):
    evals, evecs=np.linalg.eig(matrix)
    ei = sorted(evals) #sorted eigenvalues
    n = len(evals)
    a = 0
    b = 1
    while a < n and b < n:
        if ei[a] - ei[b] > -0.005:
            print("Degeneracy of orbitals", a+1, "and", b+1, "is 2")
        elif ei[a-1] - ei[a] < -0.005 or a == 0:
            print("Degeneracy of orbital", a+1, "is 1")
        a += 1
        b += 1
    if ei[n-2] - ei[n-1] < -0.005:
         print("Degeneracy of orbital", n, "is 1")

#Degeneracies upto 5
def degeneracy2(matrix):
    evals, evecs=np.linalg.eig(matrix)
    ei = sorted(evals) #sorted eigenvalues
    n = len(evals)
    if ei[0] -ei[1] < -0.005:
        print("Degeneracy of orbital 1 is 1")
    a = 0
    b = 1
    c = 2
    d = 3
    e = 4
    while e < n:
        if ei[a] - ei[e] > -0.005:
            print("Degeneracy of orbitals", a+1, ",", b+1, ",", c+1, ",", d+1, "and", e+1, "is 5")
        elif ei[a] - ei[d] > -0.005 and ei[a-1] - ei[a] < -0.005:
            print("Degeneracy of orbitals", a+1, ",", b+1, ",", c+1, "and", d+1, "is 4")
        elif ei[a] - ei[c] > -0.005 and ei[a-1] - ei[a] < -0.005:
            print("Degeneracy of orbitals", a+1, ",", b+1, "and", c+1, "is 3")
        elif ei[a] - ei[b] > -0.005 and ei[a-1] - ei[a] < -0.005:
            print("Degeneracy of orbitals", a+1, "and", b+1, "is 2")
            escape = 1
        a += 1
        b += 1
        c += 1
        d += 1
        e += 1   
    if ei[n-5] - ei[n-1] > -0.005 and ei[n-6] - ei[n-5] < -0.005:
        print("Degeneracy of orbitals", n-4, ",", n-3, ",", n-2, ",", n-1, "and", n, "is 5")
    if ei[n-5] - ei[n-2] > -0.005 and ei[n-6] - ei[n-5] < -0.005 and ei[n-2] - ei[n-1] < -0.005:
        print("Degeneracy of orbitals", n-4, ",", n-3, ",", n-2, "and", n-1, "is 4")
    if ei[n-5] - ei[n-3] > -0.005 and ei[n-6] - ei[n-5] < -0.005 and ei[n-3] - ei[n-2] < -0.005:
        print("Degeneracy of orbitals", n-4, ",", n-3, "and", n-2, "is 3")
    if ei[n-5] - ei[n-4] > -0.005 and ei[n-6] - ei[n-5] < -0.005 and ei[n-4] - ei[n-3] < -0.005 and escape != 1:
        print("Degeneracy of orbitals", n-4, "and", n-3, "is 2")
    if ei[n-4] - ei[n-1] > -0.005:
        print("Degeneracy of orbitals", n-3, ",", n-2, ",", n-1, "and", n, "is 4")
    if ei[n-4] - ei[n-2] > -0.005 and ei[n-5] - ei[n-4] < -0.005 and ei[n-2] - ei[n-1] < -0.005:
        print("Degeneracy of orbitals", n-3, ",", n-2, "and", n-1, "is 3")
    if ei[n-4] - ei[n-3] > -0.005 and ei[n-5] - ei[n-4] < -0.005 and ei[n-3] - ei[n-2] < -0.005:
        print("Degeneracy of orbitals", n-3, "and", n-2, "is 2")
    if ei[n-3] - ei[n-1] > -0.005:
        print("Degeneracy of orbitals", n-2, ",", n-1, "and", n, "is 3")
    if ei[n-3] - ei[n-2] > -0.005 and ei[n-4] - ei[n-3] < -0.005 and ei[n-2] - ei[n-1] < -0.005:
        print("Degeneracy of orbitals", n-2, "and", n-1, "is 2")
    if ei[n-2] - ei[n-1] < -0.005:
        print("Degeneracy of orbital", n, "is 1")

    
#Degeneracies of sp2-hybridised Platonic solids
#Tetrahedron
tet = np.zeros((4,4))
tet[0][1]=tet[1][0]=tet[2][0]=tet[0][2]=tet[3][0]=tet[0][3]=tet[1][2]=tet[2][1]=tet[1][3]=tet[3][1]=tet[2][3]=tet[3][2]=-1


#Cube
cube = np.zeros((8,8))
cube[0][1]=cube[1][0]=cube[0][3]=cube[3][0]=cube[0][4]=cube[4][0] = -1
cube[1][2]=cube[2][1]=cube[1][5]=cube[5][1] = -1
cube[2][3]=cube[3][2]=cube[2][6]=cube[6][2] = -1
cube[3][7]=cube[7][3] = -1
cube[4][5]=cube[5][4]=cube[4][7]=cube[7][4] = -1
cube[5][6]=cube[6][5] = -1
cube[6][7]=cube[7][6] = -1


#Dodecahedron
dod = np.zeros((20,20))
dod[0][1]=dod[1][0]=dod[0][4]=dod[4][0]=dod[0][5]=dod[5][0] = -1
dod[1][2]=dod[2][1]=dod[1][7]=dod[7][1] = -1
dod[2][3]=dod[3][2]=dod[2][9]=dod[9][2] = -1
dod[3][4]=dod[4][3]=dod[3][11]=dod[11][3] = -1
dod[4][13]=dod[13][4] = -1
dod[5][6]=dod[6][5]=dod[5][14]=dod[14][5] = -1
dod[6][7]=dod[7][6]=dod[6][16]=dod[16][6] = -1
dod[7][8]=dod[8][7] = -1
dod[8][9]=dod[9][8]=dod[8][17]=dod[17][8] = -1
dod[9][10]=dod[10][9] = -1
dod[10][11]=dod[11][10]=dod[10][18]=dod[18][10] = -1
dod[11][12]=dod[12][11] = -1
dod[12][13]=dod[13][12]=dod[12][19]=dod[19][12] = -1
dod[13][14]=dod[14][13] = -1
dod[14][15]=dod[15][14] = -1
dod[15][16]=dod[16][15]=dod[15][19]=dod[19][15] = -1
dod[16][17]=dod[17][16] = -1
dod[17][18]=dod[18][17] = -1
dod[18][19]=dod[19][18] = -1


#User interface
system = input("Type of system? ie. linear, cyclic or Platonic \n")
number = int(input("How many carbons? Please use numbers \n"))

if system.lower() == "linear":
    print(lin(number))
    get_evals(lin(number))
    degeneracy(lin(number))

elif system.lower() == "cyclic":
    print(cyc(number))
    get_evals(cyc(number))
    degeneracy(cyc(number))

elif system.lower() == "platonic" and number == 4:
    print(tet)
    get_evals(tet)
    degeneracy2(tet)

elif system.lower() == "platonic" and number == 8:
    print(cube)
    get_evals(cube)
    degeneracy2(cube)

elif system.lower() == "platonic" and number == 20:
    print(dod)
    get_evals(dod)
    degeneracy2(dod)
    
else:
    print("That is not a defined configuration")
               






