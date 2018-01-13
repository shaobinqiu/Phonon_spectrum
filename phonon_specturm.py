from numpy import *
import numpy as np
import itertools
import os

def from_string(content):
    # move empty line
    lines = [l for l in content.split('\n') if l.rstrip()]
    comment = lines[0]
    zoom = float(lines[1])
    lattice = np.around(np.array([[float(i) for i in line.split()]
                                                 for line in lines[2:5]]),
                                 decimals=6)
    if zoom < 0:
        # In vasp, a negative scale factor is treated as a volume. We need
        # to translate this to a proper lattice vector scaling.
        vol = abs(np.linalg.det(lattice))
        lattice *= (-zoom / vol) ** (1 / 3)
    else:
        lattice *= zoom
    #nsymbols = [Specie(s).Z for s in lines[5].split()]
    natoms = [int(i) for i in lines[6].split()]
    positions = np.around(np.array([[float(i) for i in line.split()[0:3]]
                                                 for line in lines[8:]]),
                          decimals=6)
    return lattice, positions, natoms


def dist_matrix(latt, pos, sum_atom):#计算距离矩阵
    dis_m = np.zeros((sum_atom, sum_atom))
    for i in range(0,sum_atom):
        for j in range(0,sum_atom):

            radiu = mat(pos[i]-pos[j])*mat(latt)
            diss =math.sqrt((np.array(radiu)**2).sum())
            dis_m[i][j] = diss
    return dis_m

def main():
    M_Z = np.array([28 , 72.5 , 1])#相对质量

    R = 3#截断半径
    with open('CONTCAR.vasp') as f:
        content = f.read()
        latt, pos, numbers = from_string(content)
        pos_cc = np.array(mat(pos)*mat(latt))
    #print(latt, pos, numbers)
    sum_atom = sum(numbers)
    M = np.zeros(sum_atom)#质量
    for m in range(0,sum_atom):
        M[m]=M_Z[0]*bool(m<numbers[0])+M_Z[1]*bool(numbers[0]-1<m<numbers[0]
                 +numbers[1])+M_Z[2]*bool(m>numbers[0]+numbers[1]-1)

    print(M)
    dist_m = dist_matrix(latt, pos, sum_atom)
    #print(dist_m)
    for ax in [0,1,2]:#XYZ方向
        a = np.zeros((sum_atom, sum_atom))
        #print(a)
        for i in range(0,sum_atom):
            for j in range(0,sum_atom):
                if dist_m[i][j]<R and i!=j:#计算一定距离内原子相互作用
                    print(abs(pos_cc[i][ax]-pos_cc[j][ax])/dist_m[i][j])
                    K = 1*abs(pos_cc[i][ax]-pos_cc[j][ax])/dist_m[i][j]/M[i]
                    a[i][j] = a[i][j] + K
                    a[i][i] = a[i][i] - K
    print(a)


main()
