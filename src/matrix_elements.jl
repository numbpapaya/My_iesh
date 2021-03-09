"""
Computes Independent Matrix Elements, i.e. Neutral, Ionic and
Coupling surfaces. Following arguments are needed

- "N": Number of Au Atom
- "x": 2d array of atom positions, shape: (N+2, 3)
- "nn": 2d array of nearest neighbor Au atoms; shape: (N, 12)
- "r0": 2d arrayequilibrium position of 12 nearest neighbor atoms, shape: (3, 12)
- "aPBC": 1d array of box dimensions, shape(3, )
- "mass": 1d array of masses of the nuclei, shape(N+2, )
- "Hp": 1d array of neutral, ion and coupling energies of 3(N+2) nuclei, shape: (3,)
- "dHp": 2d array of neutral, ion and coupling forces of 3(N+2) nuclei, shape: (3*(N+2), 3)
"""

function bond_length_NO(x)
    """
    Computes the bond length, i.e. the distance between the N and O Nuclei
    """
    r_N = x[:, 1]
    r_O = x[:, 2]
    return Euclidean()(r_N, r_O)
end

function zcom_NO(x, mass)
    """
    Computes the z component of the center of mass of the NO projectile
    """
    z_N = x[3, 1]
    z_O = x[3, 2]
    m_N = mass[1]
    m_O = mass[2]
    zcom = (z_N * m_N + z_O*m_O)/(m_N + m_O)
    return zcom
end

function V_au_au(x)
    
end
