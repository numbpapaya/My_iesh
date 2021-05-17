

function bond_length_NO(x)
    """
    Computes the bond length, i.e. the distance between the N and O Nuclei
    """
    r_N = x[1, :]
    r_O = x[2, :]
    return Euclidean()(r_N, r_O)
end

function zcom_NO(x)
    """
    Computes the z component of the center of mass of the NO projectile
    """
    z_N = x[1, 3]
    z_O = x[2, 3]
    zcom = (z_N * m_N + z_O*m_O)/(m_N + m_O)
    return zcom
end




function generate_therm_surf!(surfp, occnum, surfh, surfpinit)
    global eigval_H
    #generate thermalized populations for electronic states
    for j in 1:Ne
        surfp[j] = j
    end

    for j in 1:Ne
        occnum[surfp[j]] = 1
    end
    k = 1
    for j in 1:Ms+1
        if occnum[j] == 0
            surfh[k] = j
            k = k+1
        end
    end
    for t in 1:Ne
        for k in 1:Ne*(Ms + 1 - Ne)
            hoprand = rand(Float64, 3)
            jp = Int(floor(hoprand[1]*Ne)) + 1
            jh = Int(floor(hoprand[2] * (Ms + 1 - Ne))) + 1
            rtemp = exp(-(eigval_H[surfh[jh]] - eigval_H[surfp[jp]])) / (kb * tsurf)
            if hoprand[3] < rtemp
                m = surfh[jh]
                surfh[jh] = surfp[jp]
                surfp[jp] = m
            end
        end
    end
    surfpinit = copy(surfp)
    return surfp, occnum, surfh, surfpinit
end
