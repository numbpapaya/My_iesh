export zcom_NO
function zcom_NO(x)::Float64
    """
    Computes the z component of the center of mass of the NO projectile
    """
    z_N = x[1, 3]
    z_O = x[2, 3]
    zcom = (z_N * m_N + z_O*m_O)/(m_N + m_O)
    return zcom
end
