function create_cartesian_grid(N, M, Lx, Ly, Nfields)
    dx = Lx/N
    dy = Ly/M
    return dx, dy, zeros(N, M, Nfields)
end