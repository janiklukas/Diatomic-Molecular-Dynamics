"""
The function `interaction_force` calculates the force between two atoms located at `loc_i` and `loc_j` for given potential parameters `a`, `b` and `c`.
"""
function interaction_force(loc_i::Vector{Float64}, loc_j::Vector{Float64}, a::Float64, b::Float64, c::Float64)
    r_ij = norm(loc_j - loc_i)

    return 1/a^2 * (c * (r_ij/a)^2 - (b+c)) * exp(-(r_ij/a)^2) * (loc_j - loc_i)
end

"""
The function `wall_force` calculates the force on an atom located at `loc_i` due to the boundary wall potential.
It should only be called if the molecule center of mass is near or outside the boundary.
"""
function wall_force(loc_i::Vector{Float64}, L::Float64, k_wall::Float64)
    out = zeros(2)

    for n in 1:2
        if loc_i[n] < 0 
            out[n] += - k_wall * loc_i[n]
        elseif loc_i[n] > L
            out[n] += - k_wall * (loc_i[n] - L)
        end
    end

    return out
end

"""
The function `assign_grid_locations` takes in the state vector and returns a matrix whose entry (i,j) contains the indices of all molecules with center of mass in cell (i,j).
Particles outside the system boundary are assigned to the outermost cells and thus flagged for wall forces.
"""
function assign_grid_locations(state::Matrix{Float64}, L::Float64, N_cells::Int64) 
    L_cell = L / N_cells

    inds = ceil.(Int, state[:, 1:2] / L_cell)
    inds = min.(inds, N_cells)
    inds = max.(inds, 1)

    out = fill(Vector{Int64}([]), (N_cells, N_cells))
    out = map(copy, out)

    for i in 1:N
        push!(out[inds[i, :]...], i)
    end

    out = permutedims(out)[end:-1:1, :]

    return out
end

"""
The function `find_neighbours` takes in the output of `assign_grid_locations` and returns a list whose i-th entry is a list of all molecules with center of mass in grid cells adjacent to molecule i.
"""
function find_neighbours(grid_locations::Matrix{Vector{Int64}}, N::Int64, N_cells::Int64)
    cells = CartesianIndices((N_cells, N_cells))
    kernel = CartesianIndex(-1, -1):CartesianIndex(1, 1)
    
    out = fill(Vector{Int64}([]), N)
    out = map(copy, out)

    for i in 1:N
        ind = findall(x -> i in x, grid_locations)[1]
        for k in kernel
            if ind+k in cells
                append!(out[i], filter(x -> x != i, grid_locations[ind+k]))
            end
        end
    end

    return sort.(out)
end

"""
The function `accceleration` takes in the full state vector and all necessary parameters and returns a matrix whose i-th row contains the total accceleration for molecule i.
"""
function accceleration(state::Matrix{Float64}, N::Int64, d::Float64, L::Float64, a::Float64, b::Float64, c::Float64, k_wall::Float64, N_cells::Int64)
    grid_locations = assign_grid_locations(state, L, N_cells)
    neighbours = find_neighbours(grid_locations, N, N_cells)

    out = zeros(size(state))

    #Implement the monatomic limit:
    if d > 0
        signs = [1, -1]
    else
        signs = [1]
    end

    for i in 1:N
        #position of one atom relative to the molecule center of mass:
        pos_1 = d/2 * [cos(state[i, 3]), sin(state[i, 3])]

        for σ in signs
            #atom position in global coordinates:
            loc_1 = state[i, 1:2] + σ * pos_1

            force = zeros(2)

            #interaction force on each atom
            for j in neighbours[i]
                pos_2 = d/2 * [cos(state[j, 3]), sin(state[j, 3])]
                    
                for ρ in signs
                    loc_2 = state[j, 1:2] + ρ * pos_2
                    force += interaction_force(loc_1, loc_2, a, b, c)
                end
            end
            
            #wall force on each atom
            ind = findall(x -> i in x, grid_locations)[1]

            if ind[1] in [1, N_cells] || ind[2] in [1, N_cells]
                force += wall_force(loc_1, L, k_wall)
            end

            out[i, 1:2] += force

            #total z-torque on each atom
            if d > 0
                out[i, 3] += 4/d^2 * σ * (pos_1[1] * force[2] - pos_1[2] * force[1])
            end
        end
    end

    return out
end

"""
The function `verlet_step` takes in the previous and current state vectors as well as all necessary parameters and returns the updated state vector
"""
function verlet_step(state_prev::Matrix{Float64}, state::Matrix{Float64}, N::Int64, d::Float64, L::Float64, a::Float64, b::Float64, c::Float64, k_Wall::Float64, N_cells::Int64, Δt::Float64)
    state_next = 2 * state - state_prev + Δt^2 * accceleration(state, N, d, L, a, b, c, k_Wall, N_cells)

    return state_next
end