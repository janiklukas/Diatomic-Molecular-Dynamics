"""
The function `motion_energy` calculates the total kinetic and rotational energy from the current and previous system state.
"""
function motion_energy(state_prev::Matrix{Float64}, state::Matrix{Float64}, N::Int64, d::Float64, Δt::Float64)
    kinetic_x::Float64 = 0
    kinetic_y::Float64 = 0
    rotational::Float64 = 0

    for i in 1:N
        v_x = (state[i, 1] - state_prev[i, 1]) / Δt
        v_y = (state[i, 2] - state_prev[i, 2]) / Δt
        kinetic_x += 1/2 * v_x^2
        kinetic_y += 1/2 * v_y^2

        ω = (state[i, 3] - state_prev[i, 3]) / Δt
        rotational += d^2/8 * ω^2
    end

    return kinetic_x, kinetic_y, rotational
end

"""
The function `interaction_potential` calculates the potential energy due to the interaction between two atoms located at `loc_i` and `loc_j`.
"""
function interaction_potential(loc_i::Vector{Float64}, loc_j::Vector{Float64}, a::Float64, b::Float64, c::Float64)
    r_ij = norm(loc_j - loc_i)

    return 1/2 * (b - c * (r_ij/a)^2) * exp(-(r_ij/a)^2)
end

"""
The function `plot_interaction` creates a plot of the atom interaction potential.
"""
function plot_interaction(r_max::Float64, a::Float64, b::Float64, c::Float64)
    r_array = range(start = 0.0, stop = r_max, length = 100)
    V_array = @. 1/2 * (b - c * (r_array/a)^2) * exp(-(r_array/a)^2)

    plot(r_array, V_array, label = nothing, size = (400,300))
    xlabel!(L"r/r_0"); ylabel!(L"V/V_0")
    xlims!(0, r_max)
end

"""
The function `wall_potential` calculates the wall potential energy for an atom located at `loc_i`.
"""
function wall_potential(loc_i::Vector{Float64}, L::Float64, k_wall::Float64)
    out::Float64 = 0

    for n in 1:2
        if loc_i[n] < 0 
            out += k_wall/2 * loc_i[n]^2
        elseif loc_i[n] > L
            out += k_wall/2 * (loc_i[n] - L)^2
        end
    end

    return out
end

"""
The function `potential_energy` calculates the total interaction and wall energy for the current system state.
"""
function potential_energy(state::Matrix{Float64}, N::Int64, d::Float64, a::Float64, b::Float64, c::Float64, k_wall::Float64)
    grid_locations = assign_grid_locations(state, L, N_cells)
    neighbours = find_neighbours(grid_locations, N, N_cells)

    if d > 0
        signs = [1, -1]
    else
        signs = [1]
    end

    interaction_energy::Float64 = 0
    wall_energy::Float64 = 0

    for i in 1:N
        pos_1 = d/2 * [cos(state[i, 3]), sin(state[i, 3])]

        for σ in signs
            loc_1 = state[i, 1:2] + σ * pos_1

            #interaction potential
            for j in neighbours[i]
                pos_2 = d/2 * [cos(state[j, 3]), sin(state[j, 3])]
                
                for ρ in signs
                    loc_2 = state[j, 1:2] + ρ * pos_2
                    interaction_energy += 1/2 * interaction_potential(loc_1, loc_2, a, b, c)
                    #factor of 1/2 to avoid double counting!
                end
            end
            
            #wall potential
            ind = findall(x -> i in x, grid_locations)[1]

            if ind[1] in [1, N_cells] || ind[2] in [1, N_cells]
                wall_energy += wall_potential(loc_1, L, k_wall)
            end
        end
    end

    return interaction_energy, wall_energy
end

"""
The function `plot_energy` creates a plot of the various forms of energy in the system over time.
"""
function plot_energy(times::Vector{Float64}, E_kin_x::Vector{Float64}, E_kin_y::Vector{Float64}, E_rot::Vector{Float64}, E_int::Vector{Float64}, E_wall::Vector{Float64})
    total = E_kin_x + E_kin_y + E_rot + E_int + E_wall
    p = plot(times, total, c = 5, label = "total")
    xlabel!(L"t/t_0"); ylabel!(L"E/V_0")

    plot!(times, E_kin_x, c = 1, label = "kinetic (x)")
    plot!(times, E_kin_y, c = 7, label = "kinetic (y)")
    plot!(times, E_rot, c = 2, label = "rotational")
    plot!(times, E_int, c = 3, label = "interaction")
    plot!(times, E_wall, c = 4, label = "wall")

    display(p)
end

"""
The function `evolve` solves the equations of motion for a given initial state and returns an animation of the result.
# Arguments:
- `state_init::Matrix{Float64}` : the initial state vector
- `vel_init` : the initial velocity vector
- `N::Int64` : the number of molecules
- `d::Float64` : the distance between molecule atoms
- `L::Float64` : the box side length
- `a::Float64`, `b::Float64`, `c::Float64` : the interaction potential parameters
- `k_wall::Float64` : the wall spring constant
- `N_cells::Int64` : the number of cell divisions in each dimension used for interaction calculations
- `Δt::Float64`, `t_max::Float64` : the time step and total evolution time
- `fps::Int64`, `len::Int64` : the animation frames per second and length
- `name::String` : the filename for the animation
"""
function evolve(state_init::Matrix{Float64}, vel_init::Matrix{Float64}, N::Int64, d::Float64, L::Float64, a::Float64, b::Float64, c::Float64, k_wall::Float64, N_cells::Int64, Δt::Float64, t_max::Float64, fps::Int64, len::Int64, name::String; create_plot::Bool=true)
    @assert d >= 0 "`d` must be positive"
    @assert a != 0 "`a` cannot be zero"
    @assert k_wall >= 0 "`k_wall` must be greater than zero"
    @assert N_cells > 0 "`N_cells`must be a positive integer"
    @assert L / N_cells >= 3*a "`N_cells` should be chosen such that the cell size is larger than 3*`a`"
    
    state_prev = state_init
    state = state_init + Δt * vel_init + Δt^2/2 * accceleration(state_init, N, d, L, a, b, c, k_wall, N_cells)

    N_steps = ceil(t_max/Δt)

    if create_plot
        times = Vector{Float64}([])
        E_kin_x = Vector{Float64}([])
        E_kin_y = Vector{Float64}([])
        E_rot = Vector{Float64}([])
        E_int = Vector{Float64}([])
        E_wall = Vector{Float64}([])
    end

    anim = @animate for k in 1:N_steps
        scatter(leg = false, aspect_ratio = :equal, xticks = [0, L], yticks = [0, L], size = (400,400))
        scatter!(minorticks = N_cells, minorgrid = true, minorgridalpha = 0.3, minorgridcolor = "gray")
        hline!([0, L], c = 1); vline!([0, L], c = 1)
        xlabel!(L"x/r_0"); ylabel!(L"y/r_0")
        xlims!(0, L); ylims!(0, L)
        title!("Step " * string(round(Int, (k-1))))
        
        scatter!(state[:, 1] + d/2 * cos.(state[:, 3]), state[:, 2] + d/2 * sin.(state[:, 3]), c = 2)
        if d > 0 
            scatter!(state[:, 1] - d/2 * cos.(state[:, 3]), state[:, 2] - d/2 * sin.(state[:, 3]), c = 2)
        end

        if create_plot
            push!(times, (k-1)*Δt)

            e1,e2,e3 = motion_energy(state_prev, state, N, d, Δt)
            push!(E_kin_x, e1)
            push!(E_kin_y, e2)
            push!(E_rot, e3)

            e1, e2 = potential_energy(state, N, d, a, b, c, k_wall)
            push!(E_int, e1)
            push!(E_wall, e2)
        end
        
        state_next = verlet_step(state_prev, state, N, d, L, a, b, c, k_wall, N_cells, Δt)
        state_prev = state
        state = state_next
    end every round(Int, 1/(len*fps) * N_steps)

    display(gif(anim, name * ".gif", fps = fps))

    if create_plot
        plot_energy(times, E_kin_x, E_kin_y, E_rot, E_int, E_wall)
    end

end