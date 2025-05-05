using Random, Plots, Statistics

#set seed for reproducibility
# Random.seed!(1234)

# --- Parameters ---
L = 30
T = 1000
ρ = 0.5
v_0 = 0.5
cell_length = 1
interaction_radius = 1
dt = 1
N = Int(ρ * L^2)
noise_lims = (-pi/2, pi/2) # limits for the white noise
η = 0.1
num_cells = Int(div(L, cell_length)) #number of cells in one direction
total_cells = num_cells * num_cells
frameskip = 2 #print every 200 time steps
# --- Functions ---
# Define a custom mutable particle struct
mutable struct Particle
    id::Int
    x::Float64
    y::Float64
    θ::Float64
    v::Float64
    cell_index::Tuple{Int,Int}  # (i, j) cell coordinates (1-indexed)
    x_p::Float64
    y_p::Float64
    θ_p::Float64
    cell_index_p::Tuple{Int,Int}  # pre-update cell coordinates (1-indexed)
end
# Given x, y coordinates, return the cell indices (1-indexed)
function cell_index(x, y, L, cell_length)
    i = Int(floor(x/cell_length))+1
    j = Int(floor(y/cell_length))+1
    # in a single time step, the particles should not move more than one cell
    # at boundary conditions, wrap around, also x, y are always >=0 so i, j >= 1
    if i >= L+1
        i = 1
    end
    if j >= L+1
        j = 1
    end
    return (i, j)
end
# Initialize particles
function initialize_particles(N, L, v_0, cell_length)
    particles = Vector{Particle}(undef, N)
    for i in 1:N
        x = rand() * L
        y = rand() * L
        θ = 2π * rand()
        cell = cell_index(x, y, L, cell_length)
        x_p = x
        y_p = y
        θ_p = θ
        cell_p = cell
        #p is the pre-update quantity
        particles[i] = Particle(i, x, y, θ, v_0, cell, x_p, y_p, θ_p, cell_p)
    end
    return particles
end
# Convert cell (i, j) to 1D index
function cell_to_index(cell::Tuple{Int,Int})
    i, j = cell
    return i + (j - 1) * num_cells
end

function build_cell_lookup(particles, cell_lookup, L, cell_length)

    # clear all cells
    for i in eachindex(cell_lookup)
        empty!(cell_lookup[i])
    end

    for particle in particles
        index = cell_to_index(particle.cell_index)
        push!(cell_lookup[index], particle.id)
    end
    return cell_lookup
end


# Periodic distance function
function periodic_distance(x1, y1, x2, y2, L)
    dx = abs(x1 - x2)
    dy = abs(y1 - y2)
    dx = min(dx, L - dx)
    dy = min(dy, L - dy)
    return sqrt(dx^2 + dy^2)
end
#function breaks for 2*2 cells because of overlapping cells, i-1 is the same as i+1
function find_neighbors(particle::Particle, particles, reference_cell_lookup, L, cell_length, interaction_radius)
    i, j = particle.cell_index
    neighbors = Int[]

    for di in -1:1, dj in -1:1
        i_neighbor = mod1(i + di, num_cells)
        j_neighbor = mod1(j + dj, num_cells)
        idx = cell_to_index((i_neighbor, j_neighbor))
    
        for neighbor_id in reference_cell_lookup[idx]
            neighbor = particles[neighbor_id]
            if periodic_distance(particle.x_p, particle.y_p, neighbor.x_p, neighbor.y_p, L) < interaction_radius #use the pre-update positions
                push!(neighbors, neighbor_id)
            end
        end
    end
    return neighbors
end
# Velocity vector
function velocity_vector(v, θ)
    return (v*cos(θ), v*sin(θ))
end
# White noise generator
function WhiteNoise(η, noise_lims)
    #generate uniform white White
    #η is the Control Parameter Noise Amplitude
    lower_limit, upper_limit = noise_lims
    return η*(rand()*(upper_limit - lower_limit) + lower_limit)
end

function OrderParameters(particles)
    # Calculate the polar order parameter
    sum_sin = sum(sin(p.θ) for p in particles)
    sum_cos = sum(cos(p.θ) for p in particles)
    P_t = sqrt(sum_sin^2 + sum_cos^2) / length(particles)  
    # Calculate the nematic order parameter
    sum_sin2theta = sum(sin(2 * p.θ) for p in particles)
    sum_cos2theta = sum(cos(2 * p.θ) for p in particles)
    S_t = sqrt(sum_sin2theta^2 + sum_cos2theta^2) / length(particles)
    return P_t, S_t
end

function SignOfCosine(θ)
    if cos(θ) > 0
        return 1
    else
        return -1
    end
end

# Update a single particle
function update_particle!(particle::Particle, particles, reference_cell_lookup, L, cell_length, interaction_radius, η, noise_lims, dt)
    # Find neighbors
    neighbors = find_neighbors(particle, particles, reference_cell_lookup, L, cell_length, interaction_radius)
   
    # Update angle using Viscek Rules
    sum_sin = 0
    sum_cos = 0
    for n in neighbors
        neighbor = particles[n]
        sign = SignOfCosine(neighbor.θ_p - particle.θ_p)
        sum_sin += sign * sin(neighbor.θ_p)
        sum_cos += sign * cos(neighbor.θ_p)
    end
    θ_new = atan(sum_sin, sum_cos) + WhiteNoise(η, noise_lims)

    # Update position
    v_x, v_y = velocity_vector(particle.v, θ_new)
    x_new = mod(particle.x + v_x*dt, L)
    y_new = mod(particle.y + v_y*dt, L)
    new_cell = cell_index(x_new, y_new, L, cell_length)

    # Update particle data
    particle.x, particle.y, particle.θ, particle.cell_index = x_new, y_new, θ_new, new_cell
end
#Update pre-update values of all particles

function update_pre_update_values!(particles)
    Threads.@threads for i in eachindex(particles)
        particle = particles[i]
        particle.x_p, particle.y_p, particle.θ_p, particle.cell_index_p = particle.x, particle.y, particle.θ, particle.cell_index
    end
end

# Parallelized Monte Carlo step
function monte_carlo_step!(particles, reference_cell_lookup, L, cell_length, interaction_radius, η, noise_lims, dt)
    # Update pre-update values
    update_pre_update_values!(particles)
    Threads.@threads for i in eachindex(particles)
        update_particle!(particles[i], particles, reference_cell_lookup, L, cell_length, interaction_radius, η, noise_lims, dt)
    end
end

# Run simulation
function run_simulation(initial_particles, L, cell_length, interaction_radius, η, noise_lims, dt, num_time_steps, output_filename, frameskip=1)

    particles = copy(initial_particles)
    reference_cell_lookup = [Int[] for _ in 1:total_cells] # Initialize empty lists for each cell
    println("t = 0")

    # Open files for particle data and order parameters
    open(output_filename, "w") do output_file
        open("order_parameters_L$(L)_T$(T)_N$(N)_eta$(η).csv", "w") do order_file
            # Write headers
            write(output_file, "time_step,particle_id,x,y,theta\n")
            write(order_file, "time_step,P_t,S_t\n")
            
            for t in 0:num_time_steps
                if t % frameskip == 0
                    for particle in particles
                        write(output_file, "$t, $(particle.id), $(particle.x), $(particle.y), $(particle.θ)\n")
                    end
                    # Calculate and write order parameters
                    P_t, S_t = OrderParameters(particles)
                    write(order_file, "$t, $P_t, $S_t\n")
                end
                reference_cell_lookup = build_cell_lookup(particles, reference_cell_lookup, L, cell_length) # Build the cell lookup
                monte_carlo_step!(particles, reference_cell_lookup, L, cell_length, interaction_radius, η, noise_lims, dt)
                println("t = $t")
            end
        end
    end

    println("Simulation completed.")
    return particles
end

# Initialize and run simulation
initial_particles = initialize_particles(N, L, v_0, cell_length)
output_filename = "particles_L$(L)_T$(T)_N$(N)_eta$(η).csv"
final_particles = run_simulation(initial_particles, L, cell_length, interaction_radius, η, noise_lims, dt, T, output_filename, frameskip)
