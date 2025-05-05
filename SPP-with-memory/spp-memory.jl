using Random, Plots, Statistics

#set seed for reproducibility
# Random.seed!(1234)

# --- Parameters ---
L = 64
T = 1200
ρ = 5
v_0 = 1/2
cell_length = 1
interaction_radius = 1
dt = 1
N = Int(ρ * L^2)
num_cells = Int(div(L, cell_length)) #number of cells in one direction
total_cells = num_cells * num_cells
ω_0 = -0.1
σ_noise = 0.003
α = 0.1
τ = 100
frameskip = 5 #print every 200 time steps

# --- Functions ---
# Define a custom mutable particle struct
mutable struct Particle
    id::Int
    x::Float64
    y::Float64
    θ::Float64
    ω::Float64
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
        ω = 0.0
        cell = cell_index(x, y, L, cell_length)
        x_p = x
        y_p = y
        θ_p = θ
        cell_p = cell
        #p is the pre-update quantity
        particles[i] = Particle(i, x, y, θ, ω, v_0, cell, x_p, y_p, θ_p, cell_p)
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
function white_noise(σ, μ=0.0)
    return μ + σ*randn()
end
# Update a single particle
function update_particle!(particle::Particle, particles, reference_cell_lookup, L, cell_length, interaction_radius, ω_0, σ_noise, α, τ, dt)

    # Find neighbors
    neighbors = find_neighbors(particle, particles, reference_cell_lookup, L, cell_length, interaction_radius)
    
    # Update angular velocity
    ω_new = particle.ω * (1 - dt/τ) + ω_0 * dt/τ + white_noise(σ_noise) * dt

    # Update angle using nematic alignment
    sum_sin = sum(sin(2*(particles[n].θ_p - particle.θ_p)) for n in neighbors)
    θ_new = particle.θ + particle.ω * dt + α * (sum_sin / max(length(neighbors), 1)) * dt

    # Update position
    v_x, v_y = velocity_vector(particle.v, θ_new)
    x_new = mod(particle.x + v_x*dt, L)
    y_new = mod(particle.y + v_y*dt, L)
    new_cell = cell_index(x_new, y_new, L, cell_length)

    # Update particle data
    particle.x, particle.y, particle.θ, particle.ω, particle.cell_index = x_new, y_new, θ_new, ω_new, new_cell
end
function update_pre_update_values!(particles)
    Threads.@threads for i in eachindex(particles)
        particle = particles[i]
        particle.x_p, particle.y_p, particle.θ_p, particle.cell_index_p = particle.x, particle.y, particle.θ, particle.cell_index
    end
end

# Parallelized Monte Carlo step
function monte_carlo_step!(particles, reference_cell_lookup, L, cell_length, interaction_radius, ω_0, σ_noise, α, τ, dt)
    update_pre_update_values!(particles)
    Threads.@threads for i in eachindex(particles)
        update_particle!(particles[i], particles, reference_cell_lookup, L, cell_length, interaction_radius, ω_0, σ_noise, α, τ, dt)
    end
end
# Run simulation
function run_simulation(
    initial_particles, L, cell_length, interaction_radius, ω_0, σ_noise, α, τ, dt,
    num_time_steps, output_filename, frameskip=1
)
    particles = copy(initial_particles)
    reference_cell_lookup = [Int[] for _ in 1:total_cells] # Initialize empty lists for each cell
    println("t = 0")
    # println("Particle 1: x = $(round(particles[1].x, digits=3)), x_p = $(round(particles[1].x_p, digits=3)), y = $(round(particles[1].y, digits=3)), y_p = $(round(particles[1].y_p, digits=3)), cell_index = $(particles[1].cell_index), cell_index_p = $(particles[1].cell_index_p)")
    # println("Initial Particles", particles)

    open(output_filename, "w") do output_file
        write(output_file, "time_step,particle_id,x,y,theta,omega\n")
        
        for t in 0:num_time_steps
            if t % frameskip == 0
                for particle in particles
                    write(output_file, "$t, $(particle.id), $(particle.x), $(particle.y), $(particle.θ), $(particle.ω)\n")
                end
            end
            if t % 200 == 0
                scatter([particle.x for particle in particles], [particle.y for particle in particles], 
                    title="Particle Positions at t=$t",
                    xlabel="x",
                    ylabel="y",
                    xlims=(0, L),
                    ylims=(0, L),
                    legend=false,
                    markersize=0.5,
                    color=:blue)
                savefig("particles_t$t.png")
            end
            reference_cell_lookup = build_cell_lookup(particles, reference_cell_lookup, L, cell_length) #build the cell lookup
            monte_carlo_step!(particles, reference_cell_lookup, L, cell_length, interaction_radius, ω_0, σ_noise, α, τ, dt)
	    println("t = $t")
            #print the first particle's x, x_p, y, y_p, cell_index, cell_index_p truncated to 3 decimal places
            # println("Particle 1: x = $(round(particles[1].x, digits=3)), x_p = $(round(particles[1].x_p, digits=3)), y = $(round(particles[1].y, digits=3)), y_p = $(round(particles[1].y_p, digits=3)), cell_index = $(particles[1].cell_index), cell_index_p = $(particles[1].cell_index_p)")
        end
    end

    println("Simulation completed.")
    return particles
end

# Initialize and run simulation
initial_particles = initialize_particles(N, L, v_0, cell_length)
output_filename = "L$(L)_N$(N)_T$(T)_ρ$(ρ)_w0$(ω_0)_s$(σ_noise)_a$(α)_t$(τ).csv"
final_particles = run_simulation(initial_particles, L, cell_length, interaction_radius, ω_0, σ_noise, α, τ, dt, T, output_filename, frameskip)
