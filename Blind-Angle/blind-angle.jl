using Random, Statistics

# --- Parameters ---
L = 20            # System size
T = 1000          # Total simulation time steps
ρ = 2.0           # Particle density
v_0 = 0.1         # Particle speed
cell_length = 1.0 # Cell size for neighbor list
interaction_radius = 1.0 # Interaction radius
dt = 1.0          # Time step
N = Int(ρ * L^2)  # Number of particles
num_cells = Int(div(L, cell_length)) # Number of cells in one direction
total_cells = num_cells^2            # Total number of cells
frameskip = 2    # Output data every frameskip steps
r = interaction_radius
ω = v_0/(r*2.03)  # Maximum angular velocity
η =  0.1*ω           # Noise strength
noise_bounds = (-0.5, 0.5) # Noise bounds
blind_angle = [deg2rad(10), deg2rad(350)] # Blind angle range in radians

# --- Define Particle struct ---
mutable struct Particle
    id::Int
    x::Float64
    y::Float64
    θ::Float64
    v::Float64
    cell_index::Tuple{Int,Int}  # (i, j) cell coordinates (1-indexed)
    x_p::Float64  # previous x position
    y_p::Float64  # previous y position
    θ_p::Float64  # previous orientation
    cell_index_p::Tuple{Int,Int}  # previous cell coordinates
end

# --- Cell and Neighbor Finding Functions ---
# Convert (x,y) coordinates to cell indices (1-indexed)
function cell_index(x, y, L, cell_length)
    i = mod(floor(Int, x/cell_length), num_cells) + 1
    j = mod(floor(Int, y/cell_length), num_cells) + 1
    return (i, j)
end

# Convert cell (i,j) to 1D index
function cell_to_index(cell::Tuple{Int,Int})
    i, j = cell
    return i + (j - 1) * num_cells
end

# Build cell lookup table
function build_cell_lookup(particles, cell_lookup)
    # Clear all cells
    for i in eachindex(cell_lookup)
        empty!(cell_lookup[i])
    end
    
    # Add particles to cells
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

# Check if an angle is within the blind spot
function is_in_blind_spot(θ_self, θ_relative, blind_range)
    # Shift the blind angle by the particle's orientation
    θ_start = mod2pi(blind_range[1] + θ_self)
    θ_end = mod2pi(blind_range[2] + θ_self)
    
    # Check if the relative angle is within the blind spot
    if θ_start <= θ_end
        # Blind angle doesn't wrap around 0
        return θ_start <= θ_relative <= θ_end
    else
        # Blind angle wraps around 0
        return θ_relative >= θ_start || θ_relative <= θ_end
    end
end

# Find neighbors for a particle
function find_neighbors(particle::Particle, particles, reference_cell_lookup, L, cell_length, interaction_radius, blind_angle)
    i, j = particle.cell_index
    neighbors = Int[]
    
    # Check neighboring cells
    for di in -1:1, dj in -1:1
        i_neighbor = mod1(i + di, num_cells)
        j_neighbor = mod1(j + dj, num_cells)
        idx = cell_to_index((i_neighbor, j_neighbor))
        
        for neighbor_id in reference_cell_lookup[idx]
            if neighbor_id == particle.id
                push!(neighbors, neighbor_id)  # Include self
                continue
            end
            
            neighbor = particles[neighbor_id]
            dist = periodic_distance(particle.x_p, particle.y_p, neighbor.x_p, neighbor.y_p, L)
            
            if dist < interaction_radius
                # Calculate relative angle to the neighbor
                dx = neighbor.x_p - particle.x_p
                dy = neighbor.y_p - particle.y_p
                
                # Adjust for periodic boundaries
                if dx > L/2
                    dx -= L
                elseif dx < -L/2
                    dx += L
                end
                
                if dy > L/2
                    dy -= L
                elseif dy < -L/2
                    dy += L
                end
                
                # Calculate the angle from particle to neighbor
                θ_relative = mod2pi(atan(dy, dx))
                
                # Check if neighbor is in blind spot
                if !is_in_blind_spot(particle.θ_p, θ_relative, blind_angle)
                    push!(neighbors, neighbor_id)
                end
            end
        end
    end
    
    return neighbors
end

# --- Order Parameter Calculation ---
function calculate_order_parameters(particles)
    N = length(particles)
    
    # Pre-calculate sin and cos values for each particle
    sin_theta = [sin(p.θ) for p in particles]
    cos_theta = [cos(p.θ) for p in particles]
    sin_2theta = [sin(2 * p.θ) for p in particles]
    cos_2theta = [cos(2 * p.θ) for p in particles]
    
    # Calculate polar order parameter
    sum_sin = sum(sin_theta)
    sum_cos = sum(cos_theta)
    P_t = sqrt(sum_sin^2 + sum_cos^2) / N
    
    # Calculate nematic order parameter
    sum_sin2theta = sum(sin_2theta)
    sum_cos2theta = sum(cos_2theta)
    S_t = sqrt(sum_sin2theta^2 + sum_cos2theta^2) / N
    
    # Calculate center of mass
    xcom = mean(p.x for p in particles)
    ycom = mean(p.y for p in particles)
    
    # Calculate average absolute value of the normalized angular momentum
    m_a = 0.0
    for i in 1:N
        p = particles[i]
        # Adjust for periodic boundaries
        dx = p.x - xcom
        dy = p.y - ycom
        
        # Apply periodic boundary conditions
        if dx > L/2
            dx -= L
        elseif dx < -L/2
            dx += L
        end
        
        if dy > L/2
            dy -= L
        elseif dy < -L/2
            dy += L
        end
        
        # Calculate distance from center of mass
        r_cm_i = sqrt(dx^2 + dy^2)
        
        # Skip particles at the center of mass (to avoid division by zero)
        if r_cm_i < 1e-10
            continue
        end
        
        # Use pre-calculated values
        u_x = cos_theta[i]
        u_y = sin_theta[i]
        
        # Cross product of position and unit velocity vector
        cross_product = abs(dx * u_y - dy * u_x)
        
        # Normalize by distance from center of mass
        normalized_angular_momentum = cross_product / r_cm_i
        
        # Add to sum
        m_a += normalized_angular_momentum
    end
    
    # Divide by number of particles
    m_a /= N
    
    return P_t, S_t, m_a
end

# --- Initialization ---
function initialize_particles(N, L, v_0, cell_length)
    particles = Vector{Particle}(undef, N)
    for i in 1:N
        x = rand() * L
        y = rand() * L
        θ = 2π * rand()
        cell = cell_index(x, y, L, cell_length)
        particles[i] = Particle(i, x, y, θ, v_0, cell, x, y, θ, cell)
    end
    return particles
end

# --- Update Functions ---
# Generate white noise
function white_noise(η, bounds)
    lower_limit, upper_limit = bounds
    return η * (rand() * (upper_limit - lower_limit) + lower_limit)
end

# Update a single particle
function update_particle!(particle::Particle, particles, reference_cell_lookup, L, cell_length, interaction_radius, ω, η, noise_bounds, dt, blind_angle)
    # Find neighbors
    neighbors = find_neighbors(particle, particles, reference_cell_lookup, L, cell_length, interaction_radius, blind_angle)
    
    # Calculate average direction of neighbors
    sum_sin = 0.0
    sum_cos = 0.0
    for neighbor_id in neighbors
        neighbor = particles[neighbor_id]
        sum_sin += sin(neighbor.θ_p)
        sum_cos += cos(neighbor.θ_p)
    end
    
    # Calculate average direction
    θ_avg = atan(sum_sin, sum_cos)
    θ_avg = mod2pi(θ_avg)  # Ensure angle is in [0, 2π)
    
    # Calculate angular difference
    δθ = θ_avg - particle.θ_p
    # Normalize to [-π, π]
    if δθ > π
        δθ -= 2π
    elseif δθ < -π
        δθ += 2π
    end
    
    # Apply maximum angular velocity constraint
    max_turn = ω * dt
    θ_new = particle.θ_p
    
    if abs(δθ) <= max_turn
        θ_new = θ_avg
    elseif δθ < 0
        θ_new = particle.θ_p - max_turn
    else
        θ_new = particle.θ_p + max_turn
    end
    
    # Add noise
    θ_new = mod2pi(θ_new + white_noise(η, noise_bounds))
    
    # Update position
    vx = particle.v * cos(θ_new)
    vy = particle.v * sin(θ_new)
    x_new = mod(particle.x + vx * dt, L)
    y_new = mod(particle.y + vy * dt, L)
    
    # Update cell index
    new_cell = cell_index(x_new, y_new, L, cell_length)
    
    # Update particle data
    particle.x, particle.y, particle.θ, particle.cell_index = x_new, y_new, θ_new, new_cell
end

# Update pre-update values for all particles
function update_pre_update_values!(particles)
    Threads.@threads for i in eachindex(particles)
        particle = particles[i]
        particle.x_p = particle.x
        particle.y_p = particle.y
        particle.θ_p = particle.θ
        particle.cell_index_p = particle.cell_index
    end
end

# Perform one Monte Carlo step
function monte_carlo_step!(particles, reference_cell_lookup, L, cell_length, interaction_radius, ω, η, noise_bounds, dt, blind_angle)
    # Update pre-update values
    update_pre_update_values!(particles)
    
    # Build cell lookup
    reference_cell_lookup = build_cell_lookup(particles, reference_cell_lookup)
    
    # Update particles in parallel
    Threads.@threads for i in eachindex(particles)
        update_particle!(particles[i], particles, reference_cell_lookup, L, cell_length, interaction_radius, ω, η, noise_bounds, dt, blind_angle)
    end
end

# --- Main Simulation Function ---
function run_simulation(L, T, ρ, v_0, cell_length, interaction_radius, ω, η, noise_bounds, dt, frameskip, blind_angle)
    N = Int(ρ * L^2)
    
    # Initialize particles
    particles = initialize_particles(N, L, v_0, cell_length)
    
    # Initialize cell lookup
    reference_cell_lookup = [Int[] for _ in 1:total_cells]
    
    # File names for output
    particles_filename = "modified_vicsek_particles_L$(L)_T$(T)_N$(N)_eta$(η)_omega$(round(ω, digits=2)).csv"
    order_filename = "modified_vicsek_order_L$(L)_T$(T)_N$(N)_eta$(η)_omega$(round(ω, digits=2)).csv"
    
    # Open files for output
    open(particles_filename, "w") do particle_file
        open(order_filename, "w") do order_file
            # Write headers
            write(particle_file, "time_step,particle_id,x,y,theta\n")
            # write(order_file, "time_step,P_t,S_t,avg_vx,avg_vy,angular_momentum\n")
            write(order_file, "time_step,P_t,S_t,m_a\n")
            
            # Initial state
            for particle in particles
                write(particle_file, "0,$(particle.id),$(particle.x),$(particle.y),$(particle.θ)\n")
            end
            
            # Calculate and write initial order parameters
            P_t, S_t, m_a = calculate_order_parameters(particles)
            write(order_file, "0,$P_t,$S_t,$m_a\n")            
            # Main simulation loop
            for t in 1:T
                # Update particles
                monte_carlo_step!(particles, reference_cell_lookup, L, cell_length, interaction_radius, ω, η, noise_bounds, dt, blind_angle)
                
                # Write data every frameskip steps
                if t % frameskip == 0
                    for particle in particles
                        write(particle_file, "$t,$(particle.id),$(particle.x),$(particle.y),$(particle.θ)\n")
                    end
                    
                    # Calculate and write order parameters
                    P_t, S_t, m_a = calculate_order_parameters(particles)
                    write(order_file, "$t,$P_t,$S_t,$m_a\n")
                    
                end
                
                # Print progress
                if t % 50 == 0
                    println("t = $t")
                end
            end
        end
    end
    
    println("Simulation completed. Data saved to $particles_filename and $order_filename")
    return particles
end

# Run the simulation
@time run_simulation(L, T, ρ, v_0, cell_length, interaction_radius, ω, η, noise_bounds, dt, frameskip, blind_angle)
