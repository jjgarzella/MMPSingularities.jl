
using Plots
plotly()

function plot_simplex()

    # Generate a lattice of integer points within a specified range (e.g., -5 to 5 for each axis)
    range = 0:3  # Define the range for the lattice
    lattice_points = []
    
    # Loop over all integer combinations in the range
    for i in range
        for j in range
            for k in range
                if i + j + k ≤ 3
                    push!(lattice_points, [i, j, k])
                end
            end
        end
    end
    
    # Extract x, y, and z coordinates for the lattice points
    lattice_x = [p[1] for p in lattice_points]
    lattice_y = [p[2] for p in lattice_points]
    lattice_z = [p[3] for p in lattice_points]
    
    # Add the lattice points to the plot
    p = scatter(lattice_x, lattice_y, lattice_z, marker=:o, label="Lattice Points", color=:red)

    # Define the vertices of the simplex
    vertices = [
        [0, 0, 0],   # Origin
        [3, 0, 0],   # (3, 0, 0)
        [0, 3, 0],   # (0, 3, 0)
        [0, 0, 3]    # (0, 0, 3)
    ]
    
    # Extract x, y, and z coordinates of the vertices
    x = [v[1] for v in vertices]
    y = [v[2] for v in vertices]
    z = [v[3] for v in vertices]
    
    p = plot3d!(p, x, y, z, marker=:o, label="Vertices", color=:blue)
    
    # Add lines connecting every vertex to every other vertex
    for i in 1:length(vertices)
        for j in i+1:length(vertices)
            plot3d!(p, [x[i], x[j]], [y[i], y[j]], [z[i], z[j]], color=:black, linewidth=10, label="")
        end
    end

    
    # Show the plot
    display(p)
end
