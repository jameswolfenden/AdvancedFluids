using HDF5
using Plots

# Function to read h5 file and extract pressure data
function extract_pressure_data(h5_file::String)
    pressures_over_time = []
    time_steps = []

    # Open the h5 file
    h5open(h5_file, "r") do file
        # Assuming a structure where time steps are indexed (adjust as necessary)
        for timestep in keys(file)
            push!(time_steps, timestep)

            # Access each domain for the current time step
            domain_pressures = []
            for domain in keys(file[timestep])
                if domain == "time"
                    continue
                end
                # Read pressure data
                pressure = file[timestep][domain]["p"]
                # Check the dimensions of the dataset
                dims = size(pressure)
                println("Dimensions of pressure dataset: ", dims)
                # Extract values along the 1D portion (line_indices)
                line_pressures = pressure[:, :, 3]
                push!(domain_pressures, line_pressures)
            end
            push!(pressures_over_time, domain_pressures)
        end
    end

    return pressures_over_time, time_steps
end

# Visualization
function plot_pressures_over_time(pressures_over_time, time_steps, domain_labels::Vector{String})
    # plot the pressure data for the first domain using a 2d histogram
    # convert the data from an array of arrays of arrays to a 3d matrix

    a = size(pressures_over_time, 1) # number of time steps
    b = size(pressures_over_time[1], 1) # number of domains
    c = size(pressures_over_time[1][1], 1)
    f = zeros(a, b, c, c)
    for i in 1:a
        for j in 1:b
            f[i, j, :, :] .= pressures_over_time[i][j]
        end
    end
    #println("domain 0 data: ", f[:, 1, :,3])
    #println("domain 1 data: ", f[:, 2, :,3])
    #p1 = heatmap(f[:, 1, :,3], c=:viridis, title="Domain 0", xlabel="Distance", ylabel="Time")
    #p2 = heatmap(f[:, 2, :,3], c=:viridis, title="Domain 1", xlabel="Distance", ylabel="Time")
    #plot(p1, p2, layout=(1, 2), size=(1000, 500))
    println("time 0 data: ", f[1, 1, :, :])
    println("time 1 data: ", f[2, 1, :, :])
    println("time 2 data: ", f[3, 1, :, :])
    #println("time 3 data: ", f[4, 1, :, :])

    #p3 = heatmap(f[1,1,:,:], c=:viridis, title="D 0, T 0")
    #p4 = heatmap(f[2,1,:,:], c=:viridis, title="D 0, T 1")
    #p5 = heatmap(f[3,1,:,:], c=:viridis, title="D 0, T 2")
    #p6 = heatmap(f[4,1,:,:], c=:viridis, title="D 0, T 3")
    #plot(p3, p4, p5, p6, layout=(2, 2), size=(1000, 1000))

    p4 = heatmap(f[1,1,:,:], c=:viridis, title="D 0, T 0")
    p5 = heatmap(f[1,2,:,:], c=:viridis, title="D 1, T 0")
    #p6 = heatmap(f[1,3,:,:], c=:viridis, title="D 2, T 0")
    p7 = heatmap(f[2,1,:,:], c=:viridis, title="D 0, T 1")
    p8 = heatmap(f[2,2,:,:], c=:viridis, title="D 1, T 1")
    #p9 = heatmap(f[2,3,:,:], c=:viridis, title="D 2, T 1")
    p10 = heatmap(f[3,1,:,:], c=:viridis, title="D 0, T 2")
    p11 = heatmap(f[3,2,:,:], c=:viridis, title="D 1, T 2")
    #p12 = heatmap(f[3,3,:,:], c=:viridis, title="D 2, T 2")
    #plot(p4, p5, p6, p7, p8, p9, p10, p11, p12, layout=(3, 3), size=(1000, 1000))
    plot(p4, p5, p7, p8, p10, p11, layout=(3, 2), size=(1000, 1000))

end

# Define parameters
h5_file = "test.h5"  # Path to your h5 file
domain_labels = ["Domain 0", "Domain 1"]  # Adjust based on your data

# Extract data and plot
pressures, time_steps = extract_pressure_data(h5_file)
plot_pressures_over_time(pressures, time_steps, domain_labels)
