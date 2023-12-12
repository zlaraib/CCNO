using CUDA

# Get the number of available CUDA devices
num_devices = CUDA.device_count()

println("Number of CUDA devices found: ", num_devices)

# Print information about each device
for i in 1:num_devices
    device = CUDA.device(i)
    println("Device $i:")
    println("  Name: ", CUDA.name(device))
    println("  Compute capability: ", CUDA.compute_capability(device))
    println("  Total memory: ", CUDA.total_memory(device))
    println("  Driver version: ", CUDA.driver_version(device))
    println("  CUDA version: ", CUDA.runtime_version(device))
end