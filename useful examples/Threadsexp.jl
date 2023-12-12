# Import the Threads module
using Base.Threads

# Function to be executed by each thread
function myfunction(thread_id)
    println("Thread $thread_id started")
    # Your threaded code here
    for i in 1:5
        println("Thread $thread_id: iteration $i")
        sleep(1)  # Simulate some work
    end
    println("Thread $thread_id finished")
end

# Number of threads to create
num_threads = 3

# Create and start the threads
for i in 1:num_threads
    t = Threads.@spawn myfunction(i)  # Create a thread calling myfunction with argument i
    wait(t)  # Wait for the thread to finish (optional)
end

println("All threads have finished")
