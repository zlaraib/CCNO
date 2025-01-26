module CCNO

for file in readdir(@__DIR__)
    if endswith(file, ".jl") && file != "CCNO.jl"
        println("Included $file")
        include(file)
    end
end

end # module