include("Annealing_updated.jl")


function Annealing(filename::String,minc::Int64,maxc::Int64,cooling::Float64)
    
    comm, q = simulator.annealing(filename*".gml",minc,maxc,100.0, 1e-8, cooling,10.0,filename*"_"*string(minc)*"_"*string(maxc)*"_"*string(cooling))
    println(comm,q)
end


Annealing("lesmis",1,30, 0.98)