using JLD
using Plots
using PrettyTables


function table(name,certification)
    data = [0 0 0 0 0 0 0]
    for deg in 2:3
        dictionnary = load("logs/"*name*"_"*string(deg)*"_"*string(certification)*"withScip.jld")["data"]
        data = [data; Integer(deg) maximum(dictionnary["array_lb"]) dictionnary["lb"] dictionnary["ub"] dictionnary["time"] dictionnary["iter"] dictionnary["certif_time"]]
    end
    println("certified gap", 100*(minimum(data[2:end,4])-maximum(data[2:end,3]))/minimum(data[2:end,4]));
    println("estimated gap", 100*(minimum(data[2:end,4])-maximum(data[2:end,2]))/minimum(data[2:end,4]));
    return pretty_table(data[2:end,:]; header = ["Degree", "Estimated LB","Certified LB","""Value feas. control (CL_V)""","Solution time (in s)","Iterations number","Certification time (in s)"],backend = Val(:latex))
end

