using JLD
using Plots
using PrettyTables


function table(name,certification)
    data = [0 0 0 0 0 0 0]
    for deg in 2:8
        dictionnary = load("logs/"*name*"_"*string(deg)*"_"*string(certification)*".jld")["data"]
        data = [data; Integer(deg) maximum(dictionnary["array_lb"]) dictionnary["lb"] dictionnary["ub"] dictionnary["time"] dictionnary["iter"] -1 ];
        #### A terme mettre la ligne ci-dessous
        #row = [deg;maximum(dictionnary["array_lb"]);dictionnary["lb"];dictionnary["ub"];dictionnary["time_heur"]+dictionnary["time"];dictionnary["iter"];dictionnary["certif_time"] ]
        #data   = [ row]
            #data = vcat(data,row)
    end
    println("certified gap", 100*(minimum(data[2:end,4])-maximum(data[2:end,3]))/minimum(data[2:end,4]));
    println("estimated gap", 100*(minimum(data[2:end,4])-maximum(data[2:end,2]))/minimum(data[2:end,4]));
    return pretty_table(data; header = ["Degree", "Estimated LB","Certified LB","""Value feas. control (CL_V)""","Solution time (in s)","Iterations number","Certification time (in s)"],backend = Val(:latex))
end

