module Get_Lambda

export get_Lambda

using PyCall 
using NamedArrays
using CSV


path_to_GetLambda_dir = "home/john/constraints/Get_Lambda"
path_to_darkkrome_dir = "home/john/darkkrome"
density_of_nT_grid = 10

cd("..")
cd("..")
cd("..")
cd("..")
cd("..")
cd(path_to_GetLambda_dir) # Path to constraints dir
import FromFile: @from
@from "DK_mods.jl" using DK_modded_files
cd("..")
cd("..")
cd("..")
cd("..")

### Run DarkKROME w/o compiling ###
function run_DK_wo_comp(n, T, r_m, r_M, r_α, ξ, ϵ, z)
    cd("..")
    cd("..")
    cd("..")
    cd("..")
    cd(path_to_darkkrome_dir)
    r_rotE = ((r_α^2) * (r_m^2))/(r_M)   # Molecular Rotational Energy
    r_binE = ((r_α^2) * (r_m))           # Atomic Binding Energy
    generate_react(r_m, r_M, r_α, ξ)                         # Generate new react file
    generate_test(n,T,ξ,ϵ,z,r_m,r_M,r_α, r_rotE, r_binE)     # Generate new test file
    cd("build")
    Lambda = replace(readchomp(pipeline(`./test $T $n`)), "E" => "e")
    if Lambda == ""
        println("rstep value in the test file is not large enough!")    # If the rstep value is not large enough to reach t_lim
    end
    cd("..")  
    Lambda   # Cooling rate at this n and T
end
function get_Lambda_raw_data(n::Vector, T::Vector, r_m, r_M, r_α, ξ, ϵ, z)
    cd("..")
    cd("..")
    cd("..")
    cd("..")
    cd("..")
    cd(path_to_darkkrome_dir)
    r_rotE = ((r_α^2) * (r_m^2))/(r_M)   # Molecular Rotational Energy
    r_binE = ((r_α^2) * (r_m))           # Atomic Binding Energy
    generate_react(r_m, r_M, r_α, ξ)                         # Generate new react file
    generate_test(n,T,ξ,ϵ,z,r_m,r_M,r_α, r_rotE, r_binE)     # Generate new test file
    run(`./darkkrome --test darkEquil -t -C`)                                    # Compile DarkKROME
    cooling_table = parse.(Float64,run_DK_wo_comp.(n, T', r_m, r_M, r_α, ξ, ϵ, z))      # Read in the values of Lambda over n and T
    cd("..")
    cd("..")
    cd("..")
    cd("..")
    cd("..")
    cooling_table # Return the raw cooling table
end

function get_Lambda(r_m, r_M, r_α, ξ, ϵ, z)
    n=10 .^ range(0, stop=10, length=density_of_nT_grid)                 # Create n part of the n/T grid
    T=10 .^ range(0, stop=10, length=density_of_nT_grid)                 # Create n part of the n/T grid
    data = get_Lambda_raw_data(n, T, r_m, r_M, r_α, ξ, ϵ, z)             # Make raw cooling table
    nameT = string.(T)                                                   # Make n column for the named array
    namen = string.(n)                                                   # Make T row for the named array
    cooling_array = (NamedArray(data, (namen, nameT), ("n", "T")))       # Make the cooling table with n and T values labled
end
    


end