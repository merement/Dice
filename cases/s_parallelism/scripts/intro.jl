using DrWatson
@quickactivate "s_parallelism"

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())
"""
)

sim_output_dir = datadir("sims/");
data_input_dir = datadir("data_sets/");

modules_path = projectdir() * "/src"
push!(LOAD_PATH, modules_path)
