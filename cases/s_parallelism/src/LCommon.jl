# project: logic
# module: LCommon

# Some common service functions for the logic project

module LCommon

using Dice

function save_trajectories(traj_collection::Vector{Vector{Hybrid}},
                           prefix::String,
                           suffix::String)
    # Save the trajectories in the Spin Reader format
    # traj_collection is a collection of trajectories collected after
    # agitations. We save each agitation in a separate file.
    for (ind, traj) in enumerate(traj_collection)
        out_file_name = prefix * "_$(ind)_" * suffix * ".dat"
        open(out_file_name, "w") do outf
            for (t, state) in enumerate(traj)
                print(outf, "$t")
                for (s, x) in zip(state[1], state[2])
                    print(outf, " $s $x")
                end
                println(outf, "")
            end
        end
        println("File $out_file_name is generated")
    end
    
end

# # is the data type defined in Dice
# traj_collection::Vector{Vector{Dice.Hybrid}} = []
# out_traj_name_prefix = "and-trajectories"
# regular_suffix = "regular"
# fail_suffix = "failed"
# # traj = Dice.trajectories_pinned(model, 200, state, pinned_list)
# traj = Dice.trajectories(model, 300, state)
# # push!(traj_collection, traj)
# @show state = traj[end]
# p1 = show_state(state)
# plot(p1, size=(200, 200))

end # module ends here
