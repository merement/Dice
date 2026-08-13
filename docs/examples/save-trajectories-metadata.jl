using Dice

traj = Hybrid[
    (Int8[1, -1, 1], Float64[0.12, 0.44, 0.91]),
    (Int8[1, -1, 1], Float64[0.15, 0.43, 0.90]),
]

time = [0.0, 0.1]

header = [
    Dice.comment("Simulation metadata"),
    Dice.format_record(),
    Dice.created("2026-01-26T18:44:00Z"),
    # or Dice.created() for inserting the current time
    Dice.title("Trajectory for test run 42"),
    Dice.node_indexing(1),
    Dice.node(1, "CTRL_1"),
    Dice.node(3, "pass_a"),
    Dice.record(:node, index = 3, name = "pass_a"),
    Dice.record(:solver; name = "my-relaxer", tolerance = 1e-8),
    Dice.metadata((git_commit = "abc123", seed = 42)),
]

Dice.save_trajectories(stdout, traj, time, header)
# trajectories are automatically numerated
# Dice.save_trajectories(stdout, traj, header)
