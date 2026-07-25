using Dice
using Graphs
using SimpleWeightedGraphs
using Test

@testset "Dice.jl Test Suite" begin

    @testset "Distance Metrics" begin
        s1 = Int8[1, -1, 1]
        s2 = Int8[-1, -1, 1]
        @test hamming_distance(s1, s2) == 1
        @test_throws DimensionMismatch hamming_distance(s1, Int8[1, 1])

        v1 = [1.0, 2.0, 3.0]
        v2 = [1.0, 2.0, 4.0]
        @test euclid_distance(v1, v2) == 1.0
        @test_throws DimensionMismatch euclid_distance(v1, [1.0, 2.0])
    end

    @testset "Model Construction & Propagation" begin
        # 1. Unweighted (binary) graph
        g = Graphs.path_graph(4)
        model = Model(g, Dice.model_2_hybrid_coupling, 0.1)
        @test model.graph === g
        @test model.scale == 0.1
        @test model.silence == 3 # default silence

        state = get_random_hybrid(nv(g), 1.0)
        @test length(state[1]) == 4
        @test length(state[2]) == 4

        state_tuple = get_random_hybrid(nv(g), (-0.5, 0.5))
        @test length(state_tuple[1]) == 4
        @test length(state_tuple[2]) == 4
        @test all(x -> -0.5 <= x <= 0.5, state_tuple[2])

        # Run propagation
        final_state = propagate(model, 10, state)
        @test length(final_state[1]) == 4
        @test length(final_state[2]) == 4
        @test all(s -> s == 1 || s == -1, final_state[1])

        # 2. Weighted graph
        wg = SimpleWeightedGraph(4)
        add_edge!(wg, 1, 2, 0.5)
        add_edge!(wg, 2, 3, 1.5)
        wmodel = Model(wg, Dice.model_2_hybrid_coupling, 0.05)
        @test wmodel.graph === wg
        @test wmodel.scale == 0.05
    end

    @testset "Rounding Methods" begin
        g = Graphs.path_graph(3)
        # Test get_best_configuration
        V = [1.5, -0.5, 0.5]
        becu, beco = get_best_configuration(g, V)
        @test length(beco) == 3
        @test all(s -> s == 1 || s == -1, beco)
    end

    @testset "Statistical Methods" begin
        arr = [1, 2, 2, 3]
        vals, pvals = integer_distribution(arr)
        @test vals == [1, 2, 3]
        @test pvals == [0.25, 0.5, 0.25]

        cvals, cpvals = cumulative_integer_distribution(arr)
        @test cvals == [0, 1, 2, 4]
        @test cpvals == [0.0, 0.25, 0.75, 1.0]
    end

    @testset "File Operations" begin
        g = Graphs.cycle_graph(4)
        temp_file = "temp_test_graph.mtx"
        
        try
            save_mtx_graph(g, temp_file)
            @test isfile(temp_file)
            
            g_loaded = load_mtx_graph(temp_file)
            @test nv(g_loaded) == 4
            @test ne(g_loaded) == 4
        finally
            if isfile(temp_file)
                rm(temp_file)
            end
        end
    end

    @testset "Deprecated Backward Compatibility" begin
        # Verify legacy wrappers are present and throw warnings but work
        s1 = Int8[1, -1, 1]
        s2 = Int8[-1, -1, 1]
        @test_logs (:warn, "Function `HammingD` is deprecated, use `hamming_distance` instead") @test HammingD(s1, s2) == 1

        v1 = [1.0, 2.0, 3.0]
        v2 = [1.0, 2.0, 4.0]
        @test_logs (:warn, "Function `EuclidD` is deprecated, use `euclid_distance` instead") @test EuclidD(v1, v2) == 1.0

        g = Graphs.path_graph(3)
        temp_file = "temp_legacy_graph.mtx"
        try
            @test_logs (:warn, "Function `saveMTXGraph` is deprecated, use `save_mtx_graph` instead") saveMTXGraph(g, temp_file)
            @test_logs (:warn, "Function `loadMTXGraph` is deprecated, use `load_mtx_graph` instead") loadMTXGraph(temp_file)
        finally
            if isfile(temp_file)
                rm(temp_file)
            end
        end
    end

end
