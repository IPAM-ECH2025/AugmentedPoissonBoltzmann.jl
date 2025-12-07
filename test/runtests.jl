using AugmentedPoissonBoltzmann
using AugmentedPoissonBoltzmann.SolverCore
using Test: @testset, @test
using LessUnitful
using ExplicitImports, Markdown

@testset "ExplicitImports" begin
    @test ExplicitImports.check_no_implicit_imports(AugmentedPoissonBoltzmann, skip = [Base, Core, Markdown, AugmentedPoissonBoltzmann]) === nothing
end

@testset "AugmentedPoissonBoltzmann" begin
    @testset "Units" begin
        @test AugmentedPoissonBoltzmann.Units.RT(0) == 0
    end
    @testset "Grids" begin
        grid = AugmentedPoissonBoltzmann.Grid.create_full_cell(
            AugmentedPoissonBoltzmann.Grid.GeometricGrid(10.0f0, 0, 1.0f0, 1.0f0, false),
        )
        @test AugmentedPoissonBoltzmann.Grid.is_equivalent(grid, 1, 11, 10, 3)

        grid = AugmentedPoissonBoltzmann.Grid.create_half_cell(
            AugmentedPoissonBoltzmann.Grid.GeometricGrid(10.0f0, 0, 1.0f0, 1.0f0, false),
        )
        @test AugmentedPoissonBoltzmann.Grid.is_equivalent(grid, 1, 11, 10, 2)
    end
    @testset "solvercore" begin
        let
            data = AugmentedPBData()
            set_molarity!(data, 0.01)
            data.χ = 78.49 - 1
            cdl0 = SolverCore.dlcap0(data)
            @test cdl0 ≈ 22.84669184882525ufac"μF/cm^2"
            # Check with Bard/Faulkner: the value must be $(22.8)μF/cm^2")
        end
        let
            data = AugmentedPBData()
            set_molarity!(data, 0.01)
            ddata = SolverCore.DerivedData(data)
            sumyz = 0.0
            sumyv = ddata.y0_E * data.v0
            sumy = ddata.y0_E
            for α in 1:data.N
                v = (1.0 + data.κ[α]) * data.v0
                sumyz += ddata.y_E[α] * data.z[α]
                sumyv += ddata.y_E[α] * v
                sumy += ddata.y_E[α]
            end
            @test sumy ≈ 1.0
        end

    end


    @testset "mpbsolve" begin
        X, C0, CP, CM = mpbpsolve(n = 5)

        X_ref = 0.0:2.5e-11:1.0e-10

        C0_ref = [53.08541695923597, 53.24118825379715, 53.292594052218185, 53.24118825379715, 53.085416959235964]
        CP_ref = [0.15635629457868624, 0.1257104894377271, 0.10070027035371844, 0.08036330567162285, 0.06387852730895277]
        CM_ref = [0.06387852730895317, 0.08036330567162334, 0.10070027035371906, 0.12571048943772797, 0.1563562945786872]

        @test X ≈ X_ref
        @test C0 ≈ C0_ref
        @test CP ≈ CP_ref
        @test CM ≈ CM_ref
    end
    @testset "icmpbsolve" begin
        X, C0, CP, CM = icmpbpsolve(n = 5)

        X_ref = 0.0:2.5e-11:1.0e-10

        C0_ref = [33.00334734467418, 33.60728611643063, 33.81408042246456, 33.607286116430636, 33.00334734467418]
        CP_ref = [1.4492043173155262, 1.2120889663429761, 0.9860872535243383, 0.7788850230724209, 0.5966731968050023]
        CM_ref = [0.5966731968050024, 0.7788850230724208, 0.9860872535243385, 1.2120889663429766, 1.4492043173155262]

        @test X ≈ X_ref
        @test C0 ≈ C0_ref
        @test CP ≈ CP_ref
        @test CM ≈ CM_ref
    end
end
