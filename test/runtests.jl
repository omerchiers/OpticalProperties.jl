using OpticalProperties
using Test

@testset "metals" begin
    @test permittivity(Al(), 1.0e15) ≈ -493.40133059541773 + 60.31696233264097im
    @test permittivity(Au(), 1.0e15) ≈ -172.92729201753394 + 20.04870903024803im
    @test permittivity(Au_latella(1.0, 0.0), 1.0e15) ≈ -186.16029145670757 + 9.956927505496843im
    @test permittivity(Cu, 1.0e15) ≈ -162.83122496251903 + 17.1567318225186im
    @test permittivity(Ti, 1.0e15) ≈ -13.517226294883221 + 1.0446596041797966im
    @test permittivity(TiW, 1.0e15) ≈ -42.785310734463266 + 43.785310734463266im
    @test permittivity(TiN, 1.0e15) ≈ -38.478487614080834 + 26.31899174272056im
    @test permittivity(TiW_v2, 1.0e15) ≈ -27.24858757062147 + 28.24858757062147im
    @test permittivity(W, 1.0e15) ≈ -30.00599475458973 + 9.36155863619333im
end

@testset "dielectrics" begin
    @test permittivity(SiO2, 1.0e15) ≈ 2.0719132476965396 + 8.63648521966404e-5im
    @test permittivity(Vacuum, 1.0e15) ≈ 1.0 + 0.0im
end

@testset "Semi-conductors" begin
    @test permittivity(AlGaAs(0.0, 300.0), 1.0e15) ≈ 9.899045845146038 + 0.07545599286968264im
    @test permittivity(cBN, 1.0e15) ≈ 4.364022312148895 + 9.925510975982143e-5im
    @test permittivity(Si, 1.0e15) ≈ 11.953741590306542 + 0.0im
    @test permittivity(Si_cst, 1.0e15) ≈ 11.7 + 0.0im
    @test permittivity(SiC, 1.0e15) ≈ 6.62441617687685 + 6.935638015943665e-5im
    @test permittivity(SiN, 1.0e15) ≈ 3.947544978277769 + 0.0im
end

@testset "effective medium models" begin
    bm = Bruggeman(Cst(3.0 + 0.5 * im), Cst(5.0 + 1.5 * im), 1.0)
    mg = MaxwellGarnett(Cst(3.0 + 0.5 * im), Cst(5.0 + 1.5 * im), 0.0)

    @test permittivity(bm, 1e15) ≈ permittivity(Cst(5.0 + 1.5 * im), 1e15)
    @test permittivity(mg, 1e15) ≈ permittivity(Cst(5.0 + 1.5 * im), 1e15)
end
