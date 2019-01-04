using OpticalProperties
using Test

@testset "effective medium models" begin

    bm = Bruggeman(Cst(3.0 + 0.5*im),Cst(5.0 + 1.5*im),1.0)
    mg = MaxwellGarnett(Cst(3.0 + 0.5*im),Cst(5.0 + 1.5*im),0.0)

    @test permittivity(bm,1e15) == permittivity(Cst(5.0 + 1.5*im),1e15)
    @test permittivity(mg,1e15) == permittivity(Cst(5.0 + 1.5*im),1e15)

end
