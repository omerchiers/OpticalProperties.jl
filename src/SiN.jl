
struct Brendel{T} <: OptProp
    dl :: Model
    σ  :: T
end


struct SerieBrendel <: OptProp
   serie :: Vector{Brendel}
end

permittivity(material :: Brendel,w) = permittivity(material.dl, w)
susceptibility(material :: Model,w) = permittivity(material, w) - material.eps0

function susceptibility(material :: Brendel, w)
    @unpack eps0,w0,wp,gamma0 = material.dl
    @unpack σ = material
    integrand(x) = exp(-(x-w0)^2/(2.0*σ^2))*susceptibility(Model(eps0,wp,x,gamma0),w)/sqrt(2.0*pi)/σ
    val, err = quadgk(integrand,10*(σ-w0),10*(σ+w0))
    return val
end


function susceptibility(material :: SerieBrendel, w)
    χ = 0.0 + 0.0*im
    for i in material.serie
        χ += susceptibility(i, w)
    end
    return χ
end

permittivity(material :: SerieBrendel , w) = material.serie[1].dl.eps0 + susceptibility(material,w)

b1 = Brendel(Model(6.12, 1240*c0*2*pi, 826*c0*pi, 2*c0*pi), 283*c0*2*pi)
b2 = Brendel(Model(6.12,    1*c0*2*pi, 925*c0*pi, 4*c0*pi),  44*c0*2*pi)
b3 = Brendel(Model(6.12,   63*c0*2*pi,1063*c0*pi,36*c0*pi),3468*c0*2*pi)
b4 = Brendel(Model(6.12,  702*c0*2*pi,1185*c0*pi,57*c0*pi), 592*c0*2*pi)
b5 = Brendel(Model(6.12,    0*c0*2*pi,7931*c0*pi, 1*c0*pi),2987*c0*2*pi)

SN20_300 = SerieBrendel([b1 , b2, b3, b4 ,b5])
