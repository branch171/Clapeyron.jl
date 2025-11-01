using Clapeyron

# a VLLLE mixture, HELD finds all phases

components = ["hexane","aniline","1,1,2-trichloroethane","water","nitrogen"]
model = SRK(components;mixing=MHV2Rule, activity=UNIFAC)

p = 1.0e5
T = 5+273.15

z=ones(length(model))/length(model)

verbose = false
HELDflash = Clapeyron.tp_flash_impl(model,p,T,z, HELDTPFlash(add_random_guess = true, verbose = verbose))
verbose = true
GFEMflash = Clapeyron.tp_flash_impl(model,p,T,z, GFEMTPFlash(flashin = HELDflash, verbose = verbose))

beta = GFEMflash.fractions
xp = GFEMflash.compositions
vp = GFEMflash.volumes
Gsol =  GFEMflash.data.g

verbose = false
if !verbose
    for ip in eachindex(beta)
        println("Phase beta[$(ip)[] = $(beta[ip])")
    end
    println("Phase mole fraction:")
    for ip in eachindex(beta)
        println("Phase x[$(ip)[] = $(xp[ip])")
    end
    println("Phase volumes:")
    for ip in eachindex(beta)
        println("Phase volume[$(ip)[] = $(vp[ip])")
    end
    println("Minimum Gibbs Energy = $(Gsol)")
end