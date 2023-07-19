module QuantitativeMRI

using LsqFit
using Optim
using EPGsim

include("T2mapping/mese.jl")
include("T2mapping/mese_epg.jl")
include("T1mapping/mp2rage.jl")

end # module
