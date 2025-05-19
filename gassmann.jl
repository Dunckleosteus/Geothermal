### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 45e865ce-30a3-11f0-3998-6f69a1127874
begin
	using Pkg; Pkg.activate("."); Pkg.status()
	using PlutoUI, GLMakie, Statistics
end

# ╔═╡ 47888aac-a650-4856-b1af-a811b4920a21
md"""
# Gassmann equation
$$\frac{K_{sat}}{K_0-K_{sat}} = \frac{K_{dry}}{K_0-K_{dry}} + \frac{K_{fluid}}{(K_0-K_{fluid})*\phi}$$

Which basically translates to : wet rock = dry rock + water

Need to isolate Kdry
"""

# ╔═╡ a22ae92b-78d2-4a98-afda-bea18c86b487
md"""
- ρ $(@bind ρ NumberField(0:0.1:5, default=2.33)) g/cc
- ρ fluid $(@bind ρwater NumberField(0:0.1:5, default=1.038)) g/cc
- vp fluid $(@bind vp_water NumberField(0:0.1:5, default=1.7)) km/s
- μ fluid $(@bind μwater NumberField(0:0.1:5, default=0)) g/cc
"""

# ╔═╡ 426ec1ee-20dd-4891-b20a-642ed9fc58ac
global function K(ρ, vp, μ)::Real
	return (ρ * vp^2) - ((4/3) * μ)
end

# ╔═╡ 19b56726-9994-4052-ac00-0cd662d2a279
global function Vp(K, μ, ρ)
	return sqrt(
		(K+(3/4)*μ) / ρ
	)
end

# ╔═╡ 70ec516a-06c5-44dd-abca-030e48486a1d
function Ksat(K_dry, K_0, K_fluid, phi)
	C = (K_dry / (K_0 - K_dry)) + (K_fluid / ((K_0 - K_fluid) * phi))
	K_sat = (C * K_0) / (1 + C)
	return K_sat
end

# ╔═╡ 9488cd80-477e-476d-857d-aa1cac4d1ec2
global function Vs(μ, ρ)
	return sqrt(
		μ/ρ
	)
end

# ╔═╡ 65f45fe5-8894-4de5-87cd-40333cc40e09
function Kdry(K_sat, K_0, K_fluid, phi)
	  B = (K_sat / (K_0 - K_sat)) - (K_fluid / ((K_0 - K_fluid) * phi))
	  K_dry = (B * K_0) / (1 + B)
	  return K_dry
end

# ╔═╡ b8f81a7e-320b-46e7-8360-38fd3eed22c0
begin 
	# for gas mass 0.168, vp 0.590
	
	# burial = 2000 # m
	ϕ = 0.2
	dt_compressional = 83.67 	# μs/ft°, 
	dt_shear = 148.17 			# μs/ft
	# pressure = 350 				# bars 
	# temperature = 95 			# C

	K0 = 38 # K mineral
	# convert dt values into vp / vs values
	vp = 304.8/dt_compressional
	vs = 304.8/dt_shear

	# calculating ksat
	μ = ρ * vs^2

	ksat = K(ρ, vp, μ)


	# calculating kfluid
	Kwater   = K(ρwater, vp_water, 0) # should be water

	kdry = Kdry(ksat, K0, Kwater, ϕ)
	# results
	print("μ = $μ \nksat = $ksat \nVp = $vp\nVs = $vs\nKwater = $Kwater \nKdry = $kdry")
end

# ╔═╡ 7621a087-0183-4d26-8c1d-224256622d1b
md"""
## Fluid subsitution
Now that we have Kdry, we calculate a new ksat for a rock saturated with gas with the following pvt : 
- ρ  = $(@bind ρ_new NumberField(0:0.1:5, default=1.168)) g/cc
- vp = $(@bind vp_new NumberField(0:0.1:5, default=0)) km/s
- μ  = $(@bind μ_new NumberField(0:0.1:5, default=0))

Whith gas or oil, the ksat value will be close to kdry ($(round(kdry; digits=2))) because the fluid is not dense
"""

# ╔═╡ 71775f1d-eb14-4420-b9d4-3b1a56eb6bd7
begin 
	manual = 11.78
	local kdry = isnothing(manual) ? kdry : manual # remplace calculated kdry
		
	K_fluid_new = K(ρ_new, vp_new, μ_new)
	ksat_new 	= Ksat(kdry, K0, K_fluid_new, ϕ)
	vp_new_ 		= Vp(ksat_new, μ_new, ρ_new)
	vs_new 		= Vs(μ_new, ρ_new)
	print("rho = $ρ_new\nvp = $vp_new\nμ = $μ_new\nkfluid = $K_fluid_new\nksat = $ksat_new\nvp = $vp_new\nvs = $vs_new\nkdry = $kdry")
end

# ╔═╡ Cell order:
# ╟─45e865ce-30a3-11f0-3998-6f69a1127874
# ╟─47888aac-a650-4856-b1af-a811b4920a21
# ╟─a22ae92b-78d2-4a98-afda-bea18c86b487
# ╠═b8f81a7e-320b-46e7-8360-38fd3eed22c0
# ╠═426ec1ee-20dd-4891-b20a-642ed9fc58ac
# ╠═19b56726-9994-4052-ac00-0cd662d2a279
# ╠═70ec516a-06c5-44dd-abca-030e48486a1d
# ╠═9488cd80-477e-476d-857d-aa1cac4d1ec2
# ╠═65f45fe5-8894-4de5-87cd-40333cc40e09
# ╟─7621a087-0183-4d26-8c1d-224256622d1b
# ╟─71775f1d-eb14-4420-b9d4-3b1a56eb6bd7
