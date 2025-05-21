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

# ╔═╡ 0a66c416-3484-11f0-0778-a1e104318d7b
begin
	using Pkg;
	Pkg.activate(".");
	Pkg.status();
	using GLMakie, PlutoUI, GeoStats, Random;
end

# ╔═╡ 7c88cee5-b895-42b5-9df4-94bcb4425574
nt_::Int = 20

# ╔═╡ f1d283dc-9b8a-40a1-a947-f29588b98c2f
md"""
## 2D diffusion equation
Diffusion in the context of heat transfer refers to the spreading of heat due to random molecular motions, which is essentially the same as conduction.

The 2D diffusion equation can be expressed as [follows](https://scipython.com/book/chapter-7-matplotlib/examples/the-two-dimensional-diffusion-equation/) : 

$$\frac{\nabla U}{\nabla t} = D(\frac{\nabla ^2 U}{\nabla x^2}+\frac{\nabla ^2U}{\nabla y^2})$$

Where : 

|Parameter|Meaning|
|---|---|
|$D$|diffusion coefficient|
|$U$|?|

In code: 
```python
for i in range(1, nx-1):
    for j in range(1, ny-1):
        uxx = (u0[i+1,j] - 2*u0[i,j] + u0[i-1,j]) / dx2
        uyy = (u0[i,j+1] - 2*u0[i,j] + u0[i,j-1]) / dy2
        u[i,j] = u0[i,j] + dt * D * (uxx + uyy)
```
$(@bind t_view PlutoUI.Slider(1:nt_, show_value=true))
"""

# ╔═╡ 9534870e-2b83-4590-8e3b-030f63816c43
md"""
## 2D convection equation
Convection only occurs in fluids and can be expressed with [this](https://drzgan.github.io/Python_CFD/9.%202D%20nonlinear%20convection.html) equation :

X momentum equation : 

$$\frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} = -\frac{1}{\rho_0}\frac{\partial P}{\partial x}+ \nu (\frac{\partial^2u}{\partial x^2} + \frac{\partial^2u}{\partial y^2})$$

Y momentum equation :

$$\frac{\partial v}{\partial t} + u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} = -\frac{1}{\rho_0}\frac{\partial P}{\partial y}+ \nu (\frac{\partial^2v}{\partial x^2} + \frac{\partial^2v}{\partial y^2}) - g \beta (T-T_{ref})$$

[Source](https://scientiairanica.sharif.edu/article_4022_4c02163d61a0e092d1b0b0519dd23f4e.pdf)

Where : 

| | |
|---|---|
|u|x mometum|
|v|y momentum|
|P|Pressure in Pa|
|$g \beta (T-T_{ref})$|Buoyancy term, makes the fluid rise due to differences in density|
|$\nu (\frac{\partial^2v}{\partial x^2}$|Viscous advection terms|
|T|Temperature (in Kelvin)|
|t|time (seconds)|
|ν|kinematic viscosity|
|g|acceleration due to gravity $m.s²$|
|β|expansion coefficient|

And distances in meters.

## Advection diffusion for temperature

$$\frac{\partial T}{\partial t}+u\frac{\partial T}{\partial x} + v\frac{\partial T}{\partial y} = \alpha (\frac{\partial^2T}{\partial x^2}+\frac{\partial^2 T}{\partial y^2})$$
Where terms on the right is the thermal diffusion. 

[Source](https://www3.nd.edu/~powers/ame.60614/hwex.2004/hw4.pdf)
## Discretization
Uses finite differences, with forward for time and central for spatial derivatives. 

$$u_{i,j}^{n+1} = u_{i,j}^n + \Delta t \left[ - \left( u_{i,j}^n \frac{u_{i+1,j}^n - u_{i-1,j}^n}{2 \Delta x} + v_{i,j}^n \frac{u_{i,j+1}^n - u_{i,j-1}^n}{2 \Delta y} \right) + \nu \left( \frac{u_{i+1,j}^n - 2u_{i,j}^n + u_{i-1,j}^n}{(\Delta x)^2} + \frac{u_{i,j+1}^n - 2u_{i,j}^n + u_{i,j-1}^n}{(\Delta y)^2} \right) \right]$$

$v_{i,j}^{n+1} = v_{i,j}^n + \Delta t \left[ - \left( u_{i,j}^n \frac{v_{i+1,j}^n - v_{i-1,j}^n}{2 \Delta x} + v_{i,j}^n \frac{v_{i,j+1}^n - v_{i,j-1}^n}{2 \Delta y} \right) + \nu \left( \frac{v_{i+1,j}^n - 2v_{i,j}^n + v_{i-1,j}^n}{(\Delta x)^2} + \frac{v_{i,j+1}^n - 2v_{i,j}^n + v_{i,j-1}^n}{(\Delta y)^2} \right) + g \beta (T_{i,j}^n - T_{ref}) \right]$

$$T_{i,j}^{n+1} = T_{i,j}^n + \Delta t \left[ - \left( u_{i,j}^n \frac{T_{i+1,j}^n - T_{i-1,j}^n}{2 \Delta x} + v_{i,j}^n \frac{T_{i,j+1}^n - T_{i,j-1}^n}{2 \Delta y} \right) + \alpha \left( \frac{T_{i+1,j}^n - 2T_{i,j}^n + T_{i-1,j}^n}{(\Delta x)^2} + \frac{T_{i,j+1}^n - 2T_{i,j}^n + T_{i,j-1}^n}{(\Delta y)^2} \right) \right]$$
"""

# ╔═╡ 214f2271-730f-4236-9761-a574295f12bf
md"""
## Input parameters
|Parameter|Value|
|---|---|
|Nx|$(@bind nx PlutoUI.Slider(10:200, default=41, show_value=true))|
|Ny|$(@bind ny PlutoUI.Slider(10:200, default=41, show_value=true))|
|dt|$(@bind dt PlutoUI.NumberField(0.001:1, default=0.02))|
|Nt|$(@bind nt PlutoUI.NumberField(100:100:5000, default=1200))|
|lx|$(@bind lx PlutoUI.Slider(1:10:1000, default=1.0, show_value=true))|
|ly|$(@bind ly PlutoUI.Slider(1:10:200, default=1.0, show_value=true))|
|$T_{res}$|$(@bind T_res_ PlutoUI.Slider(10:200, default=120, show_value=true))|
|$T_{inj}$|$(@bind T_inj_ PlutoUI.Slider(10:200, default=90, show_value=true))|
|Current|$(@bind current PlutoUI.Slider(-1:0.1:1, default=0.1, show_value=true))|
|Use permeability Variogram|$(@bind use_vario CheckBox())|
"""

# ╔═╡ 909d507c-972e-4ce3-a8c3-3965823b0e2b
begin
	local nx::Int = 20
	local ny::Int = 20 
	
	local dx::Float64 = 10.0
	local dy::Float64 = 10.0
	local dt::Float64 = 20.0
	local D::Float64 = 0.01 # TODO
	
	ud = Array{Float64, 3}(undef, (nx, ny, nt))
	fill!(ud, 0.0)
	ud[5:10, 5:10, :] .= 10

	for coord in CartesianIndices((nx-1, ny-1, nt_-1))
		i = coord[1]
		j = coord[2]
		t = coord[3]

	    uxx = (ud[i+1, j, t] - 2*ud[i, j, t] + ud[max(1, i-1), j, t]) / dx
	    uyy = (ud[i, j+1, t] - 2*ud[i, j, t] + ud[i, max(1, j-1), t])/ dy
	    ud[i, j, t+1] = ud[i, j, t] + dt * D * (uxx + uyy)
	end
	heatmap(ud[:, :, t_view])
end

# ╔═╡ 59b99db1-7852-40eb-b02a-8d3b6f79fcca
@bind taa PlutoUI.Slider(2:nt, show_value=true)

# ╔═╡ 64dbeb7c-dfe3-4ea5-95da-6d7d336df9b2
global function CreateMask(radius::Real, center::Tuple, nx, ny)::BitMatrix
	dists  = Array{Float64, 2}(undef, (nx, ny))
	
	cart_dist(a::Tuple{Real,Real},b::Tuple{Real,Real}) = sqrt(
	    (a[1]-b[1])^2 + (a[2]-b[2])^2
	)
	 
	dists = [cart_dist(center, (x,y)) for x in 1:nx, y in 1:ny] .< radius
	
	return dists
end

# ╔═╡ 097b2ba3-c060-433d-bf10-a3d157a84c3d
md"""
## Permeability maps
Create matrices containing both horizontal and vertical permeabilities. These are applied as multiplipliers to horizontal and vertical permeability values. 
"""

# ╔═╡ 3be91d09-7b8b-4fd7-8719-6c6ce9e16a8f
md"""
## Permeability
Rayleigh number can be linked to permeability with Raleigh-Darcy number ([source](https://web.mit.edu/1.63/www/Lec-notes/chap6_porous/6-6Lapwood.pdf)) : 

$$Ra = \frac{\rho \beta \Delta T klg}{n\alpha}$$
With : 
- k the permeability of the porous medium ✓
- ρ density of the fluid ✓
- ΔT difference in temperature between T₀ and T₀+ΔT ≈ Tres-Tinj ✓
- β expansion coefficient ✓
- l characteristic length ≈ height of reservoir in meters ✓
- ν kinematic viscosity ✓
- α thermal diffusivity ✓

"""

# ╔═╡ 6300158b-af62-4863-9944-77fc4569df4f
global function Variogram(nx::Int, ny::Int, range=35, n_rand=10)
	table = (; z=rand(Float64, n_rand)) # table containing values to fit
	ra = [tuple(x[1]*nx, x[2]*ny) for x in zip(rand(Float64, n_rand), rand(Float64, n_rand))]
	coord = ra# coordinate of table values
	geotable = georef(table, coord) # georeferencing values
	grid = CartesianGrid(nx, ny)
	model = Kriging(GaussianVariogram(range=range))
	interp = geotable |> Interpolate(grid, model=model)
	print(length(interp.z))
	return reshape(interp.z, (nx, ny))
end

# ╔═╡ f32aacf3-650a-4f54-8ca0-d26a91e921d3
md"""
### Vertical permeability profile
"""

# ╔═╡ a709ce76-97f3-45aa-9bae-a32e3dedd8e8
begin 
	x = 1:nx
	fperm(x, ny, λ=10, ▽=3) = (sin.(x ./ (ny/λ)) .* 0.1  .+ x / (ny*2))
	lines(x, fperm.(x, ny))
end

# ╔═╡ cc561259-1ad6-46a6-a323-c0d12d859012
begin
	permₕ = Array{Float64, 2}(undef, (nx, ny)); fill!(permₕ, 1.0)
	permₕ = [fperm(y, ny) for y in 1:nx, y in 1:ny]
	# dists = [cart_dist(center, (x, y)) < radius for x in 1:nx, y in 1:ny]
	if use_vario 
		permₕ = permₕ .* Variogram(nx, ny, 20.)
	end
	permₕ = permₕ ./10
	permᵥ = permₕ ./5
	surface(permᵥ)
end

# ╔═╡ c30bde2d-e0df-4fb3-b5b1-11455a8364f2
begin
	T_inj = T_inj_ + 273
	T_res = T_res_ + 273
	dx = lx / (nx - 1)
	dy = ly / (ny - 1)
	Pr = 7.01  # Prandtl number (for water)
	Ra = 1e5 # Rayleigh number (controls convection strength)
		
	g = 9.81 # Acceleration due to gravity
	β = 0.75e-3  # Thermal expansion coefficient (water @ 100°C)
	# TODO: change to not use Ra
	ν_water = 2.938e-7
	ν = sqrt(g * β * abs(T_res - T_inj) * ny^3 / Ra) / sqrt(Ra * Pr)
	α = ν/ Pr 		# Thermal diffusivity ~ 0.02 
	
	ρ_water = 1000 	# kg.m³
	# array of x velocities
	uc = Array{Float64, 3}(undef, nx, ny, nt)
	fill!(uc, 0.0)
	# array of y velocities
	vc = Array{Float64, 3}(undef, nx, ny, nt)
	fill!(vc, 0.0)
	# temperature array
	T = Array{Float64, 3}(undef, nx, ny, nt)
	fill!(T, (T_inj + T_res) / 2)
	# pressure array
	p = Array{Float64, 3}(undef, nx, ny, nt)
	fill!(p, 0.0)
	ΔT = abs(T_inj-T_res)
	# Initial conditions
	x_offset, y_offset = 10, 2
	# creating bubble shape 
	radius = 5
	center = (nx÷2, (ny÷4)*3)
	mask = CreateMask(radius, center, nx, ny)
	T[:, :, 1] .= ifelse.(mask, T_inj, T[:, :, 1])
	# Permeability 
	
	for t in 1:(nt - 1) # time step 
		for coords ∈ CartesianIndices((2:(nx - 1), 2:(ny - 1)))
			i = coords[1]
			j = coords[2]
			kₕ = permₕ[i, j]
			kᵥ = permᵥ[i, j]
			
			# Ra_Da = DarcyRayleigh(ρ_water, β, ΔT, k, ly, ν, α)
			# α = alpha(β, ΔT, Ra_Da, Pr, g)
	    	# Spatial derivatives (central differences)
	        duc_dx = (uc[i+1, j, t] - uc[i-1, j, t]) / (2 * dx)
	        duc_dy = (uc[i, j+1, t] - uc[i, j-1, t]) / (2 * dy)
	        dvc_dx = (vc[i+1, j, t] - vc[i-1, j, t]) / (2 * dx)
	        dvc_dy = (vc[i, j+1, t] - vc[i, j-1, t]) / (2 * dy)

	        dT_dx = (T[i+1, j, t] - T[i-1, j, t]) / (2 * dx) # ΔT/Δy
	        dT_dy = (T[i, j+1, t] - T[i, j-1, t]) / (2 * dy) # ΔT/Δy
	
	    	d2uc_dx2 = (uc[i+1, j, t] - 2 * uc[i, j, t] + uc[i-1, j, t]) / dx^2
	        d2uc_dy2 = (uc[i, j+1, t] - 2 * uc[i, j, t] + uc[i, j-1, t]) / dy^2
	        d2vc_dx2 = (vc[i+1, j, t] - 2 * vc[i, j, t] + vc[i-1, j, t]) / dx^2
	        d2vc_dy2 = (vc[i, j+1, t] - 2 * vc[i, j, t] + vc[i, j-1, t]) / dy^2
	        d2T_dx2 = (T[i+1, j, t] - 2 * T[i, j, t] + T[i-1, j, t]) / dx^2
	        d2T_dy2 = (T[i, j+1, t] - 2 * T[i, j, t] + T[i, j-1, t]) / dy^2
	
	        # Momentum equations (x and y) with Boussinesq approximation
duc_dt = -(uc[i,j,t]*duc_dx+vc[i,j,t]*duc_dy)+ν*(d2uc_dx2+d2uc_dy2)+current*vc[i,j,t]*kₕ
dvc_dt = -(uc[i,j,t]*dvc_dx+vc[i, j, t]*dvc_dy)+ν*(d2vc_dx2 + d2vc_dy2) +g*β*(T[i,j,t]-(T_inj+T_res)/2)*kᵥ
	        dT_dt = -(uc[i,j,t]*dT_dx+vc[i,j,t]*dT_dy)+α*(d2T_dx2+d2T_dy2)
	
	        # Update values
	        uc[i, j, t+1] = uc[i, j, t] + dt * duc_dt
	        vc[i, j, t+1] = vc[i, j, t] + dt * dvc_dt
	        T[i, j, t+1] = T[i, j, t] + dt * dT_dt
	    	# end
	    end
	
	    # Boundary conditions (simplified no-slip and fixed temperature)
		# UC boundaries (no velocity on edges)
	    uc[:, 1, t+1].= 0.0
		uc[:, ny, t+1].= 0.0
		uc[1, :, t+1].= 0.0
		uc[nx, :, t+1] .= 0.0
		# VC boundaries (no velocity on edges)
	    vc[:, 1, t+1] .= 0.0
		vc[:, ny, t+1] .= 0.0
		vc[1, :, t+1] .= 0.0
		vc[nx, :, t+1] .= 0.0
		# Temperature boundary conditions
		# T[15:20, 15:20, t+1] .= T_inj
		T[:, :, t+1] .= ifelse.(mask, T_inj, T[:, :, t+1])
		
		avg = (T_res + T_inj) / 2
	    
		T[:, ny, t+1] .= avg # top boundary
		T[1, :, t+1]  .= avg # side boundaries
		T[nx, :, t+1] .= avg # side boundaries
		T[:, 1, t+1]  .= T_res # bottom boundary
	    #T[1, :, t+1] .= (T_res + T_inj) / 2 # side boundary
	    #T[nx, :, t+1] .= (T_res + T_inj) / 2 # side boundary
	end
	temperature_field=T
	print()
end

# ╔═╡ 0a72c6aa-0602-44b2-85a3-d68401aa7c31
begin
	fig = Figure(size = (800, 600))
	ax = Axis(fig[1,1], title = "Temperature Field at iteration $taa")
	hm = contourf!(ax, temperature_field[:, :, taa])
	Colorbar(fig[1,2], hm)
	fig
end

# ╔═╡ da7216cd-865f-47ee-8e7b-9086cbaf8fb5
permₕ

# ╔═╡ Cell order:
# ╠═0a66c416-3484-11f0-0778-a1e104318d7b
# ╟─f1d283dc-9b8a-40a1-a947-f29588b98c2f
# ╟─7c88cee5-b895-42b5-9df4-94bcb4425574
# ╟─909d507c-972e-4ce3-a8c3-3965823b0e2b
# ╟─9534870e-2b83-4590-8e3b-030f63816c43
# ╠═214f2271-730f-4236-9761-a574295f12bf
# ╠═c30bde2d-e0df-4fb3-b5b1-11455a8364f2
# ╠═59b99db1-7852-40eb-b02a-8d3b6f79fcca
# ╠═0a72c6aa-0602-44b2-85a3-d68401aa7c31
# ╟─64dbeb7c-dfe3-4ea5-95da-6d7d336df9b2
# ╟─097b2ba3-c060-433d-bf10-a3d157a84c3d
# ╟─3be91d09-7b8b-4fd7-8719-6c6ce9e16a8f
# ╠═6300158b-af62-4863-9944-77fc4569df4f
# ╠═cc561259-1ad6-46a6-a323-c0d12d859012
# ╠═da7216cd-865f-47ee-8e7b-9086cbaf8fb5
# ╟─f32aacf3-650a-4f54-8ca0-d26a91e921d3
# ╠═a709ce76-97f3-45aa-9bae-a32e3dedd8e8
