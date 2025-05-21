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
|Nx|$(@bind nx PlutoUI.Slider(10:200, default=25, show_value=true))|
|Ny|$(@bind ny PlutoUI.Slider(10:200, default=25, show_value=true))|
|Nt|$(@bind nt PlutoUI.NumberField(10:100:5000, default=500))|
|lx|$(@bind lx PlutoUI.Slider(10:10:1000, default=10.0, show_value=true))|
|ly|$(@bind ly PlutoUI.Slider(10:10:200, default=10.0, show_value=true))|
|$T_{res}$|$(@bind T_res PlutoUI.Slider(10:200, default=120, show_value=true))|
|$T_{inj}$|$(@bind T_inj PlutoUI.Slider(10:200, default=90, show_value=true))|
|Rayleigh|$(@bind Ra_ PlutoUI.NumberField(1:6, default=1e4))
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

# ╔═╡ 3be91d09-7b8b-4fd7-8719-6c6ce9e16a8f
md"""
## Permeability
Variations in permeability are modelled by simply applying a modifier to the x and y velocity using a permeability map as a reference.

"""

# ╔═╡ 59b99db1-7852-40eb-b02a-8d3b6f79fcca
@bind taa PlutoUI.Slider(1:nt, show_value=true)

# ╔═╡ 6f2fa035-c06f-409c-a43e-238b782aac3f
begin
	global function Variogram(nx::Int, ny::Int, range=35.)
		table = (; z=[1.,0.,1.]) # table containing values to fit
		coord = [
			(nx/4, ny/4), 
			(nx/2, 0.75*ny), 
			(0.75*nx, ny/2)
		] # coordinate of table values
		geotable = georef(table, coord) # georeferencing values
		grid = CartesianGrid(100, 100)
		model = Kriging(GaussianVariogram(range=35.))
		interp = geotable |> Interpolate(grid, model=model)
		print(length(interp.z))
		return reshape(interp.z, (100, 100))
	end

	heatmap(Variogram(100, 100))
end

# ╔═╡ c30bde2d-e0df-4fb3-b5b1-11455a8364f2
begin
	dx = lx / (nx - 1)
	dy = ly / (ny - 1)
		
	dt = 0.0001 # Time step size
	Pr = 7.01  # Prandtl number (for water) # 7.01
	Ra = 1e5 # Rayleigh number (controls convection strength)
		
	g = 9.81 # Acceleration due to gravity
	β = 1 # Thermal expansion coefficient (water @ 100°C)
	ν_water = 2.938e-7 # Water dynamic viscosity
	ν = sqrt(g * β * abs(T_res - T_inj) * ny^3 / Ra) / sqrt(Ra * Pr)
	# ν = ν_water
	α = ν_water / Pr 		# Thermal diffusivity ~ 0.02 
	
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
	x_offset, y_offset = 2, 1
	T[(nx÷2)-x_offset:(nx÷2)+x_offset, ny:(ny-(2*y_offset)), :] .= T_inj

	# permeability values 
	perm = Variogram(nx, ny)
	for t in 1:(nt - 1) # time step 
		for coords ∈ CartesianIndices((2:(nx - 1), 2:(ny - 1)))
			i = coords[1]
			j = coords[2]
			kᵢⱼ = perm[i, j]
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
	        duc_dt = -(uc[i,j,t]*duc_dx+vc[i,j,t]*duc_dy)+ν*(d2uc_dx2+d2uc_dy2)*kᵢⱼ
	        dvc_dt = -(uc[i,j,t]*dvc_dx+vc[i, j, t]*dvc_dy)+ν*(d2vc_dx2 + d2vc_dy2) +g*β*(T[i,j,t]-(T_inj+T_res)/2)*kᵢⱼ
	        dT_dt = -(uc[i,j,t]*dT_dx+vc[i,j,t]*dT_dy)+α*(d2T_dx2+d2T_dy2)
	
	        # Update values
	        uc[i, j, t+1] = uc[i, j, t] + dt * duc_dt
	        vc[i, j, t+1] = vc[i, j, t] + dt * dvc_dt
	        T[i, j, t+1]  = T[i, j, t] + dt * dT_dt
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
		if t<10
			T[(nx÷2)-x_offset:(nx÷2)+x_offset, 15:20, t+1] .= T_inj
		end
		avg = (T_res + T_inj) / 2
	    
		T[:, ny, t+1] .= avg # top boundary
		T[1, :, t+1]  .= avg # side boundaries
		T[nx, :, t+1] .= avg # side boundaries
		T[:, 1, t+1]  .= avg # bottom boundary
	    #T[1, :, t+1] .= (T_res + T_inj) / 2 # side boundary
	    #T[nx, :, t+1] .= (T_res + T_inj) / 2 # side boundary
	end
	temperature_field=T
	print()
end

# ╔═╡ 0a72c6aa-0602-44b2-85a3-d68401aa7c31
begin
	# plotting 
	c_range = (T_inj, T_res) # color range
	
	fig, ax, hm = heatmap(temperature_field[:, :, taa], colorrange=c_range)
	Colorbar(fig[:, end+1], colorrange=c_range)
	fig
end

# ╔═╡ Cell order:
# ╠═0a66c416-3484-11f0-0778-a1e104318d7b
# ╟─f1d283dc-9b8a-40a1-a947-f29588b98c2f
# ╟─7c88cee5-b895-42b5-9df4-94bcb4425574
# ╟─909d507c-972e-4ce3-a8c3-3965823b0e2b
# ╟─9534870e-2b83-4590-8e3b-030f63816c43
# ╠═214f2271-730f-4236-9761-a574295f12bf
# ╟─3be91d09-7b8b-4fd7-8719-6c6ce9e16a8f
# ╠═c30bde2d-e0df-4fb3-b5b1-11455a8364f2
# ╟─59b99db1-7852-40eb-b02a-8d3b6f79fcca
# ╟─0a72c6aa-0602-44b2-85a3-d68401aa7c31
# ╠═6f2fa035-c06f-409c-a43e-238b782aac3f
