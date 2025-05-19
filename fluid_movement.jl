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
	using GLMakie, PlutoUI;
end

# ╔═╡ e61c7792-6431-4438-83b0-a56ed1b772eb
md"""
## 2D dispersion equation

## Convection-Diffusion
A mixing of both terms is presented [here](https://louis-finegan.github.io/Convection-Diffusion-Pages/pages/Two_Dimensional_Models.html).
"""

# ╔═╡ 59b99db1-7852-40eb-b02a-8d3b6f79fcca
@bind taa PlutoUI.Slider(1:500)

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
$(@bind t_view PlutoUI.Slider(1:nt, show_value=true))
"""

# ╔═╡ 909d507c-972e-4ce3-a8c3-3965823b0e2b
begin
	ud = copy(u) # diffusion matrix
	for coord in CartesianIndices((nx-1, ny-1, nt-1))
		i = coord[1]
		j = coord[2]
		t = coord[3]

	    uxx = (ud[i+1, j, t] - 2*ud[i,j, t] + ud[i==1 ? 1 : i-1, j, t]) / dx
	    uyy = (ud[i, j+1, t] - 2*ud[i,j, t] + ud[i,j==1 ? 1 : j-1, t]) / dy
	    ud[i, j, t+1] = ud[i, j, t] + dt * D * (uxx + uyy)
	end
	heatmap(ud[:, :, t_view])
end

# ╔═╡ 9534870e-2b83-4590-8e3b-030f63816c43
md"""
## 2D convection equation
Convection only occurs in fluids and can be expressed with [this](https://drzgan.github.io/Python_CFD/9.%202D%20nonlinear%20convection.html) equation :

$$\frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} = 0$$

$$\frac{\partial v}{\partial t} + u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} = 0$$

Where are the u and v component of velocity (x velocity and y velocity)

### Discretization 

$$u_{i,j}^{n+1} = u_{i,j}^n - u_{i,j} \frac{\Delta t}{\Delta x} (u_{i,j}^n-u_{i-1,j}^n) - v_{i,j}^n \frac{\Delta t}{\Delta y} (u_{i,j}^n-u_{i,j-1}^n)$$ 

And v: 

$$v_{i,j}^{n+1} = v_{i,j}^n - u_{i,j} \frac{\Delta t}{\Delta x} (v_{i,j}^n-v_{i-1,j}^n) - v_{i,j}^n \frac{\Delta t}{\Delta y} (v_{i,j}^n-v_{i,j-1}^n)$$

### Code implementation
$(@bind t_step PlutoUI.Slider(1:1:nt, show_value=true))
"""

# ╔═╡ 99f22b76-b00d-4ef9-a11c-db774aa1f490
# ╠═╡ disabled = true
#=╠═╡
begin
	local uc = Array{Float64, 3}(undef, (nx, ny, nt)); fill!(uc, 0.0)
	local vc = Array{Float64, 3}(undef, (nx, ny, nt)); fill!(vc, 0.0)
	local ϕ  = Array{Float64, 3}(undef, (nx, ny, nt)); fill!(ϕ, 0.0)

	# boundary condi
	# initial conditions
	uc[9:11, 9:11, 1] .= 2
	vc[9:11, 9:11, 1] .= 0.2
	
	local dxdt = dx/dt # ratio of dx and dt
	local dydt = dy/dt # ratio of dy and dt
	for coord in CartesianIndices((nx-1, ny-1, nt-1))
		# unpacking indices
		i = coord[1]
		j = coord[2] 
		t = coord[3] # timestep

		# cell values
		local U = uc[coord] # value @ current cell
		local V = vc[coord] # value @ current cell
		
		# x velocity (u)
		local uᵢ = U * dxdt * (U - uc[max(i-1, 1), j, t])
		local uⱼ = V * dydt * (U - uc[i, max(j-1, 1), t])
		uc[i, j, t+1] = U - uᵢ - uⱼ 	# assigning value to next timestep

		# y velocity (v)
		local vᵢ=U*dxdt*(V-vc[max(i-1, 1), j, t])
		local vⱼ=V*dydt*(V-vc[i, max(j-1, 1), t])
		vc[i, j, t+1] = V - vᵢ - vⱼ   # assigning valeus to next timestep

		# Convection term
		ϕ[i, j, t+1] = ϕ[i, j, t] - vc[i, j, t+1] - uc[i, j, t+1]
	end
	heatmap(ϕ[:, :, t_step])
end
  ╠═╡ =#

# ╔═╡ ec938beb-8fc0-483f-b2b0-c02bd68a1c80
# ╠═╡ disabled = true
#=╠═╡
begin
	nx::Int, ny::Int, nt::Int = 20, 20, 20
	u = Array{Float64, 3}(undef, (nx, ny, nt))
	fill!(u, 0.0)
	u[5:10, 5:10, :] .= 10
	
	dx::Float64, dy::Float64, dt::Float64 = 10.0, 10.0, 20.0
	D::Float64 = 0.01 # TODO
	
end
  ╠═╡ =#

# ╔═╡ c30bde2d-e0df-4fb3-b5b1-11455a8364f2
begin

function simulate_convection(nx, ny, nt, dx, dy, dt, Pr, Ra, T_hot, T_cold)
    # Initialize arrays
	#
    uc = Array{Float64, 3}(undef, nx, ny, nt); fill!(uc, 0.0)
    vc = Array{Float64, 3}(undef, nx, ny, nt); fill!(vc, 0.0)
    T = Array{Float64, 3}(undef, nx, ny, nt); fill!(T, (T_hot + T_cold) / 2)
    p = Array{Float64, 3}(undef, nx, ny, nt); fill!(p, 0.0) # Pressure

    # Set initial conditions (example: hot bottom, cold top)
    T[15:20, 15:20, 1:10] .= T_hot

    # Parameters
    g = 9.81 # Acceleration due to gravity
    β = 1.0  # Thermal expansion coefficient (simplified)
	# kinematic viscosity
    ν = sqrt(g * β * abs(T_hot - T_cold) * ny^3 / Ra) / sqrt(Ra * Pr)
    α = ν / Pr # Thermal diffusivity

    # Time-stepping loop
    for t in 1:(nt - 1)
        for i in 2:(nx - 1)
            for j in 2:(ny - 1)
                # Spatial derivatives (central difference for now)
                duc_dx = (uc[i+1, j, t] - uc[i-1, j, t]) / (2 * dx)
                duc_dy = (uc[i, j+1, t] - uc[i, j-1, t]) / (2 * dy)
                dvc_dx = (vc[i+1, j, t] - vc[i-1, j, t]) / (2 * dx)
                dvc_dy = (vc[i, j+1, t] - vc[i, j-1, t]) / (2 * dy)
				
                dT_dx = (T[i+1, j, t] - T[i-1, j, t]) / (2 * dx)
                dT_dy = (T[i, j+1, t] - T[i, j-1, t]) / (2 * dy)

                d2uc_dx2 = (uc[i+1, j, t] - 2 * uc[i, j, t] + uc[i-1, j, t]) / dx^2
                d2uc_dy2 = (uc[i, j+1, t] - 2 * uc[i, j, t] + uc[i, j-1, t]) / dy^2
                d2vc_dx2 = (vc[i+1, j, t] - 2 * vc[i, j, t] + vc[i-1, j, t]) / dx^2
                d2vc_dy2 = (vc[i, j+1, t] - 2 * vc[i, j, t] + vc[i, j-1, t]) / dy^2
                d2T_dx2 = (T[i+1, j, t] - 2 * T[i, j, t] + T[i-1, j, t]) / dx^2
                d2T_dy2 = (T[i, j+1, t] - 2 * T[i, j, t] + T[i, j-1, t]) / dy^2

                # Momentum equations (x and y) with Boussinesq approximation
                duc_dt = - (uc[i, j, t] * duc_dx + vc[i, j, t] * duc_dy) + ν * (d2uc_dx2 + d2uc_dy2)
                dvc_dt = - (uc[i, j, t] * dvc_dx + vc[i, j, t] * dvc_dy) + ν * (d2vc_dx2 + d2vc_dy2) + g * β * (T[i, j, t] - (T_hot + T_cold) / 2)
                dT_dt = - (uc[i, j, t] * dT_dx + vc[i, j, t] * dT_dy) + α * (d2T_dx2 + d2T_dy2)

                # Update values
                uc[i, j, t+1] = uc[i, j, t] + dt * duc_dt
                vc[i, j, t+1] = vc[i, j, t] + dt * dvc_dt
                T[i, j, t+1] = T[i, j, t] + dt * dT_dt
            end
        end

        # Boundary conditions (simplified no-slip and fixed temperature)
        uc[:, 1, t+1] .= 0.0; uc[:, ny, t+1] .= 0.0; uc[1, :, t+1] .= 0.0; uc[nx, :, t+1] .= 0.0
        vc[:, 1, t+1] .= 0.0; vc[:, ny, t+1] .= 0.0; vc[1, :, t+1] .= 0.0; vc[nx, :, t+1] .= 0.0
        # T[:, 1, t+1] .= T_hot; T[:, ny, t+1] .= T_cold
        T[1, :, t+1] .= (T_hot + T_cold) / 2 # Example side boundary
        T[nx, :, t+1] .= (T_hot + T_cold) / 2 # Example side boundary
    end

    return T
end

# --- Simulation parameters ---
nx = 31  # Number of grid points in x
ny = 31  # Number of grid points in y
nt = 500 # Number of time steps
lx = 10.0  # Length in x
ly = 10.0  # Length in y
dx = lx / (nx - 1)
dy = ly / (ny - 1)
dt = 0.001 # Time step size
Pr = 0.71  # Prandtl number (for air)
Ra = 1e5   # Rayleigh number (controls convection strength)
T_hot = 120
T_cold = 90
# t_step = nt # Time step to visualize

# Run the simulation
temperature_field = simulate_convection(nx, ny, nt, dx, dy, dt, Pr, Ra, T_hot, T_cold)

# Visualize the temperature field at the last time step
heatmap(temperature_field[:, :, taa])
end

# ╔═╡ Cell order:
# ╟─0a66c416-3484-11f0-0778-a1e104318d7b
# ╠═ec938beb-8fc0-483f-b2b0-c02bd68a1c80
# ╟─f1d283dc-9b8a-40a1-a947-f29588b98c2f
# ╟─909d507c-972e-4ce3-a8c3-3965823b0e2b
# ╟─9534870e-2b83-4590-8e3b-030f63816c43
# ╠═99f22b76-b00d-4ef9-a11c-db774aa1f490
# ╟─e61c7792-6431-4438-83b0-a56ed1b772eb
# ╠═59b99db1-7852-40eb-b02a-8d3b6f79fcca
# ╠═c30bde2d-e0df-4fb3-b5b1-11455a8364f2
