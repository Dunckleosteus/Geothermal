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

# ╔═╡ 9205d372-2a74-11f0-0eb7-f13c28ba9f81
begin
	using Pkg; 
	Pkg.activate(".");
	Pkg.status();
	using GLMakie, PlutoUI
end

# ╔═╡ 582b7363-eb51-4a3c-a4e4-1b0622f57132
md"""
2D heat conduction equation

$$\frac{\partial u}{\partial t} = \alpha (\frac{\partial ^2u}{\partial x^2}+\frac{\partial ^2u}{\partial y^2})$$


Where: 
- $u(x,y,t)$ temperature at a point x, y and time t
- $\alpha$ is thermal diffusivity of the material
- $\frac{\partial u}{\partial t}$ is partial derivative of u with respect to time
- $\frac{\partial ^2u}{\partial x^2} \text{ and } \frac{\partial^2u}{\partial y^2}$ the second partial derivatives of u with respect to x and y respectively.

# Discretization
Using the finite elements method : 

$$\frac{df}{dx}\approx \frac{f(x + dx) - f(x)}{dx}$$

$$\frac{d^2f}{dx^2}\approx \frac{f(x+dx) - 2f(x) + f(x-dx)}{dx^2}$$

The thermal conductivity equation can be re-written as : 

$$u(t+dt)=\alpha \times dt \times (
u(x+dx)-2u(x)+u(x-dx) + u(y+dy)-2u(y)+u(y-dy)
) + u(t)$$

And translated into code : 
```
u[x, y, t+dt]= α * dt * (
	u[x+dx, y, t]-2*u[x, y, t]+u[x-dx, y, t] +
	u[x, y+dy, t]-2*u[x, y, t]+u[x, y+dy, t]
) + u[x, y, t]
```
"""

# ╔═╡ ff05d6e5-e0d8-4231-a186-a43d5c142988
begin 

end

# ╔═╡ dbe15ec3-5e7d-4d3b-b61c-bad08c6d97eb
begin 
	start_temperature = 100     # temperature underground
	water_temperature = 10
	nx, ny, nt = 100, 100, 100 	# array size
	u = zeros(nx, ny, nt)       # create array
	u .= start_temperature      # set all values to initial temperature
	α = 0.1 						# thermal diffusivity

	u[1:10, 1:10, begin:20] .= water_temperature

	for x ∈ 1:nx
		for y ∈ 1:ny
			for t ∈ 1:(nt-1)
				dt = 1
				u[x, y, t+1]= α * dt * (
					u[x<nx ? x+1 : x, y, t]-2*u[x, y, t]+u[x==1 ? x : x-1, y, t]+
					u[x, y<ny ? y+1 : y, t]-2*u[x, y, t]+u[x, y==1 ? y : y-1, t]
				) + u[x, y, t]
			end
		end
	end
end

# ╔═╡ 0b2c0b34-ae4a-46ca-9124-4e76c1bda032
@bind tempus PlutoUI.Slider(1:1:100)

# ╔═╡ 3679bcb7-0544-4a19-bb33-4d2f68390149
GLMakie.heatmap(u[:, :, tempus])

# ╔═╡ 930ca8ca-43bf-4105-92b9-4e905cffabb7
u[:, :, tempus]

# ╔═╡ Cell order:
# ╠═9205d372-2a74-11f0-0eb7-f13c28ba9f81
# ╟─582b7363-eb51-4a3c-a4e4-1b0622f57132
# ╠═ff05d6e5-e0d8-4231-a186-a43d5c142988
# ╠═dbe15ec3-5e7d-4d3b-b61c-bad08c6d97eb
# ╟─0b2c0b34-ae4a-46ca-9124-4e76c1bda032
# ╟─3679bcb7-0544-4a19-bb33-4d2f68390149
# ╠═930ca8ca-43bf-4105-92b9-4e905cffabb7
