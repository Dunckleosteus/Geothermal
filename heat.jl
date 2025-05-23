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
	using GLMakie, PlutoUI, LinearAlgebra, Statistics
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

# ╔═╡ 11ed96d4-047c-4886-94b6-51e24df24197
md"""
Start temperature: $(@bind start_temperature NumberField(0:200, default=120))

Constant temperature boundaries $(@bind constant_boundaries CheckBox(default=true))

Water injection temperature: $(@bind water_temperature_ NumberField(0:200, default=90))

Time steps $(@bind nt NumberField(1:5000, default=600))
"""

# ╔═╡ 9275a3c2-2583-4eff-8937-4a99df419a6f
@bind time_list_ PlutoUI.RangeSlider(1:1:nt, show_value=true)

# ╔═╡ 73d1885c-30e5-4399-a03a-fd9ef66be088
"""
# sdEllipse
Signed distance function to create an ellipse
# Inputs
- p -> point with x, y coordinates
- ab -> list of 2 itemps containing minor and major axis of the ellipse
"""
global function sdEllipse(p::Vector{A}, ab::Vector{T}) where {T <: Real, A <: Real}
    p = abs.(p)
    if p[1] > p[2]
        p = reverse(p)
        ab = reverse(ab)
    end

    l = ab[2]^2 - ab[1]^2
    m = ab[1] * p[1] / l
    m2 = m^2
    n = ab[2] * p[2] / l
    n2 = n^2
    c = (m2 + n2 - 1.0) / 3.0
    c3 = c^3
    q = c3 + m2 * n2 * 2.0
    d = c3 + m2 * n2
    g = m + m * n2
    co = 0.0

    if d < 0.0
        h = acos(q / c3) / 3.0
        s = cos(h)
        t = sin(h) * sqrt(3.0)
        rx = sqrt(-c * (s + t + 2.0) + m2)
        ry = sqrt(-c * (s - t + 2.0) + m2)
        co = (ry + sign(l) * rx + abs(g) / (rx * ry) - m) / 2.0
    else
        h = 2.0 * m * n * sqrt(d)
        s = sign(q + h) * abs(q + h)^(1.0 / 3.0)
        u = sign(q - h) * abs(q - h)^(1.0 / 3.0)
        rx = -s - u - c * 4.0 + 2.0 * m2
        ry = (s - u) * sqrt(3.0)
        rm = sqrt(rx^2 + ry^2)
        co = (ry / sqrt(rm - rx) + 2.0 * g / rm - m) / 2.0
    end

    r = ab .* [co, sqrt(1.0 - co^2)]
    return norm(r - p) * sign(p[2] - r[2])
end


# ╔═╡ da05c62f-f69f-4b03-965e-74dc2293cf8c
"""
# Bubble 
Given a 3D matrix (x, y, time), replace the values of cells under an ellipse to a given value over a list of times. 
## Inputs 
- Matrix : x, y and time array containing float values
- center : tuple representing the x, y position of center
- axis1  : axis 1 of ellipse
- axis2  : axis 2 of ellipse
- value  : the value to replace selected cells with
- time   : Vector of times at which the value should be changed. If left empty will change all values. 
## Returns
Modified matrix with same dimensions as input matrix
"""
global function Bubble(
		matrix::Array{Float64, 3}, 
		center::Tuple{A, A}, 
		axis1::A, 
		axis2::A,
		value::C;
		time::Union{Vector{T}, Nothing}=nothing
	) where {T <: Real, A<:Integer, C<:Real}
	time = isnothing(time) ? range(1, size(matrix)[end]) : time;
	for t ∈ time 
		for i in 1:size(matrix)[1]
		    for j in 1:size(matrix)[2]
			   	point = [i, j]
			   	dist = sdEllipse(point .- center, [axis1, axis2])
				matrix[i, j, t] = dist <= 0 ? value : matrix[i, j, t]
		    end
		end
	end
	return matrix
end

# ╔═╡ 0281916d-03bf-4302-b616-741c70d02061
"""
# Mask
Given a 3D matrix (x, y, time), generate a filter array with same dimensions as input matrix where : 
- False -> outside ellipse sdf
- True  -> under ellipse sdf
## Inputs 
- Matrix : x, y and time array containing float values-
- center : tuple representing the x, y position of center
- axis1 : axis 1 of ellipse
- axis2 : axis 2 of ellipse
- value : the value to replace selected cells with
- time : Vector of times at which the value should be changed. If left empty will change all values.
## Outputs 
- Boolean array of same dimension as input matrix.
"""
global function Mask(
		matrix::Array{Float64, 3}, 
		center::Tuple{A, A}, 
		axis1::A, 
		axis2::A,
		value::C;
		time::Union{Vector{T}, Nothing}=nothing,
		bounds::Bool=true
	) where {T <: Real, A<:Integer, C<:Real}
	result = Array{Bool, 3}(undef, size(matrix))
	fill!(result, false)
	time = isnothing(time) ? range(1,size(matrix)[end]) : time;
	for t ∈ time 
		for i in 1:size(matrix)[1]
		    for j in 1:size(matrix)[2]
			   	point = [i, j]
			   	dist = sdEllipse(point .- center, [axis1, axis2])
				result[i, j, t] = dist <= 0
		    end
		end
	end
	return result
end

# ╔═╡ dbe15ec3-5e7d-4d3b-b61c-bad08c6d97eb
begin 
	water_temperature = water_temperature_ == start_temperature ? start_temperature -1 : water_temperature_

	bubble_width=5
	bubble_height=3
	
	nx, ny = 30, 15	# array size
	bubble_center=(nx÷2, ny÷2)
	u = zeros(nx, ny, nt)       # create array
	u .= start_temperature      # set all values to initial temperature
	α = 1.5e-7# thermal diffusivity
	time_list = collect(time_list_)
	u = Bubble(u, (nx÷2, ny÷2), bubble_width, bubble_height, water_temperature; time=time_list)
	mask::Array{Bool, 3} = Mask(u, (nx÷2, ny÷2), bubble_width, bubble_height, water_temperature; time=time_list)
	dt = 5e7 # 6
	dx = 10
	dy = 10
	for coord in CartesianIndices(size(u))
		x = coord[1]
		y = coord[2]
		t = coord[3]
		
		u[[1, end], :, t] .= start_temperature # constant temperature bounds
		u[:, [1, end], t] .= start_temperature

			
		if t < size(u)[end]
			u[x, y, t+1] = mask[x, y, t] ? water_temperature : α * dt * (
				(u[x<nx ? x+1 : x, y, t]-2*u[x, y, t]+u[x==1 ? x : x-1, y, t])/dx^2+
				(u[x, y<ny ? y+1 : y, t]-2*u[x, y, t]+u[x, y==1 ? y : y-1, t])/dy^2
			) + u[x, y, t]
			
		end
	end
end

# ╔═╡ 0b2c0b34-ae4a-46ca-9124-4e76c1bda032
md"""
Iteration : $(@bind tempus PlutoUI.Slider(1:1:size(u)[end], show_value=true))
"""

# ╔═╡ 3679bcb7-0544-4a19-bb33-4d2f68390149
begin
	# plotting
	
	tempus_meus = round((tempus * dt) * 0.00000038026486; digits=2)
	
	fig = Figure(size=(800, 600))
	ax = Axis(
		fig[1, 1],
		title = "t = $tempus_meus months",
		xlabel = "X",
		ylabel = "Y"
	)
	c_range = (water_temperature, start_temperature) # color range
	
	hm = heatmap!(u[:, :, tempus], colorrange=c_range)
	Colorbar(fig[:, end+1], colorrange=c_range)
	fig
end

# ╔═╡ 7eab2a58-5e34-45e8-9adc-4cd84880c819
md"""
#### Animations
"""

# ╔═╡ 72565615-d29b-4f18-88e4-cadbd6f9f05c
begin
	local fig = Figure(size=(800, 600))
	local ax = Axis(
		fig[1, 1],
		title = "t = $tempus_meus months",
		xlabel = "X",
		ylabel = "Y"
	)
	local c_range = (water_temperature, start_temperature)
	local hm = heatmap!(u[:, :, tempus], colorrange=c_range)
	
	Colorbar(fig[:, end+1], colorrange=c_range)
	
	# animation parameters
	framerate = 30
	filepath  = "animation.mp4"
	iterator = 1:1:nt
	
	record(fig, filepath, iterator; framerate = framerate) do it
		hm[1] = u[:, :, it]
	    # hm = heatmap!(u[:, :, it], colorrange=c_range)
		ax.title = "iteration $it"
	end
end

# ╔═╡ aeda06d7-4917-4bb8-8d7e-d4bbcaa707b3
md"""
# Change in system
"""

# ╔═╡ ad140132-4da3-46a2-bdb4-88661bf0e1b8
global function ToMonths(tempus, dt)
	return round((tempus * dt) * 0.00000038026486; digits=2)
end

# ╔═╡ 7f4caa81-4d6e-463b-a4e2-24bf4fb3f091
begin
	f = Figure()
	
	ax1 = Axis(f[1, 1], yticklabelcolor = :blue, title = "Average temperature of system and it's derivative")
	ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right, xlabel = "time in months")
	hidespines!(ax2)
	#hidexdecorations!(ax2)
	
	temperaturae = [mean(u[:, :, i]) for i in 1:size(u)[end]]
	times = [ToMonths(i, dt) for i in 1:size(u)[end]]
	temp_line = lines!(ax1, times, temperaturae, color = :blue)
	deriv_line = lines!(ax2, times[begin:end-1], diff(temperaturae), color=:red)

	axislegend(ax1, [[temp_line], [deriv_line]], ["temperature", "derivative"], framevisible = false, valign=:center)
	
	f
end

# ╔═╡ Cell order:
# ╟─9205d372-2a74-11f0-0eb7-f13c28ba9f81
# ╟─da05c62f-f69f-4b03-965e-74dc2293cf8c
# ╟─582b7363-eb51-4a3c-a4e4-1b0622f57132
# ╟─0281916d-03bf-4302-b616-741c70d02061
# ╟─11ed96d4-047c-4886-94b6-51e24df24197
# ╟─9275a3c2-2583-4eff-8937-4a99df419a6f
# ╟─dbe15ec3-5e7d-4d3b-b61c-bad08c6d97eb
# ╟─73d1885c-30e5-4399-a03a-fd9ef66be088
# ╟─0b2c0b34-ae4a-46ca-9124-4e76c1bda032
# ╠═3679bcb7-0544-4a19-bb33-4d2f68390149
# ╟─7eab2a58-5e34-45e8-9adc-4cd84880c819
# ╠═72565615-d29b-4f18-88e4-cadbd6f9f05c
# ╟─aeda06d7-4917-4bb8-8d7e-d4bbcaa707b3
# ╟─ad140132-4da3-46a2-bdb4-88661bf0e1b8
# ╠═7f4caa81-4d6e-463b-a4e2-24bf4fb3f091
