### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 5dcbc267-cf1c-4b41-a581-0b2e4278887c
import Pkg; Pkg.add("CuArrays"); using CuArrays

# ╔═╡ 5e416548-6483-4823-929d-9611e67a18c6
using DifferentialEquations

# ╔═╡ 4abf82f1-d2ce-429f-b3aa-9c039ada2c50
using Plots; gr()

# ╔═╡ 6d5d93e7-3a20-46df-98f4-b13724fd467c
using PlutoUI

# ╔═╡ c3e6134d-f30a-445a-973b-1f4ef96dc777
md"""
# ODE Example

This is an example of solving a simple ode of exponential growth.  The ODE is:

$$u' = 0.98u$$
"""

# ╔═╡ b089a9fe-8de5-4d4a-84f5-3ea999368367
f(u,p,t) = 0.98u

# ╔═╡ f732f584-c91f-465c-8a30-f0d72513a185
u0 = 1.0

# ╔═╡ 1f291dc5-f09b-4b5b-9965-c1a609098a2d
tspan = (0.0,1.0)

# ╔═╡ 02d2da85-3dbb-4f18-a3cb-b4c064b54cdf
begin
	prob = ODEProblem(f, u0, tspan)
	sol = solve(prob)
end

# ╔═╡ 73632b20-612f-4f93-9183-9070e5137c1a
md"We have a simple solution to the ODE, now lets plot and look at it"

# ╔═╡ 35694b84-55fa-4f6e-aaef-97f4596a3057
begin
	plot(sol, linewidth=5,
		title="Solution to teh linear ODE with a thick line",
		xaxis="Time (t)", yaxis="u(t) (in μm)",
		label="My thick line ;)")
	plot!(sol.t, t->1.0*exp(0.98t), lw=3, ls=:dash, label="True sol")
end

# ╔═╡ 81d21420-0a5a-42d9-b7cf-7110b2adf5f6
[t+u for (u, t) in tuples(sol)]

# ╔═╡ accc3bce-21ad-490e-adab-e3f32a450e45
sol

# ╔═╡ 7b022ac8-36bf-4e0b-913a-929f5766ea57
sol(0.45)

# ╔═╡ 96d9a1ab-de9a-420f-b923-737723796d4b
md"# Systems of ODES

Now we're doing multiple ODEs at once. Like what we're gonna do for the oscillator

$$\frac{dx}{dt}=\sigma(y-x)$$

$$\frac{dy}{dt}=x(\rho - z)-y$$

$$\frac{dz}{dt}=xy-\beta z$$
"

# ╔═╡ d2df2642-8310-427a-960e-a347f1cf6d15
function lorentz!(du, u, p ,t)
	σ, ρ, β = p
	du[1] = σ*(u[2] - u[1])
	du[2] = u[1]*(ρ-u[3]) - u[2]
	du[3] = u[1]*u[2] - β*u[3]
end

# ╔═╡ 082dc4b3-c286-49e4-9649-6bba7f830420
u0₁ = [1,0,0]

# ╔═╡ d5bbc172-ce3f-4500-b396-edb4d8aea06d
p = (10, 28, 8/3)

# ╔═╡ 0c30464c-9b95-4b42-8bc5-be4e61754102
begin
	tspan₁ = (0.0, 100.0)
	prob₁ = ODEProblem(lorentz!, u0₁, tspan₁, p)
	sol₁ = solve(prob₁)
end

# ╔═╡ 0f55deef-c8c5-4287-bb07-a5cf9020912a
md"**an important feature of the ODEProblem solver is that it INTERPOLATES by default**. So while it displays points, it also has those points fully interpolated"

# ╔═╡ c408a6c3-a0cb-4420-a1b0-bdd8e59410b4
with_terminal() do
	println(sol₁.t[10], ",", sol₁[10])
	println(sol₁.t[20], ",", sol₁[20])
	println(sol₁.t[30], ",", sol₁[30])
	println(sol₁.t[40], ",", sol₁[40])
end

# ╔═╡ 3648a422-f88b-4f0c-a5b3-cd9ae79259e1
plot(sol₁, vars=(1,2,3), lw=1)

# ╔═╡ 17e193cc-b014-4413-b780-49fa74b6a3cd
md"Here we see the effects of NOT using the interpolation"

# ╔═╡ 58c45b88-7a23-414e-9954-5cc6b8c7d78f
plot(sol₁, vars=(1,2,3), lw=1, denseplot=false)

# ╔═╡ 047051d4-dc81-4863-b3f8-5b9359ec9d53
md"# System of ODEs, but make it GPU"

# ╔═╡ 638bf840-1d86-4953-a313-2c9f235913ab
function lorentz_2!(du, u, p ,t)
	σ, ρ, β = p
	du[1] = σ*(u[2] - u[1])
	du[2] = u[1]*(ρ-u[3]) - u[2]
	du[3] = u[1]*u[2] - β*u[3]
end

# ╔═╡ 5a23e15a-2bd0-4f6a-95b2-9290b7461f01


# ╔═╡ Cell order:
# ╟─c3e6134d-f30a-445a-973b-1f4ef96dc777
# ╠═5e416548-6483-4823-929d-9611e67a18c6
# ╠═b089a9fe-8de5-4d4a-84f5-3ea999368367
# ╠═f732f584-c91f-465c-8a30-f0d72513a185
# ╠═1f291dc5-f09b-4b5b-9965-c1a609098a2d
# ╠═02d2da85-3dbb-4f18-a3cb-b4c064b54cdf
# ╟─73632b20-612f-4f93-9183-9070e5137c1a
# ╠═4abf82f1-d2ce-429f-b3aa-9c039ada2c50
# ╠═35694b84-55fa-4f6e-aaef-97f4596a3057
# ╠═81d21420-0a5a-42d9-b7cf-7110b2adf5f6
# ╠═accc3bce-21ad-490e-adab-e3f32a450e45
# ╠═7b022ac8-36bf-4e0b-913a-929f5766ea57
# ╟─96d9a1ab-de9a-420f-b923-737723796d4b
# ╠═d2df2642-8310-427a-960e-a347f1cf6d15
# ╠═082dc4b3-c286-49e4-9649-6bba7f830420
# ╠═d5bbc172-ce3f-4500-b396-edb4d8aea06d
# ╠═0c30464c-9b95-4b42-8bc5-be4e61754102
# ╟─0f55deef-c8c5-4287-bb07-a5cf9020912a
# ╠═6d5d93e7-3a20-46df-98f4-b13724fd467c
# ╠═c408a6c3-a0cb-4420-a1b0-bdd8e59410b4
# ╠═3648a422-f88b-4f0c-a5b3-cd9ae79259e1
# ╟─17e193cc-b014-4413-b780-49fa74b6a3cd
# ╠═58c45b88-7a23-414e-9954-5cc6b8c7d78f
# ╟─047051d4-dc81-4863-b3f8-5b9359ec9d53
# ╠═5dcbc267-cf1c-4b41-a581-0b2e4278887c
# ╠═638bf840-1d86-4953-a313-2c9f235913ab
# ╠═5a23e15a-2bd0-4f6a-95b2-9290b7461f01
