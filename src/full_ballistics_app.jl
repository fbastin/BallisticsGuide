
using GLMakie
using Observables
using JSON3
using LinearAlgebra
using DelimitedFiles
using Interpolations

#=========================================
# Constants
#=========================================
g = 9.81

#=========================================
# Doppler Drag
#=========================================
mach = [0.3,0.5,0.9,1.0,1.2,1.5,2.0,2.5,3.0]
cd   = [0.115,0.120,0.140,0.250,0.300,0.270,0.240,0.230,0.220]
cd_spline = CubicSplineInterpolation(mach, cd)

function Cd(M)
    return cd_spline(M)
end

#=========================================
# Solver
#=========================================
function simulate(v0, angle; dt=0.001)

    vx, vy = v0*cos(angle), v0*sin(angle)
    x, y = 0.0, 0.0

    xs, ys, vs = Float64[], Float64[], Float64[]

    while y >= 0
        push!(xs, x)
        push!(ys, y)

        v = sqrt(vx^2 + vy^2)
        push!(vs, v)

        M = v/340
        drag = Cd(M)*v^2*0.001

        vx -= drag*(vx/v)*dt
        vy -= (g + drag*(vy/v))*dt

        x += vx*dt
        y += vy*dt
    end

    return xs, ys, vs
end

#=========================================
# DOPE
#=========================================
function generate_dope(xs, ys, vs)
    table = []
    for d in 100:100:1000
        idx = argmin(abs.(xs .- d))
        drop = ys[idx]
        vel  = vs[idx]
        mil = (drop/d)*1000
        moa = mil*3.43775
        push!(table, (d, drop, mil, moa, vel))
    end
    return table
end

#=========================================
# GUI
#=========================================
v0 = Observable(800.0)
angle = Observable(2.0)

traj = Observable((Float64[],Float64[],Float64[]))
dope = Observable([])

function update!()
    xs, ys, vs = simulate(v0[], deg2rad(angle[]))
    traj[] = (xs, ys, vs)
    dope[] = generate_dope(xs, ys, vs)
end

onany(v0, angle) do _
    update!()
end

update!()

fig = Figure(resolution=(1200,700))
tabs = TabLayout(fig[1,1])

# Trajectory
tab1 = tabs[1] = Tab("Trajectory")
ax = Axis(tab1[1,1])
lines!(ax, lift(t->t[1], traj), lift(t->t[2], traj))

# DOPE
tab2 = tabs[2] = Tab("DOPE")
label = Label(tab2[1,1], "")
on(dope) do t
    label.text[] = join([string(r) for r in t], "\n")
end

# Export
function export_csv()
    writedlm("dope.csv", dope[], ',')
end

Button(tab2[2,1], label="Export CSV") do _
    export_csv()
end

fig
