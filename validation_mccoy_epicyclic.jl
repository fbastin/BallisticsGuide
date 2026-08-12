include(joinpath(@__DIR__, "src", "CompetitionBallistics.jl"))
using Printf
import .CompetitionBallistics.SixDOF: mccoy_308_168_shot, solve_6dof, initial_spin_rate
import .CompetitionBallistics.DragModels: interp_table
import .CompetitionBallistics.SixDOF: MCCOY_308_168_CMA0

function beat()
    sp = mccoy_308_168_shot()
    m = sp.mass_grains*0.0647989e-3; d = sp.caliber_in*0.0254
    v0 = sp.muzzle_vel_fps*0.3048; tw = sp.twist_in*0.0254
    A = pi*d^2/4; rho = 1.225; a_s = 340.3
    p0 = initial_spin_rate(v0, tw); CMa = interp_table(MCCOY_308_168_CMA0, v0/a_s)
    # theorie linearisee, en radians par calibre puis par seconde
    P = (sp.Ix/sp.Iy) * (p0*d/v0)
    M = rho*A*d^3*CMa/(2*sp.Iy)
    disc = P^2 - 4M
    VD = v0/d
    println("=== theorie linearisee (McCoy ch. 10), inerties mesurees ===")
    @printf("   P = %.6f /cal   M = %.6e /cal^2   P^2-4M = %.3e\n", P, M, disc)
    @printf("   phi_F' = %.1f rad/s   phi_S' = %.1f rad/s\n",
            0.5*(P+sqrt(disc))*VD, 0.5*(P-sqrt(disc))*VD)
    battement = sqrt(disc)*VD
    @printf("   battement phi_F'-phi_S' = %.1f rad/s  -> periode %.3f ms\n",
            battement, 2pi/battement*1e3)
    K = 25.0/battement
    @printf("   amplitude pour une vitesse initiale de 25 rad/s : 2K = %.2f deg  (publie : 2.0)\n",
            2*rad2deg(K))

    println("\n=== solveur : intervalle entre maxima d'incidence (30 premieres ms) ===")
    tr = solve_6dof(mccoy_308_168_shot(record_every=2, target_range_yd=40.0))
    tops = Float64[]
    for i in 2:length(tr)-1
        if tr[i].alpha_t > tr[i-1].alpha_t && tr[i].alpha_t >= tr[i+1].alpha_t
            push!(tops, tr[i].t)
        end
    end
    if length(tops) >= 3
        dts = diff(tops)
        @printf("   %d maxima, intervalle median = %.3f ms  -> battement %.1f rad/s\n",
                length(tops), 1e3*sort(dts)[div(end,2)], 2pi/sort(dts)[div(end,2)])
        @printf("   rapport solveur/theorie = %.2f\n", (2pi/sort(dts)[div(end,2)])/battement)
    else
        println("   pas assez de maxima detectes")
    end
end
beat()
