include(joinpath(@__DIR__, "src", "CompetitionBallistics.jl"))
using Printf
import .CompetitionBallistics.SixDOF: mccoy_308_168_shot, mccoy_308_168_aero, solve_6dof,
       moments_of_inertia, initial_spin_rate, gyroscopic_stability_6dof
import .CompetitionBallistics.DragModels: interp_table
import .CompetitionBallistics.SixDOF: MCCOY_308_168_CMA0

function valider()
    sp = mccoy_308_168_shot()
    m = sp.mass_grains*0.0647989e-3; d = sp.caliber_in*0.0254
    v0 = sp.muzzle_vel_fps*0.3048; tw = sp.twist_in*0.0254
    A = pi*d^2/4; rho = 1.225; a_s = 340.3
    p0 = initial_spin_rate(v0, tw)
    CMa = interp_table(MCCOY_308_168_CMA0, v0/a_s)
    Sg  = gyroscopic_stability_6dof(sp.Ix, p0, rho, A, d, sp.Iy, v0, CMa)
    Sgf = let (ixf,iyf) = moments_of_inertia(m,d,sp.bullet_length_in*0.0254)
        gyroscopic_stability_6dof(ixf, p0, rho, A, d, iyf, v0, CMa) end
    println("=== 1. Facteur de stabilite gyroscopique a la bouche ===")
    @printf("   inerties mesurees   : Sg = %.3f     <- publie : 1.70\n", Sg)
    @printf("   estimateur geometrique : Sg = %.3f  (ecart %+.0f %%)\n", Sgf, 100*(Sgf/Sg-1))
    println("""
     Cet ecart est ASSUME, ne pas le rattraper en retouchant les defauts.
     Les defauts de BulletGeometry sont des medianes de populations mesurees ;
     celle-ci porte un bateau de 0,51 cal, le 2e percentile des balles match
     modernes (mediane 0,782), et de 13 deg la ou la mediane est 7,5. Un defaut
     qui la viserait serait faux partout ailleurs.
     L'estimateur se juge sur la ligne < contour publie > de validation_inerties.jl,
     qui lui donne la vraie cote : -0,2 % sur Ix et +2,4 % sur Iy.""")

    println("\n=== 2. Amplitude du mouvement de tangage-lacet ===")
    tr = solve_6dof(mccoy_308_168_shot(record_every=5))
    function env(a, b)
        v = [p.alpha_t for p in tr if a <= p.x/0.9144 <= b]
        isempty(v) ? NaN : rad2deg(maximum(v))
    end
    @printf("   premier maximum (0-30 yd) : %.2f deg    <- publie : 2.0 deg\n", env(0,30))
    @printf("   180-200 yd                : %.2f deg    <- figure 9.3 : ~1.75 deg\n", env(180,200))
    @printf("   580-600 yd                : %.2f deg    <- figure 9.4 : ~2.2 deg\n", env(580,600))
    @printf("   870-900 yd                : %.2f deg    <- figure 9.5\n", env(870,900))
    @printf("   maximum sur tout le parcours : %.2f deg <- publie : jamais > 5 deg\n",
            rad2deg(maximum(p.alpha_t for p in tr)))
    e = tr[end]
    @printf("\n   arrivee : x=%.1f m  chute=%.3f m  derive=%.3f m  v=%.1f m/s (%.0f fps)\n",
            e.x, e.y, e.z, e.v_total, e.v_total/0.3048)
end
valider()
