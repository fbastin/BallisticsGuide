# Point 3 : l'ecart sur les inerties, propage sur la DERIVE INTEGREE et non
# plus seulement sur l'angle de repos a la bouche.
#
# Chaine physique (McCoy ch. 10 et 12) :
#   angle de repos      aR = 2 Ix p g cos(phi) / (rho S d V^3 C_Ma)
#   force laterale      Fz = (rho/2) V^2 S C_La aR
#   -> acceleration     az = g (C_La/C_Ma) (Ix p)/(m d V)
# Le rho et la surface se simplifient : az ne depend que de Ix, du rapport
# C_La/C_Ma, et de p/V.
#
# Si C_Ma est ESTIME en inversant Miller (c'est le plan PMM), alors
#   C_Ma ~ Ix^2 /(Iy Sg)   donc   az ~ (Iy/Ix) * Sg
# et l'erreur sur la derive est celle du RAPPORT Iy/Ix, multipliee par celle
# de Sg. C'est ce que ce script verifie par integration, sur tout le parcours.

include(joinpath(@__DIR__, "src", "CompetitionBallistics.jl"))
using Printf
import .CompetitionBallistics.SixDOF: mccoy_308_168_shot, initial_spin_rate,
       moments_of_inertia, bullet_inertia, BulletGeometry,
       MCCOY_308_168_CD0, MCCOY_308_168_CLA, MCCOY_308_168_CMA0, MCCOY_308_168_CLP
import .CompetitionBallistics.DragModels: interp_table
import .CompetitionBallistics.BallisticUtils: miller_stability
import .CompetitionBallistics.ExteriorBallistics: spin_drift_inches

const G = 9.80665

"""Derive par integration de l'angle de repos, jusqu'a `portee_yd`.
Rend la derive (m) aux jalons demandes."""
function derive(Ix, Iy, k_cma; jalons = [200.0, 400.0, 600.0, 800.0, 1000.0])
    sp  = mccoy_308_168_shot()
    m   = sp.mass_grains * 0.0647989e-3
    d   = sp.caliber_in * 0.0254
    S   = pi * d^2 / 4
    v0  = sp.muzzle_vel_fps * 0.3048
    tw  = sp.twist_in * 0.0254
    rho = 1.225; a_s = 340.3
    p   = initial_spin_rate(v0, tw)

    x = 0.0; z = 0.0; vz = 0.0; v = v0; t = 0.0
    dt = 1e-4
    out = Dict{Float64,Float64}()
    j = 1
    while j <= length(jalons)
        M   = v / a_s
        CD  = interp_table(MCCOY_308_168_CD0,  M)
        CLa = interp_table(MCCOY_308_168_CLA,  M)
        CMa = interp_table(MCCOY_308_168_CMA0, M) * k_cma
        Clp = interp_table(MCCOY_308_168_CLP,  M)

        # traînee (tir tendu : on neglige l'inclinaison du vecteur vitesse)
        v -= (rho * S * CD * v^2 / (2m)) * dt
        # decroissance du spin
        p += (rho * S * d^2 * v * Clp / (2 * Ix)) * p * dt
        # angle de repos et acceleration laterale
        aR  = 2 * Ix * p * G / (rho * S * d * v^3 * CMa)
        az  = (rho * S * CLa * v^2 / (2m)) * aR
        vz += az * dt
        z  += vz * dt
        x  += v * dt
        t  += dt
        while j <= length(jalons) && x / 0.9144 >= jalons[j]
            out[jalons[j]] = z; j += 1
        end
        t > 5.0 && break
    end
    return out
end

function principal()
    sp = mccoy_308_168_shot()
    m  = sp.mass_grains * 0.0647989e-3
    d  = sp.caliber_in * 0.0254
    L  = sp.bullet_length_in * 0.0254
    v0 = sp.muzzle_vel_fps * 0.3048
    S  = pi * d^2 / 4
    p0 = initial_spin_rate(v0, sp.twist_in * 0.0254)
    CMa0 = interp_table(MCCOY_308_168_CMA0, v0 / 340.3)

    Ixm, Iym = sp.Ix, sp.Iy
    Ixc, Iyc = 0.5m * (d / 2)^2, m * (L^2 / 12 + (d / 2)^2 / 4)
    Ixg, Iyg = moments_of_inertia(m, d, L)

    Sg_mil = miller_stability(mass_gr = sp.mass_grains, caliber_in = sp.caliber_in,
                              bullet_length_in = sp.bullet_length_in,
                              twist_in = sp.twist_in, muzzle_vel_fps = sp.muzzle_vel_fps)
    # C_Ma estime par inversion de Miller, rapporte a la valeur mesuree
    inv(Ix, Iy) = Ix^2 * p0^2 / (2 * 1.225 * Iy * S * d * v0^2 * Sg_mil)

    cas = [("inerties mesurees, C_Ma mesure ", Ixm, Iym, 1.0),
           ("inerties mesurees, C_Ma Miller ", Ixm, Iym, inv(Ixm, Iym) / CMa0),
           ("cylindre plein,    C_Ma Miller ", Ixc, Iyc, inv(Ixc, Iyc) / CMa0),
           ("estimateur geom.,  C_Ma Miller ", Ixg, Iyg, inv(Ixg, Iyg) / CMa0)]

    ref = derive(Ixm, Iym, 1.0)
    jalons = [200.0, 400.0, 600.0, 800.0, 1000.0]
    println("  derive par angle de repos, en cm, et ecart relatif par jalon\n")
    @printf("  %-31s %s\n", "cas", join([@sprintf("%9.0f yd", y) for y in jalons]))
    for (nom, Ix, Iy, k) in cas
        r = derive(Ix, Iy, k)
        @printf("  %-31s %s\n", nom,
                join([@sprintf("%6.1f cm", 100 * r[y]) for y in jalons]))
        if nom != cas[1][1]
            @printf("  %-31s %s\n", "",
                    join([@sprintf("  %+6.1f %%", 100 * (r[y] / ref[y] - 1)) for y in jalons]))
        end
    end

    println("\n  reperes :")
    @printf("    Iy/Ix  mesure %.2f   cylindre %.2f (%+.0f %%)   estimateur %.2f (%+.0f %%)\n",
            Iym/Ixm, Iyc/Ixc, 100*((Iyc/Ixc)/(Iym/Ixm)-1), Iyg/Ixg, 100*((Iyg/Ixg)/(Iym/Ixm)-1))
    @printf("    Sg Miller / Sg vrai : %+.1f %%\n",
            100*(Sg_mil/(Ixm^2*p0^2/(2*1.225*Iym*S*d*v0^2*CMa0)) - 1))
    # formule empirique de Litz, celle que sert le site
    t1000 = 0.0
    let v = v0, x = 0.0, dt = 1e-4
        while x/0.9144 < 1000.0
            v -= (1.225*S*interp_table(MCCOY_308_168_CD0, v/340.3)*v^2/(2m))*dt
            x += v*dt; t1000 += dt
        end
    end
    @printf("    Litz (site) a 1000 yd, Sg Miller %.2f : %.1f cm   (t = %.2f s)\n",
            Sg_mil, 2.54*spin_drift_inches(t1000, Sg_mil, 1), t1000)
end

principal()

println("""
  Lecture :

  1. L'ecart relatif est CONSTANT sur tout le parcours -- identique au dixieme
     de pourcent de 200 a 1000 yd. L'angle de repos croit bien en aval (la
     derive passe de 0,6 a 31,9 cm), mais az ~ (Iy/Ix) * Sg et aucun des deux
     facteurs ne depend de la vitesse : l'erreur est un facteur multiplicatif,
     pas une derive cumulative. Les chiffres calcules a la bouche valaient donc
     deja pour toute la trajectoire.

  2. Le cylindre plein ajoutait 15 cm de derive parasite a 1000 yd (46,6 au lieu
     de 31,9), ce qui compte a cette distance. L'estimateur geometrique tombe a
     1,5 %, soit 0,5 cm -- sous le bruit de tout le reste.

  3. Controle croise independant : la formule empirique de Litz, celle que sert
     pas-de-rayure.php, donne 29,4 cm la ou l'integration de l'angle de repos
     donne 31,9 cm. Deux chemins sans rien de commun se rejoignent a 8 %.
""")
