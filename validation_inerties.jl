# Validation de l'estimateur d'inerties depuis la geometrie de balle.
#
# Deux cas publies, independants l'un de l'autre :
#
#   1. .308" 168 gr Sierra International -- McCoy, table 9.2 pour les inerties
#      et le centre de gravite, annexe A du chapitre 9 pour le croquis cote.
#      Balle chemisee a noyau de plomb : eprouve le modele de repartition interne.
#
#   2. Modele 20 mm cone-cylindre en acier -- McCoy, exemple 12.1 (p. 257).
#      Corps homogene de geometrie exacte : eprouve l'integrateur seul, sans
#      aucune hypothese sur la repartition de la masse.
#
# Les deux ont servi a caler l'estimateur : le premier fixe l'epaisseur de
# chemise par defaut, le second ne fixe rien du tout.

include(joinpath(@__DIR__, "src", "CompetitionBallistics.jl"))
using Printf
import .CompetitionBallistics.SixDOF: bullet_inertia, BulletGeometry,
       MCCOY_308_168_GEOMETRY, mccoy_308_168_shot, initial_spin_rate
import .CompetitionBallistics.DragModels: interp_table
import .CompetitionBallistics.SixDOF: MCCOY_308_168_CMA0

const LB_IN2 = 0.45359237 * 6.4516e-4      # lb.in^2 -> kg.m^2
ecart(a, b) = 100 * (a / b - 1)

function cas_sierra()
    println("=== 1. .308\" 168 gr Sierra International (McCoy table 9.2) ===")
    sp = mccoy_308_168_shot()
    m  = sp.mass_grains * 0.0647989e-3
    d  = sp.caliber_in * 0.0254
    L  = sp.bullet_length_in * 0.0254
    Ixr, Iyr, cgr = sp.Ix, sp.Iy, 0.474 * 0.0254

    # l'ancienne approximation, conservee ici comme point de comparaison
    Ixc = 0.5 * m * (d / 2)^2
    Iyc = m * (L^2 / 12 + (d / 2)^2 / 4)

    println("   modele                        Ix        Iy        Iy/Ix     CG")
    @printf("   cylindre plein (ancien)    %+6.1f%%   %+6.1f%%   %+6.1f%%       --\n",
            ecart(Ixc, Ixr), ecart(Iyc, Iyr), ecart(Iyc / Ixc, Iyr / Ixr))
    for (nom, g) in [("contour publie    ", MCCOY_308_168_GEOMETRY),
                     ("contour par defaut", BulletGeometry())]
        Ix, Iy, cg, rho = bullet_inertia(m, d, L; geometry = g)
        @printf("   %s         %+6.1f%%   %+6.1f%%   %+6.1f%%   %+6.1f%%   (rho %.0f kg/m3)\n",
                nom, ecart(Ix, Ixr), ecart(Iy, Iyr),
                ecart(Iy / Ix, Iyr / Ixr), ecart(cg, cgr), rho)
    end

    # ce que l'ecart devient une fois propage sur les grandeurs utiles
    v0 = sp.muzzle_vel_fps * 0.3048
    tw = sp.twist_in * 0.0254
    S  = pi * d^2 / 4
    p0 = initial_spin_rate(v0, tw)
    CMa = interp_table(MCCOY_308_168_CMA0, v0 / 340.3)
    Sg(Ix, Iy) = Ix^2 * p0^2 / (2 * 1.225 * Iy * S * d * v0^2 * CMa)
    aR(Ix, Iy) = 4 * Iy * Sg(Ix, Iy) * 9.80665 / (Ix * p0 * v0)

    println("\n   propagation (Sg publie 1.70) :")
    for (nom, Ix, Iy) in [("inerties mesurees ", Ixr, Iyr),
                          ("cylindre plein    ", Ixc, Iyc),
                          ("contour par defaut",
                           bullet_inertia(m, d, L)[1:2]...)]
        @printf("     %s  Sg = %.3f (%+5.1f%%)   angle de repos = %.4f mrad (%+5.1f%%)\n",
                nom, Sg(Ix, Iy), ecart(Sg(Ix, Iy), Sg(Ixr, Iyr)),
                1e3 * aR(Ix, Iy), ecart(aR(Ix, Iy), aR(Ixr, Iyr)))
    end
end

function cas_cone_cylindre()
    println("\n=== 2. modele 20 mm cone-cylindre acier (McCoy exemple 12.1) ===")
    # d = 0.786", L = 3.93" : nez conique 2 cal + cylindre 3 cal, acier homogene
    d = 0.786 * 0.0254
    L = 3.93 * 0.0254
    m = 0.392 * 0.45359237
    Ixr = 0.0282 * LB_IN2
    Iyr = 0.3263 * LB_IN2
    cgr = 1.47 * 0.0254

    # un cone est une ogive de rayon infini : on l'obtient en integrant
    # directement, l'estimateur ne modelisant que les ogives tangentes.
    r(z) = z <= 3d ? d / 2 : (d / 2) * (1 - (z - 3d) / (2d))
    n = 200_000; h = L / n
    w(i) = (i == 0 || i == n) ? 1.0 : (isodd(i) ? 4.0 : 2.0)
    v = m1 = q = 0.0
    for i in 0:n
        z = i * h; y2 = r(z)^2; c = w(i)
        v += c * y2; m1 += c * y2 * z; q += c * y2 * y2
    end
    cg = m1 / v
    rho = m / (pi * (h / 3) * v)
    Ix = rho * pi * (h / 3) * q / 2
    t = 0.0
    for i in 0:n
        z = i * h; y2 = r(z)^2
        t += w(i) * (y2 * y2 / 4 + y2 * (z - cg)^2)
    end
    Iy = rho * pi * (h / 3) * t

    @printf("   densite implicite %.0f kg/m3   (acier 7850 ; McCoy retire 0.007 lb en percant)\n", rho)
    @printf("   CG   %+6.1f%%      Ix  %+6.1f%%      Iy  %+6.1f%%      Iy/Ix  %+6.1f%%\n",
            ecart(cg, cgr), ecart(Ix, Ixr), ecart(Iy, Iyr), ecart(Iy / Ix, Iyr / Ixr))
    println("   -> l'integrateur est juste sur un corps homogene de forme connue ;")
    println("      ce qui restait faux dans l'ancienne estimation, c'etait la forme.")
end

cas_sierra()
cas_cone_cylindre()
