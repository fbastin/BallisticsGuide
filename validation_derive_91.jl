# Dérive du cas 9.1 de McCoy (.308, 168 gr Sierra, 2 600 fps, pas de 12", 1 000 yd) : d'où vient l'écart entre les
# 9,3 po publiés (Modern Exterior Ballistics, 2e éd., p. 198, figure 9.8) et les 12,1 po du 6-DOF sans perturbation ?
#
#     julia validation_derive_91.jl
#
# Trois constats (2026-10-17), chacun calculé ici :
#  1. La perturbation initiale du livre (« an initial yaw rate, or tip-off rate, of 25 radians per second, that started
#     the bullet's nose moving to the left of the line of fire », p. 196) est celle du solveur : `initial_pitch_rate`
#     tourne autour de l'axe 2 du corps, vertical au départ, et positive elle porte le nez à GAUCHE. Sur les
#     orientations de la perturbation, l'écart au tir non perturbé se décompose en une partie TOURNANTE (le saut
#     aérodynamique, opposé à la vitesse initiale du nez, comme le veut la théorie linéaire pour C_Mα > 0) et une partie
#     FIXE, surtout verticale : la chute supplémentaire due à la traînée du lacet de 2 à 5° qu'entretient la perturbation.
#     L'élévation de McCoy (14,68 mrad contre 14,46 sans perturbation) ne désigne donc pas un saut vers le bas.
#  2. La dérive de McCoy ne contient pas le saut : 11,1 / 9,3 / 8,0 po pour 10, 12 et 14" suivent exactement l'inverse
#     du pas, ce qu'un saut, à peu près indépendant de la rotation, casserait ; et les courbes de la figure 9.8 partent
#     à plat de la bouche.
#  3. Intégrée avec la portance de l'annexe A, la courbe de lacet de repos du livre (figure 9.7, relevée au pixel)
#     donne ~9,8 po : figures 9.7 et 9.8 sont cohérentes entre elles. L'écart tient au lacet de repos lui-même : la
#     figure le fait partir de zéro à la bouche et ne rejoindre l'équilibre linéaire (0,0022° dès la bouche) que vers
#     600 yd. Le 6-DOF sans perturbation (12,1 po), qui n'impose pas l'équilibre mais le laisse s'établir, et la formule
#     empirique de Litz (11,6 po) donnent l'équilibre. Comment le livre a obtenu sa courbe n'est pas établi.
#
# Le constat 1 tourne sur une COPIE du module dont l'état initial porte la perturbation sur q et r selon un angle φ
# (q = 25 cos φ, r = 25 sin φ) ; rien n'est modifié dans src/.

using Printf
const PHI = Ref(0.0)
let src = read(joinpath(@__DIR__, "src", "CompetitionBallistics.jl"), String),
    old = "p0, sp.initial_pitch_rate, 0.0)"
    count(old, src) == 1 || error("état initial du 6-DOF introuvable : ce script est à revoir")
    include_string(Main, replace(src, old => "p0, sp.initial_pitch_rate*cos(Main.PHI[]), sp.initial_pitch_rate*sin(Main.PHI[]))"),
                   "CompetitionBallistics (copie, perturbation orientable)")
end
import .CompetitionBallistics.SixDOF: mccoy_308_168_shot, solve_6dof, initial_spin_rate,
       MCCOY_308_168_CD0, MCCOY_308_168_CLA, MCCOY_308_168_CMA0, MCCOY_308_168_CLP
import .CompetitionBallistics.DragModels: interp_table

const PO = 0.0254

function constat1()
    println("1. Écart au tir non perturbé à 1 000 yd, selon l'orientation φ des 25 rad/s (φ = 0 : nez vers la gauche)")
    fin(rate) = solve_6dof(mccoy_308_168_shot(initial_pitch_rate = rate, record_every = 50))[end]
    e0 = fin(0.0)
    @printf("   sans perturbation : dérive %.2f po\n", e0.z / PO)
    dy, dz = Float64[], Float64[]
    for deg in (0, 90, 180, 270)
        PHI[] = deg2rad(deg); e = fin(25.0)
        push!(dy, e.y - e0.y); push!(dz, e.z - e0.z)
        @printf("   φ %3d° : hauteur %+6.2f po, dérive %+6.2f po (dérive totale %5.2f po)\n", deg, dy[end] / PO, dz[end] / PO, e.z / PO)
    end
    PHI[] = 0.0
    fy, fz = sum(dy) / 4, sum(dz) / 4
    @printf("   partie fixe : hauteur %+.2f po, dérive %+.2f po ; partie tournante : %.2f po d'amplitude\n",
            fy / PO, fz / PO, sum(hypot.(dy .- fy, dz .- fz)) / 4 / PO)
end

function constat2()
    println("\n2. Dérives publiées (p. 198) contre l'inverse du pas, depuis 9,3 po au pas de 12\"")
    for (pas, pub) in ((10, 11.1), (12, 9.3), (14, 8.0))
        @printf("   %2d\" : publié %4.1f po, 9,3 × 12/%d = %5.2f po\n", pas, pub, pas, 9.3 * 12 / pas)
    end
end

# Figure 9.7, trait plein (12") : relevé au pixel sur un rendu à 300 ppp (graduations : 350 px pour 0,01°, 211,7 px
# pour 100 yd ; trait du milieu là où les trois courbes se séparent). Yards -> degrés.
const FIG97 = [(0.0, 0.0), (25.0, 0.0004), (75.0, 0.0008), (125.0, 0.0013), (175.0, 0.0022), (225.0, 0.0025),
               (275.0, 0.0029), (325.0, 0.0035), (375.0, 0.0042), (425.0, 0.0051), (475.0, 0.0060), (525.0, 0.0070),
               (575.0, 0.0080), (625.0, 0.0093), (675.0, 0.0104), (725.0, 0.0125), (775.0, 0.0138), (825.0, 0.0153),
               (875.0, 0.0169), (925.0, 0.0194), (975.0, 0.0215), (1000.0, 0.0226)]
interp(tab, x) = (i = findlast(p -> p[1] <= x, tab); i == length(tab) ? tab[end][2] :
                  tab[i][2] + (tab[i+1][2] - tab[i][2]) * (x - tab[i][1]) / (tab[i+1][1] - tab[i][1]))

function constat3()
    println("\n3. Dérive du lacet de repos : figure 9.7 du livre, ou équilibre linéaire (point matériel, annexe A)")
    G = 9.80665; sp = mccoy_308_168_shot()
    m = sp.mass_grains * 0.06479891e-3; d = sp.caliber_in * PO; S = pi * d^2 / 4; Ix = sp.Ix
    v = sp.muzzle_vel_fps * 0.3048; p = initial_spin_rate(v, sp.twist_in * PO); rho = 1.225; a_s = 340.3
    x = 0.0; zF = 0.0; vF = 0.0; zE = 0.0; vE = 0.0; dt = 1e-4
    jal = [10.0, 100.0, 200.0, 400.0, 600.0, 800.0, 1000.0]; j = 1
    while j <= length(jal)
        M = v / a_s
        CD = interp_table(MCCOY_308_168_CD0, M); CLa = interp_table(MCCOY_308_168_CLA, M)
        CMa = interp_table(MCCOY_308_168_CMA0, M); Clp = interp_table(MCCOY_308_168_CLP, M)
        v -= (rho * S * CD * v^2 / (2m)) * dt
        p += (rho * S * d^2 * v * Clp / (2 * Ix)) * p * dt
        bE = 2 * Ix * p * G / (rho * S * d * v^3 * CMa)              # équilibre linéaire (McCoy ch. 10)
        bF = deg2rad(interp(FIG97, x / 0.9144))
        q = rho * S * CLa * v^2 / (2m)
        vF += q * bF * dt; zF += vF * dt; vE += q * bE * dt; zE += vE * dt
        x += v * dt
        if x / 0.9144 >= jal[j]
            @printf("   %5.0f yd : lacet figure %.4f°, équilibre %.4f°  |  dérive figure %5.2f po, équilibre %5.2f po\n",
                    jal[j], rad2deg(bF), rad2deg(bE), zF / PO, zE / PO)
            j += 1
        end
    end
    println("   publié : 9,3 po à 1 000 yd")
end

constat1(); constat2(); constat3()
