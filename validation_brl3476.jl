# Validation du jeu 5,56 mm NATO mesure (BRL-MR-3476, McCoy 1985) contre les
# parametres de mouvement que le meme rapport publie en Table 5.
#
# C'est le TROISIEME cas de validation du 6-DOF, et le premier sur arme legere :
# la .308 168 gr est le cas de manuel de McCoy, le 105 mm HE M1 est de l'artillerie,
# et ceci est une balle de service reellement produite. C'est aussi le premier jeu
# de coefficients MESURES pour une arme legere que la bibliotheque possede.
#
# Ce qui est confronte : les tables 3 et 4 donnent les coefficients tir par tir,
# la table 5 donne S_g, S_d, lambda_F, lambda_S, phi'_F et phi'_S deduits des memes
# tirs. Nourrir notre theorie lineaire des premiers doit reproduire les seconds.
#
# LE PAS DU CANON N'EST PAS PUBLIE. Il se deduit : le S_g de 2,49 a Mach 2,730
# impose p = 32150 rad/s, soit 7,15 pouces par tour -- le 1:7" du M16A2 et du M249
# dont le rapport discute la dispersion. C'est une verification croisee du jeu,
# pas une entree libre.

include(joinpath(@__DIR__, "src", "CompetitionBallistics.jl"))
using Printf
import .CompetitionBallistics.SixDOF: brl3476_aero, brl3476_curve, brl3476_geometry,
       bullet_inertia, gyroscopic_stability_6dof, dynamic_stability_6dof,
       BRL3476_SS109, BRL3476_M855, BRL3476_PHYSICAL

const D    = 5.69e-3
const A    = pi * D^2 / 4
const RHO  = 1.225
const ASND = 340.3

# Table 5, colonnes Rd / Mach / S_g / S_d / lambda_F*1e3 / lambda_S*1e3 (1/cal).
# NaN = tiret. Les tirs 144xx viennent d'un canon Obermeyer 7" different.
const T5_SS109 = [
    2.638  2.53   0.79  -0.106  -0.063
    2.622  2.56   0.51  -0.179  -0.044
    1.860  2.36   0.73  -0.125  -0.060
    1.179  2.24  -0.26  -0.167   0.040
    0.736  2.06   0.19  -0.036   0.001
    2.629  2.59   0.74  -0.154  -0.077
    2.625  2.54   0.68  -0.172  -0.073 ]
const T5_M855 = [
    2.730  2.49   0.81  -0.118  -0.072
    2.714  2.50   0.44  -0.155  -0.027
    1.869  2.29   0.77  -0.115  -0.062
    1.072  2.19 -11.60  -0.109   0.096
    0.674  2.18  -2.87  -0.113   0.076 ]

ecart(x, r) = 100 * (x - r) / abs(r)
marque(e, tol) = abs(e) <= tol ? "ok" : "ECART"

# Pas deduit du S_g publie au Mach le plus eleve de chaque serie.
function pas_deduit(which, t5)
    p = BRL3476_PHYSICAL[which]
    Ix, Iy = p.Ix * 1e-7, p.Iy * 1e-7
    ma, sg = t5[1, 1], t5[1, 2]
    v = ma * ASND
    cma = brl3476_aero(which).CMa(ma)
    spin = sqrt(sg * 2 * RHO * A * D * Iy * v^2 * cma / Ix^2)
    return 2pi * v / spin, spin, cma
end

function cas(nom, which, t5)
    p  = BRL3476_PHYSICAL[which]
    m  = p.mass_g * 1e-3
    Ix, Iy = p.Ix * 1e-7, p.Iy * 1e-7
    aero = brl3476_aero(which)
    twist, _, _ = pas_deduit(which, t5)

    println("\n=== $nom ===")
    @printf("   pas deduit du S_g publie : %.2f pouces  (le rapport traite le 1:7\" du M16A2/M249)\n",
            twist / 0.0254)
    @printf("   k_x^-2 = %.3f   k_y^-2 = %.3f\n", m * D^2 / Ix, m * D^2 / Iy)
    println("\n   Mach     S_g pub   S_g nous   ecart      S_d pub   S_d nous   ecart")
    esg = Float64[]
    for i in axes(t5, 1)
        ma, sgr, sdr = t5[i, 1], t5[i, 2], t5[i, 3]
        v  = ma * ASND
        sp = 2pi * v / twist
        sg = gyroscopic_stability_6dof(Ix, sp, RHO, A, D, Iy, v, aero.CMa(ma))
        sd = dynamic_stability_6dof(aero.CNa(ma) - aero.Cd0_func(ma), aero.Cd0_func(ma),
                                    aero.CMpa(ma), aero.CMq_CMad(ma),
                                    m * D^2 / Ix, m * D^2 / Iy)
        push!(esg, ecart(sg, sgr))
        @printf("   %5.3f   %7.2f   %8.2f  %+6.1f%%   %8.2f   %8.2f  %+7.1f%%  %s\n",
                ma, sgr, sg, ecart(sg, sgr), sdr, sd, ecart(sd, sdr), marque(ecart(sg, sgr), 6))
    end
    @printf("\n   S_g : ecart median %+.1f %%, maximum %.1f %%\n",
            sort(esg)[cld(length(esg), 2)], maximum(abs.(esg)))
    return esg
end

function inerties()
    println("\n=== estimateur d'inerties contre la Table 1 (cote de la figure 2) ===")
    println("   balle     Ix        Iy        CG      <- ecart de bullet_inertia")
    for which in (:ss109, :m855)
        p = BRL3476_PHYSICAL[which]
        m, L = p.mass_g * 1e-3, p.len_cal * D
        Ix, Iy, cg, rho = bullet_inertia(m, D, L; geometry = brl3476_geometry(which))
        @printf("   %-8s %+6.1f%%   %+6.1f%%   %+6.1f%%   (rho implicite %.0f kg/m3)\n",
                which, ecart(Ix, p.Ix * 1e-7), ecart(Iy, p.Iy * 1e-7),
                ecart(cg, p.cg_cal * D), rho)
    end
    println("""
   Ces ecarts sont ATTENDUS et ne sont pas un defaut d'integration : le SS-109 et
   le M855 portent un PENETRATEUR D'ACIER devant un noyau de plomb, quand
   BulletGeometry ne connait qu'un noyau unique. C'est la limite du modele sur
   du ball militaire, mesuree ici pour la premiere fois.
   Ne PAS retoucher les defauts la-dessus : leur meilleur score apparent sur ces
   deux balles vient d'une compensation d'erreurs, comme sur le cas McCoy.""")
end

function couverture()
    println("\n=== bande de Mach reellement couverte par les tirs ===")
    for (nom, rows) in (("SS-109", BRL3476_SS109), ("M855", BRL3476_M855))
        for (c, n) in ((4, "C_D"), (5, "C_Ma"), (6, "C_La"), (7, "C_Mpa"), (8, "C_Mq+C_Mad"))
            t = brl3476_curve(rows, c)
            @printf("   %-7s %-11s %d points   Mach %.2f a %.2f\n", nom, n, size(t, 1), t[1, 1], t[end, 1])
        end
    end
    println("""
   Hors de ces bandes, _brl_interp bloque au bord ET previent une fois. Les
   traceurs L110 et M856 sont bien plus maigres (deux points de moment pour le
   L110) : utilisables en supersonique seulement.""")
end

function main()
    println("Validation du jeu 5,56 mm NATO mesure -- BRL-MR-3476 (McCoy, octobre 1985)")
    cas("SS-109 (ball)", :ss109, T5_SS109)
    cas("M855 (ball)",   :m855,  T5_M855)
    inerties()
    couverture()
    println("""

=== ce que ce cas etablit, et ce qu'il n'etablit pas ===
   ETABLI : les coefficients mesures des tables 3 et 4, passes dans notre theorie
   lineaire, reproduisent le S_g publie en table 5. Les entrees et les sorties du
   meme rapport sont coherentes chez nous.
   NON ETABLI : le S_d ne se reproduit pas et ne le peut pas -- il depend du C_Mpa
   MOYEN SUR LE CYCLE, alors que la table 3 donne un tir isole a un lacet isole.
   Meme piege qu'au §10.11 de McCoy, ou la valeur a lacet nul (-0,33) et la moyenne
   de cycle (-0,22) donnent des verdicts opposes. Les colonnes S_d sont imprimees
   ici pour memoire, PAS comme un controle.
   NON ETABLI NON PLUS : aucun cycle limite. Le C_Mpa de ce jeu ne depend que du
   Mach, faute d'assez de tirs pour separer l'effet du lacet.""")
end

main()
