# Validation de l'estimateur d'inerties sur les trois balles de match 7,62 mm de
# BRL-MR-3733 (McCoy, decembre 1988).
#
# POURQUOI CE CAS EST LE MEILLEUR DU LOT. Partout ailleurs le contour et l'inertie
# viennent de deux sources : la 190 Sierra appariait un contour de 2015 avec une
# mesure des annees 1980, et Sierra a re-outille entre les deux. Ici le rapport dit
# lui-meme que CINQ exemplaires de chacun des trois types ont ete mesures, que la
# table 1 en donne la moyenne, et que les figures 3, 4 et 5 en donnent le croquis
# cote. Meme lot, meme campagne, un seul document.
#
# CE QUE LE CROQUIS PORTE DE PLUS. Chaque figure cote le rayon d'ogive ET, quand
# l'arc n'est pas tangent, la position axiale de son centre. Deux des trois ogives
# sont SECANTES. BulletGeometry n'a pas de champ pour un rayon : `_profile` impose
# la tangente. Ce fichier mesure ce que cette approximation coute, en integrant
# l'arc reel a cote -- comme validation_inerties.jl integre deja directement le
# cone-cylindre que le module ne sait pas representer.
#
# VERDICT (2026-08-18) : l'ecart est reel sur la 190 mais ne justifie PAS d'ajouter
# le parametre. Voir le bloc final, et la note de BulletGeometry.

include(joinpath(@__DIR__, "src", "CompetitionBallistics.jl"))
using Printf
import .CompetitionBallistics.SixDOF: bullet_inertia, BulletGeometry,
       BRL3733_PHYSICAL, brl3733_geometry

const D = 7.82e-3                       # calibre de reference du rapport
const ORDRE = (:m118, :sierra190, :sierra168)
const NOMS = Dict(:m118 => "M118", :sierra190 => "190 Sierra", :sierra168 => "168 Sierra")

ecart(x, r) = 100 * (x / r - 1)

# --- ogive secante : arc de rayon R passant par (0, r_meplat) et (L_nez, 1/2) ----
# Rendu (R, x_centre, y_centre) en calibres, x compte depuis la POINTE -- c'est le
# repere du croquis, celui dans lequel BRL cote la position du centre.
#   xc^2 + (rm - yc)^2 = R^2  et  (nose - xc)^2 + (0.5 - yc)^2 = R^2
# La difference des deux donne xc affine en yc, qu'on reinjecte dans la premiere.
function arc_secant(nose, rm, R)
    p = (rm - 0.5) / nose
    q = (nose^2 - rm^2 + 0.25) / (2nose)          # xc = q + p*yc
    a = p^2 + 1
    b = 2 * p * q - 2 * rm
    c = q^2 + rm^2 - R^2
    disc = b^2 - 4a * c
    disc < 0 && return nothing
    yc = (-b - sqrt(disc)) / (2a)       # racine basse : centre sous le profil
    return (R, q + p * yc, yc)
end

rayon_tangent(nose, rm) = (nose^2 + (0.5 - rm)^2) / (2 * (0.5 - rm))

# --- integrateur local, seul endroit qui sait dessiner une ogive secante --------
# Meme modele interne que bullet_inertia : chemise d'epaisseur radiale constante,
# noyau de plomb montant du culot, cavite au-dessus, hauteur de noyau resolue pour
# rendre la masse. Rien n'est ajuste sur les inerties.
const RHO_PB, RHO_CH = 11340.0, 8860.0

function inertie_secante(m_kg, d_m, L_cal, nose, mep, btl, btd, R)
    rm = mep / 2
    arc = arc_secant(nose, rm, R)
    arc === nothing && error("arc impossible")
    _, xc, yc = arc
    zsh = L_cal - nose                              # z compte depuis le CULOT
    rbase = 0.5 - btl * tand(btd)
    function rayon(z)
        if btl > 0 && z <= btl
            return rbase + (0.5 - rbase) * (z / btl)
        elseif z <= zsh
            return 0.5
        else
            x = L_cal - z                           # distance a la pointe
            return yc + sqrt(max(R^2 - (x - xc)^2, 0.0))
        end
    end
    simpson(f, a, b, n) = b <= a ? 0.0 : begin
        h = (b - a) / n; s = f(a) + f(b)
        for i in 1:n-1; s += (isodd(i) ? 4.0 : 2.0) * f(a + i * h); end
        s * h / 3
    end
    integ(f, brk) = begin
        bn = sort(unique(clamp.(vcat(0.0, brk, L_cal), 0.0, L_cal)))
        sum(simpson(f, bn[i], bn[i+1], 400) for i in 1:length(bn)-1)
    end
    t = 0.095
    d3 = d_m^3; d5 = d_m^5
    lin(z, zc) = begin
        ro = rayon(z); ri = max(ro - t, 0.0)
        (ro^2 - ri^2) * RHO_CH + ri^2 * (z <= zc ? RHO_PB : 0.0)
    end
    qua(z, zc) = begin
        ro = rayon(z); ri = max(ro - t, 0.0)
        (ro^4 - ri^4) * RHO_CH + ri^4 * (z <= zc ? RHO_PB : 0.0)
    end
    # hauteur de noyau par dichotomie, pour rendre exactement la masse donnee
    lo, hi = 0.0, L_cal
    for _ in 1:80
        mid = (lo + hi) / 2
        (pi * d3 * integ(z -> lin(z, mid), [btl, zsh, mid]) < m_kg) ? (lo = mid) : (hi = mid)
    end
    zc = (lo + hi) / 2
    br = [btl, zsh, zc]
    mtot = pi * d3 * integ(z -> lin(z, zc), br)
    cg = integ(z -> z * lin(z, zc), br) / integ(z -> lin(z, zc), br)
    Ix = (pi / 2) * d5 * integ(z -> qua(z, zc), br)
    Iy = (pi / 4) * d5 * integ(z -> qua(z, zc), br) +
         pi * d5 * integ(z -> (z - cg)^2 * lin(z, zc), br)
    return (Ix, Iy, cg * d_m, mtot)
end

# ------------------------------------------------------------------------------
function integrateur_conforme()
    println("=== 0. L'integrateur local est-il le meme modele que le module ? ===")
    println("   Nourri du rayon TANGENT, il doit reproduire bullet_inertia. Sinon la")
    println("   colonne (a) ci-dessous ne mesurerait pas l'arc mais deux codes.")
    ok = true
    for w in ORDRE
        p = BRL3733_PHYSICAL[w]
        m = p.mass_g * 1e-3
        Rt = rayon_tangent(p.nose_cal, p.meplat_cal / 2)
        a = inertie_secante(m, D, p.len_cal, p.nose_cal, p.meplat_cal, p.bt_cal, p.bt_deg, Rt)
        b = bullet_inertia(m, D, p.len_cal * D; geometry = brl3733_geometry(w))
        e = (ecart(a[1], b[1]), ecart(a[2], b[2]), ecart(a[3], b[3]))
        ok &= all(abs.(e) .< 0.05)
        @printf("   %-12s  Ix %+.3f %%   Iy %+.3f %%   CG %+.3f %%\n", NOMS[w], e...)
    end
    println(ok ? "   -> conforme.\n" : "   -> ECART : les deux integrateurs divergent, tout ce qui suit est suspect.\n")
end

function croquis_coherent()
    println("=== 1. Le croquis se verifie-t-il lui-meme ? ===")
    println("   Le rayon imprime doit replacer le centre d'arc a la cote imprimee.")
    println("   balle         R cote   R tangent   verdict    centre calcule / cote")
    for w in ORDRE
        p = BRL3733_PHYSICAL[w]
        rm = p.meplat_cal / 2
        Rt = rayon_tangent(p.nose_cal, rm)
        sec = p.ogive_R_cal > 1.01 * Rt
        @printf("   %-12s  %5.2f     %5.2f     %-9s", NOMS[w], p.ogive_R_cal, Rt,
                sec ? "SECANTE" : "tangente")
        if isnan(p.ogive_xc_cal)
            println("  -- (tangente : le centre est au raccord, non cote)")
        else
            _, xc, _ = arc_secant(p.nose_cal, rm, p.ogive_R_cal)
            @printf("  %.3f / %.2f  (%+.1f %%)\n", xc, p.ogive_xc_cal,
                    ecart(xc, p.ogive_xc_cal))
        end
    end
    println("""
   -> Les cotes sont redondantes et concordantes. Le croquis n'est pas lu de
      travers : la 190 tombe a 0,1 %, le M118 a 1,3 % (arrondi du dessin).""")
end

function inerties()
    println("\n=== 2. Inerties estimees contre la Table 1 du meme rapport ===")
    println("   Trois colonnes, de la plus informee a la plus generique :")
    println("     (a) contour vrai + arc reel   -- integrateur local, secante comprise")
    println("     (b) contour vrai + tangente   -- ce que le module sait faire")
    println("     (c) contour par defaut        -- ce que 97,7 % de la base utilise")
    println()
    println("   balle              (a) Ix      Iy   |  (b) Ix      Iy   |  (c) Ix      Iy")
    acc = Dict(:a => Float64[], :b => Float64[], :c => Float64[])
    for w in ORDRE
        p = BRL3733_PHYSICAL[w]
        m = p.mass_g * 1e-3
        L = p.len_cal * D
        Ixr, Iyr = p.Ix * 1e-7, p.Iy * 1e-7
        a = inertie_secante(m, D, p.len_cal, p.nose_cal, p.meplat_cal,
                            p.bt_cal, p.bt_deg, p.ogive_R_cal)
        b = bullet_inertia(m, D, L; geometry = brl3733_geometry(w))
        c = bullet_inertia(m, D, L)
        for (k, v) in ((:a, a), (:b, b), (:c, c))
            push!(acc[k], ecart(v[1], Ixr), ecart(v[2], Iyr))
        end
        @printf("   %-12s %+7.1f%% %+7.1f%% | %+7.1f%% %+7.1f%% | %+7.1f%% %+7.1f%%\n",
                NOMS[w], ecart(a[1], Ixr), ecart(a[2], Iyr),
                ecart(b[1], Ixr), ecart(b[2], Iyr),
                ecart(c[1], Ixr), ecart(c[2], Iyr))
    end
    rms(v) = sqrt(sum(v .^ 2) / length(v))
    @printf("\n   RMS sur les six ecarts :  (a) %.2f   (b) %.2f   (c) %.2f\n",
            rms(acc[:a]), rms(acc[:b]), rms(acc[:c]))
    println("""
   -> La 168, seule ogive reellement tangente, donne la MEME valeur en (a) et en
      (b) : c'est le temoin integre, et il dit que l'ecart entre ces deux colonnes
      est bien l'arc et rien d'autre.
   -> (b) est PIRE que (c). Donner la vraie longueur d'ogive en imposant une
      tangente est moins bon que de tout laisser generique : la geometrie a moitie
      juste perd contre la geometrie franchement generique.""")
end

function cg_et_masse()
    println("\n=== 3. Centre de gravite, et le residu du M118 ===")
    println("   balle           CG estime   CG mesure     ecart     masse/enveloppe")
    ref = nothing
    for w in ORDRE
        p = BRL3733_PHYSICAL[w]
        m = p.mass_g * 1e-3
        a = inertie_secante(m, D, p.len_cal, p.nose_cal, p.meplat_cal,
                            p.bt_cal, p.bt_deg, p.ogive_R_cal)
        cgr = p.cg_cal * D
        @printf("   %-12s  %7.2f mm   %7.2f mm   %+6.1f%%      %.3f g/cal\n",
                NOMS[w], a[3] * 1e3, cgr * 1e3, ecart(a[3], cgr),
                p.mass_g / p.len_cal)
    end
    # la masse que la longueur du M118 laisserait attendre, interpolee 168 -> 190
    p1, p2, pm = BRL3733_PHYSICAL[:sierra168], BRL3733_PHYSICAL[:sierra190], BRL3733_PHYSICAL[:m118]
    att = p1.mass_g + (p2.mass_g - p1.mass_g) *
          (pm.len_cal - p1.len_cal) / (p2.len_cal - p1.len_cal)
    dIy = ecart(inertie_secante(pm.mass_g * 1e-3, D, pm.len_cal, pm.nose_cal,
                                pm.meplat_cal, pm.bt_cal, pm.bt_deg, pm.ogive_R_cal)[2],
                pm.Iy * 1e-7)
    @printf("""

   Le M118 pese %.2f g quand son enveloppe en laisse attendre %.2f g (interpolation
   lineaire 168 -> 190 sur la longueur), soit %.2f g de moins. Notre modele evacue
   cette masse en descendant le noyau, donc il alourdit le culot : c'est exactement
   le signe de l'ecart de CG ci-dessus, et de son %+.0f %% sur Iy.
   RESIDU OUVERT : il tient a la construction interne, que MR-3733 ne decrit pas.
   Ne PAS le rattraper en retouchant jacket_cal -- balaye sur les trois balles, ce
   parametre ne fait qu'echanger le M118 contre les deux autres.
""", pm.mass_g, att, att - pm.mass_g, dIy)
end

function verdict()
    println("""
=== ce que ce cas etablit, et ce qu'il n'etablit pas ===
   ETABLI : deux des trois ogives sont secantes, les croquis le disent en cotant le
   centre d'arc, et les lire correctement fait passer la 190 de -3,5 % a +0,7 % sur
   Iy. L'ogive n'est donc PAS hors de cause -- la conclusion inverse du 2026-08-15
   reposait sur une longueur de nez fausse (2,69 cal au lieu de 2,09) tiree d'une
   base commerciale.

   ETABLI AUSSI, et c'est ce qui decide : le parametre ne sera PAS ajoute.
     - S_g est exactement proportionnel a 1/Iy, donc le gain vaut ~4 % de S_g sur
       une seule balle, dans une estimation qui porte deja +20/-16 % venus de la
       constante de famille masse<->longueur.
     - Le contour par defaut fait aussi bien que le contour vrai sur ces trois
       balles : Iy est fixe pour l'essentiel par la masse et la longueur.
     - Le rayon est introuvable en production. Sur 3528 balles : 2,3 % portent une
       longueur d'ogive, 0,1 % un meplat, TROIS le contour complet. Aucune source
       exploitable ne publie le rayon -- pas meme la table Berger, la plus riche
       dont nous disposions.
     - Et un champ que personne ne peut remplir finit rempli depuis un curseur de
       forme qui n'est pas un rayon.

   NON ETABLI : le residu du M118 (-10 % sur Iy, -9 % sur le CG). Ni le contour ni
   la chemise ne l'atteignent. C'est aujourd'hui le plus gros ecart du lot.

   MORT POUR DE BON : la chemise fuselee. Elle survivait a n = 2 et paraissait
   physique ; balayee sur les trois balles, l'optimum est tau = 1,0, c'est-a-dire
   pas de fuselage. C'etait un artefact de deux points.""")
end

function main()
    println("Validation sur BRL-MR-3733 -- McCoy, The Aerodynamic Characteristics of")
    println("7.62mm Match Bullets, decembre 1988 -- DTIC ADA205633.")
    println("Oeuvre du gouvernement federal americain : sans copyright US (17 USC 105),")
    println("scan DTIC et non Google. La mention <public release> du rapport est un")
    println("controle de DIFFUSION, pas un enonce de droit d'auteur -- ne pas confondre.\n")
    integrateur_conforme()
    croquis_coherent()
    inerties()
    cg_et_masse()
    verdict()
end

main()
