# Validation du 6-DOF sur un SECOND cas publie, delibérément dissemblable du
# premier : McCoy exemple 9.2, obus de 105 mm HE M1.
#
#   14,97 kg au lieu de 11 g ; tir a 45 et 70 degres de site au lieu du tir
#   tendu ; subsonique (charge 1, Mach 0,602) autant que supersonique (charge 7,
#   Mach 1,449) ; rapport I_x/I_y de 0,101 contre 0,134 ; et un moment de Magnus
#   majoritairement POSITIF la ou celui de la balle est negatif a petit lacet.
#
# Rien ici n'a servi a calibrer quoi que ce soit : les trois corrections du
# 2026-08-13 (signe de la force normale, signe du Magnus, convention BRL) ont
# ete etablies sur le cas 9.1 seul.
#
# Sorties publiees (p. 201, table 9.4, figure 9.11) :
#   charge 1, 45 deg : portee 3760 m,  duree > 28,5 s,  sommet ~1000 m
#   charge 1, 70 deg : portee ~2320 m, duree ~38,5 s,   sommet ~1750 m
#   charge 7, 45 deg : portee ~11500 m, duree > 52,5 s, sommet 3460 m
#   charge 7, 70 deg : portee 7300 m,  duree ~70,5 s,   sommet > 6000 m
#   Sg de bouche 3,1 (1/18) et 1,6 (1/25) ; premier maximum de lacet 3,0 deg
#   figure 9.11 : angle de repos ~2,1 deg a l'apogee, cone lent K_S = 1,9 deg,
#   cycle limite du mode lent ~3 deg sur la fin de la branche descendante,
#   mode rapide (nutation) pratiquement eteint a l'impact.

include(joinpath(@__DIR__, "src", "CompetitionBallistics.jl"))
using Printf
import .CompetitionBallistics.SixDOF: mccoy_105mm_m1_shot, mccoy_105mm_m1_aero,
       solve_6dof, gyroscopic_stability_6dof, initial_spin_rate
import .CompetitionBallistics.DragModels: interp_table
import .CompetitionBallistics.SixDOF: MCCOY_105_CMA0

const D = 0.1048; const IX = 0.02326; const IY = 0.23118
const S_ = pi * D^2 / 4

ec(x, r) = 100 * (x / r - 1)

function stabilite()
    println("=== 1. Facteur de stabilite gyroscopique (table 9.4) ===")
    for (ch, v) in [("charge 1", 205.0), ("charge 7", 493.0)]
        CMa = interp_table(MCCOY_105_CMA0, v / 340.3)
        for (tw, pub) in [(18.0, 3.1), (25.0, 1.6)]
            p = initial_spin_rate(v, tw * D)
            s = gyroscopic_stability_6dof(IX, p, 1.225, S_, D, IY, v, CMa)
            @printf("   %s  1/%.0f : Sg = %.2f   <- publie %.1f  (%+.0f %%)\n",
                    ch, tw, s, pub, ec(s, pub))
        end
    end
end

function trajectoires()
    println("\n=== 2. Trajectoires, pas 1/18 (p. 201) ===")
    println("   cas                portee              duree             sommet          1er lacet")
    cas = [(1, 45, 3760.0, 28.5, 1000.0), (1, 70, 2320.0, 38.5, 1750.0),
           (7, 45, 11500.0, 52.5, 3460.0), (7, 70, 7300.0, 70.5, 6000.0)]
    for (ch, qe, xr, tr_, hr) in cas
        tj = solve_6dof(mccoy_105mm_m1_shot(twist_cal = 18, charge = ch, qe_deg = qe))
        e = tj[end]
        h = maximum(p.y for p in tj)
        y1 = rad2deg(maximum(p.alpha_t for p in tj[1:min(400, end)]))
        @printf("   Ch%d %2d deg   %6.0f m (%+5.1f%%)  %5.1f s (%+5.1f%%)  %5.0f m (%+5.1f%%)  %.2f deg\n",
                ch, qe, e.x, ec(e.x, xr), e.t, ec(e.t, tr_), h, ec(h, hr), y1)
    end
    println("   publie       3760/2320/11500/7300 m ; 28,5/38,5/52,5/70,5 s ;")
    println("                1000/1750/3460/6000 m ; premier maximum 3,0 deg partout")
end

function lacet()
    println("\n=== 3. Mouvement de lacet, charge 1 a 45 deg (figure 9.11) ===")
    tj = solve_6dof(mccoy_105mm_m1_shot(twist_cal = 18, charge = 1, qe_deg = 45))
    h = [p.y for p in tj]
    isom = argmax(h)
    fen(a, b) = (v = [rad2deg(p.alpha_t) for p in tj[a:b]]; (minimum(v), maximum(v)))
    n = length(tj)
    lo_ap, hi_ap = fen(max(1, isom - n ÷ 20), min(n, isom + n ÷ 20))
    lo_fin, hi_fin = fen(max(1, n - n ÷ 10), n)
    lo_deb, hi_deb = fen(1, n ÷ 20)
    @printf("   pres de l'apogee   : alpha_t de %.2f a %.2f deg\n", lo_ap, hi_ap)
    println("      figure 9.11 : cercle de rayon K_S = 1,9 autour de beta_R = 2,1,")
    println("      soit une oscillation attendue entre 0,2 et 4,0 deg")
    @printf("   fin de descente    : alpha_t de %.2f a %.2f deg   <- cycle limite ~3 deg\n",
            lo_fin, hi_fin)
    @printf("   maximum sur tout le vol : %.2f deg\n", maximum(rad2deg(p.alpha_t) for p in tj))
    @printf("   depart (nutation)  : alpha_t de %.2f a %.2f deg\n", lo_deb, hi_deb)
    println("      la nutation doit s'amortir et avoir pratiquement disparu a l'impact")
end

function sommet()
    println("\n=== 4. D'ou vient le residu de portee a 70 degres ===")
    println("   cas         lacet a l'apogee   surcroit de trainee   ecart de portee")
    for (ch, qe, xr) in [(1,45,3760.0), (1,70,2320.0), (7,45,11500.0), (7,70,7300.0)]
        tj = solve_6dof(mccoy_105mm_m1_shot(twist_cal = 18, charge = ch, qe_deg = qe))
        i = argmax([p.y for p in tj]); n = length(tj)
        w = tj[max(1, i - n ÷ 25):min(n, i + n ÷ 25)]
        a = maximum(rad2deg(p.alpha_t) for p in w)
        @printf("   Ch%d %2d deg     %6.2f deg          +%4.0f %% sur C_D0        %+5.1f %%\n",
                ch, qe, a, 3.2 * sind(a)^2 / 0.124 * 100, ec(tj[end].x, xr))
    end
    println("""
   Ce n'est pas la trainee de base : aucun facteur unique sur C_D0 n'aligne les
   quatre cas (a x0,90 le Ch7/70 tombe juste mais le Ch7/45 part a +3,6 %).
   Le residu est propre au site eleve, et la lecture est directe : a 45 degres
   McCoy note que << le lacet de repos au sommet n'est pas grand >>, on trouve
   1,5 et 4,0 degres, et la portee tombe a 0,5 %. A 70 degres la montee
   sommitale porte le lacet a 21-26 degres, ou la trainee de lacet vaut trois a
   cinq fois C_D0 : la portee y depend de l'exactitude d'un lacet de 20 degres
   bien plus que de l'integration. McCoy choisit d'ailleurs ce site comme cas le
   plus defavorable, quelques degres sous l'angle de basculement de l'obus.
   L'effet Eotvos, absent du solveur, ne vaut lui que ~13 m sur 7300.""")
end

stabilite(); trajectoires(); lacet(); sommet()
