# Preakable (a) du 4-DOF : OU la correlation de derive gyroscopique de Litz
# decroche-t-elle du calcul par equations de moment ?
#
# Litz :  dz = 1.25 (Sg + 1.2) t^1.83   pouces
#
# La question n'est pas « est-elle approchee » -- elle l'est par construction --
# mais « existe-t-il un domaine ou elle s'ecarte NOTABLEMENT du calcul physique,
# assez pour qu'un modele de plus se justifie ».
#
# ─────────────────────────────────────────────────────────────────────────────
# POURQUOI CE BANC N'A BESOIN D'AUCUNE DONNEE NEUVE
#
# C_La et C_Ma sont fonctions du Mach SEUL : elles ne dependent pas de la
# rotation. On peut donc faire varier le PAS DE RAYURE d'une balle dont les
# coefficients sont MESURES, et obtenir la derive physique exacte -- sans
# estimateur, sans IntLift, sans les 7 % d'ecart median qui vont avec.
#
# C'est le seul axe du probleme ou l'on puisse balayer Sg sans rien estimer.
#
# ─────────────────────────────────────────────────────────────────────────────
# CE QUE LA THEORIE PREDIT AVANT TOUT CALCUL (McCoy ch. 12)
#
#   angle de repos    aR = 2 Ix p g / (rho S d V^3 C_Ma)
#   accel. laterale   az = g (C_La/C_Ma) (Ix p) / (m d V)
#
# az est LINEAIRE en p. Le taux de decroissance du spin, (rho S d^2 V Clp)/(2 Ix),
# ne depend pas de p : donc p(t) proportionnel a p0 sur toute la trajectoire, et
# V(t) est inchange (le pas ne touche pas la trainee a incidence nulle). Donc
#
#           z(t) proportionnel a p0   EXACTEMENT,  et  Sg proportionnel a p0^2
#           => z proportionnel a sqrt(Sg)   le long de l'axe du pas.
#
# Litz, lui, donne z proportionnel a (Sg + 1.2). Le rapport des deux vaut
#
#           (Sg + 1.2) / sqrt(Sg)
#
# qui est MINIMAL en Sg = 1.2 et croit des deux cotes. Il est plat a 3 % pres
# entre Sg = 0.8 et 2.0 -- la ou Litz a ajuste sa correlation -- et diverge
# au-dela. La prediction est donc : Litz tient tres bien dans sa zone, et
# SUR-ESTIME la derive des balles SUR-STABILISEES.
#
# Ce script verifie cette prediction par integration, sur deux balles mesurees.

include(joinpath(@__DIR__, "src", "CompetitionBallistics.jl"))
using Printf
import .CompetitionBallistics.SixDOF: mccoy_308_168_shot, initial_spin_rate,
       brl3476_aero, BRL3476_PHYSICAL, BRL3476_CALIBER_M,
       MCCOY_308_168_CD0, MCCOY_308_168_CLA, MCCOY_308_168_CMA0, MCCOY_308_168_CLP
import .CompetitionBallistics.DragModels: interp_table
import .CompetitionBallistics.BallisticUtils: miller_stability
import .CompetitionBallistics.ExteriorBallistics: spin_drift_inches

const G   = 9.80665
const RHO = 1.225
const A_S = 340.3

"""Une balle a coefficients MESURES, reduite a ce que la derive demande."""
struct Balle
    nom::String
    m::Float64          # kg
    d::Float64          # m
    L::Float64          # m
    Ix::Float64         # kg m^2
    v0::Float64         # m/s
    cd0::Function       # Mach -> C_D0
    cla::Function       # Mach -> C_La
    cma::Function       # Mach -> C_Ma
    clp::Function       # Mach -> C_lp
    source::String
end

function balle_308()
    sp = mccoy_308_168_shot()
    Balle("«.308 Win 168 gr match»",
          sp.mass_grains * 0.0647989e-3, sp.caliber_in * 0.0254,
          sp.bullet_length_in * 0.0254, sp.Ix, sp.muzzle_vel_fps * 0.3048,
          M -> interp_table(MCCOY_308_168_CD0, M),
          M -> interp_table(MCCOY_308_168_CLA, M),
          M -> interp_table(MCCOY_308_168_CMA0, M),
          M -> interp_table(MCCOY_308_168_CLP, M),
          "McCoy, annexe A du ch. 9 (couloir a etincelles)")
end

function balle_556(which::Symbol = :ss109)
    p = BRL3476_PHYSICAL[which]
    a = brl3476_aero(which)
    d = BRL3476_CALIBER_M
    # C_La = C_Na - C_D0, l'inverse de ce que brl3476_aero assemble.
    # Vitesse initiale prise au Mach le plus haut REELLEMENT tire pour ce
    # projectile (2,63 pour le SS109, soit 2936 fps) : au-dela, brl3476_aero
    # bloque C_D0 au bord de la bande et previent qu'il extrapole. On reste
    # dedans plutot que de faire tourner le banc sur une valeur bloquee.
    Balle("«5,56 mm OTAN $(uppercase(String(which)))»",
          p.mass_g * 1e-3, d, get(p, :len_cal, 4.06) * d, p.Ix * 1e-7,
          2.63 * A_S,
          M -> a.Cd0_func(M),
          M -> a.CNa(M) - a.Cd0_func(M),
          M -> a.CMa(M),
          # ⚠️ C_lp n'est PAS mesure par le BRL-MR-3476 : il n'y figure que dans
          # la liste des symboles. Valeur de service. Sans consequence ici : la
          # decroissance du spin est un facteur commun a tout le balayage.
          _ -> -0.010,
          "BRL-MR-3476 (5,56 mm OTAN, couloir a etincelles d'Aberdeen)")
end

"""Derive par integration de l'angle de repos jusqu'a `portee_m`.
Rend (derive_m, temps_de_vol_s, mach_final)."""
function derive_repos(b::Balle, twist_m::Float64, portee_m::Float64)
    S  = pi * b.d^2 / 4
    p  = initial_spin_rate(b.v0, twist_m)
    x = 0.0; z = 0.0; vz = 0.0; v = b.v0; t = 0.0
    dt = 2e-5
    while x < portee_m && t < 6.0
        M   = v / A_S
        v  -= (RHO * S * b.cd0(M) * v^2 / (2b.m)) * dt
        p  += (RHO * S * b.d^2 * v * b.clp(M) / (2 * b.Ix)) * p * dt
        aR  = 2 * b.Ix * p * G / (RHO * S * b.d * v^3 * b.cma(M))
        az  = (RHO * S * b.cla(M) * v^2 / (2b.m)) * aR
        vz += az * dt
        z  += vz * dt
        x  += v * dt
        t  += dt
    end
    return z, t, v / A_S
end

sg_miller(b::Balle, twist_m) = miller_stability(
    mass_gr = b.m / 0.0647989e-3, caliber_in = b.d / 0.0254,
    bullet_length_in = b.L / 0.0254, twist_in = twist_m / 0.0254,
    muzzle_vel_fps = b.v0 / 0.3048)

litz_m(t, sg) = 0.0254 * spin_drift_inches(t, sg, 1)

# ─────────────────────────────────────────────────────────────────────────────

"""Balayage du pas de rayure a balle et vitesse constantes."""
function balayage(b::Balle, twists_in, portee_m)
    @printf("\n%s — %s\n", b.nom, b.source)
    @printf("  portee %.0f m ; balayage du PAS SEUL (balle, vitesse, air inchanges)\n\n", portee_m)
    println("    pas      Sg       derive physique   Litz      Litz/phys   z/sqrt(Sg)")
    println("   ------  ------   ----------------  --------   ---------   ----------")
    ref = nothing
    lignes = Tuple{Float64,Float64,Float64,Float64}[]
    for tw_in in twists_in
        tw = tw_in * 0.0254
        z, t, _ = derive_repos(b, tw, portee_m)
        sg = sg_miller(b, tw)
        zl = litz_m(t, sg)
        # invariant attendu : z / sqrt(Sg) constant le long de cet axe
        inv = 100z / sqrt(sg)
        ref === nothing && (ref = zl / z)
        @printf("   1:%-4.1f  %5.2f      %7.1f cm      %6.1f cm    %5.2f      %8.2f\n",
                tw_in, sg, 100z, 100zl, zl / z, inv)
        push!(lignes, (sg, z, zl, t))
    end
    println()
    # ecart de Litz rapporte a sa zone d'ajustement (Sg ~ 1.5-2.0)
    i0 = argmin([abs(l[1] - 1.75) for l in lignes])
    k  = lignes[i0][3] / lignes[i0][2]
    @printf("  Recale sur Sg = %.2f (zone d'ajustement de Litz), l'ecart de Litz devient :\n",
            lignes[i0][1])
    print("   ")
    for (sg, z, zl, _) in lignes
        @printf(" Sg %.1f: %+5.1f %% ", sg, 100 * ((zl / z) / k - 1))
    end
    println("\n")
    return lignes
end

function principal()
    println(repeat("=", 78))
    println("  OU LA CORRELATION DE LITZ DECROCHE — banc sur coefficients MESURES")
    println(repeat("=", 78))

    b308 = balle_308()
    b556 = balle_556(:ss109)

    # 1. Verification de l'invariant structurel annonce plus haut.
    println("\n1. CONTROLE : la derive est-elle bien proportionnelle a sqrt(Sg) ?")
    println("   Si oui, la colonne z/sqrt(Sg) est constante — et l'ecart de Litz")
    println("   est alors entierement du a la forme (Sg + 1.2) de sa correlation.")

    l1 = balayage(b308, [16.0, 14.0, 12.0, 11.0, 10.0, 9.0, 8.0, 7.0], 914.4)
    l2 = balayage(b556, [12.0, 10.0, 9.0, 8.0, 7.0, 6.5, 6.0], 550.0)

    for (nom, l) in (("308", l1), ("5,56", l2))
        inv = [l[i][2] / sqrt(l[i][1]) for i in eachindex(l)]
        @printf("  %-5s : z/sqrt(Sg) varie de %+.2f %% sur tout le balayage\n",
                nom, 100 * (maximum(inv) / minimum(inv) - 1))
    end

    # 2. L'AUTRE AXE : changer de BALLE a pas constant.
    #
    # Sur l'axe du pas, Sg bouge par p et la derive suit sqrt(Sg). Sur l'axe de
    # la balle, Sg bouge surtout par C_Ma, et az ~ 1/C_Ma ~ Sg : la pente est
    # DEUX FOIS plus raide. Rien ne dit a priori qu'une correlation a un seul
    # parametre puisse tenir les deux. Litz a ajuste sur des balles, pas sur des
    # pas -- c'est donc l'axe ou son +1.2 a une chance d'etre le bon compromis.
    #
    # Les quatre projectiles du BRL-MR-3476 sont le seul jeu ou l'on dispose de
    # plusieurs balles MESUREES au meme calibre : deux balles ordinaires et deux
    # tracantes, dont les inerties different de 60 % sur Iy.
    println("\n2. L'AUTRE AXE — quatre balles mesurees, pas et vitesse identiques")
    println("   (5,56 mm OTAN, 1:7, portee 550 m)\n")
    #
    # ⚠️ Les deux tracantes n'ont PAS de cote de longueur dans le rapport, et le
    # Sg de Miller en depend au premier ordre. On ne s'en remet donc pas a
    # Miller ici : Ix, Iy et C_Ma etant tous mesures, on calcule le Sg VRAI par
    # sa definition,  Sg = Ix^2 p^2 / (2 rho Iy S d V^2 C_Ma),  et la comparaison
    # ne repose plus sur aucune cote supposee. Le Sg de Miller est garde en
    # regard, parce que c'est lui que le site sert a la formule de Litz.
    println("    balle     Sg vrai  Sg Miller   Ix     Iy      derive phys.  Litz(vrai)  rapport")
    println("   --------  --------  ---------  ------ ------  ------------  ----------  -------")
    for w in (:ss109, :m855, :l110, :m856)
        b  = balle_556(w)
        p  = BRL3476_PHYSICAL[w]
        tw = 7.0 * 0.0254
        z, t, _ = derive_repos(b, tw, 550.0)
        S   = pi * b.d^2 / 4
        p0  = initial_spin_rate(b.v0, tw)
        Iy  = p.Iy * 1e-7
        sgv = b.Ix^2 * p0^2 / (2 * RHO * Iy * S * b.d * b.v0^2 * b.cma(b.v0 / A_S))
        sgm = sg_miller(b, tw)
        zl  = litz_m(t, sgv)
        @printf("   %-8s   %5.2f     %5.2f     %6.4f %6.3f    %6.1f cm     %5.1f cm     %5.2f\n",
                uppercase(String(w)), sgv, sgm, p.Ix, p.Iy, 100z, 100zl, zl / z)
    end

    # 3. L'exposant du temps.
    println("\n3. L'EXPOSANT DU TEMPS — Litz pose t^1.83 ; que rend l'integration ?")
    println("   (.308, pas 1:12, portees croissantes)\n")
    println("    portee    t (s)    derive physique    exposant local")
    println("   --------  -------  ----------------   ---------------")
    prev = nothing
    for r in [200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0]
        z, t, _ = derive_repos(b308, 12.0 * 0.0254, r * 0.9144)
        e = prev === nothing ? NaN : log(z / prev[1]) / log(t / prev[2])
        @printf("   %6.0f yd  %6.3f      %7.1f cm         %s\n", r, t, 100z,
                isnan(e) ? "—" : @sprintf("%.2f", e))
        prev = (z, t)
    end

    # 4. Le meme balayage vu comme Litz le voit.
    println("\n4. LE POINT DE BASCULE")
    println("   Rapport (Sg + 1.2)/sqrt(Sg), normalise a sa valeur en Sg = 1.75 :")
    print("   ")
    for sg in [1.0, 1.2, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]
        r = ((sg + 1.2) / sqrt(sg)) / ((1.75 + 1.2) / sqrt(1.75))
        @printf(" %.1f:%+5.1f%%", sg, 100 * (r - 1))
    end
    println()

    # 5. Ce qu'une correlation de forme correcte donnerait, SANS donnee neuve.
    #
    #   Litz     : 1.25 (Sg + 1.2) t^1.83
    #   candidat : 1.25 (Sg0 + 1.2) sqrt(Sg/Sg0) t^1.83  = 2.79 sqrt(Sg) t^1.83
    #
    # avec Sg0 = 1.75. Meme valeur que Litz dans la zone ou Litz est ajustee,
    # et la dependance en Sg que la theorie impose sur l'axe du pas. C'est un
    # changement d'une ligne : ni modele de plus, ni coefficient a mesurer.
    #
    # ⚠️ Ce n'est PAS applique au site. La formule de Litz est ajustee sur des
    # tirs reels, la ligne ci-dessous sur une theorie du premier ordre : le
    # candidat corrige une forme, il n'apporte pas d'etalonnage neuf.
    litz_corrige(t, sg) = 0.0254 * 1.25 * (1.75 + 1.2) * sqrt(sg / 1.75) * t^1.83
    println("\n5. CE QU'UNE FORME EN sqrt(Sg) RATTRAPERAIT (candidat, non applique)")
    println("   .308 168 gr a 1000 yd, balayage du pas\n")
    println("    pas      Sg      physique    Litz  ecart     candidat  ecart")
    println("   ------  ------   ---------  ------ -------   --------- -------")
    for tw_in in [16.0, 12.0, 10.0, 8.0, 7.0]
        tw = tw_in * 0.0254
        z, t, _ = derive_repos(b308, tw, 914.4)
        sg = sg_miller(b308, tw)
        zl, zc = litz_m(t, sg), litz_corrige(t, sg)
        @printf("   1:%-4.1f  %5.2f    %6.1f cm  %5.1f cm %+5.1f %%   %5.1f cm %+5.1f %%\n",
                tw_in, sg, 100z, 100zl, 100 * (zl / z - 1), 100 * zc, 100 * (zc / z - 1))
    end
end

principal()

println("""

  ─────────────────────────────────────────────────────────────────────────────
  LECTURE

  1. LA PREDICTION STRUCTURELLE EST VERIFIEE, ET EXACTEMENT.  z/sqrt(Sg) ne
     bouge pas de 0,00 % sur les deux balayages, deux calibres, deux jeux de
     coefficients mesures independants. Sur l'axe du pas de rayure, la derive
     par angle de repos est proportionnelle a sqrt(Sg) -- ce n'est pas un
     ajustement, c'est une identite : az est lineaire en p, le taux de
     decroissance du spin ne depend pas de p, et Sg ~ p^2.

  2. LITZ DECROCHE PAR SUR-STABILISATION, PAS PAR SOUS-STABILISATION.  C'est
     l'inverse de ce que le chantier supposait. Le rapport (Sg+1.2)/sqrt(Sg)
     est MINIMAL en Sg = 1,2 et croit des deux cotes ; il est plat a 2 % pres
     de Sg 1,0 a 2,0 -- soit toute la zone des balles de match longue distance
     sur lesquelles Litz a ajuste -- puis diverge : +9 % a Sg 3, +17 % a Sg 4,
     +24 % a Sg 5. La correlation est donc EXCELLENTE la ou elle a ete faite,
     et il n'y a rien a lui reprocher de ce cote.

  3. TOUT L'ECART EN Sg TIENT DANS LA FORME (Sg + 1.2).  La section 5 le montre
     sans ambiguite : en remplacant (Sg+1.2) par 2,95 sqrt(Sg/1,75), le residu
     devient un biais CONSTANT de -7,5 % sur tout le balayage. Un biais constant
     n'est pas une erreur de forme, c'est un ecart d'etalonnage entre une
     correlation ajustee sur des tirs reels et une theorie du premier ordre --
     et 7,5 % tombe dans la dispersion de 10 % deja constatee entre les trois
     methodes (validation_derive.jl).

  4. L'EXPOSANT 1,83 EST TROP BAS POUR LE DOMAINE COURANT.  L'integration rend
     2,00 jusque vers 1000 yd, et ne descend a 1,90 qu'au-dela. Rapporte a un
     etalonnage a 1000 yd, cela fait environ +10 % de derive de trop a 600 yd.
     C'est le meme ordre de grandeur que le defaut en Sg, et il porte, lui, sur
     la zone ou le tir se pratique.

  5. LE VRAI POINT AVEUGLE N'EST PAS Sg, C'EST QUE Sg NE SUFFIT PAS.  Section 2 :
     quatre projectiles du meme calibre, meme pas, meme vitesse. Les deux balles
     ordinaires (SS109, M855) donnent 0,95 et 0,94 -- Litz les traite tres bien.
     Les deux tracantes, a Sg VRAI de 1,98 et 2,07, c'est-a-dire en plein dans la
     zone d'ajustement de Litz, derivent de 14,7 et 12,0 cm quand Litz en annonce
     8,6 et 8,4 : elle en sous-estime une de 42 %. Rien dans Sg ne pouvait le
     voir. La derive va comme (Iy/Ix)*Sg*C_La, et le rapport Iy/Ix passe de 7,8
     sur la SS109 a 11,9 sur la L110 -- une balle a centre de gravite tres
     recule, que Sg ne distingue pas de sa voisine.

     ⚠️ Deux reserves sur ce point, a ne pas gommer. Les tracantes n'ont pas de
     cote de longueur dans le rapport, ce qui a impose de calculer leur Sg par
     sa definition plutot que par Miller -- avec Miller le site les aurait vues
     a Sg 2,47, plus faux encore. Et une tracante n'est pas une balle de sport :
     elle brule, son centre de gravite se deplace en vol, et le BRL la mesure
     precisement parce qu'elle est atypique.

  ─────────────────────────────────────────────────────────────────────────────
  CONSEQUENCE POUR LE 4-DOF

  Le prealable (a) est CLOS, et il ne rouvre pas le chantier.

  Les defauts 2 et 4 sont des defauts d'ALGEBRE de la correlation : ils se
  corrigent dans la formule elle-meme, sans modele de plus, sans coefficient a
  mesurer, sans donnee neuve. Un 4-DOF ne servirait ici a rien qu'une ligne ne
  fasse.

  Le defaut 5 est le seul qui demanderait vraiment un 4-DOF, puisqu'il faut
  connaitre Ix, Iy, C_La et C_Ma par balle pour le voir. Mais c'est aussi celui
  dont la portee pratique est la plus etroite : il se manifeste sur une
  construction atypique, pas entre deux balles de match. Et l'alimenter
  supposerait IntLift, avec ses 7 % d'ecart median et son biais systematique de
  +7,6 % sur C_Ma -- du meme ordre que le defaut qu'on chercherait a corriger.

  ─────────────────────────────────────────────────────────────────────────────
  CE QUE CE BANC NE DIT PAS

  - Il ne compare pas a des TIRS. Les deux membres sont des modeles ; ce qui est
    etabli est un desaccord de forme entre eux, pas laquelle des deux formes
    colle au reel. La theorie de l'angle de repos est du premier ordre en
    incidence, et son etalonnage absolu vaut les ~10 % du controle a trois voies.
  - Le resultat de la section 1 est structurel et ne depend d'aucun etalonnage :
    il tiendrait meme si les deux methodes etaient decalees de 30 %. C'est
    pourquoi les ecarts sont donnes RECALES sur la zone d'ajustement de Litz.
  - Sg > 3 et temps de vol long vont rarement ensemble : une balle sur-stabilisee
    est courte, donc de faible CB, donc rarement tiree loin. L'ecart en % y est
    maximal quand la derive en cm y est petite. Sur la SS109 a 1:7 -- la
    combinaison de service -- l'ecart a 550 m vaut 0,9 cm.
""")
