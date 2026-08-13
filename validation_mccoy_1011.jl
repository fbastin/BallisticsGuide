# Validation de la theorie linearisee contre McCoy §10.11, qui publie TOUS ses
# intermediaires pour la .308" 168 gr Sierra International a la bouche :
#
#   k_x^-2 = 9.224 ; k_y^-2 = 1.238 ; C_D0 = 0.331 ; C_Ma = 2.63 ; C_La = 2.75 ;
#   C_Mq + C_Mad = -8.2 ; Average C_Mpa = -0.22 ; (pd/V)_0 = 0.16127 rad/cal
#   rho S d / 2m = 0.000021158 ; P = 0.021645
#
#   -> M = 0.000068889   H = 0.00026597   T = 0.000015249   Sg = 1.70
#      Sd = 0.1147       Sg requis > 4.6
#      phi'_F = 0.017768  phi'_S = 0.003877  rad/calibre
#      lambda_F = -0.00031645   lambda_S = +0.00005048  par calibre
#
# ATTENTION au C_Mpa : McCoy prend la MOYENNE SUR LE CYCLE, -0.22, et non la
# valeur a lacet nul de la table, -0.33. Le Magnus est le seul coefficient
# nettement non lineaire a ces incidences, et Sd est une propriete du mouvement.
# Avec -0.33 on obtient Sd = -0.047 : verdict qualitativement oppose.

include(joinpath(@__DIR__, "src", "CompetitionBallistics.jl"))
using Printf
import .CompetitionBallistics.SixDOF: dynamic_stability_6dof, dynamically_stable_6dof

const KX2, KY2 = 9.224, 1.238
const CD0, CMA, CLA, CMQ = 0.331, 2.63, 2.75, -8.2
const PREF, P, M = 0.000021158, 0.021645, 0.000068889

verdict(x, ref, tol) = abs(x/ref - 1) <= tol ? "ok" : "ECART"

function epicyclique(CMpa)
    H = PREF * (CLA - CD0 - KY2 * CMQ)
    T = PREF * (CLA + KX2 * CMpa)
    disc = sqrt(P^2 - 4M)
    phiF, phiS = 0.5*(P + disc), 0.5*(P - disc)
    lF = (-H*phiF + P*T) / (phiF - phiS)
    lS = ( H*phiS - P*T) / (phiF - phiS)
    Sd = dynamic_stability_6dof(CLA, CD0, CMpa, CMQ, KX2, KY2)
    return (; H, T, Sd, phiF, phiS, lF, lS, Sg = P^2/(4M))
end

function principal()
    r = epicyclique(-0.22)
    println("=== McCoy §10.11, C_Mpa moyenne sur le cycle = -0.22 ===")
    println("   grandeur      publie          calcule         ")
    for (nom, calc, pub, tol) in [
            ("H       ", r.H,    0.00026597,   1e-4),
            ("T       ", r.T,    0.000015249,  1e-3),
            ("Sg      ", r.Sg,   1.70,         1e-3),
            ("Sd      ", r.Sd,   0.1147,       1e-3),
            ("phi'_F  ", r.phiF, 0.017768,     1e-3),
            ("phi'_S  ", r.phiS, 0.003877,     1e-3),
            ("lambda_F", r.lF,  -0.00031645,   1e-3),
            ("lambda_S", r.lS,   0.00005048,   1e-3)]
        @printf("   %s  %+13.8f  %+13.8f   %s\n", nom, pub, calc, verdict(calc, pub, tol))
    end
    Sg_req = 1 / (r.Sd * (2 - r.Sd))
    @printf("\n   Sg requis pour la stabilite dynamique : %.2f   <- publie : > 4.6\n", Sg_req)
    @printf("   a 1:12 (Sg %.2f) : %s      a 1:7 (Sg %.2f) : %s\n",
            r.Sg, dynamically_stable_6dof(r.Sg, r.Sd) ? "stable" : "INSTABLE",
            r.Sg*(12/7)^2, dynamically_stable_6dof(r.Sg*(12/7)^2, r.Sd) ? "stable" : "INSTABLE")

    z = epicyclique(-0.33)
    println("\n=== le meme calcul avec la valeur a LACET NUL, -0.33 ===")
    @printf("   Sd = %+.4f  (contre %+.4f)   lambda_S = %+.3e  (contre %+.3e)\n",
            z.Sd, r.Sd, z.lS, r.lS)
    @printf("   Sd < 0 : aucune rotation ne stabiliserait jamais -- verdict oppose.\n")
    println("   C'est pourtant CETTE valeur qu'attendent les formules de cycle limite")
    println("   du §13.7, ou lambda_S = lambda_S0 + lambda_S2 delta^2 porte le lacet")
    println("   explicitement. Les deux usages sont justes, les entrees different.")
end

principal()
