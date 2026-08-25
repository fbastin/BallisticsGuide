#!/usr/bin/env python3
"""Contrôle de parité entre le modèle Julia de référence et son portage JavaScript.

    python3 scripts/check_js_parity.py           # sort en 1 si divergence
    python3 scripts/check_js_parity.py --verbose # affiche chaque cas

POURQUOI. `src/CompetitionBallistics.jl` est la référence ;
`CompetitionBallistics.js` est ce que le site SERT. Les deux sont écrits à la main
et rien ne les confrontait. Le 2026-08-13, `miller_stability` s'est révélée fausse
**du côté Julia** — correction de vitesse en rapport direct au lieu de racine
cubique, et correction atmosphérique inversée — pendant que le JS, lui, était
juste. Personne ne l'a vu parce que rien ne regardait. Un an de divergence
possible sur la formule la plus consultée de la bibliothèque.

CE QUI EST COMPARÉ, ET POURQUOI PAS DES CONSTANTES. `check_model_drift.py` compare
les constantes du modèle de vibrations, ce qui convient à deux simulations. Ici on
compare des **SORTIES sur une grille d'entrées** : une erreur de formule ne se voit
pas dans les constantes, et c'est précisément une erreur de formule qui avait
échappé. Les deux implémentations sont exécutées une fois chacune, sur les mêmes
cas, et les résultats confrontés.

LA CORRESPONDANCE EST DÉCLARÉE, PAS DEVINÉE. `CAS` nomme explicitement la fonction
Julia, la fonction JS et les arguments. Deviner l'appariement sur les noms
reviendrait à tester que deux conventions de nommage coïncident, pas que deux
physiques coïncident — et un nom qui ne s'apparie pas passerait silencieusement
pour « rien à comparer », le pire des résultats.

TOLÉRANCE. 1e-9 en relatif. Les deux langages font du flottant 64 bits sur les
mêmes expressions : au-delà de l'ordre des opérations, tout écart est un écart de
formule. Une tolérance lâche laisserait passer exactement ce qu'on cherche.

ET LES DONNÉES, DEPUIS LE 2026-08-25. `COMMON_BULLETS` et `CARTRIDGE_MAX_PSI`
sont des tables de valeurs, pas des formules : la parité des fonctions ne les
voyait pas. Le CB faux de la Sierra 175 MK (0.259 contre 0.243) vivait des deux
côtés de ce contrôle sans jamais le déclencher — il a été corrigé à la main le
2026-08-25, par discipline et non par garde-fou. Les profils et les pressions
max sont désormais confrontés champ à champ, et `--selftest` perturbe un CB de
1 % pour prouver que le garde mord.
"""
import argparse
import json
import math
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
TOL = 1e-9

# Chaque cas : (libellé, appel Julia, appel JS, dict d'arguments).
# Les appels sont des GABARITS où {nom} est remplacé par la valeur de l'argument.
# Julia prend des mots-clés là où le JS prend des positionnels : d'où deux gabarits
# et non une signature commune.
CAS = [
    ("miller_stability",
     "BU.miller_stability(mass_gr={m}, caliber_in={d}, bullet_length_in={L}, "
     "twist_in={t}, muzzle_vel_fps={v}, temp_F={T}, pressure_inhg={P})",
     "BU.millerStability({{massGr:{m}, caliberIn:{d}, bulletLengthIn:{L}, "
     "twistIn:{t}, muzzleVelFps:{v}, tempF:{T}, pressureInhg:{P}}})",
     [dict(m=168, d=0.308, L=1.215, t=10, v=2600, T=59, P=29.92),
      dict(m=55,  d=0.224, L=0.760, t=9,  v=3200, T=59, P=29.92),
      dict(m=300, d=0.338, L=1.819, t=9.4, v=2750, T=59, P=29.92),
      # Hors des conditions standard : c'est là que vivaient les deux facteurs faux.
      dict(m=168, d=0.308, L=1.215, t=10, v=4000, T=59, P=29.92),
      dict(m=168, d=0.308, L=1.215, t=10, v=1800, T=59, P=29.92),
      dict(m=168, d=0.308, L=1.215, t=10, v=2600, T=100, P=29.92),
      dict(m=168, d=0.308, L=1.215, t=10, v=2600, T=-20, P=29.92),
      dict(m=168, d=0.308, L=1.215, t=10, v=2600, T=59, P=20.0),
      dict(m=168, d=0.308, L=1.215, t=10, v=2600, T=59, P=31.0)]),

    ("greenhill_twist",
     "BU.greenhill_twist({d}, {L}; C={C})",
     "BU.greenhillTwist({d}, {L}, {C})",
     [dict(d=0.308, L=1.215, C=150), dict(d=0.224, L=0.760, C=180),
      dict(d=0.338, L=1.819, C=150)]),

    ("spin_drift_inches",
     "EB.spin_drift_inches({t}, {sg}, {dir})",
     "EB.spinDriftInches({t}, {sg}, {dir})",
     [dict(t=1.0, sg=1.5, dir=1), dict(t=2.5, sg=2.4, dir=1),
      dict(t=1.8, sg=1.9, dir=-1)]),

    ("coriolis_horizontal",
     "EB.coriolis_horizontal({r}, {tof}, {lat})",
     "EB.coriolisHorizontal({r}, {tof}, {lat})",
     [dict(r=1000, tof=1.5, lat=50), dict(r=500, tof=0.7, lat=-33),
      dict(r=1800, tof=3.1, lat=0)]),

    ("eotvos_vertical",
     "EB.eotvos_vertical({r}, {tof}, {lat}, {az})",
     "EB.eotvosVertical({r}, {tof}, {lat}, {az})",
     [dict(r=1000, tof=1.5, lat=50, az=90), dict(r=1000, tof=1.5, lat=50, az=270),
      dict(r=800, tof=1.1, lat=-33, az=45)]),

    ("wind_deflection_lag_rule",
     "EB.wind_deflection_lag_rule({w}, {tof}, {r}, {v})",
     "EB.windDeflectionLagRule({w}, {tof}, {r}, {v})",
     [dict(w=4.5, tof=1.5, r=1000, v=800), dict(w=-2.0, tof=0.9, r=500, v=850)]),

    ("aerodynamic_jump",
     "EB.aerodynamic_jump({c}, {k}, {a}, {p}, {r})",
     "EB.aerodynamicJump({c}, {k}, {a}, {p}, {r})",
     [dict(c=2.8, k=0.9, a=0.002, p=0.02, r=0.004)]),

    ("form_factor",
     "BU.form_factor({m}, {d}, {bc})",
     "BU.formFactor({m}, {d}, {bc})",
     [dict(m=168, d=0.308, bc=0.462), dict(m=140, d=0.264, bc=0.610)]),

    ("drop_to_moa",
     "BU.drop_to_moa({drop}, {r})",
     "BU.dropToMoa({drop}, {r})",
     [dict(drop=-12.5, r=300), dict(drop=-88.0, r=1000)]),

    ("air_density",
     "Atm.air_density(P={P}, T={T}, H={H})",
     "Atm.airDensity({{P:{P}, T:{T}, H:{H}}})",
     [dict(P=101325, T=288.15, H=0), dict(P=101325, T=288.15, H=80),
      dict(P=85000, T=303.15, H=50)]),

    ("density_altitude",
     "Atm.density_altitude({rho})",
     "Atm.densityAltitude({rho})",
     [dict(rho=1.225), dict(rho=1.05), dict(rho=1.35)]),

    ("cd_g1",
     "DM.cd_g1({mach})",
     "DM.cdG1({mach})",
     [dict(mach=m) for m in (0.5, 0.9, 1.0, 1.05, 1.2, 2.0, 3.0)]),

    ("cd_g7",
     "DM.cd_g7({mach})",
     "DM.cdG7({mach})",
     [dict(mach=m) for m in (0.5, 0.9, 1.0, 1.05, 1.2, 2.0, 3.0)]),

    # Les splines, ENTRE les nœuds. Échantillonner sur les nœuds ne teste rien : c'est
    # exactement l'erreur qui m'a fait conclure « la traînée concorde » alors que le
    # Julia interpolait linéairement et le JS par spline — aux nœuds les deux méthodes
    # coïncident par construction.
    ("cd_g7_spline (entre nœuds)",
     "DM.cd_g7_spline({mach})",
     "DM.cdG7Spline({mach})",
     [dict(mach=m) for m in (0.412, 0.873, 1.017, 1.133, 1.677, 2.317, 2.866, 4.049)]),

    ("cd_g1_spline (entre nœuds)",
     "DM.cd_g1_spline({mach})",
     "DM.cdG1Spline({mach})",
     [dict(mach=m) for m in (0.412, 0.873, 1.017, 1.133, 1.677, 2.317, 2.866, 4.049)]),

    ("drag_coefficient (défaut)",
     "DM.drag_coefficient(:G7, {mach})",
     "DM.dragCoefficient('G7', {mach})",
     [dict(mach=m) for m in (0.412, 0.873, 1.017, 2.866)]),

    # Le solveur complet, sur ses quatre sorties d'arrivée. C'est le cas qui compte
    # le plus : il fait travailler ensemble atmosphère, traînée, angle de zérotage,
    # Coriolis et vent, et une divergence n'importe où s'y voit. Il a d'ailleurs
    # découvert que le Julia ne s'exécutait pas du tout — `g_along` et les trois
    # composantes `oe_*` n'étaient jamais affectées, depuis 4c2343c.
    #
    # ⚠️ ON COMPARE À INDICE FIXE, PAS AU DERNIER POINT. Les deux boucles s'arrêtent
    # sur `x <= target_range`, et un écart d'intégration infime décide de faire ou non
    # un pas de plus : sur une .224 à 500 yd, Julia enregistre 1063 points et le JS
    # 1062, soit 0,5310 s contre 0,5305 s. Comparer « le dernier point » comparait donc
    # DEUX INSTANTS DIFFÉRENTS, et rendait 0,83 m/s d'écart — exactement la
    # décélération d'un pas. Ce n'était pas une divergence de physique mais un artefact
    # de bord, et il aurait fait passer une vraie divergence pour du bruit connu.
    # Les deux enregistrent à chaque pas depuis t=0 : l'indice 500 du JS et le 501 du Julia (indexation à 0 contre 1) sont donc le même
    # instant des deux côtés, et la comparaison redevient exacte.
    ("solve_trajectory · chute @500",
     "EB.solve_trajectory(EB.ShotParameters(mass_grains={m}, caliber_in={d}, bc={bc}, "
     "muzzle_vel_fps={v}, zero_range_yd={zr}, target_range_yd={tr}, latitude_deg={lat}, "
     "azimuth_deg={az}, wind_speed_mph={w}, incline_deg={inc}, humidity_pct={hum}, enable_coriolis=true))[501].drop_m",
     "EB.solveTrajectory({{...EB.DEFAULT_SHOT, massGrains:{m}, caliberIn:{d}, bc:{bc}, "
     "muzzleVelFps:{v}, zeroRangeYd:{zr}, targetRangeYd:{tr}, latitudeDeg:{lat}, "
     "azimuthDeg:{az}, windSpeedMph:{w}, inclineDeg:{inc}, humidityPct:{hum}, enableCoriolis:true}})[500].dropM",
     [dict(m=168, d=0.308, bc=0.462, v=2600, zr=100, tr=1000, lat=50, az=90, w=0, inc=0, hum=0),
      dict(m=168, d=0.308, bc=0.462, v=2600, zr=100, tr=1000, lat=50, az=270, w=0, inc=0, hum=0),
      dict(m=140, d=0.264, bc=0.610, v=2700, zr=200, tr=800, lat=-33, az=45, w=10, inc=0, hum=0),
      dict(m=175, d=0.308, bc=0.505, v=2600, zr=100, tr=600, lat=45, az=0, w=5, inc=30, hum=0),
      # Humidité non nulle : `_to_si` la passe telle quelle des deux côtés, mais rien
      # ne le vérifiait — l'écriture diffère (`si.H` contre `p.humidityPct`) et seul
      # un cas à humidité réelle peut établir qu'elles désignent bien la même valeur.
      dict(m=168, d=0.308, bc=0.462, v=2600, zr=100, tr=1000, lat=50, az=90, w=0, inc=0, hum=80),
      dict(m=168, d=0.308, bc=0.462, v=2600, zr=100, tr=1000, lat=50, az=90, w=0, inc=0, hum=100),
      dict(m=140, d=0.264, bc=0.610, v=2700, zr=200, tr=800, lat=-33, az=45, w=10, inc=0, hum=55)]),

    ("solve_trajectory · dérive @500",
     "EB.solve_trajectory(EB.ShotParameters(mass_grains={m}, caliber_in={d}, bc={bc}, "
     "muzzle_vel_fps={v}, zero_range_yd={zr}, target_range_yd={tr}, latitude_deg={lat}, "
     "azimuth_deg={az}, wind_speed_mph={w}, incline_deg={inc}, humidity_pct={hum}, enable_coriolis=true))[501].windage_m",
     "EB.solveTrajectory({{...EB.DEFAULT_SHOT, massGrains:{m}, caliberIn:{d}, bc:{bc}, "
     "muzzleVelFps:{v}, zeroRangeYd:{zr}, targetRangeYd:{tr}, latitudeDeg:{lat}, "
     "azimuthDeg:{az}, windSpeedMph:{w}, inclineDeg:{inc}, humidityPct:{hum}, enableCoriolis:true}})[500].windageM",
     [dict(m=168, d=0.308, bc=0.462, v=2600, zr=100, tr=1000, lat=50, az=90, w=0, inc=0, hum=0),
      dict(m=140, d=0.264, bc=0.610, v=2700, zr=200, tr=800, lat=-33, az=45, w=10, inc=0, hum=0)]),

    ("solve_trajectory · vitesse @500",
     "EB.solve_trajectory(EB.ShotParameters(mass_grains={m}, caliber_in={d}, bc={bc}, "
     "muzzle_vel_fps={v}, zero_range_yd={zr}, target_range_yd={tr}, humidity_pct={hum}))[501].v_total",
     "EB.solveTrajectory({{...EB.DEFAULT_SHOT, massGrains:{m}, caliberIn:{d}, bc:{bc}, "
     "muzzleVelFps:{v}, zeroRangeYd:{zr}, targetRangeYd:{tr}, humidityPct:{hum}}})[500].vTotal",
     [dict(m=168, d=0.308, bc=0.462, v=2600, zr=100, tr=1000, hum=0),
      dict(m=55, d=0.224, bc=0.255, v=3200, zr=100, tr=500, hum=0),
      dict(m=55, d=0.224, bc=0.255, v=3200, zr=100, tr=500, hum=95)]),

    # La dérive gyroscopique, isolée du reste du windage. Elle est la seule des trois
    # composantes que le tireur applique séparément — le vent change d'un coup à
    # l'autre, elle non — et c'est le seul point où le 4DOF de Hornady nous devançait
    # réellement : de la présentation, pas de la physique.
    ("solve_trajectory · dérive gyro @500",
     "EB.solve_trajectory(EB.ShotParameters(mass_grains={m}, caliber_in={d}, bc={bc}, "
     "muzzle_vel_fps={v}, zero_range_yd={zr}, target_range_yd={tr}, twist_in={tw}, "
     "wind_speed_mph={w}, humidity_pct={hum}))[501].spin_drift_m",
     "EB.solveTrajectory({{...EB.DEFAULT_SHOT, massGrains:{m}, caliberIn:{d}, bc:{bc}, "
     "muzzleVelFps:{v}, zeroRangeYd:{zr}, targetRangeYd:{tr}, twistIn:{tw}, "
     "windSpeedMph:{w}, humidityPct:{hum}}})[500].spinDriftM",
     [dict(m=175, d=0.308, bc=0.505, v=2600, zr=100, tr=1100, tw=10, w=0, hum=0),
      dict(m=175, d=0.308, bc=0.505, v=2600, zr=100, tr=1100, tw=8, w=10, hum=0),
      dict(m=140, d=0.264, bc=0.610, v=2700, zr=200, tr=1000, tw=8, w=10, hum=60)]),

    ("find_zero_angle",
     "EB.find_zero_angle(EB.ShotParameters(mass_grains={m}, caliber_in={d}, bc={bc}, "
     "muzzle_vel_fps={v}, zero_range_yd={zr}, humidity_pct={hum}))",
     "EB.findZeroAngle({{...EB.DEFAULT_SHOT, massGrains:{m}, caliberIn:{d}, bc:{bc}, "
     "muzzleVelFps:{v}, zeroRangeYd:{zr}, humidityPct:{hum}}})",
     [dict(m=168, d=0.308, bc=0.462, v=2600, zr=100, hum=0),
      dict(m=168, d=0.308, bc=0.462, v=2600, zr=600, hum=0),
      dict(m=168, d=0.308, bc=0.462, v=2600, zr=600, hum=90)]),
]


def julia_source():
    lignes = ['include("%s")' % os.path.join(SRC, "CompetitionBallistics.jl"),
              "using .CompetitionBallistics",
              "using .CompetitionBallistics.BallisticUtils",
              "using .CompetitionBallistics.Atmosphere",
              "import .CompetitionBallistics: BallisticUtils, Atmosphere, DragModels, ExteriorBallistics",
              "using .CompetitionBallistics.DragModels",
              "const BU = BallisticUtils; const Atm = Atmosphere; const DM = DragModels; const EB = ExteriorBallistics",
              "import .CompetitionBallistics.ReferenceData: COMMON_BULLETS, CARTRIDGE_MAX_PSI",
              "res = Any[]"]
    for nom, gab_jl, _, grille in CAS:
        for args in grille:
            lignes.append('push!(res, try Float64(%s) catch e; "ERREUR: " * string(e) end)'
                          % gab_jl.format(**args))
    lignes.append('print("[", join([x isa String ? "\\"" * x * "\\"" : '
                  '(isfinite(x) ? string(x) : "null") for x in res], ","), "]")')
    # Données de référence, sur une seconde ligne. JSON écrit à la main : pas de
    # dépendance JSON.jl ici, et les noms de balles ne contiennent ni guillemets
    # ni antislash — l'échappement complet serait de la défense contre du vent.
    lignes.append('println()')
    lignes.append('bullets = String[]')
    lignes.append('for k in sort(collect(keys(COMMON_BULLETS)))')
    lignes.append('    b = COMMON_BULLETS[k]')
    lignes.append('    push!(bullets, "{\\"name\\":\\"" * b.name * "\\",\\"massGr\\":" * string(b.mass_gr) * '
                  '",\\"caliberIn\\":" * string(b.caliber_in) * ",\\"bcG7\\":" * string(b.bc_g7) * '
                  '",\\"lengthIn\\":" * string(b.length_in) * ",\\"twistIn\\":" * string(b.twist_in) * "}")')
    lignes.append('end')
    lignes.append('psi = String[]')
    lignes.append('for k in sort(collect(keys(CARTRIDGE_MAX_PSI)))')
    lignes.append('    push!(psi, "\\"" * k * "\\":" * string(CARTRIDGE_MAX_PSI[k]))')
    lignes.append('end')
    lignes.append('print("{\\"bullets\\":[" * join(bullets, ",") * "],\\"psi\\":{" * join(psi, ",") * "}}")')
    return "\n".join(lignes)


def js_source(selftest=False):
    lignes = ['const CB = require("%s");' % os.path.join(SRC, "CompetitionBallistics.js"),
              "const BU = CB.BallisticUtils, Atm = CB.Atmosphere, DM = CB.DragModels, EB = CB.ExteriorBallistics;",
              "const res = [];"]
    for nom, _, gab_js, grille in CAS:
        for args in grille:
            lignes.append('try { const v = %s; res.push(Number.isFinite(v) ? v : null); }'
                          ' catch (e) { res.push("ERREUR: " + e.message); }'
                          % gab_js.format(**args))
    lignes.append('const RD = CB.ReferenceData;')
    lignes.append('const bullets = Object.keys(RD.COMMON_BULLETS).sort().map(k => {')
    lignes.append('    const b = RD.COMMON_BULLETS[k];')
    # --selftest : perturber un CB de 1 % pour vérifier que le contrôle détecte.
    lignes.append('    const pert = %s && k === "Sierra 175 MK (.308)" ? 1.01 : 1;'
                  % ('true' if selftest else 'false'))
    lignes.append('    return { name: b.name, massGr: b.massGr, caliberIn: b.caliberIn, '
                  'bcG7: b.bcG7 * pert, lengthIn: b.lengthIn, twistIn: b.twistIn };')
    lignes.append('});')
    lignes.append('const psi = {};')
    lignes.append('Object.keys(RD.CARTRIDGE_MAX_PSI).sort().forEach(k => { psi[k] = RD.CARTRIDGE_MAX_PSI[k]; });')
    lignes.append('process.stdout.write(JSON.stringify(res) + "\\n" + JSON.stringify({ bullets, psi }));')
    return "\n".join(lignes)


def executer(cmd, source, suffixe):
    chemin = os.path.join("/tmp", "parity_%d%s" % (os.getpid(), suffixe))
    with open(chemin, "w", encoding="utf-8") as fh:
        fh.write(source)
    try:
        p = subprocess.run(cmd + [chemin], capture_output=True, text=True, timeout=300)
        if p.returncode != 0:
            sys.exit("%s a échoué :\n%s" % (cmd[0], (p.stderr or p.stdout)[-2000:]))
        lignes = p.stdout.strip().split("\n")
        return json.loads(lignes[0]), json.loads(lignes[-1])
    finally:
        os.unlink(chemin)


def comparer_donnees(data_jl, data_js, ecarts):
    """Confronte COMMON_BULLETS et CARTRIDGE_MAX_PSI champ à champ.

    Retourne le nombre de valeurs comparées. Les clés doivent coïncider
    exactement : une balle présente d'un seul côté est une divergence, pas un
    « rien à comparer ».
    """
    n = 0
    kj = {b["name"]: b for b in data_jl.get("bullets", [])}
    ks = {b["name"]: b for b in data_js.get("bullets", [])}
    for nom in sorted(set(kj) | set(ks)):
        if nom not in kj or nom not in ks:
            ecarts.append(("COMMON_BULLETS", {"balle": nom}, "absent", "présent", "un seul côté"))
            continue
        for champ in ("massGr", "caliberIn", "bcG7", "lengthIn", "twistIn"):
            n += 1
            v_jl, v_js = kj[nom].get(champ), ks[nom].get(champ)
            if v_jl is None or v_js is None:
                ecarts.append(("COMMON_BULLETS." + champ, {"balle": nom}, v_jl, v_js, "champ manquant"))
                continue
            if abs(v_jl - v_js) / max(abs(v_jl), abs(v_js), 1e-300) > TOL:
                ecarts.append(("COMMON_BULLETS." + champ, {"balle": nom}, v_jl, v_js, "écart de donnée"))
    jp, jsp = data_jl.get("psi", {}), data_js.get("psi", {})
    for nom in sorted(set(jp) | set(jsp)):
        n += 1
        if nom not in jp or nom not in jsp:
            ecarts.append(("CARTRIDGE_MAX_PSI", {"cartouche": nom}, "absent", "présent", "un seul côté"))
            continue
        if abs(jp[nom] - jsp[nom]) / max(abs(jp[nom]), abs(jsp[nom]), 1e-300) > TOL:
            ecarts.append(("CARTRIDGE_MAX_PSI", {"cartouche": nom}, jp[nom], jsp[nom], "écart de donnée"))
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--selftest", action="store_true",
                    help="perturbe un CB côté JS et vérifie que le contrôle détecte")
    args = ap.parse_args()

    jl, data_jl = executer(["julia"], julia_source(), ".jl")
    js, data_js = executer(["node"], js_source(selftest=args.selftest), ".js")

    total = sum(len(g) for _, _, _, g in CAS)
    if not (len(jl) == len(js) == total):
        sys.exit("comptes incohérents : %d cas déclarés, %d Julia, %d JS — "
                 "un gabarit ne s'est pas développé" % (total, len(jl), len(js)))

    ecarts, i = [], 0
    for nom, _, _, grille in CAS:
        for a in grille:
            v_jl, v_js = jl[i], js[i]
            i += 1
            if isinstance(v_jl, str) or isinstance(v_js, str):
                ecarts.append((nom, a, v_jl, v_js, "erreur d'exécution"))
                continue
            if v_jl is None and v_js is None:
                continue
            if v_jl is None or v_js is None:
                ecarts.append((nom, a, v_jl, v_js, "un côté non fini"))
                continue
            denom = max(abs(v_jl), abs(v_js), 1e-300)
            rel = abs(v_jl - v_js) / denom
            if rel > TOL:
                ecarts.append((nom, a, v_jl, v_js, "écart relatif %.3g" % rel))
            elif args.verbose:
                print("   ok  %-26s %s -> %.10g" % (nom, a, v_jl))

    n_data = comparer_donnees(data_jl, data_js, ecarts)

    if args.selftest:
        cible = [e for e in ecarts if e[0] == "COMMON_BULLETS.bcG7"
                 and e[1].get("balle") == "Sierra 175 MK (.308)"]
        if len(cible) == 1 and len(ecarts) == 1:
            print("Autotest : CB perturbé de 1 % → divergence détectée sur la bonne balle. OK.")
            return 0
        sys.exit("Autotest raté : %d divergence(s) au lieu d'une seule sur Sierra 175 MK (.308).\n%s"
                 % (len(ecarts), ecarts))

    print("Parité Julia ↔ JS : %d cas sur %d fonctions, %d données de référence."
          % (total, len(CAS), n_data))
    if not ecarts:
        print("  OK — les deux implémentations rendent les mêmes valeurs.")
        return 0
    print("  /!\\ %d divergence(s) :" % len(ecarts))
    for nom, a, v_jl, v_js, why in ecarts:
        print("     %-26s %s" % (nom, a))
        print("        Julia %-22s JS %-22s (%s)" % (v_jl, v_js, why))
    print("\n  -> La référence est le Julia, mais c'est le JS que le site SERT :")
    print("     établir lequel a raison AVANT de corriger. En 2026-08-13 c'était")
    print("     le Julia qui était faux ; le 2026-08-25, une donnée (CB Sierra)")
    print("     était fausse des deux côtés — la parité n'y pouvait rien, la")
    print("     référence tierce a tranché.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
