#!/usr/bin/env python3
"""Confronte NOTRE solveur 3-DOF à une référence tierce gelée.

    python3 scripts/check_ballistics_reference.py             # sort en 1 si divergence
    python3 scripts/check_ballistics_reference.py --verbose   # affiche chaque cas
    python3 scripts/check_ballistics_reference.py --selftest  # prouve qu'il sait crier

CE QUE CE CONTRÔLE AJOUTE À `check_js_parity.py`. La parité compare le Julia au JS :
deux rédactions de la même physique, par la même main. Elle voit une divergence, jamais
une erreur COMMUNE. Le 2026-08-22, les deux portaient la même table de traînée G7
fausse au-dessus de Mach 1,03 : la parité était parfaite, et les deux faux. Ici on
compare à une implémentation écrite par d'autres, à partir des mêmes tables publiées —
ce que la parité ne peut pas faire par construction.

CE QUE CE CONTRÔLE N'EST PAS. Un oracle. `py-ballisticcalc` est un deuxième avis. Quand
il diverge, la question « lequel a raison » reste ENTIÈRE et s'arbitre contre les tables
publiées (Litz, McCoy) ou la DOPE mesurée — jamais en s'alignant sur le tiers parce
qu'il est tiers.

LA RÉFÉRENCE EST GELÉE, pas recalculée : `data/ballistics/pointmass_reference.json`,
produit délibérément par `scripts/gen_ballistics_reference.py`. Le contrôle nocturne
n'a donc besoin ni de la bibliothèque tierce, ni du réseau, et une mise à jour de
py-ballisticcalc ne peut pas faire crier le garde-fou un matin sans que rien n'ait bougé
chez nous.

LES TOLÉRANCES SONT PHYSIQUES, pas numériques. Contrairement à la parité (1e-9 : mêmes
formules, mêmes flottants), deux implémentations indépendantes diffèrent légitimement —
ordre des opérations, pas d'intégration, interpolation de la table. On tolère donc un
écart de l'ordre de ce que le tireur ne peut pas mesurer, et on double la tolérance
sous Mach 1,2 : en transsonique, la mise à l'échelle par facteur de forme d'un modèle
G7 diverge honnêtement d'une implémentation à l'autre.
"""
import argparse
import json
import math
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JS = os.path.join(ROOT, "src", "CompetitionBallistics.js")
REFERENCE = os.path.join(ROOT, "data", "ballistics", "pointmass_reference.json")

# (clé, libellé, unité, tolérance relative, plancher absolu)
#
# Le plancher absolu existe parce qu'une tolérance purement relative n'a pas de sens
# près de zéro : à 100 yd la chute vaut 0 par construction (c'est la distance de
# zérotage) et le moindre écart y serait infiniment grand en relatif.
GRANDEURS = [
    ("drop_in",    "chute",   "in", 0.010, 0.30),
    ("windage_in", "dérive",  "in", 0.020, 0.20),
    ("vel_fps",    "vitesse", "fps", 0.005, 1.00),
    ("tof_s",      "ToF",     "s",   0.005, 0.002),
]

MACH_TRANSSONIQUE_FPS = 1400.0   # ~Mach 1,25 en atmosphère standard


def source_js(cas_liste, perturbation):
    """Fabrique le pilote node : notre solveur, sur la grille de la référence.

    ON INTERPOLE, on ne prend pas « le point le plus proche ». Le pas d'intégration
    vaut 0,5 ms, soit ~0,4 m à 800 m/s : prendre le premier point au-delà de la portée
    voulue introduirait jusqu'à 40 cm de portée en trop, et donc un écart de chute qui
    ne vient pas de la physique. C'est exactement le piège qui avait fait comparer deux
    instants différents dans `check_js_parity.py`.
    """
    lignes = [
        'const CB = require(%s);' % json.dumps(JS),
        'const EB = CB.ExteriorBallistics, BU = CB.BallisticUtils;',
        'const CAS = %s;' % json.dumps(cas_liste),
        'const PERTURBATION = %s;' % json.dumps(perturbation),
        '''
const YD = 0.9144, IN = 0.0254;
const sortie = [];
for (const c of CAS) {
  const p = {
    massGrains: c.massGrains, caliberIn: c.caliberIn, bc: c.bc * PERTURBATION,
    dragModel: c.dragModel, muzzleVelFps: c.muzzleVelFps,
    bulletLengthIn: c.bulletLengthIn, sightHeightIn: c.sightHeightIn,
    zeroRangeYd: c.zeroRangeYd, targetRangeYd: c.maxRangeYd,
    tempF: c.tempF, pressureInhg: c.pressureInhg, humidityPct: c.humidityPct,
    pressureIsSeaLevel: false, altitudeFt: 0,
    windSpeedMph: c.windSpeedMph, windAngleDeg: c.windAngleDeg,
    twistIn: c.twistIn, twistDirection: 1,
    latitudeDeg: c.latitudeDeg, azimuthDeg: c.azimuthDeg,
    enableCoriolis: c.enableCoriolis, enableSpinDrift: c.enableSpinDrift,
    inclineDeg: 0,
  };
  let traj;
  try { traj = EB.solveTrajectory(p); }
  catch (e) { sortie.push({ id: c.id, erreur: e.message }); continue; }

  const lignes = [];
  for (const yd of c.portees) {
    const cible = yd * YD;
    let i = 1;
    while (i < traj.length && traj[i].rangeM < cible) i++;
    if (i >= traj.length) { lignes.push({ yd, absent: true }); continue; }
    const a = traj[i - 1], b = traj[i];
    const f = (cible - a.rangeM) / (b.rangeM - a.rangeM);
    const lerp = (x, y) => x + f * (y - x);
    lignes.push({
      yd,
      drop_in:    lerp(a.dropM, b.dropM) / IN,
      windage_in: lerp(a.windageM, b.windageM) / IN,
      vel_fps:    lerp(a.vTotal, b.vTotal) / 0.3048,
      tof_s:      lerp(a.time, b.time),
    });
  }
  sortie.push({ id: c.id, lignes });
}
process.stdout.write(JSON.stringify(sortie));
''']
    return "\n".join(lignes)


def executer_js(source):
    chemin = os.path.join("/tmp", "ballref_%d.js" % os.getpid())
    with open(chemin, "w", encoding="utf-8") as fh:
        fh.write(source)
    try:
        p = subprocess.run(["node", chemin], capture_output=True, text=True, timeout=600)
        if p.returncode != 0:
            sys.exit("node a échoué :\n%s" % (p.stderr or p.stdout)[-2000:])
        return json.loads(p.stdout)
    finally:
        os.unlink(chemin)


def comparer(reference, perturbation=1.0):
    """Rend (nb_points, liste d'écarts). Un écart = (id, yd, libellé, ref, nous, pourquoi)."""
    defauts = dict(tempF=59.0, pressureInhg=29.92, humidityPct=0.0,
                   windSpeedMph=0.0, windAngleDeg=90.0, twistIn=10.0,
                   latitudeDeg=45.0, azimuthDeg=0.0,
                   enableSpinDrift=False, enableCoriolis=False)
    cas_js = []
    for c in reference["cas"]:
        d = dict(defauts)
        d.update({k: v for k, v in c.items() if k != "lignes"})
        d["portees"] = [l["yd"] for l in c["lignes"]]
        cas_js.append(d)

    nous = {o["id"]: o for o in executer_js(source_js(cas_js, perturbation))}

    ecarts, n = [], 0
    for c in reference["cas"]:
        o = nous.get(c["id"])
        if o is None:
            ecarts.append((c["id"], None, "cas", None, None, "absent de notre sortie"))
            continue
        if "erreur" in o:
            ecarts.append((c["id"], None, "cas", None, None, "notre solveur a levé : " + o["erreur"]))
            continue
        par_yd = {l["yd"]: l for l in o["lignes"]}
        for ref in c["lignes"]:
            mien = par_yd.get(ref["yd"])
            if mien is None or mien.get("absent"):
                ecarts.append((c["id"], ref["yd"], "portée", None, None,
                               "notre trajectoire n'atteint pas cette portée"))
                continue
            n += 1
            # Le transsonique est le seul régime où deux implémentations honnêtes
            # divergent visiblement : on y double la tolérance plutôt que d'y crier.
            mou = 2.0 if ref["vel_fps"] < MACH_TRANSSONIQUE_FPS else 1.0
            for cle, libelle, unite, tol_rel, plancher in GRANDEURS:
                a, b = ref[cle], mien[cle]
                if not (math.isfinite(a) and math.isfinite(b)):
                    ecarts.append((c["id"], ref["yd"], libelle, a, b, "valeur non finie"))
                    continue
                limite = max(abs(a) * tol_rel, plancher) * mou
                if abs(a - b) > limite:
                    ecarts.append((c["id"], ref["yd"], libelle, a, b,
                                   "écart %.3g %s > %.3g %s" % (abs(a - b), unite, limite, unite)))
    return n, ecarts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--selftest", action="store_true",
                    help="perturbe le CB de 1 %% et exige que le contrôle le voie")
    args = ap.parse_args()

    if not os.path.exists(REFERENCE):
        print("Référence absente : %s" % os.path.relpath(REFERENCE, ROOT))
        print("  La produire (délibérément) avec :")
        print("    gemini-env/bin/python scripts/gen_ballistics_reference.py")
        return 2
    with open(REFERENCE, encoding="utf-8") as fh:
        reference = json.load(fh)

    # L'AUTOTEST D'ABORD, ET IL EST BLOQUANT. Un contrôle muet ne signale rien et
    # rassure à tort — c'est ce qui est arrivé au détecteur de dérive de contenu
    # jusqu'au 2026-07-29. On perturbe le coefficient balistique de 1 %, ce qui déplace
    # la chute à 1000 yd de plusieurs pouces : si le contrôle ne le voit pas, il ne
    # verra rien, et son silence ne vaut rien.
    if args.selftest:
        n, ecarts = comparer(reference, perturbation=1.01)
        if not ecarts:
            print("AUTOTEST ÉCHOUÉ : un CB perturbé de 1 %% n'a produit aucun écart "
                  "sur %d points. Le contrôle est muet." % n)
            return 1
        print("Autotest : CB perturbé de 1 %% → %d écart(s) détecté(s) sur %d points. OK."
              % (len(ecarts), n))
        return 0

    prov = reference.get("_provenance", {})
    n, ecarts = comparer(reference)
    print("Référence point-masse (%s, moteur %s, gelée le %s) : %d cas, %d points."
          % (prov.get("bibliotheque", "?"), prov.get("moteur", "?"),
             prov.get("date", "?"), len(reference["cas"]), n))
    if not ecarts:
        print("  OK — notre solveur et la référence tierce concordent.")
        return 0

    # Regroupé par cas : trois cents lignes d'écarts individuels ne se lisent pas, et
    # une alerte illisible s'apprend à s'ignorer.
    par_cas = {}
    for e in ecarts:
        par_cas.setdefault(e[0], []).append(e)
    print("  /!\\ %d divergence(s) sur %d cas :" % (len(ecarts), len(par_cas)))
    for cid, liste in par_cas.items():
        desc = next((c.get("desc", cid) for c in reference["cas"] if c["id"] == cid), cid)
        print("     %s — %s" % (cid, desc))
        montrees = liste if args.verbose else liste[:3]
        for _, yd, libelle, a, b, why in montrees:
            if a is None:
                print("        %s : %s" % (libelle, why))
            else:
                print("        %4s yd  %-8s référence %12.3f   nous %12.3f   (%s)"
                      % (yd, libelle, a, b, why))
        if len(liste) > len(montrees):
            print("        … %d autres (--verbose pour tout voir)" % (len(liste) - len(montrees)))
    print("\n  -> La référence n'est PAS un oracle : établir lequel des deux a raison")
    print("     contre les tables publiées (Litz, McCoy) ou une DOPE mesurée AVANT")
    print("     de corriger quoi que ce soit. Et toute correction va dans les DEUX")
    print("     sources — CompetitionBallistics.js ET CompetitionBallistics.jl.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
