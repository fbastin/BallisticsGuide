#!/usr/bin/env python3
"""Fabrique la RÉFÉRENCE POINT-MASSE : une grille de trajectoires calculées par une
implémentation TIERCE, gelée dans un fichier que le contrôle nocturne relit.

    gemini-env/bin/python scripts/gen_ballistics_reference.py          # écrit la référence
    gemini-env/bin/python scripts/gen_ballistics_reference.py --stdout # sans écrire

POURQUOI UNE RÉFÉRENCE TIERCE. `check_js_parity.py` confronte le Julia au JS : deux
rédactions de la MÊME physique, écrites par la même main. Il attrape une divergence,
jamais une erreur commune. Le 2026-08-22, les deux portaient exactement la même table
de traînée G7 fausse au-dessus de Mach 1,03 — la parité était parfaite et le résultat
faux des deux côtés. Il fallait un tiers.

POURQUOI PAS UN SITE. La première idée était d'interroger JBM Ballistics. Il a retiré
ses calculateurs (constat déjà porté au journal de maintenance le 2026-08-14), et le
calculateur de Londero Sports est derrière un pare-feu qui refuse tout accès
automatisé. Une référence hébergée chez autrui n'est pas une fixture : c'est un
emprunt, révocable sans préavis. py-ballisticcalc s'installe, s'exécute hors ligne, et
son code se LIT — un désaccord devient diagnosticable au lieu d'être une énigme.

CE QUI EST GELÉ, ET POURQUOI. Le contrôle nocturne ne relance pas la bibliothèque
tierce : il relit ce fichier. Deux raisons. D'abord une mise à jour de py-ballisticcalc
ne doit pas pouvoir faire crier notre garde-fou un matin sans que rien n'ait changé
chez nous. Ensuite une référence qu'on régénère à chaque exécution finit par suivre le
code qu'elle est censée surveiller. La régénération est donc DÉLIBÉRÉE, et le fichier
porte sa provenance : version de la bibliothèque, moteur d'intégration, date.

CE N'EST PAS LA VÉRITÉ. py-ballisticcalc est un deuxième avis, pas un oracle. La
vérité reste les tables publiées (Litz, McCoy) et la DOPE mesurée. Quand les deux
divergent, il faut ARBITRER — pas s'aligner.
"""
import argparse
import datetime
import json
import os
import sys

try:
    from py_ballisticcalc import (Ammo, Atmo, Calculator, DragModel, Distance, Pressure,
                                  Shot, TableG1, TableG7, Temperature, Velocity, Weapon,
                                  Wind, logger)
except ImportError:
    sys.exit("py-ballisticcalc absent. Ce script tourne dans le venv :\n"
             "  gemini-env/bin/pip install py-ballisticcalc\n"
             "  gemini-env/bin/python scripts/gen_ballistics_reference.py")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SORTIE = os.path.join(ROOT, "data", "ballistics", "pointmass_reference.json")

TABLES = {"G1": TableG1, "G7": TableG7}

# --------------------------------------------------------------------------------
# LA GRILLE.
#
# Chaque cas est décrit dans NOS conventions (celles de `ExteriorBallistics`) ; la
# traduction vers py-ballisticcalc est faite plus bas, en un seul endroit. C'est
# volontaire : une convention traduite à deux endroits finit par l'être de deux façons.
#
# Les cas « nu » (ni vent, ni Coriolis, ni dérive gyroscopique) viennent d'abord et
# comptent le plus : ils isolent le cœur — table de traînée, atmosphère, angle de
# zérotage, intégration. Les effets sont ensuite ajoutés UN PAR UN, jamais ensemble :
# deux effets actifs à la fois rendent un écart impossible à attribuer.
# --------------------------------------------------------------------------------
GRILLE = [
    # -- Le cœur, atmosphère standard ICAO --------------------------------------
    dict(id="308-175-g7-nu", desc=".308 Win 175 gr G7 0,243 @ 2600 fps — nu",
         massGrains=175, caliberIn=0.308, bc=0.243, dragModel="G7", muzzleVelFps=2600,
         bulletLengthIn=1.240, sightHeightIn=1.5, zeroRangeYd=100, maxRangeYd=1100),
    dict(id="65cm-140-g7-nu", desc="6.5 Creedmoor 140 gr G7 0,311 @ 2700 fps — nu",
         massGrains=140, caliberIn=0.264, bc=0.311, dragModel="G7", muzzleVelFps=2700,
         bulletLengthIn=1.340, sightHeightIn=1.5, zeroRangeYd=100, maxRangeYd=1200),
    dict(id="338lm-300-g7-nu", desc=".338 Lapua Magnum 300 gr G7 0,419 @ 2750 fps — nu",
         massGrains=300, caliberIn=0.338, bc=0.419, dragModel="G7", muzzleVelFps=2750,
         bulletLengthIn=1.820, sightHeightIn=2.2, zeroRangeYd=100, maxRangeYd=1500),
    # G1 : le modèle des catalogues, et celui que saisit le tireur qui n'a que la
    # valeur du fabricant. Sa branche supersonique diverge davantage encore que G7.
    dict(id="308-168-g1-nu", desc=".308 Win 168 gr G1 0,462 @ 2650 fps — nu",
         massGrains=168, caliberIn=0.308, bc=0.462, dragModel="G1", muzzleVelFps=2650,
         bulletLengthIn=1.215, sightHeightIn=1.5, zeroRangeYd=100, maxRangeYd=1000),
    # Balle légère et rapide : Mach 3 au départ, donc le haut de la table de traînée,
    # là où aucun des cas précédents ne va.
    dict(id="223-55-g1-nu", desc=".223 Rem 55 gr G1 0,255 @ 3240 fps — nu",
         massGrains=55, caliberIn=0.224, bc=0.255, dragModel="G1", muzzleVelFps=3240,
         bulletLengthIn=0.760, sightHeightIn=2.6, zeroRangeYd=100, maxRangeYd=700),

    # -- L'atmosphère ------------------------------------------------------------
    # Chaud et haut : densité ~0,85 de la standard. Une erreur d'échelle sur la
    # pression ou la température sort ici, pas au niveau de la mer.
    dict(id="65cm-140-g7-chaud-haut", desc="6.5 CM 140 gr — 30 °C, 850 hPa station, 50 % HR",
         massGrains=140, caliberIn=0.264, bc=0.311, dragModel="G7", muzzleVelFps=2700,
         bulletLengthIn=1.340, sightHeightIn=1.5, zeroRangeYd=100, maxRangeYd=1200,
         tempF=86.0, pressureInhg=25.10, humidityPct=50),
    dict(id="65cm-140-g7-froid", desc="6.5 CM 140 gr — −15 °C, 1030 hPa, 0 % HR",
         massGrains=140, caliberIn=0.264, bc=0.311, dragModel="G7", muzzleVelFps=2700,
         bulletLengthIn=1.340, sightHeightIn=1.5, zeroRangeYd=100, maxRangeYd=1200,
         tempF=5.0, pressureInhg=30.42, humidityPct=0),

    # -- Le vent -----------------------------------------------------------------
    # 90° : le travers pur. 45° : la composante AXIALE, celle où les conventions se
    # séparent — chez nous 0° est un vent debout, chez py-ballisticcalc 0° est un vent
    # arrière. Un cas à 90° seul n'aurait jamais montré l'inversion.
    dict(id="65cm-140-g7-vent90", desc="6.5 CM 140 gr — vent 10 mph à 90°",
         massGrains=140, caliberIn=0.264, bc=0.311, dragModel="G7", muzzleVelFps=2700,
         bulletLengthIn=1.340, sightHeightIn=1.5, zeroRangeYd=100, maxRangeYd=1200,
         windSpeedMph=10, windAngleDeg=90),
    dict(id="65cm-140-g7-vent45", desc="6.5 CM 140 gr — vent 10 mph à 45°",
         massGrains=140, caliberIn=0.264, bc=0.311, dragModel="G7", muzzleVelFps=2700,
         bulletLengthIn=1.340, sightHeightIn=1.5, zeroRangeYd=100, maxRangeYd=1200,
         windSpeedMph=10, windAngleDeg=45),

    # -- La dérive gyroscopique ---------------------------------------------------
    # Même formule des deux côtés (Litz, 1,25·(Sg+1,2)·t^1,83) : ce cas ne teste pas la
    # formule mais le Sg de Miller qui l'alimente, et le temps de vol qui l'alimente
    # aussi. Un Sg faux se voit ici et nulle part ailleurs.
    dict(id="308-175-g7-spin", desc=".308 175 gr — dérive gyroscopique, pas 1:10\"",
         massGrains=175, caliberIn=0.308, bc=0.243, dragModel="G7", muzzleVelFps=2600,
         bulletLengthIn=1.240, sightHeightIn=1.5, zeroRangeYd=100, maxRangeYd=1100,
         twistIn=10, enableSpinDrift=True),
    dict(id="65cm-140-g7-spin", desc="6.5 CM 140 gr — dérive gyroscopique, pas 1:8\"",
         massGrains=140, caliberIn=0.264, bc=0.311, dragModel="G7", muzzleVelFps=2700,
         bulletLengthIn=1.340, sightHeightIn=1.5, zeroRangeYd=100, maxRangeYd=1200,
         twistIn=8, enableSpinDrift=True),

    # -- Coriolis -----------------------------------------------------------------
    # Est et Ouest au même endroit : l'effet Eötvös change de signe entre les deux, et
    # c'est la seule paire qui puisse établir que le signe est bon. Le Julia l'avait
    # inversé sans que rien ne le voie.
    dict(id="65cm-140-g7-coriolis-est", desc="6.5 CM 140 gr — Coriolis, 45° N, azimut 90° (Est)",
         massGrains=140, caliberIn=0.264, bc=0.311, dragModel="G7", muzzleVelFps=2700,
         bulletLengthIn=1.340, sightHeightIn=1.5, zeroRangeYd=100, maxRangeYd=1200,
         latitudeDeg=45, azimuthDeg=90, enableCoriolis=True),
    dict(id="65cm-140-g7-coriolis-ouest", desc="6.5 CM 140 gr — Coriolis, 45° N, azimut 270° (Ouest)",
         massGrains=140, caliberIn=0.264, bc=0.311, dragModel="G7", muzzleVelFps=2700,
         bulletLengthIn=1.340, sightHeightIn=1.5, zeroRangeYd=100, maxRangeYd=1200,
         latitudeDeg=45, azimuthDeg=270, enableCoriolis=True),
    # Hémisphère sud : la dérive horizontale change de côté.
    dict(id="65cm-140-g7-coriolis-sud", desc="6.5 CM 140 gr — Coriolis, 33° S, azimut 0° (Nord)",
         massGrains=140, caliberIn=0.264, bc=0.311, dragModel="G7", muzzleVelFps=2700,
         bulletLengthIn=1.340, sightHeightIn=1.5, zeroRangeYd=100, maxRangeYd=1200,
         latitudeDeg=-33, azimuthDeg=0, enableCoriolis=True),
]

DEFAUTS = dict(tempF=59.0, pressureInhg=29.92, humidityPct=0.0,
               windSpeedMph=0.0, windAngleDeg=90.0, twistIn=10.0,
               enableSpinDrift=False, enableCoriolis=False)


def distances(cas):
    """Les portées échantillonnées : de 100 yd au maximum du cas, par pas de 100."""
    return list(range(100, int(cas["maxRangeYd"]) + 1, 100))


def calculer(cas):
    """Traduit UN cas de nos conventions vers py-ballisticcalc, et l'exécute.

    Les trois traductions qui comptent, toutes vérifiées et non devinées :

    - VENT. Chez nous, `windAngleDeg` = 0 est un vent DEBOUT (`wx = −w·cos θ`) ; chez
      py-ballisticcalc, `direction_from` = 0 est un vent ARRIÈRE. Les deux vecteurs
      coïncident pour φ = 180 − θ. À 90° les deux conventions se confondent, ce qui
      rend l'inversion invisible sur un cas de travers pur — d'où le cas à 45°.
    - DÉRIVE GYROSCOPIQUE. Elle ne s'éteint pas par un drapeau chez eux : c'est
      `twist = 0` qui l'annule.
    - CORIOLIS. Il ne s'éteint pas non plus : ce sont `latitude` et `azimuth` laissées
      à None.
    """
    p = dict(DEFAUTS)
    p.update(cas)

    dm = DragModel(p["bc"], TABLES[p["dragModel"]],
                   weight=p["massGrains"], diameter=p["caliberIn"],
                   length=p["bulletLengthIn"])
    ammo = Ammo(dm, Velocity.FPS(p["muzzleVelFps"]))
    arme = Weapon(sight_height=Distance.Inch(p["sightHeightIn"]),
                  twist=Distance.Inch(p["twistIn"]) if p["enableSpinDrift"] else 0)
    atmo = Atmo(altitude=Distance.Foot(0),
                pressure=Pressure.InHg(p["pressureInhg"]),
                temperature=Temperature.Fahrenheit(p["tempF"]),
                humidity=p["humidityPct"])
    vents = []
    if p["windSpeedMph"]:
        vents.append(Wind(velocity=Velocity.MPH(p["windSpeedMph"]),
                          direction_from=180.0 - p["windAngleDeg"]))

    tir = dict(ammo=ammo, weapon=arme, atmo=atmo, winds=vents)
    if p["enableCoriolis"]:
        tir["latitude"] = p["latitudeDeg"]
        tir["azimuth"] = p["azimuthDeg"]
    shot = Shot(**tir)

    calc = Calculator(engine="rk4_engine")   # même intégrateur que le nôtre
    calc.set_weapon_zero(shot, Distance.Yard(p["zeroRangeYd"]))

    portees = distances(p)
    res = calc.fire(shot, trajectory_range=Distance.Yard(portees[-1]),
                    trajectory_step=Distance.Yard(100))

    par_yd = {}
    for r in res.trajectory:
        yd = round(r.distance >> Distance.Yard)
        par_yd[yd] = dict(
            yd=yd,
            drop_in=round(r.height >> Distance.Inch, 4),
            windage_in=round(r.windage >> Distance.Inch, 4),
            vel_fps=round(r.velocity >> Velocity.FPS, 4),
            tof_s=round(r.time, 6),
        )
    lignes = [par_yd[yd] for yd in portees if yd in par_yd]
    if len(lignes) != len(portees):
        manquantes = [yd for yd in portees if yd not in par_yd]
        raise SystemExit("cas %s : portées absentes de la sortie tierce : %s"
                         % (p["id"], manquantes))
    return lignes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stdout", action="store_true",
                    help="écrit sur la sortie standard au lieu du fichier de référence")
    args = ap.parse_args()

    logger.setLevel(40)   # la bibliothèque annonce son moteur en INFO : on se tait

    import importlib.metadata as md
    doc = {
        "_provenance": {
            "produit_par": "scripts/gen_ballistics_reference.py",
            "bibliotheque": "py-ballisticcalc " + md.version("py-ballisticcalc"),
            "moteur": "rk4_engine",
            "date": datetime.date.today().isoformat(),
            "avertissement": "Deuxième avis, pas oracle. En cas de divergence, ARBITRER "
                             "contre les tables publiées (Litz, McCoy) ou la DOPE mesurée.",
        },
        "cas": [],
    }
    for cas in GRILLE:
        entree = {k: v for k, v in cas.items()}
        entree["lignes"] = calculer(cas)
        doc["cas"].append(entree)

    texte = json.dumps(doc, indent=1, ensure_ascii=False) + "\n"
    if args.stdout:
        sys.stdout.write(texte)
        return 0
    os.makedirs(os.path.dirname(SORTIE), exist_ok=True)
    with open(SORTIE, "w", encoding="utf-8") as fh:
        fh.write(texte)
    print("Référence écrite : %s" % os.path.relpath(SORTIE, ROOT))
    print("  %d cas, %d points." % (len(doc["cas"]), sum(len(c["lignes"]) for c in doc["cas"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
