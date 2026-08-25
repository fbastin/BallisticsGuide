// ==========================================================================
// CompetitionBallistics.js
// JavaScript port of the core ballistic computation modules from
// CompetitionBallistics.jl — ready for browser or Node.js usage.
//
// Modules ported:
//   - BallisticUtils      : Unit conversions, sectional density, stability
//   - ReferenceData       : Common bullet profiles, cartridge pressure limits
//   - Atmosphere          : ICAO standard atmosphere, humidity, density altitude
//   - InteriorBallistics  : Le Duc model, pressure, burn rate, temp sensitivity
//   - DragModels          : G1/G7 BRL tabular data with interpolation
//   - ExteriorBallistics  : Full 3-DOF solver (RK4), Coriolis, spin drift, wind
//   - ParameterMetadata   : Physical explanations and limits (NEW)
// ==========================================================================

// --------------------------------------------------------------------------
// MODULE 1: BallisticUtils
// --------------------------------------------------------------------------
const BallisticUtils = (() => {

  // -- Mass --
  const grainsToKg    = (gr) => gr * 6.47989e-5;
  const kgToGrains    = (kg) => kg / 6.47989e-5;
  const grainsToGrams = (gr) => gr * 0.06480;
  const gramsToGrains = (g)  => g / 0.06480;

  // -- Velocity --
  const fpsToMs = (fps) => fps * 0.3048;
  const msToFps = (ms)  => ms / 0.3048;

  // -- Length --
  const inchesToM  = (x) => x * 0.0254;
  const mToInches  = (x) => x / 0.0254;
  const inchesToMm = (x) => x * 25.4;
  const mmToInches = (x) => x / 25.4;
  const yardsToM   = (x) => x * 0.9144;
  const mToYards   = (x) => x / 0.9144;

  // -- Temperature --
  const fahrenheitToKelvin  = (f) => (f - 32.0) * 5.0 / 9.0 + 273.15;
  const kelvinToFahrenheit  = (k) => (k - 273.15) * 9.0 / 5.0 + 32.0;
  const fahrenheitToRankine = (f) => f + 459.67;
  const celsiusToKelvin     = (c) => c + 273.15;
  const kelvinToCelsius     = (k) => k - 273.15;

  // -- Pressure --
  const inhgToPa  = (x) => x * 3386.389;
  const paToInhg  = (x) => x / 3386.389;
  const inhgToHpa = (x) => x * 33.8639;
  const hpaToInhg = (x) => x / 33.8639;

  // -- Angles --
  const moaToRad = (moa) => moa * (Math.PI / (180.0 * 60.0));
  const radToMoa = (rad) => rad * (180.0 * 60.0 / Math.PI);
  const milToRad = (mil) => mil * 0.001;
  const radToMil = (rad) => rad * 1000.0;
  const moaToMil = (moa) => moa * 0.29089;
  const milToMoa = (mil) => mil * 3.43775;

  const dropToMoa = (dropInches, rangeYards) =>
    dropInches / (rangeYards * 1.04720 / 100.0);

  const dropToMil = (dropInches, rangeYards) => {
    const rangeM = yardsToM(rangeYards);
    const dropM  = inchesToM(dropInches);
    return dropM / rangeM * 1000.0;
  };

  // -- Energy --
  const kineticEnergyJ     = (massKg, velMs) => 0.5 * massKg * velMs * velMs;
  const kineticEnergyFtlbs = (massGr, velFps) => massGr * velFps * velFps / 450436.0;

  // -- Bullet descriptors --
  const sectionalDensity = (massGr, caliberIn) =>
    (massGr / 7000.0) / (caliberIn * caliberIn);

  const formFactor = (massGr, caliberIn, bc) =>
    sectionalDensity(massGr, caliberIn) / bc;

  const millerStability = ({
    massGr, caliberIn, bulletLengthIn, twistIn,
    muzzleVelFps = 2800.0, tempF = 59.0, pressureInhg = 29.92
  }) => {
    const d  = caliberIn;
    const l  = bulletLengthIn / d;
    const tw = twistIn / d;
    const tR = fahrenheitToRankine(tempF);
    let sg = 30.0 * massGr / (tw * tw * d * d * d * l * (1.0 + l * l));
    // Correction de vitesse de Miller : racine cubique du rapport à 2800 fps,
    // et non rapport direct. L'écart est faible près de 2800 fps mais atteint
    // 27 % à 4000 fps.
    sg *= Math.cbrt(muzzleVelFps / 2800.0);
    // Correction atmosphérique : Sg varie comme 1/rho, et rho ~ P/T.
    // Air chaud ou raréfié (altitude) = Sg plus élevé.
    sg *= (tR / 518.67) * (29.92 / pressureInhg);
    return sg;
  };

  // Formule amendée de Courtney & Miller (Precision Shooting, janvier 2012) pour
  // les balles à pointe plastique. Le plastique est ~10 fois moins dense que le
  // métal : il allonge la balle sans rien ajouter aux moments d'inertie. Seul le
  // terme (1 + l²), qui vient des moments d'inertie, prend la longueur métallique.
  // L'autre l vient du centre de pression, propriété aérodynamique de forme, et
  // garde la longueur totale. Appliquer le Miller ordinaire à ces balles
  // SOUS-ESTIME leur stabilité.
  const millerStabilityTipped = ({
    massGr, caliberIn, bulletLengthIn, metalLengthIn, twistIn,
    muzzleVelFps = 2800.0, tempF = 59.0, pressureInhg = 29.92
  }) => {
    const d  = caliberIn;
    const l  = bulletLengthIn / d;
    const lm = metalLengthIn / d;
    const tw = twistIn / d;
    let sg = 30.0 * massGr / (tw * tw * d * d * d * l * (1.0 + lm * lm));
    sg *= Math.cbrt(muzzleVelFps / 2800.0);
    sg *= (fahrenheitToRankine(tempF) / 518.67) * (29.92 / pressureInhg);
    return sg;
  };

  // --- Passage masse <-> longueur ------------------------------------------
  //
  // Miller exige une LONGUEUR ; le tireur connaît un POIDS. Pour une balle de
  // forme et de densité moyenne données, m ~ rho * d^2 * L, donc la longueur en
  // calibres l = L/d varie comme m/d^3. Le coefficient ci-dessous a été calibré
  // le 2026-07-27 sur neuf balles match chemisées dont la longueur est publiée
  // par le fabricant (Sierra 55/77/155/175, Berger 105/140/180/185/300), du .224
  // au .338 : écart maximal 6 % sur la longueur.
  //
  // ATTENTION : c'est une constante de FAMILLE, pas une loi. 6 % d'écart sur la
  // longueur valent +20/-16 % sur Sg. Et la relation se dégrade hors du match
  // chemisé — un monobloc cuivre mesuré (Barnes TSX 168 gr .308, 1,318") est 8 %
  // plus long que la relation ne le prédit, soit 20 % de Sg en moins. Dès que la
  // longueur publiée est disponible, elle prime sur toute estimation.
  const BULLET_SHAPE_K = 6.9024e-4;

  const bulletLengthFromWeight = (massGr, caliberIn, k = BULLET_SHAPE_K) =>
    k * massGr / (caliberIn * caliberIn);          // = (k*m/d^3) * d

  const bulletWeightFromLength = (bulletLengthIn, caliberIn, k = BULLET_SHAPE_K) =>
    bulletLengthIn * caliberIn * caliberIn / k;

  // --- Sens inverse : quel pas faut-il ? -----------------------------------
  //
  // PIÈGE : ce n'est pas l'inverse littéral de millerStability. Miller applique
  // au PAS des exposants deux fois plus petits qu'au facteur de stabilité —
  // « multiply the calculated s by the 1/3 root of (v/2800) and calculated t by
  // the 1/6 root », et fTP « multiplies the calculated s and t², and its square
  // root multiplies t itself ». Inverser naïvement le code de Sg donne un pas faux.
  const millerTwistForSg = ({
    massGr, caliberIn, bulletLengthIn, sgTarget = 2.0,
    muzzleVelFps = 2800.0, tempF = 59.0, pressureInhg = 29.92
  }) => {
    const d = caliberIn;
    const l = bulletLengthIn / d;
    const twSq = 30.0 * massGr / (sgTarget * d * d * d * l * (1.0 + l * l));
    const fTP = (fahrenheitToRankine(tempF) / 518.67) * (29.92 / pressureInhg);
    const tw = Math.sqrt(twSq) * Math.pow(muzzleVelFps / 2800.0, 1.0 / 6.0) * Math.sqrt(fTP);
    return tw * d;                                  // pas en pouces par tour
  };

  // Poids maximal admissible pour un pas donné, via la relation masse <-> longueur.
  // En y substituant m = l*d^3/k, d et l se simplifient de Miller sauf dans
  // (1 + l²) : Sg = 30 * fv * fTP / (k * t² * (1 + l²)). Renvoie null si le pas
  // ne peut atteindre la cible pour aucune longueur (balle plus courte impossible).
  const maxBulletForTwist = ({
    caliberIn, twistIn, sgTarget = 2.0,
    muzzleVelFps = 2800.0, tempF = 59.0, pressureInhg = 29.92, k = BULLET_SHAPE_K
  }) => {
    const d = caliberIn;
    const tw = twistIn / d;
    const fv = Math.cbrt(muzzleVelFps / 2800.0);
    const fTP = (fahrenheitToRankine(tempF) / 518.67) * (29.92 / pressureInhg);
    const lSq = 30.0 * fv * fTP / (k * tw * tw * sgTarget) - 1.0;
    if (lSq <= 0) return null;
    const l = Math.sqrt(lSq);
    return { bulletLengthIn: l * d, massGr: l * d * d * d / k };
  };

  const greenhillTwist = (caliberIn, bulletLengthIn, C = 150.0) =>
    C * caliberIn * caliberIn / bulletLengthIn;

  const loadConsistencyScore = (esFps, sdFps) => {
    let score = 100.0 - 5.0 * sdFps - 2.0 * esFps;
    return Math.max(0.0, Math.min(score, 100.0));
  };

  const bcFromTwoChronographs = (v1Fps, v2Fps, distFt, calIn, massGr, model = "G7", tempF = 59.0, presInhg = 29.92, rho0 = 1.225) => {
    const v1 = fpsToMs(v1Fps);
    const v2 = fpsToMs(v2Fps);
    const dist = distFt * 0.3048;
    const T = fahrenheitToKelvin(tempF);
    
    // Inline simplified atmosphere to avoid cyclic dependencies
    const P_Pa = inhgToPa(presInhg);
    const R_air = 287.0528;
    const rho = P_Pa / (R_air * T);
    const aS = Math.sqrt(1.4 * R_air * T);

    const vAvg = (v1 + v2) / 2.0;
    const Ma = vAvg / aS;
    
    // We assume DragModels is available when this function is called
    // If called directly before DragModels is initialized, it will throw.
    const cdStd = CompetitionBallistics.DragModels.dragCoefficient(model, Ma);

    // Décroissance de vitesse sous une traînée en v^2, intégrée sur la DISTANCE :
    //     dv/dx = -k v   =>   ln(v1/v2) = k x
    // La forme en (1/v2 - 1/v1) qui figurait ici est celle de l'intégration sur le
    // TEMPS ; appliquée à une distance, elle rendait un BC ~10^5 fois trop grand.
    //
    // k est la constante de retardation du solveur, k = (rho/2) Cd C0 / BC, où C0
    // convertit la densité sectionnelle (lb/in^2, définie sur d^2) vers le SI :
    //     C0 = 0.0254^2 / 0.45359237
    // C'est le même appariement SD <-> table de traînée que dans solveTrajectory.
    const C0 = 0.0254 * 0.0254 / 0.45359237;
    const bc = (rho * cdStd * C0 * dist) / (2.0 * Math.log(v1 / v2));
    // Pas de renormalisation par rho0 : rho est déjà celui de la mesure, donc le BC
    // obtenu est le vrai BC du projectile, indépendant des conditions.
    return Math.abs(bc);
  };

  return {
    grainsToKg, kgToGrains, grainsToGrams, gramsToGrains,
    fpsToMs, msToFps,
    inchesToM, mToInches, inchesToMm, mmToInches, yardsToM, mToYards,
    fahrenheitToKelvin, kelvinToFahrenheit, fahrenheitToRankine,
    celsiusToKelvin, kelvinToCelsius,
    inhgToPa, paToInhg, inhgToHpa, hpaToInhg,
    moaToRad, radToMoa, milToRad, radToMil, moaToMil, milToMoa,
    dropToMoa, dropToMil,
    kinetic_energy_J: kineticEnergyJ, kinetic_energy_ftlbs: kineticEnergyFtlbs, // matching Julia names
    kineticEnergyJ, kineticEnergyFtlbs,
    sectionalDensity, formFactor,
    millerStability, millerStabilityTipped, greenhillTwist,
    millerTwistForSg, maxBulletForTwist,
    bulletLengthFromWeight, bulletWeightFromLength, BULLET_SHAPE_K,
    loadConsistencyScore, bcFromTwoChronographs,
  };
})();


// --------------------------------------------------------------------------
// MODULE 1.1: ReferenceData
// --------------------------------------------------------------------------
const ReferenceData = (() => {

  const COMMON_BULLETS = {
    "Berger 105 Hybrid (6mm)":       { name: "Berger 105 Hybrid (6mm)",       massGr: 105.0, caliberIn: 0.243, bcG7: 0.275, lengthIn: 1.10, twistIn: 8.0 },
    "Berger 140 Hybrid (6.5mm)":     { name: "Berger 140 Hybrid (6.5mm)",     massGr: 140.0, caliberIn: 0.264, bcG7: 0.311, lengthIn: 1.34, twistIn: 8.0 },
    "Hornady 147 ELD-M (6.5mm)":     { name: "Hornady 147 ELD-M (6.5mm)",     massGr: 147.0, caliberIn: 0.264, bcG7: 0.351, lengthIn: 1.42, twistIn: 7.5 },
    "Berger 180 Hybrid (7mm)":       { name: "Berger 180 Hybrid (7mm)",       massGr: 180.0, caliberIn: 0.284, bcG7: 0.350, lengthIn: 1.50, twistIn: 8.5 },
    "Sierra 175 MK (.308)":          { name: "Sierra 175 MK (.308)",          massGr: 175.0, caliberIn: 0.308, bcG7: 0.243, lengthIn: 1.24, twistIn: 10.0 },
    "Berger 185 Juggernaut (.308)":  { name: "Berger 185 Juggernaut (.308)",  massGr: 185.0, caliberIn: 0.308, bcG7: 0.283, lengthIn: 1.30, twistIn: 10.0 },
    "Hornady 178 ELD-M (.308)":      { name: "Hornady 178 ELD-M (.308)",      massGr: 178.0, caliberIn: 0.308, bcG7: 0.274, lengthIn: 1.31, twistIn: 10.0 },
    "Berger 300 Hybrid (.338)":      { name: "Berger 300 Hybrid (.338)",      massGr: 300.0, caliberIn: 0.338, bcG7: 0.419, lengthIn: 1.82, twistIn: 9.4 },
  };

  const CARTRIDGE_MAX_PSI = {
    ".223 Remington":        55000.0,
    "6.5 Creedmoor":         62000.0,
    ".308 Winchester":       62000.0,
    ".300 Winchester Magnum": 64000.0,
    ".338 Lapua Magnum":     60916.0,
  };

  return { COMMON_BULLETS, CARTRIDGE_MAX_PSI };
})();


// --------------------------------------------------------------------------
// MODULE 1.2: ParameterMetadata (NEW)
// --------------------------------------------------------------------------
const ParameterMetadata = (() => {
  const DATA = {
    mass: {
      label: "Masse",
      help: "L'inertie dépend de la masse. À calibre égal, une balle lourde garde mieux sa vitesse.",
      physics: "Une masse plus élevée augmente la densité sectionnelle et réduit la décélération pour une force de traînée donnée (F=ma).",
      min: 1, max: 1000
    },
    bc: {
      label: "Coefficient Balistique",
      help: "Plus le CB est élevé, mieux la balle fend l'air. Diminue généralement avec la masse.",
      physics: "Le CB (G7) compare la traînée du projectile à celle d'un projectile de référence profilé.",
      min: 0.01, max: 1.2
    },
    v0: {
      label: "Vitesse Initiale",
      help: "Vitesse réelle mesurée. Une balle plus légère part plus vite avec la même charge.",
      physics: "Énergie cinétique E = 1/2 mv². Une diminution de m à E constant augmente v de façon quadratique.",
      min: 100, max: 5000
    }
  };
  return { DATA };
})();


// --------------------------------------------------------------------------
// MODULE 2: Atmosphere
// --------------------------------------------------------------------------
const Atmosphere = (() => {

  const T0        = 288.15;     // K
  const P0        = 101325.0;   // Pa
  const rho0      = 1.2250;     // kg/m³
  const L         = 0.0065;     // K/m
  const g0        = 9.80665;    // m/s²
  const R_air     = 287.0528;   // J/(kg·K) specific gas constant for dry air (ISA 1976)
  const gamma_air = 1.4;

  const stdTemperature = (h) => T0 - L * h;

  const stdPressure = (h) => {
    const T = stdTemperature(h);
    return P0 * Math.pow(T / T0, g0 / (L * R_air));
  };

  const saturationVaporPressure = (tK) => {
    const tC = tK - 273.15;
    return 610.78 * Math.pow(10.0, 7.5 * tC / (tC + 237.3));
  };

  const airDensity = ({ P = P0, T = T0, H = 0.0 } = {}) => {
    const es = saturationVaporPressure(T);
    const e  = (H / 100.0) * es;
    return (P - 0.37802 * e) / (R_air * T);
  };

  const speedOfSound = (T = T0) => Math.sqrt(gamma_air * R_air * T);

  const densityRatio = ({ P = P0, T = T0, H = 0.0 } = {}) =>
    airDensity({ P, T, H }) / rho0;

  const densityAltitude = (rho) =>
    (T0 / L) * (1.0 - Math.pow(rho / rho0, 0.234969));

  const densityAltitudeFt = (rho) =>
    145442.16 * (1.0 - Math.pow(rho / rho0, 0.234969));

  return {
    T0, P0, rho0, L, g0, R_air, gamma_air,
    stdTemperature, stdPressure,
    saturationVaporPressure, airDensity,
    speedOfSound, densityRatio,
    densityAltitude, densityAltitudeFt,
  };
})();


// --------------------------------------------------------------------------
// MODULE 3: InteriorBallistics
// --------------------------------------------------------------------------
const InteriorBallistics = (() => {

  const leducVelocity = (x, { a, b }) => a * x / (b + x);

  const leducPressure = (x, { a, b, mBullet, mCharge, dBore }) => {
    const mE = mBullet + mCharge / 3.0;
    const A  = Math.PI * dBore * dBore / 4.0;
    return mE * a * a * b * x / (A * Math.pow(b + x, 3));
  };

  const peakPressure = ({ a, b, mBullet, mCharge, dBore }) => {
    const mE = mBullet + mCharge / 3.0;
    const A  = Math.PI * dBore * dBore / 4.0;
    return 4.0 * mE * a * a / (27.0 * A * b);
  };

  const effectiveMass = (mBullet, mCharge) => mBullet + mCharge / 3.0;

  const burnFraction = (xi, theta = 0.0) => xi * (1.0 + theta * xi);

  const propellantEnergy = (chargeMassKg, impetus, z) =>
    chargeMassKg * impetus * z;

  const chamberPressureClosedBomb = ({ C, f, z, Vch, rhoP, etaSp }) => {
    const Vgas = Vch - C * (1.0 - z) / rhoP - C * z * etaSp;
    return C * f * z / Vgas;
  };

  const barrelLengthCorrection = (vRef, lBarrel, lRef, exponent = 0.27) =>
    vRef * Math.pow(lBarrel / lRef, exponent);

  const muzzleVelocityTempCorrection = (v0, dT, sigmaT = 1.0) =>
    v0 + sigmaT * dT;

  return {
    leducVelocity, leducPressure, peakPressure,
    effective_mass: effectiveMass, burn_fraction: burnFraction, // matching Julia names
    effectiveMass, burnFraction, propellantEnergy,
    chamberPressureClosedBomb,
    barrelLengthCorrection, muzzleVelocityTempCorrection,
  };
})();


// --------------------------------------------------------------------------
// MODULE 4: DragModels
// --------------------------------------------------------------------------
const DragModels = (() => {

  // Fonctions de traînée étalons BRL G7 et G1 (C_D en fonction du nombre de Mach).
  //
  // ⚠️ CES VALEURS SE RECOPIENT, ELLES NE SE PROLONGENT PAS. Jusqu'au 2026-08-22 les
  // deux tables étaient exactes dans leur partie basse puis INVENTÉES au-delà — G7 à
  // partir de Mach 1,05, G1 dès Mach 0,75, où la fin du domaine subsonique reprenait
  // en miroir les valeurs du début. L'écart atteignait −33 % (G7) et −45 % (G1) à
  // Mach 2. Le solveur portait par-dessus un facteur 4/π calé pour qu'UNE .308 arrive
  // à 1000 yd sur une vitesse plausible : l'erreur de données était donc masquée à un
  // point de fonctionnement et fausse partout ailleurs.
  //
  // La parité Julia↔JS ne pouvait rien y voir : les deux copies portaient la même
  // faute. C'est `scripts/check_ballistics_reference.py`, qui confronte à une
  // implémentation tierce, qui l'a sortie — et l'arbitrage s'est fait sur McCoy,
  // « Modern Exterior Ballistics », figure 4.6 p. 56 : le BRL-1, décrit comme *very
  // low drag* et plus fin que l'étalon G7, y tient C_D ≈ 0,26 à Mach 3, quand notre
  // table en annonçait 0,140.

  const G7_RAW = [
    0.000,0.1198,  0.050,0.1197,  0.100,0.1196,
    0.150,0.1194,  0.200,0.1193,  0.250,0.1194,
    0.300,0.1194,  0.350,0.1194,  0.400,0.1193,
    0.450,0.1193,  0.500,0.1194,  0.550,0.1193,
    0.600,0.1194,  0.650,0.1197,  0.700,0.1202,
    0.725,0.1207,  0.750,0.1215,  0.775,0.1226,
    0.800,0.1242,  0.825,0.1266,  0.850,0.1306,
    0.875,0.1368,  0.900,0.1464,  0.925,0.1660,
    0.950,0.2054,  0.975,0.2993,  1.000,0.3803,
    1.025,0.4015,  1.050,0.4043,  1.075,0.4034,
    1.100,0.4014,  1.125,0.3987,  1.150,0.3955,
    1.200,0.3884,  1.250,0.3810,  1.300,0.3732,
    1.350,0.3657,  1.400,0.3580,  1.500,0.3440,
    1.550,0.3376,  1.600,0.3315,  1.650,0.3260,
    1.700,0.3209,  1.750,0.3160,  1.800,0.3117,
    1.850,0.3078,  1.900,0.3042,  1.950,0.3010,
    2.000,0.2980,  2.050,0.2951,  2.100,0.2922,
    2.150,0.2892,  2.200,0.2864,  2.250,0.2835,
    2.300,0.2807,  2.350,0.2779,  2.400,0.2752,
    2.450,0.2725,  2.500,0.2697,  2.550,0.2670,
    2.600,0.2643,  2.650,0.2615,  2.700,0.2588,
    2.750,0.2561,  2.800,0.2533,  2.850,0.2506,
    2.900,0.2479,  2.950,0.2451,  3.000,0.2424,
    3.100,0.2368,  3.200,0.2313,  3.300,0.2258,
    3.400,0.2205,  3.500,0.2154,  3.600,0.2106,
    3.700,0.2060,  3.800,0.2017,  3.900,0.1975,
    4.000,0.1935,  4.200,0.1861,  4.400,0.1793,
    4.600,0.1730,  4.800,0.1672,  5.000,0.1618,
  ];

  const G1_RAW = [
    0.000,0.2629,  0.050,0.2558,  0.100,0.2487,
    0.150,0.2413,  0.200,0.2344,  0.250,0.2278,
    0.300,0.2214,  0.350,0.2155,  0.400,0.2104,
    0.450,0.2061,  0.500,0.2032,  0.550,0.2020,
    0.600,0.2034,  0.700,0.2165,  0.725,0.2230,
    0.750,0.2313,  0.775,0.2417,  0.800,0.2546,
    0.825,0.2706,  0.850,0.2901,  0.875,0.3136,
    0.900,0.3415,  0.925,0.3734,  0.950,0.4084,
    0.975,0.4448,  1.000,0.4805,  1.025,0.5136,
    1.050,0.5427,  1.075,0.5677,  1.100,0.5883,
    1.125,0.6053,  1.150,0.6191,  1.200,0.6393,
    1.250,0.6518,  1.300,0.6589,  1.350,0.6621,
    1.400,0.6625,  1.450,0.6607,  1.500,0.6573,
    1.550,0.6528,  1.600,0.6474,  1.650,0.6413,
    1.700,0.6347,  1.750,0.6280,  1.800,0.6210,
    1.850,0.6141,  1.900,0.6072,  1.950,0.6003,
    2.000,0.5934,  2.050,0.5867,  2.100,0.5804,
    2.150,0.5743,  2.200,0.5685,  2.250,0.5630,
    2.300,0.5577,  2.350,0.5527,  2.400,0.5481,
    2.450,0.5438,  2.500,0.5397,  2.600,0.5325,
    2.700,0.5264,  2.800,0.5211,  2.900,0.5168,
    3.000,0.5133,  3.100,0.5105,  3.200,0.5084,
    3.300,0.5067,  3.400,0.5054,  3.500,0.5040,
    3.600,0.5030,  3.700,0.5022,  3.800,0.5016,
    3.900,0.5010,  4.000,0.5006,  4.200,0.4998,
    4.400,0.4995,  4.600,0.4992,  4.800,0.4990,
    5.000,0.4988,
  ];

  function parseTable(raw) {
    const mach = [];
    const cd   = [];
    for (let i = 0; i < raw.length; i += 2) {
      mach.push(raw[i]);
      cd.push(raw[i + 1]);
    }
    return { mach, cd };
  }

  const G7_TABLE = parseTable(G7_RAW);
  const G1_TABLE = parseTable(G1_RAW);

  // -- Linear interpolation --
  function searchSortedLast(arr, val) {
    let lo = 0, hi = arr.length - 1;
    if (val < arr[0]) return -1;
    if (val >= arr[hi]) return hi;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (arr[mid] <= val) lo = mid; else hi = mid - 1;
    }
    return lo;
  }

  function interpTable(table, mach) {
    const M = table.mach;
    const C = table.cd;
    const n = M.length;
    if (mach <= M[0])     return C[0];
    if (mach >= M[n - 1]) return C[n - 1];
    let idx = searchSortedLast(M, mach);
    if (idx < 0) idx = 0;
    if (idx >= n - 1) idx = n - 2;
    const frac = (mach - M[idx]) / (M[idx + 1] - M[idx]);
    return C[idx] + frac * (C[idx + 1] - C[idx]);
  }

  // -- Cubic spline (C² continuous) --
  function buildCubicSpline(xArr, yArr) {
    const n = xArr.length;
    const h = new Float64Array(n - 1);
    const a = Float64Array.from(yArr);
    for (let i = 0; i < n - 1; i++) h[i] = xArr[i + 1] - xArr[i];

    const alpha = new Float64Array(n);
    for (let i = 1; i < n - 1; i++) {
      alpha[i] = 3.0 * ((a[i + 1] - a[i]) / h[i] - (a[i] - a[i - 1]) / h[i - 1]);
    }

    const l  = new Float64Array(n); l[0] = 1.0;
    const mu = new Float64Array(n);
    const z  = new Float64Array(n);

    for (let i = 1; i < n - 1; i++) {
      l[i]  = 2.0 * (xArr[i + 1] - xArr[i - 1]) - h[i - 1] * mu[i - 1];
      mu[i] = h[i] / l[i];
      z[i]  = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    const c = new Float64Array(n);
    const b = new Float64Array(n - 1);
    const d = new Float64Array(n - 1);

    for (let j = n - 2; j >= 0; j--) {
      c[j] = z[j] - mu[j] * c[j + 1];
      b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
      d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
    }

    return {
      x: Float64Array.from(xArr),
      a: a.slice(0, n - 1),
      b, c: c.slice(0, n - 1), d,
      n,
    };
  }

  function evalSpline(sp, xq) {
    if (xq <= sp.x[0]) return sp.a[0];
    if (xq >= sp.x[sp.n - 1]) {
      const k  = sp.a.length - 1;
      const dx = sp.x[sp.n - 1] - sp.x[k];
      return sp.a[k] + sp.b[k] * dx + sp.c[k] * dx * dx + sp.d[k] * dx * dx * dx;
    }
    let k = searchSortedLast(sp.x, xq);
    if (k < 0) k = 0;
    if (k >= sp.a.length) k = sp.a.length - 1;
    const dx = xq - sp.x[k];
    return sp.a[k] + sp.b[k] * dx + sp.c[k] * dx * dx + sp.d[k] * dx * dx * dx;
  }

  const G7_SPLINE = buildCubicSpline(G7_TABLE.mach, G7_TABLE.cd);
  const G1_SPLINE = buildCubicSpline(G1_TABLE.mach, G1_TABLE.cd);

  const cdG1       = (mach) => interpTable(G1_TABLE, mach);
  const cdG7       = (mach) => interpTable(G7_TABLE, mach);
  const cdG1Spline = (mach) => evalSpline(G1_SPLINE, mach);
  const cdG7Spline = (mach) => evalSpline(G7_SPLINE, mach);

  const dragCoefficient = (model, mach, spline = true) => {
    if (model === "G1") return spline ? cdG1Spline(mach) : cdG1(mach);
    if (model === "G7") return spline ? cdG7Spline(mach) : cdG7(mach);
    throw new Error(`Unknown drag model: ${model}. Supported: "G1", "G7"`);
  };

  // -- Custom drag model from Doppler radar or CFD --
  function buildCustomDrag(machVec, cdVec, name = "Custom") {
    const indices = machVec.map((_, i) => i).sort((a, b) => machVec[a] - machVec[b]);
    const mach = indices.map(i => machVec[i]);
    const cd   = indices.map(i => cdVec[i]);
    const spline = buildCubicSpline(mach, cd);
    return { mach, cd, spline, name };
  }

  const cdmEval       = (cdm, mach) => evalSpline(cdm.spline, mach);
  const cdmEvalLinear = (cdm, mach) => interpTable({ mach: cdm.mach, cd: cdm.cd }, mach);

  const cdFromRadar = (velocity, dvDt, massKg, rho, A) =>
    -2.0 * massKg * dvDt / (rho * A * velocity * velocity);

  // -- Single-velocity BC conversion between reference models --
  // A bullet with BC_from has form factor i_from = SD/BC_from relative to the
  // `from` reference, so its absolute drag is Cd = i_from·Cd_from(M). Expressed
  // against the `to` reference, i_to = Cd/Cd_to(M) and BC_to = SD/i_to, hence:
  //     BC_to = BC_from · Cd_to(M) / Cd_from(M).
  // The ratio varies with Mach, so a representative velocity must be supplied —
  // this is exactly why catalogue G1 BCs are velocity-banded. Returns the
  // converted BC and the drag-curve ratio used.
  function convertBC(bc, fromModel, toModel, mach) {
    if (fromModel === toModel) return { bc, ratio: 1.0, mach };
    const cdFrom = dragCoefficient(fromModel, mach);
    const cdTo   = dragCoefficient(toModel, mach);
    const ratio  = cdTo / cdFrom;
    return { bc: bc * ratio, ratio, mach };
  }

  return {
    G7_TABLE, G1_TABLE, G7_SPLINE, G1_SPLINE,
    cdG1, cdG7, cdG1Spline, cdG7Spline,
    dragCoefficient,
    buildCubicSpline, evalSpline,
    buildCustomDrag, cdmEval, cdmEvalLinear,
    cdFromRadar, convertBC,
  };
})();


// --------------------------------------------------------------------------
// MODULE 5: ExteriorBallistics
// --------------------------------------------------------------------------
const ExteriorBallistics = (() => {

  const {
    grainsToKg, inchesToM, fpsToMs, msToFps, mToInches, mToYards,
    yardsToM, fahrenheitToKelvin, inhgToPa, moaToRad,
    millerStability, dropToMoa,
  } = BallisticUtils;

  const { airDensity, speedOfSound, g0 } = Atmosphere;
  const { dragCoefficient } = DragModels;

  const DEFAULT_SHOT = {
    massGrains:      175.0,
    caliberIn:       0.308,
    bc:              0.275,
    dragModel:       "G7",
    muzzleVelFps:    2600.0,
    sightHeightIn:   1.5,
    zeroRangeYd:     100.0,
    elevationMoa:    0.0,
    windageMoa:      0.0,
    tempF:           59.0,
    pressureInhg:    29.92,
    humidityPct:     0.0,
    altitudeFt:      0.0,
    pressureIsSeaLevel: false,   // true = la pression saisie est un QNH, à ramener via altitudeFt
    windSpeedMph:    0.0,
    windAngleDeg:    90.0,
    targetRangeYd:   1000.0,
    inclineDeg:      0.0,
    latitudeDeg:     45.0,
    azimuthDeg:      0.0,
    enableCoriolis:  true,
    twistIn:         10.0,
    twistDirection:  1,
    bulletLengthIn:  1.24,
    enableSpinDrift: true,
    dt:              0.0005,
  };

  function makeShotParams(overrides = {}) {
    return { ...DEFAULT_SHOT, ...overrides };
  }

  // Pression réellement vue par le projectile.
  //
  // `altitudeFt` a longtemps été déclaré dans DEFAULT_SHOT et transmis par le
  // formulaire sans être lu NULLE PART : seule la pression saisie agissait, si
  // bien qu'indiquer 1500 m d'altitude en laissant 1013 hPa rendait des
  // résultats de niveau de la mer, sans le moindre signalement.
  //
  // Deux conventions coexistent chez les tireurs : la pression *station*, lue
  // sur place, et la pression ramenée au niveau de la mer (QNH), celle des
  // bulletins météo. La première se suffit à elle-même — l'altitude n'apporte
  // alors rien de plus. La seconde doit être ramenée au point de tir, et c'est
  // là que l'altitude sert.
  //
  // Ce helper est le SEUL point de vérité : _toSi et millerSg passent tous deux
  // par lui, précisément pour qu'ils ne puissent plus diverger.
  function stationPressureInhg(p) {
    if (!p.pressureIsSeaLevel) return p.pressureInhg;
    const hM = (p.altitudeFt || 0) * 0.3048;
    return p.pressureInhg * (Atmosphere.stdPressure(hM) / Atmosphere.P0);
  }

  function _toSi(p) {
    return {
      m:   grainsToKg(p.massGrains),
      d:   inchesToM(p.caliberIn),
      v0:  fpsToMs(p.muzzleVelFps),
      T:   fahrenheitToKelvin(p.tempF),
      P:   inhgToPa(stationPressureInhg(p)),
      H:   p.humidityPct,
      w:   p.windSpeedMph * 0.44704,
      sh:  inchesToM(p.sightHeightIn),
      zr:  yardsToM(p.zeroRangeYd),
      tr:  yardsToM(p.targetRangeYd),
      L:   inchesToM(p.bulletLengthIn),
      tw:  inchesToM(p.twistIn),
    };
  }

  function millerSg(p) {
    return millerStability({
      massGr:        p.massGrains,
      caliberIn:     p.caliberIn,
      bulletLengthIn: p.bulletLengthIn,
      twistIn:       p.twistIn,
      muzzleVelFps:  p.muzzleVelFps,
      tempF:         p.tempF,
      pressureInhg:  stationPressureInhg(p),
    });
  }

  function spinDriftInches(t, sg, direction = 1) {
    return direction * 1.25 * (sg + 1.2) * Math.pow(t, 1.83);
  }

  function coriolisHorizontal(rangeM, tof, latitudeDeg) {
    const omega = 7.2921e-5;
    return omega * Math.sin(latitudeDeg * Math.PI / 180.0) * rangeM * tof;
  }

  function eotvosVertical(rangeM, tof, latitudeDeg, azimuthDeg) {
    const omega = 7.2921e-5;
    return omega * Math.cos(latitudeDeg * Math.PI / 180.0) *
           Math.sin(azimuthDeg * Math.PI / 180.0) * rangeM * tof;
  }

  function windDeflectionLagRule(windCrossMs, tof, rangeM, v0Ms) {
    const tVac = rangeM / v0Ms;
    return windCrossMs * (tof - tVac);
  }

  function aerodynamicJump(cLAlpha, kT, alphaTrim, pS, rS) {
    return (cLAlpha / (2.0 * kT * kT)) * (alphaTrim / (pS - rS));
  }

  function inclinedFireEffectiveRange(slantRange, angleDeg) {
    return slantRange * Math.cos(angleDeg * Math.PI / 180.0);
  }

  // -- Zero angle finder --
  function findZeroAngle(p, tol = 1e-6, maxIter = 50) {
    const si  = _toSi(p);
    const rho = airDensity({ P: si.P, T: si.T, H: si.H });
    const aS  = speedOfSound(si.T);
    const A   = Math.PI * si.d * si.d / 4.0;
    // Form factor i = SD/BC scales the reference drag to the actual bullet.
    const ff  = BallisticUtils.formFactor(p.massGrains, p.caliberIn, p.bc);
    let theta = Math.atan(g0 * si.zr / (2.0 * si.v0 * si.v0));
    const dt  = p.dt;

    for (let iter = 0; iter < maxIter; iter++) {
      let x = 0.0, y = -si.sh, z = 0.0;
      let vx = si.v0 * Math.cos(theta);
      let vy = si.v0 * Math.sin(theta);
      let vz = 0.0;
      let t  = 0.0;

      function derivsZero(sv) {
        const [_x, _y, _z, _vx, _vy, _vz] = sv;
        const v = Math.sqrt(_vx * _vx + _vy * _vy);
        const Ma = v / aS;
        const cd = dragCoefficient(p.dragModel, Ma);
        const df = (rho * cd * A * ff) / (2.0 * si.m);   // cf. `derivs` : plus de 4/π
        return [_vx, _vy, 0, -df * v * _vx, -df * v * _vy - g0, 0];
      }

      while (x < si.zr && t < 15.0) {
        const s1 = [x, y, z, vx, vy, vz];
        const k1 = derivsZero(s1);
        const s2 = s1.map((v, i) => v + 0.5 * dt * k1[i]);
        const k2 = derivsZero(s2);
        const s3 = s1.map((v, i) => v + 0.5 * dt * k2[i]);
        const k3 = derivsZero(s3);
        const s4 = s1.map((v, i) => v + dt * k3[i]);
        const k4 = derivsZero(s4);
        const nextS = s1.map((v, i) => v + (dt / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
        [x, y, z, vx, vy, vz] = nextS;
        t += dt;
      }

      if (Math.abs(y) < tol) return theta;
      theta -= y / si.zr * Math.cos(theta);
    }
    return theta;
  }

  // -- Full 3-DOF trajectory solver (RK4) --
  function solveTrajectory(params) {
    const p   = makeShotParams(params);
    const si  = _toSi(p);
    const rho = airDensity({ P: si.P, T: si.T, H: p.humidityPct });
    const aS  = speedOfSound(si.T);
    const A   = Math.PI * si.d * si.d / 4.0;
    // Form factor i = SD/BC scales the reference drag to the actual bullet.
    const ff  = BallisticUtils.formFactor(p.massGrains, p.caliberIn, p.bc);
    const omegaE = 7.2921e-5;
    const lat = p.latitudeDeg * Math.PI / 180.0;
    const azm = p.azimuthDeg  * Math.PI / 180.0;   // firing azimuth, clockwise from North
    const dt  = p.dt;

    // Earth-rotation vector in the shot frame (x = downrange/LOS, y = up, z = right).
    // Ω in (East, North, Up) = ω(0, cos lat, sin lat); projected onto a shot frame
    // whose downrange axis points at azimuth `azm` from North:
    //   downrange = (sin azm, cos azm, 0), up = (0,0,1), right = (cos azm, -sin azm, 0).
    const oeX = omegaE * Math.cos(lat) * Math.cos(azm);
    const oeY = omegaE * Math.sin(lat);
    const oeZ = -omegaE * Math.cos(lat) * Math.sin(azm);

    // Wind components
    const wAng = p.windAngleDeg * Math.PI / 180.0;
    const wx = -si.w * Math.cos(wAng);
    const wz =  si.w * Math.sin(wAng);
    const wy = 0.0;

    // Incline: resolve gravity in the LOS frame. The component perpendicular to the
    // line of sight (g·cos α) drives the drop — the physical basis of the
    // "rifleman's rule" — while g·sin α acts along the flight path. α = 0 ⇒ flat fire.
    const inc    = p.inclineDeg * Math.PI / 180.0;
    const gAlong = -g0 * Math.sin(inc);
    const gPerp  = -g0 * Math.cos(inc);

    // Zero angle + dialed corrections
    let theta0 = findZeroAngle(p);
    theta0 += moaToRad(p.elevationMoa);
    const psi0 = moaToRad(p.windageMoa);

    const sg = millerSg(p);

    // Initial state
    let x  = 0.0, y  = -si.sh, z  = 0.0;
    let vx = si.v0 * Math.cos(theta0) * Math.cos(psi0);
    let vy = si.v0 * Math.sin(theta0);
    let vz = si.v0 * Math.cos(theta0) * Math.sin(psi0);
    let t  = 0.0;

    const results = [];

    function derivs(sv) {
      const [_x, _y, _z, _vx, _vy, _vz] = sv;
      const vrx = _vx - wx, vry = _vy - wy, vrz = _vz - wz;
      const vrel = Math.sqrt(vrx * vrx + vry * vry + vrz * vrz);
      const Ma = vrel / aS;
      const cd = dragCoefficient(p.dragModel, Ma);
      // Équation de traînée point-masse : a = ρ·C_D·i·A·v²/(2m), avec A = πd²/4 l'aire
      // frontale et i = SD/BC le facteur de forme. Le d² et la masse se simplifient —
      // c'est le coefficient balistique seul qui pilote la trajectoire, comme il se doit.
      //
      // ⚠️ IL Y AVAIT ICI UN FACTEUR 4/π, retiré le 2026-08-22. Son commentaire le
      // défendait comme l'appariement de la densité sectionnelle (définie sur d²) à la
      // table BRL (référée à πd²/4), et interdisait de le « simplifier ». La dérivation
      // ne le donne pas. Ce qu'il faisait réellement, c'était compenser des tables de
      // traînée supersoniques fausses, de sorte qu'une .308 175 gr arrive à 1000 yd sur
      // une vitesse plausible. Un ajustement global calé sur un point de fonctionnement
      // rend juste ce point-là et faux tous les autres.
      const D  = (rho * cd * A * ff) / (2.0 * si.m);

      let ax = -D * vrel * vrx + gAlong;
      let ay = -D * vrel * vry + gPerp;
      let az = -D * vrel * vrz;

      if (p.enableCoriolis) {
        ax += -2.0 * (oeY * _vz - oeZ * _vy);
        ay += -2.0 * (oeZ * _vx - oeX * _vz);
        az += -2.0 * (oeX * _vy - oeY * _vx);
      }
      return [_vx, _vy, _vz, ax, ay, az];
    }

    // ⚠️ LA BOUCLE ENREGISTRE AVANT D'INTÉGRER, et doit donc dépasser la cible d'UN pas.
    // Écrite `while (x <= si.tr)`, elle s'arrêtait au dernier point situé AVANT la
    // portée demandée : `trajectoryTable` n'atteignait jamais sa dernière ligne, et une
    // demande à 1000 m rendait une table finissant à 900. Corrigé le 2026-08-22 ; c'est
    // la référence tierce qui l'a signalé, en réclamant un point que nous n'avions pas.
    while (t < 15.0) {
      const vTot = Math.sqrt(vx * vx + vy * vy + vz * vz);
      const mach = vTot / aS;
      const ek   = 0.5 * si.m * vTot * vTot;

      results.push({
        time:     t,
        rangeM:   x,
        dropM:    y,
        windageM: z,          // TOTAL : vent + Coriolis + dérive gyroscopique
        spinDriftM: 0.0,      // rempli après l'intégration
        vx, vy, vz,
        vTotal:   vTot,
        mach,
        energyJ:  ek,
      });

      if (x > si.tr) break;   // un point au-delà de la cible est enregistré, pas deux

      // RK4 integration for state vector [x, y, z, vx, vy, vz]
      const s1 = [x, y, z, vx, vy, vz];
      const k1 = derivs(s1);

      const s2 = s1.map((v, i) => v + 0.5 * dt * k1[i]);
      const k2 = derivs(s2);

      const s3 = s1.map((v, i) => v + 0.5 * dt * k2[i]);
      const k3 = derivs(s3);

      const s4 = s1.map((v, i) => v + dt * k3[i]);
      const k4 = derivs(s4);

      const nextS = s1.map((v, i) => v + (dt / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));

      [x, y, z, vx, vy, vz] = nextS;
      t += dt;
    }

    // Post-hoc spin drift
    if (p.enableSpinDrift) {
      for (let i = 0; i < results.length; i++) {
        const pt = results[i];
        const sdM = spinDriftInches(pt.time, sg, p.twistDirection) * 0.0254;
        pt.windageM += sdM;
        pt.spinDriftM = sdM;
      }
    }

    return results;
  }

  // -- Trajectory table at regular range intervals --
  // Chaque ligne est INTERPOLÉE à la distance demandée, elle n'est pas le point
  // d'intégration le plus proche.
  //
  // La table retenait auparavant le premier point vérifiant `rangeM >= nextR`, donc un
  // point qui dépasse toujours un peu — d'un pas d'intégration, soit ~0,8 m à 800 m/s.
  // Une demande à 25 m d'intervalle rendait « 275, 300, 326, 350, 375, 401, 425 » :
  // l'écart est minuscule, mais une DOPE se recopie à la main et ces 326 et 401 n'ont
  // aucun sens pour le tireur qui a demandé un pas de 25. Le dépassement se cumulait
  // en outre avec un second arrondi côté affichage, l'aller-retour mètres → yards.
  //
  // Entre deux points d'intégration séparés de moins d'un mètre, l'interpolation
  // linéaire est sans conséquence physique : la RK4 elle-même n'est pas plus fine.
  function trajectoryTable(params, stepYd = 100.0) {
    const p    = makeShotParams(params);
    const traj = solveTrajectory(p);
    const stepM = yardsToM(stepYd);
    const rows  = [];
    let nextR = stepM;

    for (let i = 1; i < traj.length; i++) {
      const a = traj[i - 1], b = traj[i];
      // `while` et non `if` : près de la bouche, un pas d'intégration peut franchir
      // deux distances demandées d'un coup si le pas de table est très serré.
      while (b.rangeM >= nextR && a.rangeM < nextR) {
        const span = b.rangeM - a.rangeM;
        const f    = span > 0 ? (nextR - a.rangeM) / span : 0;
        const lerp = (u, v) => u + (v - u) * f;

        const rYd    = mToYards(nextR);
        const dropM  = lerp(a.dropM, b.dropM);
        const dropIn = mToInches(dropM);
        const windIn = mToInches(lerp(a.windageM, b.windageM));
        const spinIn = mToInches(lerp(a.spinDriftM || 0, b.spinDriftM || 0));
        let   elev   = dropToMoa(Math.abs(dropIn), rYd);
        elev = dropM < 0 ? elev : -elev;
        const vFps   = msToFps(lerp(a.vTotal, b.vTotal));
        const ekFtlb = lerp(a.energyJ, b.energyJ) / 1.35582;

        rows.push({
          rangeM:     +nextR.toFixed(2),   // exact : la distance DEMANDÉE
          rangeYd:    Math.round(rYd),
          dropIn:     +dropIn.toFixed(1),
          elevMoa:    +elev.toFixed(1),
          windIn:     +windIn.toFixed(1),
          spinDriftIn: +spinIn.toFixed(1),
          velFps:     Math.round(vFps),
          tofS:       +lerp(a.time, b.time).toFixed(3),
          energyFtlb: Math.round(ekFtlb),
        });
        nextR += stepM;
      }
    }
    return rows;
  }

  return {
    DEFAULT_SHOT, makeShotParams,
    millerSg, spinDriftInches, stationPressureInhg,
    coriolisHorizontal, eotvosVertical,
    windDeflectionLagRule,
    aerodynamicJump, inclinedFireEffectiveRange,
    findZeroAngle, solveTrajectory, trajectoryTable,
  };
})();


// --------------------------------------------------------------------------
// Export
// --------------------------------------------------------------------------
// UMD-style: works as ES module, CommonJS, or browser global.
const CompetitionBallistics = {
  BallisticUtils,
  ReferenceData,
  ParameterMetadata,
  Atmosphere,
  InteriorBallistics,
  DragModels,
  ExteriorBallistics,
};

/* eslint-disable no-undef */
if (typeof exports === "object" && typeof module === "object") {
  module.exports = CompetitionBallistics;
} else if (typeof define === "function" && define.amd) {
  define([], () => CompetitionBallistics);
} else if (typeof globalThis !== "undefined") {
  globalThis.CompetitionBallistics = CompetitionBallistics;
}
/* eslint-enable no-undef */
