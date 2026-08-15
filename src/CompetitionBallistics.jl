# ============================================================================
# CompetitionBallistics.jl
# A Complete Ballistic Computation Package for Competition Rifle Shooting
#
# Companion code to: "Competition Rifle Ballistics: A Comprehensive
#                     Mathematical Manual"
#
# Modules:
#   - BallisticUtils      : Unit conversions, sectional density, stability
#   - Atmosphere          : ICAO standard atmosphere, humidity, density altitude
#   - InteriorBallistics  : Le Duc model, pressure, burn rate, temp sensitivity
#   - DragModels          : G1/G7 BRL tabular data with interpolation
#   - ExteriorBallistics  : Full 3-DOF solver (RK4), Coriolis, spin drift, wind
#   - ReloadingAnalysis   : Ladder test analysis, ES/SD, vertical dispersion,
#                           charge weight optimization, seating depth sensitivity
#
# Usage:
#   include("CompetitionBallistics.jl")
#   using .CompetitionBallistics
# ============================================================================

module CompetitionBallistics

# ──────────────────────────────────────────────────────────────────────────────
# MODULE 1: BallisticUtils
# ──────────────────────────────────────────────────────────────────────────────
module BallisticUtils

export grains_to_kg, kg_to_grains, grains_to_grams,
       fps_to_ms, ms_to_fps,
       inches_to_m, m_to_inches, inches_to_mm, yards_to_m, m_to_yards,
       fahrenheit_to_kelvin, kelvin_to_fahrenheit, fahrenheit_to_rankine,
       celsius_to_kelvin, kelvin_to_celsius,
       inhg_to_pa, pa_to_inhg, inhg_to_hpa,
       moa_to_rad, rad_to_moa, mil_to_rad, rad_to_mil,
       moa_to_mil, mil_to_moa,
       drop_to_moa, drop_to_mil,
       kinetic_energy_J, kinetic_energy_ftlbs,
       sectional_density, form_factor,
       miller_stability, greenhill_twist

# ── Mass conversions ──
grains_to_kg(gr::Real)    = gr * 6.47989e-5
kg_to_grains(kg::Real)    = kg / 6.47989e-5
grains_to_grams(gr::Real) = gr * 0.06480

# ── Velocity conversions ──
fps_to_ms(fps::Real)  = fps * 0.3048
ms_to_fps(ms::Real)   = ms / 0.3048

# ── Length conversions ──
inches_to_m(x::Real)   = x * 0.0254
m_to_inches(x::Real)   = x / 0.0254
inches_to_mm(x::Real)  = x * 25.4
yards_to_m(x::Real)    = x * 0.9144
m_to_yards(x::Real)    = x / 0.9144

# ── Temperature conversions ──
fahrenheit_to_kelvin(f::Real)  = (f - 32.0) * 5.0/9.0 + 273.15
kelvin_to_fahrenheit(k::Real)  = (k - 273.15) * 9.0/5.0 + 32.0
fahrenheit_to_rankine(f::Real) = f + 459.67
celsius_to_kelvin(c::Real)     = c + 273.15
kelvin_to_celsius(k::Real)     = k - 273.15

# ── Pressure conversions ──
inhg_to_pa(x::Real)   = x * 3386.389
pa_to_inhg(x::Real)   = x / 3386.389
inhg_to_hpa(x::Real)  = x * 33.8639

# ── Angle conversions ──
moa_to_rad(moa::Real)  = moa * (π / (180.0 * 60.0))
rad_to_moa(rad::Real)  = rad * (180.0 * 60.0 / π)
mil_to_rad(mil::Real)  = mil * 0.001
rad_to_mil(rad::Real)  = rad * 1000.0
moa_to_mil(moa::Real)  = moa * 0.29089
mil_to_moa(mil::Real)  = mil * 3.43775

"""
    drop_to_moa(drop_inches, range_yards)

Convert a bullet drop [inches] at a given range [yards] to MOA.
"""
function drop_to_moa(drop_inches::Real, range_yards::Real)
    return drop_inches / (range_yards * 1.04720 / 100.0)
end

"""
    drop_to_mil(drop_inches, range_yards)

Convert a bullet drop [inches] at a given range [yards] to milliradians.
"""
function drop_to_mil(drop_inches::Real, range_yards::Real)
    range_m = yards_to_m(range_yards)
    drop_m  = inches_to_m(drop_inches)
    return drop_m / range_m * 1000.0
end

"""
    kinetic_energy_J(mass_kg, vel_ms)

Kinetic energy [Joules] from mass [kg] and velocity [m/s].
"""
kinetic_energy_J(mass_kg::Real, vel_ms::Real) = 0.5 * mass_kg * vel_ms^2

"""
    kinetic_energy_ftlbs(mass_gr, vel_fps)

Kinetic energy [ft-lbs] from mass [grains] and velocity [fps].
"""
kinetic_energy_ftlbs(mass_gr::Real, vel_fps::Real) = mass_gr * vel_fps^2 / 450436.0

"""
    sectional_density(mass_gr, caliber_in)

Sectional density [lb/in²]: (mass in lb) / caliber² .
"""
function sectional_density(mass_gr::Real, caliber_in::Real)
    return (mass_gr / 7000.0) / caliber_in^2
end

"""
    form_factor(mass_gr, caliber_in, bc)

Form factor i = SD / BC.
"""
function form_factor(mass_gr::Real, caliber_in::Real, bc::Real)
    return sectional_density(mass_gr, caliber_in) / bc
end

"""
    miller_stability(; mass_gr, caliber_in, bullet_length_in, twist_in,
                       muzzle_vel_fps=2800, temp_F=59, pressure_inhg=29.92)

Gyroscopic stability factor Sg using the Miller twist rule.
Returns Sg; values > 1.3 indicate adequate stability.

The base formula is quoted at 2800 fps in standard atmosphere; the velocity and
atmospheric corrections bring it to the actual conditions. Both are Miller's own
and neither is optional at the edges of the envelope — see the comments below.
"""
function miller_stability(;
    mass_gr::Real,
    caliber_in::Real,
    bullet_length_in::Real,
    twist_in::Real,
    muzzle_vel_fps::Real = 2800.0,
    temp_F::Real = 59.0,
    pressure_inhg::Real = 29.92
)
    d  = caliber_in
    l  = bullet_length_in / d       # length in calibers
    tw = twist_in / d               # twist in calibers per turn
    T_R = fahrenheit_to_rankine(temp_F)
    sg = 30.0 * mass_gr / (tw^2 * d^3 * l * (1.0 + l^2))
    # Miller's velocity correction is the CUBE ROOT of the ratio to 2800 fps,
    # not the ratio itself. The two agree near 2800 fps and part company at the
    # edges: 27 % apart at 4000 fps.
    sg *= cbrt(muzzle_vel_fps / 2800.0)
    # Atmospheric correction: Sg varies as 1/ρ, and ρ ~ P/T. Warm or thin air
    # (altitude) therefore RAISES Sg — the ratios go this way up, not the other.
    sg *= (T_R / 518.67) * (29.92 / pressure_inhg)
    return sg
end

"""
    greenhill_twist(caliber_in, bullet_length_in; C=150)

Greenhill's formula: recommended twist rate [inches/turn].
Use C=150 for v < 2800 fps, C=180 for v > 2800 fps.
"""
function greenhill_twist(caliber_in::Real, bullet_length_in::Real; C::Real=150.0)
    return C * caliber_in^2 / bullet_length_in
end

end # module BallisticUtils


# ──────────────────────────────────────────────────────────────────────────────
# MODULE 1.1: ReferenceData
# ──────────────────────────────────────────────────────────────────────────────
module ReferenceData

export BulletProfile, COMMON_BULLETS, CARTRIDGE_MAX_PSI

"""Data profile for a specific bullet."""
struct BulletProfile
    name::String
    mass_gr::Float64
    caliber_in::Float64
    bc_g7::Float64
    length_in::Float64
    twist_in::Float64 # Recommended twist
end

# Data from Chapter 11 of the manual and Tireur.org integration
const COMMON_BULLETS = Dict{String, BulletProfile}(
    "Berger 105 Hybrid (6mm)" => BulletProfile("Berger 105 Hybrid (6mm)", 105.0, 0.243, 0.275, 1.10, 8.0),
    "Berger 140 Hybrid (6.5mm)" => BulletProfile("Berger 140 Hybrid (6.5mm)", 140.0, 0.264, 0.311, 1.34, 8.0),
    "Hornady 147 ELD-M (6.5mm)" => BulletProfile("Hornady 147 ELD-M (6.5mm)", 147.0, 0.264, 0.351, 1.42, 7.5),
    "Berger 180 Hybrid (7mm)" => BulletProfile("Berger 180 Hybrid (7mm)", 180.0, 0.284, 0.350, 1.50, 8.5),
    "Sierra 175 MK (.308)" => BulletProfile("Sierra 175 MK (.308)", 175.0, 0.308, 0.259, 1.24, 10.0),
    "Berger 185 Juggernaut (.308)" => BulletProfile("Berger 185 Juggernaut (.308)", 185.0, 0.308, 0.283, 1.30, 10.0),
    "Hornady 178 ELD-M (.308)" => BulletProfile("Hornady 178 ELD-M (.308)", 178.0, 0.308, 0.274, 1.31, 10.0),
    "Berger 300 Hybrid (.338)" => BulletProfile("Berger 300 Hybrid (.338)", 300.0, 0.338, 0.419, 1.82, 9.4)
)

# Data from Chapter 9 of the manual
const CARTRIDGE_MAX_PSI = Dict{String, Float64}(
    ".223 Remington" => 55000.0,
    "6.5 Creedmoor" => 62000.0,
    ".308 Winchester" => 62000.0,
    ".300 Winchester Magnum" => 64000.0,
    ".338 Lapua Magnum" => 60916.0
)

end # module ReferenceData


# ──────────────────────────────────────────────────────────────────────────────
# MODULE 2: Atmosphere
# ──────────────────────────────────────────────────────────────────────────────
module Atmosphere

export T0, P0, rho0, L, g0, R_air, gamma_air,
       std_temperature, std_pressure,
       saturation_vapor_pressure, air_density,
       speed_of_sound, density_altitude, density_ratio

# ICAO Standard Atmosphere sea-level constants
const T0        = 288.15      # K
const P0        = 101325.0    # Pa
const rho0      = 1.2250      # kg/m³
const L         = 0.0065      # K/m  (tropospheric lapse rate)
const g0        = 9.80665     # m/s²
const R_air     = 287.0528    # J/(kg·K) specific gas constant for dry air (ISA 1976)
const gamma_air = 1.4         # ratio of specific heats (Cp/Cv)

"""
    std_temperature(h) -> T [K]

Temperature at altitude `h` [m] in the ICAO standard atmosphere (troposphere).
"""
std_temperature(h::Real) = T0 - L * h

"""
    std_pressure(h) -> P [Pa]

Pressure at altitude `h` [m] in the ICAO standard atmosphere (troposphere).
Barometric formula: P = P0 * (T/T0)^(g0 / (L·R_air))
"""
function std_pressure(h::Real)
    T = std_temperature(h)
    return P0 * (T / T0)^(g0 / (L * R_air))
end

"""
    saturation_vapor_pressure(T_K) -> e_s [Pa]

August–Roche–Magnus formula for saturation vapor pressure.
"""
function saturation_vapor_pressure(T_K::Real)
    T_C = T_K - 273.15
    return 610.78 * 10.0^(7.5 * T_C / (T_C + 237.3))
end

"""
    air_density(; P=101325, T=288.15, H=0) -> ρ [kg/m³]

Air density with humidity correction (virtual temperature method).
- `P`: ambient pressure [Pa]
- `T`: temperature [K]
- `H`: relative humidity [0–100 %]
"""
function air_density(; P::Real=P0, T::Real=T0, H::Real=0.0)
    e_s = saturation_vapor_pressure(T)
    e   = (H / 100.0) * e_s
    return (P - 0.37802 * e) / (R_air * T)
end

"""
    speed_of_sound(T=288.15) -> a [m/s]

Speed of sound in air at temperature T [K].
a = √(γ · R_air · T)
"""
speed_of_sound(T::Real=T0) = sqrt(gamma_air * R_air * T)

"""
    density_ratio(; P, T, H) -> ρ/ρ₀

Ratio of actual air density to ICAO standard density.
"""
function density_ratio(; P::Real=P0, T::Real=T0, H::Real=0.0)
    return air_density(P=P, T=T, H=H) / rho0
end

"""
    density_altitude(ρ) -> DA [m]

Density altitude: altitude in the standard atmosphere with density ρ [kg/m³].
Formula: h = (T0/L) * (1 - (rho/rho0)^(1 / (g/(L*R) - 1)))
Exponent ≈ 1 / (5.2559 - 1) ≈ 0.234969
"""
function density_altitude(rho::Real)
    return (T0 / L) * (1.0 - (rho / rho0)^0.234969)
end

"""
    density_altitude_ft(ρ) -> DA [ft]

Density altitude in feet (practical approximation).
"""
function density_altitude_ft(rho::Real)
    return 145442.16 * (1.0 - (rho / rho0)^0.234969)
end

end # module Atmosphere


# ──────────────────────────────────────────────────────────────────────────────
# MODULE 3: InteriorBallistics
# ──────────────────────────────────────────────────────────────────────────────
module InteriorBallistics

using ..BallisticUtils

export leduc_velocity, leduc_pressure, peak_pressure,
       barrel_length_correction, muzzle_velocity_temp_correction,
       burn_fraction, propellant_energy, effective_mass,
       chamber_pressure_closed_bomb

"""
    leduc_velocity(x; a, b) -> v [m/s]

Le Duc velocity model: v(x) = a·x / (b + x).
- `x`: bullet travel distance in barrel [m]
- `a`: theoretical maximum velocity [m/s]
- `b`: shape parameter [m]
"""
function leduc_velocity(x::Real; a::Real, b::Real)
    return a * x / (b + x)
end

"""
    leduc_pressure(x; a, b, m_bullet, m_charge, d_bore) -> P [Pa]

Pressure from Le Duc velocity model.
- `x`: bullet travel [m]
- `a`, `b`: Le Duc parameters
- `m_bullet`, `m_charge`: bullet and charge mass [kg]
- `d_bore`: bore diameter [m]
"""
function leduc_pressure(x::Real; a::Real, b::Real,
                        m_bullet::Real, m_charge::Real,
                        d_bore::Real)
    m_e = m_bullet + m_charge / 3.0
    A = π * d_bore^2 / 4.0
    return m_e * a^2 * b * x / (A * (b + x)^3)
end

"""
    peak_pressure(; a, b, m_bullet, m_charge, d_bore) -> P_max [Pa]

Peak pressure from Le Duc model, occurs at x = b/2.
P_max = 4·m_e·a² / (27·A·b)
"""
function peak_pressure(; a::Real, b::Real,
                       m_bullet::Real, m_charge::Real,
                       d_bore::Real)
    m_e = m_bullet + m_charge / 3.0
    A = π * d_bore^2 / 4.0
    return 4.0 * m_e * a^2 / (27.0 * A * b)
end

"""
    effective_mass(m_bullet, m_charge) -> m_e [kg]

Effective mass = bullet mass + 1/3 charge mass.
Accounts for the kinetic energy carried by accelerating propellant gas.
"""
effective_mass(m_bullet::Real, m_charge::Real) = m_bullet + m_charge / 3.0

"""
    burn_fraction(ξ; θ=0.0) -> z

Geometric form function: z(ξ) = ξ·(1 + θ·ξ).
- `ξ`: linear fraction burned (e/e₀), 0 ≤ ξ ≤ 1
- `θ`: form factor (0 = slab, >0 progressive, <0 degressive)
"""
burn_fraction(ξ::Real; θ::Real=0.0) = ξ * (1.0 + θ * ξ)

"""
    propellant_energy(charge_mass_kg, impetus, z) -> Q [J]

Total chemical energy released: Q = C·f·z.
- `charge_mass_kg`: propellant charge mass [kg]
- `impetus`: force constant f [J/kg]
- `z`: fraction burned (0–1)
"""
propellant_energy(charge_mass_kg::Real, impetus::Real, z::Real) =
    charge_mass_kg * impetus * z

"""
    chamber_pressure_closed_bomb(; C, f, z, V_ch, rho_p, eta_sp) -> P [Pa]

Closed-bomb pressure model.
P = C·f·z / (V_ch - C·(1-z)/ρ_p - C·z·η_sp)
- `C`: charge mass [kg]
- `f`: propellant impetus [J/kg]
- `z`: fraction burned
- `V_ch`: chamber volume [m³]
- `rho_p`: solid propellant density [kg/m³]
- `eta_sp`: specific co-volume [m³/kg]
"""
function chamber_pressure_closed_bomb(; C::Real, f::Real, z::Real,
                                       V_ch::Real, rho_p::Real,
                                       eta_sp::Real)
    V_gas = V_ch - C * (1.0 - z) / rho_p - C * z * eta_sp
    return C * f * z / V_gas
end

"""
    barrel_length_correction(v_ref, L_barrel, L_ref; exponent=0.27) -> v [m/s]

Correct muzzle velocity for barrel length difference.
v = v_ref · (L_barrel / L_ref)^exponent
Typical exponent: 0.25–0.30.
"""
function barrel_length_correction(v_ref::Real, L_barrel::Real,
                                  L_ref::Real; exponent::Real=0.27)
    return v_ref * (L_barrel / L_ref)^exponent
end

"""
    muzzle_velocity_temp_correction(v0, ΔT; σ_T=1.0) -> v [fps or m/s]

Temperature correction: v = v0 + σ_T · ΔT.
- `v0`: base muzzle velocity
- `ΔT`: temperature deviation from reference
- `σ_T`: temperature sensitivity (default 1.0 fps/°F; ~0.5 for Extreme powders)
"""
function muzzle_velocity_temp_correction(v0::Real, dT::Real;
                                         sigma_T::Real=1.0)
    return v0 + sigma_T * dT
end

end # module InteriorBallistics


# ──────────────────────────────────────────────────────────────────────────────
# MODULE 4: DragModels
# ──────────────────────────────────────────────────────────────────────────────
module DragModels

# Exports are at the bottom of this module after all definitions

# G7 BRL tabular data: [Mach  Cd]
const G7_TABLE = Float64[
    0.000 0.1198;  0.050 0.1197;  0.100 0.1196;
    0.150 0.1194;  0.200 0.1193;  0.250 0.1194;
    0.300 0.1194;  0.350 0.1194;  0.400 0.1193;
    0.450 0.1193;  0.500 0.1194;  0.550 0.1193;
    0.600 0.1194;  0.650 0.1197;  0.700 0.1202;
    0.725 0.1207;  0.750 0.1215;  0.775 0.1226;
    0.800 0.1242;  0.825 0.1266;  0.850 0.1306;
    0.875 0.1368;  0.900 0.1464;  0.925 0.1660;
    0.950 0.2054;  0.975 0.2993;  1.000 0.3803;
    1.025 0.4015;  1.050 0.3845;  1.075 0.3710;
    1.100 0.3597;  1.125 0.3497;  1.150 0.3405;
    1.200 0.3250;  1.250 0.3131;  1.300 0.2992;
    1.350 0.2880;  1.400 0.2778;  1.450 0.2686;
    1.500 0.2602;  1.550 0.2525;  1.600 0.2452;
    1.650 0.2383;  1.700 0.2321;  1.750 0.2261;
    1.800 0.2204;  1.850 0.2147;  1.900 0.2093;
    1.950 0.2042;  2.000 0.1993;  2.050 0.1947;
    2.100 0.1905;  2.150 0.1866;  2.200 0.1830;
    2.250 0.1796;  2.300 0.1763;  2.350 0.1729;
    2.400 0.1695;  2.450 0.1660;  2.500 0.1629;
    3.000 0.1400;  3.500 0.1225;  4.000 0.1090;
    4.500 0.0980;  5.000 0.0893
]

# G1 BRL tabular data: [Mach  Cd]
const G1_TABLE = Float64[
    0.000 0.2629;  0.050 0.2558;  0.100 0.2487;
    0.150 0.2413;  0.200 0.2344;  0.250 0.2278;
    0.300 0.2214;  0.350 0.2155;  0.400 0.2104;
    0.450 0.2061;  0.500 0.2032;  0.550 0.2020;
    0.600 0.2034;  0.650 0.2085;  0.700 0.2165;
    0.750 0.2230;  0.800 0.2313;  0.850 0.2417;
    0.875 0.2487;  0.900 0.2558;  0.925 0.2705;
    0.950 0.2939;  0.975 0.3200;  1.000 0.4528;
    1.025 0.4748;  1.050 0.4888;  1.075 0.4951;
    1.100 0.4992;  1.125 0.4973;  1.150 0.4950;
    1.200 0.4790;  1.250 0.4621;  1.300 0.4493;
    1.350 0.4369;  1.400 0.4253;  1.450 0.4145;
    1.500 0.4042;  1.550 0.3945;  1.600 0.3855;
    1.650 0.3769;  1.700 0.3687;  1.750 0.3608;
    1.800 0.3532;  1.850 0.3460;  1.900 0.3392;
    1.950 0.3326;  2.000 0.3264;  2.500 0.2756;
    3.000 0.2394;  3.500 0.2126;  4.000 0.1923;
    4.500 0.1765;  5.000 0.1637
]

"""
    interp_table(table, mach) -> Cd

Linear interpolation of a [Mach Cd] drag table.
Clamps at table boundaries.
"""
function interp_table(table::Matrix{Float64}, mach::Real)
    M = @view table[:, 1]
    C = @view table[:, 2]
    n = length(M)
    if mach <= M[1]
        return C[1]
    elseif mach >= M[n]
        return C[n]
    end
    idx = searchsortedlast(M, Float64(mach))
    idx = clamp(idx, 1, n - 1)
    frac = (mach - M[idx]) / (M[idx+1] - M[idx])
    return C[idx] + frac * (C[idx+1] - C[idx])
end

# ── Cubic spline interpolation (C²-continuous) ──

"""
    CubicSpline

Precomputed natural cubic spline through (x, y) data.
Coefficients a, b, c, d for each interval:
  S_k(x) = a_k + b_k*(x-x_k) + c_k*(x-x_k)^2 + d_k*(x-x_k)^3
"""
struct CubicSpline
    x::Vector{Float64}
    a::Vector{Float64}   # = y values
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
end

"""
    build_cubic_spline(x, y) -> CubicSpline

Construct a natural cubic spline (S''=0 at endpoints) through
the data points (x_i, y_i). Solves the tridiagonal system in O(n).
"""
function build_cubic_spline(x::AbstractVector{<:Real},
                             y::AbstractVector{<:Real})
    n = length(x)
    @assert n == length(y) && n >= 3
    h = diff(Float64.(x))
    a = Float64.(y)

    # Build tridiagonal system for c coefficients
    # Interior equations: h_{k-1} c_{k-1} + 2(h_{k-1}+h_k) c_k + h_k c_{k+1}
    #                     = 3[(a_{k+1}-a_k)/h_k - (a_k-a_{k-1})/h_{k-1}]
    alpha = zeros(n)
    for i in 2:n-1
        alpha[i] = 3.0 * ((a[i+1] - a[i]) / h[i] - (a[i] - a[i-1]) / h[i-1])
    end

    # Solve with Thomas algorithm (tridiagonal)
    l = ones(n)
    mu = zeros(n)
    z = zeros(n)

    for i in 2:n-1
        l[i] = 2.0 * (x[i+1] - x[i-1]) - h[i-1] * mu[i-1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i]
    end

    c = zeros(n)
    b = zeros(n - 1)
    d = zeros(n - 1)

    for j in (n-1):-1:1
        c[j] = z[j] - mu[j] * c[j+1]
        b[j] = (a[j+1] - a[j]) / h[j] - h[j] * (c[j+1] + 2.0 * c[j]) / 3.0
        d[j] = (c[j+1] - c[j]) / (3.0 * h[j])
    end

    return CubicSpline(Float64.(x), a[1:end-1], b, c[1:end-1], d)
end

"""
    eval_spline(sp::CubicSpline, xq) -> y

Evaluate the cubic spline at query point xq.
"""
function eval_spline(sp::CubicSpline, xq::Real)
    if xq <= sp.x[1]
        return sp.a[1]
    elseif xq >= sp.x[end]
        # Evaluate last segment at its endpoint
        k = length(sp.a)
        dx = sp.x[end] - sp.x[k]
        return sp.a[k] + sp.b[k]*dx + sp.c[k]*dx^2 + sp.d[k]*dx^3
    end
    # Find interval
    k = searchsortedlast(sp.x, Float64(xq))
    k = clamp(k, 1, length(sp.a))
    dx = xq - sp.x[k]
    return sp.a[k] + sp.b[k]*dx + sp.c[k]*dx^2 + sp.d[k]*dx^3
end

# Precompute splines at module load time
const G7_SPLINE = build_cubic_spline(G7_TABLE[:, 1], G7_TABLE[:, 2])
const G1_SPLINE = build_cubic_spline(G1_TABLE[:, 1], G1_TABLE[:, 2])

"""Drag coefficient from G1 table, linear interpolation."""
cd_g1(mach::Real) = interp_table(G1_TABLE, mach)

"""Drag coefficient from G7 table, linear interpolation."""
cd_g7(mach::Real) = interp_table(G7_TABLE, mach)

"""Drag coefficient from G1 table, cubic spline interpolation (C²)."""
cd_g1_spline(mach::Real) = eval_spline(G1_SPLINE, mach)

"""Drag coefficient from G7 table, cubic spline interpolation (C²)."""
cd_g7_spline(mach::Real) = eval_spline(G7_SPLINE, mach)

"""
    drag_coefficient(model, mach; spline=true) -> Cd

Return drag coefficient for `:G1` or `:G7` at `mach`, by C²-continuous cubic
spline. Pass `spline=false` for raw linear interpolation of the table.

!!! warning "Le défaut était `false` jusqu'au 2026-08-14"
    Le portage JS a `spline = true` par défaut, et le README de la bibliothèque
    annonce l'interpolation « C²-continuous cubic spline ». Le Julia interpolait
    donc LINÉAIREMENT, en contradiction avec sa propre documentation et avec le
    code que le site sert. Aux nœuds de la table les deux coïncident — c'est
    pourquoi un balayage grossier n'y voyait rien —, mais entre deux nœuds l'écart
    atteint 0,6 % : à Mach 2,866 le linéaire rend 0,14612 et la spline 0,14524.

    C'était la source unique de TOUT le résidu observé entre les deux solveurs 3-DOF,
    dès le premier pas d'intégration. Trouvé par le contrôle de parité Julia ↔ JS.
"""
function drag_coefficient(model::Symbol, mach::Real; spline::Bool=true)
    if model == :G1
        return spline ? cd_g1_spline(mach) : cd_g1(mach)
    elseif model == :G7
        return spline ? cd_g7_spline(mach) : cd_g7(mach)
    else
        error("Unknown drag model: $model. Supported: :G1, :G7")
    end
end

"""
    CustomDragModel

A custom drag model from Doppler radar data or CFD.
Stores the raw (Mach, Cd) table and a precomputed cubic spline.
"""
struct CustomDragModel
    table::Matrix{Float64}
    spline::CubicSpline
    name::String
end

"""
    build_custom_drag(mach_vec, cd_vec; name="Custom") -> CustomDragModel

Build a custom drag model from user-supplied (Mach, Cd) data pairs.
The data is sorted by Mach number and a cubic spline is fitted.

# Example: Doppler-radar measured data for a specific bullet
```julia
cdm = build_custom_drag(
    [0.0, 0.5, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.5, 2.0, 2.5],
    [0.11, 0.11, 0.12, 0.15, 0.22, 0.38, 0.37, 0.34, 0.31, 0.27, 0.23, 0.20],
    name="Berger 140 Hybrid CDM"
)
cd_at_mach = cdm_eval(cdm, 1.8)
```
"""
function build_custom_drag(mach_vec::AbstractVector{<:Real},
                            cd_vec::AbstractVector{<:Real};
                            name::String="Custom")
    @assert length(mach_vec) == length(cd_vec) >= 3
    idx = sortperm(mach_vec)
    M = Float64.(mach_vec[idx])
    C = Float64.(cd_vec[idx])
    table = hcat(M, C)
    sp = build_cubic_spline(M, C)
    return CustomDragModel(table, sp, name)
end

"""Evaluate a custom drag model at a given Mach number (spline)."""
cdm_eval(cdm::CustomDragModel, mach::Real) = eval_spline(cdm.spline, mach)

"""Evaluate a custom drag model at a given Mach number (linear fallback)."""
cdm_eval_linear(cdm::CustomDragModel, mach::Real) = interp_table(cdm.table, mach)

"""
    cd_from_radar(velocity, dv_dt, mass_kg, rho, A) -> Cd

Extract drag coefficient from Doppler radar velocity and deceleration data.
Cd = -2·m·(dv/dt) / (ρ·A·v²)
"""
function cd_from_radar(velocity::Real, dv_dt::Real,
                       mass_kg::Real, rho::Real, A::Real)
    return -2.0 * mass_kg * dv_dt / (rho * A * velocity^2)
end

export cd_g1, cd_g7, cd_g1_spline, cd_g7_spline,
       drag_coefficient, G1_TABLE, G7_TABLE,
       CubicSpline, build_cubic_spline, eval_spline,
       CustomDragModel, build_custom_drag, cdm_eval, cdm_eval_linear,
       cd_from_radar

end # module DragModels


# ──────────────────────────────────────────────────────────────────────────────
# MODULE 5: ExteriorBallistics
# ──────────────────────────────────────────────────────────────────────────────
module ExteriorBallistics

using ..Atmosphere
using ..DragModels
using ..BallisticUtils
using Printf

export TrajectoryPoint, ShotParameters,
       solve_trajectory, find_zero_angle,
       miller_sg, spin_drift_inches,
       coriolis_horizontal, eotvos_vertical,
       wind_deflection_lag_rule,
       aerodynamic_jump, inclined_fire_effective_range,
       trajectory_table

# ── Data types ──

"""Result at a single point along the trajectory."""
struct TrajectoryPoint
    time::Float64           # s
    range_m::Float64        # m (downrange)
    drop_m::Float64         # m (vertical, relative to line of sight)
    windage_m::Float64      # m (cross-range, positive = right) — TOTAL
    spin_drift_m::Float64   # m, part of windage_m due to spin drift alone
    vx::Float64             # m/s
    vy::Float64             # m/s
    vz::Float64             # m/s
    v_total::Float64        # m/s
    mach::Float64
    energy_J::Float64       # J
end

"""All parameters defining a shot."""
Base.@kwdef struct ShotParameters
    # ── Bullet ──
    mass_grains::Float64        = 175.0     # bullet mass [gr]
    caliber_in::Float64         = 0.308     # bullet diameter [in]
    bc::Float64                 = 0.275     # ballistic coefficient
    drag_model::Symbol          = :G7       # :G1 or :G7

    # ── Launch ──
    muzzle_vel_fps::Float64     = 2600.0    # muzzle velocity [fps]
    sight_height_in::Float64    = 1.5       # scope center above bore [in]
    zero_range_yd::Float64      = 100.0     # zero range [yd]
    elevation_moa::Float64      = 0.0       # additional dial-up [MOA]
    windage_moa::Float64        = 0.0       # windage dial [MOA]

    # ── Atmosphere ──
    temp_F::Float64             = 59.0      # ambient temperature [°F]
    pressure_inhg::Float64      = 29.92     # barometric pressure [inHg]
    humidity_pct::Float64       = 0.0       # relative humidity [%]
    altitude_ft::Float64        = 0.0       # station altitude [ft]

    # ── Wind ──
    wind_speed_mph::Float64     = 0.0       # wind speed [mph]
    wind_angle_deg::Float64     = 90.0      # 0°=head, 90°=full right cross

    # ── Geometry ──
    target_range_yd::Float64    = 1000.0    # maximum range to compute [yd]
    incline_deg::Float64        = 0.0       # incline angle (+ = uphill)

    # ── Earth rotation ──
    latitude_deg::Float64       = 45.0      # firing latitude [°]
    azimuth_deg::Float64        = 0.0       # firing direction from N [°]
    enable_coriolis::Bool       = true

    # ── Spin ──
    twist_in::Float64           = 10.0      # twist rate [in/turn]
    twist_direction::Int        = 1         # +1=right, -1=left
    bullet_length_in::Float64   = 1.24      # bullet length [in]
    enable_spin_drift::Bool     = true

    # ── Solver ──
    dt::Float64                 = 0.0005    # integration time step [s]
end

# ── Helper: convert to SI ──
function _to_si(p::ShotParameters)
    (
        m   = grains_to_kg(p.mass_grains),
        d   = inches_to_m(p.caliber_in),
        v0  = fps_to_ms(p.muzzle_vel_fps),
        T   = fahrenheit_to_kelvin(p.temp_F),
        P   = inhg_to_pa(p.pressure_inhg),
        H   = p.humidity_pct,
        w   = p.wind_speed_mph * 0.44704,
        sh  = inches_to_m(p.sight_height_in),
        zr  = yards_to_m(p.zero_range_yd),
        tr  = yards_to_m(p.target_range_yd),
        L   = inches_to_m(p.bullet_length_in),
        tw  = inches_to_m(p.twist_in),
    )
end

"""
    miller_sg(p::ShotParameters) -> Sg

Gyroscopic stability factor (Miller twist rule).
"""
function miller_sg(p::ShotParameters)
    return miller_stability(
        mass_gr         = p.mass_grains,
        caliber_in      = p.caliber_in,
        bullet_length_in = p.bullet_length_in,
        twist_in        = p.twist_in,
        muzzle_vel_fps  = p.muzzle_vel_fps,
        temp_F          = p.temp_F,
        pressure_inhg   = p.pressure_inhg,
    )
end

"""
    spin_drift_inches(t, Sg, direction) -> Δz [inches]

Litz empirical spin-drift formula: Δz = 1.25·(Sg+1.2)·t^1.83.
Direction: +1 for right-hand twist, -1 for left-hand twist.
"""
function spin_drift_inches(t::Real, sg::Real, direction::Int=1)
    return direction * 1.25 * (sg + 1.2) * t^1.83
end

"""
    coriolis_horizontal(range_m, tof, latitude_deg) -> Δz [m]

Horizontal Coriolis deflection. Positive = right (Northern Hemisphere).
Δz ≈ ω_E · sin(φ) · R · t
"""
function coriolis_horizontal(range_m::Real, tof::Real, latitude_deg::Real)
    ω = 7.2921e-5
    return ω * sin(deg2rad(latitude_deg)) * range_m * tof
end

"""
    eotvos_vertical(range_m, tof, latitude_deg, azimuth_deg) -> Δy [m]

Eötvös vertical deflection, **positive up** — the same convention as the solver's
`drop_m`.

    Δy ≈ +ω_E · cos(φ) · sin(ψ) · R · t

ψ is the firing azimuth clockwise from North, so ψ = 90° is due East: firing east
lands **high**, firing west lands low. That is the sense the full solver produces
when its Coriolis term is switched on — a .308 at 1000 yd and 50° latitude comes
in 0.058 m high to the east and 0.052 m low to the west.

!!! warning "This carried the opposite sign until 2026-08-14"
    Both the code and the line above read `-ω_E`, and were self-consistent, which
    is why nothing looked wrong. The JS port had the correct sign all along; the
    Julia ↔ JS parity check is what surfaced it. Neither version is called by
    anything, so the error was latent — but it was the reference implementation
    that was wrong, and it would have been copied the first time someone used it.
"""
function eotvos_vertical(range_m::Real, tof::Real,
                         latitude_deg::Real, azimuth_deg::Real)
    ω = 7.2921e-5
    return ω * cos(deg2rad(latitude_deg)) * sin(deg2rad(azimuth_deg)) *
           range_m * tof
end

"""
    wind_deflection_lag_rule(wind_cross_ms, tof, range_m, v0_ms) -> Δz [m]

Lag-rule (Didion) approximation for crosswind deflection.
Δz ≈ w_cross · (t_actual - R/v0)
"""
function wind_deflection_lag_rule(wind_cross_ms::Real, tof::Real,
                                  range_m::Real, v0_ms::Real)
    t_vac = range_m / v0_ms
    return wind_cross_ms * (tof - t_vac)
end

"""
    aerodynamic_jump(C_L_alpha, k_t, alpha_trim, p_s, r_s) -> θ_jump [rad]

Aerodynamic jump angle (Eq. 6.13 in manual).
"""
function aerodynamic_jump(C_L_alpha::Real, k_t::Real, alpha_trim::Real,
                          p_s::Real, r_s::Real)
    return (C_L_alpha / (2.0 * k_t^2)) * (alpha_trim / (p_s - r_s))
end

"""
    inclined_fire_effective_range(slant_range, angle_deg) -> R_eff

Rifleman's Rule: effective horizontal range for inclined fire (Eq. 12.1).
"""
function inclined_fire_effective_range(slant_range::Real, angle_deg::Real)
    return slant_range * cos(deg2rad(angle_deg))
end

"""
    find_zero_angle(p::ShotParameters; tol=1e-6, max_iter=50) -> θ [rad]

Iteratively find the bore elevation angle that zeros the rifle
at `p.zero_range_yd`.
"""
function find_zero_angle(p::ShotParameters; tol::Real=1e-6, max_iter::Int=50)
    si = _to_si(p)
    ρ = air_density(P=si.P, T=si.T, H=si.H)
    a_s = speed_of_sound(si.T)
    A = π * si.d^2 / 4.0
    ff = form_factor(p.mass_grains, p.caliber_in, p.bc)

    θ = atan(Atmosphere.g0 * si.zr / (2.0 * si.v0^2))  # initial guess

    for _ in 1:max_iter
        x, y, z = 0.0, -si.sh, 0.0
        vx = si.v0 * cos(θ)
        vy = si.v0 * sin(θ)
        vz = 0.0
        t  = 0.0
        dt = p.dt

        function derivs_zero(sv)
            _x, _y, _z, _vx, _vy, _vz = sv
            v = sqrt(_vx^2 + _vy^2)
            Ma = v / a_s
            cd = drag_coefficient(p.drag_model, Ma)
            df = (ρ * cd * A * ff) / (2.0 * si.m) * (4.0 / π)
            return (_vx, _vy, 0.0, -df * v * _vx, -df * v * _vy - Atmosphere.g0, 0.0)
        end

        while x < si.zr && t < 15.0
            s1 = (x, y, z, vx, vy, vz)
            k1 = derivs_zero(s1)
            s2 = s1 .+ 0.5 .* dt .* k1
            k2 = derivs_zero(s2)
            s3 = s1 .+ 0.5 .* dt .* k2
            k3 = derivs_zero(s3)
            s4 = s1 .+ dt .* k3
            k4 = derivs_zero(s4)
            next_s = s1 .+ (dt / 6.0) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
            x, y, z, vx, vy, vz = next_s
            t += dt
        end

        if abs(y) < tol
            return θ
        end
        θ -= y / si.zr * cos(θ)
    end
    return θ
end

"""
    solve_trajectory(p::ShotParameters) -> Vector{TrajectoryPoint}

Full 3-DOF trajectory solver with RK4 integration, Coriolis, and spin drift.
"""
function solve_trajectory(p::ShotParameters)
    si = _to_si(p)
    ρ  = air_density(P=si.P, T=si.T, H=si.H)
    a_s = speed_of_sound(si.T)
    A  = π * si.d^2 / 4.0
    ω_E = 7.2921e-5
    lat = deg2rad(p.latitude_deg)
    azi = deg2rad(p.azimuth_deg)

    # Earth-rotation vector in the shot frame (x = downrange/LOS, y = up, z = right).
    # Ω in (East, North, Up) = ω(0, cos lat, sin lat); projected onto a shot frame
    # whose downrange axis points at azimuth `azi` from North:
    #   downrange = (sin azi, cos azi, 0), up = (0,0,1), right = (cos azi, -sin azi, 0).
    #
    # Missing alongside `g_along`/`g_perp` since 4c2343c — same cause, same silence.
    oe_x =  ω_E * cos(lat) * cos(azi)
    oe_y =  ω_E * sin(lat)
    oe_z = -ω_E * cos(lat) * sin(azi)

    # Wind components
    w_ang = deg2rad(p.wind_angle_deg)
    wx = -si.w * cos(w_ang)
    wz =  si.w * sin(w_ang)
    wy = 0.0

    # Incline: resolve gravity in the LOS frame. The component perpendicular to the
    # line of sight (g·cos α) drives the drop — the physical basis of the
    # "rifleman's rule" — while g·sin α acts along the flight path. α = 0 ⇒ flat fire.
    #
    # These two were MISSING and `solve_trajectory` threw `UndefVarError: g_along`
    # on every call, from commit 4c2343c until the Julia ↔ JS parity check found it
    # on 2026-08-14. The JS port carried them all along; nothing compared the two.
    inc    = deg2rad(p.incline_deg)
    g_along = -Atmosphere.g0 * sin(inc)
    g_perp  = -Atmosphere.g0 * cos(inc)

    # Zero angle + dialed elevation/windage
    θ0 = find_zero_angle(p)
    θ0 += moa_to_rad(p.elevation_moa)
    ψ0  = moa_to_rad(p.windage_moa)

    sg = miller_sg(p)
    dt = p.dt
    ff = form_factor(p.mass_grains, p.caliber_in, p.bc)

    # Initial state
    x, y, z = 0.0, -si.sh, 0.0
    # Pas de cos(inc) ici : la vitesse initiale est portée par l'axe du canon, qui
    # pointe θ0 au-dessus de la ligne de visée quelle que soit l'inclinaison de
    # l'ensemble. Incliner l'arme ne ralentit pas la balle — le dévers n'entre que
    # par la décomposition de la gravité (g_along / g_perp) posée plus haut.
    # Le facteur cos(inc) était présent jusqu'au 2026-08-14 : à 30° il retirait 13 %
    # de vitesse initiale et 30 % sur la chute à 600 yd. Trouvé par la parité avec le JS.
    vx = si.v0 * cos(θ0) * cos(ψ0)
    vy = si.v0 * sin(θ0)
    vz = si.v0 * cos(θ0) * sin(ψ0)
    t  = 0.0

    results = TrajectoryPoint[]

    # ── RK4 derivatives ──
    function derivs(sv)
        _x, _y, _z, _vx, _vy, _vz = sv
        vrx = _vx - wx;  vry = _vy - wy;  vrz = _vz - wz
        vrel = sqrt(vrx^2 + vry^2 + vrz^2)
        Ma = vrel / a_s
        cd = drag_coefficient(p.drag_model, Ma)
        D = (ρ * cd * A * ff) / (2.0 * si.m) * (4.0 / π)

        ax = -D * vrel * vrx + g_along
        ay = -D * vrel * vry + g_perp
        az = -D * vrel * vrz

        if p.enable_coriolis
            ax += -2.0 * (oe_y * _vz - oe_z * _vy)
            ay += -2.0 * (oe_z * _vx - oe_x * _vz)
            az += -2.0 * (oe_x * _vy - oe_y * _vx)
        end
        return (_vx, _vy, _vz, ax, ay, az)
    end

    while x <= si.tr && t < 15.0
        v_tot = sqrt(vx^2 + vy^2 + vz^2)
        mach  = v_tot / a_s
        ek    = 0.5 * si.m * v_tot^2

        push!(results, TrajectoryPoint(
            t, x, y, z, 0.0,          # spin_drift_m rempli après l'intégration
            vx, vy, vz, v_tot, mach, ek
        ))

        # RK4 integration for state vector [x, y, z, vx, vy, vz]
        s1 = (x, y, z, vx, vy, vz)
        k1 = derivs(s1)

        s2 = s1 .+ 0.5 .* dt .* k1
        k2 = derivs(s2)

        s3 = s1 .+ 0.5 .* dt .* k2
        k3 = derivs(s3)

        s4 = s1 .+ dt .* k3
        k4 = derivs(s4)

        next_s = s1 .+ (dt / 6.0) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)

        x, y, z, vx, vy, vz = next_s
        t += dt
    end

    # Post-hoc spin drift
    if p.enable_spin_drift
        for i in eachindex(results)
            pt = results[i]
            sd_m = spin_drift_inches(pt.time, sg, p.twist_direction) * 0.0254
            results[i] = TrajectoryPoint(
                pt.time, pt.range_m, pt.drop_m,
                pt.windage_m + sd_m, sd_m,
                pt.vx, pt.vy, pt.vz,
                pt.v_total, pt.mach, pt.energy_J
            )
        end
    end

    return results
end

"""
    trajectory_table(p::ShotParameters; step_yd=100) -> formatted output

Solve trajectory and print a table at regular range intervals.
"""
function trajectory_table(p::ShotParameters; step_yd::Real=100.0)
    traj = solve_trajectory(p)
    step_m = yards_to_m(step_yd)

    println("┌────────┬──────────┬──────────┬──────────┬──────────┬──────────┬────────┬──────────┐")
    println("│Rng (yd)│ Drop (in)│ Elev(MOA)│ Wind (in)│ Spin (in)│ Vel(fps) │ ToF(s) │ Eng(ftlb)│")
    println("├────────┼──────────┼──────────┼──────────┼──────────┼──────────┼────────┼──────────┤")

    next_r = step_m
    for pt in traj
        if pt.range_m >= next_r
            r_yd   = m_to_yards(pt.range_m)
            drop_in = m_to_inches(pt.drop_m)
            wind_in = m_to_inches(pt.windage_m)
            spin_in = m_to_inches(pt.spin_drift_m)
            elev    = drop_to_moa(abs(drop_in), r_yd)
            elev    = pt.drop_m < 0 ? elev : -elev
            v_fps   = ms_to_fps(pt.v_total)
            ek_ftlb = pt.energy_J / 1.35582

            @printf("│ %6.0f │ %+8.1f │ %+8.1f │ %+8.1f │ %+8.1f │ %8.0f │ %6.3f │ %8.0f │\n",
                    r_yd, drop_in, elev, wind_in, spin_in, v_fps, pt.time, ek_ftlb)

            next_r += step_m
        end
    end
    println("└────────┴──────────┴──────────┴──────────┴──────────┴──────────┴────────┴──────────┘")
end

end # module ExteriorBallistics


# ──────────────────────────────────────────────────────────────────────────────
# MODULE 6: ReloadingAnalysis
# ──────────────────────────────────────────────────────────────────────────────
module ReloadingAnalysis

using ..BallisticUtils
using ..Atmosphere
using ..ExteriorBallistics
using Statistics

export chronograph_stats, velocity_vertical_dispersion,
       ladder_test_analysis, seating_depth_sensitivity,
       charge_weight_node, temperature_shift_table,
       bc_from_two_chronographs, optimal_barrel_time,
       load_consistency_score, pressure_estimate_psi

"""
    chronograph_stats(velocities) -> NamedTuple

Compute standard statistics from a chronograph string.
Returns (mean, es, sd, min, max, n, cv_pct).
"""
function chronograph_stats(v::AbstractVector{<:Real})
    n   = length(v)
    μ   = mean(v)
    es  = maximum(v) - minimum(v)
    σ   = n > 1 ? std(v) : 0.0
    cv  = μ > 0 ? 100.0 * σ / μ : 0.0
    return (mean=μ, es=es, sd=σ, min=minimum(v), max=maximum(v),
            n=n, cv_pct=cv)
end

"""
    velocity_vertical_dispersion(p::ShotParameters, σ_v_fps;
                                 range_yd=1000, Δv=5.0) -> σ_y [inches]

Estimate vertical group dispersion at range due to velocity SD.
Uses numerical differentiation: ∂y/∂v₀ ≈ [y(v+Δv) - y(v-Δv)] / (2Δv).
"""
function velocity_vertical_dispersion(p::ShotParameters, sigma_v_fps::Real;
                                       range_yd::Real=0.0,
                                       dv::Real=5.0)
    R = range_yd > 0 ? range_yd : p.target_range_yd
    range_m = yards_to_m(R)

    # Trajectory at v0 + Δv
    p_hi = ShotParameters(p; muzzle_vel_fps = p.muzzle_vel_fps + dv,
                          target_range_yd = R)
    traj_hi = solve_trajectory(p_hi)
    y_hi = _drop_at_range(traj_hi, range_m)

    # Trajectory at v0 - Δv
    p_lo = ShotParameters(p; muzzle_vel_fps = p.muzzle_vel_fps - dv,
                          target_range_yd = R)
    traj_lo = solve_trajectory(p_lo)
    y_lo = _drop_at_range(traj_lo, range_m)

    # Sensitivity dy/dv [m per m/s]
    dv_ms = fps_to_ms(dv)
    dy_dv = (y_hi - y_lo) / (2.0 * dv_ms)

    # Vertical σ in meters, then inches
    sigma_v_ms = fps_to_ms(sigma_v_fps)
    sigma_y_m  = abs(dy_dv) * sigma_v_ms
    return m_to_inches(sigma_y_m)
end

function _drop_at_range(traj::Vector{TrajectoryPoint}, range_m::Real)
    for i in 2:length(traj)
        if traj[i].range_m >= range_m
            # Linear interpolation
            f = (range_m - traj[i-1].range_m) / (traj[i].range_m - traj[i-1].range_m)
            return traj[i-1].drop_m + f * (traj[i].drop_m - traj[i-1].drop_m)
        end
    end
    return traj[end].drop_m
end

"""
    ladder_test_analysis(charges, velocities) -> NamedTuple

Analyse a Satterlee/ladder test: identify the flattest velocity node.
Returns charge–velocity pairs sorted by charge, plus node identification.
"""
function ladder_test_analysis(charges::AbstractVector{<:Real},
                               velocities::AbstractVector{<:Real})
    @assert length(charges) == length(velocities)
    idx = sortperm(charges)
    ch = charges[idx]
    ve = velocities[idx]
    n = length(ch)

    # Compute Δv/Δcharge slopes
    slopes = Float64[]
    for i in 2:n
        push!(slopes, (ve[i] - ve[i-1]) / (ch[i] - ch[i-1]))
    end

    # Find the minimum-slope region (velocity node)
    abs_slopes = abs.(slopes)
    best_idx = argmin(abs_slopes)
    node_charge = (ch[best_idx] + ch[best_idx+1]) / 2.0
    node_vel    = (ve[best_idx] + ve[best_idx+1]) / 2.0

    return (charges=ch, velocities=ve, slopes=slopes,
            node_charge=node_charge, node_velocity=node_vel,
            node_slope=slopes[best_idx])
end

"""
    seating_depth_sensitivity(p::ShotParameters, jump_offsets_in;
                              range_yd=100) -> Vector{NamedTuple}

Estimate the effect of seating depth changes on group vertical via the
velocity sensitivity they produce. Assumes each 0.001" of jump change
produces approximately 3–8 fps change (user-adjustable).
"""
function seating_depth_sensitivity(p::ShotParameters,
                                    jump_offsets_in::AbstractVector{<:Real};
                                    range_yd::Real=100.0,
                                    fps_per_thou::Real=5.0)
    results = []
    for offset in jump_offsets_in
        # Velocity change per offset in thousandths
        dv = offset * 1000.0 * fps_per_thou
        p_mod = ShotParameters(p; muzzle_vel_fps=p.muzzle_vel_fps + dv,
                               target_range_yd=range_yd)
        traj = solve_trajectory(p_mod)
        range_m = yards_to_m(range_yd)
        drop = _drop_at_range(traj, range_m)
        push!(results, (offset_in=offset, delta_v_fps=dv,
                        drop_in=m_to_inches(drop)))
    end
    return results
end

"""
    charge_weight_node(charges, group_sizes) -> NamedTuple

Given a vector of charge weights and corresponding group sizes,
identify the OCW (Optimal Charge Weight) node as the charge at which
group size is minimised and the point-of-impact vertical shift is flattest.
"""
function charge_weight_node(charges::AbstractVector{<:Real},
                             group_sizes::AbstractVector{<:Real})
    @assert length(charges) == length(group_sizes)
    idx = sortperm(charges)
    ch = charges[idx]
    gs = group_sizes[idx]
    best = argmin(gs)
    return (optimal_charge=ch[best], best_group=gs[best],
            charges=ch, groups=gs)
end

"""
    temperature_shift_table(p::ShotParameters; temps_F, range_yd=1000,
                            σ_T=1.0) -> Vector{NamedTuple}

Compute trajectory shift at `range_yd` across a range of temperatures,
accounting for both velocity change (propellant) and air density change.
"""
function temperature_shift_table(p::ShotParameters;
                                  temps_F::AbstractVector{<:Real}=
                                      collect(0.0:20.0:120.0),
                                  range_yd::Real=0.0,
                                  sigma_T::Real=1.0)
    R = range_yd > 0 ? range_yd : p.target_range_yd
    range_m = yards_to_m(R)
    ref_T = p.temp_F
    results = []

    for T in temps_F
        dT = T - ref_T
        v_new = p.muzzle_vel_fps + sigma_T * dT
        p_mod = ShotParameters(p;
            muzzle_vel_fps = v_new,
            temp_F         = T,
            target_range_yd = R
        )
        traj = solve_trajectory(p_mod)
        drop = _drop_at_range(traj, range_m)
        push!(results, (temp_F=T, velocity_fps=v_new,
                        drop_in=m_to_inches(drop),
                        shift_moa=drop_to_moa(m_to_inches(drop), R)))
    end
    return results
end

"""
    bc_from_two_chronographs(v1_fps, v2_fps, distance_ft;
                             caliber_in, mass_gr, drag_model=:G7,
                             temp_F=59, pressure_inhg=29.92) -> BC

Estimate BC from two chronograph readings separated by `distance_ft`.
Uses the drag equation integrated over the interval.
"""
function bc_from_two_chronographs(v1_fps::Real, v2_fps::Real,
                                   distance_ft::Real;
                                   caliber_in::Real,
                                   mass_gr::Real,
                                   drag_model::Symbol=:G7,
                                   temp_F::Real=59.0,
                                   pressure_inhg::Real=29.92)
    v1 = fps_to_ms(v1_fps)
    v2 = fps_to_ms(v2_fps)
    dist = distance_ft * 0.3048
    d = inches_to_m(caliber_in)
    m = grains_to_kg(mass_gr)
    T = fahrenheit_to_kelvin(temp_F)
    P = inhg_to_pa(pressure_inhg)
    ρ = air_density(P=P, T=T)
    A = π * d^2 / 4.0
    a_s = speed_of_sound(T)

    # Average velocity & Mach
    v_avg = (v1 + v2) / 2.0
    Ma = v_avg / a_s
    cd_std = drag_coefficient(drag_model, Ma)

    # From drag equation: dv/dx = -ρ·Cd·A·v / (2m)
    # Integrating: ln(v2/v1) = -ρ·Cd_std·A·dist / (2·BC·d²) ... via BC definition
    # BC = -ρ·Cd_std·π·dist / (8·ln(v2/v1))  ... but BC is in the denominator
    # Actually: using BC form: decel = ρ·Cd_std·π·v² / (8·BC)
    # ∫ dv/v² = -ρ·Cd_std·π/(8·BC) · ∫dx
    # 1/v1 - 1/v2 = ρ·Cd_std·π·dist / (8·BC)

    bc = ρ * cd_std * π * dist / (8.0 * (1.0/v2 - 1.0/v1))
    # Convert to standard conditions (ρ₀)
    bc *= Atmosphere.rho0 / ρ

    return abs(bc)  # BC should be positive
end

"""
    optimal_barrel_time(muzzle_vel_fps, barrel_length_in) -> t_barrel [ms]

Barrel time estimate: t ≈ 2·L / v₀ (simplified uniform acceleration).
Useful for barrel harmonic analysis.
"""
function optimal_barrel_time(muzzle_vel_fps::Real, barrel_length_in::Real)
    v0 = fps_to_ms(muzzle_vel_fps)
    L  = inches_to_m(barrel_length_in)
    # Assuming average velocity = v0/2 during barrel travel
    t = 2.0 * L / v0
    return t * 1000.0  # ms
end

"""
    load_consistency_score(es_fps, sd_fps) -> score [0–100]

A simple quality score for a load based on ES and SD.
100 = perfect; benchrest-grade loads score > 90.
"""
function load_consistency_score(es_fps::Real, sd_fps::Real)
    # Penalty: each fps of SD costs 5 points, each fps of ES costs 2 points
    score = 100.0 - 5.0 * sd_fps - 2.0 * es_fps
    return clamp(score, 0.0, 100.0)
end

"""
    pressure_estimate_psi(; charge_gr, case_capacity_gr_h2o,
                           bullet_weight_gr, bore_area_in2,
                           barrel_length_in, muzzle_vel_fps) -> P_est [psi]

Rough chamber pressure estimate from observed muzzle velocity using
energy balance: P_avg ≈ m_e · v₀² / (2 · A · L_travel).
This is an approximation; actual peak pressure is ~1.5–2× the average.
"""
function pressure_estimate_psi(;
    charge_gr::Real,
    bullet_weight_gr::Real,
    bore_area_in2::Real,
    barrel_length_in::Real,
    muzzle_vel_fps::Real
)
    m_e = (bullet_weight_gr + charge_gr / 3.0) / 7000.0  # lb
    v0  = muzzle_vel_fps  # fps
    A   = bore_area_in2   # in²
    L   = barrel_length_in # in (approximate travel)

    # Work = F·d = P_avg·A·L = ½·m_e·v²
    # But in imperial: ½·(m/g)·v² = P·A·L, with g = 32.174 ft/s² = 386.09 in/s²
    P_avg = m_e * v0^2 / (2.0 * 386.09 * A * L)

    # Peak pressure ≈ 1.6 × average for typical rifle cartridges
    P_peak = 1.6 * P_avg
    return P_peak
end

end # module ReloadingAnalysis


# ──────────────────────────────────────────────────────────────────────────────
# MODULE 7: SixDOF — Full 6-Degree-of-Freedom Rigid-Body Solver
# ──────────────────────────────────────────────────────────────────────────────
module SixDOF

using ..Atmosphere
using ..DragModels
using ..BallisticUtils
using LinearAlgebra

export AeroCoefficients6DOF, ShotParameters6DOF, State6DOF,
       solve_6dof, initial_spin_rate, moments_of_inertia,
       BulletGeometry, bullet_inertia, MCCOY_308_168_GEOMETRY, FLAT_BASE_GEOMETRY,
       gyroscopic_stability_6dof, dynamic_stability_6dof, dynamically_stable_6dof,
       mccoy_308_168_aero, mccoy_308_168_shot, mccoy_308_168_cmpa,
       mccoy_105mm_m1_aero, mccoy_105mm_m1_shot,
       brl3476_aero, brl3476_shot, brl3476_geometry, brl3476_curve,
       BRL3476_SS109, BRL3476_M855, BRL3476_L110, BRL3476_M856, BRL3476_PHYSICAL

"""
Aerodynamic coefficient set for 6-DOF simulation.

Every coefficient may be given either as a **number** — held constant — or as a
**function of Mach number**, which is how measured data actually comes. See
[`mccoy_308_168_aero`](@ref) for a set built from published spark-range
measurements, and note that the defaults below are none of that: they are
placeholders that describe no particular bullet.

!!! danger "A right solver fed placeholders is worse than a wrong one"
    Now that [`solve_6dof`](@ref) reproduces both published cases, its output
    *looks* credible whatever it is fed — and a wrong answer that looks credible
    is the more dangerous failure. Set `provenance = :measured` on any set built
    from real data; leaving the default makes the solver say once, per session,
    that it is running on invented numbers.
"""
Base.@kwdef struct AeroCoefficients6DOF
    # Cd0 as a function of Mach (zero-yaw drag); defaults to G7 × form factor
    # Spline, comme `drag_coefficient` et comme le portage JS : sur un test de
    # décimation (un nœud sur deux retiré puis reconstruit), elle retrouve les valeurs
    # enlevées avec deux fois moins d'erreur médiane que l'interpolation linéaire sur
    # G7, six fois moins sur G1. Ce défaut-ci reste de toute façon un bouche-trou,
    # signalé par `provenance` — les jeux mesurés fournissent leur propre table.
    Cd0_func::Function                = (Ma) -> cd_g7_spline(Ma)
    Cda2::Union{Float64,Function}     = 3.5     # yaw drag coefficient C_{D,δ²}
    CNa::Union{Float64,Function}      = 2.8     # normal force derivative C_{N,α} [per rad]
    CMa::Union{Float64,Function}      = 3.2     # overturning moment coeff C_{M,α} [per rad]
    Clp::Union{Float64,Function}      = -0.010  # spin-damping moment C_{l,p}
    CMq_CMad::Union{Float64,Function} = -8.0    # pitch damping sum (C_{M,q} + C_{M,α̇})
    CNpa::Union{Float64,Function}     = 0.1     # Magnus force coefficient C_{N,pα}
    CMpa::Union{Float64,Function}     = -0.2    # Magnus moment coefficient C_{M,pα}

    # `:measured` for a set built from spark-range or wind-tunnel data, whatever
    # its source. Anything else is taken as invented and warned about.
    provenance::Symbol                = :placeholder
end

"""
Normalize a coefficient to a two-argument function `(mach, alpha_t) -> value`.

Accepts a constant, a function of Mach alone, or a function of Mach and total
angle of attack (in radians). The arity is resolved **once**, when the solver
starts, so the inner loop never pays for the dispatch.
"""
_as_coef(c::Float64) = (_ma, _al) -> c
_as_coef(c::Function) = applicable(c, 1.0, 1.0) ? c : (ma, _al) -> c(ma)

"""
    BulletGeometry(; nose_frac, meplat_cal, bt_len_cal, bt_deg, jacket_cal, construction)

External shape of a spitzer bullet, plus how its mass is laid out inside. Every
length is in **calibers**, so one geometry describes a family across bore sizes.

  - `nose_frac`   — ogive length as a fraction of total bullet length.
  - `meplat_cal`  — diameter of the nose-tip flat.
  - `bt_len_cal`  — boat-tail length (0 for a flat base).
  - `bt_deg`      — boat-tail half-angle, measured from the axis.
  - `jacket_cal`  — radial jacket thickness (`:jacketed_lead` only).
  - `construction` — `:jacketed_lead` or `:homogeneous`.

The ogive is taken **tangent**: given the nose length and the meplat, the radius
follows from the tangency condition, `R = (Lₙ² + a²)/2a` with `a = (d − meplat)/2`.
A secant ogive of the same nose length carries slightly more volume forward, so a
VLD comes out marginally light in the nose here.

Defaults describe a modern boat-tailed match bullet. They are *shape* defaults:
supply the real numbers whenever a drawing or a caliper is at hand, because the
estimate is only as good as the contour it is given.

`nose_frac`, `bt_len_cal` and `bt_deg` are all medians of measured populations
(see "Where the defaults come from" below). **`meplat_cal` alone is uncalibrated** —
no source publishes it. Sweeping it over its plausible range now moves Ix by 1.9 %,
Iy by 3.2 % and S_g by only 0.9 %: that is the floor left for a default-driven
estimate. It used to be 6.3 / 4.2 / **9.7 %** when `bt_deg` was a guess too — pinning
the boat-tail angle is what removed almost all of it.

Flat-base bullets are a different family, not a variation on this one: use
[`FLAT_BASE_GEOMETRY`](@ref). Applying these boat-tail defaults to a flat base
costs 9.7 % on S_g, and always in the same direction.

# Where the defaults come from

`nose_frac` and `bt_len_cal`: medians of a published manufacturer's dimension
table, 94 boat-tailed and 13 flat-base jacketed match bullets from .17 to .375,
read 2026-08-13. Residual spread on S_g, against each bullet's own measured
contour: median absolute error 3.8 % (boat tail) and 4.1 % (flat base), with a
third of the population outside ±5 %. That spread is **irreducible by a better
constant** — a linear rule in L/d was fitted and rejected, buying 0.4 points of
median error while widening the tails.

`bt_deg`: median of 83 boat-tail cones measured by handloaders and pooled in a
community database, each recorded as two diameters and a height, so the half-angle
follows as `atan((d_major − d_minor) / 2h)`. Restricted to cones that start at full
bullet diameter, on rifle calibers, with a physically plausible angle. The value is
stable across every slice of that population — 7.5° over all 83, 7.7° over the 51
from match makers, 7.0° over the 14 Berger — so 7.5° is taken rather than the
narrower Berger-only figure, which would over-fit fourteen bullets.

Neither source is redistributed; only these three numbers are derived from them.

!!! note "Why `bt_deg` cannot be validated the way the other two were"
    The other defaults were checked by holding a bullet's own measured nose and
    boat-tail length against them. No public table publishes the boat-tail
    **angle**, so there is no population of measured inertias to score it on: it
    rests on the community measurements alone. Sanity check that it passes — with
    `bt_len_cal = 0.782` it puts the base at 0.795 calibers, which is where match
    boat tails actually sit; the previous 9.0° default implied 0.753, visibly too
    narrow.

!!! warning "The reference case is not representative"
    The .308" 168 gr Sierra International of [`MCCOY_308_168_GEOMETRY`] carries a
    0.51 caliber boat tail, which is the **2nd percentile** of the modern match
    population (median 0.782). The previous default of 0.50 matched it almost
    exactly and therefore scored well on the one case with published inertias,
    while being wrong for the bullets this library is actually used on. Do not
    re-tune these defaults on that case.

!!! note "Where the residual +2.4 % on `Iy` comes from — measured, not fixed"
    Investigated 2026-08-15. It is **not** an unmodelled tip cavity: the cavity is
    part of the `:jacketed_lead` construction, and the core height is solved from
    the mass. The residual tracks `jacket_cal`, and that parameter trades the two
    inertias against each other — a thicker shell loads the nose (`Iy` up) and
    lightens the axis (`Ix` down), so no single scalar zeroes both:

    | `jacket_cal` | | ΔIx | ΔIy | Iy/Ix |
    |---|---|---:|---:|---:|
    | 0.080 | 0.63 mm | +0.9 % | −0.5 % | 7.35 |
    | **0.085** | **0.66 mm** | **+0.6 %** | **+0.5 %** | **7.45** |
    | 0.095 *(current)* | 0.74 mm | −0.1 % | +2.4 % | 7.64 |

    Measured `Iy/Ix` is **7.45**. So 0.085 would hit exactly the ratio this
    docstring calls the column that matters, and both inertias inside 0.6 %.
    **It is deliberately NOT applied.** Gilding-metal jackets run 0.025–0.030″
    (0.081–0.098 cal), so physics does not choose between the two values, and the
    calibration set is **one bullet** — the very trap the warning above describes.
    What would settle it: a second bullet carrying both a measured inertia and a
    contour. The Sierra 190 gr HPBT MatchKing (#2210) already has the BRL inertias
    and a JBM length; only its contour is missing.
"""
Base.@kwdef struct BulletGeometry
    nose_frac::Float64    = 0.541
    meplat_cal::Float64   = 0.22
    bt_len_cal::Float64   = 0.782
    bt_deg::Float64       = 7.5
    jacket_cal::Float64   = 0.095
    construction::Symbol  = :jacketed_lead
end

"""
Shape defaults for a **flat-base** jacketed match or varmint bullet, medians of a
13-bullet population from the same table as [`BulletGeometry`](@ref)'s own
defaults.

Flat-base bullets are short (L/d median 3.34 against 4.87 for boat tails) and
carry proportionally longer noses (0.638 against 0.541), so they are a separate
family rather than a boat tail with `bt_len_cal = 0`. Reaching for the boat-tail
defaults instead biases S_g by −9.7 %.
"""
const FLAT_BASE_GEOMETRY = BulletGeometry(
    nose_frac  = 0.638,
    bt_len_cal = 0.0,
)

"""
Contour of the .308", 168 gr Sierra International, read off the dimensioned
sketch in Appendix A of McCoy chapter 9 (all cotes in calibers, 1 cal = 7.82 mm):
total length 3.98, tangent ogive of 7.00 caliber radius over 2.26, 0.25 meplat,
0.51 of boat tail at 13°. This is the reference case the estimator is checked
against — see [`bullet_inertia`](@ref).
"""
const MCCOY_308_168_GEOMETRY = BulletGeometry(
    nose_frac  = 2.26 / 3.98,
    meplat_cal = 0.25,
    bt_len_cal = 0.51,
    bt_deg     = 13.0,
    jacket_cal = 0.095,
)

"""Parameters for a 6-DOF shot."""
Base.@kwdef struct ShotParameters6DOF
    # Bullet physical
    mass_grains::Float64        = 140.0
    caliber_in::Float64         = 0.264
    bullet_length_in::Float64   = 1.34

    # Launch
    muzzle_vel_fps::Float64     = 2700.0
    twist_in::Float64           = 8.0       # inches per turn, RH positive
    twist_direction::Int        = 1         # +1 = right, -1 = left
    initial_yaw_deg::Float64    = 0.5       # initial yaw angle [deg]
    initial_pitch_rate::Float64 = 0.0       # initial q [rad/s]

    # Atmosphere
    temp_F::Float64             = 59.0
    pressure_inhg::Float64      = 29.92
    humidity_pct::Float64       = 0.0
    altitude_ft::Float64        = 0.0

    # Wind
    wind_speed_mph::Float64     = 0.0
    wind_angle_deg::Float64     = 90.0

    # Target
    target_range_yd::Float64    = 1000.0
    zero_range_yd::Float64      = 100.0

    # Aerodynamic coefficients
    aero::AeroCoefficients6DOF  = AeroCoefficients6DOF()

    # Solver
    #
    # `dt` must resolve the FAST body-frame mode, not the epicyclic period seen
    # from the ground. The gyroscopic coupling in the rate equations spins (q,r)
    # at |(Ix-Iy)/Iy|·p — about 23 000 rad/s for a match rifle bullet, i.e. a
    # 270 µs period — while the epicyclic motion an observer sees is ten to
    # twenty times slower. Sizing the step on the slow one is what left the
    # legacy default at 1e-4 s, where the solution blew up (see `solve_6dof`).
    # 5e-6 s gives ~50 steps per fast cycle, the rule of thumb McCoy reports for
    # 6-DOF work, and RK4 is then both stable and converged.
    dt::Float64                 = 5.0e-6
    integrator::Symbol          = :rk4      # :rk4, or :euler for the legacy scheme
    record_every::Int           = 20        # sample the output every N steps

    # Measured moments of inertia [kg·m²], axial and transverse. Leave at
    # `nothing` to have them estimated by `bullet_inertia` from `geometry`, which
    # holds the published contour to a few percent — but measured beats estimated,
    # so pass them when they exist. Supply both or neither.
    Ix::Union{Float64,Nothing}  = nothing
    Iy::Union{Float64,Nothing}  = nothing

    # Bullet contour used for that estimate, ignored when Ix/Iy are given.
    geometry::BulletGeometry    = BulletGeometry()

    # Quadrant elevation [deg]. When given it replaces the zeroing: the piece is
    # laid on this angle and `zero_range_yd` is ignored. Artillery, in short.
    launch_angle_deg::Union{Float64,Nothing} = nothing
    # Stop when the trajectory comes back to the muzzle horizontal, rather than at
    # `target_range_yd`. The only sensible end for a lofted shot.
    stop_on_impact::Bool        = false
    # Hard cap on flight time [s]. 15 s covers any small-arms range; a 70° shell
    # needs 80.
    max_time_s::Float64         = 15.0
end

"""Snapshot of the 6-DOF state at one time step."""
struct State6DOF
    t::Float64              # time [s]
    x::Float64; y::Float64; z::Float64    # position [m]
    u::Float64; v::Float64; w::Float64    # velocity [m/s]
    q0::Float64; q1::Float64; q2::Float64; q3::Float64  # quaternion
    p_spin::Float64         # spin rate [rad/s]
    q_pitch::Float64        # pitch rate [rad/s]
    r_yaw::Float64          # yaw rate [rad/s]
    alpha_t::Float64        # total angle of attack [rad]
    mach::Float64
    v_total::Float64
end


const _RHO_LEAD   = 11340.0   # kg/m³
const _RHO_JACKET = 8860.0    # kg/m³, gilding metal 95/5

"""Outer radius in calibers at axial station `z` (calibers from the base)."""
function _profile(g::BulletGeometry, L::Real, z::Real)
    r_body   = 0.5
    r_base   = r_body - g.bt_len_cal * tand(g.bt_deg)
    nose     = g.nose_frac * L
    z_should = L - nose
    if g.bt_len_cal > 0 && z <= g.bt_len_cal
        return r_base + (r_body - r_base) * (z / g.bt_len_cal)
    elseif z <= z_should
        return r_body
    else
        a = r_body - g.meplat_cal / 2          # radius drop over the ogive
        R = (nose^2 + a^2) / (2a)              # tangency condition
        return -(R - r_body) + sqrt(max(R^2 - (z - z_should)^2, 0.0))
    end
end

"""Simpson's rule for `f` on `[a, b]` with `n` (even) subintervals."""
function _simpson(f, a::Real, b::Real, n::Int)
    b <= a && return 0.0
    isodd(n) && (n += 1)
    h = (b - a) / n
    s = f(a) + f(b)
    for i in 1:(n - 1)
        s += (isodd(i) ? 4.0 : 2.0) * f(a + i * h)
    end
    return s * h / 3
end

"""
Integrate `f(z)` over the whole bullet, splitting at every slope or density
break so Simpson never straddles a kink.
"""
function _integrate(f, g::BulletGeometry, L::Real, breaks::Vector{Float64})
    bounds = sort(unique(clamp.(vcat(0.0, g.bt_len_cal, (1 - g.nose_frac) * L,
                                     L, breaks), 0.0, L)))
    total = 0.0
    for i in 1:(length(bounds) - 1)
        total += _simpson(f, bounds[i], bounds[i + 1], 400)
    end
    return total
end

"""
    bullet_inertia(mass_kg, caliber_m, length_m; geometry) -> (Ix, Iy, cg_m, density)

Moments of inertia of a bullet from its **shape**, by integrating the solid of
revolution rather than pretending it is a cylinder. `Ix` is axial, `Iy`
transverse about the centre of gravity, `cg_m` is measured from the base, and
`density` is the mean density the model implies — worth reading as a check, since
a value far below the constituent materials means the real projectile carries a
cavity the contour does not describe.

For `:jacketed_lead` the bullet is a gilding-metal shell of constant radial
thickness with an open base, a lead core rising from that base, and a cavity
above it — the ordinary hollow-point match construction. **The core height is not
a free parameter**: it is solved so the total comes to the given mass. Shape and
mass together therefore determine the internal layout, and nothing is fitted to
the inertias themselves.

For `:homogeneous` — a monolithic turned-copper or solid-steel projectile — the
density is simply mass over volume.

If the mass falls outside what the jacketed model can reach (denser than a
solid-lead fill, or lighter than the bare jacket), it silently falls back to a
uniform density, which reproduces the mass by construction.

# Accuracy
Against the .308", 168 gr Sierra International of McCoy Table 9.2 — the one
bullet here with measured values — using its published contour
([`MCCOY_308_168_GEOMETRY`](@ref)):

| | measured | this function | solid cylinder |
|---|---:|---:|---:|
| Ix | 7.23e-8 kg·m² | +0.4 % | +15 % |
| Iy | 5.38e-7 kg·m² | +2 % | **+71 %** |
| Iy/Ix | 7.44 | +2 % | +49 % |
| CG from base | 0.474 in | −1 % | — |

The `Iy/Ix` ratio is what governs both the gyroscopic coupling in the rate
equations and the yaw of repose, so that is the column that matters. The
integrator itself was checked separately on McCoy's Example 12.1, a homogeneous
20 mm steel cone-cylinder of exactly known geometry: CG −0.4 %, Ix −0.5 %,
Iy −4 %.

Measured values still beat any estimate — pass them through
`ShotParameters6DOF(Ix=…, Iy=…)` when they exist.
"""
function bullet_inertia(mass_kg::Real, caliber_m::Real, length_m::Real;
                        geometry::BulletGeometry = BulletGeometry())
    g = geometry
    L = length_m / caliber_m                    # total length, in calibers
    d3 = caliber_m^3                            # calibers³ → m³
    r(z) = _profile(g, L, z)

    # Volume-based quantities, all in caliber units; density restores the scale.
    vol_of(rho) = π * _integrate(z -> rho(z) * r(z)^2, g, L, Float64[])

    rho_body, cg_cal, Ix_cal, Iy_cal = if g.construction === :jacketed_lead
        t = g.jacket_cal
        # Mass carried below `zc`, as a function of core height.
        mass_at(zc) = π * d3 * _integrate(g, L, [zc]) do z
            ro = r(z); ri = max(ro - t, 0.0)
            (ro^2 - ri^2) * _RHO_JACKET + ri^2 * (z <= zc ? _RHO_LEAD : 0.0)
        end
        m_hollow, m_full = mass_at(0.0), mass_at(L)
        if mass_kg <= m_hollow || mass_kg >= m_full
            (nothing, nothing, nothing, nothing)          # → uniform fallback
        else
            lo, hi = 0.0, L
            for _ in 1:60
                mid = (lo + hi) / 2
                mass_at(mid) < mass_kg ? (lo = mid) : (hi = mid)
            end
            zc = (lo + hi) / 2
            # Linear mass density and dIx/dz, both per caliber of length.
            lam(z) = (ro = r(z); ri = max(ro - t, 0.0);
                      π * ((ro^2 - ri^2) * _RHO_JACKET +
                           ri^2 * (z <= zc ? _RHO_LEAD : 0.0)))
            axi(z) = (ro = r(z); ri = max(ro - t, 0.0);
                      (π / 2) * ((ro^4 - ri^4) * _RHO_JACKET +
                                 ri^4 * (z <= zc ? _RHO_LEAD : 0.0)))
            m  = _integrate(lam, g, L, [zc])
            cg = _integrate(z -> lam(z) * z, g, L, [zc]) / m
            Ia = _integrate(axi, g, L, [zc])
            It = _integrate(z -> axi(z) / 2 + lam(z) * (z - cg)^2, g, L, [zc])
            (m, cg, Ia, It)
        end
    else
        (nothing, nothing, nothing, nothing)
    end

    if rho_body === nothing                      # uniform density
        v  = _integrate(z -> r(z)^2, g, L, Float64[])
        cg = _integrate(z -> r(z)^2 * z, g, L, Float64[]) / v
        Ia = _integrate(z -> r(z)^4 / 2, g, L, Float64[])
        It = _integrate(z -> r(z)^4 / 4 + r(z)^2 * (z - cg)^2, g, L, Float64[])
        rho = mass_kg / (π * v * d3)
        return (π * rho * Ia * d3 * caliber_m^2, π * rho * It * d3 * caliber_m^2,
                cg * caliber_m, rho)
    end

    scale = d3 * caliber_m^2                     # calibers⁵ → m⁵
    return (Ix_cal * scale, Iy_cal * scale, cg_cal * caliber_m,
            rho_body * d3 / (π * _integrate(z -> r(z)^2, g, L, Float64[]) * d3))
end

"""
    moments_of_inertia(mass_kg, caliber_m, length_m; geometry) -> (Ix, Iy)

Axial and transverse moments of inertia of a bullet, from
[`bullet_inertia`](@ref) with the default match-bullet contour. Kept as a
two-value shorthand for the call sites that only need the pair.

!!! note "The default contour is a guess, the mass is not"
    Without a drawing, `nose_frac`, the boat tail and the jacket thickness are
    assumptions; only mass, caliber and length are known. Passing a measured
    contour through `geometry` tightens the answer considerably, and measured
    inertias via `ShotParameters6DOF(Ix=…, Iy=…)` settle it outright.
"""
function moments_of_inertia(mass_kg::Real, caliber_m::Real, length_m::Real;
                            geometry::BulletGeometry = BulletGeometry())
    Ix, Iy, _, _ = bullet_inertia(mass_kg, caliber_m, length_m; geometry = geometry)
    return (Ix, Iy)
end

"""
    initial_spin_rate(v0_ms, twist_m) -> p [rad/s]

Spin rate from muzzle velocity and twist rate.
p = 2π·v₀ / twist
"""
function initial_spin_rate(v0_ms::Real, twist_m::Real)
    return 2π * v0_ms / twist_m
end

"""
    gyroscopic_stability_6dof(Ix, p, rho, A, d, Iy, v, CMa) -> Sg

Exact gyroscopic stability factor from 6-DOF theory.
"""
function gyroscopic_stability_6dof(Ix::Real, p::Real, rho::Real, A::Real,
                                    d::Real, Iy::Real, v::Real, CMa::Real)
    return Ix^2 * p^2 / (2.0 * rho * A * d * Iy * v^2 * CMa)
end

"""
    dynamic_stability_6dof(CLa, CD, CMpa, CMq_CMad, kx_inv2, ky_inv2) -> Sd

Dynamic stability factor of the linearized pitching and yawing motion,

    Sd = 2T/H,   T = C_Lα + k_x⁻² C_Mpα,   H = C_Lα − C_D − k_y⁻² (C_Mq + C_Mα̇)

where `kx_inv2` = m d²/Ix and `ky_inv2` = m d²/Iy are the inverse squared radii
of gyration. McCoy states T, H and Sd in *starred* coefficients — each multiplied
by the relative density factor ρSd/2m — but that factor is common to numerator and
denominator and cancels, so plain coefficients are used here.

`Sd` is of order one for a real projectile. Combine it with the gyroscopic
stability factor through [`dynamically_stable_6dof`](@ref); neither number decides
alone, and a gyroscopically stable projectile may still be dynamically unstable.

!!! warning "Feed it the cycle-averaged `C_Mpα`, not the zero-yaw value"
    The Magnus moment coefficient is the one coefficient that goes markedly
    nonlinear within the first couple of degrees of yaw, and `Sd` is a property of
    a *motion*, not of an instant — so it wants `C_Mpα` averaged over the yaw
    cycle being described. McCoy works the reference case at an average of −0.22
    where the zero-yaw table value is −0.33, and the difference is not cosmetic:
    the average gives `Sd = +0.1147`, the zero-yaw value `Sd = −0.047`. Those are
    two different verdicts. A positive `Sd` says enough spin can stabilize the
    motion — here `Sg > 4.63`, a 1:7 twist — while a negative one says no amount
    of spin ever will. Reproduced in `validation_mccoy_1011.jl`.

    The zero-yaw value is the right input for the *limit-cycle* formulas of §13.7
    instead, where `λ_S = λ_S0 + λ_S2 δ²` carries the yaw dependence explicitly.

!!! note "What this replaced"
    Until 2026-08-12 this function returned an expression with no counterpart in the
    literature: it added `2m·kt2/(ρAd)` — a quantity of order 10⁵ — to `C_Lα`, and
    divided by `Sg(Sg−1)`, which belongs to the stability *criterion*, not to the
    factor. With the module's own defaults it returned 7.2e8 where the answer is
    0.31. Nothing in the repository called it, so no published figure depended on it.
"""
function dynamic_stability_6dof(CLa::Real, CD::Real, CMpa::Real,
                                 CMq_CMad::Real, kx_inv2::Real, ky_inv2::Real)
    T = CLa + kx_inv2 * CMpa
    H = CLa - CD - ky_inv2 * CMq_CMad
    return 2T / H
end

"""
    dynamic_stability_6dof(aero, mach, m, d, Ix, Iy) -> Sd

Convenience method taking the module's own coefficient set. The lift-curve slope
is recovered from the normal-force slope as `C_Lα = C_Nα − C_D`, and the drag is
read from `aero.Cd0_func` at the given Mach number.
"""
function dynamic_stability_6dof(aero::AeroCoefficients6DOF, mach::Real,
                                 m::Real, d::Real, Ix::Real, Iy::Real)
    CD = aero.Cd0_func(mach)
    return dynamic_stability_6dof(aero.CNa - CD, CD, aero.CMpa, aero.CMq_CMad,
                                  m * d^2 / Ix, m * d^2 / Iy)
end

"""
    dynamically_stable_6dof(Sg, Sd) -> Bool

The generalized dynamic stability criterion, `1/Sg < Sd(2 − Sd)`. It is symmetric
about `Sd = 1`: a projectile is dynamically unstable when `Sd` strays far enough
from unity in either direction, and a larger `Sg` widens the tolerated band. Note
that `Sd ≤ 0` or `Sd ≥ 2` cannot be rescued by any amount of spin.
"""
dynamically_stable_6dof(Sg::Real, Sd::Real) = 1.0 / Sg < Sd * (2.0 - Sd)

# ─────────────────────────────────────────────────────────────────────────────
# Reference case: .308", 168 gr Sierra International — McCoy, Example 9.1
#
# Spark-range measurements tabulated in Appendix A of the 6-DOF chapter, and
# physical characteristics from Table 9.2 of the same chapter. This is the one
# fully documented case available to validate the module end to end: the source
# publishes both the inputs and the resulting 6-DOF flight, so a run here can be
# checked rather than merely admired.
#
# Not transcribed, and therefore left at library defaults:
#   * C_Npα (Magnus FORCE) — not tabulated in the appendix for this bullet;
#   * nothing else — the yaw dependence of C_D, C_Mα and C_Mpα is carried since
#     2026-08-12, the last being the mechanism that sets the amplitude at all.
# ─────────────────────────────────────────────────────────────────────────────

const MCCOY_308_168_CD0 = [
    0.0 .140; 0.8 .140; 0.85 .142; 0.90 .160; 0.95 .240; 1.00 .430; 1.05 .449;
    1.1 .447; 1.2 .434; 1.4 .410; 1.6 .385; 1.8 .365; 2.0 .350; 2.2 .339; 2.5 .320]

const MCCOY_308_168_CDD2 = [
    0.0 2.9; 0.95 2.9; 1.0 3.0; 1.05 3.1; 1.1 3.6; 1.2 6.5; 1.4 7.6; 1.6 7.3;
    1.8 6.8; 2.0 6.1; 2.2 5.4; 2.5 4.4]

const MCCOY_308_168_CLP = [
    0.0 -.0150; 0.5 -.0125; 0.8 -.0108; 0.85 -.0107; 0.90 -.0105; 0.95 -.0103;
    1.00 -.0100; 1.05 -.0099; 1.1 -.0098; 1.2 -.0095; 1.4 -.0088; 1.6 -.0083;
    1.8 -.0080; 2.0 -.0075; 2.2 -.0073; 2.5 -.0068]

const MCCOY_308_168_CLA = [
    0.0 1.75; 0.5 1.63; 0.8 1.45; 0.85 1.40; 0.90 1.35; 0.95 1.30; 1.0 1.35;
    1.05 1.55; 1.1 1.70; 1.2 1.90; 1.4 2.15; 1.6 2.32; 1.8 2.45; 2.0 2.58;
    2.2 2.68; 2.5 2.85]

const MCCOY_308_168_CMA0 = [
    0.0 3.05; 0.5 3.26; 0.8 3.38; 0.85 3.40; 0.90 3.43; 0.95 3.45; 1.0 3.24;
    1.05 3.17; 1.1 3.15; 1.2 3.12; 1.4 3.06; 1.6 2.98; 1.8 2.88; 2.0 2.79;
    2.2 2.69; 2.5 2.56]

# Pitch damping sum. Positive — i.e. ANTI-damping — below Mach 1.05, as measured.
const MCCOY_308_168_CMQ = [
    0.0 1.2; 1.05 1.2; 1.1 0.5; 1.2 -3.6; 1.4 -7.3; 1.6 -8.2; 2.5 -8.2]

# Magnus moment. The appendix tabulates it against BOTH Mach and α_t², in square
# DEGREES, and the yaw dependence is the whole story for this bullet: the
# coefficient changes sign a little above two degrees of yaw. Each row is
# (Mach, α_t² breakpoints, values); the breakpoint itself moves with Mach.
const MCCOY_308_168_CMPA_TABLE = [
    (0.00, [0.0, 29.2, 400.0], [-2.60, 0.06, 0.06]),
    (0.90, [0.0, 29.2, 400.0], [-2.60, 0.06, 0.06]),
    (1.10, [0.0, 18.4, 400.0], [-1.35, 0.05, 0.05]),
    (1.40, [0.0,  9.9, 400.0], [-0.51, 0.24, 0.24]),
    (1.70, [0.0,  5.6, 400.0], [-0.33, 0.10, 0.10]),
    (2.50, [0.0,  5.6, 400.0], [-0.33, 0.10, 0.10]),
]

"""Piecewise-linear lookup in a vector of breakpoints, clamped at both ends."""
function _interp1(xs::Vector{Float64}, ys::Vector{Float64}, x::Real)
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, Float64(x))
    i = clamp(i, 1, length(xs) - 1)
    return ys[i] + (x - xs[i]) / (xs[i+1] - xs[i]) * (ys[i+1] - ys[i])
end

"""
    mccoy_308_168_cmpa(mach, alpha_t) -> C_Mpα

Magnus moment coefficient of the reference bullet, interpolated in Mach and in
total angle of attack. `alpha_t` is in radians; the table's yaw axis is in square
degrees, so it is converted here.

Interpolation is done along yaw first, inside each of the two bracketing Mach
rows, then between them — the yaw breakpoints are not the same from one Mach to
the next, so a plain rectangular bilinear scheme would not apply.
"""
function _interp_yaw_table(rows, mach::Real, alpha_t::Real)
    a2 = rad2deg(alpha_t)^2
    mach <= rows[1][1] && return _interp1(rows[1][2], rows[1][3], a2)
    mach >= rows[end][1] && return _interp1(rows[end][2], rows[end][3], a2)
    i = findlast(r -> r[1] <= mach, rows)
    i = clamp(i, 1, length(rows) - 1)
    lo = _interp1(rows[i][2],   rows[i][3],   a2)
    hi = _interp1(rows[i+1][2], rows[i+1][3], a2)
    f  = (mach - rows[i][1]) / (rows[i+1][1] - rows[i][1])
    return lo + f * (hi - lo)
end

mccoy_308_168_cmpa(mach::Real, alpha_t::Real) =
    _interp_yaw_table(MCCOY_308_168_CMPA_TABLE, mach, alpha_t)

# Yaw dependence of the overturning moment: C_Mα = C_Mα0 + C_Mα2 sin²α_t.
const MCCOY_308_168_CMA2 = [
    0.0 -4.3; 0.95 -4.3; 1.0 -4.35; 1.05 -4.4; 2.5 -4.4]

"""
    mccoy_308_168_aero() -> AeroCoefficients6DOF

Measured coefficient set for the .308", 168 gr Sierra International bullet.

The normal-force slope the solver wants is recovered from the tabulated lift-curve
slope as `C_Nα = C_Lα + C_D`, both read at the same Mach number.

!!! danger "Ne pas « améliorer » ces tables en spline"
    Toutes les lectures ci-dessous passent par `interp_table`, c'est-à-dire une
    interpolation **linéaire**, alors que `drag_coefficient` et le G7 standard sont
    passés à la spline le 2026-08-14. Ce n'est pas un oubli.

    Un test de décimation — un nœud sur deux retiré puis reconstruit — donne des
    verdicts **opposés** selon la table :

    | table | nœuds | linéaire | spline |
    |---|---|---|---|
    | G7 standard | 65 | 0,099 % | **0,049 %** |
    | G1 standard | 53 | 0,424 % | **0,069 %** |
    | McCoy .308 C_D0 | 15 | **0,68 %** | 8,89 % |
    | McCoy 105 C_D0 | 12 | **2,86 %** | 6,93 % |

    Le critère n'est pas « spline contre linéaire » mais **la densité
    d'échantillonnage rapportée à la courbure**. Les tables standard sont denses et
    régulières : la spline y suit la courbe. Les tables de McCoy sont creuses et
    **volontairement resserrées sur le genou transsonique** — ses points sont à
    0,90 / 0,95 / 1,00 / 1,05 / 1,10 / 1,20 —, si bien que la montée brutale impose
    des dérivées que la spline propage en oscillant dans les intervalles voisins.
    Elle y est jusqu'à **treize fois pire** que la ligne droite.

    Les validations publiées ont par ailleurs été calées avec ces lectures
    linéaires : les changer obligerait à tout revalider pour un résultat moins bon.
"""
function mccoy_308_168_aero()
    cd0(Ma) = DragModels.interp_table(MCCOY_308_168_CD0, Ma)
    return AeroCoefficients6DOF(
        Cd0_func = cd0,
        Cda2     = Ma -> DragModels.interp_table(MCCOY_308_168_CDD2, Ma),
        CNa      = Ma -> DragModels.interp_table(MCCOY_308_168_CLA, Ma) + cd0(Ma),
        CMa      = (Ma, al) -> DragModels.interp_table(MCCOY_308_168_CMA0, Ma) +
                               DragModels.interp_table(MCCOY_308_168_CMA2, Ma) * sin(al)^2,
        Clp      = Ma -> DragModels.interp_table(MCCOY_308_168_CLP, Ma),
        CMq_CMad = Ma -> DragModels.interp_table(MCCOY_308_168_CMQ, Ma),
        CMpa     = mccoy_308_168_cmpa,
        provenance = :measured,
    )
end

#  ────────────────────────────────────────────────────────────────────────────
#  105 mm HE M1 — McCoy, Appendix B of chapter 9 (BRL aeroballistic nomenclature)
#
#  The second published 6-DOF case, and a deliberately unlike one: a 14.97 kg
#  artillery shell instead of an 11 g bullet, fired subsonic (Charge 1, Mach
#  0.602) as well as supersonic (Charge 7, Mach 1.449), at 45° and 70° quadrant
#  elevation. Physical properties are Table 9.3, the eight runs Table 9.4.
#
#  Two differences from the bullet set are worth noting before use: here the
#  LIFT varies with yaw too (C_Lα = C_Lα0 + C_Lα2 sin²α), and the Magnus moment
#  is mostly POSITIVE where the bullet's is negative at small yaw.
#  ────────────────────────────────────────────────────────────────────────────

const MCCOY_105_CD0 = [
    0.0 0.124; 0.875 0.124; 0.925 0.150; 0.965 0.200; 0.990 0.350;
    1.025 0.375; 1.085 0.415; 1.19 0.415; 1.35 0.385; 1.80 0.335;
    2.0 0.318; 2.5 0.276]

const MCCOY_105_CDD2 = [
    0.0 3.2; 0.88 3.2; 0.97 6.3; 0.99 4.0; 1.15 5.0; 1.25 5.4; 1.3 5.5; 2.5 5.5]

const MCCOY_105_CLP = [
    0.0 -0.0178; 0.43 -0.0149; 0.70 -0.0135; 0.91 -0.0126; 1.4 -0.0110;
    1.75 -0.0101; 2.1 -0.0094; 2.5 -0.0087]

const MCCOY_105_CLA0 = [
    0.0 1.63; 0.4 1.63; 0.7 1.41; 0.89 1.22; 0.99 1.73; 1.09 1.57;
    1.5 1.97; 2.0 2.25; 2.5 2.50]

const MCCOY_105_CLA2 = [
    0.0 0.1; 0.2 0.1; 0.6 3.5; 0.8 6.6; 0.985 9.2; 1.09 8.8; 1.3 12.0;
    1.5 13.7; 2.0 16.0; 2.5 17.0]

const MCCOY_105_CMA0 = [
    0.0 3.55; 0.46 3.55; 0.61 3.76; 0.78 3.92; 0.87 3.96; 0.925 4.85;
    0.97 4.0; 1.09 3.83; 1.5 3.75; 2.5 3.75]

const MCCOY_105_CMA2 = [
    0.0 -2.9; 0.4 -2.9; 0.45 -3.1; 0.65 -4.4; 0.78 -3.45; 0.885 -1.78;
    0.98 -3.0; 1.075 -2.1; 1.25 -3.325; 1.5 -4.45; 2.0 -4.6; 2.5 -4.6]

const MCCOY_105_CMQ = [0.0 -3.15; 0.79 -3.15; 1.15 -9.1; 1.55 -9.5]

# Yaw-dependent tables: (Mach, α² breakpoints in deg², coefficient values).
# Repeated Mach rows are McCoy's way of holding a curve flat between breakpoints.
const MCCOY_105_CNPA_TABLE = [
    (0.000, [0.0, 632.0, 908.0, 1316.0],          [-0.34, -0.91, -1.42, -2.63]),
    (0.220, [0.0, 632.0, 908.0, 1316.0],          [-0.34, -0.91, -1.42, -2.63]),
    (0.310, [0.0, 21.4, 364.5, 638.0, 1316.0],    [-0.125, -0.465, -0.503, -1.015, -2.92]),
    (0.480, [0.0, 348.5, 1316.0],                 [-0.34, -0.591, -2.45]),
    (0.999, [0.0, 348.5, 1316.0],                 [-0.34, -0.591, -2.45]),
    (1.001, [0.0, 706.0],                         [-0.36, -1.68]),
    (1.550, [0.0, 706.0],                         [-0.36, -1.68]),
]

const MCCOY_105_CMPA_TABLE = [
    (0.000, [0.0, 403.6, 630.2, 1316.0],          [0.10, 0.173, 0.345, 2.35]),
    (0.220, [0.0, 403.6, 630.2, 1316.0],          [0.10, 0.173, 0.345, 2.35]),
    (0.310, [0.0, 410.8, 637.7, 915.9, 1316.0],   [0.10, 0.133, 0.471, 1.276, 2.35]),
    (0.480, [0.0, 27.5, 375.2, 1316.0],           [-0.46, 0.08, 0.022, 0.94]),
    (0.810, [0.0, 27.5, 375.2, 1316.0],           [-0.46, 0.08, 0.022, 0.94]),
    (0.870, [0.0, 315.3, 743.9],                  [0.4175, 0.053, 0.285]),
    (0.920, [0.0, 315.3, 743.9],                  [0.4175, 0.053, 0.285]),
    (0.960, [0.0, 322.2, 1316.0],                 [0.3747, 0.05, 0.665]),
    (0.995, [0.0, 322.2, 1316.0],                 [0.3747, 0.05, 0.665]),
    (1.020, [0.0, 375.2],                         [0.20, 0.301]),
    (1.100, [0.0, 375.2],                         [0.20, 0.301]),
    (1.210, [0.0, 403.6, 705.7],                  [0.193, 0.50, 0.445]),
    (1.280, [0.0, 403.6, 705.7],                  [0.193, 0.50, 0.445]),
    (1.460, [0.0, 410.8],                         [0.215, 0.495]),
    (1.550, [0.0, 410.8],                         [0.215, 0.495]),
]

"""
    mccoy_105mm_m1_aero() -> AeroCoefficients6DOF

Measured coefficient set for the 105 mm HE M1 shell, Appendix B of chapter 9.
"""
function mccoy_105mm_m1_aero()
    cd0(Ma) = DragModels.interp_table(MCCOY_105_CD0, Ma)
    return AeroCoefficients6DOF(
        Cd0_func = cd0,
        Cda2     = Ma -> DragModels.interp_table(MCCOY_105_CDD2, Ma),
        # normal force from the tabulated lift-curve slope, which here is itself
        # yaw-dependent: C_Nα = C_Lα0 + C_Lα2 sin²α + C_D
        CNa      = (Ma, al) -> DragModels.interp_table(MCCOY_105_CLA0, Ma) +
                               DragModels.interp_table(MCCOY_105_CLA2, Ma) * sin(al)^2 +
                               cd0(Ma),
        CMa      = (Ma, al) -> DragModels.interp_table(MCCOY_105_CMA0, Ma) +
                               DragModels.interp_table(MCCOY_105_CMA2, Ma) * sin(al)^2,
        Clp      = Ma -> DragModels.interp_table(MCCOY_105_CLP, Ma),
        CMq_CMad = Ma -> DragModels.interp_table(MCCOY_105_CMQ, Ma),
        CNpa     = (Ma, al) -> _interp_yaw_table(MCCOY_105_CNPA_TABLE, Ma, al),
        CMpa     = (Ma, al) -> _interp_yaw_table(MCCOY_105_CMPA_TABLE, Ma, al),
        provenance = :measured,
    )
end

"""
    mccoy_105mm_m1_shot(; twist_cal=18.0, charge=7, qe_deg=70.0, kwargs...)

One of the eight runs of Table 9.4. `twist_cal` is 18 or 25 calibers per turn,
`charge` 1 (205 m/s) or 7 (493 m/s), `qe_deg` 45 or 70.

The muzzle yaw *rate* is the published one for that run — McCoy chose each to
produce a first maximum yaw of 3.0°, the value observed in field firings of the
M103 howitzer.

Published results to check a run against: muzzle gyroscopic stability factor
**3.1** at 1/18 and **1.6** at 1/25, at either charge; first maximum yaw **3.0°**
in all eight runs; and, for Charge 7 at 70° QE, a level-ground range of **7300 m**,
a time of flight of about **70.5 s** and a summit **slightly over 6000 m**.
"""
function mccoy_105mm_m1_shot(; twist_cal::Real = 18.0, charge::Int = 7,
                               qe_deg::Real = 70.0, kwargs...)
    d_in = 104.8 / 25.4
    v_ms = charge == 1 ? 205.0 : 493.0
    yaw_rate = Dict((18, 1, 45) => 1.44, (18, 1, 70) => 1.47,
                    (18, 7, 45) => 3.61, (18, 7, 70) => 3.64,
                    (25, 1, 45) => 0.76, (25, 1, 70) => 0.79,
                    (25, 7, 45) => 1.97, (25, 7, 70) => 1.98
                   )[(round(Int, twist_cal), charge, round(Int, qe_deg))]
    return ShotParameters6DOF(;
        mass_grains        = 14.97 / 0.0647989e-3,   # 14.97 kg
        caliber_in         = d_in,
        bullet_length_in   = 49.47 / 2.54,           # 49.47 cm
        muzzle_vel_fps     = v_ms / 0.3048,
        twist_in           = twist_cal * d_in,
        initial_yaw_deg    = 0.0,
        initial_pitch_rate = yaw_rate,
        aero               = mccoy_105mm_m1_aero(),
        Ix                 = 0.02326,                # Table 9.3
        Iy                 = 0.23118,
        launch_angle_deg   = Float64(qe_deg),
        stop_on_impact     = true,
        max_time_s         = 120.0,
        # The fast body-frame mode here turns at |(Ix−Iy)/Iy|·p ≈ 1500 rad/s, so
        # 1e-4 s still gives ~40 steps per cycle — the bullet's 5e-6 would cost
        # 14 million steps for a 70 s flight and buy nothing.
        dt                 = 1.0e-4,
        record_every       = 100,
        kwargs...)
end

#  ────────────────────────────────────────────────────────────────────────────
#  5.56 mm NATO — McCoy, BRL-MR-3476 (October 1985)
#
#  "Aerodynamic and Flight Dynamic Characteristics of the New Family of 5.56mm
#  NATO Ammunition", spark-range firings at Aberdeen. US Government work,
#  approved for public release, distribution unlimited.
#
#  The first MEASURED small-arms Magnus and pitch-damping data this library has.
#  The .308 set above is McCoy's textbook case and the 105 mm is artillery; this
#  is a production service bullet, and the one regime the site's readers shoot.
#
#  Tables 3 and 4 are transcribed BELOW AS FIRED — one row per round, exactly as
#  printed, with NaN where the reduction produced no value. Every curve the
#  library uses is derived from these rows at load time, so there is no second
#  hand-copied copy to drift out of step.
#  ────────────────────────────────────────────────────────────────────────────

"""
Spark-range rounds from Table 3 (ball) and Table 4 (tracer) of BRL-MR-3476.

Columns: `round · Mach · α_t (deg) · C_D · C_Mα · C_Lα · C_Mpα · (C_Mq+C_Mα̇) · CP_N`.
`NaN` marks a dashed cell — the reduction of that round did not yield that
coefficient. Round numbers ending in `.1` / `.2` are McCoy's `(a)` / `(b)` split
reductions of one round; the three `SS-109` rounds numbered 144xx were previously
fired from an Obermeyer 7" twist barrel, which is why they sit apart in Mach.

⚠️ **`C_D` is measured AT the listed `α_t`, not at zero yaw.** Only the low-yaw
rounds are used for the drag curve — see [`brl3476_aero`](@ref).

⚠️ One cell is deliberately `NaN` that is not dashed in the report: `C_Mα` of round
14415, printed over a smudge. Its neighbour 14416 reads 2.54 at the same Mach, so
the value is almost certainly 2.55 — which is exactly why it is left out rather
than guessed.
"""
const BRL3476_SS109 = [
#   round     Mach   α_t     C_D    C_Mα   C_Lα   C_Mpα  C_Mq+C_Mα̇  CP_N
    16221.0   2.638  1.30   0.2898  2.52   NaN   -0.04   -3.77      NaN
    16220.0   2.622  1.28   0.2857  2.55   2.82  -0.03   -6.42      2.34
    16231.0   1.964  0.88   0.3346  2.76   NaN    NaN     NaN       NaN
    16230.0   1.860  3.78   0.3730  2.77   2.67   0.04   -5.02      2.43
    16238.0   1.179  2.90   0.4614  2.97   2.15  -0.31   -3.15      2.65
    16239.0   1.138  1.47   0.4488  NaN    NaN    NaN     NaN       NaN
    16245.1   0.757  2.57   0.1753  NaN    NaN    NaN     NaN       NaN
    16245.2   0.746  4.06   0.1832  NaN    1.81   NaN     NaN       NaN
    16246.0   0.736  4.82   0.2204  3.05   1.61  -0.16   -0.22      3.19
    14414.0   2.645  0.40   0.2834  NaN    NaN    NaN     NaN       NaN
    14416.0   2.629  1.70   0.2919  2.54   NaN    0.16   -7.05      NaN
    14415.0   2.625  1.84   0.2956  NaN    NaN    0.13   -7.38      NaN ]

"@ref BRL3476_SS109 — Table 3, M855 ball."
const BRL3476_M855 = [
    16222.0   2.730  1.90   0.3039  2.40   2.61   0.10   -5.52      2.36
    16223.0   2.714  1.22   0.3007  2.35   NaN   -0.03   -6.51      NaN
    16228.0   1.875  1.62   0.3741  2.66   NaN    NaN     NaN       NaN
    16229.0   1.869  3.65   0.3851  2.64   2.59   0.06   -4.99      2.42
    16237.0   1.137  2.27   0.4874  2.82   2.16   NaN     NaN       2.60
    16238.0   1.072  2.87   0.4826  2.85   2.14  -0.59    1.27      2.62
    16243.0   0.674  6.04   0.2918  2.67   1.85  -0.47    0.00      2.79 ]

"@ref BRL3476_SS109 — Table 4, L110 tracer."
const BRL3476_L110 = [
    16224.0   2.543  0.78   0.2983  2.44   NaN    0.39  -10.4       NaN
    16225.0   2.533  2.76   0.3083  2.26   2.96   0.44  -10.1       3.21
    16234.0   1.890  6.14   0.4267  2.50   2.89   0.89  -18.6       3.27
    16235.0   1.822  0.52   0.3737  NaN    NaN    NaN     NaN       NaN
    16242.0   1.121  5.21   0.5398  NaN    NaN    NaN     NaN       NaN
    16249.1   0.854  4.31   0.2470  NaN    NaN    NaN     NaN       NaN
    16249.2   0.841  5.29   0.2074  NaN    NaN    NaN     NaN       NaN ]

"@ref BRL3476_SS109 — Table 4, M856 tracer."
const BRL3476_M856 = [
    16226.0   2.575  1.64   0.3044  2.28   2.73   0.45  -14.6       3.32
    16227.0   2.566  1.70   0.3140  2.29   2.81   0.56  -12.5       3.30
    16232.0   1.995  2.77   0.3308  2.76   2.87   0.59  -14.9       3.43
    16233.0   1.870  4.01   0.4025  2.70   NaN    0.95  -19.8       NaN
    16240.1   1.184  4.03   0.5031  2.88   2.10   NaN     NaN       3.67
    16240.2   1.140  5.95   0.5334  2.86   2.11   NaN     NaN       3.65
    16241.0   1.117  4.67   0.4850  2.87   2.05  -0.30   -5.15      3.70
    16247.0   0.758  7.59   0.3174  2.93   1.92   NaN     NaN       3.87
    16248.0   0.738 10.04   0.2847  NaN    NaN    NaN     NaN       NaN ]

"""
Average physical characteristics, Table 1 of BRL-MR-3476, plus the dimensioned
sketches of Figures 2 and 3. One caliber is 5.69 mm throughout.

Fields: `mass_g · cg_cal · Ix · Iy` (both g·cm², as printed) then the contour
`len_cal · nose_cal · ogive_R_cal · bt_cal · bt_deg`.

⚠️ The `2.00` / `1.90` cote on Figure 2 is **not** the ogive length: read that way
the tangency condition returns a meplat of 0.52 caliber, which is absurd. The
ogive is the `2.76` / `2.71` cote, and the meplat then comes out at 0.067 / 0.041
caliber — see [`brl3476_geometry`](@ref).
"""
const BRL3476_PHYSICAL = Dict(
    :ss109 => (mass_g = 4.03, cg_cal = 1.52, Ix = 0.1425, Iy = 1.112,
               len_cal = 4.07, nose_cal = 2.76, ogive_R_cal = 8.4, bt_cal = 0.45, bt_deg = 9.75),
    :m855  => (mass_g = 4.05, cg_cal = 1.54, Ix = 0.1426, Iy = 1.150,
               len_cal = 4.05, nose_cal = 2.71, ogive_R_cal = 7.9, bt_cal = 0.40, bt_deg = 8.50),
    :l110  => (mass_g = 4.09, cg_cal = 2.52, Ix = 0.1573, Iy = 1.874),
    :m856  => (mass_g = 4.19, cg_cal = 2.57, Ix = 0.1634, Iy = 1.987))

const BRL3476_CALIBER_M = 5.69e-3

"""
    brl3476_curve(rows, col; max_yaw, tol) -> [mach value]

Collapse one coefficient column of a [`BRL3476_SS109`](@ref)-shaped table into a
strictly increasing Mach curve fit for `interp_table`.

Rounds whose cell is `NaN` are dropped; rounds fired above `max_yaw` degrees are
dropped too (used to keep the drag curve near zero yaw). Rounds closer together
than `tol` in Mach are averaged — the test series fired several rounds within a
few thousandths of each other, and `interp_table` needs a strictly increasing
abscissa.
"""
function brl3476_curve(rows::AbstractMatrix, col::Integer;
                       max_yaw::Real = Inf, tol::Real = 0.05)
    pts = sort([(rows[i, 2], rows[i, col]) for i in axes(rows, 1)
                if !isnan(rows[i, col]) && rows[i, 3] <= max_yaw], by = first)
    isempty(pts) && error("BRL-MR-3476 : la colonne $col n'a aucune mesure exploitable")
    groups = Vector{Vector{Tuple{Float64,Float64}}}()
    for p in pts
        if !isempty(groups) && p[1] - groups[end][end][1] < tol
            push!(groups[end], p)
        else
            push!(groups, [p])
        end
    end
    return hcat([sum(first, g) / length(g) for g in groups],
                [sum(last,  g) / length(g) for g in groups])
end

"""
    _brl_interp(tab, name, which) -> (Mach -> value)

Interpolate a BRL-MR-3476 curve, and **say so when asked outside the Mach band
that was actually fired**.

`interp_table` clamps at both ends, which for a dense table is harmless and for
these is not: the L110 tracer has three `C_Mpα` rounds, all above Mach 1.89, so a
bare lookup answers a confident +0.89 at Mach 1.0 — a value nobody measured. Nine
rounds cannot be made to cover a flight; the honest response is to hand back the
clamped value **and warn once**, so an out-of-band run is visible in the log
instead of silently plausible. Same reasoning as the `provenance` guard on
[`AeroCoefficients6DOF`](@ref).
"""
function _brl_interp(tab::AbstractMatrix, name::AbstractString, which::Symbol)
    lo, hi = tab[1, 1], tab[end, 1]
    # 2 % of slack at each end. Several rounds sit within a few thousandths of
    # each other and get averaged into one abscissa, so asking at the Mach of an
    # individual round can fall a hair outside the merged band. Warning on that
    # would train the reader to ignore the warning, which is worse than silence.
    lo_w, hi_w = lo * 0.98, hi * 1.02
    return function (Ma)
        if Ma < lo_w || Ma > hi_w
            @warn "BRL-MR-3476 : $name du $which extrapolé hors de la bande de tir" *
                  " (Mach $(round(Ma, digits=3)) hors de [$lo, $hi]) — valeur bloquée au bord" maxlog = 1
        end
        return DragModels.interp_table(tab, Ma)
    end
end

"""
    brl3476_geometry(which) -> BulletGeometry

Contour of a 5.56 mm NATO ball projectile from Figure 2, with the meplat solved
from the tangent-ogive condition `r_m = √(R² − Lₙ²) − (R − ½)` in calibers.

`which` is `:ss109` or `:m855`.

!!! warning "The core is not what `bullet_inertia` models"
    Both projectiles carry a **steel penetrator** ahead of a lead rear core, so
    their mass is laid out in two materials. `BulletGeometry` knows a single lead
    core inside a gilding-metal jacket. Fed this contour and the true mass, the
    estimator lands 4–5 % low on Ix, 6–9 % low on Iy and 5–7 % low on the CG
    against Table 1. That is the composite core showing, not an integration
    error, and it is the honest limit of the estimator on military ball.
"""
function brl3476_geometry(which::Symbol)
    p = BRL3476_PHYSICAL[which]
    meplat = 2 * (sqrt(p.ogive_R_cal^2 - p.nose_cal^2) - (p.ogive_R_cal - 0.5))
    return BulletGeometry(nose_frac  = p.nose_cal / p.len_cal,
                          meplat_cal = meplat,
                          bt_len_cal = p.bt_cal,
                          bt_deg     = p.bt_deg)
end

"""
    brl3476_aero(which) -> AeroCoefficients6DOF

Measured coefficient set for one 5.56 mm NATO projectile: `:ss109`, `:m855`,
`:l110` or `:m856`. Marked `provenance = :measured`, so `solve_6dof` runs it
without the placeholder warning.

# What is measured here, and what is not

Taken from the range reductions: **C_Mα**, **C_Lα** (entered as `C_Nα = C_Lα + C_D0`),
**C_Mpα** and **(C_Mq + C_Mα̇)**. The drag curve uses only rounds fired below 2° of
yaw, where the yaw-drag term is under 1 % of C_D — comfortably inside the round-to-
round scatter, which is itself about 1.4 % (rounds 16221 and 16220 differ by that
much at the same Mach and yaw).

**Not measured by this report, left at the struct defaults**: `Clp` (roll damping —
it appears only in the list of symbols), `Cda2` (yaw drag) and `CNpa` (Magnus
force). The first two matter little here; `C_Npα` our own sensitivity study found
negligible.

!!! warning "This set cannot produce a limit cycle"
    Each `C_Mpα` value is a single round at a single angle of attack, so the curve
    below is a function of **Mach only** — the yaw dependence is not resolved. It is
    precisely the sign change of `C_Mpα` with yaw that fixes limit-cycle amplitude,
    so a run on this set will not reproduce one. Compare `mccoy_308_168_cmpa`,
    which does carry a yaw axis because McCoy published the fit rather than the
    individual rounds. The scatter here hints at the same behaviour — M855 goes
    from +0.10 at 1.90° to −0.59 at 2.87° — but three points cannot separate a
    Mach effect from a yaw effect, and pretending otherwise would be inventing.
"""
function brl3476_aero(which::Symbol; yaw_drag_coef::Real = 3.5)
    rows = which === :ss109 ? BRL3476_SS109 :
           which === :m855  ? BRL3476_M855  :
           which === :l110  ? BRL3476_L110  :
           which === :m856  ? BRL3476_M856  :
           error("BRL-MR-3476 ne couvre que :ss109, :m855, :l110 et :m856")

    # C_D is measured at the round's own yaw. Reduce every round to zero yaw with
    # C_D0 = C_D − C_Dδ²·sin²α_t rather than keeping only the near-zero-yaw ones:
    # on the M855 that filter leaves TWO points, both supersonic, and the curve
    # then clamps a supersonic value straight across the transonic peak.
    corrected = copy(rows)
    for i in axes(corrected, 1)
        corrected[i, 4] -= yaw_drag_coef * sind(corrected[i, 3])^2
    end

    cd0_tab = brl3476_curve(corrected, 4)
    cd0     = _brl_interp(cd0_tab, "C_D0", which)
    cla     = _brl_interp(brl3476_curve(rows, 6), "C_Lα", which)
    return AeroCoefficients6DOF(
        Cd0_func = cd0,
        CNa      = Ma -> cla(Ma) + cd0(Ma),
        CMa      = _brl_interp(brl3476_curve(rows, 5), "C_Mα", which),
        CMpa     = _brl_interp(brl3476_curve(rows, 7), "C_Mpα", which),
        CMq_CMad = _brl_interp(brl3476_curve(rows, 8), "C_Mq+C_Mα̇", which),
        provenance = :measured,
    )
end

"""
    brl3476_shot(which = :m855; kwargs...) -> ShotParameters6DOF

A 5.56 mm NATO shot with the **measured** moments of inertia of Table 1 and the
measured aerodynamics of [`brl3476_aero`](@ref).

⚠️ Two inputs are *not* from the report and are set to service values: the twist,
taken as the 1:7″ of the M16A2 and M249 that the report discusses in its
dispersion section, and the muzzle velocity, taken as the highest Mach of the test
series converted at ICAO sea level (Mach 2.730 → 3048 fps for the M855). Override
them for any specific rifle.
"""
function brl3476_shot(which::Symbol = :m855; kwargs...)
    p = BRL3476_PHYSICAL[which]
    d = BRL3476_CALIBER_M
    v0_fps = maximum((which === :ss109 ? BRL3476_SS109 :
                      which === :m855  ? BRL3476_M855  :
                      which === :l110  ? BRL3476_L110  : BRL3476_M856)[:, 2]) *
             340.3 / 0.3048
    return ShotParameters6DOF(;
        mass_grains        = p.mass_g / 0.06479891,
        caliber_in         = d / 0.0254,
        bullet_length_in   = get(p, :len_cal, 4.06) * d / 0.0254,
        muzzle_vel_fps     = v0_fps,
        twist_in           = 7.0,
        twist_direction    = 1,
        initial_yaw_deg    = 0.0,
        initial_pitch_rate = 25.0,
        zero_range_yd      = 600.0,
        target_range_yd    = 600.0,
        aero               = brl3476_aero(which),
        Ix                 = p.Ix * 1.0e-7,        # g·cm² -> kg·m²
        Iy                 = p.Iy * 1.0e-7,
        kwargs...)
end

"""
    mccoy_308_168_shot(; kwargs...) -> ShotParameters6DOF

Example 9.1 as published: 2600 fps, 12" right-hand twist, sea-level ICAO, zeroed
at 1000 yards, with the measured moments of inertia of Table 9.2 and a muzzle yaw
*rate* of 25 rad/s — which is what produces the published first maximum yaw of
2.0 degrees, rather than an initial yaw angle. Any keyword overrides a default.

Published results to check a run against: muzzle gyroscopic stability factor
**1.70**; pitch-yaw amplitude of about **1.75°** over 180–200 yards (Figure 9.3)
and about **2.2°** over 580–600 yards (Figure 9.4); and, over the whole 1000
yards, an amplitude that never exceeds 5 degrees.
"""
function mccoy_308_168_shot(; kwargs...)
    lb_in2 = 0.45359237 * 6.4516e-4        # 1 lb·in² -> kg·m²
    return ShotParameters6DOF(;
        mass_grains        = 168.0,
        caliber_in         = 0.308,
        bullet_length_in   = 1.226,         # Table 9.2
        muzzle_vel_fps     = 2600.0,
        twist_in           = 12.0,
        twist_direction    = 1,
        initial_yaw_deg    = 0.0,
        initial_pitch_rate = 25.0,          # rad/s, muzzle yaw rate
        zero_range_yd      = 1000.0,
        target_range_yd    = 1000.0,
        aero               = mccoy_308_168_aero(),
        Ix                 = 0.000247 * lb_in2,
        Iy                 = 0.001838 * lb_in2,
        kwargs...)
end

"""Quaternion multiplication."""
function quat_mult(a, b)
    return (
        a[1]*b[1] - a[2]*b[2] - a[3]*b[3] - a[4]*b[4],
        a[1]*b[2] + a[2]*b[1] + a[3]*b[4] - a[4]*b[3],
        a[1]*b[3] - a[2]*b[4] + a[3]*b[1] + a[4]*b[2],
        a[1]*b[4] + a[2]*b[3] - a[3]*b[2] + a[4]*b[1],
    )
end

"""Rotation matrix from unit quaternion."""
function quat_to_dcm(q0, q1, q2, q3)
    return [
        1-2(q2^2+q3^2)    2(q1*q2-q0*q3)    2(q1*q3+q0*q2);
        2(q1*q2+q0*q3)    1-2(q1^2+q3^2)    2(q2*q3-q0*q1);
        2(q1*q3-q0*q2)    2(q2*q3+q0*q1)    1-2(q1^2+q2^2)
    ]
end

"""
    solve_6dof(p::ShotParameters6DOF) -> Vector{State6DOF}

Full 6-DOF rigid-body trajectory solver, RK4 on the whole 13-component state
(position, velocity, quaternion, body rates). Quaternion kinematics are
singularity-free; the quaternion is renormalized after every step.

Set `integrator = :euler` to fall back on the first-order scheme this module used
until 2026-08-12. It is kept as a reference point, not as a working alternative —
see below.

# Choosing `dt`

The step must resolve the **fast body-frame mode**, not the epicyclic motion an
observer sees. The gyroscopic coupling in the rate equations rotates (q, r) at
|(Ix−Iy)/Iy|·p — 24 000 rad/s for the default 6.5 mm bullet, a 262 µs period —
while the epicyclic period is a tenfold longer. The default `dt = 5e-6 s` puts
~52 steps in the fast cycle, the rule of thumb McCoy reports for 6-DOF work.

Measured convergence at 100 yd, against a run at `dt = 6.25e-7 s`:

| `dt` (s) | error on drop | error on total yaw |
|---------:|--------------:|-------------------:|
| 4e-5     | 1.9e-3 m      | 8.5e-3 rad         |
| 1e-5     | 7.4e-6 m      | 5.5e-5 rad         |
| **5e-6** | **2.0e-6 m**  | **2.9e-6 rad**     |

The yaw column falls by a factor of ~16 per halving, i.e. the fourth order is
actually attained; the drop column reaches its round-off floor near 1e-6 m.

!!! danger "Why `:euler` is a fallback and not an option"
    Explicit Euler is **unconditionally unstable** on the undamped gyroscopic
    oscillator these equations contain: its amplification factor is
    √(1+(ωh)²) > 1 for every step size. With the historical default of `dt = 1e-4 s`
    the body rates reached 1e159 rad/s and the quaternion went `NaN` within 39 ms of
    flight — the solver never produced a usable 1000 yd trajectory. Shrinking the
    step only slows the blow-up: at `dt = 1e-6 s` the pitch rate still reaches
    1e18 rad/s. RK4 is stable here because its region of absolute stability covers
    the imaginary axis up to |ωh| ≈ 2.83, and the default step sits at ωh ≈ 0.12.

!!! note "Placeholder aerodynamic coefficients"
    With the library's default [`AeroCoefficients6DOF`](@ref) the computed total
    angle of attack **grows** along the trajectory (0.46° at the muzzle, 1.26° at
    270 m, 18.7° at 1000 yd). Those coefficients are placeholders, not measurements
    for any particular bullet. Supply a measured set before reading anything
    physical into a run.

!!! note "Validated against the published reference case"
    Checked on 2026-08-13 against [`mccoy_308_168_shot`](@ref):

    | | published | here |
    |---|---:|---:|
    | muzzle Sg | 1.70 | **1.701** |
    | first maximum yaw | 2.0° | **1.97°** |
    | 180–200 yd (Fig. 9.3) | ~1.75° | **1.53°** |
    | 580–600 yd (Fig. 9.4) | ~2.2° | **2.09°** |
    | over 1000 yd | never > 5° | **5.19°** |
    | drop at the 1000 yd zero | ≈ 0 | **−0.24 m** |

    Getting here took **three** corrections, and no two of them were enough:

    1. the normal force was applied along `+v_perp` — that is, along the bullet's
       own crossflow rather than the air's, so the lift pointed the wrong way;
    2. the Magnus moment carried the opposite sign, which reversed the roles of
       the two branches of `C_Mpα` and turned the limit cycle into a runaway;
    3. the rate-normalized coefficients were divided by `2V` (NACA convention)
       instead of the `V` of the BRL system the source's tables are quoted in.

    All three spare the epicyclic frequency, which is set by `C_Mα` and the spin,
    and land squarely on the damping exponents — which is why the beat validated
    at 2 % throughout while the amplitude ran to 179°. An earlier pass tried each
    sign **on its own** and drew the wrong conclusion from each: flipping only the
    normal force over-damps to extinction (0.02° by 600 yd), which reads as a
    refutation and is not one.

    Cross-checks that now agree: the growth rate in the linear regime matches
    McCoy's λ_S, and the amplitude settles inside the 1.90° limit cycle that eq.
    (13.62)–(13.65) predict from the same coefficients. §10.11 states the mechanism
    outright — the slow arm is undamped at small yaw, and `C_Mpα` becoming less
    negative at larger yaw is what bounds it.

!!! note "Second published case, nothing tuned to it"
    Example 9.2, the 105 mm HE M1 shell ([`mccoy_105mm_m1_shot`](@ref)), tests the
    same three corrections on a deliberately unlike projectile: 14.97 kg against
    11 g, 45° and 70° quadrant elevation against flat fire, subsonic as well as
    supersonic, and a Magnus moment that is mostly *positive* where the bullet's is
    negative. Against the four published trajectories of p. 201, at 1/18 twist:

    | | range | time of flight | summit |
    |---|---:|---:|---:|
    | Charge 1, 45° | −0.3 % | +0.2 % | −0.3 % |
    | Charge 1, 70° | −2.2 % | +0.2 % | −0.4 % |
    | Charge 7, 45° | −0.5 % | −0.1 % | −0.6 % |
    | Charge 7, 70° | −4.3 % | +0.0 % | −1.4 % |

    Muzzle `Sg` comes out 3.13 and 1.62 against the published 3.1 and 1.6. The
    sharpest check is Figure 9.11: at apogee the solver swings between **0.25° and
    4.02°**, where the figure's inset gives a coning circle of `K_S` = 1.9° about a
    yaw of repose of 2.1° — 0.2° to 4.0°.

    The two 70° ranges are the weakest numbers here, and the cause is identified
    rather than hidden: it is not base drag — no single factor on `C_D0` aligns all
    four cases — but the summital yaw surge, which reaches 21–26° there against
    1.5–4° at 45°, and where yaw drag runs three to five times `C_D0`. Range in that
    regime is governed by the accuracy of a 20° yaw, not by the integration. Note
    also that the solver still carries no Coriolis term; over these 70 s flights
    the Eötvös contribution is only about 13 m in 7300."""
function solve_6dof(sp::ShotParameters6DOF)
    # Convert to SI
    m   = grains_to_kg(sp.mass_grains)
    d   = inches_to_m(sp.caliber_in)
    L   = inches_to_m(sp.bullet_length_in)
    v0  = fps_to_ms(sp.muzzle_vel_fps)
    tw  = inches_to_m(sp.twist_in)
    T_K = fahrenheit_to_kelvin(sp.temp_F)
    P_Pa = inhg_to_pa(sp.pressure_inhg)
    rho = air_density(P=P_Pa, T=T_K, H=sp.humidity_pct)
    a_s = speed_of_sound(T_K)

    # Air properties at height `h` above the muzzle. The given temperature and
    # pressure describe the launch point; the ICAO lapse carries them upward from
    # there, so this reduces exactly to the constants above at h = 0. Holding them
    # constant is harmless over the few metres a flat-fire trajectory rises, and
    # badly wrong for artillery: a 70° shot tops out near 6 km, where the density
    # is half its sea-level value.
    @inline function air_at(h::Float64)
        h <= 0.0 && return (rho, a_s)
        T = T_K - Atmosphere.L * h
        T <= 1.0 && return (rho, a_s)
        P = P_Pa * (T / T_K)^(Atmosphere.g0 / (Atmosphere.L * Atmosphere.R_air))
        return (air_density(P=P, T=T, H=sp.humidity_pct), speed_of_sound(T))
    end
    A   = π * d^2 / 4.0
    tr  = yards_to_m(sp.target_range_yd)

    if (sp.Ix === nothing) != (sp.Iy === nothing)
        throw(ArgumentError("supply both Ix and Iy, or neither"))
    end
    Ix, Iy = sp.Ix === nothing ?
             moments_of_inertia(m, d, L; geometry = sp.geometry) : (sp.Ix, sp.Iy)
    aero = sp.aero
    if aero.provenance !== :measured
        @warn """
        6-DOF run on PLACEHOLDER aerodynamic coefficients — they describe no real
        projectile. The solver is validated against two published cases, so its
        output will look entirely plausible and mean nothing. Supply a measured
        set (see `mccoy_308_168_aero`, `mccoy_105mm_m1_aero`) and mark it
        `provenance = :measured`.""" maxlog=1
    end

    # Resolve every coefficient to a (Mach, alpha_t) function, once.
    f_Cda2 = _as_coef(aero.Cda2);  f_CNa  = _as_coef(aero.CNa)
    f_CMa  = _as_coef(aero.CMa);   f_Clp  = _as_coef(aero.Clp)
    f_CMqa = _as_coef(aero.CMq_CMad)
    f_CNpa = _as_coef(aero.CNpa);  f_CMpa = _as_coef(aero.CMpa)
    dt = sp.dt

    # Wind
    w_ang = deg2rad(sp.wind_angle_deg)
    wx = -sp.wind_speed_mph * 0.44704 * cos(w_ang)
    wz =  sp.wind_speed_mph * 0.44704 * sin(w_ang)

    # Initial conditions
    alpha0 = deg2rad(sp.initial_yaw_deg)
    p0 = sp.twist_direction * initial_spin_rate(v0, tw)

    # Find zero angle (simplified 2D iteration). Still a first-order step, and
    # deliberately so: this is a point-mass translation with no fast rotational
    # mode, where Euler at `dt` is accurate and the loop only has to land a
    # launch angle to 1e-6 m of drop at the zero range.
    # A quadrant elevation, when given, replaces the zeroing entirely: artillery
    # is laid on an angle, not sighted in on a target.
    zr = yards_to_m(sp.zero_range_yd)
    theta0 = sp.launch_angle_deg === nothing ?
             atan(Atmosphere.g0 * zr / (2 * v0^2)) : deg2rad(sp.launch_angle_deg)
    for _ in 1:(sp.launch_angle_deg === nothing ? 30 : 0)
        xt, yt = 0.0, 0.0
        vxt = v0 * cos(theta0); vyt = v0 * sin(theta0)
        while xt < zr
            vt = sqrt(vxt^2 + vyt^2)
            Ma = vt / a_s
            cd = aero.Cd0_func(Ma)
            D = rho * cd * A / (2m)
            vxt += (-D * vt * vxt) * dt
            vyt += (-D * vt * vyt - Atmosphere.g0) * dt
            xt += vxt * dt; yt += vyt * dt
        end
        abs(yt) < 1e-6 && break
        theta0 -= yt / zr * cos(theta0)
    end

    # ── State vector, 13 components ──────────────────────────────────────
    #   1-3    position     x, y, z                      [m]
    #   4-6    velocity     u, v, w                      [m/s]
    #   7-10   quaternion   q0, q1, q2, q3
    #   11-13  body rates   p (spin), q (pitch), r (yaw) [rad/s]

    # Derivative of the state, as a pure function of the state alone. Writing it
    # this way is what lets an integrator sample the slope more than once inside
    # a step; everything else it needs is captured from the enclosing scope.
    function derivs6(s::NTuple{13,Float64})
        y_pos = s[2]                        # height above the muzzle, for `air_at`
        u, v_y, w_z = s[4], s[5], s[6]
        q0_q, q1_q, q2_q, q3_q = s[7], s[8], s[9], s[10]
        p_s, q_p, r_y = s[11], s[12], s[13]

        R = quat_to_dcm(q0_q, q1_q, q2_q, q3_q)
        b1 = R[:, 1]                        # body x-axis in the inertial frame

        vrel = [u - wx, v_y, w_z - wz]      # velocity relative to the air
        vrel_mag = norm(vrel)
        cos_alpha = clamp(dot(vrel, b1) / (vrel_mag + 1e-20), -1.0, 1.0)
        alpha_t = acos(cos_alpha)           # total angle of attack

        rho_h, a_h = air_at(y_pos)          # thins out with height; see `air_at`
        qbar = 0.5 * rho_h * vrel_mag^2
        Ma = vrel_mag / a_h

        # Coefficients at this Mach number and angle of attack. Constants and
        # Mach-only functions were wrapped once, before the loop.
        c_Cda2 = f_Cda2(Ma, alpha_t);  c_CNa  = f_CNa(Ma, alpha_t)
        c_CMa  = f_CMa(Ma, alpha_t);   c_Clp  = f_Clp(Ma, alpha_t)
        c_CMqa = f_CMqa(Ma, alpha_t)
        c_CNpa = f_CNpa(Ma, alpha_t);  c_CMpa = f_CMpa(Ma, alpha_t)

        # --- Forces (inertial frame) ---
        # Yaw drag follows the source's form, C_D = C_D0 + C_Dδ² sin²α_t.
        Cd_total = aero.Cd0_func(Ma) + c_Cda2 * sin(alpha_t)^2
        F_drag = -qbar * A * Cd_total * vrel / (vrel_mag + 1e-20)

        # `vrel` is the bullet's velocity THROUGH the air, so the crossflow the
        # body actually meets is its negative — hence the leading minus. With the
        # nose pitched up, `v_perp` points down while the normal force must point
        # up. Taking `+v_perp` inverts the lift, and that inversion hides well: it
        # leaves the epicyclic frequency untouched (P and M are set by C_Mα and
        # the spin, not by the forces) and only flips the sign of the C_Lα term in
        # McCoy's damping exponents — the yaw then grows where it should decay.
        # `n_hat` also orients the Magnus force below, so both follow from this
        # one direction.
        v_perp = vrel - dot(vrel, b1) * b1
        v_perp_mag = norm(v_perp)
        if v_perp_mag > 1e-10
            n_hat = -v_perp / v_perp_mag
            F_lift = qbar * A * c_CNa * alpha_t * n_hat
        else
            n_hat = zeros(3)
            F_lift = zeros(3)
        end

        F_magnus = zeros(3)
        if v_perp_mag > 1e-10 && abs(p_s) > 1e-3
            magnus_dir = cross(b1, n_hat)
            magnus_dir_mag = norm(magnus_dir)
            if magnus_dir_mag > 1e-10
                F_magnus = qbar * A * c_CNpa * (p_s * d / vrel_mag) *
                           alpha_t * (magnus_dir / magnus_dir_mag)
            end
        end

        F_grav = [0.0, -m * Atmosphere.g0, 0.0]
        F_total = F_drag + F_lift + F_magnus + F_grav

        # --- Moments (body frame) ---
        # BRL aeroballistic normalization: angular rates are made dimensionless
        # with (pd/V) and (q_t d/V), NOT the (pd/2V) of the NACA/aircraft
        # convention. Every coefficient a spark range publishes — C_lp, C_Mpα,
        # C_Mq+C_Mα̇, C_Npα — is quoted against that scale, so dividing by 2V
        # silently halves all four. Like the lift sign above, the error spares the
        # epicyclic frequency and lands entirely on the damping exponents.
        pdv = p_s * d / (vrel_mag + 1e-20)
        v_body = R' * vrel
        alpha_p = -v_body[3] / (v_body[1] + 1e-20)   # pitch plane
        beta_y  =  v_body[2] / (v_body[1] + 1e-20)   # yaw plane

        # `alpha_p` and `beta_y` above are the angles of the VELOCITY relative to
        # the body; the body's own displacement from the flow is their negative.
        # A spin-stabilized bullet is statically UNSTABLE (C_Mα > 0 overturns), so
        # the moment has to amplify that displacement — hence the leading minus.
        # Without it the solver treated the projectile as statically stable: with
        # spin removed the yaw oscillated instead of diverging, and with spin the
        # epicyclic beat came out at √(P²+4M) instead of √(P²−4M), 1.96× too fast.
        #
        # The Magnus moment carries the opposite sign to the one first written
        # here. Its own sign decides the limit cycle: C_Mpα is negative at small
        # yaw and positive past about 2.4° (Appendix A), so it must FEED the slow
        # arm near zero and DAMP it once the yaw opens — that reversal is what
        # caps the motion. Written the other way round the two roles swap and the
        # yaw runs away, which is exactly what the solver used to do.
        Mx = qbar * A * d^2 * c_Clp * pdv
        Mq_body = qbar * A * d * (-c_CMa * alpha_p +
                  c_CMpa * pdv * beta_y +
                  c_CMqa * q_p * d / (vrel_mag + 1e-20))
        Mr_body = qbar * A * d * (-c_CMa * beta_y -
                  c_CMpa * pdv * alpha_p +
                  c_CMqa * r_y * d / (vrel_mag + 1e-20))

        return (
            u, v_y, w_z,                                    # ẋ = v
            F_total[1] / m, F_total[2] / m, F_total[3] / m, # v̇ = F/m
            0.5 * (-p_s*q1_q - q_p*q2_q - r_y*q3_q),        # quaternion kinematics
            0.5 * ( p_s*q0_q + r_y*q2_q - q_p*q3_q),
            0.5 * ( q_p*q0_q - r_y*q1_q + p_s*q3_q),
            0.5 * ( r_y*q0_q + q_p*q1_q - p_s*q2_q),
            Mx / Ix,                                        # Euler's rotational
            (Mq_body - (Ix - Iy) * p_s * r_y) / Iy,         #   equations
            (Mr_body + (Ix - Iy) * p_s * q_p) / Iy,
        )
    end

    # Quantities that are reported but never integrated.
    function diagnostics6(s::NTuple{13,Float64})
        u, v_y, w_z = s[4], s[5], s[6]
        v_total = sqrt(u^2 + v_y^2 + w_z^2)
        b1 = quat_to_dcm(s[7], s[8], s[9], s[10])[:, 1]
        vrel = [u - wx, v_y, w_z - wz]
        cos_alpha = clamp(dot(vrel, b1) / (norm(vrel) + 1e-20), -1.0, 1.0)
        return (acos(cos_alpha), v_total / a_s, v_total)
    end

    # Classical fourth-order Runge-Kutta on the whole 13-component state.
    function step_rk4(s::NTuple{13,Float64})
        k1 = derivs6(s)
        k2 = derivs6(ntuple(i -> s[i] + 0.5dt * k1[i], Val(13)))
        k3 = derivs6(ntuple(i -> s[i] + 0.5dt * k2[i], Val(13)))
        k4 = derivs6(ntuple(i -> s[i] + dt * k3[i], Val(13)))
        return ntuple(i -> s[i] + (dt / 6.0) * (k1[i] + 2k2[i] + 2k3[i] + k4[i]), Val(13))
    end

    # Legacy first-order scheme, kept as a fallback and as a reference point.
    # Reproduces the historical behaviour exactly: velocities advance first and
    # positions ride the NEW velocities (semi-implicit Euler for translation),
    # while attitude and body rates advance on the old state (explicit Euler).
    # Do not expect it to hold up — on the gyroscopic terms it is unstable at
    # every step size, see the note on `solve_6dof`.
    function step_euler(s::NTuple{13,Float64})
        k = derivs6(s)
        u  = s[4] + k[4] * dt
        v_ = s[5] + k[5] * dt
        w_ = s[6] + k[6] * dt
        return (s[1] + u * dt, s[2] + v_ * dt, s[3] + w_ * dt,
                u, v_, w_,
                s[7] + k[7] * dt, s[8] + k[8] * dt, s[9] + k[9] * dt, s[10] + k[10] * dt,
                s[11] + k[11] * dt, s[12] + k[12] * dt, s[13] + k[13] * dt)
    end

    step = if sp.integrator === :rk4
        step_rk4
    elseif sp.integrator === :euler
        step_euler
    else
        throw(ArgumentError("integrator must be :rk4 or :euler, got :$(sp.integrator)"))
    end
    sp.record_every >= 1 || throw(ArgumentError("record_every must be >= 1"))

    # Initial state. The bullet leaves the bore pointing along its velocity, so the
    # attitude must carry the launch elevation `theta0` — the previous version left
    # the axis on the inertial x-axis and so injected a spurious initial angle of
    # attack equal to the launch angle (0.83° for a 1000-yard zero).
    # `alpha0` is added on top: being a rotation about z, in this frame (y up,
    # z right) it offsets the nose vertically, i.e. it is a pitch offset.
    tilt = theta0 + alpha0
    s = (0.0, 0.0, 0.0,
         v0 * cos(theta0), v0 * sin(theta0), 0.0,
         cos(tilt/2), 0.0, 0.0, sin(tilt/2),
         p0, sp.initial_pitch_rate, 0.0)

    results = State6DOF[]
    t = 0.0
    n = 0

    # Flat fire stops at the target range; a lofted shot stops when it comes back
    # down, which is the only end it has. `t > 0.1` keeps the muzzle out of it.
    stop = sp.stop_on_impact ?
           ((st, tt) -> st[2] < 0.0 && tt > 0.1) :
           ((st, tt) -> st[1] > tr)
    while !stop(s, t) && t < sp.max_time_s
        if n % sp.record_every == 0
            alpha_t, mach, v_total = diagnostics6(s)
            push!(results, State6DOF(t, s[1], s[2], s[3], s[4], s[5], s[6],
                  s[7], s[8], s[9], s[10], s[11], s[12], s[13],
                  alpha_t, mach, v_total))
        end

        s = step(s)

        # Renormalize the quaternion: any integrator drifts off the unit sphere.
        # How fast it drifts is also the cheapest accuracy check available — the
        # rule McCoy applies to the direction-cosine formulation is that a
        # deviation past ~1e-5 means the step is too long.
        qn = sqrt(s[7]^2 + s[8]^2 + s[9]^2 + s[10]^2)
        s = (s[1], s[2], s[3], s[4], s[5], s[6],
             s[7]/qn, s[8]/qn, s[9]/qn, s[10]/qn,
             s[11], s[12], s[13])

        t += dt
        n += 1
    end

    # With decimated output the last sample would otherwise fall short of the end.
    if n % sp.record_every != 0
        alpha_t, mach, v_total = diagnostics6(s)
        push!(results, State6DOF(t, s[1], s[2], s[3], s[4], s[5], s[6],
              s[7], s[8], s[9], s[10], s[11], s[12], s[13],
              alpha_t, mach, v_total))
    end

    return results
end

end # module SixDOF


# ──────────────────────────────────────────────────────────────────────────────
# Re-export all sub-modules
# ──────────────────────────────────────────────────────────────────────────────
using .BallisticUtils
using .ReferenceData
using .Atmosphere
using .InteriorBallistics
using .DragModels
using .ExteriorBallistics
using .ReloadingAnalysis
using .SixDOF

export BallisticUtils, ReferenceData, Atmosphere, InteriorBallistics,
       DragModels, ExteriorBallistics, ReloadingAnalysis, SixDOF

end # module CompetitionBallistics


# ═══════════════════════════════════════════════════════════════════════════════
# EXAMPLE USAGE (uncomment to run)
# ═══════════════════════════════════════════════════════════════════════════════
#=
include("CompetitionBallistics.jl")
using .CompetitionBallistics
using .CompetitionBallistics.ExteriorBallistics
using .CompetitionBallistics.ReloadingAnalysis
using .CompetitionBallistics.BallisticUtils

# ── Example 1: 6.5 Creedmoor trajectory to 1000 yd ──
params = ShotParameters(
    mass_grains     = 140.0,
    caliber_in      = 0.264,
    bc              = 0.311,
    drag_model      = :G7,
    muzzle_vel_fps  = 2700.0,
    sight_height_in = 1.5,
    zero_range_yd   = 100.0,
    temp_F          = 70.0,
    pressure_inhg   = 29.92,
    humidity_pct    = 50.0,
    wind_speed_mph  = 10.0,
    wind_angle_deg  = 90.0,
    target_range_yd = 1000.0,
    twist_in        = 8.0,
    bullet_length_in = 1.34,
    latitude_deg    = 45.0,
    azimuth_deg     = 90.0,
)

println("=== 6.5 Creedmoor Trajectory ===")
trajectory_table(params, step_yd=100)

# ── Example 2: Stability check ──
sg = miller_sg(params)
println("\nGyroscopic stability Sg = $(round(sg, digits=2))")

# ── Example 3: Chronograph analysis ──
velocities = [2695, 2702, 2688, 2710, 2699, 2693, 2701, 2707, 2694, 2698]
stats = chronograph_stats(velocities)
println("\nChronograph: mean=$(round(stats.mean,digits=1)) fps, " *
        "ES=$(round(stats.es,digits=1)), SD=$(round(stats.sd,digits=1))")

# ── Example 4: Velocity SD impact ──
vert = velocity_vertical_dispersion(params, stats.sd, range_yd=1000.0)
println("Vertical σ at 1000 yd from velocity SD: $(round(vert,digits=2)) inches")

# ── Example 5: Temperature shift table ──
println("\n=== Temperature Shift (6.5 CM, 1000 yd, σ_T=1.0 fps/°F) ===")
temp_table = temperature_shift_table(params, sigma_T=1.0)
for r in temp_table
    println("  T=$(r.temp_F)°F  v=$(round(r.velocity_fps,digits=0)) fps  " *
            "drop=$(round(r.drop_in,digits=1))\"  shift=$(round(r.shift_moa,digits=1)) MOA")
end

# ── Example 6: BC from two chronographs ──
bc_est = bc_from_two_chronographs(2700.0, 2475.0, 300.0,
            caliber_in=0.264, mass_gr=140.0, drag_model=:G7)
println("\nEstimated BC (G7) from two chronos: $(round(bc_est, digits=3))")

# ── Example 7: Load quality score ──
score = load_consistency_score(stats.es, stats.sd)
println("Load consistency score: $(round(score, digits=1)) / 100")
=#
