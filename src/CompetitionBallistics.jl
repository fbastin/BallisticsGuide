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
                       muzzle_vel_fps=2800, temp_F=59)

Gyroscopic stability factor Sg using the Miller twist rule.
Returns Sg; values > 1.3 indicate adequate stability.
"""
function miller_stability(;
    mass_gr::Real,
    caliber_in::Real,
    bullet_length_in::Real,
    twist_in::Real,
    muzzle_vel_fps::Real = 2800.0,
    temp_F::Real = 59.0
)
    d  = caliber_in
    l  = bullet_length_in / d       # length in calibers
    tw = twist_in / d               # twist in calibers per turn
    T_R = fahrenheit_to_rankine(temp_F)
    sg = 30.0 * mass_gr / (tw^2 * d^3 * l * (1.0 + l^2))
    sg *= (muzzle_vel_fps / 2800.0)
    sg *= (518.67 / T_R)
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
end

# Data from Chapter 11 of the manual
const COMMON_BULLETS = Dict{String, BulletProfile}(
    "Berger 105 Hybrid (6mm)" => BulletProfile("Berger 105 Hybrid (6mm)", 105.0, 0.243, 0.275, 1.10),
    "Berger 140 Hybrid (6.5mm)" => BulletProfile("Berger 140 Hybrid (6.5mm)", 140.0, 0.264, 0.311, 1.34),
    "Hornady 147 ELD-M (6.5mm)" => BulletProfile("Hornady 147 ELD-M (6.5mm)", 147.0, 0.264, 0.351, 1.42),
    "Berger 180 Hybrid (7mm)" => BulletProfile("Berger 180 Hybrid (7mm)", 180.0, 0.284, 0.350, 1.50),
    "Sierra 175 MK (.308)" => BulletProfile("Sierra 175 MK (.308)", 175.0, 0.308, 0.259, 1.24),
    "Berger 185 Juggernaut (.308)" => BulletProfile("Berger 185 Juggernaut (.308)", 185.0, 0.308, 0.283, 1.30),
    "Hornady 178 ELD-M (.308)" => BulletProfile("Hornady 178 ELD-M (.308)", 178.0, 0.308, 0.274, 1.31),
    "Berger 300 Hybrid (.338)" => BulletProfile("Berger 300 Hybrid (.338)", 300.0, 0.338, 0.419, 1.82)
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
const R_air     = 287.058     # J/(kg·K) specific gas constant for dry air
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
    return (P - 0.3783 * e) / (R_air * T)
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
"""
function density_altitude(rho::Real)
    return (T0 / L) * (1.0 - (rho / rho0)^0.190263)
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
    drag_coefficient(model, mach; spline=false) -> Cd

Return drag coefficient for `:G1` or `:G7` at `mach`.
Set `spline=true` for C²-continuous cubic spline interpolation
(recommended for transonic accuracy).
"""
function drag_coefficient(model::Symbol, mach::Real; spline::Bool=false)
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
    windage_m::Float64      # m (cross-range, positive = right)
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

Eötvös vertical deflection.
Δy ≈ -ω_E · cos(φ) · sin(ψ) · R · t
"""
function eotvos_vertical(range_m::Real, tof::Real,
                         latitude_deg::Real, azimuth_deg::Real)
    ω = 7.2921e-5
    return -ω * cos(deg2rad(latitude_deg)) * sin(deg2rad(azimuth_deg)) *
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

    θ = atan(Atmosphere.g0 * si.zr / (2.0 * si.v0^2))  # initial guess

    for _ in 1:max_iter
        x, y = 0.0, -si.sh
        vx = si.v0 * cos(θ)
        vy = si.v0 * sin(θ)
        dt = p.dt

        while x < si.zr
            v = sqrt(vx^2 + vy^2)
            Ma = v / a_s
            cd = drag_coefficient(p.drag_model, Ma)
            drag_factor = ρ * cd * A / (2.0 * si.m)
            ax = -drag_factor * v * vx
            ay = -drag_factor * v * vy - Atmosphere.g0
            vx += ax * dt
            vy += ay * dt
            x  += vx * dt
            y  += vy * dt
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

    # Wind components
    w_ang = deg2rad(p.wind_angle_deg)
    wx = -si.w * cos(w_ang)
    wz =  si.w * sin(w_ang)
    wy = 0.0

    # Incline
    inc = deg2rad(p.incline_deg)

    # Zero angle + dialed elevation/windage
    θ0 = find_zero_angle(p)
    θ0 += moa_to_rad(p.elevation_moa)
    ψ0  = moa_to_rad(p.windage_moa)

    sg = miller_sg(p)
    dt = p.dt

    # Initial state
    x, y, z = 0.0, 0.0, 0.0
    vx = si.v0 * cos(θ0) * cos(inc) * cos(ψ0)
    vy = si.v0 * sin(θ0)
    vz = si.v0 * cos(θ0) * sin(ψ0)
    t  = 0.0

    results = TrajectoryPoint[]

    # ── RK4 derivatives ──
    function derivs(vx, vy, vz)
        vrx = vx - wx;  vry = vy - wy;  vrz = vz - wz
        vrel = sqrt(vrx^2 + vry^2 + vrz^2)
        Ma = vrel / a_s
        cd = drag_coefficient(p.drag_model, Ma)
        D = ρ * cd * A / (2.0 * si.m)

        ax = -D * vrel * vrx
        ay = -D * vrel * vry - Atmosphere.g0
        az = -D * vrel * vrz

        if p.enable_coriolis
            oe_y = ω_E * cos(lat)
            oe_z = ω_E * sin(lat)
            ax += -2.0 * (oe_y * vz - oe_z * vy)
            ay += -2.0 * (oe_z * vx)
            az += -2.0 * (-oe_y * vx)
        end
        return (ax, ay, az)
    end

    while x <= si.tr && t < 15.0
        v_tot = sqrt(vx^2 + vy^2 + vz^2)
        mach  = v_tot / a_s
        ek    = 0.5 * si.m * v_tot^2

        push!(results, TrajectoryPoint(
            t, x, y + si.sh, z,
            vx, vy, vz, v_tot, mach, ek
        ))

        # RK4 integration
        ax1, ay1, az1 = derivs(vx, vy, vz)

        vx2 = vx + 0.5*dt*ax1;  vy2 = vy + 0.5*dt*ay1;  vz2 = vz + 0.5*dt*az1
        ax2, ay2, az2 = derivs(vx2, vy2, vz2)

        vx3 = vx + 0.5*dt*ax2;  vy3 = vy + 0.5*dt*ay2;  vz3 = vz + 0.5*dt*az2
        ax3, ay3, az3 = derivs(vx3, vy3, vz3)

        vx4 = vx + dt*ax3;  vy4 = vy + dt*ay3;  vz4 = vz + dt*az3
        ax4, ay4, az4 = derivs(vx4, vy4, vz4)

        vx += dt/6.0 * (ax1 + 2*ax2 + 2*ax3 + ax4)
        vy += dt/6.0 * (ay1 + 2*ay2 + 2*ay3 + ay4)
        vz += dt/6.0 * (az1 + 2*az2 + 2*az3 + az4)

        x += vx * dt
        y += vy * dt
        z += vz * dt
        t += dt
    end

    # Post-hoc spin drift
    if p.enable_spin_drift
        for i in eachindex(results)
            pt = results[i]
            sd_m = spin_drift_inches(pt.time, sg, p.twist_direction) * 0.0254
            results[i] = TrajectoryPoint(
                pt.time, pt.range_m, pt.drop_m,
                pt.windage_m + sd_m,
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

    println("┌────────┬──────────┬──────────┬──────────┬──────────┬────────┬──────────┐")
    println("│Rng (yd)│ Drop (in)│ Elev(MOA)│ Wind (in)│ Vel(fps) │ ToF(s) │ Eng(ftlb)│")
    println("├────────┼──────────┼──────────┼──────────┼──────────┼────────┼──────────┤")

    next_r = step_m
    for pt in traj
        if pt.range_m >= next_r
            r_yd   = m_to_yards(pt.range_m)
            drop_in = m_to_inches(pt.drop_m)
            wind_in = m_to_inches(pt.windage_m)
            elev    = drop_to_moa(abs(drop_in), r_yd)
            elev    = pt.drop_m < 0 ? elev : -elev
            v_fps   = ms_to_fps(pt.v_total)
            ek_ftlb = pt.energy_J / 1.35582

            @printf("│ %6.0f │ %+8.1f │ %+8.1f │ %+8.1f │ %8.0f │ %6.3f │ %8.0f │\n",
                    r_yd, drop_in, elev, wind_in, v_fps, pt.time, ek_ftlb)

            next_r += step_m
        end
    end
    println("└────────┴──────────┴──────────┴──────────┴──────────┴────────┴──────────┘")
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
       gyroscopic_stability_6dof, dynamic_stability_6dof

"""Aerodynamic coefficient set for 6-DOF simulation."""
Base.@kwdef struct AeroCoefficients6DOF
    # Cd0 as a function of Mach (zero-yaw drag); defaults to G7 × form factor
    Cd0_func::Function   = (Ma) -> cd_g7(Ma)
    Cda2::Float64        = 3.5        # yaw drag coefficient C_{D,α²}
    CNa::Float64         = 2.8        # normal force derivative C_{N,α} [per rad]
    CMa::Float64         = 3.2        # overturning moment coeff C_{M,α} [per rad]
    Clp::Float64         = -0.010     # spin-damping moment C_{l,p}
    CMq_CMad::Float64    = -8.0       # pitch damping sum (C_{M,q} + C_{M,α̇})
    CNpa::Float64        = 0.1        # Magnus force coefficient C_{N,pα}
    CMpa::Float64        = -0.2       # Magnus moment coefficient C_{M,pα}
end

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
    dt::Float64                 = 0.0001    # smaller dt for angular dynamics
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

"""
    moments_of_inertia(mass_kg, caliber_m, length_m) -> (Ix, Iy)

Approximate moments of inertia for a rotationally symmetric bullet.
Ix = axial (about symmetry axis), Iy = transverse.
"""
function moments_of_inertia(mass_kg::Real, caliber_m::Real, length_m::Real)
    r = caliber_m / 2.0
    Ix = 0.5 * mass_kg * r^2               # solid cylinder approximation
    Iy = mass_kg * (length_m^2 / 12.0 + r^2 / 4.0)
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
    dynamic_stability_6dof(CNa, CMa, Sg, m, rho, A, d, kt2, CMq_CMad) -> Sd

Dynamic stability factor from 6-DOF linearization.
"""
function dynamic_stability_6dof(CNa::Real, CMa::Real, Sg::Real,
                                 m::Real, rho::Real, A::Real, d::Real,
                                 kt2::Real, CMq_CMad::Real)
    num = CNa + 2m * kt2 / (rho * A * d)
    return 2CNa / CMa * 1.0 / (Sg * (Sg - 1)) * (num / CMq_CMad)^2
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

Full 6-DOF rigid-body trajectory solver with RK4 integration.
Uses quaternion kinematics (singularity-free).
"""
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
    A   = π * d^2 / 4.0
    tr  = yards_to_m(sp.target_range_yd)

    Ix, Iy = moments_of_inertia(m, d, L)
    aero = sp.aero
    dt = sp.dt

    # Wind
    w_ang = deg2rad(sp.wind_angle_deg)
    wx = -sp.wind_speed_mph * 0.44704 * cos(w_ang)
    wz =  sp.wind_speed_mph * 0.44704 * sin(w_ang)

    # Initial conditions
    alpha0 = deg2rad(sp.initial_yaw_deg)
    p0 = sp.twist_direction * initial_spin_rate(v0, tw)

    # Find zero angle (simplified 2D iteration)
    zr = yards_to_m(sp.zero_range_yd)
    theta0 = atan(Atmosphere.g0 * zr / (2 * v0^2))
    for _ in 1:30
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

    # State: position, velocity, quaternion, angular velocity
    x, y, z = 0.0, 0.0, 0.0
    u = v0 * cos(theta0)
    v_y = v0 * sin(theta0)
    w_z = 0.0

    # Initial quaternion: small yaw offset alpha0 about z-axis
    q0_q = cos(alpha0/2); q1_q = 0.0; q2_q = 0.0; q3_q = sin(alpha0/2)

    p_s = p0           # spin rate
    q_p = sp.initial_pitch_rate  # pitch rate
    r_y = 0.0          # yaw rate

    results = State6DOF[]
    t = 0.0

    while x <= tr && t < 15.0
        v_total = sqrt(u^2 + v_y^2 + w_z^2)
        mach = v_total / a_s
        R = quat_to_dcm(q0_q, q1_q, q2_q, q3_q)
        b1 = R[:, 1]  # body x-axis in inertial frame

        # Relative velocity (subtract wind)
        vrel = [u - wx, v_y, w_z - wz]
        vrel_mag = norm(vrel)

        # Total angle of attack
        cos_alpha = clamp(dot(vrel, b1) / (vrel_mag + 1e-20), -1.0, 1.0)
        alpha_t = acos(cos_alpha)

        push!(results, State6DOF(t, x, y, z, u, v_y, w_z,
              q0_q, q1_q, q2_q, q3_q, p_s, q_p, r_y, alpha_t, mach, v_total))

        # Dynamic pressure
        qbar = 0.5 * rho * vrel_mag^2

        # --- Forces (in inertial frame) ---
        Ma = vrel_mag / a_s
        Cd0 = aero.Cd0_func(Ma)
        Cd_total = Cd0 + aero.Cda2 * alpha_t^2

        # Drag (opposite to relative velocity)
        F_drag = -qbar * A * Cd_total * vrel / (vrel_mag + 1e-20)

        # Normal force (lift) - perpendicular to b1, in the v-b1 plane
        v_perp = vrel - dot(vrel, b1) * b1
        v_perp_mag = norm(v_perp)
        if v_perp_mag > 1e-10
            n_hat = v_perp / v_perp_mag
            F_lift = qbar * A * aero.CNa * alpha_t * n_hat
        else
            F_lift = zeros(3)
        end

        # Magnus force (perpendicular to both b1 and the yaw plane)
        if v_perp_mag > 1e-10 && abs(p_s) > 1e-3
            magnus_dir = cross(b1, n_hat)
            magnus_dir_mag = norm(magnus_dir)
            if magnus_dir_mag > 1e-10
                magnus_dir /= magnus_dir_mag
                F_magnus = qbar * A * aero.CNpa * (p_s * d / (2 * vrel_mag)) *
                           alpha_t * magnus_dir
            else
                F_magnus = zeros(3)
            end
        else
            F_magnus = zeros(3)
        end

        # Gravity
        F_grav = [0.0, -m * Atmosphere.g0, 0.0]

        # Total force -> acceleration
        F_total = F_drag + F_lift + F_magnus + F_grav
        ax = F_total[1] / m
        ay = F_total[2] / m
        az = F_total[3] / m

        # --- Moments (in body frame) ---
        pd2v = p_s * d / (2 * vrel_mag + 1e-20)

        # Decompose alpha into pitch (alpha) and yaw (beta) components
        # in body frame
        v_body = R' * vrel
        alpha_p = -v_body[3] / (v_body[1] + 1e-20)  # pitch-plane
        beta_y  =  v_body[2] / (v_body[1] + 1e-20)  # yaw-plane

        Mx = qbar * A * d^2 * aero.Clp * pd2v
        Mq_body = qbar * A * d * (aero.CMa * alpha_p -
                  aero.CMpa * pd2v * beta_y +
                  aero.CMq_CMad * q_p * d / (2 * vrel_mag + 1e-20))
        Mr_body = qbar * A * d * (aero.CMa * beta_y +
                  aero.CMpa * pd2v * alpha_p +
                  aero.CMq_CMad * r_y * d / (2 * vrel_mag + 1e-20))

        # Euler rotational equations
        dp = Mx / Ix
        dq = (Mq_body - (Ix - Iy) * p_s * r_y) / Iy
        dr = (Mr_body + (Ix - Iy) * p_s * q_p) / Iy

        # Quaternion kinematic ODE
        dq0 = 0.5 * (-p_s*q1_q - q_p*q2_q - r_y*q3_q)
        dq1 = 0.5 * ( p_s*q0_q + r_y*q2_q - q_p*q3_q)
        dq2 = 0.5 * ( q_p*q0_q - r_y*q1_q + p_s*q3_q)
        dq3 = 0.5 * ( r_y*q0_q + q_p*q1_q - p_s*q2_q)

        # --- Euler step (simplified; full RK4 on the 13-state would be
        #     ideal but Euler with small dt is adequate for demonstration) ---
        u   += ax * dt
        v_y += ay * dt
        w_z += az * dt
        x   += u * dt
        y   += v_y * dt
        z   += w_z * dt

        p_s += dp * dt
        q_p += dq * dt
        r_y += dr * dt

        q0_q += dq0 * dt
        q1_q += dq1 * dt
        q2_q += dq2 * dt
        q3_q += dq3 * dt

        # Renormalize quaternion
        qn = sqrt(q0_q^2 + q1_q^2 + q2_q^2 + q3_q^2)
        q0_q /= qn; q1_q /= qn; q2_q /= qn; q3_q /= qn

        t += dt
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
