### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 00000001-0000-0000-0000-000000000001
md"""
# Competition Rifle Ballistics — Step-by-Step Tutorial

**Companion notebook for *Competition Rifle Ballistics: A Comprehensive Mathematical Manual***

This notebook walks through every module of the `CompetitionBallistics.jl` library, from basic unit conversions to full 6-DOF trajectory simulation. Each section builds on the previous one, mirroring the workflow a competitive shooter would follow: select components → model atmosphere → compute trajectory → analyse loads → refine.

> **How to run:** Place `CompetitionBallistics.jl` in the same directory as this notebook, then execute cells in order.
"""

# ╔═╡ 00000002-0000-0000-0000-000000000001
md"""
---
## 0. Load the Library
"""

# ╔═╡ 00000003-0000-0000-0000-000000000001
include("CompetitionBallistics.jl")

using .CompetitionBallistics.BallisticUtils
using .CompetitionBallistics.Atmosphere
using .CompetitionBallistics.InteriorBallistics
using .CompetitionBallistics.DragModels
using .CompetitionBallistics.ExteriorBallistics
using .CompetitionBallistics.ReloadingAnalysis
using .CompetitionBallistics.SixDOF

using Printf
using Statistics

# ╔═╡ 00000004-0000-0000-0000-000000000001
md"""
---
## 1. Unit Conversions (`BallisticUtils`)

The library works in SI internally but accepts and returns Imperial units, which are standard in the shooting world.
"""

# ╔═╡ 00000005-0000-0000-0000-000000000001
md"""
### 1.1 Mass, velocity, and length
"""

# ╔═╡ 00000006-0000-0000-0000-000000000001
begin
    # A 140 grain 6.5mm bullet
    m_gr = 140.0
    m_kg = grains_to_kg(m_gr)
    m_g  = grains_to_grams(m_gr)
    println("140 grains = $(round(m_kg*1e3, digits=3)) grams = $(round(m_kg, digits=6)) kg")

    # Muzzle velocity
    v_fps = 2700.0
    v_ms  = fps_to_ms(v_fps)
    println("2700 fps = $(round(v_ms, digits=1)) m/s")

    # Barrel length
    L_in = 24.0
    L_m  = inches_to_m(L_in)
    println("24 inches = $(round(L_m, digits=4)) m = $(round(inches_to_mm(L_in), digits=1)) mm")
end

# ╔═╡ 00000007-0000-0000-0000-000000000001
md"""
### 1.2 Angular conversions
"""

# ╔═╡ 00000008-0000-0000-0000-000000000001
begin
    println("1 MOA = $(round(moa_to_mil(1.0), digits=4)) mil")
    println("1 mil = $(round(mil_to_moa(1.0), digits=3)) MOA")
    println("1 MOA = $(round(moa_to_rad(1.0)*1e6, digits=2)) μrad")

    # A 10-inch drop at 500 yards is how many MOA?
    drop_moa = drop_to_moa(10.0, 500.0)
    drop_mil = drop_to_mil(10.0, 500.0)
    println("\n10 inches of drop at 500 yd = $(round(drop_moa, digits=2)) MOA = $(round(drop_mil, digits=2)) mil")
end

# ╔═╡ 00000009-0000-0000-0000-000000000001
md"""
### 1.3 Sectional density, form factor, and kinetic energy
"""

# ╔═╡ 00000010-0000-0000-0000-000000000001
begin
    cal = 0.264  # 6.5mm in inches
    sd_val = sectional_density(m_gr, cal)
    println("Sectional density (140 gr / .264 cal) = $(round(sd_val, digits=3)) lb/in²")

    bc_g7 = 0.311
    ff = form_factor(m_gr, cal, bc_g7)
    println("Form factor (G7 BC=0.311) = $(round(ff, digits=3))")
    println("  → form factor near 1.0 means the bullet closely matches the G7 reference shape")

    ke_ftlb = kinetic_energy_ftlbs(m_gr, v_fps)
    ke_J    = kinetic_energy_J(m_kg, v_ms)
    println("\nKinetic energy at $(v_fps) fps: $(round(ke_ftlb, digits=0)) ft·lbf = $(round(ke_J, digits=0)) J")
end

# ╔═╡ 00000011-0000-0000-0000-000000000001
md"""
### 1.4 Gyroscopic stability (Miller formula)
"""

# ╔═╡ 00000012-0000-0000-0000-000000000001
begin
    sg = miller_stability(
        mass_gr = 140.0,
        caliber_in = 0.264,
        bullet_length_in = 1.34,
        twist_in = 8.0,
        muzzle_vel_fps = 2700.0,
        temp_F = 59.0,
    )
    println("Miller Sg = $(round(sg, digits=2))")
    println(sg > 1.5 ? "  ✓ Well stabilized" :
            sg > 1.0 ? "  ⚠ Marginal — consider faster twist" :
                        "  ✗ UNSTABLE — bullet will tumble")

    # Try a heavier bullet in the same barrel
    sg_147 = miller_stability(
        mass_gr=147.0, caliber_in=0.264, bullet_length_in=1.42,
        twist_in=8.0, muzzle_vel_fps=2600.0, temp_F=59.0,
    )
    println("\n147 gr ELD-M in 1:8 twist: Sg = $(round(sg_147, digits=2))")

    # Greenhill's formula for comparison
    twist_gh = greenhill_twist(0.264, 1.34, C=150.0)
    println("Greenhill recommended twist: 1:$(round(twist_gh, digits=1)) inches")
end

# ╔═╡ 00000013-0000-0000-0000-000000000001
md"""
---
## 2. Atmospheric Model (`Atmosphere`)

Accurate air density is critical: it directly scales the drag force. A 10% density error translates to ≈10% drag error.
"""

# ╔═╡ 00000014-0000-0000-0000-000000000001
md"""
### 2.1 Standard atmosphere
"""

# ╔═╡ 00000015-0000-0000-0000-000000000001
begin
    println("=== ICAO Standard Sea-Level Conditions ===")
    println("T₀ = $(Atmosphere.T0) K ($(round(kelvin_to_fahrenheit(Atmosphere.T0), digits=1)) °F)")
    println("P₀ = $(Atmosphere.P0) Pa ($(round(pa_to_inhg(Atmosphere.P0), digits=2)) inHg)")
    println("ρ₀ = $(Atmosphere.rho0) kg/m³")
    println("a₀ = $(round(speed_of_sound(Atmosphere.T0), digits=1)) m/s ($(round(ms_to_fps(speed_of_sound(Atmosphere.T0)), digits=0)) fps)")

    println("\n=== At 5000 ft altitude ===")
    h = 5000 * 0.3048  # feet to meters
    T_5k = std_temperature(h)
    P_5k = std_pressure(h)
    rho_5k = air_density(P=P_5k, T=T_5k)
    println("T = $(round(T_5k, digits=1)) K ($(round(kelvin_to_fahrenheit(T_5k), digits=1)) °F)")
    println("P = $(round(P_5k, digits=0)) Pa ($(round(pa_to_inhg(P_5k), digits=2)) inHg)")
    println("ρ = $(round(rho_5k, digits=4)) kg/m³ ($(round(100*rho_5k/Atmosphere.rho0, digits=1))% of sea level)")
end

# ╔═╡ 00000016-0000-0000-0000-000000000001
md"""
### 2.2 Humidity correction and density altitude
"""

# ╔═╡ 00000017-0000-0000-0000-000000000001
begin
    # Summer day in Phoenix: 105°F, 29.80 inHg, 20% humidity
    T_phx = fahrenheit_to_kelvin(105.0)
    P_phx = inhg_to_pa(29.80)
    rho_dry = air_density(P=P_phx, T=T_phx, H=0.0)
    rho_humid = air_density(P=P_phx, T=T_phx, H=20.0)
    da = density_altitude(rho_humid)

    println("=== Phoenix summer day: 105°F, 29.80 inHg, 20% RH ===")
    println("Dry air density:   $(round(rho_dry, digits=4)) kg/m³")
    println("Humid air density: $(round(rho_humid, digits=4)) kg/m³")
    println("Humidity effect:   $(round((rho_humid - rho_dry)/rho_dry * 100, digits=2))% (humid air is lighter!)")
    println("Density altitude:  $(round(da, digits=0)) m ($(round(da/0.3048, digits=0)) ft)")
    println("  → Bullets will shoot FLATTER than at sea level standard")

    # Compare: Winter in Quebec, -20°F, 30.10 inHg, 40% RH
    T_qc = fahrenheit_to_kelvin(-20.0)
    P_qc = inhg_to_pa(30.10)
    rho_qc = air_density(P=P_qc, T=T_qc, H=40.0)
    da_qc = density_altitude(rho_qc)
    println("\n=== Quebec winter: -20°F, 30.10 inHg, 40% RH ===")
    println("Air density: $(round(rho_qc, digits=4)) kg/m³ ($(round(100*rho_qc/Atmosphere.rho0, digits=1))% of std)")
    println("Density altitude: $(round(da_qc, digits=0)) m ($(round(da_qc/0.3048, digits=0)) ft)")
    println("Drag ratio Phoenix/Quebec: $(round(rho_humid/rho_qc, digits=3))")
end

# ╔═╡ 00000018-0000-0000-0000-000000000001
md"""
---
## 3. Drag Models (`DragModels`)

The G1 and G7 BRL tables define the standard reference projectiles. G7 is preferred for modern VLD competition bullets.
"""

# ╔═╡ 00000019-0000-0000-0000-000000000001
begin
    println("=== Drag coefficients at selected Mach numbers ===")
    println("  Mach  │   G1 Cd  │   G7 Cd  │  G1/G7 ratio")
    println("────────┼──────────┼──────────┼──────────────")
    for Ma in [0.5, 0.8, 0.95, 1.0, 1.1, 1.5, 2.0, 2.5, 3.0]
        g1 = cd_g1(Ma)
        g7 = cd_g7(Ma)
        @printf("  %4.2f  │  %6.4f  │  %6.4f  │    %5.2f\n", Ma, g1, g7, g1/g7)
    end
    println("\nNote: G1 drag is ~2× higher than G7 at all Mach numbers,")
    println("which is why G1 BC values are ~2× larger than G7 BC values")
    println("for the SAME bullet.")
end

# ╔═╡ 00000020-0000-0000-0000-000000000001
md"""
---
## 4. Interior Ballistics (`InteriorBallistics`)

Before the bullet flies, we need to understand what happens inside the barrel.
"""

# ╔═╡ 00000021-0000-0000-0000-000000000001
md"""
### 4.1 Le Duc velocity model
"""

# ╔═╡ 00000022-0000-0000-0000-000000000001
begin
    # Model a .308 Win: 175 gr bullet, 44 gr Varget, 24" barrel
    # Le Duc parameters (fitted to typical .308 data)
    a_leduc = 1100.0   # theoretical max velocity [m/s] ≈ 3609 fps
    b_leduc = 0.10     # shape parameter [m]
    barrel_travel = 0.55  # ~21.5 inches of bullet travel

    println("=== Le Duc Velocity Profile (.308 Win 175 gr / Varget) ===")
    println("  Travel (in) │ Velocity (fps) │ Pressure (psi)")
    println("──────────────┼────────────────┼────────────────")

    m_308 = grains_to_kg(175.0)
    m_chg = grains_to_kg(44.0)
    d_308 = inches_to_m(0.308)

    for x_in in [2, 4, 6, 8, 10, 14, 18, 22]
        x_m = inches_to_m(x_in)
        v = leduc_velocity(x_m, a=a_leduc, b=b_leduc)
        P = leduc_pressure(x_m, a=a_leduc, b=b_leduc,
                          m_bullet=m_308, m_charge=m_chg, d_bore=d_308)
        @printf("    %5.1f     │    %7.0f     │    %7.0f\n",
                Float64(x_in), ms_to_fps(v), P / 6894.76)  # Pa to psi
    end

    P_max = peak_pressure(a=a_leduc, b=b_leduc,
                         m_bullet=m_308, m_charge=m_chg, d_bore=d_308)
    @printf("\nPeak pressure: %.0f psi (at x = %.1f inches)\n",
            P_max / 6894.76, m_to_inches(b_leduc/2))
end

# ╔═╡ 00000023-0000-0000-0000-000000000001
md"""
### 4.2 Temperature sensitivity and barrel length correction
"""

# ╔═╡ 00000024-0000-0000-0000-000000000001
begin
    v0_base = 2700.0  # fps at 70°F
    println("=== Temperature Sensitivity ===")
    println("  Temp (°F)  │  Δv (Varget, σ=0.5)  │  Δv (IMR 4064, σ=1.8)")
    println("─────────────┼──────────────────────┼──────────────────────")
    for T_F in [0, 20, 40, 59, 70, 90, 110]
        dT = T_F - 70.0
        v_stable   = muzzle_velocity_temp_correction(v0_base, dT, sigma_T=0.5)
        v_unstable = muzzle_velocity_temp_correction(v0_base, dT, sigma_T=1.8)
        @printf("    %4.0f     │      %+6.0f fps       │      %+6.0f fps\n",
                Float64(T_F), v_stable - v0_base, v_unstable - v0_base)
    end

    println("\n=== Barrel Length Correction ===")
    v_ref = 2700.0  # fps from a 24" barrel
    for L in [20, 22, 24, 26, 28, 30]
        v_new = barrel_length_correction(v_ref, Float64(L), 24.0, exponent=0.27)
        @printf("  %d\" barrel → %.0f fps (Δ = %+.0f)\n", L, v_new, v_new - v_ref)
    end
end

# ╔═╡ 00000025-0000-0000-0000-000000000001
md"""
---
## 5. Exterior Ballistics — 3-DOF Trajectory (`ExteriorBallistics`)

This is the core of the library. We define a shot and solve the full trajectory.
"""

# ╔═╡ 00000026-0000-0000-0000-000000000001
md"""
### 5.1 Define a shot: 6.5 Creedmoor at 1000 yards
"""

# ╔═╡ 00000027-0000-0000-0000-000000000001
params_65cm = ShotParameters(
    # Bullet: Berger 140 Hybrid Target
    mass_grains     = 140.0,
    caliber_in      = 0.264,
    bc              = 0.311,
    drag_model      = :G7,

    # Launch
    muzzle_vel_fps  = 2700.0,
    sight_height_in = 1.5,
    zero_range_yd   = 100.0,

    # Atmosphere: standard day
    temp_F          = 59.0,
    pressure_inhg   = 29.92,
    humidity_pct    = 0.0,

    # Wind: 10 mph full-value crosswind from the right
    wind_speed_mph  = 10.0,
    wind_angle_deg  = 90.0,

    # Target
    target_range_yd = 1000.0,

    # Spin: 1:8 RH twist
    twist_in        = 8.0,
    twist_direction = 1,
    bullet_length_in = 1.34,

    # Location: mid-latitude
    latitude_deg    = 45.0,
    azimuth_deg     = 90.0,    # firing East
    enable_coriolis = true,
    enable_spin_drift = true,
)

# ╔═╡ 00000028-0000-0000-0000-000000000001
md"""
### 5.2 Solve and print the trajectory table
"""

# ╔═╡ 00000029-0000-0000-0000-000000000001
trajectory_table(params_65cm, step_yd=100)

# ╔═╡ 00000030-0000-0000-0000-000000000001
md"""
### 5.3 Extract specific data from the trajectory
"""

# ╔═╡ 00000031-0000-0000-0000-000000000001
begin
    traj = solve_trajectory(params_65cm)
    println("Total trajectory points: $(length(traj))")
    println("Time of flight to $(round(m_to_yards(traj[end].range_m), digits=0)) yd: $(round(traj[end].time, digits=3)) s")

    # Find the point closest to 600 yards
    target_m = yards_to_m(600.0)
    pt_600 = traj[findfirst(p -> p.range_m >= target_m, traj)]

    println("\n=== At 600 yards ===")
    println("  Velocity:  $(round(ms_to_fps(pt_600.v_total), digits=0)) fps")
    println("  Mach:      $(round(pt_600.mach, digits=3))")
    println("  Drop:      $(round(m_to_inches(pt_600.drop_m), digits=1)) inches")
    println("  Windage:   $(round(m_to_inches(pt_600.windage_m), digits=1)) inches (includes spin drift)")
    println("  Energy:    $(round(pt_600.energy_J / 1.35582, digits=0)) ft·lbf")
    println("  ToF:       $(round(pt_600.time, digits=3)) s")
end

# ╔═╡ 00000032-0000-0000-0000-000000000001
md"""
### 5.4 Quick estimates: wind lag rule, Coriolis, spin drift
"""

# ╔═╡ 00000033-0000-0000-0000-000000000001
begin
    R_m = yards_to_m(1000.0)
    tof = traj[end].time
    v0_ms = fps_to_ms(2700.0)
    w_cross = 10.0 * 0.44704  # 10 mph in m/s

    # Wind deflection via lag rule
    wind_lag = wind_deflection_lag_rule(w_cross, tof, R_m, v0_ms)
    println("=== Quick Estimates at 1000 yd ===")
    println("Wind deflection (lag rule): $(round(m_to_inches(wind_lag), digits=1)) inches")

    # Coriolis
    cor_hz = coriolis_horizontal(R_m, tof, 45.0)
    println("Coriolis horizontal:        $(round(m_to_inches(cor_hz), digits=1)) inches (right in N. hemisphere)")

    # Eötvös (firing East)
    eot_vt = eotvos_vertical(R_m, tof, 45.0, 90.0)
    println("Eötvös vertical:            $(round(m_to_inches(eot_vt), digits=1)) inches (up when firing East)")

    # Spin drift
    sg_val = miller_sg(params_65cm)
    sd_in = spin_drift_inches(tof, sg_val, 1)
    println("Spin drift (1:8 RH):        $(round(sd_in, digits=1)) inches right")
end

# ╔═╡ 00000034-0000-0000-0000-000000000001
md"""
### 5.5 Compare cartridges: .308 Win vs 6.5 CM vs .223 Rem
"""

# ╔═╡ 00000035-0000-0000-0000-000000000001
begin
    cartridges = [
        ("6.5 CM 140 Hybrid",  ShotParameters(mass_grains=140.0, caliber_in=0.264, bc=0.311,
            drag_model=:G7, muzzle_vel_fps=2700.0, twist_in=8.0, bullet_length_in=1.34,
            wind_speed_mph=10.0, wind_angle_deg=90.0, target_range_yd=1000.0)),
        (".308 Win 175 SMK",   ShotParameters(mass_grains=175.0, caliber_in=0.308, bc=0.259,
            drag_model=:G7, muzzle_vel_fps=2600.0, twist_in=10.0, bullet_length_in=1.24,
            wind_speed_mph=10.0, wind_angle_deg=90.0, target_range_yd=1000.0)),
        (".223 Rem 77 SMK",    ShotParameters(mass_grains=77.0, caliber_in=0.224, bc=0.190,
            drag_model=:G7, muzzle_vel_fps=2750.0, twist_in=7.0, bullet_length_in=1.00,
            wind_speed_mph=10.0, wind_angle_deg=90.0, target_range_yd=1000.0)),
    ]

    println("=== Cartridge Comparison at 1000 yd, 10 mph crosswind ===")
    println("  Cartridge            │ Vel(fps)│ Drop(in)│ Wind(in)│ ToF(s) │ Eng(ft·lb)")
    println("───────────────────────┼─────────┼─────────┼─────────┼────────┼───────────")

    for (name, p) in cartridges
        t = solve_trajectory(p)
        pt = t[end]
        v = round(ms_to_fps(pt.v_total), digits=0)
        d = round(m_to_inches(pt.drop_m), digits=1)
        w = round(m_to_inches(pt.windage_m), digits=1)
        tf = round(pt.time, digits=3)
        e = round(pt.energy_J / 1.35582, digits=0)
        @printf("  %-21s │ %7.0f │ %+7.1f │ %+7.1f │ %6.3f │ %9.0f\n", name, v, d, w, tf, e)
    end
    println("\n  → The 6.5 CM has 30% less wind deflection than .308 and 50% less than .223")
end

# ╔═╡ 00000036-0000-0000-0000-000000000001
md"""
---
## 6. Reloading Analysis (`ReloadingAnalysis`)

Now we analyse handloading data using the ballistic models.
"""

# ╔═╡ 00000037-0000-0000-0000-000000000001
md"""
### 6.1 Chronograph data analysis
"""

# ╔═╡ 00000038-0000-0000-0000-000000000001
begin
    # Simulated 10-round string from the range
    velocities = [2695, 2702, 2688, 2710, 2699, 2693, 2701, 2707, 2694, 2698]

    stats = chronograph_stats(velocities)
    println("=== Chronograph Results (10 rounds, 6.5 CM 140 Hybrid / H4350) ===")
    println("  Mean:    $(round(stats.mean, digits=1)) fps")
    println("  ES:      $(round(stats.es, digits=1)) fps")
    println("  SD:      $(round(stats.sd, digits=1)) fps")
    println("  CV:      $(round(stats.cv_pct, digits=3))%")
    println("  Min/Max: $(round(stats.min, digits=0)) / $(round(stats.max, digits=0)) fps")

    score = load_consistency_score(stats.es, stats.sd)
    println("\n  Load quality score: $(round(score, digits=1)) / 100")
    println(score > 85 ? "  ✓ Excellent competition load" :
            score > 70 ? "  ⚠ Acceptable, but room for improvement" :
                         "  ✗ Needs work — check charge consistency")
end

# ╔═╡ 00000039-0000-0000-0000-000000000001
md"""
### 6.2 Velocity SD impact on groups at distance
"""

# ╔═╡ 00000040-0000-0000-0000-000000000001
begin
    println("=== Vertical Dispersion from Velocity SD ===")
    println("  Range (yd) │ σ_vertical (in) │ MOA equivalent")
    println("─────────────┼─────────────────┼───────────────")

    for R in [300, 500, 600, 800, 1000, 1200]
        sigma_y = velocity_vertical_dispersion(params_65cm, stats.sd, range_yd=Float64(R))
        moa_eq = drop_to_moa(sigma_y, Float64(R))
        @printf("    %5d    │     %6.2f       │    %5.2f\n", R, sigma_y, moa_eq)
    end
    println("\n  With SD=$(round(stats.sd, digits=1)) fps:")
    println("  • At 600 yd (PRS): vertical σ < 1\" — negligible on target")
    println("  • At 1000 yd (F-Class): vertical σ ≈ 1.5\" — matters on the X-ring")
end

# ╔═╡ 00000041-0000-0000-0000-000000000001
md"""
### 6.3 Satterlee ladder test analysis
"""

# ╔═╡ 00000042-0000-0000-0000-000000000001
begin
    # Ladder test data: 1 round per charge, 0.3 gr increments
    charges = [39.5, 39.8, 40.1, 40.4, 40.7, 41.0, 41.3, 41.6, 41.9, 42.2]
    ladder_v = [2580, 2610, 2638, 2665, 2680, 2690, 2695, 2708, 2730, 2758]

    result = ladder_test_analysis(charges, ladder_v)

    println("=== Satterlee Ladder Test (6.5 CM / H4350 / 140 Berger) ===")
    println("  Charge (gr) │ Velocity (fps) │ Slope (fps/gr)")
    println("──────────────┼────────────────┼───────────────")
    for i in eachindex(result.charges)
        slope_str = i < length(result.charges) ?
            @sprintf("%+7.1f", result.slopes[i]) : "    ---"
        @printf("    %5.1f     │     %5.0f      │  %s\n",
                result.charges[i], result.velocities[i], slope_str)
    end

    println("\n  ★ Velocity node found at $(round(result.node_charge, digits=2)) grains")
    println("    (slope = $(round(result.node_slope, digits=1)) fps/gr — flattest region)")
    println("    Expected velocity: $(round(result.node_velocity, digits=0)) fps")
    println("    → This charge weight is least sensitive to small charge variations")
end

# ╔═╡ 00000043-0000-0000-0000-000000000001
md"""
### 6.4 Temperature shift table
"""

# ╔═╡ 00000044-0000-0000-0000-000000000001
begin
    println("=== Season Temperature Shift at 1000 yd ===")
    println("  (6.5 CM, 140 Berger, H4350 σ_T=0.5 fps/°F)")
    println("  Temp (°F)  │  MV (fps)  │  Drop Shift (MOA)")
    println("─────────────┼────────────┼──────────────────")

    shifts = temperature_shift_table(params_65cm,
        temps_F=collect(0.0:20.0:120.0),
        range_yd=1000.0,
        sigma_T=0.5)

    ref_moa = shifts[4].shift_moa  # 59°F as reference
    for s in shifts
        @printf("    %5.0f    │   %6.0f   │     %+5.1f\n",
                s.temp_F, s.velocity_fps, s.shift_moa - ref_moa)
    end
    println("\n  Total swing 0°F to 120°F: $(round(abs(shifts[end].shift_moa - shifts[1].shift_moa), digits=1)) MOA")
    println("  → With a temp-stable powder, you need ≈ 2-3 MOA of DOPE adjustment across seasons")
end

# ╔═╡ 00000045-0000-0000-0000-000000000001
md"""
### 6.5 BC estimation from two chronographs
"""

# ╔═╡ 00000046-0000-0000-0000-000000000001
begin
    # Muzzle chrono reads 2700 fps, downrange chrono at 100 yd reads 2580 fps
    bc_est = bc_from_two_chronographs(
        2700.0, 2580.0, 300.0,   # v1, v2, distance in feet
        caliber_in=0.264, mass_gr=140.0, drag_model=:G7,
        temp_F=70.0, pressure_inhg=29.85)

    println("=== BC Estimation from Two Chronographs ===")
    println("  Muzzle:    2700 fps")
    println("  At 100 yd: 2580 fps")
    println("  Estimated BC (G7): $(round(bc_est, digits=3))")
    println("  Published BC (G7): 0.311 (Berger 140 Hybrid)")
    println("  Difference: $(round(abs(bc_est - 0.311)/0.311*100, digits=1))%")
end

# ╔═╡ 00000047-0000-0000-0000-000000000001
md"""
### 6.6 Barrel time and pressure estimation
"""

# ╔═╡ 00000048-0000-0000-0000-000000000001
begin
    bt = optimal_barrel_time(2700.0, 24.0)
    println("=== Barrel Time ===")
    println("  Estimated barrel time (24\" at 2700 fps): $(round(bt, digits=3)) ms")
    println("  (At ~15 kHz barrel resonance, this is $(round(bt/1000 * 15000, digits=1)) oscillation cycles)")

    bore_area = π * (0.264/2)^2  # in²
    p_est = pressure_estimate_psi(
        charge_gr=41.0, bullet_weight_gr=140.0,
        bore_area_in2=bore_area, barrel_length_in=22.0,
        muzzle_vel_fps=2700.0)
    println("\n=== Rough Pressure Estimate (6.5 CM) ===")
    println("  Estimated peak pressure: $(round(p_est, digits=0)) psi")
    println("  SAAMI MAP: 62,000 psi")
    println("  Margin: $(round((62000 - p_est)/62000*100, digits=1))%")
end

# ╔═╡ 00000049-0000-0000-0000-000000000001
md"""
---
## 7. Six-Degree-of-Freedom Model (`SixDOF`)

The 6-DOF model resolves the bullet's angular motion: yaw, precession, nutation, and spin decay. This predicts spin drift from first principles rather than using the Litz empirical formula.

!!! note "Why the 3-DOF solver is enough for the shooting this library targets"
    Reach for this module for stability and flight-dynamic work — not because the 3-DOF
    trajectory is suspect. McCoy (*Modern Exterior Ballistics*, 2nd ed., §9.6) is explicit:
    6-DOF trajectories **are not required for routine work in exterior ballistics**. If the
    total angle of attack stays small everywhere along the flight path, a point-mass
    trajectory is "often sufficiently accurate for all practical purposes". The criterion is
    the **yaw level, not the range**: 6-DOF becomes necessary for large-yaw flight — an
    artillery or mortar shell at high quadrant elevation, or a projectile launched sidewise
    into a several-hundred-mph crosswind. McCoy's own worked example, a .308", 168 gr Sierra
    International match bullet fired flat to 1000 yards, never exceeds 5° of pitch and yaw,
    and he notes it "can actually be done by simpler methods". Sport shooting, long range
    included, sits entirely on the small-yaw side of that line.

!!! warning "Garbage in, garbage out"
    A 6-DOF run is only as good as its aerodynamic coefficients, and the coefficients set up
    in §7.1 below are **library defaults, not measurements for this bullet**. The full set —
    normal force, overturning moment, pitch and roll damping, Magnus force and moment — comes
    from spark-range or wind-tunnel reduction and is published for almost no commercial
    bullet; MCDRAG returns the zero-yaw drag coefficient alone. A 6-DOF trajectory fed with
    placeholder coefficients is worth less than a 3-DOF trajectory fed with a measured BC.
    McCoy's own summary of the method carries the practitioners' motto: *GI-GO*.
"""

# ╔═╡ 00000050-0000-0000-0000-000000000001
md"""
### 7.1 Set up a 6-DOF simulation
"""

# ╔═╡ 00000051-0000-0000-0000-000000000001
begin
    # Define aerodynamic coefficients for a 6.5mm 140gr match bullet
    aero = AeroCoefficients6DOF(
        Cd0_func  = (Ma) -> cd_g7(Ma) * 1.02,  # slight form factor correction
        Cda2      = 3.5,
        CNa       = 2.8,
        CMa       = 3.2,
        Clp       = -0.010,
        CMq_CMad  = -8.0,
        CNpa      = 0.1,
        CMpa      = -0.2,
    )

    params6 = ShotParameters6DOF(
        mass_grains     = 140.0,
        caliber_in      = 0.264,
        bullet_length_in = 1.34,
        muzzle_vel_fps  = 2700.0,
        twist_in        = 8.0,
        twist_direction = 1,
        initial_yaw_deg = 1.0,    # 1° initial disturbance
        temp_F          = 59.0,
        pressure_inhg   = 29.92,
        target_range_yd = 600.0,  # shorter range for 6-DOF demo
        aero            = aero,
        dt              = 0.00005, # fine time step for angular dynamics
    )

    println("=== 6-DOF Parameters ===")
    m_6 = grains_to_kg(140.0)
    d_6 = inches_to_m(0.264)
    L_6 = inches_to_m(1.34)
    Ix, Iy = moments_of_inertia(m_6, d_6, L_6)
    p0 = initial_spin_rate(fps_to_ms(2700.0), inches_to_m(8.0))

    println("  Axial MOI (Ix):     $(round(Ix*1e6, digits=3)) × 10⁻⁶ kg·m²")
    println("  Transverse MOI (Iy): $(round(Iy*1e6, digits=3)) × 10⁻⁶ kg·m²")
    println("  Iy/Ix ratio:        $(round(Iy/Ix, digits=2))")
    println("  Initial spin rate:  $(round(p0, digits=0)) rad/s ($(round(p0/(2π)*60, digits=0)) RPM)")
end

# ╔═╡ 00000052-0000-0000-0000-000000000001
md"""
### 7.2 Run the 6-DOF solver and examine yaw history
"""

# ╔═╡ 00000053-0000-0000-0000-000000000001
begin
    println("Running 6-DOF solver (this may take a few seconds)...")
    traj6 = solve_6dof(params6)
    println("  Done. $(length(traj6)) time steps computed.")
    println("  Final range: $(round(m_to_yards(traj6[end].x), digits=0)) yd")

    # Sample every 2000th point for display
    println("\n=== 6-DOF Trajectory Summary ===")
    println("  Range(yd) │ Vel(fps) │ Mach  │ AoA(deg)  │ Spin(RPM)  │ Drop(in)")
    println("────────────┼──────────┼───────┼───────────┼────────────┼─────────")

    step = max(1, length(traj6) ÷ 12)
    for i in 1:step:length(traj6)
        s = traj6[i]
        r_yd = m_to_yards(s.x)
        v_fps = ms_to_fps(s.v_total)
        aoa_deg = rad2deg(s.alpha_t)
        rpm = abs(s.p_spin) / (2π) * 60
        drop_in = m_to_inches(s.y)
        @printf("   %7.0f  │ %7.0f  │ %5.3f │  %7.4f  │ %10.0f │ %+7.1f\n",
                r_yd, v_fps, s.mach, aoa_deg, rpm, drop_in)
    end
end

# ╔═╡ 00000054-0000-0000-0000-000000000001
md"""
### 7.3 Stability analysis from 6-DOF theory
"""

# ╔═╡ 00000055-0000-0000-0000-000000000001
begin
    rho_std = Atmosphere.rho0
    A_6 = π * d_6^2 / 4
    CMa = 3.2

    Sg_exact = gyroscopic_stability_6dof(Ix, p0, rho_std, A_6, d_6, Iy,
                                          fps_to_ms(2700.0), CMa)
    Sg_miller = miller_stability(mass_gr=140.0, caliber_in=0.264,
                                  bullet_length_in=1.34, twist_in=8.0,
                                  muzzle_vel_fps=2700.0)

    println("=== Stability Comparison ===")
    println("  Exact Sg (6-DOF):    $(round(Sg_exact, digits=3))")
    println("  Miller Sg (approx):  $(round(Sg_miller, digits=3))")
    println("  Difference:          $(round(abs(Sg_exact-Sg_miller)/Sg_exact*100, digits=1))%")
    println("\n  The Miller formula is a good approximation when")
    println("  the aerodynamic coefficients are not precisely known.")
end

# ╔═╡ 00000056-0000-0000-0000-000000000001
md"""
### 7.4 Spin drift comparison: 6-DOF vs empirical
"""

# ╔═╡ 00000057-0000-0000-0000-000000000001
begin
    # Extract lateral drift from 6-DOF at ~600 yd
    pt6_end = traj6[end]
    drift_6dof_in = m_to_inches(pt6_end.z)

    # Empirical Litz formula
    drift_litz_in = spin_drift_inches(pt6_end.t, Sg_miller, 1)

    println("=== Spin Drift at $(round(m_to_yards(pt6_end.x), digits=0)) yd ===")
    println("  6-DOF model:     $(round(drift_6dof_in, digits=2)) inches")
    println("  Litz empirical:  $(round(drift_litz_in, digits=2)) inches")
    println("\n  Note: the 6-DOF drift includes yaw-induced lateral forces")
    println("  and Magnus effects that the empirical formula averages out.")
    println("  Agreement is typically within 20% for well-stabilized bullets.")
end

# ╔═╡ 00000058-0000-0000-0000-000000000001
md"""
---
## 8. Putting It All Together: Match-Day Workflow

A complete workflow from component selection to DOPE card generation.
"""

# ╔═╡ 00000059-0000-0000-0000-000000000001
begin
    println("╔═══════════════════════════════════════════════════════╗")
    println("║       MATCH-DAY DOPE CARD — 6.5 Creedmoor           ║")
    println("║  140 gr Berger Hybrid / H4350 / Lapua brass          ║")
    println("║  MV: 2700 fps  SD: 6 fps  Zero: 100 yd              ║")
    println("╠═══════════════════════════════════════════════════════╣")

    # Match conditions: 75°F, 29.65 inHg, 45% RH, 3500 ft altitude
    match_params = ShotParameters(
        mass_grains=140.0, caliber_in=0.264, bc=0.311, drag_model=:G7,
        muzzle_vel_fps=2700.0, sight_height_in=1.5, zero_range_yd=100.0,
        temp_F=75.0, pressure_inhg=29.65, humidity_pct=45.0, altitude_ft=3500.0,
        wind_speed_mph=0.0,  # calm — we'll add wind separately
        target_range_yd=1200.0,
        twist_in=8.0, twist_direction=1, bullet_length_in=1.34,
        latitude_deg=39.0,   # Colorado
        azimuth_deg=180.0,   # firing South
    )

    match_traj = solve_trajectory(match_params)
    sg_match = miller_sg(match_params)

    println("║  Conditions: 75°F, 29.65 inHg, 45%RH, 3500 ft      ║")
    println("║  Lat: 39°N, Azimuth: 180° (South)                   ║")
    println("╠═════════╤═════════╤═════════╤════════╤═══════════════╣")
    println("║ Rng(yd) │Elev(MOA)│Wind/mph │ ToF(s) │ Rem. Vel(fps)║")
    println("╠═════════╪═════════╪═════════╪════════╪═══════════════╣")

    # Also compute wind deflection for 1 mph crosswind at each range
    wind1_params = ShotParameters(match_params;
        wind_speed_mph=1.0, wind_angle_deg=90.0)
    wind_traj = solve_trajectory(wind1_params)

    step_m = yards_to_m(100.0)
    next_r = step_m
    for i in eachindex(match_traj)
        pt = match_traj[i]
        if pt.range_m >= next_r
            r_yd = m_to_yards(pt.range_m)
            drop_in = m_to_inches(pt.drop_m)
            elev = drop_to_moa(abs(drop_in), r_yd)
            if pt.drop_m < 0; elev = elev; else elev = -elev; end
            v_fps = ms_to_fps(pt.v_total)

            # Find matching wind point
            wpt = wind_traj[findfirst(p -> p.range_m >= pt.range_m, wind_traj)]
            wind_per_mph = m_to_inches(wpt.windage_m)

            @printf("║  %5.0f  │ %+6.1f  │ %5.1f   │ %5.3f  │    %7.0f    ║\n",
                    r_yd, elev, wind_per_mph, pt.time, v_fps)
            next_r += step_m
            r_yd > 1200 && break
        end
    end
    println("╚═════════╧═════════╧═════════╧════════╧═══════════════╝")
    println("  Wind column = inches per 1 mph full-value crosswind")
    println("  Multiply by actual wind to get total deflection")
end

# ╔═╡ 00000060-0000-0000-0000-000000000001
md"""
---
## Summary

This notebook demonstrated the full `CompetitionBallistics.jl` library:

| Module | What it does | Key functions |
|--------|-------------|---------------|
| `BallisticUtils` | Unit conversions, SD, stability | `miller_stability`, `sectional_density` |
| `Atmosphere` | ICAO model, humidity, DA | `air_density`, `density_altitude` |
| `InteriorBallistics` | Barrel physics | `leduc_velocity`, `barrel_length_correction` |
| `DragModels` | G1/G7 BRL tables | `cd_g1`, `cd_g7`, `drag_coefficient` |
| `ExteriorBallistics` | 3-DOF trajectory solver | `solve_trajectory`, `trajectory_table` |
| `ReloadingAnalysis` | Load development tools | `chronograph_stats`, `ladder_test_analysis` |
| `SixDOF` | Full rigid-body dynamics | `solve_6dof`, `gyroscopic_stability_6dof` |

The complete mathematical derivations are in the companion PDF manual.
"""

# ╔═╡ Cell order:
# ╟─00000001-0000-0000-0000-000000000001
# ╟─00000002-0000-0000-0000-000000000001
# ╠═00000003-0000-0000-0000-000000000001
# ╟─00000004-0000-0000-0000-000000000001
# ╟─00000005-0000-0000-0000-000000000001
# ╠═00000006-0000-0000-0000-000000000001
# ╟─00000007-0000-0000-0000-000000000001
# ╠═00000008-0000-0000-0000-000000000001
# ╟─00000009-0000-0000-0000-000000000001
# ╠═00000010-0000-0000-0000-000000000001
# ╟─00000011-0000-0000-0000-000000000001
# ╠═00000012-0000-0000-0000-000000000001
# ╟─00000013-0000-0000-0000-000000000001
# ╟─00000014-0000-0000-0000-000000000001
# ╠═00000015-0000-0000-0000-000000000001
# ╟─00000016-0000-0000-0000-000000000001
# ╠═00000017-0000-0000-0000-000000000001
# ╟─00000018-0000-0000-0000-000000000001
# ╠═00000019-0000-0000-0000-000000000001
# ╟─00000020-0000-0000-0000-000000000001
# ╟─00000021-0000-0000-0000-000000000001
# ╠═00000022-0000-0000-0000-000000000001
# ╟─00000023-0000-0000-0000-000000000001
# ╠═00000024-0000-0000-0000-000000000001
# ╟─00000025-0000-0000-0000-000000000001
# ╟─00000026-0000-0000-0000-000000000001
# ╠═00000027-0000-0000-0000-000000000001
# ╟─00000028-0000-0000-0000-000000000001
# ╠═00000029-0000-0000-0000-000000000001
# ╟─00000030-0000-0000-0000-000000000001
# ╠═00000031-0000-0000-0000-000000000001
# ╟─00000032-0000-0000-0000-000000000001
# ╠═00000033-0000-0000-0000-000000000001
# ╟─00000034-0000-0000-0000-000000000001
# ╠═00000035-0000-0000-0000-000000000001
# ╟─00000036-0000-0000-0000-000000000001
# ╟─00000037-0000-0000-0000-000000000001
# ╠═00000038-0000-0000-0000-000000000001
# ╟─00000039-0000-0000-0000-000000000001
# ╠═00000040-0000-0000-0000-000000000001
# ╟─00000041-0000-0000-0000-000000000001
# ╠═00000042-0000-0000-0000-000000000001
# ╟─00000043-0000-0000-0000-000000000001
# ╠═00000044-0000-0000-0000-000000000001
# ╟─00000045-0000-0000-0000-000000000001
# ╠═00000046-0000-0000-0000-000000000001
# ╟─00000047-0000-0000-0000-000000000001
# ╠═00000048-0000-0000-0000-000000000001
# ╟─00000049-0000-0000-0000-000000000001
# ╟─00000050-0000-0000-0000-000000000001
# ╠═00000051-0000-0000-0000-000000000001
# ╟─00000052-0000-0000-0000-000000000001
# ╠═00000053-0000-0000-0000-000000000001
# ╟─00000054-0000-0000-0000-000000000001
# ╠═00000055-0000-0000-0000-000000000001
# ╟─00000056-0000-0000-0000-000000000001
# ╠═00000057-0000-0000-0000-000000000001
# ╟─00000058-0000-0000-0000-000000000001
# ╠═00000059-0000-0000-0000-000000000001
# ╟─00000060-0000-0000-0000-000000000001
