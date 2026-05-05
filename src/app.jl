# ============================================================================
# Competition Ballistics Web App
# Built with Dash.jl
# ============================================================================

using Dash
using PlotlyJS
using Printf

# Load the library
include("CompetitionBallistics.jl")
using .CompetitionBallistics
using .CompetitionBallistics.ReferenceData
using .CompetitionBallistics.ExteriorBallistics
using .CompetitionBallistics.BallisticUtils
using .CompetitionBallistics.Atmosphere

app = dash(external_stylesheets=["https://codepen.io/chriddyp/pen/bWLwgP.css"])

# Prepare bullet options for dropdown
bullet_options = [Dict("label" => k, "value" => k) for k in sort(collect(keys(COMMON_BULLETS)))]

app.layout = html_div() do
    html_h1("Competition Ballistics Calculator", style=Dict("textAlign" => "center")),
    
    html_div(style=Dict("display" => "flex", "flex-direction" => "row")) do
        # Sidebar for inputs
        html_div(style=Dict("width" => "30%", "padding" => "20px", "background-color" => "#f9f9f9", "border-radius" => "10px")) do
            html_h3("Input Parameters"),
            
            html_label("Select Bullet Profile:"),
            dcc_dropdown(id="bullet-dropdown", options=bullet_options, placeholder="Select a bullet..."),
            
            html_hr(),
            
            html_label("Bullet Mass (grains):"),
            dcc_input(id="mass", type="number", value=140.0, step=0.1),
            
            html_label("Caliber (inches):"),
            dcc_input(id="caliber", type="number", value=0.264, step=0.001),
            
            html_label("Ballistic Coefficient (G7):"),
            dcc_input(id="bc", type="number", value=0.311, step=0.001),
            
            html_label("Muzzle Velocity (fps):"),
            dcc_input(id="mv", type="number", value=2700.0, step=1),
            
            html_label("Sight Height (inches):"),
            dcc_input(id="sh", type="number", value=1.5, step=0.1),
            
            html_label("Zero Range (yards):"),
            dcc_input(id="zero", type="number", value=100.0, step=10),
            
            html_hr(),
            
            html_h4("Environment"),
            html_label("Temperature (°F):"),
            dcc_input(id="temp", type="number", value=59.0, step=1),
            
            html_label("Pressure (inHg):"),
            dcc_input(id="pressure", type="number", value=29.92, step=0.01),
            
            html_label("Wind Speed (mph):"),
            dcc_input(id="wind-speed", type="number", value=10.0, step=1),
            
            html_label("Wind Angle (deg):"),
            dcc_input(id="wind-angle", type="number", value=90.0, step=5),
            
            html_hr(),
            
            html_label("Max Range (yards):"),
            dcc_input(id="max-range", type="number", value=1000.0, step=100)
        end,
        
        # Main panel for plots and table
        html_div(style=Dict("width" => "70%", "padding" => "20px")) do
            dcc_tabs(id="tabs", value="tab-plots") do
                dcc_tab(label="Trajectory Plots", value="tab-plots") do
                    html_div() do
                        dcc_graph(id="drop-graph"),
                        dcc_graph(id="windage-graph")
                    end
                end,
                dcc_tab(label="Trajectory Table", value="tab-table") do
                    dash_datatable(
                        id="traj-table",
                        columns=[
                            Dict("name" => "Range (yd)", "id" => "range"),
                            Dict("name" => "Drop (in)", "id" => "drop"),
                            Dict("name" => "Elev (MOA)", "id" => "elev"),
                            Dict("name" => "Wind (in)", "id" => "wind"),
                            Dict("name" => "Velocity (fps)", "id" => "vel"),
                            Dict("name" => "Energy (ft-lb)", "id" => "energy")
                        ],
                        style_table=Dict("overflowX" => "auto"),
                        style_cell=Dict("textAlign" => "center")
                    )
                end
            end
        end
    end
end

# Callback to update inputs from dropdown
callback!(app,
    Output("mass", "value"),
    Output("caliber", "value"),
    Output("bc", "value"),
    Input("bullet-dropdown", "value")
) do bullet_name
    if isnothing(bullet_name) || !haskey(COMMON_BULLETS, bullet_name)
        return 140.0, 0.264, 0.311
    end
    profile = COMMON_BULLETS[bullet_name]
    return profile.mass_gr, profile.caliber_in, profile.bc_g7
end

# Callback to update graphs and table
callback!(app,
    Output("drop-graph", "figure"),
    Output("windage-graph", "figure"),
    Output("traj-table", "data"),
    Input("mass", "value"),
    Input("caliber", "value"),
    Input("bc", "value"),
    Input("mv", "value"),
    Input("sh", "value"),
    Input("zero", "value"),
    Input("temp", "value"),
    Input("pressure", "value"),
    Input("wind-speed", "value"),
    Input("wind-angle", "value"),
    Input("max-range", "value")
) do mass, caliber, bc, mv, sh, zero, temp, pressure, wind_speed, wind_angle, max_range
    # Ensure all inputs are numeric
    mass = Float64(mass)
    caliber = Float64(caliber)
    bc = Float64(bc)
    mv = Float64(mv)
    sh = Float64(sh)
    zero = Float64(zero)
    temp = Float64(temp)
    pressure = Float64(pressure)
    wind_speed = Float64(wind_speed)
    wind_angle = Float64(wind_angle)
    max_range = Float64(max_range)

    params = ShotParameters(
        mass_grains = mass,
        caliber_in = caliber,
        bc = bc,
        drag_model = :G7,
        muzzle_vel_fps = mv,
        sight_height_in = sh,
        zero_range_yd = zero,
        temp_F = temp,
        pressure_inhg = pressure,
        wind_speed_mph = wind_speed,
        wind_angle_deg = wind_angle,
        target_range_yd = max_range
    )

    traj = solve_trajectory(params)
    
    # Filter trajectory for the table (e.g., every 50 or 100 yards)
    step_yd = 50.0
    table_data = []
    next_r_m = 0.0
    for pt in traj
        if pt.range_m >= next_r_m
            r_yd = m_to_yards(pt.range_m)
            d_in = m_to_inches(pt.drop_m)
            w_in = m_to_inches(pt.windage_m)
            elev = r_yd > 0 ? drop_to_moa(abs(d_in), r_yd) : 0.0
            if pt.drop_m > 0; elev = -elev; end
            
            push!(table_data, Dict(
                "range" => round(r_yd, digits=0),
                "drop" => round(d_in, digits=1),
                "elev" => round(elev, digits=1),
                "wind" => round(w_in, digits=1),
                "vel" => round(ms_to_fps(pt.v_total), digits=0),
                "energy" => round(pt.energy_J / 1.35582, digits=0)
            ))
            next_r_m += yards_to_m(step_yd)
        end
    end

    ranges = [m_to_yards(pt.range_m) for pt in traj]
    drops = [m_to_inches(pt.drop_m) for pt in traj]
    windages = [m_to_inches(pt.windage_m) for pt in traj]

    fig_drop = Plot(
        ranges, drops,
        Layout(
            title="Bullet Drop",
            xaxis_title="Range (yards)",
            yaxis_title="Drop (inches)",
            hovermode="x"
        )
    )

    fig_wind = Plot(
        ranges, windages,
        Layout(
            title="Wind Deflection",
            xaxis_title="Range (yards)",
            yaxis_title="Deflection (inches)",
            hovermode="x"
        )
    )

    return fig_drop, fig_wind, table_data
end

run_server(app, "0.0.0.0", 8050, debug=true)
