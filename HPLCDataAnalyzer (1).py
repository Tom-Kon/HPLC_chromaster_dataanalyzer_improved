import numpy as np
import pandas as pd
from dash import dash
from dash import Dash, dcc, html, Input, Output, State, no_update
from dash import dash_table
import plotly.graph_objects as go
from scipy.integrate import simpson
import base64
from dash.exceptions import PreventUpdate
from dash import callback_context
from sklearn.linear_model import LinearRegression
import dash_bootstrap_components as dbc
from dash import html
import io
import numpy as np
import statsmodels.api as sm
from scipy.stats import t
import plotly.io as pio
import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend suitable for servers
import matplotlib.pyplot as plt
import re
import webbrowser
import threading
import time
from datetime import datetime, timedelta
import struct
import os

last_activity_time = datetime.now()
EXIT_TIMEOUT_MINUTES = 30

def inactivity_monitor():
    global last_activity_time
    while True:
        time.sleep(60)  # Check every 60 seconds
        if datetime.now() - last_activity_time > timedelta(minutes=EXIT_TIMEOUT_MINUTES):
            print("No activity detected. Exiting app.")
            os._exit(0)  # Use os._exit to ensure it works even when frozen

threading.Thread(target=inactivity_monitor, daemon=True).start()


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.CERULEAN])
app.title = "HPLC Data Analyzer"

# ---------- Function to read RW1 data ----------
def read_rw1_24bit_from_bytes(file_bytes):
    boundary_index = file_bytes.find(b"Boundary")
    if boundary_index == -1:
        raise ValueError("Boundary marker not found in file.")
    start_index = boundary_index + len(b"Boundary") + 19
    raw_data = file_bytes[start_index:]

    signal = []
    for i in range(0, len(raw_data) - 2, 3):
        b = raw_data[i:i+3]
        val = (b[0] << 16) | (b[1] << 8) | b[2]
        if val & 0x800000:
            val -= 1 << 24
        signal.append(val)

    signal = np.array(signal)
    signalNorm = signal
    sample_interval = 0.2
    time = np.arange(0, len(signal) * sample_interval, sample_interval)
    time_minutes = time / 60

    df = pd.DataFrame({'Time': time_minutes, 'Signal': signalNorm})
    return df

# ---------- App Layout with Tabs ----------
app.layout = dbc.Container([
    dbc.Row(
        dbc.Col(html.H1("HPLC Data Analyzer"), className="text-center my-4")
    ),

    dcc.Tabs([
        dcc.Tab(label='Data upload', children=[
            html.Div([
                dbc.Row(
                    dbc.Col(
                        dcc.Upload(
                            id='upload-data',
                            children=html.Div(['Drag and Drop or ', html.A('Select one or more .RW1 Files for analysis')]),
                            style={
                                'width': '60%', 'height': '60px', 'lineHeight': '60px',
                                'borderWidth': '2px', 'borderStyle': 'dashed', 'borderRadius': '5px',
                                'textAlign': 'center', 'margin': '20px auto'
                            },
                            multiple=True
                        ), width={"size": 8, "offset": 2}
                    )
                ),
            ], style={'paddingTop': '50px'})
        ]),

        dcc.Tab(label='Chromatogram Viewer', children=[
            html.Div([ 
                # Stores here if needed globally (same as before)
                dcc.Store(id='stored-data-list', data=[]),
                dcc.Store(id='current-plot-index', data=0),
                dcc.Store(id='rw1-names'),
                dcc.Store(id='selected-region'),
                dcc.Store(id='auc-storage', data=[]),
                dcc.Store(id='selection-store'),

                dcc.Graph(id='chromatogram-plot', config={"scrollZoom": True}),

                dbc.Row([
                    dbc.Col(dbc.Button("Previous Plot", id="prev-plot", n_clicks=0, color="primary"), width="auto"),
                    dbc.Col(dbc.Button("Next Plot", id="next-plot", n_clicks=0, color="primary"), width="auto"),
                    dbc.Col(
                        dbc.Button("Add to Table", id="add-to-table", n_clicks=0, color="success"), 
                        width="auto"
                    ),
                    dbc.Col(
                        dbc.Checklist(
                            options=[{'label': 'This is calibration curve data', 'value': 'enabled'}],
                            value=[],
                            id='calibration-data-checkbox',
                            switch=True,
                            inline=True,
                        ),
                        width="auto", className="d-flex align-items-center"
                    ),
                    dbc.Col(
                        dbc.Checklist(
                            options=[{'label': 'This is the blank', 'value': 'enabled'}],
                            value=[],
                            id='blank-checkbox',
                            switch=True,
                            inline=True,
                        ),
                        width="auto", className="d-flex align-items-center"
                    ),
                    dbc.Col(dbc.Button("Export All Chromatograms as an Excel", id="export-all-data-excel", n_clicks=0, color="info"), width="auto"),
                    dcc.Download(id="download-all-data-excel")
                ], justify="center", className="my-3", style={"gap": "10px"}),

                html.Div(id='checkbox-output'),
            ], style={'paddingBottom': '150px', 'paddingTop': '50px'})  # <-- padding added here
        ]),

        dcc.Tab(label='Calibration data analysis', children=[
            html.Div([
                    html.H3("Blank sample"),
                    dash.dash_table.DataTable(
                        id='blank-sample',
                        columns=[
                            {"name": "Area", "id": "area", "editable": False},
                            {"name": "Retention Time", "id": "retention_time", "editable": False},
                        ],
                        data=[],
                        style_table={
                            'overflowX': 'auto',
                            'width': '21%',
                            'marginLeft': '0',        # Aligns left
                            'marginRight': 'auto',    # Prevents horizontal centering
                            'marginBottom': '40px',   # Adds space below for separation
                        },
                        style_cell={
                            'textAlign': 'center',
                            'padding': '8px',
                            'fontFamily': 'Helvetica, Arial, sans-serif',
                            'fontSize': '14px',
                        },
                        style_header={
                            'backgroundColor': '#5e9fc2',  # Bootstrap Cerulean primary
                            'color': 'white',
                            'fontWeight': 'bold',
                            'border': '1px solid #ccc',
                        },
                        style_data={
                            'border': '1px solid #ccc',
                            'backgroundColor': '#f9f9f9',
                        },
                        style_data_conditional=[
                            {
                                'if': {'row_index': 'odd'},
                                'backgroundColor': '#d5edfa'  # Light Cerulean-ish tint for odd rows
                            }
                        ],
                        row_deletable=True,
                    ),
                    html.H3("Calibration"),
                    dash.dash_table.DataTable(
                        id='auc-calibration-regression',
                        columns=[
                            {"name": "Sample name", "id": "sample_name", "editable": True},
                            {"name": "Retention time (min)", "id": "retention_time", "editable": False},
                            {"name": "Concentration (µg/mL)", "id": "concentration", "editable": True},
                            {"name": "Area (AU·min)", "id": "area", "editable": False},
                        ],
                        data=[],
                        style_table={'overflowX': 'auto'},
                        style_cell={
                            'textAlign': 'center',
                            'padding': '8px',
                            'fontFamily': 'Helvetica, Arial, sans-serif',
                            'fontSize': '14px',
                        },
                        style_header={
                            'backgroundColor': '#5e9fc2',  # Bootstrap Cerulean primary
                            'color': 'white',
                            'fontWeight': 'bold',
                            'border': '1px solid #ccc',
                        },
                        style_data={
                            'border': '1px solid #ccc',
                            'backgroundColor': '#f9f9f9',
                        },
                        style_data_conditional=[
                            {
                                'if': {'row_index': 'odd'},
                                'backgroundColor': '#d5edfa'  # Light Cerulean-ish tint for odd rows
                            }
                        ],
                        row_deletable=True,
                    ),
                    dbc.Row([
                        dbc.Col([
                            dbc.Checklist(
                                options=[{'label': 'Autofill concentrations', 'value': 'enabled'}],
                                value=[],
                                id='concentration-autofill',
                                switch=True,
                                inline=True,
                            ),
                        ], width="auto", className="d-flex align-items-center"),
                        dbc.Col([
                            dbc.Button("Perform regression", id="perform-regression", n_clicks=0, color="primary", className="my-2"),
                        ], width="auto"),
                    ], className="my-4", style={"gap": "10px"}),
                    dcc.Graph(id="regression-plot", style={"marginTop": "20px"}),

                    html.H3("Regression Summary"),
                    dash_table.DataTable(
                        id='regression-summary-table',
                        columns=[
                            {"name": "Parameter", "id": "parameter"},
                            {"name": "Value", "id": "value", "type": "numeric", "format": {"specifier": ".6f"}},
                            {"name": "95% CI Lower", "id": "CI Lower (95%)", "type": "numeric", "format": {"specifier": ".6f"}},
                            {"name": "95% CI Upper", "id": "CI Upper (95%)", "type": "numeric", "format": {"specifier": ".6f"}},
                        ],
                        data=[],  # initially empty
                        style_table={'overflowX': 'auto'},
                        style_cell={
                            'textAlign': 'center',
                            'padding': '8px',
                            'fontFamily': 'Helvetica, Arial, sans-serif',
                            'fontSize': '14px',
                        },
                        style_header={
                            'backgroundColor': "#5e9fc2",  # Bootstrap Cerulean primary
                            'color': 'white',
                            'fontWeight': 'bold',
                            'border': '1px solid #ccc',
                        },
                        style_data={
                            'border': '1px solid #ccc',
                            'backgroundColor': '#f9f9f9',
                        },
                        style_data_conditional=[
                            {
                                'if': {'row_index': 'odd'},
                                'backgroundColor': '#d5edfa'  # Light Cerulean-ish tint for odd rows
                            }
                        ],
                    ),
                ], style={'paddingTop': '50px', 'paddingBottom': '150px'})
            ]),

            dcc.Tab(label='Sample data analysis', children=[
                html.Div([
                    dash.dash_table.DataTable(
                        id='final-sample-table',
                        columns=[
                            {"name": "Sample name", "id": "sample_name", "editable": True},
                            {"name": "Retention time (min)", "id": "retention_time", "editable": False},
                            {"name": "Area (AU·min)", "id": "area","editable": False},
                            {"name": "Concentration (µg/mL)", "id": "concentration", "editable": False},
                        ],
                        data=[],
                        style_table={'overflowX': 'auto'},
                        style_cell={
                            'textAlign': 'center',
                            'padding': '8px',
                            'fontFamily': 'Helvetica, Arial, sans-serif',
                            'fontSize': '14px',
                        },
                        style_header={
                            'backgroundColor': "#5e9fc2",  # Bootstrap Cerulean primary
                            'color': 'white',
                            'fontWeight': 'bold',
                            'border': '1px solid #ccc',
                        },
                        style_data={
                            'border': '1px solid #ccc',
                            'backgroundColor': '#f9f9f9',
                        },
                        style_data_conditional=[
                            {
                                'if': {'row_index': 'odd'},
                                'backgroundColor': '#d5edfa'  # Light Cerulean-ish tint for odd rows
                            }
                        ],
                        row_deletable=True,
                    ),
                    dbc.Row([
                        dbc.Col(dbc.Button("Perform Concentration determination", id="perform-sample-conc-det", n_clicks=0, color="primary"), width="auto"),
                        dbc.Col(dbc.Button("Export to Excel", id="export-excel-analyzed-sample-table", n_clicks=0, color="info"), width="auto"),
                        dbc.Col([
                            dbc.Checklist(
                                options=[{'label': 'I performed a dilution', 'value': 'enabled'}],
                                value=[],
                                id='dilution-checkbox',
                                switch=True,
                                inline=True,
                            ),
                        ], width="auto", style={"marginTop": "10px"}),
                        dbc.Col([
                            dcc.Input(
                                id='dilution-factor',
                                type='number',
                                placeholder='Only numbers allowed',
                                min=0,  # optional
                                step=1,  # optional
                                style={'width': '50%'}
                            )
                        ], width="auto"),    
                        dcc.Download(id="download-excel-analyzed-sample-table"),
                    ], className="my-4"),
                    html.Div([
                        dash.dash_table.DataTable(
                                    id='dilution-corrected-sample-table',
                                    columns=[
                                        {"name": "Sample name", "id": "sample_name", "editable": True},
                                        {"name": "Dilution factor", "id": "dilution_factor", "editable": True},
                                        {"name": "Concentration before dilution (µg/mL)", "id": "concentration_bef_dil", "editable": False},
                                    ],
                                    data=[],
                                    style_table={'overflowX': 'auto'},
                                    style_cell={
                                        'textAlign': 'center',
                                        'padding': '8px',
                                        'fontFamily': 'Helvetica, Arial, sans-serif',
                                        'fontSize': '14px',
                                    },
                                    style_header={
                                        'backgroundColor': "#5e9fc2",  # Bootstrap Cerulean primary
                                        'color': 'white',
                                        'fontWeight': 'bold',
                                        'border': '1px solid #ccc',
                                    },
                                    style_data={
                                        'border': '1px solid #ccc',
                                        'backgroundColor': '#f9f9f9',
                                    },
                                    style_data_conditional=[
                                        {
                                            'if': {'row_index': 'odd'},
                                            'backgroundColor': '#d5edfa'  # Light Cerulean-ish tint for odd rows
                                        }
                                    ],
                                    row_deletable=True,
                                ),                    
                    ], id='dilution-corrected-wrapper', style={'display': 'none'})
                ], style={'paddingTop': '50px'}),
            ]),      # ✅ closes dcc.Tab for Sample data analysis

            dcc.Tab(label='Pressure as a function of time', children=[
                html.Div([
                        dcc.Store(id='stored-bytes-list', data=[]),
                        dcc.Store(id='current-plot-index-pressure', data=0),


                        dcc.Graph(id='pressure-plot', config={"scrollZoom": False}),
                        dbc.Row([
                            dbc.Col(dbc.Button("Previous Plot", id="prev-plot-pressure", n_clicks=0, color="primary"), width="auto"),
                            dbc.Col(dbc.Button("Next Plot", id="next-plot-pressure", n_clicks=0, color="primary"), width="auto"),
                            
                        ], justify="center", className="my-3", style={"gap": "10px", 'paddingBottom': '50px'}), 
                        dcc.Graph(id='pressure-summary-plot', config={"scrollZoom": False}),
                ], style={'paddingTop': '50px', 'paddingBottom': '150px'})
            ])
    ]), 

    dcc.Interval(id="keepalive", interval=10*1000, n_intervals=0),  # every 10 seconds

], fluid=True)


# ---------- Fix: Chromatogram Upload + Navigation ----------
@app.callback(
    Output('stored-data-list', 'data'),
    Output('rw1-names', 'data'),
    Output('current-plot-index', 'data'),
    Output('stored-bytes-list', 'data'),
    Input('upload-data', 'contents'),
    Input('upload-data', 'filename'),
    Input('prev-plot', 'n_clicks'),
    Input('next-plot', 'n_clicks'),
    State('stored-data-list', 'data'),
    State('rw1-names', 'data'),
    State('current-plot-index', 'data'),
    prevent_initial_call=True
)
def upload_and_navigate(contents, filenames, prev_clicks, next_clicks, stored_list, names_list, current_index):
    ctx = callback_context
    if not ctx.triggered:
        raise PreventUpdate

    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]

    def extract_sample_name_from_rw1_bytes_v2(data: bytes):
        key = b'VialName'
        idx = data.find(key)
        if idx == -1:
            return None
        name_start = idx + len(key) + 7
        end_idx = data.find(b'\x00', name_start)
        if end_idx == -1:
            return None
        name_bytes = data[name_start:end_idx]
        try:
            sample_name = name_bytes.decode('latin-1').strip()
        except Exception:
            return None
        return sample_name

    def extract_numeric(filename):
        match = re.search(r'(\d+)', filename)
        return int(match.group(1)) if match else float('inf')

    if triggered_id == 'upload-data':
        if not contents or not filenames:
            raise PreventUpdate

        sorted_files = sorted(zip(contents, filenames), key=lambda x: extract_numeric(x[1]))
        data_list, extracted_names, raw_bytes_list = [], [], []

        for content, filename in sorted_files:
            content_type, content_string = content.split(',')
            decoded = base64.b64decode(content_string)

            try:
                df = read_rw1_24bit_from_bytes(decoded)
                data_list.append(df.to_json(date_format='iso', orient='split'))
                raw_bytes_list.append(content_string)
            except Exception as e:
                print(f"Upload error parsing RW1 data: {e}")
                continue

            name = extract_sample_name_from_rw1_bytes_v2(decoded) or "Unknown Sample"
            extracted_names.append(name)

        return data_list, extracted_names, 0, raw_bytes_list

    if not stored_list:
        raise PreventUpdate

    max_index = len(stored_list) - 1
    if current_index is None:
        current_index = 0

    if triggered_id == 'prev-plot':
        new_index = max(0, current_index - 1)
    elif triggered_id == 'next-plot':
        new_index = min(max_index, current_index + 1)
    else:
        new_index = current_index

    return stored_list, names_list or [], new_index, dash.no_update


# ---------- NEW Callback: Pressure Navigation Index ----------
@app.callback(
    Output('current-plot-index-pressure', 'data'),
    Input('prev-plot-pressure', 'n_clicks'),
    Input('next-plot-pressure', 'n_clicks'),
    State('stored-bytes-list', 'data'),
    State('current-plot-index-pressure', 'data'),
    prevent_initial_call=True
)
def navigate_pressure_plot(prev_clicks, next_clicks, byte_list, current_index):
    ctx = callback_context
    if not ctx.triggered or not byte_list:
        raise PreventUpdate

    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
    max_index = len(byte_list) - 1
    if current_index is None:
        current_index = 0

    if triggered_id == 'prev-plot-pressure':
        new_index = max(0, current_index - 1)
    elif triggered_id == 'next-plot-pressure':
        new_index = min(max_index, current_index + 1)
    else:
        new_index = current_index

    return new_index

# ---------- Callback: Pressure Summary Plot ----------
@app.callback(
    Output('pressure-summary-plot', 'figure'),
    Input('stored-bytes-list', 'data'),
    State('rw1-names', 'data'),
    prevent_initial_call=True
)
def plot_pressure_summary(byte_data_list, names):
    if not byte_data_list:
        raise PreventUpdate

    avg_pressures = []
    std_pressures = []
    hover_labels = []


    for i, content_string in enumerate(byte_data_list):
        try:
            decoded = base64.b64decode(content_string)
            pressure_values, _ = extract_pressure_profile_16bit(decoded)
            avg = np.mean(pressure_values)
            std = np.std(pressure_values)
            avg_pressures.append(avg)
            std_pressures.append(std)
            name = names[i] if names and i < len(names) else f"Sample {i}"
            hover_labels.append(f"{name}<br>Avg: {avg:.2f} bar<br>Std: {std:.2f} bar")
        except Exception as e:
            print(f"Error processing pressure profile: {e}")
            avg_pressures.append(None)
            std_pressures.append(None)
            hover_labels.append("Error")

    indices = list(range(len(avg_pressures)))

    # Use actual names if provided
    if names and len(names) == len(avg_pressures):
        hover_text = [f"Sample: {n}" for n in names]
    else:
        hover_text = [f"Sample: {i}" for i in indices]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=indices,
        y=avg_pressures,
        error_y=dict(
            type='data',
            array=std_pressures,
            visible=True
        ),
        mode='markers+lines',
        name='Avg Pressure',
        marker=dict(color='darkblue', size=10),
        text=hover_labels,
        hovertemplate='%{text}<extra></extra>'
    ))

    fig.update_layout(
        title=dict(
            text="Average Pressure per Sample",
            font=dict(size=26)  # Correct placement of size
        ), 
        xaxis=dict(
            title="Sample index",
            title_font=dict(size=22),
            tickfont=dict(size=18),
            tickmode = 'linear',
            tick0 = 0,
            dtick=1,           # major ticks every 1 unit (1 min)
            minor=dict(
                ticklen=6,     # length of minor ticks (smaller)
                tickcolor="black",
                dtick=0.5      # minor ticks every 0.2 units (0.2 min)
            ),
            ticks="outside",   # major ticks outside
            minor_ticks="outside", # minor ticks outside too
            linecolor="black",  # Force axis line color
            linewidth=1,        # Make it thicker
        ),
        yaxis=dict(
            title="Average Pressure ± SD (bar)",
            title_font=dict(size=22),
            tickfont=dict(size=18), 
            linecolor="black",  # Force axis line color
            linewidth=1,        # Make it thicker
        ),
        margin=dict(t=50, b=40, l=60, r=20),
        font=dict(size=14),
        height=600,
        hovermode='x unified',
    )

    return fig





def extract_pressure_profile_16bit(data: bytes):
    key_start = b'PressureProfileA'
    key_end = b'PressureProfileB'

    idx_start = data.find(key_start)
    idx_end = data.find(key_end)

    if idx_start == -1 or idx_end == -1 or idx_end <= idx_start:
        print("PressureProfile markers not found or invalid order.")
        return None

    # Offset +5 to skip unknown header bytes (based on your earlier function)
    start = idx_start + len(key_start) + 7
    end = idx_end - 21
    raw_profile_data = data[start:end]
    sample_interval = 2
    time = np.arange(0, len(raw_profile_data) * sample_interval, sample_interval)
    time_minutes_pressure = time / 60  # Convert to minutes

    # Make sure the length is even to unpack 16-bit integers
    if len(raw_profile_data) % 2 != 0:
        raw_profile_data = raw_profile_data[:-1]  # Trim last byte if odd

    # Unpack as little-endian signed 16-bit integers
    pressure_values = struct.unpack("<" + "h" * (len(raw_profile_data) // 2), raw_profile_data)

    
    return pressure_values, time_minutes_pressure

@app.callback(
    Output('pressure-plot', 'figure'),
    Input('stored-bytes-list', 'data'),
    Input('current-plot-index-pressure', 'data'),
    State('rw1-names', 'data'),
    prevent_initial_call=True
)
def update_pressure_plot(byte_data_list, current_index, current_name):
    if not byte_data_list or current_index is None or current_index >= len(byte_data_list):
        return go.Figure(layout=dict(title="No pressure profile to display"))

    try:
        decoded = base64.b64decode(byte_data_list[current_index])
        data = extract_pressure_profile_16bit(decoded)
        pressure = data[0]
        time_minutes_pressure = data[1]

    except Exception as e:
        print(f"Error decoding pressure profile: {e}")
        return go.Figure(layout=dict(title="Error reading pressure profile"))

    if not pressure:
        return go.Figure(layout=dict(title="No pressure data found"))

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        y=pressure,
        x = time_minutes_pressure,
        mode='lines',
        line=dict(color='blue'),
        name='Pressure',
    ))

    fig.update_layout(
        title=dict(
            text=f"Pressure as a function of time - {current_name[current_index]}",
            font=dict(size=26)  # Title font size
        ),
        xaxis=dict(
            title="Time (minutes)",
            title_font=dict(size=22),
            tickfont=dict(size=18),
            dtick=1,           # major ticks every 1 unit (1 min)
            minor=dict(
                ticklen=6,     # length of minor ticks (smaller)
                tickcolor="black",
                dtick=0.5      # minor ticks every 0.2 units (0.2 min)
            ),
            ticks="outside",   # major ticks outside
            minor_ticks="outside", # minor ticks outside too
            linecolor="black",  # Force axis line color
            linewidth=1,        # Make it thicker
        ),
        yaxis=dict(
            title="Pressure (bar)",
            title_font=dict(size=22),
            tickfont=dict(size=18), 
            linecolor="black",  # Force axis line color
            linewidth=1,        # Make it thicker
        ),
        height=600,
        margin=dict(t=50, b=40, l=60, r=20),
        font=dict(size=14),
    )

    return fig


# ---------- Callback: Store selectedData separately ----------
@app.callback(
    Output('selection-store', 'data'),
    Input('chromatogram-plot', 'selectedData'),
    prevent_initial_call=True
)
def store_selection(selectedData):
    if selectedData and 'range' in selectedData and 'x' in selectedData['range']:
        return selectedData
    return no_update

# ---------- Callback: Create plot ----------
@app.callback(
    Output('chromatogram-plot', 'figure'),
    Input('stored-data-list', 'data'),
    Input('current-plot-index', 'data'),
    Input('selection-store', 'data'),
    State('rw1-names', 'data'),  # <- Add this
    State('blank-sample', 'data'),  # <- Add this
    prevent_initial_call=True
)

def create_plot_multi(data_list, current_index, selectedData, name_list, blank_sample):
    if not data_list or current_index is None or current_index >= len(data_list):
        return go.Figure(layout=dict(title="No data to display"))

    df = pd.read_json(io.StringIO(data_list[current_index]), orient='split')

    fig = go.Figure()
    # Full chromatogram line
    
    fig.add_trace(go.Scatter(
        x=df["Time"],
        y=df["Signal"],
        mode='lines',
        hovertemplate='Time: %{x:.2f} min<br>Signal: %{y:.3f}',  # no <extra></extra>
        name=name_list[current_index] if name_list and current_index < len(name_list) else f"RW1 File {current_index+1}"
    ))
    auc = None  # initialize

    if selectedData and 'range' in selectedData and 'x' in selectedData['range']:
        x0, x1 = selectedData['range']['x']
        mask = (df["Time"] >= x0) & (df["Time"] <= x1)
        selected_time = df["Time"][mask]
        selected_signal = df["Signal"][mask]

        if len(selected_time) > 1:
            # Linear baseline between first and last points in selection
            y0 = selected_signal.iloc[0]
            y1 = selected_signal.iloc[-1]
            baseline = np.linspace(y0, y1, len(selected_signal))
            corrected_signal = selected_signal - baseline

            # Fill the area under the corrected curve (shifted by baseline to be aligned)
            fig.add_trace(go.Scatter(
                x=np.concatenate(([selected_time.iloc[0]], selected_time, [selected_time.iloc[-1]])),
                y=np.concatenate(([baseline[0]], selected_signal, [baseline[-1]])),
                fill='toself',
                fillcolor='rgba(255, 0, 0, 0.3)',  # semi-transparent red fill
                line=dict(color='rgba(255, 0, 0, 0)'),  # no line
                hoverinfo='skip',
                showlegend=False,
                name='Area under curve'
            ))

            # Calculate AUC here
            if blank_sample and 'area' in blank_sample[0] and 'retention_time' in blank_sample[0]:
                auc = simpson(corrected_signal, selected_time) - blank_sample[0]['area']

            else:
                auc = simpson(corrected_signal, selected_time)

    # Plot title and annotations as before
    plot_title = name_list[current_index] if name_list and current_index < len(name_list) else f"RW1 File {current_index+1}"

    fig.update_layout(
        title=dict(
            text=f"Chromatogram - {plot_title}",
            font=dict(size=26)  # Title font size
        ),
        xaxis=dict(
            title="Time (minutes)",
            title_font=dict(size=22),
            tickfont=dict(size=18),
            dtick=1,           # major ticks every 1 unit (1 min)
            minor=dict(
                ticklen=6,     # length of minor ticks (smaller)
                tickcolor="black",
                dtick=0.5      # minor ticks every 0.2 units (0.2 min)
            ),
            ticks="outside",   # major ticks outside
            minor_ticks="outside", # minor ticks outside too
            linecolor="black",  # Force axis line color
            linewidth=1,        # Make it thicker
        ),
        yaxis=dict(
            title="Signal (a.u.)",
            title_font=dict(size=22),
            tickfont=dict(size=18), 
            linecolor="black",  # Force axis line color
            linewidth=1,        # Make it thicker
        ),
        font=dict(
            size=16  # Default font size for axis labels, legend, hover labels, etc.
        ),
        dragmode="select",
        selectdirection="h",
        showlegend=False,
        height=650
    )

    # Add vertical lines at selection start and end
    if selectedData and 'range' in selectedData and 'x' in selectedData['range']:
        fig.add_vline(
            x=x0,
            line=dict(color='red', dash='dash'),
            annotation_text="Start",
            annotation_position="top right",
            line_width=1
        )
        fig.add_vline(
            x=x1,
            line=dict(color='red', dash='dash'),
            annotation_text="End",
            annotation_position="top left",
            line_width=1
        )

        # Add annotation if auc was computed
        if auc is not None:
            # Calculate retention time as time at max signal in selected region
            retention_idx = selected_signal.idxmax()
            retention_time = selected_time.loc[retention_idx]

            # Place annotation at the middle of the selected range, y at max signal in selection + some offset
            line_color = fig.data[0].line.color  # this gets the line color of the first trace


            y_pos = selected_signal.max() + 0.05 * (df["Signal"].max() - df["Signal"].min())
            fig.add_annotation(
                x=x0-0.02*x0,
                xanchor='right',  # This makes the annotation box's right edge align with x=x_pos
                y=y_pos-0.025*y_pos,
                text=f"AUC: {auc:.2f} AU·min<br>Retention Time: {retention_time:.2f} min",
                showarrow=False,
                bordercolor= line_color,
                font=dict(color=line_color)
            )

    return fig

# ---------- Callback: Export all data as Excel ----------

def sanitize_sheet_name(name):
    # Excel sheet name rules:
    # max 31 chars, no : \ / ? * [ ] 
    # Also remove leading/trailing whitespace
    name = str(name)
    name = re.sub(r'[:\\/?*\[\]]', '_', name)
    name = name.strip()
    if len(name) > 31:
        name = name[:31]
    if len(name) == 0:
        name = "Sheet"
    return name


@app.callback(
    Output("download-all-data-excel", "data"),
    Input("export-all-data-excel", "n_clicks"),
    State("stored-data-list", "data"),
    State("rw1-names", "data"),
    prevent_initial_call=True
)
def export_all_data_excel(n_clicks, stored_data_list, name_list):
    if not stored_data_list:
        raise PreventUpdate
    
    buffer = io.BytesIO()
    
    with pd.ExcelWriter(buffer, engine='xlsxwriter') as writer:
        for idx, data_json in enumerate(stored_data_list):
            df = pd.read_json(io.StringIO(data_json), orient='split')
            
            # Use original name_list unchanged; only sanitize for Excel sheet name here:
            sheet_name = "Sheet" if not name_list or idx >= len(name_list) else name_list[idx]
            sheet_name = sanitize_sheet_name(sheet_name)  # Only modify for Excel
            
            df.to_excel(writer, sheet_name=sheet_name, index=False)
            
            worksheet = writer.sheets[sheet_name]
            
            # Autofit columns
            for i, col in enumerate(df.columns):
                max_len = max(df[col].astype(str).map(len).max(), len(col)) + 2
                worksheet.set_column(i, i, max_len)

    buffer.seek(0)
    return dcc.send_bytes(buffer.read(), "all_data.xlsx")


@app.callback(
    Output('auc-calibration-regression', 'data'),          # regression table
    Output('final-sample-table', 'data'),                  # sample table
    Output('final-sample-table', 'style_data_conditional'),# sample table style
    Output('final-sample-table', 'tooltip_data'),          # sample table tooltips
    Output('blank-sample', 'data'),                        # blank sample storage
    Input('add-to-table', 'n_clicks'),
    Input('perform-sample-conc-det', 'n_clicks'),
    Input('auc-calibration-regression', 'data_timestamp'),
    State('selection-store', 'data'),
    State('stored-data-list', 'data'),
    State('current-plot-index', 'data'),
    State('auc-calibration-regression', 'data'),
    State('final-sample-table', 'data'),
    State('calibration-data-checkbox', 'value'),
    State('rw1-names', 'data'),
    State('concentration-autofill', 'value'),
    State('blank-checkbox', 'value'),
    State('blank-sample', 'data'),
    prevent_initial_call=True
)
def update_tables(add_clicks, calc_clicks, data_timestamp, selectedData, data_list, current_index,
                  regression_data, sample_data, checkbox_value, name_list, autofill_value, blank_checkbox, blank_sample):

    ctx = callback_context
    if not ctx.triggered:
        raise PreventUpdate

    triggered_id = ctx.triggered_id
    regression_data = regression_data or []
    sample_data = sample_data or []
    blank_sample = blank_sample or []

    # --- Handle autofill concentration update ---
    if triggered_id == 'auc-calibration-regression':
        if not regression_data or 'enabled' not in autofill_value:
            raise PreventUpdate

        try:
            base_conc = float(regression_data[0].get('concentration'))
        except (TypeError, ValueError):
            raise PreventUpdate

        for i, row in enumerate(regression_data):
            row['concentration'] = round(base_conc * (2 ** i), 6)

        return regression_data, dash.no_update, dash.no_update, dash.no_update, dash.no_update

    # --- Handle ADD AUC button ---
    if triggered_id == 'add-to-table':
        if not selectedData or not data_list or current_index is None or current_index >= len(data_list):
            raise PreventUpdate
        if 'range' not in selectedData or 'x' not in selectedData['range']:
            raise PreventUpdate

        df = pd.read_json(io.StringIO(data_list[current_index]), orient='split')
        x0, x1 = selectedData['range']['x']
        mask = (df["Time"] >= x0) & (df["Time"] <= x1)
        selected_df = df[mask]

        if selected_df.empty:
            raise PreventUpdate

        if blank_sample and 'area' in blank_sample[0] and 'retention_time' in blank_sample[0]:
            auc = simpson(selected_df["Signal"], selected_df["Time"]) - blank_sample[0]['area']
        else:
            auc = simpson(selected_df["Signal"], selected_df["Time"])

        peak_idx = selected_df["Signal"].idxmax()
        retention_time = df.loc[peak_idx, "Time"]

        name = name_list[current_index] if name_list and current_index < len(name_list) else f"RW1 File {current_index + 1}"

        new_entry = {
            "sample_name": name,
            "retention_time": round(retention_time, 2),
            "area": round(auc, 2),
            "concentration": ""
        }

        # Add to appropriate table
        if 'enabled' in (blank_checkbox or []):
            blank_sample = [{
                "retention_time": round(retention_time, 6),
                "area": round(auc, 6)
            }]
            return regression_data, sample_data, [], [], blank_sample


        elif 'enabled' in (checkbox_value or []):
            regression_data.append(new_entry)

            if 'enabled' in (autofill_value or []) and regression_data:
                try:
                    base_conc = float(regression_data[0].get('concentration'))
                    for i, row in enumerate(regression_data):
                        row['concentration'] = round(base_conc * (2 ** i), 6)
                except Exception:
                    pass

            return regression_data, sample_data, [], [], blank_sample

        else:
            sample_data.append(new_entry)
            return regression_data, sample_data, [], [], blank_sample

    # --- Handle CALCULATE SAMPLE CONCENTRATIONS ---
    if triggered_id == 'perform-sample-conc-det':
        if not sample_data or not regression_data:
            raise PreventUpdate

        valid_cal = [r for r in regression_data if r.get('concentration') not in [None, ""] and r.get('area') not in [None, ""]]
        if len(valid_cal) < 2:
            raise PreventUpdate

        try:
            x = np.array([float(r['concentration']) for r in valid_cal]).reshape(-1, 1)
            y = np.array([float(r['area']) for r in valid_cal])
        except Exception:
            raise PreventUpdate

        model = LinearRegression().fit(x, y)
        slope, intercept = model.coef_[0], model.intercept_
        min_area, max_area = np.min(y), np.max(y)

        updated_samples, style_cond, tooltips = [], [], []

        for i, row in enumerate(sample_data):
            tooltip_row = {}
            try:
                area = float(row['area'])
                conc = (area - intercept) / slope
                row['concentration'] = round(conc, 4)

                if area < min_area or area > max_area:
                    style_cond.append({
                        'if': {'row_index': i},
                        'backgroundColor': '#ff1a1a',
                        'color': 'white'
                    })
                    tooltip_row = {
                        'area': {
                            'value': f"Area {area:.3f} outside calibration range ({min_area:.3f}–{max_area:.3f})",
                            'type': 'text'
                        }
                    }
                else:
                    tooltip_row = {'area': {'value': '', 'type': 'text'}}
            except Exception:
                row['concentration'] = ""
                tooltip_row = {'area': {'value': '', 'type': 'text'}}

            updated_samples.append(row)
            tooltips.append(tooltip_row)

        style_cond = [{'if': {'row_index': 'odd'}, 'backgroundColor': '#d5edfa'}] + style_cond
        return regression_data, updated_samples, style_cond, tooltips, blank_sample

    raise PreventUpdate



@app.callback(
    Output("regression-plot", "figure"),
    Output("regression-summary-table", "data"),
    Input("perform-regression", "n_clicks"),
    State("auc-calibration-regression", "data"),
    prevent_initial_call=True
)
def perform_linear_regression(n_clicks, regression_data):
    if not regression_data or len(regression_data) < 2:
        return go.Figure(), []

    try:
        concentrations = np.array([
            float(row["concentration"]) for row in regression_data
            if row["concentration"] != "" and row["area"] != ""
        ]).reshape(-1, 1)

        areas = np.array([
            float(row["area"]) for row in regression_data
            if row["concentration"] != "" and row["area"] != ""
        ])
    except ValueError:
        return go.Figure(), []

    if len(concentrations) != len(areas) or len(concentrations) < 2:
        return go.Figure(), []

    # Fit sklearn model for prediction line
    model = LinearRegression()
    model.fit(concentrations, areas)
    slope = model.coef_[0]
    intercept = model.intercept_

    # Statsmodels for stats and CI
    X_sm = sm.add_constant(concentrations)
    model_sm = sm.OLS(areas, X_sm).fit()

    r_squared = model_sm.rsquared
    n = len(areas)
    p = 1
    r_squared_adj = 1 - (1 - r_squared) * (n - 1) / (n - p - 1)

    se_intercept = model_sm.bse[0]
    se_slope = model_sm.bse[1]

    df_resid = n - p - 1
    t_val = t.ppf(1 - 0.025, df_resid)  # two-sided 95% CI

    intercept_ci_lower = intercept - t_val * se_intercept
    intercept_ci_upper = intercept + t_val * se_intercept
    slope_ci_lower = slope - t_val * se_slope
    slope_ci_upper = slope + t_val * se_slope

    y_pred = model_sm.fittedvalues
    ss_total = np.sum((areas - np.mean(areas))**2)
    ss_regression = np.sum((y_pred - np.mean(areas))**2)
    ss_residual = np.sum((areas - y_pred)**2)

    # Prepare regression summary table data
    regression_summary = [
        {"parameter": "Slope", "value": slope, "CI Lower (95%)": slope_ci_lower, "CI Upper (95%)": slope_ci_upper},
        {"parameter": "Intercept", "value": intercept, "CI Lower (95%)": intercept_ci_lower, "CI Upper (95%)": intercept_ci_upper},
        {"parameter": "R squared", "value": r_squared, "CI Lower (95%)": None, "CI Upper (95%)": None},
        {"parameter": "Adjusted R squared", "value": r_squared_adj, "CI Lower (95%)": None, "CI Upper (95%)": None},
        {"parameter": "SS Total", "value": ss_total, "CI Lower (95%)": None, "CI Upper (95%)": None},
        {"parameter": "SS Regression", "value": ss_regression, "CI Lower (95%)": None, "CI Upper (95%)": None},
        {"parameter": "SS Residual", "value": ss_residual, "CI Lower (95%)": None, "CI Upper (95%)": None},
    ]

    # Plotting part (same as before)
    x_min = concentrations.min()
    x_max = concentrations.max()
    num_ticks = 12
    tick_vals = np.linspace(x_min, x_max, num_ticks)
    tick_text = [f"{val:.2f}" for val in tick_vals]

    x_range = np.linspace(x_min, x_max, 100).reshape(-1, 1)
    y_fit = model.predict(x_range)

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=concentrations.flatten(),
        y=areas,
        mode='markers',
        name='Data',
        marker=dict(size=8, color='blue')
    ))
    fig.add_trace(go.Scatter(
        x=x_range.flatten(),
        y=y_fit,
        mode='lines',
        name='Regression Line',
        line=dict(color='red', width=2)
    ))

    fig.update_layout(
        title=dict(
            text="Calibration Curve",
            font=dict(size=26)
        ),
        xaxis=dict(
            title="Concentration (µg/mL)",
            title_font=dict(size=22),
            tickfont=dict(size=18),
            tickvals=tick_vals,
            ticktext=tick_text,
            ticks="outside",
            linecolor="black",
            linewidth=1,
        ),
        yaxis=dict(
            title="Area under the peak (a.u.)",
            title_font=dict(size=22),
            tickfont=dict(size=18),
            linecolor="black",
            linewidth=1,
        ),
        font=dict(size=16),
        dragmode="select",
        selectdirection="h",
        showlegend=False,
        height=650,
        template="plotly_white"
    )

    sign = "+" if intercept >= 0 else "-"
    summary_text = f"y = {slope:.1f} X {sign} {abs(intercept):.1f}<br>R² = {r_squared:.5f}"

    fig.add_annotation(
        text=summary_text,
        xref="paper", yref="paper",
        x=0.05, y=1,
        xanchor="left", yanchor="top",
        showarrow=False,
        align="right",
        font=dict(size=18, color="black"),
        bordercolor="black",
        borderwidth=1,
        borderpad=4,
        bgcolor="white",
        opacity=0.65
    )

    return fig, regression_summary

# ---------- NEW Callback: Show/Hide and Compute Dilution Table ----------
@app.callback(
    Output('dilution-factor', 'style'),
    Input('dilution-checkbox', 'value')
)
def show_dilution_input(checkbox_val):
    if 'enabled' in checkbox_val:
        return {'width': '100%', 'marginTop': '10px'}
    return {'display': 'none'}


@app.callback(
    Output('dilution-corrected-sample-table', 'data'),
    Input('dilution-factor', 'value'),
    State('final-sample-table', 'data'),
    prevent_initial_call=True
)
def compute_concentration_before_dilution(factor, sample_data):
    if factor is None or not sample_data:
        raise PreventUpdate
    return [
        {
            'sample_name': row['sample_name'],
            'dilution_factor': factor,
            'concentration_bef_dil': round(row['concentration'] * factor, 2) if row['concentration'] is not None else None
        } for row in sample_data
    ]

@app.callback(
    Output('dilution-corrected-wrapper', 'style'),
    Input('dilution-checkbox', 'value')
)
def toggle_dilution_table_visibility(checkbox_val):
    if 'enabled' in checkbox_val:
        return {'display': 'block', 'marginTop': '20px'}
    return {'display': 'none'}


@app.callback(
    Output('dilution-corrected-sample-table', 'data', allow_duplicate=True),
    Input('dilution-corrected-sample-table', 'data_timestamp'),
    State('dilution-corrected-sample-table', 'data'),
    State('final-sample-table', 'data'),
    prevent_initial_call='initial_duplicate'
)
def recalculate_on_edit(timestamp, table_data, final_sample_data):
    # Create a lookup of concentrations from final sample table
    sample_to_conc = {row['sample_name']: row['concentration'] for row in final_sample_data}

    updated_data = []
    for row in table_data:
        try:
            factor = float(row['dilution_factor'])
            conc = float(sample_to_conc.get(row['sample_name'], 0))
            row['concentration_bef_dil'] = round(conc * factor, 2)
        except Exception:
            row['concentration_bef_dil'] = None
        updated_data.append(row)

    return updated_data

# ---------- Callback: Export to Excel with autofit and image ----------
@app.callback(
    Output("download-excel-analyzed-sample-table", "data"),
    Input("export-excel-analyzed-sample-table", "n_clicks"),
    State("final-sample-table", "data"),
    State("auc-calibration-regression", "data"),
    State("regression-summary-table", "data"),
    State("blank-sample", "data"),
    State("dilution-corrected-sample-table", "data"),
    prevent_initial_call=True,
)
def export_to_excel(n_clicks, sample_data, calibration_data, regression_summary, blank_sample, dilution_corrected_data):

    if not sample_data and not calibration_data and not regression_summary:
        raise PreventUpdate

    buffer = io.BytesIO()

    # Create Matplotlib plot image from regression data
    img_buffer = io.BytesIO()
    if calibration_data and len(calibration_data) >= 2:
        # Prepare data
        df_cal = pd.DataFrame(calibration_data)
        try:
            x = df_cal['concentration'].astype(float).values.reshape(-1, 1)
            y = df_cal['area'].astype(float).values
            model = LinearRegression()
            model.fit(x, y)
            y_pred = model.predict(x)
        except Exception:
            x = y = y_pred = None

        if x is not None and y is not None:
            plt.figure(figsize=(7, 4))
            plt.scatter(x, y, color='blue', label='Data', s=50)
            plt.plot(x, y_pred, color='red', linewidth=1, label='Regression Line')
            plt.xlabel('Concentration (µg/mL)', fontsize=12)
            plt.ylabel('Area under the peak (a.u.)', fontsize=12)
            plt.title('Calibration Curve', fontsize=14)
            plt.legend()
            plt.tight_layout()
            plt.grid(True)
            ax = plt.gca()
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.savefig(img_buffer, format='png')
            plt.close()
            img_buffer.seek(0)
        else:
            img_buffer = None
    else:
        img_buffer = None

    with pd.ExcelWriter(buffer, engine='xlsxwriter') as writer:
        workbook = writer.book

        # Calibration Data sheet

        if blank_sample:
            blank = pd.DataFrame(blank_sample)
            blank.columns = ["Retention Time (min)", "Area (AU·min)"]
            blank.to_excel(writer, sheet_name="Calibration_Data", index=False)
            worksheet = writer.sheets["Calibration_Data"]
            # Autofit columns
            for x, col in enumerate(blank.columns):
                max_len = max(blank[col].astype(str).map(len).max(), len(col)) + 2
                worksheet.set_column(x, x, max_len)
            
            add = 3

        else:
            add = 0

        if calibration_data:
            df_cal = pd.DataFrame(calibration_data)
            df_cal.columns = ["Sample name", "Retention Time (min)", "Area (AU·min)", "Concentration (µg/mL)"]
            df_cal.to_excel(writer, sheet_name="Calibration_Data", index=False, startcol = add)
            worksheet = writer.sheets["Calibration_Data"]
            # Autofit columns
            for i, col in enumerate(df_cal.columns):
                max_len = max(df_cal[col].astype(str).map(len).max(), len(col)) + 2
                worksheet.set_column(i+add, i+add, max_len)

        # Regression Summary sheet or block (optional)
        if regression_summary:
            df_reg = pd.DataFrame(regression_summary)
            df_reg.columns = ["Parameter", "Value", "Lower 95% Confidence interval", "Upper 95% Confidence interval"]
            df_reg.to_excel(writer, sheet_name="Calibration_Data", index=False, startrow=0, startcol=6+add)
            worksheet = writer.sheets["Calibration_Data"]
            for j, col in enumerate(df_reg.columns):
                max_len = max(df_reg[col].astype(str).map(len).max(), len(col)) + 2
                worksheet.set_column(add+ 6 + j, add + 6 + j, max_len)

        # Insert matplotlib image below data
        if img_buffer:
            worksheet = writer.sheets["Calibration_Data"]
            start_row = max(len(calibration_data), len(regression_summary) if regression_summary else 0) + 3
            worksheet.insert_image(start_row, 0, 'regression_plot.png', {'image_data': img_buffer, 'x_scale': 0.7, 'y_scale': 0.7})

        # Sample Data sheet
        if sample_data:
            df_sample = pd.DataFrame(sample_data)
            df_sample.columns = ["Sample name", "Retention Time (min)", "Area (AU·min)", "Concentration (µg/mL)"]
            df_sample.to_excel(writer, sheet_name="Sample_Data", index=False)
            worksheet = writer.sheets["Sample_Data"]
            for i, col in enumerate(df_sample.columns):
                max_len = max(df_sample[col].astype(str).map(len).max(), len(col)) + 4
                worksheet.set_column(i, i, max_len)


        if dilution_corrected_data:
            df_sample_corr = pd.DataFrame(dilution_corrected_data)
            df_sample_corr.columns = ["Sample name", "Dilution factor", "Concentration before dilution (µg/mL)"]
            df_sample_corr["Dilution factor"] = pd.to_numeric(df_sample_corr["Dilution factor"], errors='coerce')
            df_sample_corr.to_excel(writer, sheet_name="Sample_Data", index=False, startrow=len(df_sample) + 2)
            worksheet = writer.sheets["Sample_Data"]
            for i, col in enumerate(df_sample_corr.columns):
                max_len = max(df_sample_corr[col].astype(str).map(len).max(), len(col)) + 2
                worksheet.set_column(i, i, max_len)
        

    buffer.seek(0)
    return dcc.send_bytes(buffer.read(), "integrated_areas.xlsx")

@app.callback(Output("keepalive", "n_intervals"), Input("keepalive", "n_intervals"))
def update_last_activity(n):
    global last_activity_time
    last_activity_time = datetime.now()
    return n



if __name__ == "__main__":
    port = 8050
    url = f"http://127.0.0.1:{port}"

    webbrowser.open(url, new=2)  # Open in a new tab
    app.run(debug=False, use_reloader=False)


# pyinstaller "C:\Users\u0155764\OneDrive - KU Leuven\Tom Konings\Software\HPLC data analyzer app\Code\HPLCDataAnalyzer.py" --onefile --windowed --hidden-import dash,plotly,scipy,sklearn,pandas,numpy,xlsxwriter,matplotlib --icon="C:\Users\u0155764\OneDrive - KU Leuven\Tom Konings\Software\HPLC data analyzer app\Code\icon.ico"

