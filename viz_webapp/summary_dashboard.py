# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from functools import reduce
import dash_table
import plotly.graph_objs as go
import pandas as pd
import numpy as np
import os
import csv
import skbio
import plotly.express as px


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

SIDEBAR_STYLE = {
    'position': 'fixed',
    'top': 100,
    'left': 0,
    'bottom': 100,
    'width': '20%',
    'padding': '20px 10px',
    'background-color': '#f8f9fa'
}

# the style arguments for the main content page.
CONTENT_STYLE = {
    'margin-left': '25%',
    'margin-right': '5%',
    'padding': '20px 10p'
}

TEXT_STYLE = {
    'textAlign': 'center',
    'color': '#191970'
}

CARD_TEXT_STYLE = {
    'textAlign': 'center',
    'color': '#0074D9'
}

CARD_QC_TEXT_STYLE = {
    'textAlign': 'center',
    'color': '#228B22'
}

#TO DO
#Fetch run data from config files


#DB query
#client = MongoClient(host='mongo')
#db = client['test_1_barcodes']

#Dash webapp components
#app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

cards = []
for barcode in next(os.walk("data/"))[1]:
    qc_data = pd.read_csv("data/" + str(barcode) + "/qc_report.csv")
    bars = [go.Bar(y=[barcode],
                x=[qc_data.iloc[0, i]],
                name=str(qc_data.iloc[0, i]),
                width=0.5,
                orientation='h') 
                for i in range(1, 5)]
    graph = go.Figure(bars,layout=go.Layout(barmode='stack'))
    card = dbc.Row(
        dbc.Card(
            dbc.CardBody(
                [
                    html.H4(barcode),
                    graph,
                ]
            )
            #style={"width": "18rem"},
        )
    )
    cards.append(card)

controls = dbc.FormGroup(
    [
        # html.P('Barcode/sample select', style={
        #     'textAlign': 'center'
        # }),
        # dcc.Dropdown(
        #     id='barcode_select',
        #     value='11'
        #     ),
        html.P('Automatic update ratio', style={
            'textAlign': 'center'
        }),
        dcc.Interval(
            id='interval-component',
            interval=7*1000,
            n_intervals=0
            ),
        dbc.Card([dbc.RadioItems(
            id='update-ratio',
            options=[
            {'label': '7 seconds', 'value': "7"},
            {'label': '20 seconds', 'value': "20"},
            {'label': '1 minute', 'value': "60"}],
            value="7",
            inline = True,
            style={
                'margin': 'auto'
            }
        )]),
        html.Br(),
        dbc.Button(
            id='Update',
            n_clicks=0,
            children='update',
            color='primary',
            block=True
        ),
    ]
)

button_group = dbc.ButtonGroup(
    [
        dbc.Button("Overview", id="overview", color="success", block=True, style={'textAlign': 'center', 'width': '100%'} ),
        dbc.Button("barcode01", id="bc01",color="primary"),
        dbc.Button("barcode02", id="bc02",color="primary"),
        dbc.Button("barcode11", id="bc11",color="primary"),
    ],
    style ={'textAlign': 'center', 'width': '100%'},
    vertical=True,
    size = 'lg'
)

bc_radiolist = dbc.RadioItems(
    id="bc_radiochecklist",
    options=[{"label": barcode, "value":barcode} for barcode in next(os.walk("data/"))[1]],
    # [
    #     {"label": "barcode01", "value": 'bc01'},
    #     {"label": "barcode02", "value": 'bc02'},
    #     {"label": "barcode11", "value": 'bc11'},
    # ],
    labelCheckedStyle={"color": "green"},
    value='barcode01'
)

generate_otu = html.Div(
    [
        dbc.Button("Generate OTU file", color='primary', id="alert-toggle-auto"),
        html.Hr(),
        dbc.Alert(
            "OTU file successfully generated in viz_webapp/data",
            id="alert-auto",
            is_open=False,
            duration=4000,
        ),
    ]
)

sidebar = html.Div(
    [
        html.H3('Parameters', style=TEXT_STYLE),
        html.Hr(),
        controls,
        html.H3('Samples', style=TEXT_STYLE),
        bc_radiolist,
        html.Hr(),
        generate_otu
    ],
    style=SIDEBAR_STYLE,
)


navbar = dbc.NavbarSimple(
        children=[
            dbc.NavItem(dbc.NavLink("Home", href="/index")),
            dbc.DropdownMenu(
                children=[
                    dbc.DropdownMenuItem("Page 2", href="#"),
                    dbc.DropdownMenuItem("Page 3", href="#"),
                ],
                nav=True,
                in_navbar=True,
                label="More",
            ),
        ],
        brand="NanoRTax",
        color="primary",
        dark=True,
    )

content = html.Div(
    [
        sidebar,
        cards,
    ],
    style=CONTENT_STYLE
)

app = dash.Dash(external_stylesheets=[dbc.themes.SPACELAB])
app.layout = html.Div([navbar, content])

#UPDATE DROPDOWN
# @app.callback(Output('barcode_select', 'options'),
#               [Input('interval-component', 'n_intervals')])
# def update_barcodesel(n):
#     #barcodes = sorted(list(db.kraken2.distinct("barcode")))
#     barcodes = ["1"]
#     dropdown_options = [{'label': 'barcode' + i, 'value': i} for i in barcodes]
#     return dropdown_options

#UPDATE REFRESH RATIO
@app.callback(Output('interval-component', 'interval'),
              [Input('update-ratio', 'value')])
def update_refresh_ratio(n):
    return (int(n) * 1000)

if __name__ == '__main__':
    app.run_server(debug=True)