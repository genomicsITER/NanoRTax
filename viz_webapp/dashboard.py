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
    'top': 57,
    'left': 0,
    'bottom': 0,
    'width': '15%',
    'padding': '5px 10px',
    'background-color': '#f8f9fa'
}

# the style arguments for the main content page.
CONTENT_STYLE = {
    'margin-left': '16%',
    'margin-right': '1%',
    'padding': '0px 0px'
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


bc_radiolist = dbc.RadioItems(
    id="bc_radiochecklist",
    options=[{"label": barcode, "value":barcode} for barcode in next(os.walk("data/"))[1]],
    # [
    #     {"label": "barcode01", "value": 'bc01'},
    #     {"label": "barcode02", "value": 'bc02'},
    #     {"label": "barcode11", "value": 'bc11'},
    # ],
    labelCheckedStyle={"color": "green"},
    value='test_data'
)

generate_otu = html.Div(
    [
        dbc.Button("Generate OTU file", color='info', id="alert-toggle-auto", block=True),
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
        html.P(),
        dbc.Button("Switch view", color='info', id="sample_summary", active=False, block=True),
        html.P(),
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

#Index components:

cards = []
colors = ['green', 'crimson','red', 'darkred']
for barcode in next(os.walk("data/"))[1]:
    qc_data = pd.read_csv("data/" + str(barcode) + "/qc_report.csv")
    bars = [go.Bar(y=[barcode],
                x=[qc_data.iloc[0, i]],
                name=str(qc_data.columns[i]),
                width=0.5,
                orientation='h',
                marker_color=colors[i-1]) 
                for i in range(1, 5)]
    graph = go.Figure(bars,layout=go.Layout(barmode='stack'))
    graph_component = dcc.Graph(figure=graph, style={"height": "260px"})
    qc_data_table = dash_table.DataTable(
        data=qc_data.to_dict('records'),
        columns=[{'id': c, 'name': c} for c in qc_data.columns],
        style_cell={'textAlign': 'left', 'padding': '5px'},
        style_as_list_view=True,
        style_header={
            'fontWeight': 'bold'
        },
    )
    card = dbc.Row(
        dbc.Col(
            html.Div(
                dbc.Card(
                    [dbc.CardHeader(barcode),
                    dbc.CardBody(
                        [
                            qc_data_table,
                            graph_component,
                        ]
                    )], className="mb-3",color="primary", outline=True
                    #style={"width": "18rem"},
                )
            ), width=12
        )
    )
    cards.append(card)

index_diversity_curve = dbc.Row(
    [
        dbc.Col(
            html.Div(
                dbc.Card(
                    [dbc.CardHeader("Diversity Index for all samples"),
                    dbc.CardBody(
                        [
                            dcc.Graph(id='diversity_curve_plot_all',style={"display": "block"}),
                            dbc.Form(
                            [
                                dbc.FormGroup(
                                    [
                                        dbc.Label("Read interval", className="mr-2"),
                                        dbc.Input(type="text", placeholder="20", value="20", id="chunk_size_all_select"),
                                    ],
                                ),
                                dbc.FormGroup(
                                        [
                                            dbc.Label("Diversity Index  "),
                                            dbc.RadioItems(
                                                options=[
                                                    {"label": "Shannon", "value": 1},
                                                    {"label": "Simpson", "value": 2},
                                                ],
                                                value=1,
                                                id="diversity_index_all_checklist",
                                                inline = True,
                                            ),
                                        ],
                                ),
                                dbc.FormGroup(
                                        [
                                            dbc.Label("Tax level"),
                                            dbc.RadioItems(
                                                options=[
                                                    {"label": "Family", "value": "family"},
                                                    {"label": "Genus", "value": "genus"},
                                                    {"label": "Species", "value": "species"},
                                                ],
                                                value="genus",
                                                id="diversity_tax_level_checklist",
                                                inline = True,
                                            ),
                                        ],
                                ),
                                dbc.FormGroup(
                                        [
                                            dbc.Label("Classifier"),
                                            dbc.RadioItems(
                                                options=[
                                                    {"label": "Kraken2", "value": "kraken"},
                                                    {"label": "Centrifuge", "value": "centrifuge"},
                                                    {"label": "BLAST", "value": "blast"},
                                                ],
                                                value="kraken",
                                                id="diversity_classifier_checklist",
                                                inline = True,
                                            ),
                                        ],
                                )
                            ],
                            inline=False,
                        )
                        ]
                    )], className="mb-3",color="primary", outline=True
                    #style={"width": "18rem"},
                )
            ), width="12"
            ), 
    
    ]
)

#Sample page components

content_qc = dbc.Row(
    [
        dbc.Col([
        html.H4(children=['Quality control information']),
        html.Div(id='qc_table'),
        ],
        md=12
    )
    ]
)

content_first_row = html.Div(
    [
    dbc.Row(dbc.Col(
        html.Div(
                dbc.Card(
                    [dbc.CardHeader("Taxonomic classiffication for " + barcode),
                    dbc.CardBody(
                        [
                            dcc.Graph(id='bar_plot'),
                            dbc.Form(
                                dbc.FormGroup(
                                    [
                                    dbc.Label("\"Other\" label threshold (%)", className="mr-2"),
                                    dbc.Input(type="text", placeholder="1", value="1", id="other_threshold"),
                                    dbc.Label("  Tax level:   "),
                                    dbc.RadioItems(
                                        options=[
                                            {"label": "Family", "value": 1},
                                            {"label": "Class", "value": 2},
                                            {"label": "Order", "value": 3},
                                            {"label": "Genus", "value": 4},
                                            {"label": "Species", "value": 5},
                                        ],
                                        value=1,
                                        id="tax_level",
                                        inline=True,
                                        ),
                                    ],
                                ),inline=True
                            )
                            ,
                        ]
                    )], className="mb-3",color="primary", outline=True
                    #style={"width": "18rem"},
                )
            ), width=12
        )),
    ]
)

content_diversity_curve = dbc.Row([
    dbc.Col(html.Div(
        dbc.Card(
                    [dbc.CardHeader("Diversity Index"),
                    dbc.CardBody(
                        [
                            dcc.Graph(id='diversity_curve_plot'),
                            dbc.Form(
                            [
                                dbc.FormGroup(
                                    [
                                        dbc.Label("Read interval", className="mr-2"),
                                        dbc.Input(type="text", placeholder="20", value="20", id="chunk_size_select"),
                                    ],
                                ),
                                dbc.FormGroup(
                                        [
                                            dbc.Label("Diversity Index  "),
                                            dbc.Checklist(
                                                options=[
                                                    {"label": "Shannon", "value": 1},
                                                    {"label": "Simpson", "value": 2},
                                                    {"label": "Taxa count", "value": 3},
                                                ],
                                                value=[1],
                                                id="diversity_index_checklist",
                                                inline = True,
                                            ),
                                        ],
                                ),
                            ],
                            inline=False,
                        )
                        ]
                    )],color="primary", outline=True
                )
            ), width=7), 
    dbc.Col(html.Div(
                dbc.Card(
                    [dbc.CardHeader('Relative Abundance'),
                    dbc.CardBody(
                        [
                         html.Div(id='data_table'),
                        ]
                    )],color="primary", outline=True
                )
            ), width=5),
    ])



tax_info = html.Div(
    [
        dbc.Row([dbc.Col(
                [html.H4(children=['Diversity Index']),
                html.Div(id='diversity_table'),
                html.P(),
                content_diversity_curve],
                width=12,
                )
            ]
        )
    ]   
)

tab_content = html.Div(
    [
        html.P(),
        content_qc,
    ],

)

tabs = html.Div(
    [
        dbc.Tabs(
            [
                dbc.Tab(label="Kraken 2", tab_id="kraken"),
                dbc.Tab(label="Centrifuge", tab_id="centrifuge"),
                dbc.Tab(label="BLAST", tab_id="blast"),
            ],
            id="tabs",
            active_tab="kraken",
        ),
        
    ]
)

index_layout = html.Div(
    [
        sidebar,
        tabs,
        html.P(),
        html.Div(cards),
        index_diversity_curve,
    ],
    style=CONTENT_STYLE
)

content_layout = html.Div(
    [
        sidebar,
        tabs,
        content_qc,
        html.Div(id="content"),
        html.P(),
        content_first_row,
        html.P(),
        tax_info,
    ],
    style=CONTENT_STYLE
)

app = dash.Dash(external_stylesheets=[dbc.themes.LUMEN],suppress_callback_exceptions=True)

app.layout = html.Div([
    html.Div(id='page-content', children=[navbar, content_layout])
])

#app.layout = html.Div([navbar, content])

#UPDATE DROPDOWN
# @app.callback(Output('barcode_select', 'options'),
#               [Input('interval-component', 'n_intervals')])
# def update_barcodesel(n):
#     #barcodes = sorted(list(db.kraken2.distinct("barcode")))
#     barcodes = ["1"]
#     dropdown_options = [{'label': 'barcode' + i, 'value': i} for i in barcodes]
#     return dropdown_options
@app.callback([Output('page-content', 'children'),Output('sample_summary', 'active')],
              [Input('sample_summary', 'n_clicks')],
              [State('sample_summary', 'active')])
def display_page(n_clicks, active):
    new_active = not active
    if(active):
        return [navbar, index_layout], new_active
    else:
        return [navbar, content_layout], new_active
#UPDATE REFRESH RATIO
@app.callback(Output('interval-component', 'interval'),
              [Input('update-ratio', 'value')])
def update_refresh_ratio(n):
    return (int(n) * 1000)

#Index callbacks:

def diversity_all_curve_graph(at, level, index, chunk_size):
    fig = go.Figure()
    for barcode in next(os.walk("data/"))[1]:
        diveristy_time_data = pd.read_csv("data/" + str(barcode) + "/" + at + "_report_full.txt", delimiter="\t",  names=['seq_id', 'tax_id', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
        diversity = []
        read_count = []
        max_reads = 0
        for pos in range(int(chunk_size), diveristy_time_data.shape[0], int(chunk_size)):
            df_subset = diveristy_time_data.iloc[0:pos]
            read_count.append(df_subset.shape[0])
            tax_table = df_subset[level].value_counts()
            tax_table_df = pd.DataFrame(list(zip(tax_table.index.tolist(), tax_table.tolist())), columns =['tax_id', 'read_count'])

            if(index == 1):
                shannon = skbio.diversity.alpha.shannon(tax_table_df['read_count'])
                diversity.append(shannon)
            else:
                simpson = skbio.diversity.alpha.simpson(tax_table_df['read_count'])
                diversity.append(simpson)
        if(len(read_count) > max_reads):
            max_reads = len(read_count)
            read_count_axis = read_count
        fig.add_trace(go.Scatter(x=read_count_axis, y=diversity,
                    mode='lines',
                    name=barcode))

    fig.update_layout(title='Diversity index curve',
                   xaxis_title='Reads',
                   yaxis_title='Diversity index',
                   )
    return fig

@app.callback(Output('diversity_curve_plot_all', 'figure'), 
              Input("diversity_classifier_checklist", "value"), Input("diversity_tax_level_checklist", "value"), 
              Input("diversity_index_all_checklist", "value"), Input("chunk_size_all_select", "value"))
def update_diversity_all_graph(at, level, index, chunk_size):
    return diversity_all_curve_graph(at, level, index, chunk_size)

#Sample classification summary callbacks:

def update_assets(barcode, tool, level, other_label):
    level_dict = {1: "family", 2: "class", 3: "order", 4: "genus", 5: "species"}
    top_taxa = pd.read_csv("data/" + str(barcode) + "/" + tool + "_report_" + level_dict[level] + ".csv")
    data_table_columns = [{"name": i, "id": i} for i in sorted(top_taxa.columns)],
    if(len(top_taxa)):
        grouped_data = top_taxa.groupby("tax_id",as_index=False).sum().sort_values(by="read_count", ascending= False).reset_index(drop=True)
        #Other label for 5% data
        if(float(other_label)>0):
            total_reads = grouped_data['read_count'].sum()
            for i, row in grouped_data.iterrows():
                tax_id_label = row['tax_id']
                if (row['read_count'] < (float(total_reads)*float(other_label))/100.0):
                    tax_id_label = "Other"
                grouped_data.at[i,'tax_id'] = tax_id_label
            grouped_data = grouped_data.groupby("tax_id",as_index=False).sum().sort_values(by="read_count", ascending= False).reset_index(drop=True)

        grouped_data['Abundance'] = round((grouped_data['read_count'] / grouped_data['read_count'].sum()) * 100, 3)
        grouped_data_fn = grouped_data.rename(columns={"tax_id": "Name", "read_count": "Reads", "Abundance": "Abundance (%)"})
        data_df = grouped_data_fn.to_dict('report')
        data = dash_table.DataTable(
            data=grouped_data_fn.to_dict('records'),
            columns=[{'id': c, 'name': c} for c in grouped_data_fn.columns],
            style_cell={'textAlign': 'left'},
            page_size=20,
            style_table={'height': '600px','overflowY': 'auto'},
            style_header={
            'fontWeight': 'bold'
        }
        )
        bars = [go.Bar(y=[barcode],
                x=[grouped_data.iloc[i, 1]],
                name=str(grouped_data.iloc[i, 0]),
                width=0.5,
                orientation='h') 
                for i in range(len(grouped_data['read_count']))]
        graph = go.Figure(bars,layout=go.Layout(barmode='stack'))
        return data, graph

def diversity_curve_graph(at, barcode, level, index, chunk_size):
    level_dict = {1: "family", 2: "class", 3: "order", 4: "genus", 5: "species"}
    diveristy_time_data = pd.read_csv("data/" + str(barcode) + "/" + at + "_report_full.txt", delimiter="\t",  names=['seq_id', 'tax_id', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
    diversity_shannon = []
    diversity_simpson = []
    taxa_count= []
    read_count = []
    for pos in range(int(chunk_size), diveristy_time_data.shape[0], int(chunk_size)):
        df_subset = diveristy_time_data.iloc[0:pos]
        read_count.append(df_subset.shape[0])
        tax_table = df_subset[level_dict[level]].value_counts()
        tax_table_df = pd.DataFrame(list(zip(tax_table.index.tolist(), tax_table.tolist())), columns =['tax_id', 'read_count'])

        if(1 in index):
            shannon = skbio.diversity.alpha.shannon(tax_table_df['read_count'])
            diversity_shannon.append(shannon)
        if(2 in index):
            simpson = skbio.diversity.alpha.simpson(tax_table_df['read_count'])
            diversity_simpson.append(simpson)
        if(3 in index):
            count = tax_table_df.shape[0]
            taxa_count.append(count)
    fig = go.Figure()

    if(1 in index):
            fig.add_trace(go.Scatter(x=read_count, y=diversity_shannon,
                        mode='lines',
                        name='Shannon index'))
    if(2 in index):
            fig.add_trace(go.Scatter(x=read_count, y=diversity_simpson,
                        mode='lines',
                        name='Simpson index'))
    if(3 in index):
            fig.add_trace(go.Scatter(x=read_count, y=taxa_count,
                        mode='lines',
                        name='Taxa count'))
    fig.update_layout(title='Diversity index curve',
                   xaxis_title='Reads',
                   yaxis_title='Index')
    return fig

#@app.callback(Output('diversity_curve_plot', 'figure'), 
#              Input("tabs", "active_tab"), Input("bc_radiochecklist", "value"), Input("tax_level", "value"), 
#              Input("diversity_index_checklist", "value"), Input("chunk_size_select", "value"))
#def update_diversity_graph(at, barcode, level, index, chunk_size):
#    return diversity_curve_graph(at, barcode, level, index, chunk_size)

@app.callback(
    Output("alert-auto", "is_open"),
    [Input("alert-toggle-auto", "n_clicks")],
    [State("alert-auto", "is_open"), State("tabs", "active_tab")],
)
def generate_otu(n, is_open, at):
    samples = []
    sample_names = next(os.walk("data/"))[1]
    for barcode in sample_names:
        df = pd.read_csv("data/" + barcode + "/"+ at +"_report_full_otu.txt", delimiter="\t", names=['read_id','tax_id','otu', 'count'])
        otu_count = df['otu'].value_counts()
        otu_df = pd.DataFrame(list(zip(otu_count.index.tolist(), otu_count.tolist())), columns = ['tax_id', 'read_count']).transpose()
        otu_df.rename(columns=otu_df.iloc[0], inplace = True)
        samples.append(otu_df.iloc[1:])

    df = reduce(lambda df1,df2: pd.merge(df1,df2,how='outer'), samples)
    df.transpose().set_axis(sample_names, axis='columns').to_csv("data/" + at +"_otu_table.txt")

    if n:
        return not is_open
    return is_open


@app.callback(Output('data_table', 'children'),Output('bar_plot', 'figure'),Output('diversity_table', 'children'), Output('qc_table', 'children'), Output('diversity_curve_plot', 'figure'), 
              Input("tabs", "active_tab"), Input("bc_radiochecklist", "value"), Input("tax_level", "value"), Input("diversity_index_checklist", "value"), Input("chunk_size_select", "value"), Input("other_threshold", "value"))
def switch_tab(at, barcode, level, index, chunk_size, other_label):
    data, graph = update_assets(barcode, at, level, other_label)
    level_dict = {1: "family", 2: "class", 3: "order", 4: "genus", 5: "species"}
    diversity_data = pd.read_csv("data/" + str(barcode) + "/" + at + "_diversity_full_" + level_dict[level] + ".csv", names=["Total reads","Shannon", "Simpson"])
    diversity_data_table = dash_table.DataTable(
        data=diversity_data.to_dict('records'),
        columns=[{'id': c, 'name': c} for c in diversity_data.columns],
        #columns=[{'id': 'Shannon', 'name': 'Shannon'}, {'id': 'Simpson', 'name': 'Simpson'}],
        style_cell={'textAlign': 'left', 'padding': '5px'},
        style_as_list_view=True,
        style_table={'height': 'auto', 'overflowY': 'auto'},
        page_size=20,
        style_header={
            'fontWeight': 'bold'
        },
    )

    qc_data = pd.read_csv("data/" + str(barcode) + "/qc_report.csv")
    qc_data_table = dash_table.DataTable(
        data=qc_data.to_dict('records'),
        columns=[{'id': c, 'name': c} for c in qc_data.columns],
        style_cell={'textAlign': 'left', 'padding': '5px'},
        style_as_list_view=True,
        style_header={
            'fontWeight': 'bold'
        },
    )

    diversity_time_fig = diversity_curve_graph(at, barcode, level, index, chunk_size)

    return data, graph, diversity_data_table, qc_data_table, diversity_time_fig

if __name__ == '__main__':
    app.run_server(debug=True)