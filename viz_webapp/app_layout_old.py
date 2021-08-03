app.layout = html.Div(children=[
    html.H1(children='Dashboard'),
    html.Div(children=dcc.Graph(id='bar_plot'), className='div1'),
    html.Div(children=dash_table.DataTable(id='data_table',
        columns=[{'id': 'Taxa', 'name': 'Taxa'}, {'id': "Matches", 'name': 'Matches'}],
        style_cell_conditional=[
            {
                'if': {'column_id': c},
                'textAlign': 'left'
            } for c in ['Date', 'Region']
        ],
        style_data_conditional=[
            {
                'if': {'row_index': 'odd'},
                'backgroundColor': 'rgb(248, 248, 248)'
            }
        ],
        style_header={
            'backgroundColor': 'rgb(230, 230, 230)',
            'fontWeight': 'bold'
        }
        ), className='div2'), 

    html.Div(children=[dcc.Dropdown(
            id='barcode_select',
            value='11'
            ),
        dcc.Interval(
            id='interval-component',
            interval=7*1000,
            n_intervals=0
            ),
        dcc.RadioItems(
            id='update-ratio',
            options=[
            {'label': '7 seconds', 'value': "7"},
            {'label': '20 seconds', 'value': "20"},
            {'label': '1 minute', 'value': "60"}],
            value="7"
            )])
    ],
    className='parent'
    #style = {'display': 'inline-block', 'width': '48%'}
)
