import datetime
import dash
from dash import dcc, html
import plotly
import plotly.express as px
from dash.dependencies import Input, Output
from datetime import datetime, timedelta
import pytz
from skyfield.api import EarthSatellite
from skyfield.api import load, wgs84
from sgp4.api import Satrec, WGS72
from astropy.coordinates import CartesianDifferential, CartesianRepresentation, TEME, GCRS
from astropy import units
from astropy.time import Time, TimeDelta
from sgp4.api import jday
from astropy import units as u
from astropy.coordinates import GCRS, EarthLocation, ITRS
from IPython import embed
import plotly.graph_objects as go
from sgp4.api import SGP4_ERRORS
import pandas as pd

utc_timezone = pytz.timezone('UTC')

import requests

def load_text_file_from_norad_id(norad_id, hardcode=True):
    if hardcode:
        return [
            '1 25544U 98067A   23161.05314815  .00013937  00000-0  24614-3 0  9993',
            '2 25544  51.6422   3.8974 0005491  67.8121  67.8666 15.50633390400689'
        ]
    url = f"https://celestrak.org/NORAD/elements/gp.php?CATNR={norad_id}&FORMAT=2LE"
    response = requests.get(url)
    if response.status_code == 200:
        lines = response.text.splitlines()
        return lines
    else:
        print("Failed to load the file from the URL:", url)
        return []


class SpacecraftPropagator:
    def __init__(self, tle_file=None):
        # do tle reading
        #line1 = '1 25544U 98067A   23160.65531729  .00014326  00000+0  25287-3 0  9992'
        #line2 = '2 25544  51.6424   5.8692 0005411  66.8319   6.5873 15.50623283400624'
        #self.satellite = Satrec.twoline2rv(line1, line2)
        self.ts = load.timescale()
        self.satellite =  Satrec.twoline2rv(*load_text_file_from_norad_id(25544, hardcode=True), WGS72)

    def propagate_sgp4(self, timestamp_datetime):
        time = Time(timestamp_datetime)
        #current_datetime = datetime.now(utc_timezone)
        #jd, fr = jday(current_datetime.year, current_datetime.month, current_datetime.day, current_datetime.hour, current_datetime.minute, current_datetime.second)
        error_code, teme_p, teme_v = self.satellite.sgp4(time.jd1, time.jd2)  # in km and km/s
        if error_code != 0:
            raise RuntimeError(SGP4_ERRORS[error_code])
        teme_p = CartesianRepresentation(teme_p*u.km)
        teme_v = CartesianDifferential(teme_v*u.km/u.s)
        teme = TEME(teme_p.with_differentials(teme_v), obstime=timestamp_datetime)
        itrs_geo = teme.transform_to(ITRS(obstime=timestamp_datetime))
        geodetic = itrs_geo.earth_location.geodetic
        #cart_gcrs = teme.transform_to(GCRS(obstime=timestamp_datetime)).cartesian
        #geodetic = EarthLocation.from_geocentric(cart_gcrs.x, cart_gcrs.y, cart_gcrs.z)
        lat = geodetic.lat.deg
        lon = geodetic.lon.deg
        return lat, lon

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div(
    html.Div([
        html.H4('Droid.001 Satellite Live Feed'),
        html.Div(id='live-update-text'),
        dcc.Graph(id='live-update-graph', config={'scrollZoom': True}),
        dcc.Interval(
            id='interval-component',
            interval= 1 * 1000,  # in milliseconds
            n_intervals=0
        )
    ])
)


@app.callback(Output('live-update-text', 'children'),
              Input('interval-component', 'n_intervals'))
def update_metrics(n):
    current_datetime = datetime.now(utc_timezone)
    lat, lon = sat.propagate_sgp4(current_datetime)
    style = {'padding': '-5px', 'fontSize': '16px'}
    return [
        html.Span(f'Longitude: {lon:.3f}° ', style=style),
        html.Span(f'Latitude: {lat:.3f}° ', style=style),
    ]

# Multiple components can update every time interval gets fired.
data = pd.DataFrame(columns=['Timestamp', 'Latitude', 'Longitude'])

# Collect some data
@app.callback(Output('live-update-graph', 'figure'),
              Input('interval-component', 'n_intervals'))
def update_graph_live(n):
    current_datetime = datetime.now(utc_timezone)
    lat, lon = sat.propagate_sgp4(current_datetime)
    embed()
    data = data.append({'Timestamp': current_datetime, 'Latitude': lat, 'Longitude': lon}, ignore_index=True)

    # Create a list of dictionaries for the ground station information
    ground_stations = [
        {'name': 'Bahrain', 'lat': 26.05138, 'lon': 50.50307, 'color': "red"},
        {'name': 'Dubbo', 'lat': -32.17723, 'lon': 148.615407, 'color': "green"},
        {'name': 'Ireland', 'lat': 53.40667, 'lon': -6.22438, 'color': "orange"},
        {'name': 'CapeTown', 'lat': -34.02762, 'lon': 18.71868, 'color': "cyan"},
        {'name': 'Punta Arenas', 'lat': -52.94111, 'lon': -70.84972, 'color': "magenta"},
    ]

    fig = px.line_geo(lat="Latitude", lon="Longitude",hover_data="Timestamp", markers=True)

    # Add ground stations as scatter markers
    for station in ground_stations:
        fig.add_trace(go.Scattergeo(
            lat=[station['lat']],
            lon=[station['lon']],
            mode='markers',
            marker=dict(
                size=10,
                color=station['color'],
                opacity=0.7
            ),
            name=station['name'],
            hovertemplate=f'Name: {station["name"]}<br>Latitude: {station["lat"]:.3f}°<br>Longitude: {station["lon"]:.3f}°<extra></extra>',
        ))

    fig.update_layout(
        width=1900,  # Width of the figure (adjust to your screen resolution)
        height=950,  # Height of the figure (adjust to your screen resolution)
        uirevision=1
    )
    return fig


if __name__ == '__main__':
    sat = SpacecraftPropagator()
    app.run_server(debug=True)











