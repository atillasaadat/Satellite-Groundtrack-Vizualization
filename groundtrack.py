import datetime
import dash
from dash import dcc, html
import plotly
import plotly.express as px
from dash.dependencies import Input, Output
from datetime import timedelta
from datetime import datetime as dt
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
import numpy as np
from geopy.distance import geodesic
from geopy.point import Point

utc_timezone = pytz.UTC
la_tz = pytz.timezone('America/Los_Angeles')

import requests

def load_text_file_from_norad_id(norad_id, hardcode=True):
    if hardcode:
        return [
            '1 00001U 00000AAA 23159.93196266  .00000000  00000-0  11824-3 0  9992',
            '2 00001 097.5053 272.5491 0010051 244.4995 125.7454 15.12446858    08',
        ]
    url = f"https://celestrak.org/NORAD/elements/gp.php?CATNR={norad_id}&FORMAT=2LE"
    response = requests.get(url)
    if response.status_code == 200:
        lines = response.text.splitlines()
        return lines
    else:
        print("Failed to load the file from the URL:", url)
        return []

def JulianDayFromDate(date,calendar='standard'):

    """
creates a Julian Day from a 'datetime-like' object.  Returns the fractional
Julian Day (resolution 1 second).

if calendar='standard' or 'gregorian' (default), Julian day follows Julian 
Calendar on and before 1582-10-5, Gregorian calendar after 1582-10-15.

if calendar='proleptic_gregorian', Julian Day follows gregorian calendar.

if calendar='julian', Julian Day follows julian calendar.

Algorithm:

Meeus, Jean (1998) Astronomical Algorithms (2nd Edition). Willmann-Bell,
Virginia. p. 63
    """
    # based on redate.py by David Finlayson.
    year=date.year; month=date.month; day=date.day
    hour=date.hour; minute=date.minute; second=date.second
    # Convert time to fractions of a day
    day = day + hour/24.0 + minute/1440.0 + second/86400.0
    # Start Meeus algorithm (variables are in his notation)
    if (month < 3):
        month = month + 12
        year = year - 1
    A = int(year/100)
    jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + \
         day - 1524.5
    # optionally adjust the jd for the switch from 
    # the Julian to Gregorian Calendar
    # here assumed to have occurred the day after 1582 October 4
    if calendar in ['standard','gregorian']:
        if jd >= 2299170.5:
            # 1582 October 15 (Gregorian Calendar)
            B = 2 - A + int(A/4)
        elif jd < 2299160.5:
            # 1582 October 5 (Julian Calendar)
            B = 0
        else:
            raise ValueError('impossible date (falls in gap between end of Julian calendar and beginning of Gregorian calendar')
    elif calendar == 'proleptic_gregorian':
        B = 2 - A + int(A/4)
    elif calendar == 'julian':
        B = 0
    else:
        raise ValueError('unknown calendar, must be one of julian,standard,gregorian,proleptic_gregorian, got %s' % calendar)
    # adjust for Julian calendar if necessary
    jd = jd + B
    return jd 

def epem(date):
    """
    input: date - datetime object (assumed UTC)
    ouput: gha - Greenwich hour angle, the angle between the Greenwich
           meridian and the meridian containing the subsolar point.
           dec - solar declination.
    """
    dg2rad = np.pi/180.
    rad2dg = 1./dg2rad
    # compute julian day from UTC datetime object.
    # datetime objects use proleptic gregorian calendar.
    jday = JulianDayFromDate(date,calendar='proleptic_gregorian')
    jd = np.floor(jday) # truncate to integer.
    # utc hour.
    ut = date.hour + date.minute/60. + date.second/3600.
    # calculate number of centuries from J2000
    t = (jd + (ut/24.) - 2451545.0) / 36525.
    # mean longitude corrected for aberration
    l = (280.460 + 36000.770 * t) % 360
    # mean anomaly
    g = 357.528 + 35999.050 * t
    # ecliptic longitude
    lm = l + 1.915 * np.sin(g*dg2rad) + 0.020 * np.sin(2*g*dg2rad)
    # obliquity of the ecliptic
    ep = 23.4393 - 0.01300 * t
    # equation of time
    eqtime = -1.915*np.sin(g*dg2rad) - 0.020*np.sin(2*g*dg2rad) \
            + 2.466*np.sin(2*lm*dg2rad) - 0.053*np.sin(4*lm*dg2rad)
    # Greenwich hour angle
    gha = 15*ut - 180 + eqtime
    # declination of sun
    dec = np.arcsin(np.sin(ep*dg2rad) * np.sin(lm*dg2rad)) * rad2dg
    return gha, dec

def daynight_terminator(date, delta, lonmin, lonmax):
    """
    date is datetime object (assumed UTC).
    nlons is # of longitudes used to compute terminator."""
    dg2rad = np.pi/180.
    lons = np.arange(lonmin,lonmax+0.5*delta,delta,dtype=np.float32)
    # compute greenwich hour angle and solar declination
    # from datetime object (assumed UTC).
    tau, dec = epem(date)
    # compute day/night terminator from hour angle, declination.
    longitude = lons + tau
    lats = np.arctan(-np.cos(longitude*dg2rad)/np.tan(dec*dg2rad))/dg2rad
    return lats, lons


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
        #current_datetime = dt.now(utc_timezone)
        #jd, fr = jday(current_dt.year, current_dt.month, current_dt.day, current_dt.hour, current_dt.minute, current_dt.second)
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

    def propagate_sgp4_timeseries(self, timestamps):
        lats = []
        lons = []
        for time in timestamps:
            lat, lon = self.propagate_sgp4(time)
            lats.append(lat)
            lons.append(lon)
        return lats, lons

def get_gs_access_circle(center_lat, center_lon, radius_km=10):
    circle_coordinates = []

    # Generate 360 points around the center
    for angle_deg in range(361):
        dx = radius_km * 1000 * np.cos(np.deg2rad(angle_deg))  # Convert km to meters
        dy = radius_km * 1000 * np.sin(np.deg2rad(angle_deg))  # Convert km to meters
        new_lat = center_lat + (dy / 111320)
        new_lon = center_lon + (dx / (111320 * np.cos(np.deg2rad(center_lat))))
        circle_coordinates.append((new_lat,new_lon))
    return circle_coordinates

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
    ], style={'width': '100vw', 'height': '100vh', 'margin': '0px', 'padding': '0px'})
)


@app.callback(Output('live-update-text', 'children'),
              Input('interval-component', 'n_intervals'))
def update_metrics(n):
    current_datetime = dt.now(utc_timezone)
    lat, lon = sat.propagate_sgp4(current_datetime)
    style = {'padding': '-5px', 'fontSize': '16px'}
    return [
        html.Span(f'Timestamp: {current_datetime.strftime("%Y-%m-%d %H:%M:%S")} UTC ', style=style),
        html.Span(f'Longitude: {lon:.3f}° ', style=style),
        html.Span(f'Latitude: {lat:.3f}° ', style=style),
    ]

# Multiple components can update every time interval gets fired.

# Collect some data
@app.callback(Output('live-update-graph', 'figure'),
              Input('interval-component', 'n_intervals'))
def update_graph_live(n):
    global data
    global forward_data
    current_datetime = dt.now(utc_timezone)
    lat, lon = sat.propagate_sgp4(current_datetime)
    #data = pd.DataFrame(columns=['Timestamp', 'Latitude', 'Longitude'])
    #data = data.append({'Timestamp': current_datetime, 'Latitude': lat, 'Longitude': lon}, ignore_index=True)
    data = pd.concat([data, pd.DataFrame.from_dict({'Timestamp': [current_datetime.strftime("%Y-%m-%d %H:%M:%S")], 'Latitude': [lat], 'Longitude': [lon],
                                                    'Timestamp_LA': [current_datetime.astimezone(la_tz).strftime("%Y-%m-%d %H:%M:%S")]})], ignore_index=True)
    data = data.tail(tail_size)
    if current_datetime > dt.strptime(forward_data['Timestamp'].iloc[0],"%Y-%m-%d %H:%M:%S").replace(tzinfo=utc_timezone):
        latest_future_datetime = current_datetime + timedelta(minutes=end_time_minutes)
        last_forward_lat, last_forward_lon = sat.propagate_sgp4(latest_future_datetime)
        forward_data = forward_data.iloc[1:]
        forward_data = pd.concat([forward_data, pd.DataFrame.from_dict({'Timestamp': [latest_future_datetime.strftime("%Y-%m-%d %H:%M:%S")], 'Latitude': [last_forward_lat], 'Longitude': [last_forward_lon], 'Color': ['orange'], 
                                                                        'Timestamp_LA': [latest_future_datetime.astimezone(la_tz).strftime("%Y-%m-%d %H:%M:%S")]})], ignore_index=True)

    # Create a list of dictionaries for the ground station information
    ground_stations = [
        #{'name': 'Bahrain', 'lat': 26.05138, 'lon': 50.50307, 'color': "red"},
        #{'name': 'Dubbo', 'lat': -32.17723, 'lon': 148.615407, 'color': "green"},
        {'name': 'Ireland', 'lat': 53.40667, 'lon': -6.22438, 'color': "#15b30c"},
        {'name': 'CapeTown', 'lat': -34.02762, 'lon': 18.71868, 'color': "cyan"},
        {'name': 'Punta Arenas', 'lat': -52.94111, 'lon': -70.84972, 'color': "magenta"},
    ]

    #fig = px.line_geo(data, lat="Latitude", lon="Longitude", markers="Color", hover_data="Timestamp")
    #fig= px.line_geo(data, lat="Latitude", lon="Longitude", hover_data={"Timestamp": True})
    #fig.update_traces(mode='lines+markers', marker=dict(size=8, color='green'), line=dict(width=0.5, color='black'))

    fig = go.Figure()


    #fig = px.line_geo(lat=data["Latitude"], lon=data["Longitude"], color=data["Color"], hover_data="Timestamp", markers=True)
    #fig.add_trace(px.line_geo(forward_data, lat="Latitude", lon="Longitude", hover_data="Timestamp", markers=True)
    # fig.add_trace(
    #     go.Scattergeo(
    #         lat=[data["Latitude"].iloc[-1]], 
    #         lon=[data["Longitude"].iloc[-1]],
    #         mode = 'markers',
    #         hovertemplate = [f'Latitude: {data["Latitude"].iloc[-1]:.3f}°<br>Longitude: {data["Longitude"].iloc[-1]:.3f}°<br>Timestamp: {current_dt.strftime("%Y-%m-%d %H:%M:%S")} UTC'],
    #         marker_symbol='star-diamond',
    #         name='Droid1',
    #         #hover_text=[data["Timestamp"].iloc[-1]],
    #         marker=dict(
    #             size=15,
    #             color='red',
    #             line=dict(width=0.5, color='black'),
    #             opacity=1
    #         ),
    #     )
    # )

    # Add Scattermapbox trace for the last data point
    # Add Scattermapbox trace for the last data point
    # fig.add_trace(
    #     go.Scattermapbox(
    #         lat=[data["Latitude"].iloc[-1]],  # Latitude of the last data point
    #         lon=[data["Longitude"].iloc[-1]],  # Longitude of the last data point
    #         mode='markers',
    #         hovertemplate=f'Latitude: {{lat}}°<br>Longitude: {{lon}}°<br>Timestamp: {current_dt.strftime("%Y-%m-%d %H:%M:%S")} UTC',
    #         marker_symbol='star-diamond',
    #     )
    # )


    # # Add a single star marker for the last entry
    # fig.add_trace(
    #     go.Scattermapbox(
    #         lat=[data["Latitude"].iloc[-1]],  # Provide the latitude of the last entry
    #         lon=[data["Longitude"].iloc[-1]],  # Provide the longitude of the last entry
    #         hovertext=[data["Timestamp"].iloc[-1]],  # Set hovertext to the "Timestamp" of the last entry
    #         mode="markers",
    #         marker=dict(size=10, color='red', symbol='star'),  # Set symbol to 'star'
    #         name='Last Entry'
    #     )
    # )
    # # Add Scattergeo trace for the last data point
    # fig.add_trace(
    #     go.Scattergeo(
    #         lat=[data["Latitude"].iloc[-1]],
    #         lon=[data["Longitude"].iloc[-1]],
    #         mode='markers',
    #         hovertemplate=f'Latitude: {data["Latitude"].iloc[-1]:.3f}°<br>Longitude: {data["Longitude"].iloc[-1]:.3f}°<br>Timestamp: {current_dt.strftime("%Y-%m-%d %H:%M:%S")} UTC',
    #         marker_symbol='star-diamond',
    #         name='Droid1',
    #         marker=dict(
    #             size=15,
    #             color='red',
    #             line=dict(width=0.5, color='black'),
    #             opacity=1
    #         ),
    #     )
    # )
    
    # fig.update_geos(lataxis_showgrid=True, lonaxis_showgrid=True)
    # fig.update_geos(
    #     visible=False, resolution=50,
    #     showcountries=True, countrycolor="RebeccaPurple",
    #     showland=True, landcolor="LightGreen",
    #     showocean=True, oceancolor="LightBlue",
    #     showlakes=True, lakecolor="#3399FF",
    #     showcoastlines=True, coastlinecolor="RebeccaPurple",
    #     showsubunits=True, subunitcolor="Blue",
    #     uirevision=1,
    # )

    # Ground Station Plotting

        # Day Night Shading
    day_night_lats, day_night_lons = daynight_terminator(current_datetime, 1, -180, 180)

    # Convert lons and lats to lists and append new elements
    day_night_lons = list(day_night_lons) + [day_night_lons[-1], day_night_lons[0]]
    day_night_lats = list(day_night_lats) + [-90, -90]

    # Create a scatter plot of the terminator coordinates
    scatter = fig.add_trace(go.Scattermapbox(lon = day_night_lons,
                           lat = day_night_lats,
                           mode = 'lines', 
                           line = dict(width = 2, color = 'rgba(138,138,138,1)'),
                           fill = 'toself',
                           hoverinfo='skip',
                           fillcolor = 'rgba(138,138,138,0.2)',
                           hovertemplate=None))

    for station in ground_stations:
        circle_coordinates = get_gs_access_circle(station['lat'],station['lon'],gs_comms_radius)
        # Create the trace for the circle
        circle_trace = fig.add_trace(go.Scattermapbox(
            mode='lines',
            lat=[coord[0] for coord in circle_coordinates],
            lon=[coord[1] for coord in circle_coordinates],
            hovertext=[station['name']] * len(circle_coordinates),  # Set hover text to station['name']
            line=dict(color=station['color'], width=2),
            fill='toself',
            hovertemplate=
                "<b>Ground Station:</b> %{hovertext}<br>" +
                "<b>Latitude:</b> %{lat:.3f}°<br>" +
                "<b>Longitude:</b> %{lon:.3f}°<br>" +
                "<extra></extra>",
        ))


    # Create a color list of the same length as your data, with all 'green'
    color_list = ['blue' for _ in range(len(data['Latitude']))]

    # Change the last color to 'red'
    color_list[-1] = 'red'
    #print(forward_data["Timestamp_LA"].iloc[0].astimezone(la_tz).strftime("%Y-%m-%d %H:%M:%S"))
    fig.add_trace(
        go.Scattermapbox(
            lat=forward_data["Latitude"],  # Provide the latitude data
            lon=forward_data["Longitude"],  # Provide the longitude data
            hovertext=forward_data["Timestamp"],  # Set hovertext to "Timestamp" column
            customdata=[dt.strptime(i,"%Y-%m-%d %H:%M:%S").astimezone(la_tz).strftime("%Y-%m-%d %H:%M:%S") for i in forward_data["Timestamp_LA"]],
            hovertemplate=
                "<b>Timestamp [UTC]:</b> %{hovertext}<br>" +
                "<b>Timestamp    [PT]:</b> %{customdata}<br>" +
                "<b>Latitude:</b> %{lat:.3f}°<br>" +
                "<b>Longitude:</b> %{lon:.3f}°<br>" +
                "<extra></extra>",
            mode="lines+markers",
            marker=dict(size=8, color=forward_data['Color']),
            line=dict(width=0.5, color='black'),
            name='Droid1'
        )
    )

    # Add Scattermapbox trace
    fig.add_trace(
        go.Scattermapbox(
            lat=data["Latitude"],  # Provide the latitude data
            lon=data["Longitude"],  # Provide the longitude data
            hovertext=data["Timestamp"],  # Set hovertext to "Timestamp" column
            customdata=[dt.strptime(i,"%Y-%m-%d %H:%M:%S").astimezone(la_tz).strftime("%Y-%m-%d %H:%M:%S") for i in data["Timestamp_LA"]],
            hovertemplate=
                "<b>Timestamp [UTC]:</b> %{hovertext}<br>" +
                "<b>Timestamp    [PT]:</b> %{customdata}<br>" +
                "<b>Latitude:</b> %{lat:.3f}°<br>" +
                "<b>Longitude:</b> %{lon:.3f}°<br>" +
                "<extra></extra>",
            mode="lines+markers",
            marker=dict(size=8, color=color_list),
            line=dict(width=0.5, color='black'),
            name='Droid1',
        )
    )

    # Add ground stations as scatter markers
    #for station in ground_stations:
        # size_adjusted = 10000 / 111.32  # Approximate value for 1 degree of latitude
        # size_adjusted *= np.cos(np.deg2rad(station['lat']))
        # fig.add_trace(go.Scattergeo(
        #     lat=[station['lat']],
        #     lon=[station['lon']],
        #     mode='markers',
        #     marker=dict(
        #         size=size_adjusted,
        #         sizemode='area',  # Adjusts marker size based on zoom level
        #         color=station['color'],
        #         opacity=0.7
        #     ),
        #     name=station['name'],
        #     hovertemplate=f'Name: {station["name"]}<br>Latitude: {station["lat"]:.3f}°<br>Longitude: {station["lon"]:.3f}°<extra></extra>',
        # ))

    fig.update_layout(
        width=1900,  # Width of the figure (adjust to your screen resolution)
        height=1080,  # Height of the figure (adjust to your screen resolution)
        uirevision=1
    )
       # Set layout for the map
    fig.update_layout(
        mapbox=dict(
            #style="open-street-map",  # Choose the map style (e.g., "carto-positron", "stamen-terrain", etc.)
            #style="carto-positron",  # Choose the map style (e.g., "carto-positron", "stamen-terrain", etc.)
            #style="carto-darkmatter",  # Choose the map style (e.g., "carto-positron", "stamen-terrain", etc.)
            #style="stamen-terrain",  # Choose the map style (e.g., "carto-positron", "stamen-terrain", etc.)
            #style="stamen-toner",  # Choose the map style (e.g., "carto-positron", "stamen-terrain", etc.)
            style="stamen-watercolor",  # Choose the map style (e.g., "carto-positron", "stamen-terrain", etc.)
            center=dict(lat=0, lon=0),  # Set the initial center of the map
            zoom=1  # Set the initial zoom level
        ),
        showlegend=False  # Hide the legend
    )
    return fig


if __name__ == '__main__':
    sat = SpacecraftPropagator()
    data = {
        'Timestamp': [],
        'Latitude' : [],
        'Longitude': [],
    }

    forward_data = {
        'Timestamp': [],
        'Latitude' : [],
        'Longitude': [],
        'Color': []
    }

    data = pd.DataFrame.from_dict(data)
    forward_data = pd.DataFrame.from_dict(forward_data)
    current_datetime = dt.now(utc_timezone)
    
    interval_s = 30  # interval in seconds
    end_time_minutes = 95  # minutes into the future
    # Create a list of datetimes using list comprehension
    end_time = current_datetime + timedelta(minutes=end_time_minutes)
    future_datetimes = [(current_datetime + timedelta(seconds=i)) for i in range(0, int((end_time-current_datetime).total_seconds()), interval_s)]
    future_datetimes_str = [i.strftime("%Y-%m-%d %H:%M:%S") for i in future_datetimes]
    #future_datetimes = [current_time.strftime("%Y-%m-%d %H:%M:%S") for current_time in (current_datetime + timedelta(minutes=i) for i in range(0, (end_time - current_datetime).seconds, interval.seconds))]
    lats, lons = sat.propagate_sgp4_timeseries(future_datetimes)
    forward_data = pd.concat([forward_data, pd.DataFrame.from_dict({'Timestamp': future_datetimes_str, 'Latitude': lats, 'Longitude': lons, 'Color': ['orange']*len(future_datetimes)})], ignore_index=True)
    forward_data['Timestamp_LA'] = [i.astimezone(la_tz).strftime("%Y-%m-%d %H:%M:%S") for i in future_datetimes]


    tail_size = 50
    
    #https://file.scirp.org/pdf/IJCNS20110900004_79507510.pdf
    gs_comms_radius = 2400/2

    app.run_server(debug=True)











