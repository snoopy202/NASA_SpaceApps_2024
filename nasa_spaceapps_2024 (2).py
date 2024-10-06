import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy.table import Column
from astropy.table import join
from astropy import table
from astropy.table import QTable
import astropy.units as u

import plotly.graph_objects as go
import ipywidgets as widgets
from IPython.display import display
#https://exoplanetarchive.ipac.caltech.edu/overview/24%20Boo

star_subset = '/content/star_subset.tab'
subset = QTable.read('star_subset.tab', format='ascii.tab')
subset

x = 1.0 * u.parsec
x

subset['pl_orbper'].u = ('d')
subset['pl_orbsmax'].u = ('AU')
subset

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from astropy.table import Table

# Load data
star_subset = '/content/star_subset.tab'
subset = Table.read(star_subset, format='ascii.tab').to_pandas()  # Convert to DataFrame

# Check if required columns exist
required_columns = ['sy_dist', 'pl_orbper', 'pl_rade', 'pl_name', 'pl_orbsmax']
for col in required_columns:
    if col not in subset.columns:
        raise ValueError(f"Missing required column: {col}")

# Create 3D Visualization
def create_3d_visualization(exoplanets):
    fig = go.Figure()

    # Iterate through the DataFrame rows
    for _, exoplanet in exoplanets.iterrows():
        # Calculate coordinates within the loop
        x = exoplanet['sy_dist'] * np.cos(np.radians(exoplanet['pl_orbper']))  # X coordinate
        y = exoplanet['pl_orbper'] * np.sin(np.radians(exoplanet['pl_orbper']))  # Y coordinate
        z = exoplanet['pl_rade']  # Z coordinate

        fig.add_trace(go.Scatter3d(
            x=[x],  # X coordinate
            y=[y],  # Y coordinate
            z=[z],  # Z coordinate
            mode='markers',
            marker=dict(size=5),
            name=exoplanet['pl_name']  # Label with exoplanet name
        ))

    # Update layout for better visualization
    fig.update_layout(title='Exoplanet Observational Paths',
                      scene=dict(xaxis_title='sy_dist',
                                 yaxis_title='pl_orbper',
                                 zaxis_title='pl_rade'))
    fig.show()

# Check for empty DataFrame
if subset.empty:
    print("No data available to visualize.")
else:
    # Call the function with the subset DataFrame
    create_3d_visualization(subset)

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from astropy.table import Table
import ipywidgets as widgets
from IPython.display import display, clear_output  # Import clear_output here

# Load data
star_subset = '/content/star_subset.tab'  # Path to your data file
subset = Table.read(star_subset, format='ascii.tab').to_pandas()  # Convert to DataFrame

# Check if required columns exist
required_columns = ['sy_dist', 'pl_orbper', 'pl_rade', 'pl_name', 'pl_orbsmax']
for col in required_columns:
    if col not in subset.columns:
        raise ValueError(f"Missing required column: {col}")

# Function to create 3D visualization
def create_3d_visualization(telescope_diameter, instrument_sensitivity):
    print(f"Telescope Diameter: {telescope_diameter} m")
    print(f"Instrument Sensitivity: {instrument_sensitivity} mJy")

    # Print a preview of relevant columns
    print("Exoplanets Data Preview:")
    print(subset[['pl_name', 'sy_dist', 'pl_rade', 'pl_orbsmax']].describe())

    # Filter based on user-defined parameters
    filtered_exoplanets = subset[
        (subset['pl_rade'] <= telescope_diameter) &
        (subset['sy_dist'] <= instrument_sensitivity)
    ]

    # Check if there are any exoplanets after filtering
    if filtered_exoplanets.empty:
        print("No exoplanets found with the current parameters.")
        return

    # Create 3D visualization
    fig = go.Figure()

    for _, exoplanet in filtered_exoplanets.iterrows():
        # Calculate coordinates
        x = exoplanet['sy_dist'] * np.cos(np.radians(exoplanet['pl_orbper']))  # X coordinate
        y = exoplanet['pl_orbper'] * np.sin(np.radians(exoplanet['pl_orbper']))  # Y coordinate
        z = exoplanet['pl_rade']  # Z coordinate

        fig.add_trace(go.Scatter3d(
            x=[x],  # X coordinate
            y=[y],  # Y coordinate
            z=[z],  # Z coordinate
            mode='markers',
            marker=dict(size=5),
            name=exoplanet['pl_name']  # Label with exoplanet name
        ))

    # Update layout for better visualization
    fig.update_layout(title='Exoplanet Observational Paths',
                      scene=dict(xaxis_title='Distance (ly)',
                                 yaxis_title='Orbital Period (days)',
                                 zaxis_title='Planet Radius (R⊕)'))
    fig.show()

# Interactive widgets for user input
telescope_diameter_widget = widgets.FloatSlider(
    value=10.0,
    min=1.0,
    max=100.0,
    step=1.0,
    description='Telescope Diameter (m):',
    continuous_update=False
)

instrument_sensitivity_widget = widgets.FloatSlider(
    value=1000.0,
    min=100.0,
    max=5000.0,
    step=100.0,
    description='Instrument Sensitivity (mJy):',
    continuous_update=False
)

ui = widgets.VBox([telescope_diameter_widget, instrument_sensitivity_widget])

def on_button_click(b):
    clear_output(wait=True)
    display(ui)
    create_3d_visualization(telescope_diameter_widget.value, instrument_sensitivity_widget.value)

# Button to trigger visualization
button = widgets.Button(description="Visualize Exoplanets")
button.on_click(on_button_click)
display(ui, button)

# Check if required columns exist
required_columns = ['sy_dist', 'pl_orbper', 'pl_rade', 'pl_name', 'pl_orbsmax']
for col in required_columns:
    if col not in subset.columns:
        raise ValueError(f"Missing required column: {col}")

# Function to create 3D visualization
def create_3d_visualization(telescope_diameter, instrument_sensitivity):
    print(f"Telescope Diameter: {telescope_diameter} m")
    print(f"Instrument Sensitivity: {instrument_sensitivity} mJy")

    # Print a preview of relevant columns
    print("Exoplanets Data Preview:")
    print(subset[['pl_name', 'sy_dist', 'pl_rade', 'pl_orbsmax']].describe())

    # Filter based on user-defined parameters
    filtered_exoplanets = subset[
        (subset['pl_rade'] <= telescope_diameter) &
        (subset['sy_dist'] <= instrument_sensitivity)
    ]

    # Check if there are any exoplanets after filtering
    if filtered_exoplanets.empty:
        print("No exoplanets found with the current parameters.")
        return

    # Create 3D visualization
    fig = go.Figure()

    for _, exoplanet in filtered_exoplanets.iterrows():
        # Calculate coordinates
        x = exoplanet['sy_dist'] * np.cos(np.radians(exoplanet['pl_orbper']))  # X coordinate
        y = exoplanet['pl_orbper'] * np.sin(np.radians(exoplanet['pl_orbper']))  # Y coordinate
        z = exoplanet['pl_rade']  # Z coordinate

        fig.add_trace(go.Scatter3d(
            x=[x],  # X coordinate
            y=[y],  # Y coordinate
            z=[z],  # Z coordinate
            mode='markers',
            marker=dict(size=5),
            name=exoplanet['pl_name']  # Label with exoplanet name
        ))

    # Update layout for better visualization
    fig.update_layout(title='Exoplanet Observational Paths',
                      scene=dict(xaxis_title='Distance (ly)',
                                 yaxis_title='Orbital Period (days)',
                                 zaxis_title='Planet Radius (R⊕)'))
    fig.show()

# Interactive widgets for user input
telescope_diameter_widget = widgets.FloatSlider(
    value=10.0,
    min=1.0,
    max=100.0,
    step=1.0,
    description='Telescope Diameter (m):',
    continuous_update=False
)

instrument_sensitivity_widget = widgets.FloatSlider(
    value=1000.0,
    min=100.0,
    max=5000.0,
    step=100.0,
    description='Instrument Sensitivity (mJy):',
    continuous_update=False
)

ui = widgets.VBox([telescope_diameter_widget, instrument_sensitivity_widget])

def on_button_click(b):
    clear_output(wait=True)
    display(ui)
    create_3d_visualization(telescope_diameter_widget.value, instrument_sensitivity_widget.value)

# Button to trigger visualization
button = widgets.Button(description="Visualize Exoplanets")
button.on_click(on_button_click)
display(ui, button)

D = 6  # Telescope diameter in meters
SNR0 = 100  # Initial SNR
ES = 10  # Distance in parsecs
PS = 1.0  # Instrument sensitivity in mJy (example value)

# Add a new column for SNR calculation
def calculate_snr(row):
    R = row['st_rad']  # Stellar radius in solar radii
    RP = row['pl_rade']  # Planetary radius in Earth radii

    # Calculate SNR using the formula
    snr = SNR0 * ((R * RP * (D / 6)) / ((ES / 10) * PS)) ** 2
    return snr

# Apply the SNR calculation
subset['SNR'] = subset.apply(calculate_snr, axis=1)

# Filter exoplanets with SNR > 5
observable_exoplanets = subset[subset['SNR'] > 5]

# Display the results
print(observable_exoplanets[['pl_name', 'SNR']])

import pandas as pd
import numpy as np
import ipywidgets as widgets
from astropy.table import Table
from IPython.display import display, clear_output

# Load the exoplanet dataset
star_subset = '/content/star_subset.tab'  # Path to your data file
subset = Table.read(star_subset, format='ascii.tab').to_pandas()  # Convert to DataFrame

# Assuming constant values
SNR0 = 100  # Initial SNR
PS = 1.0  # Instrument sensitivity in mJy (example value)

# Function to calculate SNR for each exoplanet
def calculate_snr(row, D, ES):
    R = row['st_rad']  # Stellar radius in solar radii
    RP = row['pl_rade']  # Planetary radius in Earth radii
    snr = SNR0 * ((R * RP * (D / 6)) / ((ES / 10) * PS)) ** 2
    return snr

# Function to calculate maximum observable distance
def calculate_es_max(D):
    return 15 * (D / 6) / PS

# Widget for telescope diameter
telescope_diameter_widget = widgets.FloatSlider(
    value=6.0,
    min=5.0,
    max=15.0,
    step=1.0,
    description='Telescope Diameter (m):',
    continuous_update=False
)

# Function to update results based on user input
def update_results(b):
    clear_output(wait=True)
    display(telescope_diameter_widget)

    D = telescope_diameter_widget.value

    # Calculate SNR for each exoplanet
    subset['SNR'] = subset.apply(calculate_snr, axis=1, D=D, ES=10)  # Assume ES = 10 pc for filtering

    # Filter exoplanets with SNR > 5
    observable_exoplanets = subset[subset['SNR'] > 5]

    # Calculate maximum observable distance for the current telescope diameter
    es_max = calculate_es_max(D)

    # Find exoplanets within the observable distance
    characterizable_exoplanets = subset[subset['sy_dist'] <= es_max]

    print(f"Observable Exoplanets (SNR > 5) for D = {D} m:")
    print(observable_exoplanets[['pl_name', 'SNR']])

    print(f"\nCharacterizable Exoplanets within {es_max:.2f} pc:")
    print(characterizable_exoplanets[['pl_name', 'sy_dist']])

# Button to trigger results update
button = widgets.Button(description="Update Results")
button.on_click(update_results)

# Display the widgets
display(telescope_diameter_widget, button)