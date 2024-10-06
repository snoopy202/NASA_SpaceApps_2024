import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from astropy.table import Table
import astropy.units as u

# Load data
star_subset = 'star_subset.tab'  # Path to your data file
subset = Table.read(star_subset, format='ascii.tab').to_pandas()  # Convert to DataFrame

# Check if required columns exist
required_columns = ['sy_dist', 'pl_orbper', 'pl_rade', 'pl_name', 'pl_orbsmax', 'st_rad']
for col in required_columns:
    if col not in subset.columns:
        st.error(f"Missing required column: {col}")
        st.stop()

# Function to calculate SNR for each exoplanet
def calculate_snr(row, D, ES):
    R = row['st_rad']  # Stellar radius in solar radii
    RP = row['pl_rade']  # Planetary radius in Earth radii
    snr = SNR0 * ((R * RP * (D / 6)) / ((ES / 10) * PS)) ** 2
    return snr

# Constants
SNR0 = 100  # Initial SNR
PS = 1.0  # Instrument sensitivity in mJy

# Streamlit UI components
st.title("Exoplanet Visualization")
telescope_diameter = st.slider("Telescope Diameter (m)", min_value=1.0, max_value=100.0, value=10.0)
instrument_sensitivity = st.slider("Instrument Sensitivity (mJy)", min_value=100.0, max_value=5000.0, value=1000.0)

# Calculate SNR for each exoplanet
subset['SNR'] = subset.apply(calculate_snr, axis=1, D=telescope_diameter, ES=10)  # Assume ES = 10 pc for filtering
# Filter observable exoplanets
observable_exoplanets = subset[subset['SNR'] > 5]

# Check if there are observable exoplanets
if observable_exoplanets.empty:
    st.warning("No observable exoplanets found with SNR > 5.")
else:
    # Create 3D Visualization
    fig = go.Figure()

    for _, exoplanet in observable_exoplanets.iterrows():
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
    st.plotly_chart(fig)

# Show observable exoplanets
st.write("Observable Exoplanets (SNR > 5):")
st.dataframe(observable_exoplanets[['pl_name', 'SNR']])
