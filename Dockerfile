# Use an official Python runtime as a parent image
FROM python:3.12-slim

# Set the working directory in the container
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    libglu1-mesa \
    libgl1-mesa-glx \
    libgl1-mesa-dri \
    freeglut3-dev \
    libx11-dev \
    libxext-dev \
    libxrender-dev \
    libxft-dev \
    libfreetype6-dev \
    libxcursor1 \
    libxinerama1 \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip install --no-cache-dir \
    jupyterlab \
    numpy \
    pandas \
    gmsh \
    geopandas \
    pyproj \
    plotly \
    geojson\
    meshio 


# Expose port 8888 for JupyterLab
EXPOSE 8888

# Set environment variable to avoid Jupyter asking for confirmation on launch
ENV JUPYTER_ENABLE_LAB=yes

# Copy the current directory contents into the container at /app
COPY . /app

# Run JupyterLab when the container launches
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--no-browser", "--allow-root"]
#CMD ["jupyter", "lab", "Xvfb :1 -screen 0 1024x768x24 & jupyter lab --ip=0.0.0.0 --no-browser --allow-root"]

