#
# Python script for plotting meshes 
#
# Plotting the mesh using plotly 
#
import plotly.graph_objs as go
import numpy as np 
import pandas as pd

#
#
#
########################################################
#
#   
#
def plot_faults(gdf,quake_lat,quake_long,quake_name):
    # Extract the relevant latitudes, longitudes, and origin_fault 
    latitudes = []
    longitudes = []
    fault_names = []

    for _, row in gdf.iterrows():
        latitudes.append(row['lat1'])
        longitudes.append(row['lon1'])
        fault_names.append(row['origin_fault'])
        
    # Create a DataFrame associated EFSM names
    points_df = pd.DataFrame({
        'Latitude': latitudes,
        'Longitude': longitudes,
        'Fault Name': fault_names
        })


# Create a DataFrame earthquake
    quake_df = pd.DataFrame({
        'Latitude': quake_lat,
        'Longitude': quake_long,
        'Fault Name':quake_name
        })
    
    
# Create the main figure with fault nodes using go.Scattermapbox
    fig = go.Figure()

# Add fault nodes trace
    fig.add_trace(go.Scattermapbox(
        lat=points_df['Latitude'],
        lon=points_df['Longitude'],
        text=points_df['Fault Name'],
        mode='markers',
        marker=go.scattermapbox.Marker(size=10, color='blue'),  # Customize as needed
        name='Fault Nodes'
    ))



# Create the trace for earthquakes with large yellow stars
    quake_trace = go.Scattermapbox(
        lat=quake_df['Latitude'],
        lon=quake_df['Longitude'],
        text=quake_df['Fault Name'],
        mode='markers',
        marker=go.scattermapbox.Marker(size=15, color='red'),
        name='Earthquakes'
    )

# Add the earthquake trace to the figure
    fig.add_trace(quake_trace)

# Set the layout for the map
    fig.update_layout(
        mapbox_style="carto-positron",  # Choose the map style
        mapbox=dict(
            zoom=3,  # Adjust zoom level
            center=go.layout.mapbox.Center(lat=quake_df['Latitude'].mean(), lon=quake_df['Longitude'].mean())  # Center the map
        ),
        title="Picking EFSM fault that is closest to earthquake Locations"
    )

# Set the dimensions for the figure and plot
    # fig.update_layout(
    #     autosize=False,
    #     width=900,
    #     height=600)

    fig.show()
#
#
#
########################################################
#
#   This function uses plotly to plot a wireframe mesh with nodes defined in 'highlight' displayed in red
#
def plot_highlight_nodes(tri,nodes,highlight,name):
    X = nodes[:,0]  
    Y = nodes[:,1]  
    Z = nodes[:,2]  
    
# Extract edges from the mesh
    edges = []
    for i in range(tri.shape[0]):
        edges.extend([(tri[i, 0], tri[i, 1]), (tri[i, 1], tri[i, 2]), (tri[i, 2], tri[i, 0])])

# Remove duplicate edges
    unique_edges = list(set(edges))

# Extract coordinates for edges
    edge_x = []
    edge_y = []
    edge_z = []
    for edge in unique_edges:
        edge_x.extend([X[edge[0]], X[edge[1]], None])
        edge_y.extend([Y[edge[0]], Y[edge[1]], None])
        edge_z.extend([Z[edge[0]], Z[edge[1]], None])

# Create a Plotly scatter3d trace for wireframe
    trace_edges = go.Scatter3d(
        x=edge_x,
        y=edge_y,
        z=edge_z,
        mode='lines',
        line=dict(color='black', width=1),
        name = 'mesh'
    )

 # Create a Plotly  scatter3d for highlighted nodes   
    trace_highlight = go.Scatter3d(
        x=X[highlight==1],
        y=Y[highlight==1],
        z=Z[highlight==1],
        mode='markers',
        marker=dict(size=2.5,
                    line=dict(color='DarkSlateGrey', width=1)),
        name = name
    )
# Create a Plotly mesh3d trace for the triangular mesh with color based on velocity



# Create Plotly figure
    fig = go.Figure(data=[trace_edges, trace_highlight])

# Update layout and show plot
    fig.update_layout(scene=dict(aspectmode='data'))
    fig.show()
    
#
#
#
########################################################
#
#   
#
    



#
#  Plot fault mesh and location of earthquake to check earthquake is within fault plane
#
def plot_mesh_quake(tri,nodes,quake_x,quake_y,quake_z):
    X = nodes[:,0]  
    Y = nodes[:,1]  
    Z = nodes[:,2]  
    
# Extract edges from the mesh
    edges = []
    for i in range(tri.shape[0]):
        edges.extend([(tri[i, 0], tri[i, 1]), (tri[i, 1], tri[i, 2]), (tri[i, 2], tri[i, 0])])

# Remove duplicate edges
    unique_edges = list(set(edges))

# Extract coordinates for edges
    edge_x = []
    edge_y = []
    edge_z = []
    for edge in unique_edges:
        edge_x.extend([X[edge[0]], X[edge[1]], None])
        edge_y.extend([Y[edge[0]], Y[edge[1]], None])
        edge_z.extend([Z[edge[0]], Z[edge[1]], None])

# Create a Plotly scatter3d trace for wireframe
    trace = go.Scatter3d(
        x=edge_x,
        y=edge_y,
        z=edge_z,
        mode='lines',
        line=dict(color='black', width=2),
        name='Fault Mesh'
    )
# Create a Plotly scatter3d trace for earthquake epicentre
    quake_trace = go.Scatter3d(
        x= [quake_x],
        y=[quake_y],
        z=[quake_z],
        mode='markers',
        name='Earthquake Epicentre',
        )
    
    
# Create Plotly figure
    fig = go.Figure(data=[trace,quake_trace])

# Update layout and show plot
    fig.update_layout(scene=dict(aspectmode='data'))
    fig.show()    

#
#
#
########################################################
#
#   
#
def plot_mesh_node(tri,nodes,node_var,name):
    X = nodes[:,0]  
    Y = nodes[:,1]  
    Z = nodes[:,2]  
    
# Extract edges from the mesh
    edges = []
    for i in range(tri.shape[0]):
        edges.extend([(tri[i, 0], tri[i, 1]), (tri[i, 1], tri[i, 2]), (tri[i, 2], tri[i, 0])])

# Remove duplicate edges
    unique_edges = list(set(edges))

# Extract coordinates for edges
    edge_x = []
    edge_y = []
    edge_z = []
    for edge in unique_edges:
        edge_x.extend([X[edge[0]], X[edge[1]], None])
        edge_y.extend([Y[edge[0]], Y[edge[1]], None])
        edge_z.extend([Z[edge[0]], Z[edge[1]], None])

# Create a Plotly scatter3d trace for wireframe
    trace_edges = go.Scatter3d(
        x=edge_x,
        y=edge_y,
        z=edge_z,
        mode='lines',
        line=dict(color='black', width=1)
    )

# Create a Plotly mesh3d trace for the triangular mesh with color based on velocity
    trace_mesh = go.Mesh3d(
        x=X,
        y=Y,
        z=Z,
        i=tri[:, 0],
        j=tri[:, 1],
        k=tri[:, 2],
        intensity=node_var,
        colorscale='Viridis',
        colorbar=dict(title=name),
        showscale=True
    )

# Create Plotly figure
    fig = go.Figure(data=[trace_edges, trace_mesh])

# Update layout and show plot
    fig.update_layout(scene=dict(aspectmode='data'))
    fig.show()



#
#
########################################################
#
#    
def plot_mesh_cell(tri,nodes,cell_var,name):
#  plot_mesh_cell: plot mesh with colour in each cell defined by cell_var
#  
#
    X = nodes[:,0]  
    Y = nodes[:,1]  
    Z = nodes[:,2]  
    
# Extract edges from the mesh
    edges = []
    for i in range(tri.shape[0]):
        edges.extend([(tri[i, 0], tri[i, 1]), (tri[i, 1], tri[i, 2]), (tri[i, 2], tri[i, 0])])

# Remove duplicate edges
    unique_edges = list(set(edges))

# Extract coordinates for edges
    edge_x = []
    edge_y = []
    edge_z = []
    for edge in unique_edges:
        edge_x.extend([X[edge[0]], X[edge[1]], None])
        edge_y.extend([Y[edge[0]], Y[edge[1]], None])
        edge_z.extend([Z[edge[0]], Z[edge[1]], None])

# Create a Plotly scatter3d trace for wireframe
    trace_edges = go.Scatter3d(
        x=edge_x,
        y=edge_y,
        z=edge_z,
        mode='lines',
        line=dict(color='black', width=1)
    )
    var = np.zeros(len(X))   # colour is placed on vertices which is interpolated 
    for i in range(len(tri)):
        c_id = tri[i]
        var[c_id] = cell_var[i]

# Create a Plotly mesh3d trace for the triangular mesh with color based on velocity
    trace_mesh = go.Mesh3d(
        x=X,
        y=Y,
        z=Z,
        i=tri[:, 0],
        j=tri[:, 1],
        k=tri[:, 2],
        intensity=var,
        colorscale='solar',
        colorbar=dict(title=name),
        showscale=True
    )

# Create Plotly figure
    # fig = go.Figure(data=[trace_mesh])
    fig = go.Figure(data=[trace_edges, trace_mesh])

# Update layout and show plot
    fig.update_layout(scene=dict(aspectmode='data'))
    fig.show()
#
#
########################################################
#
#    
def compare_meshes(tri1,nodes1,tri2,nodes2):
#   graphically compare two meshes 
#       tri1:    no. of cellx x 3 array with connectivity between cells and nodes for first mesh
#       nodes1:  no. of nodes x 3 array with x,y,z for each node in first mesh
#
#
#       mesh1 is plotted in red and mesh2 is plotted in black 

# Extract edges from first mesh
    X1 = nodes1[:,0]  
    Y1 = nodes1[:,1]  
    Z1 = nodes1[:,2]  
  
    edges1 = []
    for i in range(tri1.shape[0]):
        edges1.extend([(tri1[i, 0], tri1[i, 1]), (tri1[i, 1], tri1[i, 2]), (tri1[i, 2], tri1[i, 0])])
# Extract edges from second mesh
    X2 = nodes2[:,0]  
    Y2 = nodes2[:,1]  
    Z2 = nodes2[:,2]  

    edges2 = []
    for i in range(tri2.shape[0]):
        edges2.extend([(tri2[i, 0], tri2[i, 1]), (tri2[i, 1], tri2[i, 2]), (tri2[i, 2], tri2[i, 0])])

# Remove duplicate edges
    unique_edges1 = list(set(edges1))
    unique_edges2 = list(set(edges2))

# Extract coordinates for edges
    edge1_x = []
    edge1_y = []
    edge1_z = []
    for edge1 in unique_edges1:
        edge1_x.extend([X1[edge1[0]], X1[edge1[1]], None])
        edge1_y.extend([Y1[edge1[0]], Y1[edge1[1]], None])
        edge1_z.extend([Z1[edge1[0]], Z1[edge1[1]], None])

    edge2_x = []
    edge2_y = []
    edge2_z = []
    for edge2 in unique_edges2:
        edge2_x.extend([X2[edge2[0]], X2[edge2[1]], None])
        edge2_y.extend([Y2[edge2[0]], Y2[edge2[1]], None])
        edge2_z.extend([Z2[edge2[0]], Z2[edge2[1]], None])

# Create a Plotly scatter3d trace for wireframe
    trace1 = go.Scatter3d(
        x=edge1_x,
        y=edge1_y,
        z=edge1_z,
        mode='lines',
        line=dict(color='red', width=2),
        name='Mesh 1'
    )
    trace2 = go.Scatter3d(
        x=edge2_x,
        y=edge2_y,
        z=edge2_z,
        mode='lines',
        line=dict(color='black', width=1),
        name='Mesh 2'
    )


# Create Plotly figure
    fig = go.Figure(data=[trace1,trace2])

# Update layout and show plot
    fig.update_layout(scene=dict(aspectmode='data'))
    fig.show()        
#
#
########################################################
#
#
def plot_mesh(tri,nodes):
    
    
    X = nodes[:,0]  
    Y = nodes[:,1]  
    Z = nodes[:,2] 
    
# Extract edges from the mesh
    edges = []
    for i in range(tri.shape[0]):
        edges.extend([(tri[i, 0], tri[i, 1]), (tri[i, 1], tri[i, 2]), (tri[i, 2], tri[i, 0])])

# Remove duplicate edges
    unique_edges = list(set(edges))

# Extract coordinates for edges
    edge_x = []
    edge_y = []
    edge_z = []
    for edge in unique_edges:
        edge_x.extend([X[edge[0]], X[edge[1]], None])
        edge_y.extend([Y[edge[0]], Y[edge[1]], None])
        edge_z.extend([Z[edge[0]], Z[edge[1]], None])

# Create a Plotly scatter3d trace for wireframe
    trace = go.Scatter3d(
        x=edge_x,
        y=edge_y,
        z=edge_z,
        mode='lines',
        line=dict(color='black', width=2)
    )

# Create Plotly figure
    fig = go.Figure(data=[trace])

# Update layout and show plot
    fig.update_layout(scene=dict(aspectmode='data'))
    fig.show()    