import plotly.graph_objects as go
import numpy as np
import scipy.io as sc

#loading in data from matlab simulations
matfile = sc.loadmat("25gcrun.mat")
imported_data = np.array(matfile["y"])

imported_data = imported_data[::20,:]

#format the data to be fed to the graph (voltages)
tbv = np.zeros((imported_data.shape[0], 125))
for i in range(125):
    tbv[:,i] = imported_data[:, 7*i]

#create central figure
fig = go.Figure()


xv, yv, zv = np.meshgrid(np.arange(5),np.arange(5),np.arange(5))
xv = np.reshape(xv,125)
yv = np.reshape(yv,125)
zv = np.reshape(zv,125)

#add traces for each timestep  
for step in np.arange(0, len(tbv)):
    fig.add_trace(
        go.Scatter3d(
            visible = False,
            x=xv,
            y=yv,
            z=zv,
            mode='markers',
            marker=dict(
            size=12,
            color=tbv[step,xv+5*yv+25*zv],
            cmin = -80,
            cmax = -20,
            colorbar = dict(title = 'Voltage', tickmode = "array"),
            colorscale='Turbo',   # choose a colorscale
            opacity=0.8
        )
     )
)

fig.data[0].visible = True

steps = []
for i in range(len(fig.data)):
    step = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)},
              {"title": "Slider switched to step: " + str(i)}],  # layout attribute
    )
    step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
    steps.append(step)

sliders = [dict(
    active=10,
    currentvalue={"prefix": "Frequency: "},
    pad={"t": 50},
    steps=steps
)]

fig.update_layout(
    sliders=sliders
)

fig.show()