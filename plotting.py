import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
import glob
import ffmpeg

# --- Configuration ---
data_dir = "build/Debug/data"
output_video = "fluid_simulation.mp4"
save_images = True  # Set to True to save individual frames as PNGs
fps = 30

# 1. Get all CSV files and sort them numerically by the tick number
files = glob.glob(os.path.join(data_dir, "frame_*.csv"))
files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))

if not files:
    print(f"No data files found in {data_dir}!")
    exit()

# Load the terain and set up the grid for plotting
terrain_location = os.path.join(data_dir, "frame_-1.csv")
if not os.path.exists(terrain_location):
    print(f"Error: Terrain data file {terrain_location} not found!")
    exit()
terrain_data = np.genfromtxt(terrain_location, delimiter=',')
grid_size = terrain_data.shape[0]
x = np.arange(grid_size)
y = np.arange(grid_size)
X, Y = np.meshgrid(x, y)

# 2. Setup the figure for 3D plotting
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Set constant axis limits based on your C++ Sim2D.h parameters
ax.set_zlim(-15, 25) 
ax.set_title("Fluid Simulation")

# Initialize the surface plot object
plot = [ax.plot_surface(X, Y, terrain_data)] #, cmap='viridis', edgecolor='none'

def update(frame_idx):
    ax.clear()
    ax.set_zlim(-15, 25)
    ax.set_title(f"Tick: {frame_idx}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    # # plot terrain
    # ax.plot_surface(X, Y, terrain_data, cmap='copper', alpha=0.6, antialiased=True) #BrBG
    
    # Load current frame water data
    water_data = np.genfromtxt(files[frame_idx], delimiter=',')
    water_surface = water_data #+ terrain_data # add terrain height to water column height to get surface height
    # water_surface[water_data < 0.01] = np.nan

    # Update surface
    surf = ax.plot_surface(X, Y, water_surface, cmap='winter', alpha=0.8, antialiased=True)
    ax.set_title(f"Tick: {frame_idx}")
    
    # Optional: Save individual frame as PNG
    if save_images:
        if not os.path.exists("frames"):
            os.makedirs("frames")
        plt.savefig(f"frames/frame_{frame_idx:04d}.png")
        
    return surf

# Create the animation
print(f"Generating animation from {len(files)} frames...")
ani = FuncAnimation(fig, update, frames=len(files[1:]), interval=1000/fps, blit=False)

# 4. Save the video
# Note: Requires 'ffmpeg' installed on your system for .mp4
try:
    ani.save(output_video, writer='ffmpeg', fps=fps)
    print(f"Success! Video saved as {output_video}")
except Exception as e:
    print("FFmpeg not found. Saving as GIF instead...")
    ani.save("fluid_simulation.gif", writer='pillow', fps=fps)
    print("Video saved as fluid_simulation.gif")
    exit()

plt.show()