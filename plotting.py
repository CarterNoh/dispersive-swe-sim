import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import cv2
import PIL

# --- Configuration ---
data_dir = "build/Debug/data"
frames_dir = "frames"
output_video = os.path.join(frames_dir, "fluid_simulation.mp4")
save_images = True  # Set to True to save individual frames as PNGs
fps = 30

# Get all CSV files and sort them numerically by the tick number
files = glob.glob(os.path.join(data_dir, "frame_*.csv"))
files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))
if not files:
    print(f"No data files found in {data_dir}!")
    exit()
# Clear the frames directory
if not os.path.exists(frames_dir):
    os.makedirs(frames_dir)
else:
    for file in glob.glob(os.path.join(frames_dir, "*.png")):
        os.remove(file)

# Load the terrain
terrain_location = os.path.join(data_dir, "terrain.csv")
if not os.path.exists(terrain_location):
    print(f"Error: Terrain data file {terrain_location} not found!")
    exit()
terrain_data = np.genfromtxt(terrain_location, delimiter=',')

# Set up the figure
# plt.switch_backend('Agg')
width, height = 1000, 800
fig = plt.figure(figsize=(width/100, height/100))
ax = fig.add_subplot(111, projection='3d')
grid_size = terrain_data.shape[0]
X, Y = np.meshgrid(np.arange(grid_size), np.arange(grid_size))


def update(frame_idx):
    ax.clear()
    ax.set_zlim(-15, 25)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title(f"Tick: {frame_idx}")
    
    # # Plot terrain
    # ax.plot_surface(X, Y, terrain_data, cmap='copper', alpha=0.6, antialiased=True) #BrBG
    
    # Load current frame water data
    water_data = np.genfromtxt(files[frame_idx], delimiter=',')
    water_surface = water_data + terrain_data # add terrain height to water column height to get surface height
    # water_surface[water_data < 0.01] = np.nan

    # Update surface
    ax.plot_surface(X, Y, water_surface, cmap='winter', alpha=0.8, antialiased=True)
    ax.set_title(f"Tick: {frame_idx}")
    if save_images:
        plt.savefig(f"{frames_dir}/frame_{frame_idx:04d}.png")
    fig.canvas.draw()
    return cv2.cvtColor(np.asarray(fig.canvas.buffer_rgba()), cv2.COLOR_RGBA2BGR)
        
# Create the video
print(f"Generating video...")
out = cv2.VideoWriter(output_video, cv2.VideoWriter_fourcc(*'mp4v'), fps, (width, height))
if not out.isOpened():
    raise RuntimeError("VideoWriter failed to open")
for i in range(len(files)):
    frame = update(i)      
    out.write(frame)
    if i % 10 == 0:
        print(f"Processed frame {i}/{len(files)}")
out.release()
plt.close(fig)
print(f"Done! Video saved to {output_video}")