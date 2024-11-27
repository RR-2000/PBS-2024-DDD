import trimesh
import numpy as np
import pyrender

# Load the mesh from an OBJ file
file_name = 'tub'
mesh = trimesh.load(f'{file_name}.obj')

# Sample points on the surface of the mesh
num_points = 20000  # Adjust this to sample more or fewer points
sampled_points, face_indices = mesh.sample(num_points, return_index=True)
with open(f'{file_name}.xyz', 'w') as f:
    f.write('x y z\n')  # Optional header
    for point in sampled_points:
        f.write(f'{point[0]:.6f},{point[1]:.6f},{point[2]:.6f}\n')

scale = 0.5
# Display or save the point cloud (optional)
point_cloud = trimesh.points.PointCloud(sampled_points*scale)
point_cloud.show()  # Opens an interactive viewer if supported