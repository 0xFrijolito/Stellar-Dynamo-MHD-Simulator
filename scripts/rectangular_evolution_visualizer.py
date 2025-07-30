import struct
import numpy as np
import pyvista as pv
import time

def read_simulation(filename):
    with open(filename, "rb") as f:
        # Leer metadata
        cells_per_dim = struct.unpack("I", f.read(4))[0]
        max_x = struct.unpack("d", f.read(8))[0]
        max_y = struct.unpack("d", f.read(8))[0]
        max_z = struct.unpack("d", f.read(8))[0]
        
        total_cells = cells_per_dim ** 3
        
        # Leer posiciones (3 f64 por celda)
        positions_size = total_cells * 3 * 8
        positions_data = f.read(positions_size)
        
        # Calcular tamaño de cada paso de tiempo
        step_size = total_cells * 3 * 8  # 3 campos escalares, f64 = 8 bytes
        
        # Leer el resto del archivo
        data = f.read()
        total_steps = len(data) // step_size

        print(f"Grid: {cells_per_dim}³ | Steps: {total_steps}")

        # Convertir datos de simulación a array numpy
        full_array = np.frombuffer(data, dtype=np.float64)
        full_array = full_array.reshape((total_steps, 3, cells_per_dim, cells_per_dim, cells_per_dim))

        return {
            "temperature": full_array[:, 0],
            "pressure":    full_array[:, 1],
            "density":     full_array[:, 2],
            "cells_per_dim": cells_per_dim,
            "max_coords": (max_x, max_y, max_z),
            "positions": np.frombuffer(positions_data, dtype=np.float64).reshape((cells_per_dim, cells_per_dim, cells_per_dim, 3))
        }
    
def visualize_voxels_animation(result, field="temperature", output_gif=None):
    # Configuración inicial
    x = result["positions"][..., 0]
    y = result["positions"][..., 1]
    z = result["positions"][..., 2]

    cells_per_dim = result["cells_per_dim"]

    max_x = result["max_coords"][0]
    max_y = result["max_coords"][1]
    max_z = result["max_coords"][2]

    dx = max_x / cells_per_dim
    dy = max_x / cells_per_dim
    dz = max_x / cells_per_dim

    # Construimos la malla estructurada
    mesh = pv.StructuredGrid(x, y, z)
    mesh[field] = result[field][0].flatten(order='F')

    # Definir región de recorte (esquina)
    corner_size = cells_per_dim / 2
    clip_bounds = (
        max_x - dx * corner_size, max_x,
        max_y - dy * corner_size, max_y,
        max_z - dz * corner_size, max_z
    )

    # Aplicar recorte conservando los IDs originales
    clipped = mesh.clip_box(bounds=clip_bounds, invert=True)

    pltr = pv.Plotter(window_size=[512, 512])
    pltr.add_mesh(
        clipped,
        scalars=field,
        cmap='inferno',
        show_scalar_bar=True,
    )

    pltr.open_gif(output_gif)
    for step in range(500):
        mesh[field] = result[field][step].flatten(order='F')
        clipped = mesh.clip_box(bounds=clip_bounds, invert=True)
        pltr.update_scalars(clipped[field], render=False)
        pltr.write_frame()

    pltr.show()

# Ejecutar simulación
result = read_simulation("output.bin")
visualize_voxels_animation(
    result, 
    field="temperature",
    output_gif="ctm.gif"  # Opcional: guardar como GIF
)