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

def visualize_sphere_animation(result, field="temperature", output_gif=None):
    # Configuración inicial
    radius = result["positions"][..., 0]
    theta  = result["positions"][..., 1]
    phi    = result["positions"][..., 2]

    # Conversión manual a cartesianas
    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)

    # Construimos la malla estructurada
    mesh = pv.StructuredGrid(x, y, z)
    mesh.plot()

    
# Ejecutar simulación
result = read_simulation("output.bin")
visualize_sphere_animation(
    result, 
    field="temperature",
    output_gif="simulation_animation.gif"  # Opcional: guardar como GIF
)