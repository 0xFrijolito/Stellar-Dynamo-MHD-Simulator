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


def visualize_voxels_animation(result, field, output_gif=None):
    cells_per_dim = result["cells_per_dim"]
    max_x, max_y, max_z = result["max_coords"]
    total_steps = result[field].shape[0]

    dx = max_x / cells_per_dim
    dy = max_y / cells_per_dim
    dz = max_z / cells_per_dim

    # Crear grid base
    grid = pv.ImageData()
    grid.dimensions = (cells_per_dim + 1, cells_per_dim + 1, cells_per_dim + 1)
    grid.origin = (0.0, 0.0, 0.0)
    grid.spacing = (dx, dy, dz)

    # Añadir IDs para rastrear celdas
    grid.cell_data["cell_ids"] = np.arange(grid.n_cells)

    # Definir región de recorte (esquina)
    corner_size = min(cells_per_dim, 64)
    clip_bounds = (
        max_x - dx * corner_size, max_x,
        max_y - dy * corner_size, max_y,
        max_z / 2, max_z
    )

    clipped = grid.clip_box(bounds=clip_bounds, invert=True)

    # Verificar que el recorte no esté vacío
    if clipped.n_cells == 0 or "cell_ids" not in clipped.cell_data:
        raise RuntimeError("Recorte vacío o sin IDs de celdas.")

    # Obtener IDs originales de las celdas recortadas
    clipped_ids = clipped.cell_data["cell_ids"]

    # Asignar datos del primer paso
    scalar_data = result[field][0].flatten(order="F")
    clipped.cell_data[field] = scalar_data[clipped_ids]

    # Crear plotter
    plotter = pv.Plotter()
    plotter.add_mesh(clipped, scalars=field, cmap="viridis", show_scalar_bar=True)

    if output_gif:
        plotter.open_gif(output_gif)

    plotter.show(interactive_update=True, auto_close=False)

    # Animación por pasos de tiempo
    for step in range(1, total_steps):
        scalar_data = result[field][step].flatten(order="F")
        clipped.cell_data[field] = scalar_data[clipped_ids]

        plotter.add_text(
            f"Step: {step}/{total_steps - 1}",
            position="lower_edge",
            name="step_text",
            font_size=10,
        )

        if output_gif:
            plotter.write_frame()
        else:
            plotter.update()
            time.sleep(0.05)

    plotter.close()


# Ejecutar simulación
result = read_simulation("output.bin")
visualize_voxels_animation(
    result, 
    field="density",
    output_gif="simulation_animation.gif"  # Opcional: guardar como GIF
)