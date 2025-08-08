import pyvista as pv
import numpy as np
from glob import glob
import re, os

# --- ordenar "naturalmente" density_00001, density_00010, ...
def _num_key(s):
    m = re.search(r"(\d+)(?=\.vtk$)", os.path.basename(s))
    return int(m.group(1)) if m else 0

files = sorted(glob("outputs/density_*.vtk"), key=_num_key)
assert files, "No se encontraron archivos density_*.vtk"

# --- rango global para una paleta estable
gmin, gmax = np.inf, -np.inf
for f in files:
    arr = pv.read(f)["density"]
    gmin = min(gmin, float(np.min(arr)))
    gmax = max(gmax, float(np.max(arr)))
clim = (gmin, gmax)

# --- primer frame y escena
mesh0 = pv.read(files[0])                      # vtkImageData (STRUCTURED_POINTS)
p = pv.Plotter(window_size=(1024, 768))
vol = p.add_volume(mesh0, scalars="density", cmap="inferno",
                   clim=clim, opacity="sigmoid", shade=True)
p.add_axes()
p.camera.zoom(1.2)

# --- elige salida
make_mp4 = True
if make_mp4:
    p.open_movie("kelvin_helmholtz.mp4", framerate=24)   # requiere imageio-ffmpeg
else:
    p.open_gif("kelvin_helmholtz.gif")                   # requiere imageio

p.show(auto_close=False, interactive=False)

# --- animar actualizando escalares IN-PLACE
arr0 = mesh0.point_data["density"]   # referencia al array interno
for f in files:
    arr_new = pv.read(f).point_data["density"]
    arr0[:] = arr_new                # actualización in-place
    p.render()
    p.write_frame()

p.close()
print("Listo: kelvin_helmholtz.mp4" if make_mp4 else "Listo: kelvin_helmholtz.gif")
