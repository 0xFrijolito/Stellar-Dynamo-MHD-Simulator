# Stellar Dynamo MHD Simulator 🔭

Un simulador simplificado de **magnetohidrodinámica (MHD)** orientado al estudio del **dínamo estelar**, basado en el método de volúmenes finitos en 3D.

---

## 🧱 Descripción del proyecto

Este proyecto busca implementar un solver MHD ideal en tres coordenadas (rectangular, esférica o cilíndrica) con enfoque en la simulación de fenómenos de dínamo en estrellas. Aunque aún no está completo, la arquitectura permite:

- Inicializar grillas estructuradas (Rectangular / Esférica / Cilíndrica)
- Almacenar métricas locales de celda para discretización en coordenadas curvilíneas
- Kernel genérico para gradiente y laplaciano escalar
- Paso numérico básico para ecuación de calor (test sencillo de difusión)

El diseño está preparado para integrar modelos más avanzados.

---

## 🚀 Características actuales

- Inicialización de grilla estructurada (grid genérico)
- Almacenamiento de posiciones, métricas, tamaños, campos escalares y vectoriales
- Funciones genéricas para gradiente y laplaciano escalar (`scalar_gradient`, `scalar_laplacian`)
- Solver de difusión térmica (`heat_equation_solver_step`) como ejemplo de evolución temporal
- Formato binario simple de exportación de datos (metadata + posiciones + temperatura/presión/densidad)
- Visualización en Python (via PyVista) de los resultados guardados

---

## 🧪 Objetivo final

Crear un simulador MHD capaz de resolver el sistema completo de ecuaciones ideales (masa, momentum, energía y campo magnético) usando el método de **volúmenes finitos** con:

- Reconstrucción de alta resolución (e.g. WENO, datos de flujo HLLD)
- Transporte magnético **sin divergencia** (`∇·B = 0`) via *Constrained Transport*
- Manejo de múltiples geometrías (especially esférica para envolturas estelares)
- Condiciones de frontera físicas (rotación estelar, convección, energía radiativa)
- Posibilidad de resolver dínamos a pequeña y gran escala

---

## 📝 Estado actual (To‑Do)

### ✅ Implementado (core básico)
- Estructura de grilla genérica con `positions`, `metrics`, `sizes`, campos
- Inicialización vía funciones callback de perfiles escalar/vectoriales
- Gradient y Laplaciano genéricos usando métricas
- Ejemplo de solver térmico (difusión) en geometría rectangular
- Output binario y script Python de lectura y visualización

### 🔧 Por implementar
- [ ] Ecuaciones completas ideal-MHD (masa, momentum, energía, inducción)
- [ ] Reconstrucción y cálculo de fluxes usando esquemas conservativos (WENO, MUSCL, HLLD/HLLC)
- [ ] Transporte magnético sin divergencia (Constrained Transport)
- [ ] Esquema temporal avanzado (Runge–Kutta o Strang splitting)
- [ ] Condiciones de frontera adecuadas (por ejemplo superficie de estrella, límites periódicos)
- [ ] Inicialización de perfiles físicos reales (estratificación, rotación, magnetización)
- [ ] Implementación de pruebas de validación (Orszag–Tang, Kelvin–Helmholtz, etc.)
- [ ] Parallelización (OpenMP o MPI)
- [ ] Soporte para geometría esférica completa (o cónica)

---

## 📦 Uso

```bash
# Clonar repositorio
git clone https://github.com/tuusuario/stellar-dynamo-mhd.git
cd stellar-dynamo-mhd

# Compilar (ejemplo con Makefile)
./build.bat

# Ejecutar ejemplo de difusión térmica
./bin/mhd_solver

# Reproducir resultados en Python
python3 viz_simulation.py output.bin
# Stellar-Dynamo-MHD-Simulator
