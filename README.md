# 2D Potential Flow Solver: Vortex & Source Panel Method for NACA 4-digit series foils

**Authors:** Bora Ugurcan Ekizoglu, Kurt Spiteri  
**Language:** Julia (Pluto.jl)  

![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/BUEkizoglu/Panel-Methods)
![GitHub last commit](https://img.shields.io/img/shields/github/last-commit/BUEkizoglu/Panel-Methods)

## 📌 Overview
This repository contains an interactive, fully parameterized interactive 2D numerical solver built in Julia (using Pluto.jl) to estimate the lift and flow field around fully submerged NACA 4-digit series hydrofoils/airfoils. 

The tool utilizes inviscid potential flow theory, allowing users to seamlessly toggle between **First-Order Source Panel** and **First-Order Vortex Panel** methodologies to evaluate continuous pressure, velocity magnitudes, and integrated force coefficients as well as the Kutta-Joukowski lift at various geometric angles of attack.

<p align="center">
  <img src="img/NACA_velocity_magnitude.png" alt="Velocity Magnitude Contour over NACA Foil" width="800">
  <br>
  <i>Figure 1: Interactive flow field visualization (Velocity Magnitude) over a parametrised NACA foil.</i>
</p>

---

## ✨ Key Features
* **Interactive Geometry Generation:** Sliders to dynamically manipulate the number of panels ($N$), chord length, thickness, camber, and angle of attack for any NACA 4-digit foil.
* **Dual Singularity Solver:** Fully interchangeable between purely Source (thickness modeling) and purely Vortex (lifting modeling) mathematical Green's functions.
* **Kutta Condition Enforcement:** Explicit matrix augmentation to enforce the trailing edge stagnation point, allowing the mathematical generation of physical circulation and lift.
* **Aerodynamic Post-Processing:** Automated integration of the continuous pressure field (via Bernoulli's Equation) to calculate integrated lift and numerical drag, alongside theoretical Kutta-Joukowski checks.
* **Flow Visualization:** Superposition of the freestream and panel-induced velocity grids to map continuous Velocity Magnitude ($V$), directional vectors ($U_x$, $U_y$), and Pressure Coefficient ($C_p$) contours.

---

## 📐 Mathematical Methodology

The solver discretizes the continuous airfoil curve into $N$ flat line segments (panels). The potential flow velocity on any panel surface is governed by the boundary condition:

$$V_{i} = V_{\infty}\hat{p} + \sum^{n}_{j=1} \lambda_{j} \int_{s_{j}} \frac{\partial G(r)}{\partial s}ds_j$$

Where $V_{\infty}$ is the freestream velocity, $\hat{p}$ is the geometric normal ($\hat{n}$) or tangent ($\vec{t}$) vector, $\lambda_{j}$ is the unknown singularity strength, and $G$ is the Green's function.

### Green's Functions
The potential induced at a point $(x,y)$ by a singularity situated at $(x_s, y_s)$ is defined by:
* **Point Source:** $G = \frac{1}{2}\ln(r)$
* **Point Vortex:** $G = \arctan\left( \frac{y-y_{s}}{x-x_{s}}\right)$

### The Matrix System
By evaluating the influence of every panel on every other panel's collocation point (the midpoint), we assemble an aerodynamic influence matrix $[A]$ and an excitation vector $[B]$ (the freestream boundary condition). 

$$\left[\begin{array}{ccc}A_{1,1} & \cdots & A_{1, N} \\ \vdots & \ddots & \vdots \\ A_{N, 1} & \cdots & A_{N, N}\end{array}\right]\left[\begin{array}{c}\lambda_1 \\ \vdots \\ \lambda_N\end{array}\right]=\left[\begin{array}{c}B_1 \\ \vdots \\ B_N\end{array}\right]$$.

### The Kutta Condition
To prevent infinite velocities around the sharp trailing edge, the system is modified to physically lock the trailing-edge stagnation point. We sacrifice the final matrix equation and replace it with the Kutta Condition, enforcing zero net vorticity at the trailing edge:

$$1 \cdot \gamma_{top} + 1 \cdot \gamma_{bottom} = 0$$

<p align="center">
  <img src="img/NACA_pressure.png" alt="Pressure Coefficient Field" width="600">
  <br>
  <i>Figure 2: Pressure Coefficient ($C_p$) contour mapping the leading-edge suction peak.</i>
</p>

---

## 📈 Validation and Grid Convergence

A rigorous validation study (`NACA_XXXX_PM_validation_convergence.jl`) accompanies the core tool. The matrix assembly and mathematical integration are validated against closed-body sanity checks (circles and ellipses) before being applied to the lifting foil.

As the panel resolution ($N$) increases, the solver demonstrates textbook spatial convergence. The integrated Lift Coefficient ($C_l$) successfully flatlines to its true grid-independent value, perfectly matching the theoretical **Kutta-Joukowski** theorem. Simultaneously, as expected from D'Alembert's Paradox, the numerical spatial truncation error (Integrated $C_d$) decays exponentially toward zero.

<p align="center">
  <img src="img/NACA_convergence.png" alt="Vortex Strength and Force Convergence" width="800">
  <br>
  <i>Figure 3: Grid convergence study demonstrating stable singularity distribution and asymptotic lift convergence.</i>
</p>

---

## 🚀 How to Use

This project is built as a reactive **Pluto.jl** notebook. 

1. Install [Julia](https://julialang.org/downloads/).
2. Open the Julia REPL and install Pluto:
   ```julia
   julia> ]
   (@v1.x) pkg> add Pluto
   julia> import Pluto
   julia> Pluto.run()
