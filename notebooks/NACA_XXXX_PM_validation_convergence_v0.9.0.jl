### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ c14bcbc0-eb5d-11ee-1bf8-8fcaaf8e34c2
begin
	using PlutoUI
	using Plots
	using BenchmarkTools
	using LinearAlgebra
	using Interpolations
	using OutMacro
	using GeometryBasics
end

# ╔═╡ 3af5e095-c9c0-4ec3-ae02-c2e060b6c153
md"""

# Implmentation of a Vortex/Source-Panel Method for a Submerged Hydrofoil 

Bora Ugurcan Ekizoglu, Kurt Spiteri

"""

# ╔═╡ cc29fa3e-cfba-4de0-ba7e-451a33d8714d
md"""
# 1. Introduction

This Pluto notebook is an extention of the numerial tool `NACA_XXXX_PM.jl` which is created to calculate the lift around a NACA $4$-digit series foil using potential theory and vortex panel method. This notebook explains in great detail the methodology used for the numerical tool and includes extensive convergence studies. This notebook will be extended in the future to include an analytical or experimental validation case as well. 

"""

# ╔═╡ 1c44dba8-ed15-4908-a3ba-968a7ca88d8b
md"""
## 1.1. Important Syntax
To make this notebook fully interactive and adaptable, the core functions rely on several key variables, keyword arguments, and operational flags. Below is a breakdown of their physical and mathematical significance.

### 1.1.1. Geometric & Domain Variables
* **`N` (Integer):** The total number of panels used to discretise the airfoil. Higher values reduce spatial truncation error and improve boundary resolution, but increase computational cost.
* **`chord`, `thickness`, `camber` (Floats):** The defining geometric parameters of the NACA 4-digit airfoil.
* **`angle_deg` (Float):** The geometric angle of attack. This physically rotates the airfoil coordinates relative to the fluid domain.
* **`XY`, `x_ran`, `y_ran`:** The spatial grid matrices and range vectors that define the fluid domain where velocities and pressures are calculated for visualization.

### 1.1.2. Solver Flags & Keyword Arguments
* **`U_free` (Tuple):** The freestream velocity vector defined as $(u_{x}, u_{y})$. It dictates the magnitude and direction of the incoming flow.
* **`G` (Symbol):** The Green's function flag (`G=source` or `G=vortex`). This dictates the fundamental mathematical solution applied to the panels:
  * `:source` creates displacement (pushing fluid away), used for capturing thickness effects.
  * `:vortex` creates circulation (spinning fluid), used for capturing lifting effects.
* **`kutta` (Boolean):** The Kutta Condition flag (`true` or `false`). When set to `true`, it overwrites the final row of the influence matrix to explicitly enforce that the circulation at the trailing edge is zero. This prevents infinite velocities at the sharp trailing edge and allows the airfoil to generate physical lift.
* **`loopdir` (Symbol):** The panel winding direction (`:CW` for Clockwise, `:CCW` for Counter-Clockwise). This determines the orientation of the geometric normal vectors ($\hat{n}$). It is critical that normals point *outward* into the fluid domain for the boundary conditions to hold true.

### 1.1.3. Post-Processing Variables
* **`γ` or `q` (Vector):** The array of solved singularity strengths for every panel. These are the primary mathematical outputs of the `solve` function and represent the local surface velocity (for vortices) or local outflow (for sources).
* **`plot_var` (Symbol):** The visualization toggle used in the flow plotting function.
  * `:magnitude`: Displays the absolute scalar velocity field ($U$).
  * `:u`: Isolates the horizontal velocity field $(u_{x})$, useful for visualizing flow acceleration over the foil.
  * `:pressure`: Applies Bernoulli's principle to convert the velocity field into a continuous Pressure Coefficient $(C_{p})$ contour, visualizing the physical forces acting on the fluid.
"""

# ╔═╡ 6de9ba9e-776c-415b-9d9b-1083e6e1540a
md"""

## 1.2. Objectives

The main objective of this project is to implement a vortex panel methodology in a Pluto notebook to develop a fast and user friendly tool to estimate the lift of a fully submerged NACA $4$-digit series hydrofoil that is not influenced by free surface effects. The methodology is to be carried out in $2D$. 

The vortex panel method is chosen as the latter method would result in zero net force for a foil at an angle, however the code is constructed to be fully interchangeable between vortex and source panel methods by simply changing the Green's function ($G$) argument of `solve( )` function.

"""

# ╔═╡ eed1913e-57fd-42c0-9430-2f676d059efc
md"""

## 1.3. Novelty

The novel feature of this methodology is that any NACA $4$-digit series foil profile under any incident velocity condition can be modelled and evaluated interactively by the user. The methodology can also be extended to other arbitrary forms. This project was undertaken as a challenge to establish the vortex /source panel method from scratch in Pluto for a fully parametrised hydrofoil since such a methodology that is thoroughly verified and validated was not found.

"""

# ╔═╡ 37cc0cb5-ee73-4791-b787-c0bc145d225c
md"""

# 2. Theoretical Overview

"""

# ╔═╡ e7624d6e-17aa-4dfd-b677-bed2116ae0e3
md"""

## 2.1. Introduction to the Vortex Panel Method

The vortex panel method is based on the philosophy of applying a vortex sheet to a hydrofoil surface. The vortex sheet must be of such strength that the foil surface becomes a streamline of the flow (Anderson, 2017). For establishing such a vortex sheet, the surface must be discretised and boundary conditions relevant to the incident velocity on the hydrofoil must be applied. The vortex strength on each panel is then computed.

The vortex panel method is governed by the following equation which is the velocity boundary condition on the panel surface:

$V_{i} = V_{\infty}\hat{p} + \sum^{n}_{j=1} \lambda_{j} \int_{s_{j}} \frac{\partial G(r)}{\partial s}ds_j$

where $V_{i}$ is the surface velocity on panel $i$, $V_{\infty}$ is the free stream velocity, $\hat{p}$ is the panel normal ($\hat{n}$) or ($\vec{t}$) tangent vector for source and vortex panels respectively, $\lambda_{j}$ is the source ($q_{j}$) or vortex ($\gamma_{j}$) strength for panel $j$, $G$ is the Green's function and $r=\sqrt{(x_{i}-x{j})^2+(y_{i}-y_{j})^2}$ is the radial distance between the collocation points of panel $i$ and panel $j$.

"""

# ╔═╡ af51eea8-1423-4720-96c3-837d12406fea
md"""

## 2.2. Implementation

To apply the Vortex/Source Panel Method, the surfaces of the foil must be discretised into a numebr of panels. Once the foil has been split up into $N$ panels constructed by $N+1$ points. Then:

- A vortex ($\gamma_{i}$) or source ($q_{i}$) distribution of constant strength along the panel length ($ds$) is applied on every panel (1 ≤ i ≤ N). This assumption limits the Vortex/Source Panel Method to the First Order. 

- At every panel, the velocity potential influence ($\phi_{ij}$) of a panel $i$ panel on itself and all the other panels $\phi$ (as a function of the vortex strength $\gamma_{i}$ ) is computed. These are assembled into the influence matrix $[A]$. The diagonal of $[A]$ holds the information for the potential influence of a vortex/source sheet on its own collocation point which o selected to be the midpoint of the panel. The diagonal is strictly set to ($\pi$) because of the non-dimensional convention selected with the Green's ($G$) function.  

- The flow velocity normal ($\partial \phi / \partial \vec{n}$) or tangential ($\partial \phi / \partial \vec{t}$) to each panel on the hydrofoil may also be computed which equals to the negative of the freestream velocity component in panel direction for potential theory. These are assembled into the Excitation Matrix B.

- Using these conditions, a matrix network may be set up to solve for the source $q_{i}$ strengths or vortex strengths $\gamma_i$ which sets up $[\Lambda]$.

$\left[\begin{array}{ccc}A_{1,1} & \cdots & A_{1, N} \\ \vdots & \ddots & \vdots \\ A_{N, 1} & \cdots & A_{N, N}\end{array}\right]\left[\begin{array}{c}\lambda_1 \\ \vdots \\ \lambda_N\end{array}\right]=\left[\begin{array}{c}B_1 \\ \vdots \\ B_N\end{array}\right]$.

"""

# ╔═╡ b92edf15-fa1f-4930-ba73-7c73e7363055
md"""

## 2.3. The Kutta Condition

For a stready state flow over a lifting surface for the flow leaves the top and bottom surfaces smoothly at the trailing edge according to potential theory. However, this cannot be achieved by simply solving the system $[A][Λ]=[B]$. This is because of the fact that vortex sheets, connecting at sharp edges like the trailing edge of a foil will induce high magnitudes of opposite tangential velocity on each other. Then, the smooth flow around the trailing edge is achieved by overriding the boundary condition for top and bottom side trailing edge panels by making sure that the stagnation point at the trailing edge physically coinsides with the trailing edge panel endpoint. This special boundary condition is called the Kutta condition (Anderson, 2017). 

The equation correspondong to zero vorticity at the trailing edge is given as follows:

$1 \cdot \gamma_{top} + 1 \cdot \gamma_{bottom} =0$

Where $\gamma_{top}$ and $\gamma_{bottom}$ are the vortex strengths for top and bottom trailing edge panels respectively. Note that the system is now overdetermined given that there N unknown vortex strengths and N+1 available equations. So, the Kutta condition has to be enforced by sacrificing one equation from the system and replacing it with the Kutta condition. This is ypically done by sacrificing the last equation. Thus, the new matrix system with the Kutta condition can be given as follows:

$\left[\begin{array}{ccc}A_{1,1} & \cdots & A_{1, N} \\ \vdots & \ddots & \vdots \\ A_{N-1, 1} & \cdots & A_{N-1, N-1} \\ 1 & 0 & 1 \end{array}\right]\left[\begin{array}{c}\gamma_1 \\ \vdots \\ \gamma_{N-1} \\ \gamma_{k} \end{array}\right]=\left[\begin{array}{c}B_1 \\ \vdots \\ B_N \\  0 \end{array}\right]$

Where $\gamma_{k}$ is the vortex strenght for top and bottom trailing edges. 

"""

# ╔═╡ a90a1cf8-2fe5-420d-9699-96bda5b00988
md"""

# 3. Hydrofoil Geometry

NACA $4$-digit foil series was considered. The hydrofoil was fully parametrised by constructing the upper surface of the foil dependin on the $x$-coordinate with the following equation, then mirroring it to close the loop:

$5 t (a_{0} \sqrt{x} + a_{1} x + a_{2} x^2 + a_{3} x^3 + a_{4} x^4)$

Where $a_{0} = 0.2969$, $a_{1}= -0.1260$, $a_{2} = -0.3516$, $a_3 = 0.2843$, $a_{4}= -0.1036$ and $t$ is the thickness of the foil.

Symmetry was used for the lower surface, and a correction was implemented to account for camber and to completely close off the trailing edge panels.

The formulations were retrieved from (Ladson & Brooks, 1975). 

"""

# ╔═╡ 0445668c-3559-4159-957c-df9e5aa40c76
md"""

## 3.1. Hydrofoil Parameters

"""

# ╔═╡ 157f9c00-4ff3-4bc6-ade1-c78886042128
begin
	Num_pointsS = @bind N Slider(5:1:256, default=128, show_value=true)
	md"""Points (N) $Num_pointsS"""
end

# ╔═╡ 740e3ad2-4397-4afa-8db1-cb2718baa1e3
num_points = Int64(round((N+2)/2))

# ╔═╡ 18dbf7f0-e90b-45e5-8de8-82db0abb6b9e
md"""

`N` is defined as a counter of the total number of panels for the whole foil loop. It takes into account both the upper and lower surfaces. The `num_points` variable is the number of points defining the upper "arch" of the foil. 

"""

# ╔═╡ c7addb5a-0a05-423b-9f73-0e237d189d11
begin
	ChordS = @bind chord Slider(0.0:0.1:5, default =2, show_value=true)	
	md"""Chord $ChordS"""
end

# ╔═╡ d155e9a4-2ae6-4229-a810-633d49df1354
begin
ThicknessS = @bind thickness Slider(0.0:1:50, default=20, show_value=true)
	md"""Thickness [% of chord] $ThicknessS"""
end

# ╔═╡ 797534ac-b52a-4e58-812d-0a7962615590
begin
	CamberS = @bind camber Slider(0.0:1:20, default=0, show_value=true)
	
	md"""Camber [% of chord] $CamberS""" 
end

# ╔═╡ 9eea78c8-3233-44d2-b375-225bce7ebc8f
begin
	Rotation_AngleS = @bind angle_deg Slider(-30:1:30, default=-10, show_value=true)
	md"""Rotation Angle [degrees] $Rotation_AngleS"""
end

# ╔═╡ f1ad33b3-6f11-4a01-a83a-9b2ebe3ca382
md"""
Freestream velociy vector is non dimensional and set to $U=(u_{x}=1,u_{y}=0)$ for all calculations. Angle of attack is determined by the rotation angle.
"""

# ╔═╡ 44686d98-0659-4157-bd34-c932b7e78d61
md"""

# 4. Results: Flow Field and Lift

"""

# ╔═╡ 6b53a7e2-692b-4093-8afb-2c07f4325e5c
md"""

# 5. Core Functions

"""

# ╔═╡ 07da2b45-6d9b-4883-8b8f-4625405018ac
md"""

## 5.1. Defining foil panels

"""

# ╔═╡ 4753ed02-0ad1-4971-8fef-219d2862b85f
begin
	function naca_airfoil(chord, thickness, camber, num_points) # Function to generate NACA airfoil coordinates
	    # Define NACA parameters
	    m = camber / 100  # Maximum camber as a fraction of chord
	    p = 0.4           # Location of maximum camber as a fraction of chord
	    t = thickness / 100  # Maximum thickness as a fraction of chord
		
	    a0 = 0.2969
	    a1 = -0.1260
	    a2 = -0.3516
	    a3 = 0.2843
		a4 = -0.1036
	
	    # Calculate coordinates
	    x = collect(range(0, stop=1, length=num_points))
	    yt = 5 * t * (a0 * sqrt.(x) .+ a1 * x .+ a2 * x.^2 .+ a3 * x.^3 .+ a4 * x.^4)
	    yc = zeros(length(x))
	    dyc_dx = zeros(length(x))
	
	    # Calculate camber line and derivative
	    if camber > 0
	        yc .= m / p^2 * (2 * p * x .- x.^2)
	        dyc_dx .= 2 * m / p^2 * (p .- x)
	    end
	
		#X coordinates of the upper surfaces
	    xu = x .- yt .* sin.(atan.(dyc_dx))
		
		#X coordinates of the lower surfaces
	    xl = x .+ yt .* sin.(atan.(dyc_dx))
		
		#Y coordinates of the upper surfaces
	    yu = yc .+ yt .* cos.(atan.(dyc_dx))
		
		#Y coordinates of the lower surfaces
	    yl = yc .- yt .* cos.(atan.(dyc_dx))
	
	    # Scale airfoil coordinates by chord
	    xu *= chord
	    xl *= chord
	    yu *= chord
	    yl *= chord
	
	    return xu, yu, xl, yl
	end
end

# ╔═╡ 18416a19-0637-4de4-9c86-daae967f443b
begin
	function rotate_and_translate_airfoil(xu, yu, xl, yl, angle_deg, rotation_point, translation) # Function to rotate and translate airfoil
	    angle_rad = deg2rad(angle_deg)
	    
	    # Translate coordinates
	    xu_trans = xu .+ translation[1]
	    yu_trans = yu .+ translation[2]
	    xl_trans = xl .+ translation[1]
	    yl_trans = yl .+ translation[2]
	    
	    # Translate coordinates to rotate about the specified point
	    xu_trans .= xu_trans .- rotation_point[1]
	    yu_trans .= yu_trans .- rotation_point[2]
	    xl_trans .= xl_trans .- rotation_point[1]
	    yl_trans .= yl_trans .- rotation_point[2]
	    
	    # Rotate upper surface
	    xu_rot = xu_trans .* cos(angle_rad) .- yu_trans .* sin(angle_rad)
	    yu_rot = xu_trans .* sin(angle_rad) .+ yu_trans .* cos(angle_rad)
	    
	    # Rotate lower surface
	    xl_rot = xl_trans .* cos(angle_rad) .- yl_trans .* sin(angle_rad)
	    yl_rot = xl_trans .* sin(angle_rad) .+ yl_trans .* cos(angle_rad)
	    
	    # Translate back to original position
	    xu_rot .+= rotation_point[1] + translation[1]
	    yu_rot .+= rotation_point[2] + translation[2]
	    xl_rot .+= rotation_point[1] + translation[1]
	    yl_rot .+= rotation_point[2] + translation[2]
	
		#The final points on the upper and lower hydrofoil surfaces.
	    return xu_rot, yu_rot, xl_rot, yl_rot
	end
end

# ╔═╡ 37a088a7-5552-460b-b01a-ad0e521d7919
begin
function rotated_and_translated_airfoil(xu_rot, yu_rot, xl_rot, yl_rot) # Function to plot rotated and translated airfoil
    # Calculate midpoints
    xu_mid = (xu_rot[1:end-1] + xu_rot[2:end]) / 2
 	xl_mid = (xl_rot[1:end-1] + xl_rot[2:end]) / 2
    yu_mid = (yu_rot[1:end-1] + yu_rot[2:end]) / 2
    yl_mid = (yl_rot[1:end-1] + yl_rot[2:end]) / 2

    # Calculate midpoint between upper and lower curves at the trailing edge
    x_end_mid = (xu_rot[end] + xl_rot[end]) / 2
    y_end_mid = (yu_rot[end] + yl_rot[end]) / 2
	return xu_mid, xl_mid, yu_mid, yl_mid, x_end_mid, y_end_mid  
end
end

# ╔═╡ 4b780945-cce8-46ea-a28f-b4c2bfd603c6
begin
	function foil_panel(xu_r,xl_r,yu_r,yl_r) # Function to arrange panel endpoints into a CW loop			
		points_x = vcat(reverse(xl_r), xu_r[2:end])
		points_y = vcat(reverse(yl_r), yu_r[2:end])
		num_points_n = length(points_x)
		return points_x, points_y, num_points_n
	end
end

# ╔═╡ 6b332adc-359f-40d7-8396-273b35d23550
begin
	rotation_point = [0, 0] # Point of rotation 
	translation_y = 0.0 # Point of translation
	translation = [0, -(translation_y/2)] # Translation vector
	function NACA(N; chord=chord, thickness=thickness, camber=camber, angle_deg=angle_deg, rotation_point=rotation_point, translation=translation) # Function to generate all panel coordinates from panel number.
		n = Int64(round((N+2)/2))
		rotation_point = [0, 0]  
		translation = [0, -(translation_y/2)] # Translation vector
		xu, yu, xl, yl = naca_airfoil(chord, thickness, camber, n)
		xu_rot, yu_rot, xl_rot, yl_rot = rotate_and_translate_airfoil(xu, yu, xl, yl, angle_deg, rotation_point, translation)
		x_foil, y_foil, N_foil = foil_panel(xu_rot,xl_rot,yu_rot,yl_rot)
		return x_foil, y_foil
	end
end

# ╔═╡ 87345c05-f3c2-4225-b44d-b68d6e34509b
md"""

## 5.2. Implementation of Source/Vortex Panel Method

The $2D$ vesion of a panel is a line segment. In general this segment will run from some point $(x_{i},y_{i})$ to $(x_{i+1},y_{i+1})$ and will have a potential flow singularity distributed over it. The potential at point $(x,y)$ for a point source situated at $(x_{s}, y_{s})$ is proportional to the following expression:

$G(x,y,x_{s},y_{s})=\frac{1}{2}ln(r)=\frac{1}{4}ln(r^2)$

Where $r^2=(x-x_{s})^2+(y-y_{s})^2$ is the square distance from a point $(x_{s},y_{s})$ on a line segment to a point $(x,y)$ in space and $G$ is the Green's function. Similarly, the potential at point $(x,y)$ for a point vortex situated at $(x_{s}, y_{s})$ is proportional to the radial distance:

$G(x,y,x_{s},y_{s})= θ = arctan\left( \frac{y-y_{s}}{x-x_{s}}\right)$

The potential of the panel is the integrated superposition of the Green's function along the length of the panel segment for all panels as in the following expression:

$\phi_{i}(x,y)=λ_{i}\int^{s_{i+1}}_{i}G(x,y,x_{s},y_{s})ds=λ_{i}F^{\phi}_{i}(x,y)$

Where $\lambda_{i}$ is the source/vortex panel strength. The function $F^{\phi}_{i}$ is the source/vortex's influence function defining the source/vortex potential induced by segment $i$ per unit strength. The potential flow velocity is the gradient of the of the potential ($\phi_{i}$) as in the following equation:

$\vec{u_{i}}=\vec{\nabla}\phi_{i}=\lambda_{i}\vec{\nabla}F^{\phi}_{i}=\lambda_{i}F^{\vec{u}}_{i}, \quad F^{\vec{u}}_{i} = \left[\frac{\partial F^{\phi}_{i}}{\partial x},\frac{\partial F^{\phi}_{i}}{\partial y}\right]$

Where $F^{\vec{u}}_{i}$ is the influence function for the velocity ($\vec{u}$).

"""

# ╔═╡ 909de60b-8906-47c6-b1fb-5f0d3c104cbc
md"""
### 5.2.1. Numerical Implementation for Single Panel

There are four general functions used thtoughout the code and the below cell defines these general functions:
- `source( )`: Defines the source's Green's function
- `vortex( )`: Defines the votex's Green's function
- `potential( )`: Defines the influence function $F^{\phi}_{i}$ by estimating the integral using a Gaussian-quadratue.
- `velocity( )`: Defines the influence function $F^{\vec{u}}_{i}$ by estimating the derivative using finite differences.

An approximate method for calculating an integral is called a quadrature. A Gaussian-quadrature the accuracy of the trapezoid rule by sampling the function at special points instead of the interval boundaries. The folloving code calculates the potential and flow around a single panel assuming unit vortex strength then plots the resulting potential contour lines and velocity vectors over a $2D$ cartesian mesh grid. The code keeps the vector components packed together. 
"""

# ╔═╡ a11a44b8-20cb-4e6b-a22a-0408b47146f6
begin
	function source(x, y, xs, ys) # Green's function for 2D point source.
		return (1/2) * log((x - xs)^2 + (y - ys)^2) 
	end
	
	function vortex(x, y, xs, ys)  # Green's function for 2D point vortex.
	    return atan(y - ys, x - xs)
	end
	
	Gauss = 0.5*(1+sqrt(1/3)) # Gaussian-quadrature sample point  
	
	function potential(x,y,x0,y0,x1,y1, G=source) # Gaussian quadrature estimate of the potential influence function.
		dG(s) = G(x,y,x0*(1-s)+x1*s,y0*(1-s)+y1*s)
		h = sqrt((x1-x0)^2+(y1-y0)^2)
		return 0.5*h*(dG(Gauss)+dG(1-Gauss)) # int_0^h G ds
	end
	
	function velocity(x,y,x0,y0,x1,y1,G=source,ϵ=1e-12) # Finite difference estimate of the velocity influence function
		phi(X,Y) = potential(X,Y,x0,y0,x1,y1,G)
		return [(phi(x+ϵ,y)-phi(x-ϵ,y))/(2*ϵ),  # dphi/dx
	            (phi(x,y+ϵ)-phi(x,y-ϵ))/(2*ϵ)]  # dphi/dy
	end
end

# ╔═╡ 00565f60-a92d-4d32-85b0-a4d5c63b5f1e
begin
	function is_inside_or_on_boundary(x, y, xb, yb; tol=1e-8) # Function to filter grid points inside or on the body.
        inside = false
        j = length(xb)
        
        for i in 1:length(xb)
            xi, yi = xb[i], yb[i]
            xj, yj = xb[j], yb[j]
            
            cross_product = (y - yi) * (xj - xi) - (x - xi) * (yj - yi)
            if abs(cross_product) < tol
                if (min(xi, xj) - tol < x < max(xi, xj) + tol) && 
                   (min(yi, yj) - tol < y < max(yi, yj) + tol)
                    return true 
                end
            end

            if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
                inside = !inside
            end
            j = i
        end
        return inside
    end

    function mask_grid(XY, xb, yb; tol=1e-8) # Function to replace grid point values 	inside or on the body with NaN.
        X_mat, Y_mat = XY
        EX = copy(X_mat)
        EY = copy(Y_mat)
        # Loop through every point in the matrices
        for i in eachindex(EX)
            if is_inside_or_on_boundary(EX[i], EY[i], xb, yb, tol=tol)
                EX[i] = NaN # Silence the X coordinate
                EY[i] = NaN # Silence the Y coordinate
            end
        end
        return (EX, EY)
    end
end

# ╔═╡ 6f76e9d9-2dfa-4b94-9d49-d5acdcc70176
begin
	function meshgrid(x::StepRangeLen, y::StepRangeLen) # Function to create a grid from a range of variables.
	    X = [i for j in y, i in x]
	    Y = [j for j in y, i in x]
	    return X, Y
	end
end

# ╔═╡ 73a9aae0-2b30-4baa-80b1-028d6bc1cebc
md"""

Then the following cell defines a source panel, plots the potential contours and velocity vectors created by the presence of the panel on a predefined bachground mesh by `meshgrid( )` function. 

"""

# ╔═╡ f4ac8540-3229-4298-8065-3b193bcc8a96
begin
    # Define line segment
    x0, y0 = -0.1, -0.5
    x1, y1 = 0.4, 0.17

	x_range_p = range(-2, 2, length=22)
	y_range_p = range(-2, 2, length=22)
	XY_p = meshgrid(x_range_p, y_range_p)

    # Calculate velocity components and potential
    UV = velocity.(XY_p..., x0, y0, x1, y1, source)
    ϕ = potential.(XY_p..., x0, y0, x1, y1, source)
	
	# Scale the plot
	magnitudes = norm.(UV)
	valid_magnitudes = filter(!isnan, magnitudes)
	global_scale = norm(valid_magnitudes) 
	
	UV = 5*(UV ./ global_scale)

    # Extract U and V matrices using Julia's vectorized getindex
	# This pulls the 1st (u) and 2nd (v) element from every tiny array in the grid
    u = getindex.(UV, 1)
    v = getindex.(UV, 2)
	
    # Decompose coordinates needed for the quiver plot
	X,Y = XY_p
    # Plot the contour lines of the potential
    p = contour(x_range_p, y_range_p, ϕ, levels=20, aspect_ratio=:equal, size= 				(700,700), legend=false)

    # Overlay the velocity vectors
    # vec() flattens the 2D matrices into 1D lists for the quiver function
    quiver!(p, vec(X), vec(Y), quiver=(vec(u), vec(v)), color=:black)
    # 3. Draw the panel line on top
    plot!(p, [x0, x1], [y0, y1], color=:blue, linewidth=3)
    # Display the plot
    p
end

# ╔═╡ db133926-612c-4fbd-99e3-835c3974a487
md"""
The resulting flow resembles point-source flow when $r\gg ds$, which is expected because the Green's function for a point vortex is proportional to the Green's function for a point-source by a constant.
"""

# ╔═╡ a31645ba-2e8c-4930-811d-ed912a2ffab7
md"""
### 5.2.2 Superposition of Multiple Panels

The next step is to add multiple plates together to form a geometry. The potential and velocity can be lineary superimposed as follows:

$\phi=\sum_{i}\phi_{i}=\sum_{i}q_{1}F^{\phi}_{i}, \quad \vec{u}=\sum_{i}\vec{u}_{i}=\sum_{i}q_{i}F^{\vec{u}}_{i}$

Then after defining a set of $N$ panels by connecting $N+1$ points to form a square, the potential around this square would simply be the sum of the potential over each panel in range $N$. The following code gives a generalized function to calculate the potential around any geometry defined by $N+1$ points and plots the resulting velocity field vectors on the cartesian mesh grid defined before.  
"""

# ╔═╡ 72a93620-14a6-4552-a44c-0040f554032d
begin
    function plot_flow(x, y, q, XY, x_ran, y_ran; G=source) # Function to plot 			potential contours and velocity vectors
        # The number of panels is the number of node points minus 1
        N = length(x) - 1
        
        # SUPERPOSITION: Calculate and sum the scaled velocity grids from all panels
        UV = sum(q[i] .* velocity.(XY..., x[i], y[i], x[i+1], y[i+1], G) for i in 			1:N)
		Gamma = sum(q[i] .* potential.(XY..., x[i], y[i], x[i+1], y[i+1], G) for i in 1:N)
				
        # NORMALIZATION: Keep arrows manageable
        # norm.(UV) gets the length of each [u, v] vector, and the outer norm() 			combines them
		magnitudes = norm.(UV)
        valid_magnitudes = filter(!isnan, magnitudes)
        global_scale = norm(valid_magnitudes) 
        UV = 5*(UV ./ global_scale)
		
		# Extract U and V matrices for the quiver plot
        u = getindex.(UV, 1)
        v = getindex.(UV, 2)
		
		# Decompose coordinates needed for the quiver plot
        X,Y = XY
        
        # Plot
        p = plot(aspect_ratio=:equal, size=(700, 700), legend=false)
        contour!(p, x_ran, y_ran, Gamma, levels=10)
        quiver!(p, vec(X), vec(Y), quiver=(vec(u), vec(v)), color=:black)
        plot!(p, x, y, color=:blue, linewidth=3, marker=:circle)
        return p
    end
	plot_flow([-1, 1, 1, -1], [-1, -1, 1, 1], [1, -1, 1], XY_p, x_range_p, y_range_p; G=source)	
end

# ╔═╡ c4372b8d-759c-48fb-a79d-31c1b0113840
md"""
### 5.2.3 Extending to Closed Geometries:

Before extending to the foil geometry there are some sanity checks that has to be carried out on closed bodies. Essentially the convergence of the code has to be checked for different geometries and for both source and vortex panels. For this purpose the folloving cell introdusec body functions for an ellipse and a circle.

"""

# ╔═╡ 55fe8cad-ccd7-444f-89e7-4ae13a2b43e8
begin
    function circle(N, R=1.0)
        # N panels require N+1 node points to close the shape
        theta = range(-π, π, length=N+1)
        
        # Calculate X and Y coordinates
        x = R .* cos.(theta)
        y = R .* sin.(theta)
        return x, y
    end
	
    function ellipse(N, a=1.0, b=1.0, theta1=π)
        # N panels require N+1 node points to close the shape
        theta = range(-π, theta1, length=N+1)
        
        # Calculate X and Y coordinates using broadcasting
        x = a .* cos.(theta)
        y = b .* sin.(theta)
        return x, y
    end
end

# ╔═╡ 57c5fd39-fe40-4901-b3b7-f1c4c9b126fb
md"""

### 5.2.4. Assembling and Solving the Matrix Equation System

Green's function is a solution to Laplace equation $\nabla^2 G=0$ meaning that any superposition of vortex strengths will result in a valid potential flow. Then the vortex strengths can be solved for usin the normal velocity boundary condition, which is given as follows on the body:

$\sum_{i}\vec{u}_{i}\cdot\hat{p}=\vec{U}\cdot\hat{p}$

Where $\hat{p}$ is the panel normal $(\hat{n})$ or tangent $(\hat{t})$ vector depending on if the body is modeled by source or vortex panels respectively, and $\vec{U}$ is the panel velocity. Substituting the equation for $\vec{u}_{i}$ and simming over each panel $j$ yields:

$\sum_{i=1}^{N}\lambda_{i}F^{\vec{u}}_{i}({c_{j}})\cdot\vec{p}=\vec{U}\cdot\vec{p}_{j}$

Where $\vec{c}_{j}$ is the collocation point (midpoint in this case) of panel $j$. defining $a_{ij}=F^{\vec{u}}_{i}(\vec{c}_{j})\cdot\vec{p}_{j}$ as the normal (for source) or tangential velocity (for vortex) influence of panel $i$ on $j$ and $b_{j}=\vec{U}\cdot\vec{p}_{j}$ yields:

$\sum_{i=1}^{N}a_{ij}\lambda_{i}=b_{j}$

This results in a linear equation system which can bw written in matrix form as $[A][\Lambda]=[B]$. Then the system can be solved for unknown source or vortex strengths $[\lambda]$. The code up to this point assumed unit source/vortex strength ($\lambda_{1}=1$). The following cell introduces a freestream velocity (`U_free`) of $U_{\infty}=(1,0)$ and then introduces the functions to calculate the individual source/vortex strength for each panel.

"""

# ╔═╡ 3d2f3a80-713a-4e04-b028-c698ef9783db
begin
    function CI_CCW(x0, y0, x1, y1) # Function to calculate panel properties for 		arrays created woth a CCW loop. 
        sx, sy = x1 .- x0, y1 .- y0          
        xc, yc = x0 .+ 0.5 .* sx, y0 .+ 0.5 .* sy  
        h = sqrt.(sx.^2 .+ sy.^2) 
        nx, ny = sy ./ h, -sx ./ h
		tx, ty = -sx ./ h, -sy ./ h 
        return xc, yc, nx, ny, tx, ty, h
    end

	function CI_CW(x0, y0, x1, y1) # Function to calculate panel properties for 	arrays created woth a CW loop. 
        sx, sy = x1 .- x0, y1 .- y0          
        xc, yc = x0 .+ 0.5 .* sx, y0 .+ 0.5 .* sy  
        h = sqrt.(sx.^2 .+ sy.^2)             
        nx, ny = -sy ./ h, sx ./ h
		tx, ty = sx ./ h, sy ./ h  
        return xc, yc, nx, ny, tx, ty, h
    end
    
    function construct_A(x, y; G=source, kutta=false, loopdir=:CCW) # Function to 		costruct [A]. 
        x0, y0 = x[1:end-1], y[1:end-1] 
        x1, y1 = x[2:end], y[2:end] 
        
        # The dorection router
        if loopdir === :CCW
            xc, yc, nx, ny, tx, ty, _ = CI_CCW(x0, y0, x1, y1) 
        elseif loopdir === :CW
            xc, yc, nx, ny, tx, ty, _ = CI_CW(x0, y0, x1, y1)
        else
            error("Invalid direction. Please use dir=:CCW or dir=:CW")
        end
        
        N = length(xc)
        A = zeros(N, N) 
        
        # Calculate panel influences
        for i in 1:N
            UV = velocity.(xc, yc, x[i], y[i], x[i+1], y[i+1], G)
            u = getindex.(UV, 1)
            v = getindex.(UV, 2)
            
            if G === source
                A[:, i] = u .* nx .+ v .* ny
            elseif G === vortex
                A[:, i] = u .* tx .+ v .* ty
            end
        end
        
        # Replace the diagonals
        if G === source
            for i in 1:N
                A[i, i] = π
            end
        elseif G === vortex
            for i in 1:N
                A[i, i] = π
            end
        end
        
        # Kutta condition
        if kutta 
            A[end, :] .= 0.0         
            A[end, 1] = 1.0          
            A[end, end] = 1.0        
        end
        
        return A, nx, ny, tx, ty, N
    end
    
    function solve(x, y, U::Tuple; G=source, kutta=false, loopdir=:CCW) # Function 		to solve the equation system
        Ux, Uy = U
        A, nx, ny, tx, ty, N = construct_A(x, y; G=G, kutta=kutta, loopdir=loopdir)
        
        # Build B Vector
        if G === source
            B = -(Ux .* nx .+ Uy .* ny)
        elseif G === vortex
            B = -(Ux .* tx .+ Uy .* ty)
        end
		# Match Kutta condition
        if kutta
            B[end] = 0.0 
        end
        return A \ B
    end
end

# ╔═╡ d7f5d921-5f0f-415c-908d-70461bf08273
md"""

### 5.2.5. Sanity Checks on the Implementation of Source/Vortex Methods

The folloving cells solves for the uniform flow around a closed circle with both source and vortex panel panel method and plots the velocity magnitude field to show that both methods are capable of obtaining the same flow field without any problems. This is the initial sanity check to determine if the matrix equation system is set up properly. 

"""

# ╔═╡ 9b91d7ca-9f35-4a9b-ba0a-2f433589ba3e
md"""

The source and vortex panels can both estimate the potential flow around a circle evident by the location of th estagnation points and the shape of the high and low velocity field contours.

"""

# ╔═╡ e7c5d436-1879-4c09-bd3c-fb3e83660d96
md"""

The folloving cells calculate the uniform flow around an ellipse with both source an vortex panels and plot the same velocity field magnitude contours to check if the code can handle different shapes.   

"""

# ╔═╡ a19b1f4b-1710-4df3-b7ec-2575005e16dd
md"""

It is evident from the velocity magnitude fields that the source and vortex panels can both estimate the potential flow around a different geometries.

"""

# ╔═╡ dd71b842-fd92-4957-a299-234b746f0c7c
md"""

## 5.3. Solution Convergence

Checking the convergence of the solution is one of the important steps of validating a numerical tool. For this case at least the source or vortex strengths has to be checked for convergence. The following code introduces the functions to check for source/vortex strength convergence for non-lifting geometries and for lifting geometries the convergence of integrated lift and Kutta-Jukowski lift is also included. The following cells contain the functions needed to plot convergence and presents the convergence plots for static circle and NACAXXX foil cases for $N=[16,32,64,256]$.

"""

# ╔═╡ eb4fad3a-2200-4d2f-9d58-6ac636c51c74
begin
    # 1. The Generalized Convergence Function
    # We pass 'body_func' (the geometry generator) and 'U' (the freestream)
    # We set default keyword arguments for the panels and the kernel!
    function plot_convergence(body_func, U; G=source, kutta=false, N_vals=[16, 32, 64, 256], loopdir=:CCW)
		if G === source
	        p = plot(
	            xlabel="Normalized Panel Number ([0] Start, [1] end of the panel loop)", 
	            ylabel="q", 
	            title="Source Strength Convergence",
	            size=(700, 500),
	            legendtitle="N",
	            legend=:bottomright )
		elseif G === vortex
			p = plot(
	            xlabel="Normalized Panel Number ([0] Start, [1] end of the panel loop)", 
	            ylabel="γ", 
	            title="Vortex Strength Convergence",
	            size=(700, 500),
	            legendtitle="N",
	            legend=:bottomright )
		end
        for N in N_vals
            # 1. Generate the shape using whatever function the user passed!
            x_nodes, y_nodes = body_func(N)
            
            # 2. Solve for strengths
            q_vals = solve(x_nodes, y_nodes, U, G=G, kutta=kutta, loopdir=loopdir)
            
            # 3. Generate the x-axis angles
            theta = (1:N) ./ (N+1)
            
            # 4. Add to plot
            plot!(p, theta, q_vals, label="$N", linewidth=2)
        end
        return p
    end
end

# ╔═╡ 5343fda3-9135-40bb-9741-e7e3ad84de1c
md"""

### 5.3.1. Source Strength Convergence for Static Circle

"""

# ╔═╡ 3a236d10-7588-4e67-86d9-e6ee41c64a24
md"""

### 5.3.2. Vortex Strength Convergence for Static Circle

"""

# ╔═╡ 035f5a34-65a8-4b5f-b9a1-4d0daf1d52c0
md"""

### 5.3.4. Vortex Strength and Force Coefficient Convergence for a NACAXXXX foil 

The last step is to extend to the NACA $4$-digit foil geometry and verify the implementation of the Kutta condition is done properly. The next cells calculates the convergence of vortex strengths and lift coefficient around the user defined NACAXXXX at the user defined flow condition from the beginning of this Pluto notebook. This part is directly connected to the user defined values so that the stability of the code can be checked by each user at ease.  

"""

# ╔═╡ 78352abf-eb80-4e3c-a8e1-863ba6f94016
md"""

## 5.4. Post-Processing

The following cells contain the functions to plot the flow field around an arbitrary body and calculationg the lift for a NACAXXXX foil during the post-processing step.

"""

# ╔═╡ af49410f-95d4-4d93-a41f-bf5b1b7b7885
begin
    function plot_flow_P(x, y, q, XY, x_ran, y_ran, U_free; G=source, plot_var=:magnitude) # Function to plot flow fields
        Ux, Uy = U_free
        X, Y = XY
        
        # The number of panels is the number of node points minus 1
        N = length(x) - 1

        # Calculate and sum the scaled velocity grids
        UV_panel = sum(q[i] .* velocity.(XY..., x[i], y[i], x[i+1], y[i+1], G) for i in 1:N)
        UV = UV_panel .+ Ref([Ux, Uy]) 

        # Extract u, v, and magnitudes
        u = getindex.(UV, 1)
        v = getindex.(UV, 2)
        magnitudes = norm.(UV)
        V_inf = norm([Ux, Uy])
        
        # Calculate the entire Pressure field using Bernoulli! <---
        Cp_field = 1.0 .- (magnitudes ./ V_inf).^2

        # Select the data and appropriate bounds
        if plot_var === :u
            plot_data = u
            plot_title = "X-Velocity Field (\$u_{x}\$)"
            c_lims = (-0.5 * V_inf, 2.0 * V_inf) 
        elseif plot_var === :magnitude
            plot_data = magnitudes
            plot_title = "Velocity Magnitude Field (\$U\$)"
            c_lims = (-0.5 * V_inf, 2.0 * V_inf) 
        elseif plot_var === :pressure
            plot_data = Cp_field
            plot_title = "Pressure Field (\$C_{p}\$)"
            c_lims = (-3.0, 1.0)  # Max Cp is 1.0, minimum can be highly negative
        end
        
        p = plot(aspect_ratio=:equal, size=(750, 750), legend=true, title=plot_title)
        
        contourf!(p, x_ran, y_ran, plot_data, levels=25, color=:bam, alpha=0.8, clims=c_lims)
        
        # Draw the physical panel segments (Nodes hidden from legend)
        plot!(p, x, y, color=:blue, linewidth=3, marker=:circle, label="")
        return p
    end
end

# ╔═╡ be57e391-4d00-4373-a938-06b6c25bae97
begin
    function forces(x, y, γ, U; loopdir=:CW) # Function to calculate forces
        Ux, Uy = U
        V_inf = hypot(Ux, Uy)
		chord = maximum(x) - minimum(x)
        
        # 1. Grab panel properties
        x0, y0 = x[1:end-1], y[1:end-1] 
        x1, y1 = x[2:end], y[2:end] 
        if loopdir === :CCW
            xc, yc, nx, ny, tx, ty, h = CI_CCW(x0, y0, x1, y1) 
        elseif loopdir === :CW
            xc, yc, nx, ny, tx, ty, h = CI_CW(x0, y0, x1, y1)
        else
            error("Invalid direction")
        end
        
        # Calculate Pressure Coefficient (Cp)
        Cp = 1.0 .- ((γ.*2π) ./ V_inf).^2
        
        # Integrate X and Y forces
        Cx = sum(-Cp .* nx .* h) / (chord * 2π)
        Cy = sum(-Cp .* ny .* h) / (chord * 2π)
        
        # Calculate Total Circulation for Kutta-Joukowski check
        Gamma_total = sum(γ .* h)
		Cl_KJ = -(2.0 * Gamma_total) / chord
		
        println("--- Force Results ---")
	    println("Integrated Lift (CL) : ", round(Cy, digits=5))
	    println("Kutta-Joukowski (CL) : ", round(Cl_KJ, digits=5))
        return Cy, Cx, Cp, Gamma_total
    end
end

# ╔═╡ 18dabe8c-cd8f-49ee-95e2-f8c12f2c5a9b
begin
	U_free = (1.0,0.0) 
	x_foil, y_foil = NACA(N)
	
	# Set grid resolution
	x_range = range(-1, 3, length=60)
	y_range = range(-2, 2, length=60)
	
	XY = meshgrid(x_range, y_range)
	XY_foil = mask_grid(XY, x_foil, y_foil)
	γ = solve(x_foil,y_foil,U_free; G=vortex, kutta=true, loopdir=:CW)
	forces(x_foil, y_foil, γ , U_free; loopdir=:CW)
	plot_flow_P(x_foil, y_foil, γ, XY_foil, x_range, y_range, U_free; G=vortex, plot_var=:magnitude)
end

# ╔═╡ d464184d-f168-4e18-b736-40b86d33dc02
begin
	x_circ, y_circ = circle(32,0.5)
	x_range_circ = range(-2, 2, length=60)
	y_range_circ = range(-2, 2, length=60)
	XY_circ = meshgrid(x_range_circ, y_range_circ)
	XY_circ = mask_grid(XY_circ, x_circ, y_circ)
	q_circ = solve(x_circ,y_circ,U_free; G=source)
	plot_flow_P(x_circ, y_circ, q_circ, XY_circ, x_range_circ, y_range_circ, U_free; G=source)
end

# ╔═╡ 11731f10-b790-420b-a9b5-abab65960824
begin
	γ_circ = solve(x_circ,y_circ,U_free; G=vortex)
	plot_flow_P(x_circ, y_circ, γ_circ, XY_circ, x_range_circ, y_range_circ, U_free; G=vortex)
end

# ╔═╡ fcc68ef4-8b34-4483-ae60-3aaf99ef1b42
begin
	x_el, y_el = ellipse(32,1,0.5)
	XY_el = meshgrid(x_range_circ, y_range_circ)
	XY_el = mask_grid(XY_el, x_el, y_el)
	q_el = solve(x_el,y_el,U_free; G=source)
	plot_flow_P(x_el, y_el, q_el, XY_el, x_range_circ, y_range_circ, U_free; G=source)
end

# ╔═╡ 6d1083ea-203f-4516-a879-2859e54a0b69
begin
	γ_el = solve(x_el,y_el,U_free; G=vortex)
	plot_flow_P(x_el, y_el, γ_el, XY_el, x_range_circ, y_range_circ, U_free; G=vortex)
end

# ╔═╡ 10ae44b6-f0b4-4210-a20a-a31940ecde41
begin
	plot_convergence(N -> circle(N), U_free, G=source, kutta=false, loopdir=:CCW)
end

# ╔═╡ 04751d3e-1b28-45af-a308-7c87a5e162ef
begin
	plot_convergence(N -> circle(N), U_free, G=vortex, kutta=false, loopdir=:CCW)
end

# ╔═╡ 0068d83c-98dd-4449-b697-0d1a5b4b097f
begin
    function plot_convergence_full(body_func, U; G=source, kutta=false, N_vals=[16, 32, 64, 128, 256], loopdir=:CCW) # Function to plot the convergence of source/vortex strengths and body foreces
        if G === source
            p_strength = plot(
                xlabel="Normalized Panel Number", 
                ylabel="q", 
                title="Source Strength Convergence",
                legendtitle="N",
                legend=:bottomright )
        elseif G === vortex
            p_strength = plot(
                xlabel="Normalized Panel Number",
                ylabel="γ", 
                title="Vortex Strength Convergence",
                legendtitle="N",
                legend=:bottomright )
        end

        # Arrays to hold our force data
        Cl_vals = Float64[]
        Cl_KJ_vals = Float64[] # New array for KJ Lift
        Cd_vals = Float64[]

        # Pre-calculate freestream velocity for KJ math
        Ux, Uy = U
        V_inf = hypot(Ux, Uy)

        # Convergence loop
        for N in N_vals
            # Generate the shape
            x_nodes, y_nodes = body_func(N)
            chord = maximum(x_nodes) - minimum(x_nodes)
            
            # Solve for strengths
            q_vals = solve(x_nodes, y_nodes, U, G=G, kutta=kutta, loopdir=loopdir)
            
            # Calculate forces (using Gamma_total to get Cl_KJ)
            Cl, Cd, Cp, Gamma_total = forces(x_nodes, y_nodes, q_vals, U; loopdir=loopdir)
            
            # Reconstruct Cl_KJ (Standard textbook formula)
            Cl_KJ = -(2.0 * Gamma_total) / (V_inf * chord)
            
            # Store the forces for plotting
            push!(Cl_vals, Cl)
            push!(Cl_KJ_vals, Cl_KJ)
            push!(Cd_vals, Cd)
            
            # Add the strength line to the first plot
            theta = (1:N) ./ (N+1)
            plot!(p_strength, theta, q_vals, label="$N", linewidth=2)
        end
        
        # Plot Integrated Lift (Cl) on the Primary (Left) Axis
        p_forces = plot(N_vals, Cl_vals, 
            xlabel="Number of Panels (N)", 
            ylabel="Lift (\$C_{L}\$)", 
            title="Force Convergence",
            marker=:circle,
            linewidth=2,
            label="Integrated \$C_{L}\$",
            color=:blue,
            yguidefontcolor=:blue, 
            legend=:topleft)

        # Add Kutta-Joukowski Lift to the SAME Primary Axis
        plot!(p_forces, N_vals, Cl_KJ_vals,
            marker=:square,      # Different marker
            linestyle=:dash,     # Dashed line to separate from Integrated Cl
            linewidth=2,
            label="K-J \$C_{L}\$",
            color=:cyan)         # Lighter blue to group visually with Lift

        # Plot Drag (Cd) on the Secondary (Right) Axis using twinx()
        p_twin = twinx(p_forces)
        plot!(p_twin, N_vals, Cd_vals, 
            ylabel="Drag (\$C_{D}\$)", 
            marker=:diamond,
            linewidth=2,
            label="Integrated \$C_{D}\$",
            color=:red,
            yguidefontcolor=:red,
            legend=:topright)

        final_plot = plot(p_strength, p_forces, layout=(1, 2), size=(800, 400), margin=6Plots.mm)
        return final_plot
    end
end

# ╔═╡ 691b0fb1-9881-47d5-be6f-3b2238a9fc17
begin
	p_full = plot_convergence_full(N -> NACA(N), U_free, G=vortex, kutta=true, loopdir=:CW)
end

# ╔═╡ 8f42b970-e661-4c11-81fb-88c0f7ea10ff
md"""

## 5.5. Validation (Future Work)

"""

# ╔═╡ 86ac31c5-7948-4c3e-bcc7-2a7e6c289850
md"""

# References

Airfoil Tools, n.d. NACA 0012 AIRFOILS. [Online] 
Available at: http://airfoiltools.com/airfoil/details?airfoil=n0012-il
[Accessed 24 03 2024].

Anderson, J., 2017. Chapter 4: Incompressible Flow over Airfoils. In: Fundamentals of Aerodynamics. s.l.:McGraw-Hill, pp. 321-423.

Bruna, P. M., 2011. Engineering the race car wing: application of the vortex panel numerical method. Sports Engineering, Volume 13, pp. 195-204.

Byrne, G., Persoons, T. & Kingston, W., 2019. Experimental validation of lift and drag forces on an asymmetrical hydrofoil for seafloor anchoring applications. Journal of Ocean and Climate, Volume 9, pp. 1-11.

Chen, Z. M., 2012. A vortex based panel method for potential flow simulation around a hydrofoil. Journal of Fluids and Structures, Volume 28, pp. 378-391.

Douglass, A., 2024. 2D Panel Methods. [Online] 
Available at: https://www.aerodynamics4students.com/subsonic-aerofoil-and-wing-theory/2d-panel-methods.php
[Accessed 28 03 2024].

Jafarimoghaddam, A. & Aberoumand, S., 2016. Introducing an Optimized Airfoil Shape Using Panel Method: A Short Report. European Journal of Advances in Engineering and Technology, , 3(7), pp. 47-52.

Ketz, J., & Plotkin, A. (1991). Low-speed aerodynamics: from wing theory to panel methods. McGraw-Hill, Incorporated.

Ladson, C. & Brooks, W., 1975. Development of a Computer Program to Obtain Oordinates for a NACA 4-digit, 4-dgit modified, 5-digit and 16-series airfoils, Hampton, Va: s.n.

Leishman, J. G. & Galbraith, R. A., 1985. G.U. Aero Report 8502 - An Algorithm For The Calculation Of The Potential Flow About An Arbitrary Two Dimensional Aerofoil, s.l.: s.n.

Louli, Z. & Winarto, H., 2007. Airfoil Analysis using First and Second Order Vortex Panel Methods with Neumann Boundary Conditions, s.l.: s.n.

Samarghandi, R., 2021. Vortex Panel Method. [Online] 
Available at: https://github.com/rezasamarghandi/vortex-panel-method
[Accessed 24 03 2024].

Weymouth, G., 2021. Solve potential flow problems using vortex panels. [Online] 
Available at: https://github.com/weymouth/MarineHydro/tree/master/pre2021/vortexpanel
[Accessed 28 03 2024].

Weymouth, G., 2024. Numerical Ship Hydro. [Online] 
Available at: https://github.com/weymouth/NumericalShipHydro
[Accessed 20 03 2024].

"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
OutMacro = "0ae4d431-9932-4135-a8f1-51ee5e017775"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
BenchmarkTools = "~1.5.0"
GeometryBasics = "~0.5.10"
Interpolations = "~0.16.2"
OutMacro = "~0.1.0"
Plots = "~1.41.1"
PlutoUI = "~0.7.58"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.6"
manifest_format = "2.0"
project_hash = "2c659b7245c232fb0c451ec1ac2c8fbe25582cff"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "0761717147821d696c9470a7a86364b2fbd22fd8"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.5.2"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "f1dff6729bc61f4d49e140da1af55dcd1ac97b2f"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.5.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d0efe2c6fdcdaa1c161d206aa8b933788397ec71"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.6+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "12177ad6b3cad7fd50c8b3825ce24a99ad61c18f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "21d088c496ea22914fe80906eb5bce65755e5ec8"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.1"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e86f4a2805f7f19bec5129bc9150c38208e5dc23"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.4"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "66381d7059b5f3f6162f28831854008040a4e905"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "70329abc09b886fd2c5d94ad2d9527639c421e3e"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.14.3+1"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "9e0fb9e54594c47f278d75063980e43066e26e20"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.1+1"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "44716a1a667cb867ee0e9ec8edc31c3e4aa5afdc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.24"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "be8a1b8065959e24fdc1b51402f39f3b6f0f6653"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.24+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "IterTools", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "1f5a80f4ed9f5a4aada88fc2db456e637676414b"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.10"

    [deps.GeometryBasics.extensions]
    GeometryBasicsGeoInterfaceExt = "GeoInterface"

    [deps.GeometryBasics.weakdeps]
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "24f6def62397474a297bfcec22384101609142ed"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.3+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "51059d23c8bb67911a2e6fd5130229113735fc7e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.11.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c0c9b76f3520863909825cbecdef58cd63de705a"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.5+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "17b94ecafcfa45e8360a4fc9ca6b583b049e4e37"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.1.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cc3ad4faf30015a3e8094c9b5b7f19e85bdf2386"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.42.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d620582b1f0cbe2c72dd1d5bd195a9ce73370ab1"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.42.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "8785729fa736197687541f7053f6d8ab7fc44f92"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.10"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ff69a2b1330bcb730b9ac1ab7dd680176f5896b8"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.1010+0"

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.OutMacro]]
deps = ["Test"]
git-tree-sha1 = "e373d6a38722569f9b374be4a3f884802bb117ed"
uuid = "0ae4d431-9932-4135-a8f1-51ee5e017775"
version = "0.1.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "12ce661880f8e309569074a61d3767e5756a199f"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.1"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "71a22244e352aa8c5f0f2adde4150f62368a3f2e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.58"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Profile]]
deps = ["StyledStrings"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "d7a4bff94f42208ce3cf6bc8e4e7d1d663e7ee8b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.10.2+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll", "Qt6Svg_jll"]
git-tree-sha1 = "d5b7dd0e226774cbd87e2790e34def09245c7eab"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.10.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "4d85eedf69d875982c46643f6b4f66919d7e157b"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.10.2+1"

[[deps.Qt6Svg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "81587ff5ff25a4e1115ce191e36285ede0334c9d"
uuid = "6de9746b-f93d-5813-b365-ba18ad4a9cf3"
version = "6.10.2+0"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "672c938b4b4e3e0169a07a5f227029d4905456f2"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.10.2+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "246a8bb2e6667f832eea063c3a56aef96429a3db"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.18"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b29c22e245d092b8b4e8d3c09ad7baa586d9f573"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.3+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "0ba01bc7396896a4ace8aab67db31403c71628f4"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.7+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c174ef70c96c76f4c3f4d3cfbe09d018bcd1b53"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "4909eb8f1cbf6bd4b1c30dd18b2ead9019ef2fad"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.18.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "ed756a03e95fff88d8f738ebc2849431bdd4fd1a"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.2.0+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "63aac0bcb0b582e11bad965cef4a689905456c03"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.125+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "45a20bd63e4fafc84ed9e4ac4ba15c8a7deff803"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.57+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "a1fc6507a40bf504527d0d4067d718f8e179b2b8"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.13.0+0"
"""

# ╔═╡ Cell order:
# ╟─3af5e095-c9c0-4ec3-ae02-c2e060b6c153
# ╟─cc29fa3e-cfba-4de0-ba7e-451a33d8714d
# ╟─1c44dba8-ed15-4908-a3ba-968a7ca88d8b
# ╟─6de9ba9e-776c-415b-9d9b-1083e6e1540a
# ╟─eed1913e-57fd-42c0-9430-2f676d059efc
# ╠═c14bcbc0-eb5d-11ee-1bf8-8fcaaf8e34c2
# ╟─37cc0cb5-ee73-4791-b787-c0bc145d225c
# ╟─e7624d6e-17aa-4dfd-b677-bed2116ae0e3
# ╟─af51eea8-1423-4720-96c3-837d12406fea
# ╟─b92edf15-fa1f-4930-ba73-7c73e7363055
# ╟─a90a1cf8-2fe5-420d-9699-96bda5b00988
# ╟─0445668c-3559-4159-957c-df9e5aa40c76
# ╟─157f9c00-4ff3-4bc6-ade1-c78886042128
# ╟─740e3ad2-4397-4afa-8db1-cb2718baa1e3
# ╟─18dbf7f0-e90b-45e5-8de8-82db0abb6b9e
# ╟─c7addb5a-0a05-423b-9f73-0e237d189d11
# ╟─d155e9a4-2ae6-4229-a810-633d49df1354
# ╟─797534ac-b52a-4e58-812d-0a7962615590
# ╟─9eea78c8-3233-44d2-b375-225bce7ebc8f
# ╟─f1ad33b3-6f11-4a01-a83a-9b2ebe3ca382
# ╟─44686d98-0659-4157-bd34-c932b7e78d61
# ╠═18dabe8c-cd8f-49ee-95e2-f8c12f2c5a9b
# ╟─6b53a7e2-692b-4093-8afb-2c07f4325e5c
# ╟─07da2b45-6d9b-4883-8b8f-4625405018ac
# ╟─4753ed02-0ad1-4971-8fef-219d2862b85f
# ╟─18416a19-0637-4de4-9c86-daae967f443b
# ╟─37a088a7-5552-460b-b01a-ad0e521d7919
# ╟─4b780945-cce8-46ea-a28f-b4c2bfd603c6
# ╟─6b332adc-359f-40d7-8396-273b35d23550
# ╟─87345c05-f3c2-4225-b44d-b68d6e34509b
# ╟─909de60b-8906-47c6-b1fb-5f0d3c104cbc
# ╠═a11a44b8-20cb-4e6b-a22a-0408b47146f6
# ╟─00565f60-a92d-4d32-85b0-a4d5c63b5f1e
# ╟─6f76e9d9-2dfa-4b94-9d49-d5acdcc70176
# ╟─73a9aae0-2b30-4baa-80b1-028d6bc1cebc
# ╠═f4ac8540-3229-4298-8065-3b193bcc8a96
# ╟─db133926-612c-4fbd-99e3-835c3974a487
# ╟─a31645ba-2e8c-4930-811d-ed912a2ffab7
# ╠═72a93620-14a6-4552-a44c-0040f554032d
# ╟─c4372b8d-759c-48fb-a79d-31c1b0113840
# ╟─55fe8cad-ccd7-444f-89e7-4ae13a2b43e8
# ╟─57c5fd39-fe40-4901-b3b7-f1c4c9b126fb
# ╟─3d2f3a80-713a-4e04-b028-c698ef9783db
# ╟─d7f5d921-5f0f-415c-908d-70461bf08273
# ╠═d464184d-f168-4e18-b736-40b86d33dc02
# ╟─9b91d7ca-9f35-4a9b-ba0a-2f433589ba3e
# ╠═11731f10-b790-420b-a9b5-abab65960824
# ╟─e7c5d436-1879-4c09-bd3c-fb3e83660d96
# ╠═fcc68ef4-8b34-4483-ae60-3aaf99ef1b42
# ╟─a19b1f4b-1710-4df3-b7ec-2575005e16dd
# ╠═6d1083ea-203f-4516-a879-2859e54a0b69
# ╟─dd71b842-fd92-4957-a299-234b746f0c7c
# ╟─eb4fad3a-2200-4d2f-9d58-6ac636c51c74
# ╟─0068d83c-98dd-4449-b697-0d1a5b4b097f
# ╟─5343fda3-9135-40bb-9741-e7e3ad84de1c
# ╠═10ae44b6-f0b4-4210-a20a-a31940ecde41
# ╟─3a236d10-7588-4e67-86d9-e6ee41c64a24
# ╠═04751d3e-1b28-45af-a308-7c87a5e162ef
# ╟─035f5a34-65a8-4b5f-b9a1-4d0daf1d52c0
# ╠═691b0fb1-9881-47d5-be6f-3b2238a9fc17
# ╟─78352abf-eb80-4e3c-a8e1-863ba6f94016
# ╟─af49410f-95d4-4d93-a41f-bf5b1b7b7885
# ╟─be57e391-4d00-4373-a938-06b6c25bae97
# ╟─8f42b970-e661-4c11-81fb-88c0f7ea10ff
# ╟─86ac31c5-7948-4c3e-bcc7-2a7e6c289850
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
