# Mini-Boson Stars: Numerical Analysis

## Project Overview

This project delves into the theoretical and computational aspects of Mini-Boson Stars, hypothetical compact astrophysical objects composed entirely of scalar bosons. Utilizing numerical methods, it solves the coupled Einstein-Klein-Gordon (EKG) equations to model spherically symmetric boson star configurations. The primary objectives are to analyze the internal structure, spacetime geometry, mass distribution, energy density, and pressure profiles of these stars. The project also investigates the relationship between the star's total mass (ADM mass) and the scalar field's oscillation frequency, and provides a robust framework for converting dimensionless simulation results into physically meaningful units for real-world astrophysical interpretation.

## Files

### `Mini-Boson Stars.ipynb`

*   **Purpose:** This Jupyter notebook serves as the central computational and analytical tool for simulating Mini-Boson Stars. It implements the numerical solution of the Einstein-Klein-Gordon equations for a self-gravitating, spherically symmetric complex scalar field and subsequently processes, visualizes, and interprets the results.

*   **Methods:**
    *   **Einstein-Klein-Gordon (EKG) System:** The notebook sets up the coupled system of four first-order ordinary differential equations (ODEs) derived from Einstein's field equations and the Klein-Gordon equation for a complex scalar field $\phi(r,t) = \phi_0(r) e^{i\omega t}$. This system describes a spherically symmetric spacetime with a metric $ds^2 = -e^{2\Phi(r)} dt^2 + e^{2\Lambda(r)} dr^2 + r^2 d\Omega^2$.
    *   **Boundary Value Problem (BVP) Solver:** It employs `scipy.integrate.solve_bvp` to numerically solve this system of ODEs. The solver enforces specific boundary conditions at the stellar center ($r=0$) and at asymptotic infinity ($r \to \infty$), ensuring physical regularity and asymptotic flatness.
    *   **Dimensional Analysis and Unit Rescaling:** The notebook includes detailed functions (`convert_to_physical_units`, `convert_to_physical_units_df`) designed to convert the dimensionless results obtained from the numerical simulations into physically meaningful SI units (e.g., kilometers, solar masses, J/m$^3$). This conversion is based on fundamental constants and the specified mass of the scalar boson.

*   **Main Logic:**
    1.  **Define ODE System (`system` function):** The four coupled first-order ODEs for the gravitational potential $\Phi(r)$, enclosed mass $m(r)$, scalar field amplitude $\phi_0(r)$, and its derivative $d\phi_0/dr$ are defined. The dimensionless scalar particle mass (denoted `m=1` in the code) is incorporated here.
    2.  **Specify Boundary Conditions (`bc` function):** Physical constraints are imposed at the origin (e.g., $\phi_0(0) = \phi_c$, $m(0)=0$, $d\phi_0/dr(0) = 0$, $\Phi(0)=0$) and at infinity (e.g., $\phi_0(\infty) = 0$, $\Phi(\infty)=0$), which are crucial for obtaining stable, physical solutions.
    3.  **Initial Guess and Grid:** An initial guess for the field profiles and a suitable radial grid are prepared to assist the BVP solver in converging to a solution.
    4.  **Solve BVP:** The `solve_bvp` function is invoked to find the numerical solution for a given scalar field oscillation frequency $\omega$.
    5.  **Calculate Derived Quantities:** From the solved profiles of $\Phi(r)$, $\Lambda(r)$, $\phi_0(r)$, and its derivative, additional physical quantities such as energy density ($\rho(r)$), radial pressure ($P_r(r)$), and tangential pressure ($P_\perp(r)$) are computed.
    6.  **Visualization:** Various plots are generated to visually represent the radial profiles of the scalar field, metric functions, mass, energy density, and pressures.
    7.  **Equation of State:** The relationships between pressure and energy density ($P_r$ vs. $\rho$ and $P_\perp$ vs. $\rho$) are plotted, offering insights into the effective equation of state within the boson star.
    8.  **ADM Mass Analysis:** The notebook demonstrates how to calculate the total Asymptotically Flat (ADM) mass of the boson star for a range of $\omega$ values and plots the characteristic $M_{ADM}$ vs. $\omega$ curve.
    9.  **Unit Conversion:** A dedicated section details the dimensional analysis and provides functions to scale the dimensionless simulation results into physical SI units, enabling real-world interpretation of boson star properties for a user-specified boson mass (e.g., in eV/c$^2$).

## Setup and Installation

### Prerequisites

Ensure you have the following software installed on your system:
*   **Python 3.x**
*   **Git**

### Clone the Repository

To obtain a local copy of this project, open your terminal or command prompt and execute the following commands:

```bash
git clone https://github.com/Hazem-74/Mini-Boson-Stars.git
cd Mini-Boson-Stars/OverView
```

### Install Dependencies

This project relies on several common Python libraries for numerical computation, data manipulation, and plotting. You can install them using `pip`:

```bash
pip install numpy scipy matplotlib seaborn pandas
```
Alternatively, if you prefer using Anaconda or Miniconda, you can create a dedicated environment:

```bash
conda create -n boson-stars python=3.9 numpy scipy matplotlib seaborn pandas
conda activate boson-stars
```

## Usage

### Running the Notebooks

1.  Navigate to the `OverView` directory where the notebook is located:
    ```bash
    cd Mini-Boson-Stars/OverView
    ```
2.  Launch Jupyter Notebook from your terminal:
    ```bash
    jupyter notebook
    ```
3.  Your default web browser will open, displaying the Jupyter interface. Click on `Mini-Boson Stars.ipynb` to open the notebook.
4.  You can then execute the cells sequentially to run the simulation, perform the analysis, and generate visualizations.

## Key Concepts and Theory

This project is grounded in the theoretical framework of General Relativity coupled with a scalar field theory. Key concepts explored include:
*   **Boson Stars:** Exotic compact objects hypothesized to form from a Bose-Einstein condensate of self-gravitating scalar particles.
*   **Einstein-Klein-Gordon (EKG) Equations:** The fundamental equations that govern the interaction between the gravitational field (described by Einstein's field equations) and a scalar field (described by the Klein-Gordon equation).
*   **Spherically Symmetric Spacetime:** The simplifying assumption of spherical symmetry allows for a one-dimensional radial analysis of the field equations.
*   **Boundary Value Problems (BVPs):** Numerical methods for solving differential equations where conditions are specified at multiple points (e.g., the center and asymptotic regions of the star).
*   **ADM Mass:** The Arnowitt-Deser-Misner (ADM) mass is a measure of the total gravitational mass of an asymptotically flat spacetime, representing the total energy of the star.
*   **Equation of State:** The relationship between pressure and energy density, which effectively describes the "matter" content (the scalar field) within the boson star.

## Results and Visualizations

The `Mini-Boson Stars.ipynb` notebook generates a series of informative plots, including:
*   Radial profiles of the scalar field amplitude ($\phi_0(r)$), gravitational potential ($\Phi(r)$), metric component ($\Lambda(r)$), and the enclosed mass ($m(r)$).
*   Radial distributions of energy density ($\rho(r)$) and the distinct radial ($P_r(r)$) and tangential ($P_\perp(r)$) pressures.
*   Implicit equation of state plots, illustrating radial and tangential pressures as functions of energy density.
*   A curve showing the total ADM mass of the boson star as a function of the scalar field's oscillation frequency ($\omega$).
*   Examples demonstrating the conversion of dimensionless simulation data into physical units (e.g., kilometers, solar masses) for a specified boson mass.

## License

This project is open-source and available under the MIT License.
