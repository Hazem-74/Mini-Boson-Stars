# Mini-Boson Stars: Numerical Analysis

This repository presents a numerical analysis of mini-boson star configurations by solving the Einsteinâ€“Kleinâ€“Gordon system of equations in a spherically symmetric spacetime. The project focuses on understanding the fundamental properties of these exotic compact objects, including their mass, density, and pressure profiles, and how they relate to the underlying scalar field parameters.

The primary goal is to simulate stable boson star solutions and convert the dimensionless results into physically interpretable units, allowing for comparison with astrophysical phenomena.

## Table of Contents

-   [Mathematical & Physical Background](#mathematical--physical-background)
    -   [Metric and Scalar Field Ansatz](#metric-and-scalar-field-ansatz)
    -   [System of ODEs](#system-of-odes)
    -   [Boundary Conditions](#boundary-conditions)
-   [Project Structure](#project-structure)
-   [Key Features & Methodology](#key-features--methodology)
-   [Requirements](#requirements)
-   [Installation](#installation)
-   [Usage](#usage)
-   [Results](#results)
-   [Dimensional Analysis and Unit Rescaling](#dimensional-analysis-and-unit-rescaling)
-   [License](#license)

## Mathematical & Physical Background

Boson stars are theoretical compact objects formed by a self-gravitating scalar field, representing a potential dark matter candidate. This project models non-rotating, spherically symmetric boson stars in general relativity.

### Metric and Scalar Field Ansatz

The spacetime is described by a static, spherically symmetric metric in Schwarzschild-like coordinates:

$$
ds^2 = -e^{2\Phi(r)} dt^2 + e^{2\Lambda(r)} dr^2 + r^2 d\Omega^2
$$

The scalar field is assumed to be complex, massive, and oscillating with a frequency $\omega$:

$$
\phi(r,t) = \phi_0(r) \, e^{i\omega t}
$$
where $\phi_0(r)$ is the real amplitude, and $\omega$ is the oscillation frequency.

### System of ODEs

By substituting the metric and scalar field ansatz into the Einstein field equations and the Klein-Gordon equation, and applying the conditions for dimensionless geometrized units ($G = c = \hbar = m_b = 1$, where $m_b$ is the boson mass), we obtain a coupled system of first-order ordinary differential equations (ODEs) for $\Phi(r)$ (gravitational potential), $m(r)$ (enclosed mass function), $\phi_0(r)$ (scalar field amplitude), and $\phi_0'(r)$ (derivative of scalar field amplitude). The $\Lambda(r)$ metric function is derived from the enclosed mass $m(r)$.

The `system` function within the notebook defines these ODEs:

-   **Mass Function**: $m'(r)$
-   **Gravitational Potential**: $\Phi'(r)$
-   **Scalar Field Equation**: $\phi_0'' + \left( \frac{2}{r} + \Lambda' - \Phi' \right)\phi_0' + \left( \omega^2 e^{2(\Lambda - \Phi)} - e^{2\Lambda} \right)\phi_0 = 0$
-   **Metric Closure**: $\Lambda(r) = -\frac{1}{2} \ln\left(1 - \frac{2m(r)}{r} \right)$

### Boundary Conditions

To obtain physically regular and asymptotically flat solutions, specific boundary conditions are imposed:

**At** $ r = 0 $:
-   $\phi_0(0) = \phi_c$ (a finite central scalar field value)
-   $\phi_0'(0) = 0$ (no cusp at the center)
-   $\Phi(0) = 0$ (regular gravitational potential at the center)
-   $m(0) = 0$ (no point mass at the origin)

**As** $ r \to \infty $:
-   $\lim_{r \to \infty} \phi_0(r) = 0$ (scalar field vanishes at infinity)
-   $\lim_{r \to \infty} \Phi(r) = 0$ (asymptotically flat spacetime)

## Project Structure

This project consists of a single Jupyter Notebook:

-   `Mini-Boson Stars.ipynb`: This notebook contains all the code for defining the system of ODEs, solving them numerically, analyzing the results, and converting them to physical units.

## Key Features & Methodology

The `Mini-Boson Stars.ipynb` notebook implements the following:

-   **Numerical Solution**: The `scipy.integrate.solve_bvp` function is used to solve the two-point boundary value problem defined by the system of ODEs and their boundary conditions.
-   **Derived Quantities**: From the primary solutions, essential physical quantities like energy density ($\rho(r)$), radial pressure ($P_r(r)$), and tangential pressure ($P_\perp(r)$) are calculated based on the energy-momentum tensor of the scalar field.
-   **Profile Visualization**: Plots illustrating the radial profiles of the scalar field amplitude, metric functions, enclosed mass, energy density, and pressures are generated.
-   **Equation of State Analysis**: The relationship between pressure and energy density ($P$ vs. $\rho$) is plotted, revealing the implicit equation of state for boson stars.
-   **Parameter Study**: The notebook includes an analysis of how the total ADM (Arnowitt-Deser-Misner) mass of the boson star varies with the scalar field's oscillation frequency ($\omega$).
-   **Dimensional Analysis**: A detailed section explains the scaling relations required to convert the dimensionless simulation results (obtained in units where $G = c = \hbar = m_b = 1$) into physical SI and astrophysical units (e.g., kilometers, solar masses, Hertz, J/m$^3$).
-   **Physical Unit Conversion**: Python functions (`convert_to_physical_units` and `convert_to_physical_units_df`) are provided to perform these conversions, allowing users to interpret the results for specific boson masses (e.g., in eV/c$^2$).

## Requirements

To run this project, you need:

-   Python 3.x
-   The following Python libraries:
    -   `numpy`
    -   `seaborn`
    -   `pandas`
    -   `matplotlib`
    -   `scipy`
    -   `jupyter` (for running the notebook)

## Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/Hazem-74/Mini-Boson-Stars.git
    ```

2.  **Navigate to the project directory:**
    ```bash
    cd Mini-Boson-Stars/Complete Project
    ```

3.  **Install the required Python packages:**
    It is recommended to use a virtual environment.
    ```bash
    pip install -r requirements.txt
    ```
    If a `requirements.txt` file is not present, you can create one with the following content:
    ```
    numpy
    seaborn
    pandas
    matplotlib
    scipy
    jupyter
    ```
    Then run `pip install -r requirements.txt`.

## Usage

1.  **Launch Jupyter Notebook:**
    ```bash
    jupyter notebook
    ```

2.  **Open the notebook:**
    In your web browser, navigate to and open `Mini-Boson Stars.ipynb`.

3.  **Run the cells:**
    Execute the notebook cells sequentially to reproduce the analysis. You can modify parameters like the central scalar field value (`phi_center`), oscillation frequency (`omega`), or the boson mass for physical unit conversion (`m_b_eV`) to explore different boson star configurations.

## Results

Upon successful execution, the notebook will generate several plots, including:

-   `basic_profiles.png`: Radial profiles of the scalar field, metric functions, and mass.
-   `pressure_vs_density.png`: The implicit equation of state for the boson star.
-   `rho_pressures_vs_radius.png`: Radial profiles of energy density and pressures.
-   `ADM_mass_vs_omega.png`: A plot showing the ADM mass of the boson star as a function of the scalar field oscillation frequency $\omega$.

These plots illustrate the characteristics of mini-boson stars and their response to varying parameters.

## Dimensional Analysis and Unit Rescaling

A key aspect of this project is the conversion of dimensionless simulation results, obtained in a natural unit system where $G = c = \hbar = m_b = 1$, into conventional physical units. This is achieved through specific scaling relations detailed within the notebook, based on the boson's mass $m_b$.

The scaling relations are provided for:
-   **Length**: Scaled by the reduced Compton wavelength $\lambda_c = \frac{\hbar}{m_b c}$.
-   **Mass**: Scaled by $M_0 = \frac{\hbar}{G m_b}$, convertible to solar masses ($M_\odot$).
-   **Frequency**: Scaled by $f_{\text{scale}} = \frac{m_b c^2}{\hbar}$, convertible to Hertz.
-   **Energy Density/Pressure**: Scaled by $\rho_0 = \frac{m_b^4 c^6}{\hbar^3}$, convertible to J/m$^3$.
-   **Scalar Field Amplitude**: Scaled to ensure consistency with energy density units.

These conversions allow the theoretical models to be compared with observational data and explore different dark matter boson candidates by varying $m_b$.

## License

This project is open-source. Please refer to the repository's license for detailed terms.

