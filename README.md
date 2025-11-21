# Mini-Boson Stars: Numerical Analysis

This repository presents a numerical study of **mini-boson stars**—hypothetical compact astrophysical objects composed of self-gravitating scalar bosons. By solving the **Einstein–Klein–Gordon (EKG) equations** under spherical symmetry, the project models the internal structure, spacetime geometry, and physical properties (mass, energy density, pressure) of these exotic objects.

A key feature is the **conversion of dimensionless simulation results into physical units** (e.g., kilometers, solar masses, eV/c²), enabling real-world astrophysical interpretation for various boson masses—potentially relevant to dark matter research.

## Structure

This repository contains two complementary implementations:

- **`/OverView`**: A self-contained Jupyter notebook with full setup instructions, theory, and visualization. Ideal for first-time users and educational purposes.  
  → See [`OverView/README.md`](OverView/README.md) for details.

- **`/Complete Project`**: A streamlined version focused on core computation and unit rescaling, with modular code structure. Suitable for parameter studies and extension.  
  → See [`Complete Project/README.md`](Complete%20Project/README.md) for details.

Both folders include the main notebook:  
**`Mini-Boson Stars.ipynb`**

## Requirements

- Python 3.x
- Libraries: `numpy`, `scipy`, `matplotlib`, `seaborn`, `pandas`, `jupyter`

## Quick Start

```bash
git clone https://github.com/Hazem-74/Mini-Boson-Stars.git
cd Mini-Boson-Stars
pip install numpy scipy matplotlib seaborn pandas jupyter
jupyter notebook
