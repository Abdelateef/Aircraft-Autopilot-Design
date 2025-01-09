# Aircraft-Autopilot-Design

A comprehensive project on designing and simulating an autopilot system for the Boeing 747 (B747-FC3). This project features flight dynamics, control systems, numerical solutions, and performance validation using MATLAB and Simulink.

## Project Structure

- **`src/`**: Contains all the MATLAB scripts and Simulink models.
  - `RK4_Solver.m`: MATLAB implementation of the Runge-Kutta 4th order method for solving ODEs.
  - `EOM_Simulator.slx`: Simulink model for the aircraft equations of motion.
  - `Control_Design/`: Subfolder with control design scripts for longitudinal and lateral controllers.

- **`docs/`**: Documentation detailing the methodology, results, and validation.
  - `Project_Report.pdf`: Comprehensive report summarizing the project.
  - `Control_Diagrams/`: Subfolder containing block diagrams and flowcharts.

- **`tests/`**: Benchmark data and test scripts.
  - `Validation_Benchmark.m`: MATLAB script for comparing simulation results with benchmark data.

- **`plots/`**: Generated plots and visualizations of system responses and controller performance.
  - `Step_Response_Analysis/`: Subfolder containing step response plots.
  - `Root_Locus_Bode/`: Subfolder with root locus and Bode plots for control systems.

- **`README.md`**: Overview of the project, including structure, setup, and usage.

## Getting Started

1. Clone the repository:
   ```bash
   git clone https://github.com/username/Aircraft-Autopilot-Design.git
   ```
2. Navigate to the project directory:
   ```bash
   cd Aircraft-Autopilot-Design
   ```
3. Open MATLAB and add the project directory to your MATLAB path.
4. Explore the `src/` folder to run simulations or view control designs.

## Key Features

- Aircraft Type: Boeing 747 (B747-FC3)
- Numerical Solutions:
  - ODE solving using RK4 method.
  - Linearization and model approximation.
- Control Systems:
  - Longitudinal and lateral control designs.
  - Velocity, pitch, yaw, roll, and altitude controllers.
- Validation:
  - MATLAB and Simulink-based simulation comparisons.
  - Benchmark testing.

---

For detailed documentation and setup instructions, refer to the `docs/` folder.
