# PendulumSim

Interactive physics simulations of single and double pendulum systems with real-time visualization using GLMakie.

## Installation

```julia
# Clone the repository
git clone https://github.com/yourusername/PendulumSim.jl.git
cd PendulumSim.jl

# Activate and instantiate the project
julia --project=.
```

```julia
# In the Julia REPL
using Pkg
Pkg.instantiate()
```

## Usage

```julia
using PendulumSim

# Launch the single pendulum simulation
run_single_pendulum_gui()

# Launch the double pendulum simulation
run_double_pendulum_gui()

# Or access submodules directly
PendulumSim.SinglePendulum.run_gui()
PendulumSim.DoublePendulum.run_double_pendulum_gui()
```

## Features

### Single Pendulum
- Interactive 2D visualization of physical pendulum (rod + bob) using RigidBodyDynamics
- Real-time phase space plots (θ vs ω or θ vs ω/ω₀)
- Comparison with Simple Harmonic Oscillator (SHO) for oscillating motion
- Accurate separatrix showing boundary between oscillation and rotation
- 4-period time series display with exact elliptic integral period calculation
- Adjustable parameters: damping coefficient, length, bob mass, rod density
- Energy tracking (KE/PE proportions) and automatic trajectory classification
- Interactive controls: play/pause, reset, randomize IC, click to set IC in phase space
- Normalized velocity toggle for circular phase space orbits in small-angle regime

### Double Pendulum
- Real-time animation of chaotic double pendulum dynamics
- Multiple visualization modes:
  - 2D trajectory tracking with color-coded trails
  - 2D phase space portraits
  - 3D phase space with selectable projections (drop θ₁, ω₁, θ₂, or ω₂)
- Interactive parameter controls:
  - Mass sliders for both bobs (m₁, m₂)
  - Length sliders for both rods (L₁, L₂)
  - Gravity slider (0-20 m/s²)
  - Damping coefficient (0-0.5)
  - Time speed multiplier (0.1-3.0x)
- Energy monitoring:
  - Real-time kinetic and potential energy display
  - Total energy tracking with drift percentage
  - Visual energy bar representation
- FPS counter for performance monitoring
- Projection dropdown for 3D phase space customization

## Project Structure

```
PendulumSim.jl/
├── src/
│   ├── PendulumSim.jl          # Main module (parent)
│   ├── SinglePendulum.jl        # Single pendulum submodule
│   └── DoublePendulum.jl        # Double pendulum submodule
├── test/                        # Test suite
├── Project.toml                 # Package dependencies
├── README.md                    # This file
└── TODO.md                      # Development roadmap
```

## Module Hierarchy

The package uses a hierarchical module structure:

- **PendulumSim** (parent module)
  - **SinglePendulum** (submodule)
    - `run_gui()` - Launch single pendulum GUI
  - **DoublePendulum** (submodule)
    - `run_double_pendulum_gui()` - Launch double pendulum GUI

Both simulation functions are re-exported at the top level for convenience.

## Dependencies

- **DifferentialEquations.jl** - ODE solving with Tsit5 integrator
- **GLMakie.jl** - GPU-accelerated interactive visualization
- **RigidBodyDynamics.jl** - 3D kinematics for single pendulum
- **StaticArrays.jl** - Efficient array operations
- **Colors.jl** - Color gradients for trail visualization

## Physics Background

### Single Pendulum
The single pendulum is governed by the nonlinear ODE:
```
θ'' + (c/m)θ' + (g/L)sin(θ) = 0
```
where θ is the angle from vertical, c is damping, m is mass, L is length, and g is gravity.

### Double Pendulum
The double pendulum exhibits chaotic dynamics governed by coupled nonlinear equations:
```
(m₁ + m₂)L₁θ₁'' + m₂L₂θ₂''cos(θ₁-θ₂) + m₂L₂θ₂'²sin(θ₁-θ₂) - (m₁+m₂)g·sin(θ₁) - cθ₁' = 0
L₂θ₂'' + L₁θ₁''cos(θ₁-θ₂) - L₁θ₁'²sin(θ₁-θ₂) - g·sin(θ₂) - cθ₂' = 0
```

The simulation uses symplectic integration (Tsit5) to preserve energy conservation in the Hamiltonian system.

## Development

### AI Agent Integration via MCPRepl

This project includes AI agent integration for development:

**For AI agents**: See [AGENTS.md](AGENTS.md) for detailed guidelines.

### Quick Start for Development

```bash
cd PendulumSim.jl
./repl  # MCP server starts automatically
```

**VS Code**: Open project, start Julia REPL
**Claude Desktop**: Config in `.mcp.json`
**Gemini**: Configured in `~/.gemini/settings.json`

## Security

**Mode**: `lax` | **Port**: `3076` | **Auth**: None (localhost only)

To change security mode:

```julia
using MCPRepl
MCPRepl.setup()
```


## Troubleshooting

### General Issues

**Simulation won't start?**
- Ensure all dependencies are installed: `Pkg.instantiate()`
- Check that GLMakie can create windows on your system

**Performance issues?**
- Reduce trail lengths in the source code
- Lower the time speed multiplier
- Close other GPU-intensive applications

### Development/MCP Issues

**Port in use?** Override with: `JULIA_MCP_PORT=3001 julia --project=.`

**Auth fails?** Check API key: `cat .env`

**Server won't start?** Restart Julia or check port: `lsof -i :3076`

## Screenshots

The package includes interactive visualizations with real-time parameter controls and multiple display modes for analyzing pendulum dynamics.

## Contributing

Contributions are welcome! Areas for improvement:

- [ ] Implement symplectic integrators for better energy conservation
- [ ] Add potential/kinetic energy trend indicators
- [ ] Fix cyclic fill/drain functionality
- [ ] Add more visualization options (Poincaré sections, Lyapunov exponents)
- [ ] Implement triple pendulum simulation
- [ ] Add export functionality for animations and data

See [TODO.md](TODO.md) for the current development roadmap.

## License

MIT License - See LICENSE file for details.

## Author

Kahli Burke <kahli@kahliburke.com>
