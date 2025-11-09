module PendulumSim

# Include submodules
include("SinglePendulum.jl")
include("DoublePendulum.jl")

# Re-export submodules so they're accessible as PendulumSim.SinglePendulum and PendulumSim.DoublePendulum
using .SinglePendulum
using .DoublePendulum

# Re-export main functions for convenience
export SinglePendulum, DoublePendulum

# Export commonly used functions from both modules
export run_single_pendulum_gui, run_double_pendulum_gui

# Make the run functions available directly
"""
    run_single_pendulum_gui()

Launch the single pendulum simulation GUI.
Alias for `SinglePendulum.run_gui()`.
"""
const run_single_pendulum_gui = SinglePendulum.run_gui

"""
    run_double_pendulum_gui()

Launch the double pendulum simulation GUI.
Alias for `DoublePendulum.run_double_pendulum_gui()`.
"""
const run_double_pendulum_gui = DoublePendulum.run_double_pendulum_gui

end # module PendulumSim
