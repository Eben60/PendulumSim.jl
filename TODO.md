# TODOS

- ✅ DONE: Run button now works, data appears, Random IC button works
- Not done: Window in floating mode at upper left using GLFW calls after display
-  Make sure phase plot is wrapping at +/- pi correctly.
Add a feature to click in the phase space and have that set up the associated ICs for theta and omega.
# Todo List

- [] Set window position to upper left - NOT DONE
  - Set window to floating mode and position at upper left of screen
- [x] Show two periods in theta vs time plot - COMPLETED
  - Adjust time axis to show exactly two periods of oscillation
NOT DONE - simulation should keep running beyond two periods, but the plot should only show the last two periods.
- [x] Add SHO comparison to theta plot - COMPLETED
  - Add SHO comparison with half-opaque line and fill between to show difference
- NOT DONE - SHO comparison is not what I want. I want the periods of the actual motion and the period of a hyptothetical SHO to match, so that you can see how the actual motion deviates from SHO. This requires calculating the exact period based on initial angle.
- create mechanism for agent to shut down and restart GUI independently. Discuss ideas with me before imlementing
- Speed up sim slightly so it runs a little faster in real time.
- Set drag value to zero.
- Add display of of potential energy and kinetic energy value at the current time, and whether the trend is for that to go up or down at the current time.
- Improve integration accuracy or use symplectic integrator.
- Add circular buffer of screenshot images, one every second, stored in a directory with timestamps, so that an AI agent can look back at recent history. Old screenshots get removed after 50 accumulate
- [x] Add omega plot with SHO comparison - COMPLETED
  - Omega plot now shows actual omega vs time with SHO comparison line (orange dashed) and fill between (yellow transparent)
