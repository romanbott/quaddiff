"""Default constants."""

# ---- Trajectory Solver constants
# Integration constants
MAX_TIME = 100000
VELOCITY_SCALE = 0.02
MAX_STEP = 0.5

# Event constants
CLOSE_2START = 0.05
CLOSE_2POLE = 0.1
CLOSE_2ZERO = 0.05
LIM = 30

# ---- Trajectory constants
# Simplification
DISTANCE_2LINE = 0.03
MIN_DISTANCE = 0.01
#Refining
MAX_DISTANCE = 0.001
# Convergence
DISTANCE_2LIMIT = CLOSE_2POLE
