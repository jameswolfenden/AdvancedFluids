# steel strength
σ = 400.0 # MPa

# localised pressure difference
p = 10000 # Pa

# pressure applicaton area
A = 0.1 # m^2

# plate size
L = 2.0 # m
W = 1.0 # m
t = 0.25e-3 # m

# stress from pressure
σ_p = 6 * p * A / (4 * pi * t^2) * log(1/A) # Pa
# convert to MPa
σ_p = σ_p / 1e6 # MPa