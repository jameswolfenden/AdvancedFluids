# thickness of gold leaf
t = 0.1e-6 # m

# Ultimate tensile strength of gold
UTS = 120.0 # MPa

# Height of gold leaf
H = 0.1 # m

# Cross sectional area of gold leaf
A = t * H # m^2

# Force required to fail gold leaf
F = UTS * A # N

println("Force required to fail gold leaf: $(F) N")