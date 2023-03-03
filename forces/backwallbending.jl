# fridge Height
H = 2.0 # m
# fridge width
W = 1.0 # m
# steel thickness
t = 0.25e-3 # m

# pressure
pATM = 1.01513
p = (pATM - 1.01325) * 10e5 # Pa
println("Pressure: $(p) Pa")

# pressure applicaton width
pwidth = 0.05 # m

# start position of load
x0 = 0.0 # m
# end position of load
x1 = 0.03 # m
# load
w = p * pwidth
# Force
F = w * (x1 - x0) # N

# Young's modulus of steel
E = 200.0e9 # Pa

# Moment of inertia
I = t^3 * W / 12 # m^4

d = (H - x0 - x1)

# Reaction at top
RT = F / (4 * H^2) * (12 * d^2 - 8 * d^3 / H + 2 * x1 * (x1 - x0)^2 / H - (x1 - x0)^3 / H - (x1 - x0)^2) # N
# Reaction at bottom
RB = F - RT # N

# Moment at top
MT = -F / (24 * H) * (24 * d^3 / H - 6 * x1 * (x1 - x0)^2 / H + 3 * (x1 - x0)^3 / H + 4 * (x1 - x0)^2 - 24 * d^2) # Nm
# Moment at bottom
MB = F / (24 * H) * (24 * d^3 / H - 6 * x1 * (x1 - x0)^2 / H + 3 * (x1 - x0)^3 / H + 2 * (x1 - x0)^2 - 48 * d^2) # Nm

# mid point of load
x = (x0 + x1) / 2 # m

# Deflection at x
y = 1 / (6 * E * I) * (RT * x^3 - 3 * MT * x^2 - w / 4 * (x - x0)^4)

println("Deflection at x = $(x) m: $(y) m")

# Deflection at something
h = 0.85 # m
y = 1 / (6 * E * I) * (RB * (H-h)^3 - 3 * MB * (H-h)^2 - w / 4 * (H - h)^4)

println("Deflection at h = $(h) m: $(y) m")

# Assume point load at mid point of load - pinned at top and bottom
xmax = H-sqrt((H^2 - (x)^2) / 3)
y = F * (x) * (H^2 - (x)^2)^(3 / 2) / (9 * sqrt(3) * H * E * I)

println("Deflection at xmax = $(xmax) m: $(y) m assuming point load at x = $(x) m")

# Assume point load at mid point of load - fixed at top and bottom
M1 = -F*x/H^2*(H-x)^2
R1 = F/H^3*(H-x)^2*(H+2*x)
h=0.85
y = M1*h^2/(2*E*I)+R1*h^3/(6*E*I)-F/(6*E*I)*(h-3*x)^3

println("Deflection at h = $(h) m: $(y) m assuming point load at x = $(x) m")