print("\n\nFridge heat loss calculator\n")
# Dimensions
fridge_height = 2.0 # m
fridge_depth = 1.0 # m
fridge_width = 1.0 # m, not stated but unit width for simplicity

# Surface areas
A_door = fridge_height * fridge_width # m^2
A_side = fridge_height * fridge_depth # m^2
A_bottom = fridge_width * fridge_depth # m^2
A_top = fridge_width * fridge_depth # m^2
A_back = fridge_height * fridge_width # m^2

# Heat transfer coefficients
q_door = 12.92 # W/m^2
q_side = 13.34 # W/m^2
q_bottom = 14.33 # W/m^2
q_top = 10.10 # W/m^2
q_back = 15.18 # W/m^2

# Temperature difference
measured_temp_delta = 31.2-4.7
new_temp_delta = 20.0-5.0

# Adjust heat transfer coefficients
q_door = q_door * (new_temp_delta/measured_temp_delta)
q_side = q_side * (new_temp_delta/measured_temp_delta)
q_bottom = q_bottom * (new_temp_delta/measured_temp_delta)
q_top = q_top * (new_temp_delta/measured_temp_delta)
q_back = q_back * (new_temp_delta/measured_temp_delta)

# Heat loss
Q_door = q_door * A_door # W
Q_side = q_side * A_side # W
Q_bottom = q_bottom * A_bottom # W
Q_top = q_top * A_top # W
Q_back = q_back * A_back # W

# Total heat loss
Q_total = Q_door + Q_side * 2 + Q_bottom + Q_top + Q_back # W
print("Total heat loss: $(Q_total) W")