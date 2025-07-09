# This script creates a doodle to explain how dmax is computed visually
import matplotlib.pyplot as plt

# Data for doodle
x = [-7, -5, -1, 0, 1, 5, 7]
y = [5, 5, 1, 1, 1, 5, 5]

# Create doodle
fig, ax = plt.subplots(figsize=(12,8))
plt.plot(x, y, color = 'blue', linewidth = 5)

plt.axhline(y=0, color='black', linestyle='-', linewidth = 2)
plt.xlim(-7, 7)
plt.xticks([-5, -1, 1, 5],
           ["-$\\Delta x_m / \\alpha$", "$-w$", "$w$", "$\\Delta x_m / \\alpha$"],
           fontsize = 20)

plt.axvline(x=0, color='black', linestyle='-', linewidth = 2)
plt.axvline(x=0.5, color='red', linestyle='--', linewidth = 3)
plt.axvline(x=-0.5, color='red', linestyle='--', linewidth = 3)
plt.ylim(-0.5, 6)
plt.yticks([1, 5],
           ["$\\alpha w$", "$\\Delta x_m$"],
           fontsize = 20)


# Display doodle
plt.tight_layout()
plt.show()