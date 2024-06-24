import matplotlib.pyplot as plt
import numpy as np

world_size = 5

cmax_values = []
for rank in range(world_size):
    filename = f"cmax_values_rank_{rank}.txt"
    with open(filename, "r") as f:
      cmaxes = [int(line.strip()) for line in f]
    plt.plot(cmaxes)
plt.show()
