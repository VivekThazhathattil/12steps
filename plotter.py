import csv
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

CSV_FILE_NAME = 'step1_out.csv'

def animate(i):
  ax.clear()
  point_data = data[i]
  line, = ax.plot(point_data, color='black')
  ax.set_xlim(0, len(point_data) - 1) 
  ax.set_ylim(min_y - min_y/25, max_y + max_y/25)
  return line,

data = []
max_y = -100000.0
min_y = 100000.0

with open(CSV_FILE_NAME, 'r') as csvfile:
  reader = csv.reader(csvfile)
  for row in reader:
    data.append([float(val) for val in row if val != ''])
    min_y = min(min(data[-1]), min_y)
    max_y = max(max(data[-1]), max_y)

fig, ax = plt.subplots(figsize=(13,3))
animation = FuncAnimation(fig, animate, frames=len(data), interval=3, blit=False)
writer = animation.save('step1_anim.gif')
plt.show()

