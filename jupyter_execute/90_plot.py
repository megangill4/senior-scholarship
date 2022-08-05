#!/usr/bin/env python
# coding: utf-8

# ## 90% Plot

# The second method I used to display the atomic orbitals was a 90% plot. The output I generated represents the region of 90% probability of finding an electron or pair of electrons. In order to create these plots, having determined my probability distribution function, I created a 50x50x50 3d grid of coordinates, and found the probability at each, and appended the coordinates and their corresponding probabilities to a list.
# 
# ```{figure} images/grid.png
# ---
# scale: 100%
# align: left
# ---
# 3D grid
# ```
# 
# The next step was to create a while loop. The loop works through the list of coordinates and probabilities, adding the probabilities together and appending the coordinates to corresponding lists. The loop stops running once the total probability reaches 0.9, or 90%. The regions with displayed points (from the lists of coordinates) therefore represent the position of the electron or pair of electrons 90% of the time. The colour of the points corresponds to the phase of the orbital in that region (red is positive, blue is negative). Below is my code for generating a 90% plot of atomic orbitals, as well as my generated plots:

# In[ ]:




