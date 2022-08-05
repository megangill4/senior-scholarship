#!/usr/bin/env python
# coding: utf-8

# # Molecular Orbitals
# 
# Now that we have displayed atomic orbitals, we can start to consider molecular orbitals. Much like atomic orbitals, molecular orbitals are described by wavefunctions, and are constructed by finding the linear combination of atomic orbitals from different atoms[<sup>1</sup>](#fn1). Similarly to atomic orbitals, we can find the probability density of an electron at each point in the molecular orbital by taking the square modulus of the wavefunction. Bonding orbitals are lower in energy than the atomic orbitals that form them, while antibonding orbitals are higher in energy. Filling an antibonding orbital effectively cancels out the bond formed by filling the bonding orbital. 
# 
# Conjugated compounds contain both sigma bonds (single bonds formed from the head on overlap of any orbital) and π bonds (double bonds formed from the side on overlap of p orbitals). For this project, I have chosen to focus on π bonds as it is π bonding that gives conjugated compounds their interesting properties. The simplest organic molecule containing a π bond is ethene:
# 
# ```{image} ../images/Picture1.png
# :alt: Ethene
# :width: 200px
# :align: center
# ```

# The π bond in ethene is constructed from the constructive side-on overlap of carbon 2p orbitals[<sup>2</sup>](#fn2).

# ```{image} ../images/pibonding.jpg
# :alt: Pi bonding
# :width: 600px
# :align: center
# ```

# An antibonding orbital is formed from the destructive overlap of orbitals.
# 
# ```{image} ../images/piantibonding.jpg
# :alt: Pi antibonding
# :width: 600px
# :align: center
# ```

# As each 2p orbital contributes one electron to the π system, there are two electrons in the π system, only the bonding MO is filled, and a π bond is formed.
# 
# ```{image} ../images/benzene.jpg
# :alt: Benzene
# :width: 150px
# :align: center
# ```
# 
# In benzene (above), six 2p orbitals contribute to the π system[<sup>3</sup>](#fn3). The number of orbitals is conserved, and we are left with six π molecular orbitals. Rather than a single bonding and a single antibonding π molecular orbital, we have six different molecular orbitals which each have different degrees of bonding and antibonding character, formed from different combinations of overlap between the carbon 2p orbitals, as demonstrated below. The three lowest energy molecular orbitals are occupied.
# 

# ```{image} ../images/benzenemo.png
# :alt: Benzene MOs
# :width: 500px
# :align: center
# ```
# 
# Using these diagrams, we can easily produce a qualitative picture of the π bonding in benzene. In order to produce a more quantative and accurate picture, we can employ Hückel Molecular Orbital theory.

# <span id="fn1"> 1: Atkins, P.W. and Julio De Paula (2010). Physical chemistry. 9th ed. New York: W.H. Freeman And Co. </span>
# 
# <span id="fn2"> 2: Atkins, P.W. and Julio De Paula (2010) </span>
# 
# <span id="fn3"> 3: Atkins, P.W. and Julio De Paula (2010) </span>

# In[ ]:




