# Single Polyanion Example
I've created a very simple example where we simulate a single polyanion chain in a box. To perform the simulation, there are two steps:

1. Initialize the structure: run `python polyanion_init.py` to generate a file called `input.data`. This file contains your initial structure.
2. Run the simulation: run `lmp -in polyanion.in`. This will run a short simulation of this single chain. You can visualize the simulation using `drop.lammpstrj` in Ovito.