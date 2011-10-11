# An attempt to recreate the algorithm for the periodic boundary condition as implemented in
# ice-Ih-gen.f.
import random

def generateGrid(rows, cols):
	"""Attempting to generate a list of oxygen atoms representing "hexagonal" ice."""
	# Characterize each atom by its position on the x-axis.
	grid = []
	for i in range(rows):
		#Each row generate cols number of atoms with random distance to the next
		#one.
		row = [random.randint(-50, 150) for x in range(cols)]
		row = sorted(row)
		grid.append(row) #make sure to sort the list before appending it.
		print grid
	return grid

def applyBoundaryCondition(grid):
	"""Attempt to apply the boundary conditions to a grid of atoms"""
	for i in range(len(grid)): #rows
		print "Row: "
		print grid[i]
		for j in range(len(grid[i])): #cols

			#Sets the next molecule -- if it's the last one in the row, it wraps.
			if j+1 > len(grid[i]) - 1:
				nextM = grid[i][0]
			else:
				nextM = grid[i][j+1]

			dist = nextM - grid[i][j]
			if dist < (-1 * size[1]/2):
				grid[i][j] = grid[i][j] + size[1]
			if dist > (size[1]/2):
				grid[i][j] = grid[i][j] - size[1]

		print "New row: "
		print grid[i]		

def applyNewBoundaryCondition(grid):
	"""Attempt to apply the new boundary conditions to a grid of atoms"""
	for i in range(len(grid)): #rows
		print "Row: "
		print grid[i]
		for j in range(len(grid[i])): #cols

			#Sets the next molecule -- if it's the last one in the row, it wraps.
			if j+1 > len(grid[i]) - 1:
				nextM = grid[i+1][0]
			else:
				nextM = grid[i][j+1]

			dist = nextM - grid[i][j]
			if dist < (-1 * size[1]/2):
				grid[i][j] = grid[i][j] + size[1]
			if dist > (size[1]/2):
				grid[i][j] = grid[i][j] - size[1]

		print "New row: "
		print grid[i]		
#Arbitrary x and y size values
size = [100,100]
grid = generateGrid(10,10)
print grid
applyBoundaryCondition(grid)
print grid
