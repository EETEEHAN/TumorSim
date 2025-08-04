import random
import json
from turtle import position
import matplotlib.pyplot as plt # might be using, not sure
import numpy as np
import argparse

class Environment():
    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.grid = [] #create a grid of unmutated cells
        #self.nutrition potentially different environemtns make a difference

    #fill grid with cells (with one cancer cell in the middle?)
    def initialize_grid(self):
        for y in range(self.height):
            row = [] # creating a list representing each row
            for x in range(self.width):
                row.append(Cell((x,y))) # filling the row with cells with correct coordinates
            self.grid.append(row)
    #cell = self.grid[y][x]  grid can be accessed like this

    #returns bool depending on if pos is within grid
    def is_valid_position(self,x,y):
        if x < self.width and y < self.height:
            return True
        return False

    #returns bool if position is taken by cancer cell (this will be good for growth)
    def is_occupied(self,x,y):
        if self.grid[y][x].cell_type == 'cancer':
            return True
        return False

    #place cell at given position (for tumor growth) if not taken by tumor cell? or depends on implementation
    def place_cell(self,cell,x,y):
        if not self.is_occupied(x,y):
            self.grid[y][x] = cell
            return True
        return False

    #returns list of neighboring cells (8 directions unless on border), this is useful for determining if crowding is occuring and if the tumor can grow further
    def get_neighbors(self,x,y):
        directions = [ # rotate around the cell starting from above and going clockwise
            (-1,0), #north
            (-1,+1),
            (0,+1), #east
            (+1,+1),
            (+1,0), #south
            (+1,-1),
            (0,-1), #west
            (-1,-1)
            ]
        neighbors = []

        for direction_row, direction_col in directions:
            coordinate_row = y + direction_row
            coordinate_col = x + direction_col

            if 0 <= coordinate_col < self.width and 0 <= coordinate_row < self.height:
                neighbors.append(self.grid[coordinate_row][coordinate_col])

        return neighbors

    #plot tumor using ascii/plt
    def visualize(self):
        grid_data = np.zeros((self.height,self.width))

        for y in range(self.height):
            for x in range(self.width):
                cell = self.grid[y][x]
                if cell.cell_type == 'cancer':
                    grid_data[y][x] = 1
                else:
                    grid_data[y][x] = 0

        plt.figure(figsize = (6,6))
        plt.imshow(grid_data, cmap = 'hot', interpolation = 'nearest')
        plt.title('Tumor Growth')
        plt.colorbar(label = 'Cell Type(0 = normal, 1 = cancer)')
        plt.show()

class Tumor():

    def __init__(self,environment):
        self.iteration_count = 0
        self.cells = [] # list of all cancerous cells
        self.environment = environment #the environment in which the tumor grows

    #this is optional for now as th env class is doing this
    def seed_initial_cancer(self, cancer_cell=None):
        """Seeds the initial cancer cell in the middle of the environment grid."""
        middle_x = self.environment.width // 2
        middle_y = self.environment.height // 2
        position = (middle_x, middle_y)

        if cancer_cell is None:
            cancer_cell = Cancer_Cell(position=position)

        self.environment.place_cell(cancer_cell, position[0], position[1])
        self.cells.append(cancer_cell)
        return self

    #run one iteration, growing cells and checking for division
    def step(self):
        self.iteration_count += 1

        # making a copy of the list so that we can shuffle it so that division is fully random, this way creates a whole list not just a reference
        shuffled_cells = self.cells[:]
        random.shuffle(shuffled_cells)

        for c in shuffled_cells:
            c.grow()

            if isinstance(c, Cancer_Cell):
                pressure = self.get_local_pressure(c)

                if c.should_divide(pressure):
                    self.divide_cell(c)
        #DEBUG
        print(f"Step {self.iteration_count}: {len(self.cells)} cancer cells")

        self.store_step()


    #this is used for crowding (higher pressure means more cancer cells around, less likely to divide)
    def get_local_pressure(self,cell):
        neighbors = self.environment.get_neighbors(cell.position[0],cell.position[1])
        total_neighbors = len(neighbors)
        cancer_neighbors = 0.0
        for c in neighbors:
            if c.cell_type == 'cancer':
                cancer_neighbors += 1.0

        if total_neighbors == 0:
            return 1.0

        return cancer_neighbors/total_neighbors



    #attempt to divide to neighboring spot, returns false if cant divide, otherwise chooses random available neighbor and creates identical cell there\
    #only cancer cells should call this function
    def divide_cell(self,cell):
        x,y = cell.position
        neighbors = self.environment.get_neighbors(cell.position[0],cell.position[1])
        pressure = self.get_local_pressure(cell)
        if pressure == 1.0:
            return False
        
        # Find empty neighbor positions
        available_positions = []
        direction_map = [
            (-1, 0), (-1, 1), (0, 1), (1, 1),
            (1, 0), (1, -1), (0, -1), (-1, -1)
        ]
        for direction_x, direction_y in direction_map:
            neighbor_x = direction_x + x
            neighbor_y = direction_y + y
                
            if 0 <= neighbor_x < self.environment.width and 0 <= neighbor_y < self.environment.height:
                if self.environment.grid[neighbor_y][neighbor_x].cell_type != 'cancer':
                    available_positions.append((neighbor_x,neighbor_y))

        if not available_positions:
            return False

        new_position = random.choice(available_positions)
        new_cell = cell.clone(new_position) #cancer cells only
        self.environment.grid[new_position[1]][new_position[0]] = new_cell
        self.cells.append(new_cell)

        return True


    #store current data on iteration as a dictionary (current idea) to be added to a df which can be converted to json for time series data
    def store_step(self):
        pass

    #next steps if time allows: treatment (bottleneck effect)



class Cell():
    def __init__(self,position):
        self.position = position
        self.age = 0
        self.cell_type = 'normal'
    
    #increases age, each increase in age is equivalent to a self-renewal
    def grow(self):
        self.age += 1
        return self


    #returns position on grid
    def get_position(self):
        return self.position

    def __repr__(self):
        return f"Cell(type={self.cell_type}, pos={self.position}, age={self.age}"

    

class Cancer_Cell(Cell):
    def __init__(self,position,mutation_rate =  0.01, proliferation_chance = 0.3,aggressiveness = 1.2):
        super().__init__(position)
        self.mutations = 0
        self.mutation_rate = mutation_rate
        self.proliferation_chance = proliferation_chance# chance to divide during iteration
        self.cell_type = 'cancer'
        self.aggressiveness = aggressiveness #scaling function to be implemented
        #add additional self.aggressiveness (this can be used to scale regular cells)

    #different grow function because cancer cells accumulate more mutations (tp53)
    def grow(self):
        self.age += 1
        self.mutate()
        return self

    def mutate(self):
        if random.random() < self.mutation_rate * self.aggressiveness:
            self.mutations += 1
            return True
        return False

    #checks if cell will divide, returns bool, takes into account pressure, which is calculated by the tumor class
    def should_divide(self,pressure = 0.0): 
        effective_chance = self.proliferation_chance * (1 - pressure)
        return random.random() < effective_chance

    #this function is useful because new cancer cells should be the same as the parent cells
    def clone(self,new_position):
        return Cancer_Cell(position = new_position,
                           mutation_rate=self.mutation_rate,
                           proliferation_chance=self.proliferation_chance,
                           aggressiveness=self.aggressiveness)

    def __repr__(self):
        return (f"Cell(type={self.cell_type}, pos={self.position}, age={self.age}, "
                f'mutations={self.mutations}, prolif={self.proliferation_chance}, aggressiveness={self.aggressiveness}')

    #next steps if time allows: resistance to treatment, metasticize

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = 'Tumor Simulator Debug',
        description = 'This is a command-line implementation of TumorSim and is for development purposes',
        epilog = 'A full web-services version will be available for public use soon')

    parser.add_argument('-W', '--width', type = int, default = 100, help = 'Width of cell matrix')
    parser.add_argument('-H', '--height', type = int, default = 100, help = 'Height of cell matrix')
    parser.add_argument('-S', '--steps', type = int, default = 50, help = 'Number of simulation iterations to run (50 default)')
    # Mutation rate override
    parser.add_argument('--mutation_rate', type=float, default=0.01, help='Mutation rate of cancer cells (default: 0.01)')
    # Proliferation chance override
    parser.add_argument('--proliferation', type=float, default=0.3, help='Base proliferation chance (default: 0.3)')
    # Aggressiveness override
    parser.add_argument('--aggressiveness', type=float, default=1.2, help='Aggressiveness multiplier for mutation/division (default: 1.2)')
    args = parser.parse_args()

    middle_position = (args.width // 2, args.height // 2)

    env = Environment(args.width,args.height)
    tumor = Tumor(env)
    tumor.environment.initialize_grid()
    tumor.seed_initial_cancer(Cancer_Cell(position = middle_position, mutation_rate=args.mutation_rate,proliferation_chance=args.proliferation,aggressiveness=args.aggressiveness))
    for i in range(args.steps):
        tumor.step()

    tumor.environment.visualize()


