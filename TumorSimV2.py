import random
import json
from turtle import position
import matplotlib.pyplot as plt # might be using, not sure
import numpy as np
import argparse
from matplotlib.colors import ListedColormap

#This could be stored in a JSON
MUTATION_TYPES = [
    # Affects proliferation_chance
    "TP53",     # ↑ proliferation_chance (loss of checkpoint)
    "KRAS",     # ↑ proliferation_chance (growth signaling)
    "EGFR",     # ↑ proliferation_chance (receptor signaling)
    "RB1",      # ↑ proliferation_chance (removes cycle block)
    "BRAF",     # ↑ proliferation_chance (MAPK pathway)
    "ALK",      # ↑ proliferation_chance (fusion kinase)
    "HER2",     # ↑ proliferation_chance (signal amplification)
    "APC",      # ↑ proliferation_chance (WNT deregulation)
    "CDKN2A",   # ↑ proliferation_chance (cycle deregulation)

    # Affects mutation_rate
    "BRCA1",    # ↑ mutation_rate (DNA repair failure)
    "BRCA2",    # ↑ mutation_rate (similar to BRCA1)
    "PIK3CA",   # ↑ mutation_rate (growth/survival pathway)
    "NOTCH1",   # ↑ mutation_rate (assumed oncogenic)

    # Affects crowding sensitivity (pressure effect)
    "PTEN",     # ↓ pressure effect (loss of inhibition)
    "SMAD4",    # ↓ pressure effect (less growth suppression)

    # Affects aggressiveness
    "MYC",      # ↑ aggressiveness (proliferation driver)
    "IDH1",     # ↓ aggressiveness slightly (metabolic shift)
]

MUTATION_EFFECTS = {
    "TP53":   {"proliferation_chance": +0.04},
    "KRAS":   {"proliferation_chance": +0.03},
    "EGFR":   {"proliferation_chance": +0.03},
    "PIK3CA": {"mutation_rate": +0.01, "resistance": +0.1},
    "PTEN":   {"pressure_sensitivity": -0.2},
    "RB1":    {"proliferation_chance": +0.02},
    "BRAF":   {"proliferation_chance": +0.03},
    "ALK":    {"proliferation_chance": +0.03},
    "MYC":    {"aggressiveness": +0.2},
    "BRCA1":  {"mutation_rate": +0.01},
    "BRCA2":  {"mutation_rate": +0.01},
    "HER2":   {"proliferation_chance": +0.03},
    "SMAD4":  {"pressure_sensitivity": -0.2},
    "APC":    {"proliferation_chance": +0.02},
    "IDH1":   {"aggressiveness": -0.1},
    "NOTCH1": {"mutation_rate": +0.01},
    "CDKN2A": {"proliferation_chance": +0.03},
}

# Define color mapping and priority for subtypes

SUBTYPE_COLORS = {
    "Proliferative": 1,
    "Mutator": 2,
    "Resistant": 3,
    "Aggressive": 4,
    "Insensitive": 5,
    "Unclassified": 6,
}
SUBTYPE_PRIORITY = [
    "Proliferative", "Mutator", "Resistant", "Aggressive", "Insensitive", "Unclassified"
]

custom_cmap = ListedColormap([
    "#ffffff",  # 0: Normal = white
    "#e41a1c",  # 1: Proliferative = red
    "#377eb8",  # 2: Mutator = blue
    "#4daf4a",  # 3: Resistant = green
    "#ff7f00",  # 4: Aggressive = orange
    "#984ea3",  # 5: Insensitive = purple
    "#999999",  # 6: Unclassified = gray
])

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
        if 0 <= x < self.width and 0 <= y < self.height:
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
        grid_data = np.zeros((self.height,self.width), dtype=int)
        
        for y in range(self.height):
            for x in range(self.width):
                cell = self.grid[y][x]
                if cell.cell_type == 'cancer':
                    for subtype in SUBTYPE_PRIORITY:
                        if subtype in cell.subtype:
                            grid_data[y][x] = SUBTYPE_COLORS[subtype]
                            break # only take highest priority here
                else:
                    grid_data[y][x] = 0 #normal cell


        plt.figure(figsize=(6, 6))
        plt.imshow(grid_data, cmap=custom_cmap, interpolation='nearest')
        colorbar_labels = ['Normal'] + SUBTYPE_PRIORITY
        cbar = plt.colorbar(ticks=range(len(colorbar_labels)))
        cbar.ax.set_yticklabels(colorbar_labels)
        plt.title('Tumor Growth by Subtype')
        plt.show()

class Tumor():

    def __init__(self,environment):
        self.iteration_count = 0
        self.cells = [] # list of all cancerous cells
        self.environment = environment#the environment in which the tumor grows

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
    
    # returns number of cancer neighbors at x,y position
    def cancer_neighbor_count(self,pos):
        x,y = pos
        neighbors = self.environment.get_neighbors(x,y)
        total = 0
        for c in neighbors:
            if c.cell_type == 'cancer':
                total += 1
        return total



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

        #prioritize gaps between cells instead of being fully random: prevents tumor for being spaced out too much
        if random.random() < 0.8:
            available_positions.sort(key=self.cancer_neighbor_count,reverse=True)
            new_position = available_positions[0]
        else:
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
        return f"Cell(type={self.cell_type}, pos={self.position}, age={self.age})"

    

class Cancer_Cell(Cell):
    def __init__(self,position,mutation_rate =  0.01, proliferation_chance = 0.3,aggressiveness = 1.2):
        super().__init__(position)
        self.mutations = set()
        self.mutation_count = 0
        self.mutation_rate = mutation_rate
        self.proliferation_chance = proliferation_chance# chance to divide during iteration
        self.cell_type = 'cancer'
        self.aggressiveness = aggressiveness #scaling function to be implemented
        #add additional self.aggressiveness (this can be used to scale regular cells)
        self.resistance = 0 #TODO: implement treatment which can be negated through increased resistance
        self.subtype = [] # dynamically add subtypes based on whether certain metrics are reached
        self.pressure_sensitivity = 1.0 # how sensitive the cell is to crowding effects (when other cells are around it, it divides more than usual cells would)
        

    #different grow function because cancer cells accumulate more mutations (tp53)
    def grow(self):
        self.age += 1
        self.mutate()
        return self

    #determining what subtype the cell is based on whether or not it meets criteria
    def determine_subtypes(self):
        self.subtype = [] #inefficient, but prevents duplicates

        # If proliferation chance exceeds ~0.35, tag as highly proliferative
        # Reference: Ki‑67 thresholds (~20–30%) used to distinguish aggressive/high‑growth tumors
        if self.proliferation_chance > 0.35:
            self.subtype.append("Proliferative")

        # Mutator phenotype often emerges early and leads to rapid mutation accrual
        # Ref: Vogelstein et al.’s model & mutator phenotype hypothesis
        if self.mutation_rate > 0.02:
            self.subtype.append("Mutator")

        # Resistance above ~0.1 implies significant survival advantage under treatment
        # (Estimates broadly align with mutation-driven resistance scales)
        if self.resistance > 0.1:
            self.subtype.append("Resistant")

        # Aggressiveness index significantly above baseline (>1.5) flags aggressive clone
        if self.aggressiveness > 1.5:
            self.subtype.append("Aggressive")

        # Pressure sensitivity below ~0.7 suggests the cell tolerates crowding better
        # (modeling assumption tied to reduced contact inhibition)
        if self.pressure_sensitivity < 0.7:
            self.subtype.append("Insensitive")

        if not self.subtype:
            self.subtype.append("Unclassified")

    def mutate(self):
        if random.random() < self.mutation_rate * self.aggressiveness:
            possible_mutations = list(set(MUTATION_TYPES) - self.mutations)
            if possible_mutations:
                new_mutation = random.choice(possible_mutations)
                self.mutations.add(new_mutation)
                self.mutation_count = len(self.mutations)
                
                #apply mutation effects to cell
                effects = MUTATION_EFFECTS.get(new_mutation, {})
                for attr, value in effects.items():
                    if attr == 'proliferation_chance':
                        self.proliferation_chance += value
                    elif attr == 'mutation_rate':
                        self.mutation_rate += value
                    elif attr == 'resistance':
                        self.resistance += value
                    elif attr == 'aggressiveness':
                        self.aggressiveness += value
                    elif attr == 'pressure_sensitivity':
                        self.pressure_sensitivity += value
                        self.pressure_sensitivity = max(0.0, self.pressure_sensitivity) #pressure sensitivity cant go below 0 that would mess up the program

                self.determine_subtypes() # recalculate subtypes
                return True
        return False

    #checks if cell will divide, returns bool, takes into account pressure, which is calculated by the tumor class
    def should_divide(self,pressure = 0.0): 
        effective_chance = self.proliferation_chance * (1 - pressure * self.pressure_sensitivity)
        return random.random() < effective_chance

    #this function is useful because new cancer cells should be the same as the parent cells
    def clone(self, new_position):
        clone = Cancer_Cell(
            position=new_position,
            mutation_rate=self.mutation_rate,
            proliferation_chance=self.proliferation_chance,
            aggressiveness=self.aggressiveness
        )
        clone.mutations = set(self.mutations)
        clone.mutation_count = self.mutation_count
        clone.resistance = self.resistance
        clone.pressure_sensitivity = self.pressure_sensitivity
        clone.subtype = list(self.subtype)
        return clone


    def __repr__(self):
        return (f"Cell(type={self.cell_type}, pos={self.position}, age={self.age}, "
                f'mutations={self.mutations}, mutationCount = {self.mutation_count}, prolif={self.proliferation_chance}, aggressiveness={self.aggressiveness}')

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
        #store_step()

    tumor.environment.visualize()


cd 