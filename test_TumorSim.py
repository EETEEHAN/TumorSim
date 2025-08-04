import pytest
from TumorSimV7 import Environment, Tumor, Cancer_Cell, Cell

@pytest.fixture
def small_env_and_tumor():
    env = Environment(width=10, height=10)
    env.initialize_grid()
    tumor = Tumor(env)
    return env, tumor

def test_environment_initialization(small_env_and_tumor):
    env, _ = small_env_and_tumor
    assert len(env.grid) == 10
    assert len(env.grid[0]) == 10
    assert env.grid[5][5].cell_type == 'normal'

def test_seed_initial_cancer(small_env_and_tumor):
    env, tumor = small_env_and_tumor
    tumor.seed_initial_cancer()
    mid_x, mid_y = env.width // 2, env.height // 2
    assert env.grid[mid_y][mid_x].cell_type == 'cancer'
    assert len(tumor.cells) == 1

def test_cancer_cell_growth_and_mutation():
    cell = Cancer_Cell(position=(5, 5))
    initial_mutation_count = cell.mutation_count
    mutated = False
    for _ in range(100):
        cell.grow()
        if cell.mutation_count > initial_mutation_count:
            mutated = True
            break
    assert mutated, "Cell did not mutate after 100 growth cycles"

def test_division_creates_new_cell(small_env_and_tumor):
    env, tumor = small_env_and_tumor
    tumor.seed_initial_cancer()
    cancer_cell = tumor.cells[0]
    cancer_cell.proliferation_chance = 1.0

    # Clear neighbors to ensure there's room
    x, y = cancer_cell.position
    for dx, dy in [(-1, 0), (1, 0)]:
        nx, ny = x + dx, y + dy
        if 0 <= nx < env.width and 0 <= ny < env.height:
            env.grid[ny][nx] = Cell(position=(nx, ny))  # Normal cell

    success = tumor.divide_cell(cancer_cell)
    assert success
    assert len(tumor.cells) > 1

def test_subtype_assignment():
    cell = Cancer_Cell(position=(0, 0))
    cell.proliferation_chance = 0.6
    cell.mutation_rate = 0.03
    cell.resistance = 0.2
    cell.aggressiveness = 1.6
    cell.pressure_sensitivity = 0.2
    cell.determine_subtypes()
    assert "Proliferative" in cell.subtype
    assert "Mutator" in cell.subtype
    assert "Resistant" in cell.subtype
    assert "Aggressive" in cell.subtype
    assert "Insensitive" in cell.subtype

def test_get_local_pressure_low_density(small_env_and_tumor):
    env, tumor = small_env_and_tumor
    tumor.seed_initial_cancer()
    pressure = tumor.get_local_pressure(tumor.cells[0])
    assert 0 <= pressure <= 1

def test_cancer_neighbor_count_accuracy(small_env_and_tumor):
    env, tumor = small_env_and_tumor
    tumor.seed_initial_cancer()
    mid_x, mid_y = env.width // 2, env.height // 2
    count = tumor.cancer_neighbor_count((mid_x, mid_y))
    assert isinstance(count, int)

def test_clone_copies_all_attributes():
    original = Cancer_Cell(position=(1, 1))
    original.mutations = {"TP53", "KRAS"}
    original.mutation_rate = 0.02
    original.proliferation_chance = 0.5
    original.resistance = 0.3
    original.aggressiveness = 1.8
    original.pressure_sensitivity = 0.5
    original.determine_subtypes()

    clone = original.clone((2, 2))
    assert clone.position == (2, 2)
    assert clone.mutations == original.mutations
    assert clone.mutation_rate == original.mutation_rate
    assert clone.proliferation_chance == original.proliferation_chance
    assert clone.resistance == original.resistance
    assert clone.aggressiveness == original.aggressiveness
    assert clone.pressure_sensitivity == original.pressure_sensitivity
    assert clone.subtype == original.subtype

def test_should_divide_true_under_low_pressure():
    cell = Cancer_Cell(position=(0, 0), proliferation_chance=1.0)
    cell.pressure_sensitivity = 0
    assert cell.should_divide(pressure=1.0) is True

def test_cell_grow_increases_age():
    cell = Cell(position=(0, 0))
    assert cell.age == 0
    cell.grow()
    assert cell.age == 1

def test_get_neighbors_returns_correct_count():
    env = Environment(5, 5)
    env.initialize_grid()
    neighbors = env.get_neighbors(2, 2)
    assert len(neighbors) == 8

    edge_neighbors = env.get_neighbors(0, 0)
    assert len(edge_neighbors) == 3  # top-left corner

def test_is_valid_position():
    env = Environment(5, 5)
    assert env.is_valid_position(2, 2) is True
    assert env.is_valid_position(-1, 0) is False
    assert env.is_valid_position(0, 5) is False
