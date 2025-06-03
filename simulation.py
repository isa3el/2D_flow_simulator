import numpy as np
import os
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
from simulation_functions import read_input_file, build_simulator, compute_well_flows, export_well_results
from plot_functions import plot_pressure_field, plot_permeability_fields, plot_pressure_map_with_values

# Para rodar o simulador, defina o diretório que contém os arquivos .txt de entrada.
# Descomente os diretório que será utilizado.
#
# Diretórios disponíveis:
# - Grid_60_60: reservatório com malha 60x60 e 9 poços (5 produtores controlados por BHP e 4 injetores controlados por vazão).
# - Grid_15_15: reservatório com malha 15x15 e 3 poços (2 produtores controlados por BHP e 1 injetor controlado por vazão).
# - Grid_5_5: reservatório com malha 5x5 e 2 poços (1 produtor controlado por BHP e 1 injetor controlado por vazão).

# Os diretórios disponíveis possuem "60_60", "15_15" e "5_5" no nome apenas para fins de identificação.
# O simulador é genérico e capaz de rodar reservatórios de qualquer tamanho, desde que os arquivos de entrada estejam no formato esperado.


# Descomente abaixo o diretório que será utilizado:
data_dir = "Grid_60_60"
#data_dir = "Grid_15_15"
#data_dir = "Grid_5_5"

#Definindo os arquivos a serem lidos a partir da 
input_path = os.path.join(data_dir, "input.txt")
grid_path = os.path.join(data_dir, "grid.txt")
permx_path = os.path.join(data_dir, "perm_x.txt")
permy_path = os.path.join(data_dir, "perm_y.txt")

NX, NY, dx, dy, h, active, kx, ky, mu, rho, wells = read_input_file(input_path, grid_path, permx_path, permy_path)

P = build_simulator(NX, NY, dx, dy, h, active, kx, ky, mu, wells)
print(P)

wells_prod = compute_well_flows(P, NX, NY, dx, dy, h, active, kx, ky, mu, wells)

plot_pressure_field(P, wells, data_dir)
plot_permeability_fields(kx, ky, wells, data_dir)
plot_pressure_map_with_values(P, wells, data_dir)

export_well_results(wells_prod, data_dir, NX, NY, dx, dy, h, mu, rho, P, kx, ky)
