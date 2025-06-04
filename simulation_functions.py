import numpy as np
import os
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def read_input_file(filename, grid_file, permx_file, permy_file):
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    def extract_block(tag):
        try:
            start = lines.index(f'[{tag}]') + 1
            end = next((i for i in range(start, len(lines)) if lines[i].startswith('[')), len(lines))
            return lines[start:end]
        except ValueError:
            raise ValueError(f"Tag [{tag}] não encontrada no arquivo.")

    NX, NY = map(int, extract_block("NX NY")[0].split())
    dx_block = extract_block("DX")[0].split('*')
    if len(dx_block) == 2:
        dx = [float(dx_block[1])] * int(dx_block[0])
    else:
        dx = list(map(float, extract_block("DX")[0].split()))
    dy_block = extract_block("DY")[0].split('*')
    if len(dy_block) == 2:
        dy = [float(dy_block[1])] * int(dy_block[0])
    else:
        dy = list(map(float, extract_block("DY")[0].split()))
    h = float(extract_block("H")[0])
    active = np.loadtxt(grid_file).reshape((NY, NX)).astype(int)
    kx = np.loadtxt(permx_file).reshape((NY, NX))
    ky = np.loadtxt(permy_file).reshape((NY, NX))
    mu, rho = map(float, extract_block("PROPS")[0].split())

    wells_block = extract_block("WELLS")
    n_wells = int(wells_block[0])
    wells = []
    for line in wells_block[1:n_wells + 1]:
        tipo_raw, i, j, controle, valor = line.split()
        tipo = "INJETOR" if tipo_raw.upper() == "INJ" else "PRODUTOR"
        well = {
            "i": int(i),
            "j": int(j),
            "tipo": tipo,
            "controle": controle.upper(),
            "valor": float(valor)
        }
        wells.append(well)

    print('Gridblocks dimensão X (NX) [número de células]:', NX)
    print('Gridblocks dimensão Y (NY) [número de células]:', NY)
    print('Tamanho das células na direção X (dx) [m]:', dx)
    print('Tamanho das células na direção Y (dy) [m]:', dy)
    print('Espessura do reservatório (h) [m]:', h)
    print('Viscosidade dinâmica (mu) [cP]:', mu)
    print('Densidade (rho) [kg/m³]:', rho)
    print('Poços:')
    for w in wells:
        print('  - Coordenadas (i, j):', (w['i'], w['j']), '| Tipo:', w['tipo'],'| Controle:',w['controle'],  '| Valor:', w['valor'], '[kPa ou m³/dia]')
    if NX*NY < 400:
        print('Mapa de células ativas/inativas:', active)
        print('Permeabilidade na direção X (kx) [mD]:', kx)
        print('Permeabilidade na direção Y (ky) [mD]:', ky)
    
    return NX, NY, dx, dy, h, active, kx, ky, mu, rho, wells


import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

def build_simulator(NX, NY, dx, dy, h, active, kx, ky, mu, wells, rw=0.1, s=0.0):
    """
    Constroi e resolve o sistema de equações de escoamento monofásico incompressível
    em duas dimensões.
    
    Parâmetros:
    - NX, NY: dimensões da malha em X e Y
    - dx, dy: tamanhos das células nas direções X e Y
    - h: espessura do reservatório
    - active: matriz (NY x NX) indicando se uma célula está ativa (1) ou não (0)
    - kx, ky: permeabilidades nas direções X e Y (mD)
    - mu: viscosidade do fluido (cP convertido para Pa.s)
    - wells: lista de dicionários com dados dos poços (posição, tipo, controle, valor)
    - rw: raio do poço (default: 0.1 m)
    - s: skin factor (default: 0.0)

    Retorno:
    - Campo de pressão (matriz NY x NX)
    """

    N = NX * NY                 # número total de células
    T = sp.lil_matrix((N, N))   # matriz de transmissibilidades (forma esparsa)
    Q = np.zeros(N)             # vetor de fontes (vazões)

    def get_index(i, j):
        # Converte coordenadas (i,j) em índice 1D
        return i * NX + j

    for i in range(NY):
        for j in range(NX):
            if active[i, j] == 0:
                continue    # pula células inativas

            idx = get_index(i, j)

            # Transmissibilidade nas 4 direções
            
            # Transmissibilidade para a esquerda
            if j > 0 and active[i, j - 1]:
                kx_eff = 2 * kx[i, j] * kx[i, j - 1] / (kx[i, j] + kx[i, j - 1])
                Tx = (kx_eff * h) / (mu * (dx[j] + dx[j - 1]) / 2)
                T[idx, idx] += Tx
                T[idx, get_index(i, j - 1)] -= Tx

            # Transmissibilidade para a direita
            if j < NX - 1 and active[i, j + 1]:
                kx_eff = 2 * kx[i, j] * kx[i, j + 1] / (kx[i, j] + kx[i, j + 1])
                Tx = (kx_eff * h) / (mu * (dx[j] + dx[j + 1]) / 2)
                T[idx, idx] += Tx
                T[idx, get_index(i, j + 1)] -= Tx

            # Transmissibilidade para cima
            if i > 0 and active[i - 1, j]:
                ky_eff = 2 * ky[i, j] * ky[i - 1, j] / (ky[i, j] + ky[i - 1, j])
                Ty = (ky_eff * h) / (mu * (dy[i] + dy[i - 1]) / 2)
                T[idx, idx] += Ty
                T[idx, get_index(i - 1, j)] -= Ty

            # Transmissibilidade para baixo
            if i < NY - 1 and active[i + 1, j]:
                ky_eff = 2 * ky[i, j] * ky[i + 1, j] / (ky[i, j] + ky[i + 1, j])
                Ty = (ky_eff * h) / (mu * (dy[i] + dy[i + 1]) / 2)
                T[idx, idx] += Ty
                T[idx, get_index(i + 1, j)] -= Ty

    for w in wells:
        i, j = w["i"], w["j"]  # posição do poço
        idx = get_index(i, j)

        if active[i, j] == 0:
            raise ValueError(f"Poço em célula inativa ({i}, {j})")

        ctrl = w["controle"]    # "PRESSAO" ou "VAZAO"
        valor = w["valor"]      # valor associado ao controle

        if ctrl == "VAZAO":
             # Controle por vazão: fonte ou sumidouro
            Q[idx] += valor

        elif ctrl == "PRESSAO":
            # Controle por pressão: WI
            kx_ij = kx[i, j]
            ky_ij = ky[i, j]

            dx_ = dx[j] if isinstance(dx, (list, np.ndarray)) else dx
            dy_ = dy[i] if isinstance(dy, (list, np.ndarray)) else dy

            if kx_ij <= 0 or ky_ij <= 0:    
                continue  

            # Cálculo do raio equivalente (req) para meio anisotrópico
            term1 = np.sqrt(ky_ij / kx_ij) * dx_**2
            term2 = np.sqrt(kx_ij / ky_ij) * dy_**2
            numerator = np.sqrt(term1 + term2)
            denominator = (ky_ij / kx_ij)**0.25 + (kx_ij / ky_ij)**0.25
            req = 0.28 * numerator / denominator

            # Termo logarítmico no denominador do WI
            log_term = np.log(req / rw) + s
            if abs(log_term) < 1e-6:
                log_term = 1e-6

            WI = (2 * np.pi * np.sqrt(kx_ij * ky_ij) * h) / log_term

            T_poço = WI / mu
            T[idx, idx] += T_poço
            Q[idx] += T_poço * valor  

    # Resolve o sistema linear T * P = Q
    T = T.tocsr()
    P = spla.spsolve(T, Q)
    return P.reshape((NY, NX))


def compute_well_flows(P, kx, ky, h, mu, wells, dx, dy, rw=0.1, s=0.0):
    """
    Calcula a vazão dos poços usando WI com anisotropia e adiciona a cada poço a chave 'q_wi'.

    Modifica os dicionários da lista wells in-place.

    Retorna:
    - Lista de tuplas: (i, j, tipo, controle, pressão [kPa], vazão [m³/dia]) para possível exportação
    """
    results = []

    for w in wells:
        i, j = w['i'], w['j']
        ctrl = w['controle']
        tipo = w.get('tipo', '')
        val = w['valor']

        if ctrl == 'PRESSAO':
            p_well = val
            p_res = P[i, j]

            kx_ij = kx[i, j]
            ky_ij = ky[i, j]

            dx_ = dx[j] if isinstance(dx, (list, np.ndarray)) else dx
            dy_ = dy[i] if isinstance(dy, (list, np.ndarray)) else dy

            # Raio efetivo com anisotropia (fórmula)
            term1 = np.sqrt(ky_ij / kx_ij) * dx_**2
            term2 = np.sqrt(kx_ij / ky_ij) * dy_**2
            numerator = np.sqrt(term1 + term2)
            denominator = (ky_ij / kx_ij)**0.25 + (kx_ij / ky_ij)**0.25
            req = 0.28 * numerator / denominator

            WI = (2 * np.pi * np.sqrt(kx_ij * ky_ij) * h) / (np.log(req / rw) + s)

            q = WI * (p_well - p_res) / mu
            w['q_wi'] = q  # adiciona ao dicionário
            w['p_wi'] = p_well

            results.append((i, j, tipo, ctrl, p_well, q))

        elif ctrl == 'VAZAO':
            q = val
            p_res = P[i, j]
            w['q_wi'] = q
            w['p_wi'] = p_res

            results.append((i, j, tipo, ctrl, p_res, q))

    return results



def export_well_results(well_results, data_dir, NX, NY, dx, dy, h, mu, rho, P, kx, ky):
    """
    Exporta os resultados dos poços e informações do reservatório.
    Se NX * NY < 400, exporta também os mapas de pressão e permeabilidade.

    Parâmetros:
    - well_results: lista de tuplas (i, j, tipo, controle, pressão, vazão)
    - output_path: caminho do arquivo .txt
    - NX, NY: dimensões da malha
    - dx, dy: listas de tamanho das células [m]
    - h: espessura do reservatório [m]
    - mu: viscosidade [cP]
    - rho: densidade [kg/m³]
    - P, kx, ky: matrizes 2D (NY x NX), opcionais
    """
    
    output_path = os.path.join(data_dir, "results.txt")
    
    with open(output_path, 'w') as f:
        f.write("=== INFORMAÇÕES DO RESERVATÓRIO ===\n")
        f.write(f"Dimensões da malha: NX = {NX}, NY = {NY}\n")
        f.write(f"Tamanho das células em X (dx) [m]: {dx}\n")
        f.write(f"Tamanho das células em Y (dy) [m]: {dy}\n")
        f.write(f"Espessura (h) [m]: {h}\n")
        f.write(f"Viscosidade do fluido (mu) [cP]: {mu}\n")
        f.write(f"Densidade do fluido (rho) [kg/m³]: {rho}\n\n")

        f.write("=== RESULTADOS DOS POÇOS ===\n")
        f.write(f"{'i':>3} {'j':>3} {'tipo':<10} {'controle':<10} {'pressao_kPa':>12} {'vazao_m3_dia':>14}\n")
        for i, j, tipo, controle, pres, q in well_results:
            f.write(f"{i:>3} {j:>3} {tipo:<10} {controle:<10} {pres:12.2f} {q:14.2f}\n")

        f.write("\n=== MAPA DE PRESSÃO [kPa] ===\n")
        for row in P:
            f.write(" ".join(f"{val:7.2f}" for val in row) + "\n")

        f.write("\n=== PERMEABILIDADE NA DIREÇÃO X (kx) [mD] ===\n")
        for row in kx:
            f.write(" ".join(f"{val:7.2f}" for val in row) + "\n")

        f.write("\n=== PERMEABILIDADE NA DIREÇÃO Y (ky) [mD] ===\n")
        for row in ky:
            f.write(" ".join(f"{val:7.2f}" for val in row) + "\n")

    print(f"Arquivo salvo em: {output_path}")


def export_well_results_comparativo(wells_prod_darcy, wells_prod_wi, data_dir, NX, NY, dx, dy, h, mu, rho, P, kx, ky):
    """
    Exporta os resultados dos poços comparando o método de Darcy e o WI.
    Inclui também informações do reservatório e, se possível, mapas.

    Parâmetros:
    - wells_prod_darcy: lista (i, j, tipo, controle, pressão, vazao_darcy)
    - wells_prod_wi: lista (i, j, tipo, controle, pressão, vazao_wi)
    - data_dir: diretório onde salvar
    - NX, NY, dx, dy, h, mu, rho: propriedades do reservatório
    - P, kx, ky: mapas 2D
    """
    
    output_path = os.path.join(data_dir, "results.txt")

    with open(output_path, 'w') as f:
        f.write("=== INFORMAÇÕES DO RESERVATÓRIO ===\n")
        f.write(f"Dimensões da malha: NX = {NX}, NY = {NY}\n")
        f.write(f"Tamanho das células em X (dx) [m]: {dx}\n")
        f.write(f"Tamanho das células em Y (dy) [m]: {dy}\n")
        f.write(f"Espessura (h) [m]: {h}\n")
        f.write(f"Viscosidade do fluido (mu) [cP]: {mu}\n")
        f.write(f"Densidade do fluido (rho) [kg/m³]: {rho}\n\n")

        f.write("=== COMPARAÇÃO DE RESULTADOS DOS POÇOS ===\n")
        f.write(f"{'i':>3} {'j':>3} {'tipo':<10} {'controle':<10} {'pressao_kPa':>12} {'q_Darcy':>12} {'q_WI':>12}\n")
        
        for wd, ww in zip(wells_prod_darcy, wells_prod_wi):
            i1, j1, tipo1, ctrl1, pres1, q_darcy = wd
            i2, j2, tipo2, ctrl2, pres2, q_wi = ww

            # Verificação básica (assume mesma ordem)
            assert (i1, j1) == (i2, j2), f"Poços em ordem diferente: ({i1},{j1}) ≠ ({i2},{j2})"

            f.write(f"{i1:>3} {j1:>3} {tipo1:<10} {ctrl1:<10} {pres1:12.2f} {q_darcy:12.2f} {q_wi:12.2f}\n")

        
        f.write("\n=== MAPA DE PRESSÃO [kPa] ===\n")
        for row in P:
            f.write(" ".join(f"{val:7.2f}" for val in row) + "\n")

        f.write("\n=== PERMEAfBILIDADE NA DIREÇÃO X (kx) [mD] ===\n")
        for row in kx:
            f.write(" ".join(f"{val:7.2f}" for val in row) + "\n")

        f.write("\n=== PERMEABILIDADE NA DIREÇÃO Y (ky) [mD] ===\n")
        for row in ky:
            f.write(" ".join(f"{val:7.2f}" for val in row) + "\n")

    print(f"Arquivo salvo em: {output_path}")
