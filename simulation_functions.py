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


def build_simulator(NX, NY, dx, dy, h, active, kx, ky, mu, wells):
    """
    Constroi e resolve o sistema de equações de escoamento monofásico incompressível
    em duas dimensões utilizando diferenças finitas e lei de Darcy.

    Parâmetros:
    - NX, NY: dimensões da malha (número de células em X e Y)
    - dx, dy: listas com os tamanhos das células em cada direção (m)
    - h: espessura do reservatório (m)
    - active: matriz (NY x NX) com 1 para células ativas e 0 para inativas
    - kx, ky: matrizes de permeabilidade nas direções X e Y (mD)
    - mu: viscosidade do fluido (cP)
    - wells: lista de dicionários com informações dos poços (posição, tipo e valor)

    Retorno:
    - Matriz (NY x NX) com o campo de pressões (kPa) em cada célula ativa
    """
    N = NX * NY
    T = sp.lil_matrix((N, N))
    Q = np.zeros(N)

    def get_index(i, j):
        return i * NX + j

    for i in range(NY):
        for j in range(NX):
            if active[i, j] == 0:
                continue

            idx = get_index(i, j)

            # Para cada vizinho:
            # Cálculo da transmissibilidade com o vizinho à esquerda:
            # - Média harmônica da permeabilidade (kx_eff)
            # - Transmissibilidade Tx pela Lei de Darcy:
            #   Tx = (kx_eff * h) / (mu * distância média entre centros)
            # - Atualização da matriz T: termo diagonal e vizinho

            # Vizinho a esquerda
            if j > 0 and active[i, j-1] == 1:
                kx_eff = 2 * kx[i, j] * kx[i, j-1] / (kx[i, j] + kx[i, j-1])
                Tx = (kx_eff * h) / (mu * (dx[j] + dx[j-1]) / 2)
                T[idx, idx] += Tx
                T[idx, get_index(i, j-1)] -= Tx

            # Vizinho à direita
            if j < NX-1 and active[i, j+1] == 1:
                kx_eff = 2 * kx[i, j] * kx[i, j+1] / (kx[i, j] + kx[i, j+1])
                Tx = (kx_eff * h) / (mu * (dx[j] + dx[j+1]) / 2)
                T[idx, idx] += Tx
                T[idx, get_index(i, j+1)] -= Tx

            # Vizinho abaixo
            if i > 0 and active[i-1, j] == 1:
                ky_eff = 2 * ky[i, j] * ky[i-1, j] / (ky[i, j] + ky[i-1, j])
                Ty = (ky_eff * h) / (mu * (dy[i] + dy[i-1]) / 2)
                T[idx, idx] += Ty
                T[idx, get_index(i-1, j)] -= Ty

            # Vizinho acima
            if i < NY-1 and active[i+1, j] == 1:
                ky_eff = 2 * ky[i, j] * ky[i+1, j] / (ky[i, j] + ky[i+1, j])
                Ty = (ky_eff * h) / (mu * (dy[i] + dy[i+1]) / 2)
                T[idx, idx] += Ty
                T[idx, get_index(i+1, j)] -= Ty

        # Aplicação dos poços com controle por pressão ou vazão
    for w in wells:
        i, j = w["i"], w["j"]
        idx = get_index(i, j)

        if active[i, j] == 0:
            raise ValueError(f"Poço posicionado em célula inativa: ({i}, {j})")

        if w["controle"] == "PRESSAO":
            # Define pressão imposta: zera linha e coluna, coloca 1 na diagonal e valor em Q
            T[idx, :] = 0
            T[idx, idx] = 1
            Q[idx] = w["valor"]

        elif w["controle"] == "VAZAO":
            # Define vazão imposta: apenas adiciona valor em Q
            Q[idx] += w["valor"]

    T = T.tocsr()
    P = spla.spsolve(T, Q)
    return P.reshape((NY, NX))


def compute_well_flows(P, NX, NY, dx, dy, h, active, kx, ky, mu, wells):
    """
    Calcula a vazão dos poços com base nas pressões simuladas.

    Parâmetros:
    - P: matriz de pressão (NY x NX) [kPa]
    - NX, NY, dx, dy, h: dimensões da malha e propriedades geométricas
    - active: matriz (NY x NX) com 1 para células ativas e 0 para inativas
    - kx, ky: permeabilidades nas direções X e Y (mD)
    - mu: viscosidade do fluido (cP)
    - wells: lista de dicionários com os poços
      (chaves: "i", "j", "tipo", "controle", "valor")

    Retorno:
    - Lista de tuplas: (i, j, tipo, controle, pressão [kPa], vazão [m³/dia])
    """

    def get_index(i, j):
        return i * NX + j

    well_results = []

    for w in wells:
        i, j = w['i'], w['j']

        if active[i, j] == 0:
            raise ValueError(f"Poço posicionado em célula inativa: ({i}, {j})")

        controle = w['controle']
        tipo = w['tipo']
        
        if controle == 'PRESSAO':
            pressure = w['valor']
            flow = 0.0

            for ni, nj in [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]:
                if 0 <= ni < NY and 0 <= nj < NX and active[ni, nj]:
                    dp = pressure - P[ni, nj]
                    if ni == i:
                        # vizinho na mesma linha → direção X
                        k_eff = 2 * kx[i, j] * kx[ni, nj] / (kx[i, j] + kx[ni, nj])
                        dist = (dx[j] + dx[nj]) / 2
                    else:
                        # vizinho na mesma coluna → direção Y
                        k_eff = 2 * ky[i, j] * ky[ni, nj] / (ky[i, j] + ky[ni, nj])
                        dist = (dy[i] + dy[ni]) / 2

                    trans = (k_eff * h) / (mu * dist)
                    flow += trans * dp

            well_results.append((i, j, tipo, 'PRESSAO', pressure, flow))

        elif controle == 'VAZAO':
            flow = w['valor']
            pressure = P[i, j]
            well_results.append((i, j, tipo, 'VAZAO', pressure, flow))

        else:
            raise ValueError(f"Controle desconhecido: {controle}")

    # Impressão dos resultados
    print("Resultados dos Poços:")
    for i, j, tipo, controle, pres, q in well_results:
        print(f"Poço em ({i}, {j}) | Tipo: {tipo} | Controle: {controle} | Pressão [kPa]: {pres:.2f} | Vazão [m³/dia]: {q:.2f}")

    return well_results


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
