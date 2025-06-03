import os 
import matplotlib.pyplot as plt



def plot_pressure_field(P, wells, data_dir): 
    '''
    Plota o campo de pressão em uma malha 2D com a localização dos poços.
    Salva o plot como imagem PNG no diretório especificado.

    Parâmetros:
    - P: matriz 2D com os valores de pressão (em kPa) por célula.
    - wells: lista de dicionários com as chaves:
        'i' (linha), 'j' (coluna), 'tipo' ('INJETOR' ou 'PRODUTOR'), 'controle' (tipo de controle)
    - data_dir: caminho do diretório onde a imagem será salva.
    '''
    plt.figure(figsize=(8, 6))
    im = plt.imshow(P, cmap="jet", origin="lower")
    plt.colorbar(im, label="Pressão [kPa]")

    for w in wells:
        i, j = w['i'], w['j']
        if w['tipo'] == 'PRODUTOR':
            plt.plot(j, i, marker='v', color='magenta', markersize=10, label='Produtor')
        elif w['tipo'] == 'INJETOR':
            plt.plot(j, i, marker='o', color='violet', markersize=10, label='Injetor')

    # Remover duplicatas da legenda
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='upper right')

    plt.title("Campo de Pressão")
    plt.xlabel("X (células)")
    plt.ylabel("Y (células)")
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(os.path.join(data_dir, "mapa_pressao.png"), dpi=600)
    plt.show()


def plot_permeability_fields(kx, ky, wells, data_dir):
    """
    Plota mapas lado a lado das permeabilidades nas direções X e Y.
    Salva o plot como imagem PNG, destacando os poços com base no tipo e controle.

    Parâmetros:
    - kx: matriz 2D com permeabilidade na direção X [mD]
    - ky: matriz 2D com permeabilidade na direção Y [mD]
    - wells: lista de dicionários com chaves: 'i', 'j', 'tipo', 'controle'
    - data_dir: caminho do diretório onde a imagem será salva.
    """

    fig, axs = plt.subplots(1, 2, figsize=(14, 6))

    def get_marker_style(well):
        if well['tipo'] == 'INJETOR':
            marker = 'o'
            color = 'violet' if well['controle'] == 'VAZAO' else 'blue'
        elif well['tipo'] == 'PRODUTOR':
            marker = 'v'
            color = 'magenta' if well['controle'] == 'PRESSAO' else 'red'
        return marker, color

    # --- Mapa kx ---
    im1 = axs[0].imshow(kx, cmap='viridis', origin='lower')
    axs[0].set_title("Permeabilidade na direção X (kx)")
    axs[0].set_xlabel("X (células)")
    axs[0].set_ylabel("Y (células)")

    for w in wells:
        i, j = w['i'], w['j']
        marker, color = get_marker_style(w)
        axs[0].plot(j, i, marker=marker, color=color, markersize=10,
                    label=f"{w['tipo']} - {w['controle']}")

    fig.colorbar(im1, ax=axs[0], label="Permeabilidade [mD]")

    # --- Mapa ky ---
    im2 = axs[1].imshow(ky, cmap='viridis', origin='lower')
    axs[1].set_title("Permeabilidade na direção Y (ky)")
    axs[1].set_xlabel("X (células)")
    axs[1].set_ylabel("Y (células)")

    for w in wells:
        i, j = w['i'], w['j']
        marker, color = get_marker_style(w)
        axs[1].plot(j, i, marker=marker, color=color, markersize=10,
                    label=f"{w['tipo']} - {w['controle']}")

    fig.colorbar(im2, ax=axs[1], label="Permeabilidade [mD]")

    # Remover duplicatas da legenda
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='upper right')

    plt.tight_layout()
    plt.savefig(os.path.join(data_dir, "mapa_permeabilidade.png"), dpi=600)
    plt.show()



def plot_pressure_map_with_values(P, wells, data_dir):
    '''
    Plota o campo de pressão com os valores numéricos em cada célula e destaca a posição dos poços.
    O resultado é salvo como imagem no diretório especificado.

    Parâmetros:
    - P: matriz 2D com os valores de pressão (em kPa).
    - wells: lista de dicionários contendo as posições e informações dos poços:
        chaves: 'i', 'j', 'tipo' ('INJETOR' ou 'PRODUTOR'), 'controle' ('PRESSAO' ou 'VAZAO')
    - data_dir: caminho do diretório onde a imagem será salva.
    '''

    def get_marker_style(well):
        if well['tipo'] == 'INJETOR':
            marker = 'o'
            color = 'violet' if well['controle'] == 'VAZAO' else 'blue'
        elif well['tipo'] == 'PRODUTOR':
            marker = 'v'
            color = 'magenta' if well['controle'] == 'PRESSAO' else 'red'
        return marker, color

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(P, cmap='jet', origin='lower')
    cbar = plt.colorbar(im, ax=ax, label="Pressão [kPa]")

    NY, NX = P.shape
    for i in range(NY):
        for j in range(NX):
            ax.text(j, i, f"{P[i, j]:.1f}", color="white", ha='center', va='center', fontsize=6)

    for w in wells:
        i, j = w['i'], w['j']
        marker, color = get_marker_style(w)
        ax.plot(j, i, marker=marker, color=color, markersize=8,
                label=f"{w['tipo']} - {w['controle']}")

    # Remover duplicatas da legenda
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper right')

    ax.set_title("Mapa de Pressão com Valores [kPa]")
    ax.set_xlabel("X (células)")
    ax.set_ylabel("Y (células)")
    plt.tight_layout()
    plt.savefig(os.path.join(data_dir, "mapa_pressao_com_valores.png"), dpi=600)
    plt.show()

