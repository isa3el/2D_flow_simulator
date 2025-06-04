# Simulador de Escoamento Monofásico 2D

Trabalho 01 da disciplina MAT2490
Aluna Isabel Goncalves - 2312237

Este simulador modela o escoamento monofásico incompressível em um meio poroso bidimensional (2D).

## Matemática por Trás do Simulador

O simulador resolve numericamente a **equação de escoamento monofásico incompressível** em um meio poroso, baseada na **Lei de Darcy** e no **princípio de conservação de massa**.

### Equação governante (forma contínua):

$$
\nabla \cdot \mathbf{q} = q_s
\quad\text{com}\quad
\mathbf{q} = -\frac{\mathbf{K}}{\mu} \nabla p
$$

Onde:

* $\mathbf{q}$: vetor de fluxo (m³/dia)
* $\mathbf{K}$: tensor de permeabilidade (kx e ky)
* $\mu$: viscosidade do fluido (cP)
* $p$: pressão (kPa)
* $q_s$: fonte ou sumidouro (vazão de poço)

### Discretização por Diferenças Finitas

A malha é estruturada, e a equação é discretizada usando a formulação de **diferenças finitas com médias harmônicas** para transmissibilidade nas interfaces entre blocos:

#### Transmissibilidade entre dois blocos adjacentes:

$$
T_{ij} = \frac{2 \cdot k_{\text{eff}} \cdot h}{\mu \cdot \Delta x}
\quad\text{com}\quad
k_{\text{eff}} = \frac{2 k_i k_j}{k_i + k_j}
$$

### Sistema Linear

A discretização resulta em um sistema linear:

$$
\mathbf{T} \cdot \mathbf{P} = \mathbf{Q}
$$

* $\mathbf{T}$: matriz de transmissibilidade (esparsa)
* $\mathbf{P}$: vetor de pressões (kPa)
* $\mathbf{Q}$: vetor de fontes (vazões de injeção ou extração)

### Condições nos Poços

Cada poço pode ser controlado por:

* **Pressão imposta**: substitui a equação da célula pela condição $p = p_{\text{poço}}$
* **Vazão imposta**: valor inserido diretamente no vetor $\mathbf{Q}$

### Observações Finais

* As médias harmônicas são usadas para preservar o comportamento físico em transições bruscas de permeabilidade.
* A matriz $\mathbf{T}$ é montada com esquema de 5 pontos (célula e seus quatro vizinhos).
* A solução é feita via solver direto do `scipy.sparse.linalg`.

### Cálculo da Vazão nos Poços

A vazão de cada poço é calculada pelo **Well Index (WI)**: baseado na permeabilidade anisotrópica (kx, ky), dimensões da célula (dx, dy), espessura do reservatório (h), raio do poço (rw) e fator de dano (skin). A vazão é então dada por:


$q = \frac{WI \cdot (p_{\text{poço}} - p_{\text{célula}})}{\mu}$

com 

$WI = \frac{2 \pi k h}{\ln\left(\frac{r_e}{r_w}\right) + s}$

---

## Funcionalidades Principais

* **Malha 2D configurável**:

  * Dimensões arbitrárias (exemplos criados: 5×5, 15×15, 60×60)
  * Tamanho individual das células em X (`dx`) e Y (`dy`)
* **Mapa de células ativas**: suporte a zonas inativas no domínio simulado
* **Poços com controle flexível**:

  * Tipos: `INJETOR` ou `PRODUTOR`
  * Modos de controle: por `PRESSAO` (kPa) ou `VAZAO` (m³/dia)
  * Possibilidade de adicionar poços em qualquer célula ativa
* **Propriedades anisotrópicas**:

  * Permeabilidades `kx` e `ky` independentes
* **Parâmetros físicos configuráveis**:

  * Espessura (`h`), viscosidade (`mu`), densidade (`rho`)
* **Cálculo automático de vazão nos poços**, com base no campo de pressão simulado

## Entradas Esperadas

* `input.txt` informações sobre os reservatórios e os poços

```
[NX NY]
dimensões X e Y da malha
[DX]
tamanho dos blocos na direção X. Entradas aceitas: 3*50 ou 50 50 50 ou 10 40 50  
[DY]
tamanho dos blocos na direção X. Entradas aceitas: 3*50 ou 50 50 50 ou 10 40 50
[H]
profundidade do reservatório em M.

[PROPS]
viscosidade (mu) e densidade (rho)

[WELLS]
número de poços
tipo do poço (INJ ou PRO), coordenada x, coordenada y, tipo de controle, valor do controle. Exemplos:
INJ 1 4 VAZAO 600.0
PRO 2 3 PRESSAO 320.0
```

* `grid.txt`: células ativas (0 ou 1)

* `perm_x.txt`: matriz de permeabilidade na direção X (kx)

* `perm_y.txt`: matriz de permeabilidade na direção Y (ky)

* Optei por ter arquivos separados para grid, perm\_x e perm\_y pois queria testar o caso 60x60. Ficaria confuso tudo em um único .txt de entrada.

## Saídas Geradas

* **Campo de pressão** (matriz 2D) em kPa (exibida apenas para malhas com `nX * nY < 400`)
* **Tabela de poços** com localização, controle, Pressão imposta ou simulada e Vazão resultante
* **Gráficos com a localização dos poços**: Campo de pressão, Mapa de permeabilidade `kx` e `ky` e Campo de pressão com valores
* **Arquivos exportados**:

  * `resultados.txt`: contém todas as informações conhecidas do reservatório, os dados de produção dos poços, mapa de pressão e mapas de permeabilidade nas direções X e Y
  * `mapa_pressao.png`, `mapa_permeabilidade.png` e `mapa_pressao_com_valores.png`

---



## Estrutura do Repositório

```
├─ README.md
├─ plot_functions.py
├─ simulation_functions.py
├─ simulation.py
├─ Grid_5_5
│  ├─ grid.txt
│  ├─ input.txt
│  ├─ perm_x.txt
│  ├─ perm_y.txt
│  ├─ mapa_permeabilidade.png
│  ├─ mapa_pressao_com_valores.png
│  ├─ mapa_pressao.png
│  ├─ results.txt
├─ Grid_15_15
│  ├─ ...
├─ Grid_60_60
│  ├─ ...
```

## Arquivos Principais

* `plot_functions.py`: funções para gerar e salvar os gráficos
* `simulation_functions.py`: funções para ler os arquivos de entrada e resolver o sistema
* `simulation.py`: script principal que executa o simulador

## Ambiente Utilizado

* Python 3.12.9
* `numpy` 1.26.4
* `matplotlib` 3.10.1
* `scipy` 1.13.1
