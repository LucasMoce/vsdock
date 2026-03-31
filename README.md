# vsdock

Pipeline de virtual screening e docking molecular por linha de comando.

## Requisitos

Antes de instalar o vsdock, você precisa ter os seguintes programas instalados:

- Python 3.9 ou superior (recomendado via [Miniconda](https://docs.conda.io/en/latest/miniconda.html))
- [AutoDock Vina 1.2+](https://github.com/ccsb-scripps/AutoDock-Vina/releases)
- [Open Babel 3.x](http://openbabel.org/wiki/Category:Installation)
- [PLIP](https://github.com/pharmai/plip)

### Instalando as dependências no Ubuntu/Debian

```bash
sudo apt install autodock-vina openbabel
pip install plip
```

### Instalando as dependências no macOS (via Homebrew)

```bash
brew install open-babel
pip install plip
# AutoDock Vina: baixe o binário em https://github.com/ccsb-scripps/AutoDock-Vina/releases
```

---

## Instalação do vsdock

```bash
pip install git+https://github.com/LucasMoce/vsdock.git
```

Verifique se funcionou:

```bash
vsdock --version
```

---

## Modos de uso

O vsdock suporta dois modos dependendo do que você tem disponível:

| | Modo com ligante de referência | Modo receptor-only |
|---|---|---|
| Tem ligante cocristalizado? | ✅ Sim | ❌ Não |
| Etapa de screening | ✅ Por similaridade | ⏭ Pulada automaticamente |
| Definição do sítio | Ligante cocristalizado | Blind docking ou resíduos |
| Score de referência no relatório | ✅ Sim | ❌ Não |

---

## Fluxo com ligante de referência

Use quando você tem um medicamento ou composto conhecido que se liga ao receptor.

```bash
mkdir meu_projeto && cd meu_projeto

vsdock init --target receptor.pdbqt --ligand maraviroc
vsdock prepare --pdbid 4MBS --ligand-code MRV
vsdock fetch --database --source chembl --max-mols 2000
vsdock screen --similarity 0.4 --max-hits 500
vsdock analyze
vsdock dock --autobox-ligand MRV1101A
vsdock plip --top-n 10
vsdock admet --atc-code J05
vsdock report --format markdown
```

## Fluxo sem ligante de referência (receptor-only)

Use quando você tem apenas a estrutura do receptor e quer descobrir ligantes do zero.

```bash
mkdir meu_projeto && cd meu_projeto

vsdock init --target receptor.pdbqt
vsdock prepare --pdbid 4MBS
vsdock fetch --database --source chembl --max-mols 2000
# etapa screen é pulada automaticamente
vsdock dock --blind                              # blind docking
# ou: vsdock dock --autobox-residues VAL100A TRP248A  # por resíduos
vsdock plip --top-n 10
vsdock admet
vsdock report --format markdown
```

---

## Referência de comandos

### `vsdock init`

Inicializa o projeto e cria a estrutura de pastas.

```bash
vsdock init --target receptor.pdbqt --ligand maraviroc   # com ligante
vsdock init --target receptor.pdbqt                      # sem ligante
```

| Parâmetro | Descrição |
|-----------|-----------|
| `--target` | Nome do arquivo `.pdbqt` do receptor (gerado pelo `prepare`) |
| `--ligand` | Nome do ligante de referência no PubChem. **Omitir para modo receptor-only** |

---

### `vsdock prepare`

Baixa e prepara o receptor a partir de um PDB ID ou arquivo local.
Separa receptor e ligante, converte para PDBQT e atualiza o estado do projeto.

```bash
vsdock prepare --pdbid 4MBS --ligand-code MRV    # com ligante
vsdock prepare --pdbid 4MBS                      # sem ligante (receptor-only)
vsdock prepare --pdb-file estrutura.pdb --ligand-code MRV  # arquivo local
```

| Parâmetro | Descrição |
|-----------|-----------|
| `--pdbid` | PDB ID para download automático (ex: `4MBS`) |
| `--pdb-file` | Arquivo PDB local (alternativa ao `--pdbid`) |
| `--ligand-code` | Código de 3 letras do ligante no PDB (ex: `MRV`). **Omitir para modo receptor-only** |

**Como descobrir o código do ligante:**
```bash
grep "HETATM" estrutura.pdb | awk '{print $4}' | sort -u
```

O comando lista os ligantes e informa automaticamente o que usar nos próximos passos:
```
[prepare] Ligantes encontrados no PDB:
  MRV1101A   ← use MRV como --ligand-code
  OLC1103A
[prepare] Pronto! Use nas próximas etapas:
  --receptor receptor.pdbqt
  --pdb      4MBS.pdb
  --autobox-ligand MRV1101A
```

---

### `vsdock fetch`

Baixa banco de moléculas para screening/docking.

```bash
vsdock fetch --database --source chembl --max-mols 2000
```

| Parâmetro | Descrição | Padrão |
|-----------|-----------|--------|
| `--source` | `chembl` (recomendado) ou `zinc` | `chembl` |
| `--max-mols` | Número de moléculas a baixar | `500` |
| `--mw-min` / `--mw-max` | Intervalo de peso molecular | `150` / `500` |
| `--logp-min` / `--logp-max` | Intervalo de logP | `-1` / `5` |

---

### `vsdock screen`

Triagem por similaridade estrutural. **Requer ligante de referência.**
No modo receptor-only, este comando informa que a etapa será pulada.

```bash
vsdock screen --similarity 0.4 --max-hits 500
```

| Parâmetro | Descrição | Padrão |
|-----------|-----------|--------|
| `--similarity` | Threshold de Tanimoto (0–1) | `0.4` |
| `--max-hits` | Máximo de hits a retornar | `500` |
| `--fingerprint` | `morgan`, `maccs` ou `rdkit` | `morgan` |

---

### `vsdock analyze`

Filtro de PAINS e regra de Lipinski. Recomendado no modo com ligante.
No modo receptor-only pode ser pulado (o banco vai direto para o docking).

```bash
vsdock analyze
vsdock analyze --no-lipinski   # desativa filtro de Lipinski (peptídeos, macrociclos)
```

---

### `vsdock dock`

Docking molecular com AutoDock Vina.

```bash
# Com ligante cocristalizado (caixa automática)
vsdock dock --autobox-ligand MRV1101A

# Por resíduos do sítio ativo
vsdock dock --autobox-residues VAL100A TRP248A LYS305A

# Blind docking (cobre a proteína toda)
vsdock dock --blind

# Coordenadas manuais
vsdock dock --center X Y Z --size 20 20 20
```

| Parâmetro | Descrição | Padrão |
|-----------|-----------|--------|
| `--autobox-ligand` | ID do ligante cocristalizado (ex: `MRV1101A`) | — |
| `--autobox-residues` | Resíduos do sítio ativo (ex: `VAL100A TRP248A`) | — |
| `--blind` | Blind docking — caixa cobre a proteína toda | — |
| `--center X Y Z` | Centro da caixa manual | — |
| `--size SX SY SZ` | Tamanho da caixa em Å | `20 20 20` |
| `--exhaustiveness` | Exaustividade do Vina (maior = mais preciso) | `8` |
| `--top` | Limita o número de moléculas a dockar | todos |

> **Nota sobre o identificador do ligante:** o formato é `NOME` + `NÚMERO` + `CADEIA`
> (ex: `MRV1101A` = ligante MRV, ID 1101, cadeia A).
> O `vsdock prepare` informa o identificador correto automaticamente.

---

### `vsdock plip`

Análise de interações proteína-ligante nas melhores poses.

```bash
vsdock plip --top-n 10
```

| Parâmetro | Descrição | Padrão |
|-----------|-----------|--------|
| `--top-n` | Número de melhores compostos a analisar | `10` |

---

### `vsdock admet`

Predição de ~100 propriedades ADMET usando ADMET-AI (offline, sem GPU necessária).

```bash
vsdock admet
vsdock admet --atc-code J05   # compara com antivirais aprovados no DrugBank
```

| Parâmetro | Descrição |
|-----------|-----------|
| `--atc-code` | Código ATC para comparação com DrugBank (ex: `J05` antivirais, `C09` cardiovascular) |
| `--weights` | YAML com pesos personalizados por endpoint |

**Personalizando os pesos** (edite `configs/default.yaml`):
```yaml
admet:
  hERG: 0.9              # cardiotoxicidade
  DILI: 0.8              # hepatotoxicidade
  AMES: 0.7              # mutagenicidade
  ClinTox: 0.8
  Bioavailability_Ma: 0.7
  Solubility_AqSolDB: 0.6
  BBB_Martini: 0.3
```

---

### `vsdock report`

Gera manuscrito estruturado consolidando todos os resultados.

```bash
vsdock report --format markdown   # → report/manuscript.md
vsdock report --format quarto     # → report/manuscript.qmd
vsdock report --format html       # → report/manuscript.html
```

A seção de métodos é preenchida automaticamente com os parâmetros usados.
As seções marcadas com `[TO COMPLETE]` precisam de redação manual
(introdução, discussão, conclusão).

---

## Estrutura de pastas gerada

```
meu_projeto/
├── vsdock_state.yaml        ← estado do projeto (atualizado automaticamente)
├── 4MBS.pdb                 ← estrutura original
├── receptor_clean.pdb       ← receptor sem água
├── receptor.pdbqt           ← receptor pronto para docking
├── mrv_crystal.pdb          ← ligante cocristalizado (se houver)
├── zinc/
│   └── chembl_database.smi  ← banco molecular
├── hits/
│   ├── hits.csv             ← hits do screening (modo com ligante)
│   ├── hits_clean.csv       ← hits após filtro PAINS/Lipinski
│   └── hits_from_bank.csv   ← banco convertido (modo receptor-only)
├── docking/
│   ├── docking_results.csv  ← scores ranqueados
│   ├── pdbqt/               ← estruturas 3D dos ligantes
│   └── poses/               ← poses geradas pelo Vina
├── plip/
│   ├── plip_interactions.csv
│   └── plip_summary.csv
├── admet/
│   ├── admet_results.csv    ← predições completas (~100 propriedades)
│   └── admet_summary.csv    ← endpoints principais
└── report/
    └── manuscript.md        ← manuscrito gerado automaticamente
```

---

## Citações

Se usar o vsdock num artigo, cite os programas subjacentes:

- **AutoDock Vina:** Eberhardt et al., *J. Chem. Inf. Model.* 2021
- **PLIP:** Salentin et al., *Nucleic Acids Res.* 2015
- **ADMET-AI:** Swanson et al., *Bioinformatics* 2024
- **RDKit:** [rdkit.org](https://www.rdkit.org)
- **Open Babel:** O'Boyle et al., *J. Cheminform.* 2011
- **autobox:** [github.com/omixlab/autobox](https://github.com/omixlab/autobox)
- **ChEMBL:** Mendez et al., *Nucleic Acids Res.* 2019

---

## Problemas comuns

**`vsdock: command not found` após instalação**
```bash
conda activate base
vsdock --version
```

**Docking com scores próximos de zero**
- Use sempre `vsdock prepare` para garantir que receptor e caixa estão corretos
- Confirme que `--autobox-ligand` corresponde ao ligante no PDB

**ADMET-AI falha com erro de GPU**
- Já tratado internamente — roda em CPU automaticamente
- Se persistir: `CUDA_VISIBLE_DEVICES="" vsdock admet`

**ZINC retorna erro SSL**
- Use `--source chembl` (padrão e recomendado)

**Ligante não encontrado no `vsdock prepare`**
```bash
grep "HETATM" estrutura.pdb | awk '{print $4}' | sort -u
```
