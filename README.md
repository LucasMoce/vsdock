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

Isso instala o vsdock e todas as dependências Python automaticamente
(RDKit, pandas, requests, ADMET-AI, etc.).

Verifique se funcionou:

```bash
vsdock --version
```

---

## Uso

O vsdock é organizado em etapas sequenciais. Cada comando lê os outputs
do anterior automaticamente — você não precisa passar caminhos de arquivo
na maioria dos casos.

### Crie uma pasta para o seu projeto

```bash
mkdir meu_projeto
cd meu_projeto
```

---

### Etapa 1 — Inicializar o projeto

```bash
vsdock init --target receptor.pdbqt --ligand "nome do ligante"
```

**O que faz:**
- Baixa o SMILES do ligante de referência via PubChem
- Cria a estrutura de pastas do projeto
- Salva o estado em `vsdock_state.yaml`

**Parâmetros:**
- `--target` : arquivo `.pdbqt` do receptor (precisa estar na pasta do projeto)
- `--ligand` : nome do ligante como aparece no PubChem (ex: `maraviroc`, `imatinib`)

**Exemplo:**
```bash
vsdock init --target receptor.pdbqt --ligand maraviroc
```

> **Nota:** O receptor `.pdbqt` precisa ser preparado antes. Veja a seção
> "Preparando o receptor" mais abaixo.

---

### Etapa 2 — Baixar banco molecular

```bash
vsdock fetch --database --source chembl --max-mols 5000
```

**O que faz:**
- Baixa moléculas do ChEMBL filtradas por peso molecular e logP
- Salva em `zinc/chembl_database.smi`

**Parâmetros principais:**
- `--source` : `chembl` (padrão, recomendado) ou `zinc`
- `--max-mols` : número de moléculas a baixar (padrão: 500)
- `--mw-min` / `--mw-max` : intervalo de peso molecular (padrão: 150–500)
- `--logp-min` / `--logp-max` : intervalo de logP (padrão: −1 a 5)

**Exemplo com filtros:**
```bash
vsdock fetch --database --source chembl --max-mols 2000 --mw-min 200 --mw-max 450
```

---

### Etapa 3 — Triagem por similaridade

```bash
vsdock screen --similarity 0.4 --max-hits 500
```

**O que faz:**
- Compara o ligante de referência com todo o banco usando fingerprints de Morgan
- Seleciona os compostos mais similares (Tanimoto ≥ threshold)
- Salva os hits em `hits/hits.csv` e `hits/hits.smi`

**Parâmetros principais:**
- `--similarity` : threshold de Tanimoto de 0 a 1 (padrão: 0.4)
- `--max-hits` : número máximo de hits (padrão: 500)
- `--fingerprint` : tipo de fingerprint — `morgan` (padrão), `maccs` ou `rdkit`

---

### Etapa 4 — Filtro PAINS e Lipinski

```bash
vsdock analyze
```

**O que faz:**
- Remove compostos com alertas PAINS (interferentes de ensaios)
- Remove compostos que violam mais de uma regra de Lipinski
- Salva os aprovados em `hits/hits_clean.csv`

**Parâmetros:**
- `--no-lipinski` : desativa o filtro de Lipinski (útil para peptídeos/macrociclos)

---

### Etapa 5 — Docking molecular

```bash
vsdock dock --pdb estrutura.pdb --autobox-ligand COD1101A --receptor receptor.pdbqt
```

**O que faz:**
- Calcula automaticamente a caixa de docking usando o ligante cocristalizado (autobox)
- Converte cada hit de SMILES para PDBQT
- Roda o AutoDock Vina para todos os hits + ligante de referência
- Salva os scores em `docking/docking_results.csv`

**Parâmetros principais:**
- `--pdb` : arquivo PDB original com o ligante cocristalizado
- `--autobox-ligand` : código do ligante no PDB (ex: `MRV1101A`)
  - Formato: `NOME` + `NÚMERO` + `CADEIA` (ex: `MRV1101A` = ligante MRV, ID 1101, cadeia A)
  - Para descobrir o código: `grep "HETATM" estrutura.pdb | awk '{print $4, $5, $6}' | sort -u`
- `--receptor` : arquivo `.pdbqt` do receptor
- `--exhaustiveness` : exaustividade do Vina (padrão: 8; aumente para resultados mais precisos)
- `--top` : limita o número de hits a dockar (útil para testes rápidos)

**Alternativa sem PDB (coordenadas manuais):**
```bash
vsdock dock --center X Y Z --size 20 20 20 --receptor receptor.pdbqt
```

---

### Etapa 6 — Análise de interações (PLIP)

```bash
vsdock plip --top-n 10
```

**O que faz:**
- Analisa as interações proteína-ligante das melhores poses
- Identifica ligações de hidrogênio, contatos hidrofóbicos, π-stacking, etc.
- Salva em `plip/plip_interactions.csv` e `plip/plip_summary.csv`

**Parâmetros:**
- `--top-n` : número de melhores compostos a analisar (padrão: 10)

---

### Etapa 7 — Predição ADMET

```bash
vsdock admet
```

**O que faz:**
- Prediz ~100 propriedades ADMET usando ADMET-AI (modelos locais, offline)
- Calcula um score ponderado por endpoint (configurável)
- Salva em `admet/admet_results.csv` e `admet/admet_summary.csv`

**Parâmetros:**
- `--atc-code` : código ATC para comparação com DrugBank (ex: `J05` para antivirais,
  `C09` para cardiovascular). Veja códigos em [whocc.no/atc](https://www.whocc.no/atc/)
- `--weights` : arquivo YAML com pesos personalizados por endpoint

**Personalizando os pesos ADMET** (edite `configs/default.yaml`):
```yaml
admet:
  hERG: 0.9              # cardiotoxicidade — peso alto
  DILI: 0.8              # hepatotoxicidade
  Hepatotoxicity_Xu: 0.8
  AMES: 0.7              # mutagenicidade
  ClinTox: 0.8
  Bioavailability_Ma: 0.7
  Solubility_AqSolDB: 0.6
  BBB_Martini: 0.3       # barreira hematoencefálica — menos relevante para antivirais
```

---

### Etapa 8 — Gerar manuscrito

```bash
vsdock report --format markdown
```

**O que faz:**
- Consolida todos os resultados em um manuscrito estruturado
- Preenche automaticamente a seção de métodos com os parâmetros usados
- Gera tabelas de docking, PLIP e ADMET prontas para o artigo
- Marca com `[TO COMPLETE]` as seções que precisam de redação manual

**Formatos:**
- `markdown` (padrão) → `report/manuscript.md`
- `quarto` → `report/manuscript.qmd` (renderizável com [Quarto](https://quarto.org))
- `html` → `report/manuscript.html`

---

## Preparando o receptor

O receptor precisa estar em formato `.pdbqt`. O jeito mais simples:

```bash
# 1. Baixe a estrutura do PDB (ex: 4MBS)
# Em https://www.rcsb.org — Download Files → PDB Format

# 2. Extraia apenas a proteína (remove água e ligantes)
grep "^ATOM" estrutura.pdb > receptor_clean.pdb

# 3. Converta para PDBQT
obabel receptor_clean.pdb -O receptor.pdbqt -xr
```

> **Recomendação:** Para resultados mais precisos, use o
> [MGLTools](https://ccsb.scripps.edu/mgltools/) para preparar o receptor,
> pois ele calcula cargas parciais de forma mais robusta.

---

## Estrutura de pastas gerada

Após rodar o pipeline completo, seu projeto terá esta estrutura:

```
meu_projeto/
├── vsdock_state.yaml       ← estado do projeto (parâmetros, SMILES, receptor)
├── receptor.pdbqt          ← receptor preparado
├── 4MBS.pdb                ← estrutura original (para autobox)
├── zinc/
│   └── chembl_database.smi ← banco molecular baixado
├── hits/
│   ├── hits.csv            ← hits do screening
│   ├── hits_clean.csv      ← hits após filtro PAINS/Lipinski
│   └── hits_pains_report.csv
├── docking/
│   ├── docking_results.csv ← scores ranqueados
│   ├── pdbqt/              ← estruturas 3D dos ligantes
│   └── poses/              ← poses geradas pelo Vina
├── plip/
│   ├── plip_interactions.csv
│   └── plip_summary.csv
├── admet/
│   ├── admet_results.csv   ← predições completas (106 propriedades)
│   └── admet_summary.csv   ← endpoints principais
└── report/
    └── manuscript.md       ← manuscrito gerado automaticamente
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

---

## Problemas comuns

**`vsdock: command not found` após instalação**
```bash
# Verifique se o ambiente conda está ativo
conda activate base
vsdock --version
```

**Docking com scores próximos de zero**
- Verifique as coordenadas da caixa — use sempre `--autobox-ligand` com o PDB original
- Confirme que o receptor está na mesma posição do ligante cocristalizado

**ADMET-AI falha com erro de GPU**
- Já tratado internamente — roda em CPU automaticamente
- Se persistir: `CUDA_VISIBLE_DEVICES="" vsdock admet`

**ZINC retorna erro SSL**
- Use `--source chembl` (padrão e recomendado)
