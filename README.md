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

## Uso

O vsdock é organizado em etapas sequenciais. Cada comando detecta os outputs
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
- `--ligand` : nome do ligante como aparece no PubChem (ex: `maraviroc`, `imatinib`)
- `--target` : nome do arquivo `.pdbqt` do receptor (será gerado na Etapa 2)

**Exemplo:**
```bash
vsdock init --target receptor.pdbqt --ligand maraviroc
```

---

### Etapa 2 — Preparar o receptor

O vsdock pode baixar a estrutura do PDB automaticamente ou usar um arquivo local.
Ele separa o receptor do ligante cocristalizado e converte tudo para PDBQT.

**A partir de um PDB ID (download automático):**
```bash
vsdock prepare --pdbid 4MBS --ligand-code MRV
```

**A partir de um arquivo PDB local:**
```bash
vsdock prepare --pdb-file estrutura.pdb --ligand-code MRV
```

**O que faz:**
- Baixa a estrutura do RCSB (se usar `--pdbid`)
- Lista todos os ligantes encontrados no PDB
- Separa receptor e ligante cocristalizado em arquivos distintos
- Converte o receptor para `.pdbqt` via Open Babel
- Atualiza o `vsdock_state.yaml` com os caminhos gerados

**Parâmetros:**
- `--pdbid` : PDB ID para download automático (ex: `4MBS`)
- `--pdb-file` : arquivo PDB local (alternativa ao `--pdbid`)
- `--ligand-code` : código de 3 letras do ligante no PDB (ex: `MRV`)

**Como descobrir o código do ligante:**
```bash
grep "HETATM" estrutura.pdb | awk '{print $4}' | sort -u
```

**Exemplo de output:**
```
[prepare] Ligantes encontrados no PDB:
  MRV1101A   ← use MRV como --ligand-code
  OLC1103A
  ...
[prepare] Pronto! Use nas próximas etapas:
  --receptor receptor.pdbqt
  --pdb      4MBS.pdb
  --autobox-ligand MRV1101A
```

---

### Etapa 3 — Baixar banco molecular

```bash
vsdock fetch --database --source chembl --max-mols 5000
```

**O que faz:**
- Baixa moléculas do ChEMBL filtradas por propriedades drug-like
- Salva em `zinc/chembl_database.smi`

**Parâmetros principais:**
- `--source` : `chembl` (padrão, recomendado) ou `zinc`
- `--max-mols` : número de moléculas a baixar (padrão: 500)
- `--mw-min` / `--mw-max` : intervalo de peso molecular (padrão: 150–500)
- `--logp-min` / `--logp-max` : intervalo de logP (padrão: −1 a 5)

---

### Etapa 4 — Triagem por similaridade

```bash
vsdock screen --similarity 0.4 --max-hits 500
```

**O que faz:**
- Compara o ligante de referência com o banco usando fingerprints de Morgan
- Retém compostos com Tanimoto ≥ threshold
- Salva os hits em `hits/hits.csv`

**Parâmetros:**
- `--similarity` : threshold de Tanimoto de 0 a 1 (padrão: 0.4)
- `--max-hits` : número máximo de hits (padrão: 500)
- `--fingerprint` : `morgan` (padrão), `maccs` ou `rdkit`

---

### Etapa 5 — Filtro PAINS e Lipinski

```bash
vsdock analyze
```

**O que faz:**
- Remove compostos com alertas PAINS
- Remove compostos que violam mais de uma regra de Lipinski
- Salva os aprovados em `hits/hits_clean.csv`

**Parâmetros:**
- `--no-lipinski` : desativa o filtro de Lipinski (útil para peptídeos/macrociclos)

---

### Etapa 6 — Docking molecular

```bash
vsdock dock --autobox-ligand MRV1101A
```

**O que faz:**
- Calcula a caixa de docking automaticamente usando o ligante cocristalizado
- Converte cada hit de SMILES para PDBQT
- Roda o AutoDock Vina para todos os hits + ligante de referência
- Salva os scores em `docking/docking_results.csv`

**Parâmetros principais:**
- `--autobox-ligand` : identificador completo do ligante no PDB (ex: `MRV1101A`)
  - Informado automaticamente pelo `vsdock prepare`
- `--exhaustiveness` : exaustividade do Vina (padrão: 8)
- `--top` : limita o número de hits a dockar (útil para testes)

**Alternativa com coordenadas manuais:**
```bash
vsdock dock --center X Y Z --size 20 20 20
```

---

### Etapa 7 — Análise de interações (PLIP)

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

### Etapa 8 — Predição ADMET

```bash
vsdock admet
```

**O que faz:**
- Prediz ~100 propriedades ADMET usando ADMET-AI (offline, sem GPU necessária)
- Calcula um score ponderado configurável por endpoint
- Salva em `admet/admet_results.csv` e `admet/admet_summary.csv`

**Parâmetros:**
- `--atc-code` : código ATC para comparação com DrugBank
  (ex: `J05` para antivirais, `C09` para cardiovascular)
- `--weights` : YAML com pesos personalizados por endpoint

**Personalizando os pesos** (edite `configs/default.yaml`):
```yaml
admet:
  hERG: 0.9              # cardiotoxicidade
  DILI: 0.8              # hepatotoxicidade
  Hepatotoxicity_Xu: 0.8
  AMES: 0.7              # mutagenicidade
  ClinTox: 0.8
  Bioavailability_Ma: 0.7
  Solubility_AqSolDB: 0.6
  BBB_Martini: 0.3       # menos relevante para antivirais periféricos
```

---

### Etapa 9 — Gerar manuscrito

```bash
vsdock report --format markdown
```

**O que faz:**
- Consolida todos os resultados em um manuscrito estruturado
- Preenche automaticamente a seção de métodos com os parâmetros reais usados
- Gera tabelas de docking, PLIP e ADMET prontas para o artigo
- Marca com `[TO COMPLETE]` as seções que precisam de redação manual

**Formatos:**
- `markdown` (padrão) → `report/manuscript.md`
- `quarto` → `report/manuscript.qmd`
- `html` → `report/manuscript.html`

---

## Fluxo completo

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

---

## Estrutura de pastas gerada

```
meu_projeto/
├── vsdock_state.yaml        ← estado do projeto (atualizado automaticamente)
├── 4MBS.pdb                 ← estrutura original
├── receptor_clean.pdb       ← receptor sem ligante
├── receptor.pdbqt           ← receptor pronto para docking
├── mrv_crystal.pdb          ← ligante cocristalizado isolado
├── zinc/
│   └── chembl_database.smi  ← banco molecular
├── hits/
│   ├── hits.csv             ← hits do screening
│   ├── hits_clean.csv       ← hits após filtro PAINS/Lipinski
│   └── hits_pains_report.csv
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
# Verifique se o ambiente conda está ativo
conda activate base
vsdock --version
```

**Docking com scores próximos de zero**
- Use sempre `vsdock prepare` para gerar o receptor e a caixa automaticamente
- Confirme que o `--autobox-ligand` corresponde ao ligante no PDB

**ADMET-AI falha com erro de GPU**
- Já tratado internamente — roda em CPU automaticamente
- Se persistir: `CUDA_VISIBLE_DEVICES="" vsdock admet`

**ZINC retorna erro SSL**
- Use `--source chembl` (padrão e recomendado)

**Ligante não encontrado no `vsdock prepare`**
- Verifique o código correto com:
  `grep "HETATM" estrutura.pdb | awk '{print $4}' | sort -u`
