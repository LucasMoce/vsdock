# vsdock

Pipeline de virtual screening e docking molecular.

## Instalação

```bash
git clone https://github.com/seu-usuario/vsdock
cd vsdock
pip install -e .
```

## Uso básico

```bash
vsdock init --target receptor.pdbqt --ligand maraviroc
vsdock screen --similarity 0.4
vsdock dock --exhaustiveness 8 --top 20
vsdock analyze
vsdock report --format quarto
```

## Módulos

| Módulo | Função |
|--------|--------|
| fetch  | Download de ligantes e bancos moleculares |
| screen | Triagem por similaridade (ChemFP/RDKit) |
| pains  | Filtro de compostos problemáticos |
| dock   | Docking com AutoDock Vina |
| plip   | Análise de interações proteína-ligante |
| admet  | Predição de propriedades ADMET |
| report | Geração de tabelas, figuras e manuscrito |
