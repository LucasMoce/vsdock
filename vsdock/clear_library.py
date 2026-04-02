"""
vsdock.clear_library
====================
Filtragem e limpeza da biblioteca de compostos.

Filtros disponíveis:
  --pains     Remove compostos PAINS (Pan-Assay Interference)
  --lipinski  Remove violações da Regra de 5 de Lipinski (>1 violação)
  --pfizer    Filtro de Pfizer (MW≤500, logP≤5, HBD≤5, HBA≤10, RotBonds≤10, TPSA≤150)
  --gsk       Filtro de GSK (MW≤400, logP≤4, TPSA≤140, sem fluor excessivo)

Referências:
  - PAINS: Baell & Holloway, J. Med. Chem. 2010
  - Pfizer: Gleeson, J. Med. Chem. 2008
  - GSK: Rishton, Drug Discov. Today 2003
"""

import pandas as pd
from pathlib import Path
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


# ---------------------------------------------------------------------------
# Filtros individuais
# ---------------------------------------------------------------------------

def _setup_pains_catalog():
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    return FilterCatalog(params)


def _lipinski_violations(mol) -> int:
    """Conta violações da Regra de 5 de Lipinski."""
    violations = 0
    if Descriptors.MolWt(mol) > 500:       violations += 1
    if Descriptors.NOCount(mol) > 10:       violations += 1
    if Descriptors.NHOHCount(mol) > 5:      violations += 1
    if Descriptors.MolLogP(mol) > 5:        violations += 1
    return violations


def _pfizer_filter(mol) -> tuple:
    """
    Filtro de Pfizer (CNS MPO-inspired, Gleeson 2008).
    Retorna (passa, motivo).
    MW≤500, logP≤5, HBD≤5, HBA≤10, RotBonds≤10, TPSA≤150
    """
    mw       = Descriptors.MolWt(mol)
    logp     = Descriptors.MolLogP(mol)
    hbd      = rdMolDescriptors.CalcNumHBD(mol)
    hba      = rdMolDescriptors.CalcNumHBA(mol)
    rot      = rdMolDescriptors.CalcNumRotatableBonds(mol)
    tpsa     = rdMolDescriptors.CalcTPSA(mol)

    failures = []
    if mw   > 500:  failures.append(f"MW={mw:.0f}>500")
    if logp > 5:    failures.append(f"logP={logp:.1f}>5")
    if hbd  > 5:    failures.append(f"HBD={hbd}>5")
    if hba  > 10:   failures.append(f"HBA={hba}>10")
    if rot  > 10:   failures.append(f"RotBonds={rot}>10")
    if tpsa > 150:  failures.append(f"TPSA={tpsa:.0f}>150")

    return (len(failures) == 0, "; ".join(failures))


def _gsk_filter(mol) -> tuple:
    """
    Filtro de GSK (Rishton 2003, refinado).
    Retorna (passa, motivo).
    MW≤400, logP≤4, TPSA≤140, fluoretos≤4, sem grupos reativos
    """
    mw   = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)

    # Conta fluor
    n_fluorine = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 9)

    failures = []
    if mw   > 400:  failures.append(f"MW={mw:.0f}>400")
    if logp > 4:    failures.append(f"logP={logp:.1f}>4")
    if tpsa > 140:  failures.append(f"TPSA={tpsa:.0f}>140")
    if n_fluorine > 4: failures.append(f"F={n_fluorine}>4")

    return (len(failures) == 0, "; ".join(failures))


# ---------------------------------------------------------------------------
# Função principal
# ---------------------------------------------------------------------------

def filter_library(
    hits_file: str,
    outdir: str = "hits",
    apply_pains: bool = True,
    apply_lipinski: bool = True,
    apply_pfizer: bool = False,
    apply_gsk: bool = False,
) -> pd.DataFrame:
    """
    Filtra a biblioteca de compostos aplicando os filtros selecionados.

    Parâmetros
    ----------
    hits_file      : CSV com colunas smiles, id (saída do similarity-search)
    outdir         : pasta de saída
    apply_pains    : remove compostos PAINS
    apply_lipinski : remove violações de Lipinski (>1 violação)
    apply_pfizer   : aplica filtro de Pfizer
    apply_gsk      : aplica filtro de GSK

    Retorna
    -------
    DataFrame com compostos aprovados
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(hits_file)
    print(f"[clear-library] {len(df)} moléculas carregadas de {hits_file}")

    active_filters = []
    if apply_pains:    active_filters.append("PAINS")
    if apply_lipinski: active_filters.append("Lipinski")
    if apply_pfizer:   active_filters.append("Pfizer")
    if apply_gsk:      active_filters.append("GSK")
    print(f"[clear-library] Filtros ativos: {', '.join(active_filters) if active_filters else 'nenhum'}")

    # Setup PAINS
    catalog = _setup_pains_catalog() if apply_pains else None

    # Colunas de resultado
    results = {
        "is_pains":         [],
        "pains_alert":      [],
        "lipinski_violation": [],
        "pfizer_fail":      [],
        "pfizer_reason":    [],
        "gsk_fail":         [],
        "gsk_reason":       [],
    }

    for smiles in tqdm(df["smiles"], desc="Filtrando", unit=" mols"):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            results["is_pains"].append(True)
            results["pains_alert"].append("invalid_smiles")
            results["lipinski_violation"].append(True)
            results["pfizer_fail"].append(True)
            results["pfizer_reason"].append("invalid_smiles")
            results["gsk_fail"].append(True)
            results["gsk_reason"].append("invalid_smiles")
            continue

        # PAINS
        if apply_pains:
            entry = catalog.GetFirstMatch(mol)
            results["is_pains"].append(bool(entry))
            results["pains_alert"].append(entry.GetDescription() if entry else "")
        else:
            results["is_pains"].append(False)
            results["pains_alert"].append("")

        # Lipinski
        results["lipinski_violation"].append(
            _lipinski_violations(mol) > 1 if apply_lipinski else False
        )

        # Pfizer
        if apply_pfizer:
            passes, reason = _pfizer_filter(mol)
            results["pfizer_fail"].append(not passes)
            results["pfizer_reason"].append(reason)
        else:
            results["pfizer_fail"].append(False)
            results["pfizer_reason"].append("")

        # GSK
        if apply_gsk:
            passes, reason = _gsk_filter(mol)
            results["gsk_fail"].append(not passes)
            results["gsk_reason"].append(reason)
        else:
            results["gsk_fail"].append(False)
            results["gsk_reason"].append("")

    for col, vals in results.items():
        df[col] = vals

    # Estatísticas
    n_pains    = df["is_pains"].sum()
    n_lipinski = df["lipinski_violation"].sum()
    n_pfizer   = df["pfizer_fail"].sum()
    n_gsk      = df["gsk_fail"].sum()

    if apply_pains:    print(f"[clear-library] PAINS removidos       : {n_pains}")
    if apply_lipinski: print(f"[clear-library] Lipinski removidos    : {n_lipinski}")
    if apply_pfizer:   print(f"[clear-library] Pfizer reprovados     : {n_pfizer}")
    if apply_gsk:      print(f"[clear-library] GSK reprovados        : {n_gsk}")

    # Salva relatório completo
    df.to_csv(outdir / "hits_filter_report.csv", index=False)

    # Filtra aprovados
    mask = (
        (~df["is_pains"]) &
        (~df["lipinski_violation"]) &
        (~df["pfizer_fail"]) &
        (~df["gsk_fail"])
    )
    df_clean = df[mask].drop(
        columns=["is_pains", "lipinski_violation", "pfizer_fail", "gsk_fail"],
        errors="ignore"
    ).reset_index(drop=True)

    print(f"[clear-library] Moléculas aprovadas    : {len(df_clean)} / {len(df)}")

    df_clean.to_csv(outdir / "hits_clean.csv", index=False)
    df_clean[["smiles", "id"]].to_csv(
        outdir / "hits_clean.smi", sep="\t", index=False, header=False
    )

    print(f"[clear-library] Resultados salvos em:")
    print(f"  {outdir}/hits_filter_report.csv  (todos com flags)")
    print(f"  {outdir}/hits_clean.csv          (aprovados)")
    print(f"  {outdir}/hits_clean.smi          (aprovados, SMILES)")

    return df_clean


# Alias para compatibilidade com código antigo
def filter_pains(hits_file, outdir="hits", also_filter_lipinski=True):
    return filter_library(
        hits_file=hits_file,
        outdir=outdir,
        apply_pains=True,
        apply_lipinski=also_filter_lipinski,
    )
