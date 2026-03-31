"""
vsdock.pains
============
Filtro de compostos PAINS (Pan-Assay Interference Compounds) usando RDKit.
Remove moléculas com subestruturas problemáticas que geram falsos positivos
em ensaios biológicos.

Referência: Baell & Holloway, J. Med. Chem. 2010, 53, 2719-2740.
"""

import pandas as pd
from pathlib import Path
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def filter_pains(
    hits_file: str,
    outdir: str = "hits",
    also_filter_lipinski: bool = True,
) -> pd.DataFrame:
    """
    Filtra compostos PAINS (e opcionalmente viola Lipinski) de um arquivo de hits.

    Parâmetros
    ----------
    hits_file          : arquivo CSV com colunas smiles, id, tanimoto (saída do screen)
    outdir             : pasta para salvar resultados
    also_filter_lipinski : remove moléculas que violam mais de 1 regra de Lipinski

    Retorna
    -------
    DataFrame filtrado com coluna adicional 'pains_alert'
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Carrega hits
    df = pd.read_csv(hits_file)
    print(f"[pains] {len(df)} moléculas carregadas de {hits_file}")

    # Configura filtro PAINS (catálogos A, B e C)
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    catalog = FilterCatalog(params)

    # Avalia cada molécula
    pains_flags = []
    pains_alerts = []
    lipinski_flags = []

    for smiles in tqdm(df["smiles"], desc="Filtro PAINS", unit=" mols"):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            pains_flags.append(True)   # descarta inválidas
            pains_alerts.append("invalid_smiles")
            lipinski_flags.append(True)
            continue

        # PAINS
        entry = catalog.GetFirstMatch(mol)
        if entry:
            pains_flags.append(True)
            pains_alerts.append(entry.GetDescription())
        else:
            pains_flags.append(False)
            pains_alerts.append("")

        # Lipinski (opcional)
        if also_filter_lipinski:
            violations = _lipinski_violations(mol)
            lipinski_flags.append(violations > 1)
        else:
            lipinski_flags.append(False)

    df["is_pains"] = pains_flags
    df["pains_alert"] = pains_alerts
    df["lipinski_violation"] = lipinski_flags

    # Estatísticas
    n_pains = sum(pains_flags)
    n_lipinski = sum(lipinski_flags)
    print(f"[pains] PAINS removidos       : {n_pains}")
    if also_filter_lipinski:
        print(f"[pains] Lipinski removidos    : {n_lipinski - n_pains} adicionais")

    # Salva relatório completo (com flags)
    df.to_csv(outdir / "hits_pains_report.csv", index=False)

    # Filtra e salva apenas os aprovados
    mask = (~df["is_pains"]) & (~df["lipinski_violation"])
    df_clean = df[mask].drop(columns=["is_pains", "lipinski_violation"]).reset_index(drop=True)

    print(f"[pains] Moléculas aprovadas    : {len(df_clean)} / {len(df)}")

    df_clean.to_csv(outdir / "hits_clean.csv", index=False)
    df_clean[["smiles", "id"]].to_csv(
        outdir / "hits_clean.smi", sep="\t", index=False, header=False
    )

    print(f"[pains] Resultados salvos em:")
    print(f"  {outdir}/hits_pains_report.csv  (todos com flags)")
    print(f"  {outdir}/hits_clean.csv         (aprovados)")
    print(f"  {outdir}/hits_clean.smi         (aprovados, formato SMILES)")

    return df_clean


def _lipinski_violations(mol) -> int:
    """
    Conta violações às regras de Lipinski (Ro5).
    MW <= 500, HBA <= 10, HBD <= 5, logP <= 5
    """
    from rdkit.Chem import Descriptors
    violations = 0
    if Descriptors.MolWt(mol) > 500:
        violations += 1
    if Descriptors.NOCount(mol) > 10:
        violations += 1
    if Descriptors.NHOHCount(mol) > 5:
        violations += 1
    if Descriptors.MolLogP(mol) > 5:
        violations += 1
    return violations
