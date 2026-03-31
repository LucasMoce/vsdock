"""
vsdock.fetch
============
Download do ligante de referência (PubChem) e de bancos moleculares (ZINC20).
"""

import os
import time
import requests
from pathlib import Path
from tqdm import tqdm


# ---------------------------------------------------------------------------
# PubChem
# ---------------------------------------------------------------------------

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def fetch_ligand(name: str, outdir: str = ".", fmt: str = "sdf") -> Path:
    """
    Baixa a estrutura 2D/3D de um ligante pelo nome via PubChem.

    Parâmetros
    ----------
    name   : nome do composto (ex: "maraviroc", "imatinib")
    outdir : pasta de destino
    fmt    : formato de saída — "sdf" ou "smiles"

    Retorna
    -------
    Path do arquivo baixado
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[fetch] Buscando '{name}' no PubChem...")

    # 1. Resolve o nome para CID
    resp = requests.get(
        f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(name)}/cids/JSON",
        timeout=15,
    )
    resp.raise_for_status()
    cid = resp.json()["IdentifierList"]["CID"][0]
    print(f"[fetch] CID encontrado: {cid}")

    # 2. Baixa a estrutura no formato desejado
    if fmt == "sdf":
        url = f"{PUBCHEM_BASE}/compound/cid/{cid}/SDF?record_type=3d"
        suffix = ".sdf"
    else:
        url = f"{PUBCHEM_BASE}/compound/cid/{cid}/property/IsomericSMILES/TXT"
        suffix = ".smi"

    resp = requests.get(url, timeout=15)
    resp.raise_for_status()

    outfile = outdir / f"{name.lower().replace(' ', '_')}{suffix}"
    outfile.write_bytes(resp.content)
    print(f"[fetch] Ligante salvo em: {outfile}")
    return outfile


def fetch_ligand_smiles(name: str) -> str:
    """
    Retorna apenas o SMILES canônico do ligante (sem salvar arquivo).
    Útil para passar direto pro módulo de screening.
    """
    resp = requests.get(
        f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(name)}/property/IsomericSMILES/TXT",
        timeout=15,
    )
    resp.raise_for_status()
    smiles = resp.text.strip()
    print(f"[fetch] SMILES de '{name}': {smiles}")
    return smiles


# ---------------------------------------------------------------------------
# Banco molecular — ChEMBL (principal) + ZINC20 (opcional)
# ---------------------------------------------------------------------------

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
ZINC20_BASE = "https://zinc20.docking.org"


def fetch_database(
    outdir: str = ".",
    source: str = "chembl",
    mw_min: float = 150,
    mw_max: float = 500,
    logp_min: float = -1,
    logp_max: float = 5,
    max_mols: int = 500,
) -> Path:
    """
    Baixa banco de moléculas pequenas para screening.

    Parâmetros
    ----------
    outdir   : pasta de destino
    source   : "chembl" (padrão, estável) | "zinc" | "local"
    mw_min/max   : intervalo de peso molecular
    logp_min/max : intervalo de logP
    max_mols : número máximo de moléculas

    Retorna
    -------
    Path do arquivo .smi gerado
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if source == "chembl":
        return _fetch_chembl(outdir, mw_min, mw_max, logp_min, logp_max, max_mols)
    elif source == "zinc":
        return _fetch_zinc(outdir, mw_min, mw_max, logp_min, logp_max, max_mols)
    else:
        raise ValueError(f"Fonte desconhecida: {source}. Use 'chembl' ou 'zinc'.")


# Mantém fetch_zinc como alias para compatibilidade
def fetch_zinc(outdir=".", availability="for-sale", mw_min=150, mw_max=500,
               logp_min=-1, logp_max=5, max_mols=500):
    return fetch_database(outdir=outdir, source="zinc", mw_min=mw_min,
                          mw_max=mw_max, logp_min=logp_min, logp_max=logp_max,
                          max_mols=max_mols)


def _fetch_chembl(
    outdir: Path,
    mw_min: float,
    mw_max: float,
    logp_min: float,
    logp_max: float,
    max_mols: int,
) -> Path:
    """
    Baixa moléculas do ChEMBL via API REST com filtros de propriedade.
    Retorna arquivo .smi com formato: SMILES<tab>ID
    """
    print(f"[fetch] Conectando ao ChEMBL...")
    print(f"  MW   : {mw_min}–{mw_max}")
    print(f"  logP : {logp_min}–{logp_max}")
    print(f"  Máx. : {max_mols} moléculas")

    params = {
        "molecule_properties__mw_freebase__lte": mw_max,
        "molecule_properties__mw_freebase__gte": mw_min,
        "molecule_properties__alogp__lte": logp_max,
        "molecule_properties__alogp__gte": logp_min,
        "molecule_type": "Small molecule",
        "limit": min(max_mols, 1000),
        "offset": 0,
    }

    outfile = outdir / "chembl_database.smi"
    all_lines = []

    with tqdm(total=max_mols, desc="Baixando ChEMBL", unit=" mols") as pbar:
        while len(all_lines) < max_mols:
            resp = requests.get(
                f"{CHEMBL_BASE}/molecule.json",
                params=params,
                timeout=30,
            )
            resp.raise_for_status()
            data = resp.json()
            mols = data.get("molecules", [])

            if not mols:
                break

            for mol in mols:
                structs = mol.get("molecule_structures") or {}
                smiles = structs.get("canonical_smiles", "")
                chembl_id = mol.get("molecule_chembl_id", "")
                if smiles:
                    all_lines.append(f"{smiles}\t{chembl_id}")
                    pbar.update(1)
                if len(all_lines) >= max_mols:
                    break

            # Paginação
            next_url = data.get("page_meta", {}).get("next")
            if not next_url or len(all_lines) >= max_mols:
                break

            params["offset"] += params["limit"]
            time.sleep(0.2)

    outfile.write_text("\n".join(all_lines))
    print(f"[fetch] {len(all_lines)} moléculas salvas em: {outfile}")
    return outfile


def _fetch_zinc(
    outdir: Path,
    mw_min: float,
    mw_max: float,
    logp_min: float,
    logp_max: float,
    max_mols: int,
) -> Path:
    """
    Tenta baixar do ZINC20. Se falhar por SSL, orienta o usuário
    a baixar manualmente e usar --source local.
    """
    print(f"[fetch] Conectando ao ZINC20...")
    print(f"  MW   : {mw_min}–{mw_max}")
    print(f"  logP : {logp_min}–{logp_max}")
    print(f"  Máx. : {max_mols} moléculas")

    url = f"{ZINC20_BASE}/substances/subsets/for-sale.smi"
    params = {
        "mw": f"{mw_min}:{mw_max}",
        "logp": f"{logp_min}:{logp_max}",
        "count": min(max_mols, 1000),
        "output_fields": "zinc_id,smiles",
    }

    try:
        resp = requests.get(url, params=params, timeout=60, stream=True, verify=False)
        resp.raise_for_status()
    except requests.exceptions.SSLError:
        print(
            "\n[fetch] ERRO: ZINC20 recusou a conexão SSL.\n"
            "  O ZINC20 tem instabilidade frequente de SSL.\n"
            "  Alternativas:\n"
            "    1) Use --source chembl (recomendado)\n"
            "    2) Baixe manualmente em https://zinc20.docking.org e use --source local\n"
        )
        raise SystemExit(1)

    outfile = outdir / "zinc_database.smi"
    lines = []
    for line in tqdm(resp.iter_lines(), desc="Baixando ZINC", unit=" mols"):
        if line:
            lines.append(line.decode("utf-8"))
        if len(lines) >= max_mols:
            break

    outfile.write_text("\n".join(lines))
    print(f"[fetch] {len(lines)} moléculas salvas em: {outfile}")
    return outfile
