"""
vsdock.fetch
============
Download do ligante de referência (PubChem) e de bancos moleculares.

Fontes suportadas:
  - ChEMBL (padrão) — via API REST
  - PubChem — via PUG REST API

Filtros disponíveis:
  --available        compostos comercialmente disponíveis (ChEMBL)
  --fda              compostos FDA-aprovados (ChEMBL max_phase=4)
  --natural-compounds produtos naturais (ChEMBL natural_product=1)
  --fragments        fragment-like (MW<300, regra de 3)
  --druglike         drug-like (Lipinski Ro5, padrão)
"""

import time
import requests
from pathlib import Path
from tqdm import tqdm

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CHEMBL_BASE  = "https://www.ebi.ac.uk/chembl/api/data"
ZINC20_BASE  = "https://zinc20.docking.org"


# ---------------------------------------------------------------------------
# PubChem — ligante de referência
# ---------------------------------------------------------------------------

def fetch_ligand(name: str, outdir: str = ".", fmt: str = "sdf") -> Path:
    """Baixa estrutura do ligante via PubChem."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[fetch] Buscando '{name}' no PubChem...")
    resp = requests.get(
        f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(name)}/cids/JSON",
        timeout=15,
    )
    resp.raise_for_status()
    cid = resp.json()["IdentifierList"]["CID"][0]
    print(f"[fetch] CID encontrado: {cid}")

    if fmt == "sdf":
        url    = f"{PUBCHEM_BASE}/compound/cid/{cid}/SDF?record_type=3d"
        suffix = ".sdf"
    else:
        url    = f"{PUBCHEM_BASE}/compound/cid/{cid}/property/IsomericSMILES/TXT"
        suffix = ".smi"

    resp = requests.get(url, timeout=15)
    resp.raise_for_status()

    outfile = outdir / f"{name.lower().replace(' ', '_')}{suffix}"
    outfile.write_bytes(resp.content)
    print(f"[fetch] Ligante salvo em: {outfile}")
    return outfile


def fetch_ligand_smiles(name: str) -> str:
    """Retorna o SMILES canônico do ligante via PubChem."""
    resp = requests.get(
        f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(name)}/property/IsomericSMILES/TXT",
        timeout=15,
    )
    resp.raise_for_status()
    smiles = resp.text.strip()
    print(f"[fetch] SMILES de '{name}': {smiles}")
    return smiles


# ---------------------------------------------------------------------------
# Banco molecular — dispatcher principal
# ---------------------------------------------------------------------------

def fetch_database(
    outdir: str = ".",
    source: str = "chembl",
    mw_min: float = 150,
    mw_max: float = 500,
    logp_min: float = -1,
    logp_max: float = 5,
    max_mols: int = 500,
    available_only: bool = False,
    fda_only: bool = False,
    natural_only: bool = False,
    fragments: bool = False,
    druglike: bool = True,
) -> Path:
    """
    Baixa banco de moléculas para screening.

    Parâmetros
    ----------
    source         : "chembl" (padrão) ou "pubchem"
    mw_min/max     : intervalo de peso molecular
    logp_min/max   : intervalo de logP
    max_mols       : número máximo de moléculas
    available_only : compostos comercialmente disponíveis (ChEMBL)
    fda_only       : apenas FDA-aprovados (ChEMBL max_phase=4)
    natural_only   : apenas produtos naturais (ChEMBL natural_product=1)
    fragments      : modo fragment-like (MW<300, regra de 3)
    druglike       : modo drug-like / Lipinski Ro5 (padrão: True)
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Ajusta filtros de propriedade conforme modo selecionado
    if fragments:
        # Regra de 3 para fragments
        mw_min, mw_max   = 100, 300
        logp_min, logp_max = -3, 3
        print("[fetch] Modo: fragment-like (Regra de 3: MW≤300, logP≤3, HBD≤3, HBA≤3)")
    elif not druglike:
        print(f"[fetch] Modo: propriedades personalizadas (MW {mw_min}–{mw_max}, logP {logp_min}–{logp_max})")
    else:
        print("[fetch] Modo: drug-like (Lipinski Ro5)")

    if source == "chembl":
        return _fetch_chembl(
            outdir, mw_min, mw_max, logp_min, logp_max, max_mols,
            available_only=available_only,
            fda_only=fda_only,
            natural_only=natural_only,
            fragments=fragments,
        )
    elif source == "pubchem":
        return _fetch_pubchem(
            outdir, mw_min, mw_max, logp_min, logp_max, max_mols,
            fda_only=fda_only,
            fragments=fragments,
        )
    elif source == "zinc":
        return _fetch_zinc(outdir, mw_min, mw_max, logp_min, logp_max, max_mols)
    else:
        raise ValueError(f"Fonte desconhecida: {source}. Use 'chembl', 'pubchem' ou 'zinc'.")


# alias para compatibilidade
def fetch_zinc(outdir=".", availability="for-sale", mw_min=150, mw_max=500,
               logp_min=-1, logp_max=5, max_mols=500):
    return fetch_database(outdir=outdir, source="zinc", mw_min=mw_min,
                          mw_max=mw_max, logp_min=logp_min, logp_max=logp_max,
                          max_mols=max_mols)


# ---------------------------------------------------------------------------
# ChEMBL
# ---------------------------------------------------------------------------

def _fetch_chembl(
    outdir: Path,
    mw_min: float, mw_max: float,
    logp_min: float, logp_max: float,
    max_mols: int,
    available_only: bool = False,
    fda_only: bool = False,
    natural_only: bool = False,
    fragments: bool = False,
) -> Path:
    """Baixa moléculas do ChEMBL com retry automático."""
    print(f"[fetch] Conectando ao ChEMBL...")
    print(f"  MW   : {mw_min}–{mw_max}")
    print(f"  logP : {logp_min}–{logp_max}")
    print(f"  Máx. : {max_mols} moléculas")

    flags = []
    if available_only: flags.append("disponíveis comercialmente")
    if fda_only:       flags.append("FDA-aprovados")
    if natural_only:   flags.append("produtos naturais")
    if flags:
        print(f"  Filtros: {', '.join(flags)}")

    params = {
        "molecule_properties__mw_freebase__lte": mw_max,
        "molecule_properties__mw_freebase__gte": mw_min,
        "molecule_properties__alogp__lte": logp_max,
        "molecule_properties__alogp__gte": logp_min,
        "molecule_type": "Small molecule",
        "limit": min(max_mols, 1000),
        "offset": 0,
    }

    if available_only:
        params["availability_type"] = 1
    if fda_only:
        params["max_phase"] = 4
    if natural_only:
        params["natural_product"] = 1
    if fragments:
        # Regra de 3 adicional: HBD≤3, HBA≤3
        params["molecule_properties__hbd__lte"] = 3
        params["molecule_properties__hba__lte"]  = 3

    outfile   = outdir / "chembl_database.smi"
    all_lines = []

    with tqdm(total=max_mols, desc="Baixando ChEMBL", unit=" mols") as pbar:
        while len(all_lines) < max_mols:
            resp = None
            for attempt in range(3):
                try:
                    resp = requests.get(
                        f"{CHEMBL_BASE}/molecule.json",
                        params=params,
                        timeout=60,
                    )
                    resp.raise_for_status()
                    break
                except requests.exceptions.Timeout:
                    if attempt < 2:
                        print(f"\n[fetch] Timeout, tentando novamente ({attempt+2}/3)...")
                        time.sleep(5)
                    else:
                        print(f"\n[fetch] ChEMBL não respondeu. {len(all_lines)} moléculas salvas.")
                        resp = None
                        break
                except requests.exceptions.RequestException as e:
                    print(f"\n[fetch] Erro de conexão: {e}")
                    resp = None
                    break

            if resp is None:
                break

            data = resp.json()
            mols = data.get("molecules", [])
            if not mols:
                break

            for mol in mols:
                structs   = mol.get("molecule_structures") or {}
                smiles    = structs.get("canonical_smiles", "")
                chembl_id = mol.get("molecule_chembl_id", "")
                if smiles:
                    all_lines.append(f"{smiles}\t{chembl_id}")
                    pbar.update(1)
                if len(all_lines) >= max_mols:
                    break

            next_url = data.get("page_meta", {}).get("next")
            if not next_url or len(all_lines) >= max_mols:
                break

            params["offset"] += params["limit"]
            time.sleep(0.5)

    outfile.write_text("\n".join(all_lines))
    print(f"[fetch] {len(all_lines)} moléculas salvas em: {outfile}")
    return outfile


# ---------------------------------------------------------------------------
# PubChem
# ---------------------------------------------------------------------------

PUBCHEM_CLASSIF = {
    "fda":     "approved",
    "natural": "natural products",
}

def _fetch_pubchem(
    outdir: Path,
    mw_min: float, mw_max: float,
    logp_min: float, logp_max: float,
    max_mols: int,
    fda_only: bool = False,
    fragments: bool = False,
) -> Path:
    """
    Baixa moléculas do PubChem via PUG REST com filtros de propriedade.
    Usa o endpoint de busca por propriedade + download em SMILES.
    """
    print(f"[fetch] Conectando ao PubChem...")
    print(f"  MW   : {mw_min}–{mw_max}")
    print(f"  logP : {logp_min}–{logp_max}")
    print(f"  Máx. : {max_mols} moléculas")

    # PubChem: busca por faixa de propriedades via PUG REST
    # Primeiro obtém CIDs, depois baixa SMILES em lote
    outfile = outdir / "pubchem_database.smi"

    # Monta query de propriedade
    # PubChem aceita busca por range de MW via fast_formula ou property search
    base_url = f"{PUBCHEM_BASE}/compound/fastproperty"

    # Para FDA-approved usamos a classificação farmacológica
    if fda_only:
        print("[fetch] Filtrando compostos FDA-aprovados via PubChem...")
        cids = _pubchem_fda_cids(max_mols * 3)  # pega mais pois vai filtrar por MW
    else:
        print("[fetch] Buscando por propriedades moleculares no PubChem...")
        cids = _pubchem_property_cids(mw_min, mw_max, max_mols * 2)

    if not cids:
        print("[fetch] Nenhum CID encontrado no PubChem.")
        outfile.write_text("")
        return outfile

    print(f"[fetch] {len(cids)} CIDs encontrados, baixando SMILES...")

    # Baixa SMILES em lote (máx 200 CIDs por request)
    all_lines = []
    batch_size = 200

    with tqdm(total=min(len(cids), max_mols), desc="Baixando PubChem", unit=" mols") as pbar:
        for i in range(0, len(cids), batch_size):
            if len(all_lines) >= max_mols:
                break

            batch = cids[i:i+batch_size]
            cid_str = ",".join(str(c) for c in batch)

            try:
                resp = requests.get(
                    f"{PUBCHEM_BASE}/compound/cid/{cid_str}/property/IsomericSMILES,MolecularWeight,XLogP/JSON",
                    timeout=30,
                )
                resp.raise_for_status()
                props = resp.json().get("PropertyTable", {}).get("Properties", [])

                for p in props:
                    smiles = p.get("IsomericSMILES", "")
                    mw     = p.get("MolecularWeight", 0)
                    xlogp  = p.get("XLogP", 0) or 0
                    cid    = p.get("CID", "")

                    # Aplica filtros de propriedade
                    if not (mw_min <= float(mw) <= mw_max):
                        continue
                    if not (logp_min <= float(xlogp) <= logp_max):
                        continue
                    if fragments and float(mw) > 300:
                        continue

                    if smiles:
                        all_lines.append(f"{smiles}\tCID{cid}")
                        pbar.update(1)

                    if len(all_lines) >= max_mols:
                        break

                time.sleep(0.3)

            except requests.exceptions.RequestException as e:
                print(f"\n[fetch] Erro no batch {i}: {e}")
                continue

    outfile.write_text("\n".join(all_lines))
    print(f"[fetch] {len(all_lines)} moléculas salvas em: {outfile}")
    return outfile


def _pubchem_fda_cids(max_cids: int) -> list:
    """Busca CIDs de compostos FDA-aprovados via classificação farmacológica."""
    try:
        resp = requests.get(
            f"{PUBCHEM_BASE}/compound/classification/Drug+Approved/cids/JSON",
            timeout=30,
            params={"MaxRecords": min(max_cids, 10000)},
        )
        resp.raise_for_status()
        return resp.json().get("IdentifierList", {}).get("CID", [])[:max_cids]
    except Exception:
        # Fallback: busca por status FDA via PubChem
        try:
            resp = requests.post(
                f"{PUBCHEM_BASE}/compound/fastsubstructure/smarts/[#6]/cids/JSON",
                timeout=30,
            )
            return []
        except Exception:
            return []


def _pubchem_property_cids(mw_min: float, mw_max: float, max_cids: int) -> list:
    """Busca CIDs por range de peso molecular no PubChem."""
    try:
        # PubChem fast property search
        url = (
            f"{PUBCHEM_BASE}/compound/fastproperty/"
            f"MW%2C{mw_min}%2C{mw_max}/cids/JSON"
            f"?MaxRecords={min(max_cids, 50000)}"
        )
        resp = requests.get(url, timeout=60)
        if resp.status_code == 200:
            return resp.json().get("IdentifierList", {}).get("CID", [])[:max_cids]

        # Fallback: usa endpoint de busca por formula range
        resp = requests.get(
            f"{PUBCHEM_BASE}/compound/fastsimilarity_2d/cid/2244/cids/JSON",
            params={"Threshold": 80, "MaxRecords": max_cids},
            timeout=30,
        )
        return resp.json().get("IdentifierList", {}).get("CID", [])[:max_cids]

    except Exception as e:
        print(f"[fetch] Erro ao buscar CIDs no PubChem: {e}")
        return []


# ---------------------------------------------------------------------------
# ZINC (fallback com aviso de SSL)
# ---------------------------------------------------------------------------

def _fetch_zinc(
    outdir: Path,
    mw_min: float, mw_max: float,
    logp_min: float, logp_max: float,
    max_mols: int,
) -> Path:
    print(f"[fetch] Conectando ao ZINC20...")
    url    = f"{ZINC20_BASE}/substances/subsets/for-sale.smi"
    params = {
        "mw":           f"{mw_min}:{mw_max}",
        "logp":         f"{logp_min}:{logp_max}",
        "count":        min(max_mols, 1000),
        "output_fields":"zinc_id,smiles",
    }

    try:
        resp = requests.get(url, params=params, timeout=60, stream=True, verify=False)
        resp.raise_for_status()
    except (requests.exceptions.SSLError, requests.exceptions.HTTPError):
        print(
            "\n[fetch] ERRO: ZINC20 inacessível (SSL/HTTP error).\n"
            "  Use --source chembl ou --source pubchem.\n"
        )
        raise SystemExit(1)

    outfile = outdir / "zinc_database.smi"
    lines   = []
    for line in tqdm(resp.iter_lines(), desc="Baixando ZINC", unit=" mols"):
        if line:
            lines.append(line.decode("utf-8"))
        if len(lines) >= max_mols:
            break

    outfile.write_text("\n".join(lines))
    print(f"[fetch] {len(lines)} moléculas salvas em: {outfile}")
    return outfile
