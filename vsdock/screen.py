"""
vsdock.screen
=============
Triagem por similaridade estrutural usando RDKit.
Compara o ligante de referência contra um banco de moléculas (SMILES)
e retorna os hits acima do threshold de Tanimoto definido.
"""

import pandas as pd
from pathlib import Path
from tqdm import tqdm

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys
from rdkit import RDLogger

# Silencia warnings do RDKit sobre moléculas inválidas
RDLogger.DisableLog("rdApp.*")


# ---------------------------------------------------------------------------
# Fingerprints disponíveis
# ---------------------------------------------------------------------------

def _get_fingerprint(mol, fp_type: str = "morgan", radius: int = 2):
    """
    Calcula o fingerprint de uma molécula.

    Tipos suportados
    ----------------
    morgan : Morgan circular (ECFP4 com radius=2)
    maccs  : MACCS 166 keys
    rdkit  : RDKit topological
    """
    if fp_type == "morgan":
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=2048)
    elif fp_type == "maccs":
        return MACCSkeys.GenMACCSKeys(mol)
    elif fp_type == "rdkit":
        return Chem.RDKFingerprint(mol)
    else:
        raise ValueError(f"Fingerprint desconhecido: {fp_type}. Use morgan, maccs ou rdkit.")


# ---------------------------------------------------------------------------
# Função principal
# ---------------------------------------------------------------------------

def screen(
    query_smiles: str,
    database_file: str,
    outdir: str = "hits",
    threshold: float = 0.4,
    fp_type: str = "morgan",
    radius: int = 2,
    max_hits: int = 500,
) -> pd.DataFrame:
    """
    Triagem por similaridade de Tanimoto.

    Parâmetros
    ----------
    query_smiles  : SMILES do ligante de referência
    database_file : arquivo .smi com o banco (SMILES<tab>ID por linha)
    outdir        : pasta para salvar os hits
    threshold     : Tanimoto mínimo (0–1)
    fp_type       : tipo de fingerprint (morgan | maccs | rdkit)
    radius        : raio para Morgan fingerprint
    max_hits      : número máximo de hits a retornar

    Retorna
    -------
    DataFrame com colunas: smiles, id, tanimoto
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # --- Processa o ligante de referência ---
    query_mol = Chem.MolFromSmiles(query_smiles)
    if query_mol is None:
        raise ValueError(f"SMILES inválido para o ligante: {query_smiles}")

    query_fp = _get_fingerprint(query_mol, fp_type, radius)
    print(f"[screen] Ligante de referência processado")
    print(f"  Fingerprint : {fp_type} (radius={radius})")
    print(f"  Threshold   : Tanimoto >= {threshold}")

    # --- Lê o banco ---
    db_path = Path(database_file)
    if not db_path.exists():
        raise FileNotFoundError(f"Banco não encontrado: {database_file}")

    lines = [l.strip() for l in db_path.read_text().splitlines() if l.strip()]
    print(f"[screen] Banco carregado: {len(lines)} moléculas")

    # --- Triagem ---
    hits = []
    invalid = 0

    for line in tqdm(lines, desc="Triagem", unit=" mols"):
        parts = line.split()
        if len(parts) < 1:
            continue

        smiles = parts[0]
        mol_id = parts[1] if len(parts) > 1 else f"mol_{len(hits)+invalid}"

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            invalid += 1
            continue

        fp = _get_fingerprint(mol, fp_type, radius)
        tanimoto = DataStructs.TanimotoSimilarity(query_fp, fp)

        if tanimoto >= threshold:
            hits.append({
                "smiles": smiles,
                "id": mol_id,
                "tanimoto": round(tanimoto, 4),
            })

    if invalid > 0:
        print(f"[screen] {invalid} moléculas ignoradas (SMILES inválido)")

    # --- Ordena por similaridade e aplica limite ---
    if not hits:
        print(f"[screen] Nenhum hit encontrado com threshold={threshold}. Tente reduzir o valor.")
        return pd.DataFrame(columns=["smiles", "id", "tanimoto"])

    df = pd.DataFrame(hits).sort_values("tanimoto", ascending=False)
    df = df.head(max_hits).reset_index(drop=True)

    print(f"[screen] {len(df)} hits encontrados (threshold={threshold})")

    # --- Salva resultados ---
    hits_smi = outdir / "hits.smi"
    hits_csv = outdir / "hits.csv"

    df[["smiles", "id"]].to_csv(hits_smi, sep="\t", index=False, header=False)
    df.to_csv(hits_csv, index=False)

    print(f"[screen] Hits salvos em:")
    print(f"  {hits_smi}")
    print(f"  {hits_csv}")

    # --- Resumo ---
    if len(df) > 0:
        print(f"\n[screen] Top 5 hits:")
        print(df.head(5).to_string(index=False))

    return df


# ---------------------------------------------------------------------------
# Utilitário: carrega SMILES do estado do projeto
# ---------------------------------------------------------------------------

def load_query_from_state(state_file: str = "vsdock_state.yaml") -> str:
    """Lê o SMILES do ligante salvo pelo comando init."""
    import yaml
    state = yaml.safe_load(Path(state_file).read_text())
    smiles = state.get("ligand_smiles")
    if not smiles:
        raise ValueError("SMILES do ligante não encontrado em vsdock_state.yaml. Rode 'vsdock init' primeiro.")
    return smiles
