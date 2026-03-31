"""
vsdock.dock
===========
Docking molecular com AutoDock Vina 1.2.x.

Fluxo por molécula:
  SMILES -> SDF (RDKit, 3D) -> PDBQT (Open Babel) -> Vina -> score parseado

Requer: obabel, vina (no PATH)
"""

import os
import re
import subprocess
import shutil
import pandas as pd
from pathlib import Path
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


# ---------------------------------------------------------------------------
# Checagem de dependências
# ---------------------------------------------------------------------------

def _check_deps():
    for tool in ["vina", "obabel"]:
        if not shutil.which(tool):
            raise EnvironmentError(
                f"'{tool}' não encontrado no PATH. "
                f"Instale antes de usar o módulo dock."
            )


# ---------------------------------------------------------------------------
# Autobox — geração automática de caixa de docking
# ---------------------------------------------------------------------------

def get_box_from_autobox(
    pdb_file: str,
    ligand: str = None,
    residues: list = None,
    padding: float = 5.0,
    blind: bool = False,
) -> dict:
    """
    Usa o autobox para calcular automaticamente as coordenadas da caixa de docking.

    Parâmetros
    ----------
    pdb_file  : arquivo PDB do receptor (com ligante cocristalizado)
    ligand    : identificador do ligante no PDB (ex: "MRV1101A")
    residues  : lista de resíduos (ex: ["VAL2A", "LYS4A"])
    padding   : padding em Angstroms ao redor do ligante/resíduos
    blind     : se True, faz docking cego (caixa = proteína toda)

    Retorna
    -------
    dict com center (x,y,z) e size (x,y,z)
    """
    if not shutil.which("autobox"):
        raise EnvironmentError(
            "autobox não encontrado. Instale com:\n"
            "pip install git+https://github.com/omixlab/autobox"
        )

    cmd = ["autobox", "--input", pdb_file, "--padding", str(int(padding))]

    if blind:
        cmd.append("--blind")
    elif ligand:
        cmd += ["--ligands", ligand]
    elif residues:
        cmd += ["--residues"] + residues
    else:
        raise ValueError("Informe --ligand, --residues ou --blind para o autobox.")

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"autobox falhou:\n{result.stderr}")

    # Parseia output: "x_center=149.85\ny_center=108.41\n..."
    box = {}
    for line in result.stdout.splitlines():
        if "=" in line:
            key, val = line.strip().split("=")
            box[key.strip()] = float(val.strip())

    center = (box["x_center"], box["y_center"], box["z_center"])
    size   = (box["x_size"],   box["y_size"],   box["z_size"])

    print(f"[dock] Caixa calculada pelo autobox:")
    print(f"  Centro : {center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}")
    print(f"  Tamanho: {size[0]:.2f}, {size[1]:.2f}, {size[2]:.2f} Å")

    return {"center": center, "size": size}


# ---------------------------------------------------------------------------
# Conversão SMILES -> PDBQT
# ---------------------------------------------------------------------------

def smiles_to_pdbqt(smiles: str, mol_id: str, outdir: Path) -> Path | None:
    """
    Converte um SMILES para PDBQT com coordenadas 3D.
    Retorna o Path do arquivo .pdbqt ou None se falhar.
    """
    sdf_path  = outdir / f"{mol_id}.sdf"
    pdbqt_path = outdir / f"{mol_id}.pdbqt"

    # 1. SMILES -> SDF com geometria 3D via RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result != 0:
        # Fallback: tenta ETKDG sem constraints
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if result != 0:
        return None

    AllChem.MMFFOptimizeMolecule(mol)

    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol)
    writer.close()

    # 2. SDF -> PDBQT via Open Babel
    cmd = [
        "obabel", str(sdf_path),
        "-O", str(pdbqt_path),
        "--partialcharge", "gasteiger",
        "-h",  # adiciona H se necessário
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    sdf_path.unlink(missing_ok=True)  # limpa SDF temporário

    if not pdbqt_path.exists():
        return None

    return pdbqt_path


# ---------------------------------------------------------------------------
# Execução do Vina
# ---------------------------------------------------------------------------

def run_vina(
    ligand_pdbqt: Path,
    receptor_pdbqt: str,
    center: tuple,
    size: tuple,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    energy_range: int = 3,
    out_dir: Path = None,
) -> float | None:
    """
    Roda o AutoDock Vina para um ligante e retorna o melhor score (kcal/mol).

    Parâmetros
    ----------
    ligand_pdbqt  : Path do ligante em PDBQT
    receptor_pdbqt: Path do receptor em PDBQT
    center        : (x, y, z) do centro da caixa de docking
    size          : (sx, sy, sz) tamanho da caixa em Angstroms
    exhaustiveness: exaustividade da busca (default 8)
    num_modes     : número de poses a gerar
    energy_range  : range de energia para poses alternativas

    Retorna
    -------
    Melhor score em kcal/mol ou None se falhar
    """
    out_pdbqt = out_dir / f"{ligand_pdbqt.stem}_out.pdbqt"

    cmd = [
        "vina",
        "--receptor", str(receptor_pdbqt),
        "--ligand",   str(ligand_pdbqt),
        "--out",      str(out_pdbqt),
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x",   str(size[0]),
        "--size_y",   str(size[1]),
        "--size_z",   str(size[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes",      str(num_modes),
        "--energy_range",   str(energy_range),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        return None

    # Parseia o score da saída do Vina
    score = _parse_vina_score(result.stdout)
    return score


def _parse_vina_score(vina_output: str) -> float | None:
    """Extrai o melhor score (modo 1) da saída do Vina."""
    for line in vina_output.splitlines():
        # Linha do melhor modo: "   1      -7.5      0.000      0.000"
        match = re.match(r"\s+1\s+([-\d.]+)", line)
        if match:
            return float(match.group(1))
    return None


# ---------------------------------------------------------------------------
# Função principal
# ---------------------------------------------------------------------------

def dock_all(
    hits_file: str,
    receptor: str,
    center: tuple,
    size: tuple = (20, 20, 20),
    outdir: str = "docking",
    exhaustiveness: int = 8,
    num_modes: int = 9,
    top_n: int = None,
    reference_smiles: str = None,
    reference_name: str = "reference",
) -> pd.DataFrame:
    """
    Roda docking de todos os hits + ligante de referência.

    Parâmetros
    ----------
    hits_file        : CSV de hits limpos (saída do pains)
    receptor         : Path do receptor .pdbqt
    center           : (x, y, z) centro da caixa
    size             : (sx, sy, sz) tamanho da caixa
    outdir           : pasta de saída
    exhaustiveness   : exaustividade Vina
    num_modes        : poses por ligante
    top_n            : limita número de hits (None = todos)
    reference_smiles : SMILES do ligante de referência (docka junto)
    reference_name   : nome do ligante de referência

    Retorna
    -------
    DataFrame ranqueado por score com colunas: id, smiles, score_kcal_mol
    """
    _check_deps()

    outdir = Path(outdir)
    pdbqt_dir = outdir / "pdbqt"
    poses_dir = outdir / "poses"
    pdbqt_dir.mkdir(parents=True, exist_ok=True)
    poses_dir.mkdir(parents=True, exist_ok=True)

    if not Path(receptor).exists():
        raise FileNotFoundError(f"Receptor não encontrado: {receptor}")

    # Carrega hits
    df_hits = pd.read_csv(hits_file)
    if top_n:
        df_hits = df_hits.head(top_n)

    print(f"[dock] {len(df_hits)} hits para dockar")
    print(f"  Receptor       : {receptor}")
    print(f"  Centro         : {center}")
    print(f"  Caixa          : {size} Å")
    print(f"  Exhaustiveness : {exhaustiveness}")

    # Adiciona ligante de referência à lista
    rows = []
    if reference_smiles:
        rows.append({"id": reference_name, "smiles": reference_smiles, "tanimoto": 1.0})

    for _, row in df_hits.iterrows():
        rows.append({"id": row["id"], "smiles": row["smiles"], "tanimoto": row.get("tanimoto", None)})

    results = []
    failed = []

    for entry in tqdm(rows, desc="Docking", unit=" mol"):
        mol_id  = str(entry["id"]).replace("/", "_").replace(" ", "_")
        smiles  = entry["smiles"]
        tanimoto = entry.get("tanimoto")

        # Converte para PDBQT
        pdbqt_path = smiles_to_pdbqt(smiles, mol_id, pdbqt_dir)
        if pdbqt_path is None:
            failed.append(mol_id)
            continue

        # Roda Vina
        score = run_vina(
            ligand_pdbqt=pdbqt_path,
            receptor_pdbqt=receptor,
            center=center,
            size=size,
            exhaustiveness=exhaustiveness,
            num_modes=num_modes,
            out_dir=poses_dir,
        )

        if score is None:
            failed.append(mol_id)
            continue

        results.append({
            "id": entry["id"],
            "smiles": smiles,
            "tanimoto": tanimoto,
            "score_kcal_mol": score,
        })

    if failed:
        print(f"[dock] {len(failed)} moléculas falharam na conversão/docking")

    df = pd.DataFrame(results).sort_values("score_kcal_mol").reset_index(drop=True)

    # Marca hits melhores que o ligante de referência
    if reference_smiles and len(df) > 0:
        ref_row = df[df["id"] == reference_name]
        if not ref_row.empty:
            ref_score = ref_row.iloc[0]["score_kcal_mol"]
            df["better_than_ref"] = df["score_kcal_mol"] < ref_score
            print(f"\n[dock] Score do ligante de referência : {ref_score} kcal/mol")
            n_better = df[df["id"] != reference_name]["better_than_ref"].sum()
            print(f"[dock] Hits melhores que referência   : {n_better}")

    # Salva resultados
    out_csv = outdir / "docking_results.csv"
    df.to_csv(out_csv, index=False)
    print(f"\n[dock] Resultados salvos em: {out_csv}")
    print(f"\n[dock] Top 10:")
    print(df.head(10)[["id", "score_kcal_mol", "tanimoto"]].to_string(index=False))

    return df
