"""
vsdock.plip
===========
Análise de interações proteína-ligante usando PLIP.

Fluxo por molécula:
  receptor.pdbqt + ligante_out.pdbqt -> PDB combinado -> PLIP -> XML parseado

Requer: plip (no PATH), obabel
"""

import subprocess
import shutil
import xml.etree.ElementTree as ET
import pandas as pd
from pathlib import Path
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Conversão PDBQT -> PDB e combinação com receptor
# ---------------------------------------------------------------------------

def _pdbqt_to_pdb(pdbqt_path: Path, pdb_path: Path):
    """Converte PDBQT para PDB via obabel."""
    cmd = ["obabel", str(pdbqt_path), "-O", str(pdb_path)]
    subprocess.run(cmd, capture_output=True)


def _build_complex_pdb(receptor_pdbqt: str, ligand_pdbqt: Path, out_pdb: Path):
    """
    Combina receptor + ligante numa PDB única para o PLIP.
    Extrai apenas o melhor modo (MODEL 1) do ligante.
    """
    # Converte receptor
    rec_pdb = out_pdb.parent / "receptor_tmp.pdb"
    _pdbqt_to_pdb(Path(receptor_pdbqt), rec_pdb)

    # Extrai modelo 1 do ligante (melhor pose)
    lig_lines = []
    in_model1 = False
    for line in Path(ligand_pdbqt).read_text().splitlines():
        if line.startswith("MODEL 1") or (not line.startswith("MODEL") and not lig_lines and line.startswith("ATOM")):
            in_model1 = True
        if in_model1:
            if line.startswith("ENDMDL") or line.startswith("MODEL 2"):
                break
            if line.startswith(("ATOM", "HETATM", "REMARK")):
                # Marca como HETATM com resname LIG
                if line.startswith("ATOM"):
                    line = "HETATM" + line[6:17] + "LIG" + line[20:]
                lig_lines.append(line)

    # Combina
    rec_content = rec_pdb.read_text()
    # Remove END do receptor antes de combinar
    rec_content = rec_content.replace("\nEND\n", "\n").replace("\nEND", "")

    combined = rec_content + "\n" + "\n".join(lig_lines) + "\nEND\n"
    out_pdb.write_text(combined)
    rec_pdb.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Execução do PLIP e parse do XML
# ---------------------------------------------------------------------------

def _run_plip(complex_pdb: Path, out_dir: Path) -> Path | None:
    """Roda PLIP em modo XML e retorna o Path do XML gerado."""
    cmd = [
        "plip",
        "-f", str(complex_pdb),
        "-o", str(out_dir),
        "-x",   # output XML
        "-q",   # quieto
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    # PLIP gera report.xml na pasta de saída
    xml_candidates = list(out_dir.glob("*.xml")) + list(out_dir.glob("report.xml"))
    if xml_candidates:
        return xml_candidates[0]
    return None


def _parse_plip_xml(xml_path: Path, mol_id: str) -> list:
    """
    Parseia o XML do PLIP e retorna lista de interações.

    Tipos de interação extraídos:
    - hydrophobic_interaction
    - hydrogen_bond
    - water_bridge
    - salt_bridge
    - pi_stack
    - pi_cation_interaction
    - halogen_bond
    """
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
    except ET.ParseError:
        return []

    interactions = []

    # Itera sobre sítios de ligação
    for bs in root.iter("bindingsite"):
        for interaction_type in [
            "hydrophobic_interaction",
            "hydrogen_bond",
            "water_bridge",
            "salt_bridge",
            "pi_stack",
            "pi_cation_interaction",
            "halogen_bond",
        ]:
            for item in bs.iter(interaction_type):
                entry = {
                    "mol_id": mol_id,
                    "interaction_type": interaction_type,
                    "residue": "",
                    "residue_nr": "",
                    "chain": "",
                    "distance": "",
                }
                # Extrai resíduo
                for tag in ["resnr", "restype", "reschain", "dist", "dist_h-a", "dist_a-w"]:
                    el = item.find(tag)
                    if el is not None and el.text:
                        if tag == "resnr":
                            entry["residue_nr"] = el.text
                        elif tag == "restype":
                            entry["residue"] = el.text
                        elif tag == "reschain":
                            entry["chain"] = el.text
                        elif tag in ["dist", "dist_h-a", "dist_a-w"]:
                            entry["distance"] = el.text

                interactions.append(entry)

    return interactions


# ---------------------------------------------------------------------------
# Função principal
# ---------------------------------------------------------------------------

def analyze_plip(
    docking_results: str,
    receptor_pdbqt: str,
    poses_dir: str = "docking/poses",
    outdir: str = "plip",
    top_n: int = 10,
) -> pd.DataFrame:
    """
    Roda PLIP nas melhores poses do docking.

    Parâmetros
    ----------
    docking_results : CSV de resultados do docking (saída do dock_all)
    receptor_pdbqt  : arquivo .pdbqt do receptor
    poses_dir       : pasta com os arquivos _out.pdbqt gerados pelo Vina
    outdir          : pasta de saída
    top_n           : número de melhores hits para analisar

    Retorna
    -------
    DataFrame com todas as interações encontradas
    """
    if not shutil.which("plip"):
        raise EnvironmentError("plip não encontrado no PATH.")

    outdir    = Path(outdir)
    poses_dir = Path(poses_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Carrega resultados do docking
    df_dock = pd.read_csv(docking_results).head(top_n)
    print(f"[plip] Analisando top {len(df_dock)} poses")

    all_interactions = []

    for _, row in tqdm(df_dock.iterrows(), total=len(df_dock), desc="PLIP", unit=" mol"):
        mol_id  = str(row["id"]).replace("/", "_").replace(" ", "_")
        pose_file = poses_dir / f"{mol_id}_out.pdbqt"

        if not pose_file.exists():
            print(f"[plip] Pose não encontrada para {mol_id}, pulando")
            continue

        mol_outdir  = outdir / mol_id
        mol_outdir.mkdir(exist_ok=True)
        complex_pdb = mol_outdir / f"{mol_id}_complex.pdb"

        # Monta PDB do complexo
        _build_complex_pdb(receptor_pdbqt, pose_file, complex_pdb)

        # Roda PLIP
        xml_path = _run_plip(complex_pdb, mol_outdir)
        if xml_path is None:
            print(f"[plip] PLIP não gerou XML para {mol_id}")
            continue

        # Parseia interações
        interactions = _parse_plip_xml(xml_path, mol_id)
        all_interactions.extend(interactions)

    if not all_interactions:
        print("[plip] Nenhuma interação encontrada.")
        return pd.DataFrame()

    df = pd.DataFrame(all_interactions)

    # Salva resultados
    out_csv = outdir / "plip_interactions.csv"
    df.to_csv(out_csv, index=False)

    # Resumo por molécula
    summary = (
        df.groupby(["mol_id", "interaction_type"])
        .size()
        .reset_index(name="count")
        .pivot_table(index="mol_id", columns="interaction_type", values="count", fill_value=0)
        .reset_index()
    )
    summary.to_csv(outdir / "plip_summary.csv", index=False)

    print(f"[plip] {len(df)} interações encontradas em {df['mol_id'].nunique()} moléculas")
    print(f"[plip] Resultados salvos em:")
    print(f"  {out_csv}")
    print(f"  {outdir}/plip_summary.csv")
    print(f"\n[plip] Resumo:")
    print(summary.to_string(index=False))

    return df
