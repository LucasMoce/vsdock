"""
vsdock.plif
===========
Protein-Ligand Interaction Fingerprint (PLIF) usando ProLIF.

O PLIF codifica as interações proteína-ligante como um vetor binário,
permitindo comparação quantitativa entre poses de docking.

Recursos:
  - Gera fingerprint de interações para cada pose
  - Calcula similaridade de Tanimoto entre fingerprints
  - Compara hits com o ligante de referência
  - Gera heatmap de interações
  - Exporta matriz de fingerprints para análise downstream

Requer: prolif, MDAnalysis, obabel
"""

import subprocess
import shutil
import warnings
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm

warnings.filterwarnings("ignore", category=DeprecationWarning)


# ---------------------------------------------------------------------------
# Conversão PDBQT -> PDB (necessário para ProLIF)
# ---------------------------------------------------------------------------

def _pdbqt_to_pdb(pdbqt_path: Path, pdb_path: Path):
    """Converte PDBQT para PDB via obabel."""
    cmd = ["obabel", str(pdbqt_path), "-O", str(pdb_path), "-h"]
    subprocess.run(cmd, capture_output=True)


def _extract_best_pose_pdb(pose_pdbqt: Path, out_pdb: Path):
    """
    Extrai o modelo 1 do PDBQT do Vina e gera PDB com elementos corretos.
    Converte tipos AutoDock (A, OA, NA, etc.) para elementos PDB padrão.
    """
    AD_TO_ELEMENT = {
        "C": "C",  "A": "C",  "N": "N",  "NA": "N", "NS": "N",
        "O": "O",  "OA": "O", "OS": "O", "S": "S",  "SA": "S",
        "H": "H",  "HD": "H", "HS": "H", "P": "P",  "F": "F",
        "Cl": "Cl","CL": "Cl","Br": "Br","BR": "Br","I": "I",
        "Mg": "Mg","Ca": "Ca","Mn": "Mn","Fe": "Fe","Zn": "Zn",
    }

    lines = []
    in_model = False

    for line in pose_pdbqt.read_text().splitlines():
        if line.startswith("MODEL") and not in_model:
            in_model = True
            continue
        if line.startswith("ENDMDL") and in_model:
            break
        if in_model and line.startswith(("ATOM", "HETATM")):
            parts = line.split()
            ad_type = parts[-1] if parts else "C"
            element = AD_TO_ELEMENT.get(ad_type, "C")
            # Reconstrói linha PDB com elemento correto na coluna 77-78
            new_line = "HETATM" + line[6:17] + "LIG" + line[20:76].ljust(56) + f"{element:>2}"
            lines.append(new_line)

    # Fallback: sem MODEL tags
    if not lines:
        for line in pose_pdbqt.read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                ad_type = parts[-1] if parts else "C"
                element = AD_TO_ELEMENT.get(ad_type, "C")
                new_line = "HETATM" + line[6:17] + "LIG" + line[20:76].ljust(56) + f"{element:>2}"
                lines.append(new_line)

    if lines:
        out_pdb.write_text("\n".join(lines) + "\nEND\n")


def _compute_plif(receptor_pdb: Path, ligand_pdb: Path) -> dict:
    """
    Calcula o PLIF para um par receptor-ligante usando ProLIF 2.x.
    Adiciona hidrogênios explícitos ao ligante antes do cálculo.
    """
    try:
        import prolif as plf
        import MDAnalysis as mda
        warnings.filterwarnings("ignore")

        # Adiciona H ao ligante via obabel antes de carregar
        ligand_pdb_h = ligand_pdb.parent / f"{ligand_pdb.stem}_H.pdb"
        cmd = ["obabel", str(ligand_pdb), "-O", str(ligand_pdb_h), "-h"]
        subprocess.run(cmd, capture_output=True)

        lig_file = ligand_pdb_h if ligand_pdb_h.exists() and ligand_pdb_h.stat().st_size > 10 \
                   else ligand_pdb

        # Carrega estruturas
        u_prot = mda.Universe(str(receptor_pdb))
        u_lig  = mda.Universe(str(lig_file))

        protein = u_prot.select_atoms("protein")
        ligand  = u_lig.select_atoms("all")

        if len(protein) == 0 or len(ligand) == 0:
            return {}

        # Converte para ProLIF Molecule — usa inferrer=None para não exigir H perfeitos
        try:
            lig_mol  = plf.Molecule.from_mda(ligand, inferrer=None)
            prot_mol = plf.Molecule.from_mda(protein, inferrer=None)
        except TypeError:
            # versão mais antiga do ProLIF sem inferrer
            try:
                lig_mol  = plf.Molecule.from_mda(ligand)
                prot_mol = plf.Molecule.from_mda(protein)
            except Exception:
                return {}
        except Exception:
            return {}

        # Calcula fingerprint
        fp = plf.Fingerprint()
        fp.run_from_iterable([lig_mol], prot_mol)
        df = fp.to_dataframe()

        if df is None or df.empty:
            return {}

        # Converte MultiIndex para dicionário plano
        result = {}
        for col in df.columns:
            try:
                if isinstance(col, tuple):
                    key = "_".join(str(c) for c in col[1:])
                else:
                    key = str(col)
                result[key] = bool(df[col].iloc[0])
            except Exception:
                pass

        # Limpa arquivo temporário
        ligand_pdb_h.unlink(missing_ok=True)

        return result

    except Exception:
        return {}


# ---------------------------------------------------------------------------
# Função principal
# ---------------------------------------------------------------------------

def compute_plif(
    docking_results: str,
    receptor_pdbqt: str,
    poses_dir: str = "docking/poses",
    outdir: str = "plif",
    top_n: int = 10,
    heatmap: bool = False,
) -> pd.DataFrame:
    """
    Calcula Protein-Ligand Interaction Fingerprints para as melhores poses.
    O ligante de referência é sempre incluído.

    Parâmetros
    ----------
    docking_results : CSV de resultados do docking
    receptor_pdbqt  : arquivo .pdbqt do receptor
    poses_dir       : pasta com os _out.pdbqt gerados pelo Vina
    outdir          : pasta de saída
    top_n           : número de melhores hits para analisar

    Retorna
    -------
    DataFrame com matriz de fingerprints (linhas=moléculas, colunas=interações)
    """
    try:
        import prolif as plf
    except ImportError:
        raise EnvironmentError(
            "ProLIF não encontrado. Instale com: pip install prolif"
        )

    if not shutil.which("obabel"):
        raise EnvironmentError("obabel não encontrado no PATH.")

    outdir    = Path(outdir)
    poses_dir = Path(poses_dir)
    pdb_dir   = outdir / "pdb_tmp"
    outdir.mkdir(parents=True, exist_ok=True)
    pdb_dir.mkdir(exist_ok=True)

    # Prepara receptor em PDB — usa receptor_clean.pdb se disponível
    receptor_pdb = pdb_dir / "receptor.pdb"
    receptor_clean = Path(".") / "receptor_clean.pdb"
    if receptor_clean.exists():
        import shutil as sh2
        sh2.copy(receptor_clean, receptor_pdb)
        print(f"[plif] Usando receptor_clean.pdb")
    else:
        _pdbqt_to_pdb(Path(receptor_pdbqt), receptor_pdb)
        print(f"[plif] Convertendo receptor.pdbqt para PDB")

    # Carrega resultados do docking
    df_dock = pd.read_csv(docking_results)

    # Separa referência dos demais
    if "is_reference" in df_dock.columns:
        ref_rows = df_dock[df_dock["is_reference"] == True]
        non_ref  = df_dock[df_dock["is_reference"] != True]
    else:
        ref_rows = pd.DataFrame()
        non_ref  = df_dock

    # Top N + referência
    top_hits   = non_ref.head(top_n)
    to_analyze = pd.concat([top_hits, ref_rows]).reset_index(drop=True) \
                 if not ref_rows.empty else top_hits

    print(f"[plif] Calculando PLIF para {len(top_hits)} hits", end="")
    if not ref_rows.empty:
        print(f" + ligante de referência ({ref_rows.iloc[0]['id']})")
    else:
        print()

    all_fps = []
    mol_ids = []
    is_refs = []

    for _, row in tqdm(to_analyze.iterrows(), total=len(to_analyze),
                       desc="PLIF", unit=" mol"):
        mol_id  = str(row["id"]).replace("/", "_").replace(" ", "_")
        is_ref  = bool(row.get("is_reference", False))
        pose_f  = poses_dir / f"{mol_id}_out.pdbqt"

        if not pose_f.exists():
            if is_ref:
                print(f"[plif] AVISO: pose do ligante de referência '{mol_id}' não encontrada.")
            else:
                print(f"[plif] Pose não encontrada para {mol_id}, pulando")
            continue

        # Converte pose para PDB
        lig_pdb = pdb_dir / f"{mol_id}_lig.pdb"
        _extract_best_pose_pdb(pose_f, lig_pdb)

        if not lig_pdb.exists() or lig_pdb.stat().st_size == 0:
            continue

        # Calcula PLIF
        fp = _compute_plif(receptor_pdb, lig_pdb)

        if fp:
            all_fps.append(fp)
            mol_ids.append(row["id"])
            is_refs.append(is_ref)

    if not all_fps:
        print("[plif] Nenhum fingerprint gerado.")
        return pd.DataFrame()

    # Monta matriz de fingerprints
    df_fp = pd.DataFrame(all_fps, index=mol_ids).fillna(False)
    df_fp.index.name = "mol_id"
    df_fp = df_fp.reset_index()
    df_fp.insert(1, "is_reference", is_refs)

    # Calcula similaridade de Tanimoto com a referência
    ref_mask = [r for r in is_refs if r]
    if ref_mask and sum(is_refs) > 0:
        ref_idx   = is_refs.index(True)
        ref_fp    = df_fp.iloc[ref_idx, 2:].astype(int).values

        similarities = []
        for i in range(len(df_fp)):
            fp_i = df_fp.iloc[i, 2:].astype(int).values
            intersection = np.logical_and(ref_fp, fp_i).sum()
            union        = np.logical_or(ref_fp, fp_i).sum()
            tanimoto     = float(intersection / union) if union > 0 else 0.0
            similarities.append(round(tanimoto, 4))

        df_fp.insert(2, "tanimoto_vs_ref", similarities)
        df_fp = df_fp.sort_values("tanimoto_vs_ref", ascending=False).reset_index(drop=True)

    # Marca referência com ★ no display
    ref_name = ref_rows.iloc[0]["id"] if not ref_rows.empty else None

    # Salva matriz completa
    out_matrix = outdir / "plif_matrix.csv"
    df_fp.to_csv(out_matrix, index=False)

    # Salva resumo — número de interações por molécula
    interaction_cols = [c for c in df_fp.columns
                        if c not in ["mol_id", "is_reference", "tanimoto_vs_ref"]]
    summary = df_fp[["mol_id", "is_reference"]].copy()
    if "tanimoto_vs_ref" in df_fp.columns:
        summary["tanimoto_vs_ref"] = df_fp["tanimoto_vs_ref"]
    summary["n_interactions"] = df_fp[interaction_cols].sum(axis=1).astype(int)
    summary["mol_id_display"] = summary.apply(
        lambda r: f"★ {r['mol_id']} (ref)" if r["is_reference"] else r["mol_id"], axis=1
    )
    summary.to_csv(outdir / "plif_summary.csv", index=False)

    print(f"\n[plif] {len(df_fp)} fingerprints gerados")
    print(f"  Interações únicas : {len(interaction_cols)}")
    print(f"\n[plif] Resultados salvos em:")
    print(f"  {out_matrix}           (matriz completa)")
    print(f"  {outdir}/plif_summary.csv     (resumo por molécula)")

    # Resumo no terminal
    display_cols = ["mol_id_display", "n_interactions"]
    if "tanimoto_vs_ref" in summary.columns:
        display_cols.append("tanimoto_vs_ref")
    print(f"\n[plif] Resumo (ordenado por similaridade com referência):")
    print(summary[display_cols].to_string(index=False))

    # Gera heatmap apenas se solicitado explicitamente
    if heatmap:
        _try_heatmap(df_fp, interaction_cols, ref_name, outdir)
    # (matplotlib pode causar segfault em ambientes headless sem --heatmap)

    # Limpa PDBs temporários
    import shutil as sh
    sh.rmtree(pdb_dir, ignore_errors=True)

    return df_fp


def _try_heatmap(df_fp: pd.DataFrame, interaction_cols: list,
                 ref_name: str, outdir: Path):
    """Gera heatmap de interações se matplotlib disponível."""
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        if len(interaction_cols) == 0 or len(df_fp) == 0:
            return

        # Limita a 30 interações mais frequentes para legibilidade
        top_cols = df_fp[interaction_cols].sum().nlargest(30).index.tolist()
        matrix   = df_fp[top_cols].astype(int)

        labels = []
        for _, row in df_fp.iterrows():
            label = f"★ {row['mol_id']} (ref)" if row["is_reference"] else str(row["mol_id"])
            labels.append(label)

        fig, ax = plt.subplots(figsize=(max(12, len(top_cols)*0.4),
                                        max(6, len(labels)*0.4)))

        im = ax.imshow(matrix.values, cmap="Blues", aspect="auto", vmin=0, vmax=1)

        ax.set_xticks(range(len(top_cols)))
        ax.set_xticklabels(top_cols, rotation=90, fontsize=7)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)

        ax.set_title("Protein-Ligand Interaction Fingerprint (PLIF)", fontsize=12, pad=15)
        ax.set_xlabel("Interação (resíduo_tipo)", fontsize=9)
        ax.set_ylabel("Composto", fontsize=9)

        plt.colorbar(im, ax=ax, label="Interação presente")
        plt.tight_layout()

        out_fig = outdir / "plif_heatmap.png"
        plt.savefig(out_fig, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  {out_fig}      (heatmap)")

    except Exception as e:
        print(f"[plif] Heatmap não gerado: {e}")
