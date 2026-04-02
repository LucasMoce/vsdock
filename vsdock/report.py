"""
vsdock.report
=============
Geração de relatório/manuscrito consolidando todos os resultados do pipeline.
O ligante de referência é sempre incluído nas tabelas, independente do ranking.
"""

import pandas as pd
from pathlib import Path
from datetime import date


def _load_results(project_dir: Path) -> dict:
    results = {}
    files = {
        "state":   project_dir / "vsdock_state.yaml",
        "docking": project_dir / "docking" / "docking_results.csv",
        "admet":   project_dir / "admet" / "admet_summary.csv",
        "plip":    project_dir / "plip" / "plip_summary.csv",
        "hits":    project_dir / "hits" / "hits_clean.csv",
        "pains":   project_dir / "hits" / "hits_pains_report.csv",
    }
    for key, path in files.items():
        if path.exists():
            if path.suffix == ".yaml":
                import yaml
                results[key] = yaml.safe_load(path.read_text())
            else:
                results[key] = pd.read_csv(path)
        else:
            results[key] = None
    return results


def _get_ref_row(df: pd.DataFrame, ligand: str) -> pd.DataFrame:
    """Retorna a linha do ligante de referência, usando is_reference ou id."""
    if "is_reference" in df.columns:
        ref = df[df["is_reference"] == True]
        if not ref.empty:
            return ref
    return df[df["id"].astype(str) == str(ligand)]


def _label_ref(val, ligand: str, is_ref_col=None) -> str:
    """Adiciona ★ ao nome se for o ligante de referência."""
    if is_ref_col is not None and is_ref_col:
        return f"★ {val} (ref)"
    if str(val) == str(ligand):
        return f"★ {val} (ref)"
    return val


def _top_with_ref(df: pd.DataFrame, ligand: str, n: int = 10,
                  sort_col: str = "score_kcal_mol", ascending: bool = True) -> pd.DataFrame:
    """
    Retorna os top N compostos garantindo que o ligante de referência está incluído.
    Se o ligante não está no top N, é adicionado ao final.
    """
    ref_rows = _get_ref_row(df, ligand)

    # Top N excluindo referência
    non_ref = df[~df.index.isin(ref_rows.index)]
    top_non_ref = non_ref.sort_values(sort_col, ascending=ascending).head(n)

    # Combina: top N + referência (se não estiver já incluída)
    ref_in_top = not ref_rows.empty and any(ref_rows.index.isin(top_non_ref.index))
    if not ref_rows.empty and not ref_in_top:
        result = pd.concat([top_non_ref, ref_rows]).reset_index(drop=True)
    else:
        result = top_non_ref.reset_index(drop=True)

    return result


def _section_title(results: dict) -> str:
    state  = results.get("state") or {}
    ligand = state.get("ligand_name")
    if ligand:
        return f"# Virtual Screening and Molecular Docking Study of {ligand.capitalize()} Analogues\n"
    return "# Virtual Screening and Molecular Docking Study\n"


def _section_abstract(results: dict) -> str:
    state    = results.get("state") or {}
    ligand   = state.get("ligand_name", "the reference ligand")
    df_dock  = results.get("docking")
    df_admet = results.get("admet")

    n_hits  = "N/A"
    n_admet = len(df_admet) if df_admet is not None else "N/A"
    ref_score = "N/A"

    if df_dock is not None:
        non_ref = df_dock[df_dock["id"].astype(str) != str(ligand)] if "is_reference" not in df_dock.columns \
                  else df_dock[df_dock["is_reference"] != True]
        n_hits = len(non_ref)
        ref_row = _get_ref_row(df_dock, ligand)
        if not ref_row.empty:
            ref_score = f"{ref_row.iloc[0]['score_kcal_mol']:.2f} kcal/mol"

    return f"""## Abstract

> **[TO COMPLETE]** Summarize the biological context, objective, key findings, and conclusions.

*Placeholder:* In this study, we performed a structure-based virtual screening campaign
to identify novel analogues of {ligand}. A total of {n_hits} compounds were screened
by molecular docking against the target receptor. The reference ligand ({ligand})
achieved a docking score of {ref_score}. ADMET profiling was performed for {n_admet}
top-ranked compounds. The most promising candidates were further characterized by
protein-ligand interaction analysis.

---
"""


def _section_introduction() -> str:
    return """## 1. Introduction

> **[TO COMPLETE]** Describe the biological target, disease context, and rationale
> for using the reference ligand as a scaffold. Include relevant references.

---
"""


def _section_methods(results: dict) -> str:
    state    = results.get("state") or {}
    ligand   = state.get("ligand_name", "N/A")
    receptor = state.get("receptor", "N/A")
    center   = state.get("box_center", ["N/A", "N/A", "N/A"])
    size     = state.get("box_size",   ["N/A", "N/A", "N/A"])

    df_hits  = results.get("hits")
    n_screened = len(df_hits) if df_hits is not None else "N/A"

    df_pains = results.get("pains")
    n_pains  = int(df_pains["is_pains"].sum()) \
               if df_pains is not None and "is_pains" in df_pains.columns else "N/A"

    df_dock  = results.get("docking")
    n_docked = "N/A"
    if df_dock is not None:
        non_ref = df_dock[df_dock["id"].astype(str) != str(ligand)] if "is_reference" not in df_dock.columns \
                  else df_dock[df_dock["is_reference"] != True]
        n_docked = len(non_ref)

    threshold = state.get("screen_threshold", 0.4)

    cx = f"{center[0]:.2f}" if isinstance(center[0], float) else center[0]
    cy = f"{center[1]:.2f}" if isinstance(center[1], float) else center[1]
    cz = f"{center[2]:.2f}" if isinstance(center[2], float) else center[2]
    sx = f"{size[0]:.1f}" if isinstance(size[0], float) else size[0]
    sy = f"{size[1]:.1f}" if isinstance(size[1], float) else size[1]
    sz = f"{size[2]:.1f}" if isinstance(size[2], float) else size[2]

    return f"""## 2. Materials and Methods

### 2.1 Reference Ligand and Database

The reference ligand **{ligand}** was retrieved from PubChem. Molecular structures
were obtained from the ChEMBL database and filtered by drug-like properties
(molecular weight 150–500 Da; logP −1 to 5).

### 2.2 Similarity-Based Virtual Screening

Structural similarity screening was performed using Morgan circular fingerprints
(radius = 2, 2048 bits) and Tanimoto coefficient as the similarity metric.
Compounds with Tanimoto similarity ≥ {threshold} to the reference ligand were retained.
A total of **{n_screened} hits** were identified.

### 2.3 PAINS and Drug-likeness Filtering

Pan-Assay Interference Compounds (PAINS) were identified and removed using the
RDKit FilterCatalog (PAINS-A, -B, -C subsets) (Baell & Holloway, 2010).
Compounds violating more than one Lipinski Rule of Five criterion were also excluded.
**{n_pains} PAINS** alerts were detected.

### 2.4 Molecular Docking

Molecular docking was performed using AutoDock Vina v1.2.3 (Eberhardt et al., 2021).
The receptor structure **{receptor}** was used as the docking target.
The docking grid box was centered at coordinates
(x = {cx}, y = {cy}, z = {cz}) Å
with dimensions {sx} × {sy} × {sz} Å.
Exhaustiveness was set to 8. A total of **{n_docked} compounds** were docked,
alongside the reference ligand **{ligand}** as a benchmark.
Three-dimensional molecular conformations were generated using RDKit (ETKDG method)
and converted to PDBQT format using Open Babel 3.1.0.
Prior to 3D embedding, salt forms were standardized by retaining only the largest
organic fragment.

### 2.5 Protein–Ligand Interaction Analysis

Protein–ligand interactions for the top-ranked compounds were analyzed using
PLIP (Protein–Ligand Interaction Profiler) (Salentin et al., 2015).
Interaction types analyzed included hydrogen bonds, hydrophobic contacts,
π-stacking, π-cation interactions, salt bridges, and halogen bonds.

### 2.6 ADMET Prediction

ADMET (Absorption, Distribution, Metabolism, Excretion, and Toxicity) properties
were predicted using ADMET-AI (Swanson et al., 2024), a machine learning framework
trained on the Therapeutics Data Commons benchmark dataset.
Key endpoints evaluated included hERG inhibition, DILI (Drug-Induced Liver Injury),
AMES mutagenicity, clinical toxicity (ClinTox), oral bioavailability, and aqueous solubility.

---
"""


def _section_results(results: dict) -> str:
    state    = results.get("state") or {}
    ligand   = state.get("ligand_name", "reference")
    df_dock  = results.get("docking")
    df_admet = results.get("admet")
    df_plip  = results.get("plip")

    sections = ["## 3. Results\n"]

    # -----------------------------------------------------------------------
    # 3.1 Docking
    # -----------------------------------------------------------------------
    sections.append("### 3.1 Molecular Docking\n")
    if df_dock is not None:
        ref_row   = _get_ref_row(df_dock, ligand)
        ref_score = ref_row.iloc[0]["score_kcal_mol"] if not ref_row.empty else None

        if ref_score is not None:
            sections.append(
                f"The reference ligand **{ligand}** achieved a docking score of "
                f"**{ref_score:.2f} kcal/mol**, which was used as the benchmark threshold.\n"
            )

        non_ref = df_dock[df_dock["id"].astype(str) != str(ligand)] if "is_reference" not in df_dock.columns \
                  else df_dock[df_dock["is_reference"] != True]

        if ref_score is not None:
            n_better = (non_ref["score_kcal_mol"] < ref_score).sum()
            sections.append(
                f"A total of **{n_better} compounds** achieved docking scores "
                f"better than the reference ligand.\n"
            )

        # Top 10 + referência garantida
        top = _top_with_ref(df_dock, ligand, n=10, sort_col="score_kcal_mol", ascending=True)

        is_ref_col = top.get("is_reference") if "is_reference" in top.columns else None
        top["Compound ID"] = top.apply(
            lambda r: f"★ {r['id']} (ref)" if (
                ("is_reference" in r and r["is_reference"]) or str(r["id"]) == str(ligand)
            ) else r["id"], axis=1
        )
        top = top.rename(columns={
            "score_kcal_mol": "Docking Score (kcal/mol)",
            "tanimoto": "Tanimoto Similarity",
        })
        cols = ["Compound ID", "Docking Score (kcal/mol)", "Tanimoto Similarity"]
        cols = [c for c in cols if c in top.columns]

        sections.append("\n**Table 1.** Top 10 compounds by docking score (★ = reference ligand).\n")
        sections.append(top[cols].to_markdown(index=False))
        sections.append("\n")
    else:
        sections.append("> *Docking results not found.*\n")

    # -----------------------------------------------------------------------
    # 3.2 PLIP
    # -----------------------------------------------------------------------
    sections.append("\n### 3.2 Protein–Ligand Interactions\n")
    if df_plip is not None:
        df_plip_display = df_plip.copy()

        # Marca referência — usa coluna is_reference se existir, senão compara por nome
        if "is_reference" in df_plip_display.columns:
            df_plip_display["mol_id"] = df_plip_display.apply(
                lambda r: f"★ {r['mol_id']} (ref)" if r["is_reference"] else r["mol_id"], axis=1
            )
            df_plip_display = df_plip_display.drop(columns=["is_reference"])
        else:
            df_plip_display["mol_id"] = df_plip_display["mol_id"].apply(
                lambda x: f"★ {x} (ref)" if str(x) == str(ligand) else x
            )

        ref_in_plip = any(str(ligand) in str(m) for m in df_plip["mol_id"].tolist())

        sections.append(
            "Key protein–ligand interactions for the top-ranked compounds "
            "are summarized in Table 2 (★ = reference ligand).\n"
        )
        if not ref_in_plip:
            sections.append(
                f"> **Note:** The reference ligand ({ligand}) was not included in the "
                f"PLIP analysis because its docking pose was not found. "
                f"Re-run `vsdock plip` after a successful docking run.\n"
            )

        sections.append("\n**Table 2.** Protein–ligand interactions (PLIP). ★ = reference ligand.\n")
        sections.append(df_plip_display.to_markdown(index=False))
        sections.append("\n")
        sections.append(
            "> **[TO COMPLETE]** Describe the key interactions observed for the "
            "best candidates compared to the reference ligand.\n"
        )
    else:
        sections.append("> *PLIP results not found.*\n")

    # -----------------------------------------------------------------------
    # 3.3 ADMET
    # -----------------------------------------------------------------------
    sections.append("\n### 3.3 ADMET Profiling\n")
    if df_admet is not None:
        admet_cols = ["id", "admet_score"] + [
            c for c in df_admet.columns
            if any(ep in c for ep in ["hERG", "DILI", "AMES", "ClinTox", "Bioavail", "Solubility"])
            and "percentile" not in c
        ]
        admet_cols = [c for c in admet_cols if c in df_admet.columns]

        # Top 10 por admet_score + referência garantida
        top_admet = _top_with_ref(df_admet, ligand, n=10,
                                   sort_col="admet_score", ascending=True)
        top_admet = top_admet[admet_cols].copy()
        top_admet["id"] = top_admet["id"].apply(
            lambda x: f"★ {x} (ref)" if str(x) == str(ligand) else x
        )

        sections.append(
            "ADMET properties for the top-ranked compounds are shown in Table 3. "
            "Lower ADMET score indicates a more favorable safety/PK profile "
            "(★ = reference ligand).\n"
        )
        sections.append("\n**Table 3.** ADMET predictions. ★ = reference ligand.\n")
        sections.append(top_admet.to_markdown(index=False))
        sections.append("\n")
        sections.append(
            "> **[TO COMPLETE]** Discuss ADMET findings, highlighting compounds "
            "with favorable profiles and flagging any toxicity concerns.\n"
        )
    else:
        sections.append("> *ADMET results not found.*\n")

    sections.append("\n---\n")
    return "\n".join(sections)


def _section_discussion() -> str:
    return """## 4. Discussion

> **[TO COMPLETE]** Integrate the docking, interaction, and ADMET findings.
> Identify the most promising candidates and justify selection based on
> score, interactions, and safety profile. Compare with literature if available.

---
"""


def _section_conclusion() -> str:
    return """## 5. Conclusion

> **[TO COMPLETE]** Summarize the key findings and propose next steps
> (e.g., MD simulations, in vitro validation).

---
"""


def _section_references() -> str:
    return """## References

- Baell, J.B. & Holloway, G.A. (2010). New Substructure Filters for Removal of Pan Assay Interference Compounds (PAINS). *J. Med. Chem.*, 53(7), 2719–2740.
- Eberhardt, J. et al. (2021). AutoDock Vina 1.2.0. *J. Chem. Inf. Model.*, 61(8), 3891–3898.
- Salentin, S. et al. (2015). PLIP: fully automated protein–ligand interaction profiler. *Nucleic Acids Res.*, 43(W1), W443–W447.
- Swanson, K. et al. (2024). ADMET-AI: A machine learning ADMET platform. *Bioinformatics*, 40(1), btad693.
- Trott, O. & Olson, A.J. (2010). AutoDock Vina. *J. Comput. Chem.*, 31(2), 455–461.

> **[TO COMPLETE]** Add target-specific references and any additional citations.
"""


def generate_report(
    outdir: str = "report",
    project_dir: str = ".",
    fmt: str = "markdown",
) -> Path:
    outdir      = Path(outdir)
    project_dir = Path(project_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("[report] Carregando resultados do pipeline...")
    results = _load_results(project_dir)

    found = [k for k, v in results.items() if v is not None and k != "state"]
    print(f"[report] Módulos encontrados: {', '.join(found)}")

    sections = [
        _section_title(results),
        f"*Generated by vsdock on {date.today().isoformat()}*\n\n---\n",
        _section_abstract(results),
        _section_introduction(),
        _section_methods(results),
        _section_results(results),
        _section_discussion(),
        _section_conclusion(),
        _section_references(),
    ]

    manuscript = "\n".join(sections)

    if fmt == "quarto":
        header = """---
title: "Virtual Screening Report"
format:
  html:
    toc: true
    theme: cosmo
  pdf:
    toc: true
execute:
  echo: false
---

"""
        manuscript = header + manuscript
        outfile = outdir / "manuscript.qmd"
    elif fmt == "html":
        try:
            import markdown
            html_body = markdown.markdown(manuscript, extensions=["tables"])
            manuscript = f"<html><body>{html_body}</body></html>"
        except ImportError:
            fmt = "markdown"
        outfile = outdir / "manuscript.html"
    else:
        outfile = outdir / "manuscript.md"

    outfile.write_text(manuscript)
    print(f"[report] Manuscrito gerado em: {outfile}")
    print(f"\n[report] Seções marcadas com [TO COMPLETE] precisam de revisão manual.")

    return outfile
