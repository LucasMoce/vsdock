"""
vsdock.report
=============
Geração de relatório/manuscrito consolidando todos os resultados do pipeline.

Gera um arquivo Markdown estruturado como manuscrito científico com:
- Seção de métodos (preenchida automaticamente com parâmetros usados)
- Tabelas de resultados (docking, ADMET, PLIP)
- Espaços marcados para o usuário completar (introdução, discussão)
"""

import pandas as pd
from pathlib import Path
from datetime import date


# ---------------------------------------------------------------------------
# Carregamento dos resultados
# ---------------------------------------------------------------------------

def _load_results(project_dir: Path) -> dict:
    """Carrega todos os CSVs gerados pelo pipeline."""
    results = {}

    files = {
        "state":    project_dir / "vsdock_state.yaml",
        "docking":  project_dir / "docking" / "docking_results.csv",
        "admet":    project_dir / "admet" / "admet_summary.csv",
        "plip":     project_dir / "plip" / "plip_summary.csv",
        "hits":     project_dir / "hits" / "hits_clean.csv",
        "pains":    project_dir / "hits" / "hits_pains_report.csv",
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


# ---------------------------------------------------------------------------
# Seções do manuscrito
# ---------------------------------------------------------------------------

def _section_title(results: dict) -> str:
    state = results.get("state") or {}
    ligand = state.get("ligand_name", "ligand").capitalize()
    return f"# Virtual Screening and Molecular Docking Study of {ligand} Analogues\n"


def _section_abstract(results: dict) -> str:
    state   = results.get("state") or {}
    ligand  = state.get("ligand_name", "the reference ligand")
    df_dock = results.get("docking")
    df_admet = results.get("admet")

    n_hits  = len(df_dock) - 1 if df_dock is not None else "N/A"  # exclui referência
    n_admet = len(df_admet) if df_admet is not None else "N/A"

    ref_score = "N/A"
    if df_dock is not None:
        ref_row = df_dock[df_dock["id"] == ligand]
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
    state = results.get("state") or {}
    ligand   = state.get("ligand_name", "N/A")
    receptor = state.get("receptor", "N/A")
    center   = state.get("box_center", ["N/A", "N/A", "N/A"])
    size     = state.get("box_size",   ["N/A", "N/A", "N/A"])

    df_hits = results.get("hits")
    n_screened = len(df_hits) if df_hits is not None else "N/A"

    df_pains = results.get("pains")
    n_pains = int(df_pains["is_pains"].sum()) if df_pains is not None and "is_pains" in df_pains.columns else "N/A"

    df_dock = results.get("docking")
    n_docked = len(df_dock) - 1 if df_dock is not None else "N/A"

    return f"""## 2. Materials and Methods

### 2.1 Reference Ligand and Database

The reference ligand **{ligand}** was retrieved from PubChem. Molecular structures
were obtained from the ChEMBL database and filtered by drug-like properties
(molecular weight 150–500 Da; logP −1 to 5).

### 2.2 Similarity-Based Virtual Screening

Structural similarity screening was performed using Morgan circular fingerprints
(radius = 2, 2048 bits) and Tanimoto coefficient as the similarity metric.
Compounds with Tanimoto similarity ≥ 0.4 to the reference ligand were retained.
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
(x = {center[0]:.2f}, y = {center[1]:.2f}, z = {center[2]:.2f}) Å
with dimensions {size[0]:.1f} × {size[1]:.1f} × {size[2]:.1f} Å.
Exhaustiveness was set to 8. A total of **{n_docked} compounds** were docked,
including the reference ligand.
Three-dimensional molecular conformations were generated using RDKit (ETKDG method)
and converted to PDBQT format using Open Babel 3.1.0.

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
    state   = results.get("state") or {}
    ligand  = state.get("ligand_name", "reference")
    df_dock = results.get("docking")
    df_admet = results.get("admet")
    df_plip  = results.get("plip")

    sections = ["## 3. Results\n"]

    # --- Docking ---
    sections.append("### 3.1 Molecular Docking\n")
    if df_dock is not None:
        ref_row = df_dock[df_dock["id"] == ligand]
        ref_score = ref_row.iloc[0]["score_kcal_mol"] if not ref_row.empty else None

        if ref_score is not None:
            sections.append(
                f"The reference ligand **{ligand}** achieved a docking score of "
                f"**{ref_score:.2f} kcal/mol**, which was used as the benchmark threshold.\n"
            )

        better = df_dock[df_dock["id"] != ligand]
        if ref_score is not None:
            better = better[better["score_kcal_mol"] < ref_score]

        sections.append(
            f"A total of **{len(better)} compounds** achieved docking scores "
            f"better than the reference ligand.\n"
        )

        # Tabela top 10
        top10 = df_dock[df_dock["id"] != ligand].head(10)[
            ["id", "score_kcal_mol", "tanimoto"]
        ].copy()
        top10.columns = ["Compound ID", "Docking Score (kcal/mol)", "Tanimoto Similarity"]
        sections.append("\n**Table 1.** Top 10 compounds by docking score.\n")
        sections.append(top10.to_markdown(index=False))
        sections.append("\n")
    else:
        sections.append("> *Docking results not found.*\n")

    # --- PLIP ---
    sections.append("\n### 3.2 Protein–Ligand Interactions\n")
    if df_plip is not None:
        sections.append(
            "Key protein–ligand interactions for the top-ranked compounds "
            "are summarized in Table 2.\n"
        )
        sections.append("\n**Table 2.** Protein–ligand interactions (PLIP analysis).\n")
        sections.append(df_plip.to_markdown(index=False))
        sections.append("\n")
        sections.append(
            "> **[TO COMPLETE]** Describe the key interactions observed for the "
            "best candidates compared to the reference ligand.\n"
        )
    else:
        sections.append("> *PLIP results not found.*\n")

    # --- ADMET ---
    sections.append("\n### 3.3 ADMET Profiling\n")
    if df_admet is not None:
        admet_cols = ["id", "admet_score"] + [
            c for c in df_admet.columns
            if any(ep in c for ep in ["hERG", "DILI", "AMES", "ClinTox", "Bioavail", "Solubility"])
            and "percentile" not in c
        ]
        admet_cols = [c for c in admet_cols if c in df_admet.columns]
        df_admet_show = df_admet[admet_cols].head(10).copy()

        sections.append(
            "ADMET properties were predicted for the top-ranked compounds. "
            "Results are shown in Table 3. Lower ADMET score indicates a "
            "more favorable overall safety and pharmacokinetic profile.\n"
        )
        sections.append("\n**Table 3.** ADMET predictions for top compounds.\n")
        sections.append(df_admet_show.to_markdown(index=False))
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


# ---------------------------------------------------------------------------
# Função principal
# ---------------------------------------------------------------------------

def generate_report(
    outdir: str = "report",
    project_dir: str = ".",
    fmt: str = "markdown",
) -> Path:
    """
    Gera o relatório/manuscrito consolidando todos os resultados.

    Parâmetros
    ----------
    outdir      : pasta de saída
    project_dir : pasta raiz do projeto vsdock
    fmt         : formato de saída ("markdown" | "quarto" | "html")

    Retorna
    -------
    Path do arquivo gerado
    """
    outdir      = Path(outdir)
    project_dir = Path(project_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("[report] Carregando resultados do pipeline...")
    results = _load_results(project_dir)

    found = [k for k, v in results.items() if v is not None and k != "state"]
    print(f"[report] Módulos encontrados: {', '.join(found)}")

    # Monta o manuscrito
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

    # Salva
    if fmt == "quarto":
        # Adiciona header YAML para Quarto
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
            print("[report] markdown não instalado, salvando como .md")
            fmt = "markdown"
        outfile = outdir / "manuscript.html"
    else:
        outfile = outdir / "manuscript.md"

    outfile.write_text(manuscript)
    print(f"[report] Manuscrito gerado em: {outfile}")
    print(f"\n[report] Seções marcadas com [TO COMPLETE] precisam de revisão manual.")

    return outfile
