"""
vsdock.admet
============
Predição de propriedades ADMET usando ADMET-AI (offline, via CLI).

Novos recursos:
  list_admet_filters()   — lista endpoints disponíveis para filtragem
  apply_admet_filters()  — filtra compostos por reference_compound ou thresholds

Referência: Swanson et al., Bioinformatics 2024
"""

import subprocess
import shutil
import pandas as pd
from pathlib import Path


# ---------------------------------------------------------------------------
# Endpoints disponíveis por categoria
# ---------------------------------------------------------------------------

ADMET_ENDPOINTS = {
    "absorption": {
        "Caco2_Wang":               "Permeabilidade Caco-2 (log cm/s). Bom: > -5.15",
        "HIA_Hou":                  "Absorção intestinal humana (0-1). Bom: > 0.3",
        "Pgp_Broccatelli":          "Substrato P-glicoproteína (0-1). Bom: < 0.5",
        "Bioavailability_Ma":       "Biodisponibilidade oral (0-1). Bom: > 0.3",
        "Lipophilicity_AstraZeneca":"Lipofilicidade logD7.4. Bom: -1 a 3",
        "Solubility_AqSolDB":       "Solubilidade aquosa (log mol/L). Bom: > -4",
    },
    "distribution": {
        "BBB_Martini":              "Barreira hematoencefálica (0-1). Bom: depende do alvo",
        "PPBR_AZ":                  "Ligação proteínas plasmáticas (%). Bom: < 90%",
        "VDss_Lombardo":            "Volume de distribuição (L/kg). Referência: 0.6-7",
    },
    "metabolism": {
        "CYP1A2_Veith":             "Inibidor CYP1A2 (0-1). Bom: < 0.5",
        "CYP2C19_Veith":            "Inibidor CYP2C19 (0-1). Bom: < 0.5",
        "CYP2C9_Veith":             "Inibidor CYP2C9 (0-1). Bom: < 0.5",
        "CYP2D6_Veith":             "Inibidor CYP2D6 (0-1). Bom: < 0.5",
        "CYP3A4_Veith":             "Inibidor CYP3A4 (0-1). Bom: < 0.5",
        "CYP2C9_Substrate_CarbonMangels": "Substrato CYP2C9 (0-1)",
        "CYP2D6_Substrate_CarbonMangels": "Substrato CYP2D6 (0-1)",
        "CYP3A4_Substrate_CarbonMangels": "Substrato CYP3A4 (0-1)",
    },
    "excretion": {
        "Half_Life_Obach":          "Meia-vida (horas). Referência: >1h",
        "Clearance_Hepatocyte_AZ":  "Clearance hepatócito (mL/min/10^6 cells)",
        "Clearance_Microsome_AZ":   "Clearance microssoma (mL/min/g)",
    },
    "toxicity": {
        "hERG":                     "Inibição hERG/cardiotoxicidade (0-1). Bom: < 0.5",
        "hERG_Karim":               "Inibição hERG alternativa (0-1). Bom: < 0.5",
        "AMES":                     "Mutagenicidade AMES (0-1). Bom: < 0.5",
        "DILI":                     "Lesão hepática induzida por droga (0-1). Bom: < 0.5",
        "Hepatotoxicity_Xu":        "Hepatotoxicidade (0-1). Bom: < 0.5",
        "Skin_Reaction":            "Reação cutânea (0-1). Bom: < 0.5",
        "Carcinogens_Lagunin":      "Carcinogenicidade (0-1). Bom: < 0.5",
        "ClinTox":                  "Toxicidade clínica (0-1). Bom: < 0.5",
        "LD50_Zhu":                 "LD50 aguda (log mol/kg). Maior = mais seguro",
    },
}

# Direção de cada endpoint: "lower_better" = menor é melhor (ex: toxicidade)
#                           "higher_better" = maior é melhor (ex: biodisponibilidade)
ENDPOINT_DIRECTION = {
    # Toxicidade — menor probabilidade = melhor
    "hERG": "lower_better", "hERG_Karim": "lower_better",
    "AMES": "lower_better", "DILI": "lower_better",
    "Hepatotoxicity_Xu": "lower_better", "ClinTox": "lower_better",
    "Skin_Reaction": "lower_better", "Carcinogens_Lagunin": "lower_better",
    "Pgp_Broccatelli": "lower_better",
    "CYP1A2_Veith": "lower_better", "CYP2C19_Veith": "lower_better",
    "CYP2C9_Veith": "lower_better", "CYP2D6_Veith": "lower_better",
    "CYP3A4_Veith": "lower_better",
    # Absorção/distribuição — maior é melhor
    "Bioavailability_Ma": "higher_better", "HIA_Hou": "higher_better",
    "Solubility_AqSolDB": "higher_better",
}

DEFAULT_WEIGHTS = {
    "hERG":              0.9,
    "DILI":              0.8,
    "Hepatotoxicity_Xu": 0.8,
    "AMES":              0.7,
    "ClinTox":           0.8,
    "Bioavailability_Ma":0.7,
    "Solubility_AqSolDB":0.6,
    "HIA_Hou":           0.6,
    "BBB_Martini":       0.3,
    "CYP3A4_Veith":      0.5,
}

TOX_KEYWORDS  = ["hERG", "AMES", "DILI", "Hepatotox", "ClinTox",
                  "Skin", "Carcino", "NR-", "SR-"]
GOOD_KEYWORDS = ["Bioavailability", "HIA", "Solubility"]


# ---------------------------------------------------------------------------
# list_admet_filters
# ---------------------------------------------------------------------------

def list_admet_filters():
    """Lista todos os endpoints ADMET disponíveis para filtragem."""
    print("\n[admet] Endpoints disponíveis para --apply-filters:\n")
    for category, endpoints in ADMET_ENDPOINTS.items():
        print(f"  {'─'*50}")
        print(f"  {category.upper()}")
        print(f"  {'─'*50}")
        for name, description in endpoints.items():
            print(f"    {name:<40} {description}")
    print()
    print("  Uso:")
    print("    vsdock admet --apply-filters hERG DILI AMES --filter-by reference_compound")
    print("    vsdock admet --apply-filters hERG DILI --filter-by thresholds --thresholds 'hERG<0.5,DILI<0.3'")
    print()


# ---------------------------------------------------------------------------
# apply_admet_filters
# ---------------------------------------------------------------------------

def apply_admet_filters(
    df: pd.DataFrame,
    filter_by: str,
    filters: list,
    outdir: str = "admet",
    reference_name: str = None,
    thresholds_str: str = None,
) -> pd.DataFrame:
    """
    Filtra compostos com base em endpoints ADMET selecionados.

    Parâmetros
    ----------
    df             : DataFrame com predições ADMET (saída de predict_admet)
    filter_by      : "reference_compound" ou "thresholds"
    filters        : lista de endpoints a considerar (ex: ["hERG", "DILI"])
    outdir         : pasta de saída
    reference_name : nome do ligante de referência (para reference_compound)
    thresholds_str : string de thresholds (ex: "hERG<0.5,DILI<0.3")

    Retorna
    -------
    DataFrame filtrado
    """
    outdir = Path(outdir)

    # Resolve colunas reais no DataFrame para cada endpoint solicitado
    resolved = {}
    for ep in filters:
        col = next((c for c in df.columns if ep in c and "percentile" not in c), None)
        if col:
            resolved[ep] = col
        else:
            print(f"[admet] AVISO: endpoint '{ep}' não encontrado nas predições. Ignorando.")

    if not resolved:
        print("[admet] ERRO: nenhum endpoint válido encontrado.")
        return df

    print(f"\n[admet] Aplicando filtros: {list(resolved.keys())}")
    print(f"  Critério: {filter_by}")

    if filter_by == "reference_compound":
        return _filter_by_reference(df, resolved, reference_name, outdir)
    elif filter_by == "thresholds":
        return _filter_by_thresholds(df, resolved, thresholds_str, outdir)
    else:
        raise ValueError(f"filter_by deve ser 'reference_compound' ou 'thresholds'.")


def _filter_by_reference(
    df: pd.DataFrame,
    resolved: dict,
    reference_name: str,
    outdir: Path,
) -> pd.DataFrame:
    """
    Mantém compostos com perfil ADMET melhor que o ligante de referência
    em todos os endpoints selecionados.
    """
    # Encontra linha da referência
    if "is_reference" in df.columns:
        ref_row = df[df["is_reference"] == True]
    elif reference_name:
        ref_row = df[df["id"].astype(str) == str(reference_name)]
    else:
        print("[admet] ERRO: ligante de referência não encontrado.")
        return df

    if ref_row.empty:
        print(f"[admet] AVISO: referência '{reference_name}' não encontrada no ADMET.")
        return df

    ref_values = {}
    for ep, col in resolved.items():
        ref_values[ep] = ref_row.iloc[0][col]
        direction = ENDPOINT_DIRECTION.get(ep, "lower_better")
        print(f"  {ep}: referência = {ref_values[ep]:.4f} ({direction})")

    # Filtra: para cada endpoint, composto deve ser melhor que referência
    mask = pd.Series([True] * len(df), index=df.index)
    for ep, col in resolved.items():
        direction = ENDPOINT_DIRECTION.get(ep, "lower_better")
        ref_val   = ref_values[ep]

        if direction == "lower_better":
            # Menor é melhor — mantém compostos com valor menor que referência
            ep_mask = df[col] <= ref_val
        else:
            # Maior é melhor — mantém compostos com valor maior que referência
            ep_mask = df[col] >= ref_val

        # Sempre mantém a referência
        if "is_reference" in df.columns:
            ep_mask = ep_mask | (df["is_reference"] == True)
        elif reference_name:
            ep_mask = ep_mask | (df["id"].astype(str) == str(reference_name))

        mask = mask & ep_mask

    df_filtered = df[mask].reset_index(drop=True)
    n_removed = len(df) - len(df_filtered)

    print(f"\n[admet] Filtro por referência:")
    print(f"  Compostos antes  : {len(df)}")
    print(f"  Removidos        : {n_removed}")
    print(f"  Aprovados        : {len(df_filtered)}")

    _save_filtered(df_filtered, outdir, "admet_filtered_by_reference.csv")
    return df_filtered


def _filter_by_thresholds(
    df: pd.DataFrame,
    resolved: dict,
    thresholds_str: str,
    outdir: Path,
) -> pd.DataFrame:
    """
    Filtra compostos por thresholds fixos definidos pelo usuário.
    Formato: "hERG<0.5,DILI<0.3,Bioavailability_Ma>0.6"
    Operadores suportados: <, <=, >, >=
    """
    if not thresholds_str:
        print("[admet] ERRO: --thresholds é obrigatório com --filter-by thresholds.")
        print("  Exemplo: --thresholds 'hERG<0.5,DILI<0.3'")
        return df

    import re
    threshold_rules = {}
    for rule in thresholds_str.split(","):
        rule = rule.strip()
        match = re.match(r"(\w+)\s*(<=|>=|<|>)\s*([\d.]+)", rule)
        if match:
            ep, op, val = match.groups()
            threshold_rules[ep] = (op, float(val))
        else:
            print(f"[admet] AVISO: regra inválida '{rule}'. Use formato 'ENDPOINT<valor'.")

    if not threshold_rules:
        print("[admet] ERRO: nenhuma regra de threshold válida.")
        return df

    print(f"  Regras:")
    for ep, (op, val) in threshold_rules.items():
        print(f"    {ep} {op} {val}")

    mask = pd.Series([True] * len(df), index=df.index)
    for ep, (op, val) in threshold_rules.items():
        # Encontra coluna correspondente
        col = next((c for c in df.columns if ep in c and "percentile" not in c), None)
        if col is None:
            print(f"[admet] AVISO: '{ep}' não encontrado. Ignorando.")
            continue

        if op == "<":    ep_mask = df[col] < val
        elif op == "<=": ep_mask = df[col] <= val
        elif op == ">":  ep_mask = df[col] > val
        elif op == ">=": ep_mask = df[col] >= val
        else:            continue

        mask = mask & ep_mask

    df_filtered = df[mask].reset_index(drop=True)
    n_removed = len(df) - len(df_filtered)

    print(f"\n[admet] Filtro por thresholds:")
    print(f"  Compostos antes  : {len(df)}")
    print(f"  Removidos        : {n_removed}")
    print(f"  Aprovados        : {len(df_filtered)}")

    _save_filtered(df_filtered, outdir, "admet_filtered_by_thresholds.csv")
    return df_filtered


def _save_filtered(df: pd.DataFrame, outdir: Path, filename: str):
    outfile = outdir / filename
    df.to_csv(outfile, index=False)
    print(f"[admet] Filtrados salvos em: {outfile}")

    if len(df) > 0:
        print(f"\n[admet] Top 5 após filtro:")
        key_cols = ["id", "admet_score"] + [
            c for c in df.columns
            if any(ep in c for ep in ["hERG", "DILI", "AMES", "ClinTox", "Bioavail"])
            and "percentile" not in c
        ]
        key_cols = [c for c in key_cols if c in df.columns]
        print(df[key_cols].head(5).to_string(index=False))


# ---------------------------------------------------------------------------
# predict_admet (função principal — sem alterações na lógica)
# ---------------------------------------------------------------------------

def _load_weights(weights_config: str = None) -> dict:
    if weights_config and Path(weights_config).exists():
        import yaml
        config = yaml.safe_load(Path(weights_config).read_text())
        return config.get("admet", DEFAULT_WEIGHTS)
    return DEFAULT_WEIGHTS


def predict_admet(
    hits_file: str,
    outdir: str = "admet",
    weights_config: str = None,
    smiles_column: str = "smiles",
    atc_code: str = None,
) -> pd.DataFrame:
    """
    Roda predições ADMET nos hits usando ADMET-AI.
    """
    if not shutil.which("admet_predict"):
        raise EnvironmentError("admet_predict não encontrado. Instale com: pip install admet-ai")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df_hits = pd.read_csv(hits_file)
    if smiles_column not in df_hits.columns:
        raise ValueError(f"Coluna '{smiles_column}' não encontrada em {hits_file}")

    admet_input  = outdir / "admet_input.csv"
    admet_output = outdir / "admet_raw.csv"

    df_input = df_hits[[smiles_column] + (["id"] if "id" in df_hits.columns else [])].copy()
    df_input = df_input.rename(columns={smiles_column: "smiles"})
    df_input.to_csv(admet_input, index=False)

    print(f"[admet] Rodando ADMET-AI em {len(df_hits)} moléculas...")

    cmd = [
        "admet_predict",
        "--data_path",     str(admet_input),
        "--save_path",     str(admet_output),
        "--smiles_column", "smiles",
    ]
    if atc_code:
        cmd += ["--atc_code", atc_code]

    env = {"CUDA_VISIBLE_DEVICES": "", **__import__("os").environ}
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)

    if not admet_output.exists():
        print(f"[admet] ERRO:\n{result.stderr}")
        raise RuntimeError("admet_predict falhou.")

    df_admet = pd.read_csv(admet_output)
    print(f"[admet] {len(df_admet)} moléculas, {len(df_admet.columns)} propriedades preditas")

    weights = _load_weights(weights_config)
    df_admet = _add_weighted_score(df_admet, weights)

    if "id" in df_hits.columns and "id" in df_admet.columns:
        df_merged = df_hits.merge(
            df_admet.drop(columns=["smiles"], errors="ignore"),
            on="id", how="left"
        )
    else:
        df_merged = pd.concat(
            [df_hits.reset_index(drop=True), df_admet.reset_index(drop=True)], axis=1
        )

    if "admet_score" in df_merged.columns:
        df_merged = df_merged.sort_values("admet_score").reset_index(drop=True)

    df_merged.to_csv(outdir / "admet_results.csv", index=False)

    key_cols = ["id", "smiles", "admet_score"] + [
        c for c in df_merged.columns
        if any(ep in c for ep in ["hERG", "DILI", "Hepatotox", "AMES",
                                   "Bioavail", "Solubility", "ClinTox", "BBB"])
        and "percentile" not in c
    ]
    key_cols = [c for c in key_cols if c in df_merged.columns]
    df_merged[key_cols].to_csv(outdir / "admet_summary.csv", index=False)

    print(f"[admet] Resultados salvos em:")
    print(f"  {outdir}/admet_results.csv")
    print(f"  {outdir}/admet_summary.csv")
    print(f"\n[admet] Top 5 (menor score = melhor perfil ADMET):")
    print(df_merged[key_cols].head(5).to_string(index=False))

    return df_merged


def _add_weighted_score(df: pd.DataFrame, weights: dict) -> pd.DataFrame:
    score = pd.Series(0.0, index=df.index)

    for endpoint, weight in weights.items():
        col = next((c for c in df.columns if endpoint in c), None)
        if col is None or weight == 0:
            continue

        vals = pd.to_numeric(df[col], errors="coerce").fillna(0.5)

        if any(k in endpoint for k in TOX_KEYWORDS):
            score += weight * vals
        elif any(k in endpoint for k in GOOD_KEYWORDS):
            score += weight * (1 - vals)
        else:
            score += weight * vals

    df["admet_score"] = score.round(4)
    return df
