"""
vsdock.admet
============
Predição de propriedades ADMET usando ADMET-AI (offline, via CLI).

Referência: Swanson et al., Bioinformatics 2024
"""

import subprocess
import shutil
import pandas as pd
from pathlib import Path


ADMET_ENDPOINTS = {
    "absorption": [
        "Caco2_Wang", "HIA_Hou", "Pgp_Broccatelli", "Bioavailability_Ma",
        "Lipophilicity_AstraZeneca", "Solubility_AqSolDB",
    ],
    "distribution": ["BBB_Martini", "PPBR_AZ", "VDss_Lombardo"],
    "metabolism": [
        "CYP1A2_Veith", "CYP2C19_Veith", "CYP2C9_Veith",
        "CYP2D6_Veith", "CYP3A4_Veith",
    ],
    "excretion": ["Half_Life_Obach", "Clearance_Hepatocyte_AZ", "Clearance_Microsome_AZ"],
    "toxicity": [
        "hERG", "AMES", "DILI", "Hepatotoxicity_Xu",
        "Skin_Reaction", "ClinTox", "LD50_Zhu",
    ],
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

    Parâmetros
    ----------
    hits_file      : CSV com coluna smiles
    outdir         : pasta de saída
    weights_config : YAML com pesos por endpoint (usa default se omitido)
    smiles_column  : nome da coluna de SMILES no CSV
    atc_code       : código ATC para comparação com DrugBank (ex: "J05" para antivirais)

    Retorna
    -------
    DataFrame com predições ADMET + score ponderado
    """
    if not shutil.which("admet_predict"):
        raise EnvironmentError("admet_predict não encontrado. Instale com: pip install admet-ai")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df_hits = pd.read_csv(hits_file)
    if smiles_column not in df_hits.columns:
        raise ValueError(f"Coluna '{smiles_column}' não encontrada em {hits_file}")

    # Prepara input para ADMET-AI
    admet_input  = outdir / "admet_input.csv"
    admet_output = outdir / "admet_raw.csv"

    df_input = df_hits[[smiles_column] + (["id"] if "id" in df_hits.columns else [])].copy()
    df_input = df_input.rename(columns={smiles_column: "smiles"})
    df_input.to_csv(admet_input, index=False)

    print(f"[admet] Rodando ADMET-AI em {len(df_hits)} moléculas...")

    cmd = [
        "admet_predict",
        "--data_path",    str(admet_input),
        "--save_path",    str(admet_output),
        "--smiles_column","smiles",
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

    # Score ponderado
    weights = _load_weights(weights_config)
    df_admet = _add_weighted_score(df_admet, weights)

    # Merge com dados originais
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

    # Salva completo
    df_merged.to_csv(outdir / "admet_results.csv", index=False)

    # Resumo com endpoints chave
    key_cols = ["id", "smiles", "admet_score"] + [
        c for c in df_merged.columns
        if any(ep in c for ep in ["hERG", "DILI", "Hepatotox", "AMES",
                                   "Bioavail", "Solubility", "ClinTox", "BBB"])
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
    """
    Score ADMET ponderado. Menor = melhor.
    Toxicidade: penaliza probabilidade alta.
    Absorção/biodisponibilidade: penaliza valores baixos.
    """
    tox_keywords  = ["hERG", "AMES", "DILI", "Hepatotox", "ClinTox",
                     "Skin", "Carcino", "NR-", "SR-"]
    good_keywords = ["Bioavailability", "HIA", "Solubility"]

    score = pd.Series(0.0, index=df.index)

    for endpoint, weight in weights.items():
        col = next((c for c in df.columns if endpoint in c), None)
        if col is None or weight == 0:
            continue

        vals = pd.to_numeric(df[col], errors="coerce").fillna(0.5)

        if any(k in endpoint for k in tox_keywords):
            score += weight * vals
        elif any(k in endpoint for k in good_keywords):
            score += weight * (1 - vals)
        else:
            score += weight * vals

    df["admet_score"] = score.round(4)
    return df
