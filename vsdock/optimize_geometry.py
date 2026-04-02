"""
vsdock.optimize_geometry
========================
Otimização de geometria e minimização de energia dos ligantes.

Usa RDKit com campo de força MMFF94 para:
  1. Gerar conformação 3D inicial (ETKDGv3)
  2. Minimizar energia com MMFF94
  3. Calcular e reportar a energia final

Isso garante que todos os ligantes estão em geometria otimizada
antes do docking, melhorando a qualidade dos resultados.
"""

import pandas as pd
from pathlib import Path
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def optimize_geometry(
    hits_file: str,
    outdir: str = "hits",
    force_field: str = "mmff94",
    max_iters: int = 2000,
    remove_salts: bool = True,
) -> pd.DataFrame:
    """
    Otimiza a geometria de todos os compostos da biblioteca.

    Parâmetros
    ----------
    hits_file    : CSV com coluna smiles (saída do clear-library)
    outdir       : pasta de saída
    force_field  : "mmff94" (padrão) ou "uff"
    max_iters    : número máximo de iterações de minimização
    remove_salts : remove sais/contrafragmentos do SMILES antes de otimizar

    Retorna
    -------
    DataFrame com colunas adicionais: smiles_3d, energy_mmff94, opt_converged
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(hits_file)
    print(f"[optimize] {len(df)} moléculas carregadas")
    print(f"  Campo de força : {force_field.upper()}")
    print(f"  Max iterações  : {max_iters}")

    energies    = []
    converged   = []
    smiles_opt  = []
    failed_ids  = []

    sdf_path = outdir / "hits_optimized.sdf"
    writer   = Chem.SDWriter(str(sdf_path))

    for _, row in tqdm(df.iterrows(), total=len(df), desc="Otimizando geometria", unit=" mol"):
        smiles = row["smiles"]
        mol_id = str(row.get("id", ""))

        # Remove sais
        if remove_salts and "." in smiles:
            frags = smiles.split(".")
            best  = max(frags, key=lambda s: len(Chem.MolFromSmiles(s).GetAtoms())
                        if Chem.MolFromSmiles(s) else 0)
            smiles = best

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            energies.append(None)
            converged.append(False)
            smiles_opt.append(smiles)
            failed_ids.append(mol_id)
            continue

        mol = Chem.AddHs(mol)

        # Gera conformação 3D
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if result != 0:
            result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        if result != 0:
            energies.append(None)
            converged.append(False)
            smiles_opt.append(smiles)
            failed_ids.append(mol_id)
            continue

        # Minimiza energia
        if force_field.lower() == "mmff94":
            ff = AllChem.MMFFGetMoleculeForceField(
                mol, AllChem.MMFFGetMoleculeProperties(mol)
            )
        else:  # UFF
            ff = AllChem.UFFGetMoleculeForceField(mol)

        if ff is None:
            energies.append(None)
            converged.append(False)
            smiles_opt.append(smiles)
            failed_ids.append(mol_id)
            continue

        conv = ff.Minimize(maxIts=max_iters)
        energy = ff.CalcEnergy()

        energies.append(round(energy, 4))
        converged.append(conv == 0)  # 0 = convergiu

        # SMILES da geometria otimizada (sem H)
        mol_noH = Chem.RemoveHs(mol)
        smiles_opt.append(Chem.MolToSmiles(mol_noH))

        # Adiciona propriedades e salva no SDF
        mol.SetProp("_Name", mol_id)
        mol.SetProp("ID", mol_id)
        mol.SetProp(f"Energy_{force_field.upper()}", str(round(energy, 4)))
        mol.SetProp("Converged", str(conv == 0))
        writer.write(mol)

    writer.close()

    df[f"energy_{force_field}"] = energies
    df["opt_converged"]         = converged
    df["smiles_optimized"]      = smiles_opt

    # Estatísticas
    n_ok     = sum(1 for e in energies if e is not None)
    n_conv   = sum(1 for c in converged if c)
    n_failed = len(failed_ids)

    print(f"\n[optimize] Resultados:")
    print(f"  Otimizadas com sucesso : {n_ok} / {len(df)}")
    print(f"  Convergiram            : {n_conv} / {n_ok}")
    if n_failed:
        print(f"  Falharam               : {n_failed}")

    if n_ok > 0:
        valid_energies = [e for e in energies if e is not None]
        print(f"  Energia média          : {sum(valid_energies)/len(valid_energies):.2f} kcal/mol")
        print(f"  Energia mínima         : {min(valid_energies):.2f} kcal/mol")
        print(f"  Energia máxima         : {max(valid_energies):.2f} kcal/mol")

    # Salva CSV atualizado
    out_csv = outdir / "hits_optimized.csv"
    df.to_csv(out_csv, index=False)

    # Salva também o .smi com SMILES otimizados
    out_smi = outdir / "hits_optimized.smi"
    df_ok = df[df[f"energy_{force_field}"].notna()]
    df_ok[["smiles_optimized", "id"]].to_csv(
        out_smi, sep="\t", index=False, header=False
    )

    print(f"\n[optimize] Arquivos salvos:")
    print(f"  {out_csv}              (CSV com energias)")
    print(f"  {sdf_path}   (SDF 3D para visualização)")
    print(f"  {out_smi}    (SMILES otimizados)")

    return df
