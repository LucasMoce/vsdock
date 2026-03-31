"""
vsdock.prepare
==============
Preparação automática do receptor a partir de um PDB ID ou arquivo PDB local.

Fluxo:
  PDB ID → download RCSB → separação receptor/ligante → PDBQT (obabel)
  arquivo PDB local      → separação receptor/ligante → PDBQT (obabel)
"""

import re
import shutil
import subprocess
import requests
from pathlib import Path


RCSB_URL = "https://files.rcsb.org/download/{pdbid}.pdb"


# ---------------------------------------------------------------------------
# Download do PDB
# ---------------------------------------------------------------------------

def _download_pdb(pdbid: str, outdir: Path) -> Path:
    """Baixa estrutura do RCSB dado um PDB ID."""
    pdbid = pdbid.upper().strip()
    url   = RCSB_URL.format(pdbid=pdbid)

    print(f"[prepare] Baixando {pdbid} do RCSB...")
    resp = requests.get(url, timeout=30)

    if resp.status_code != 200:
        raise RuntimeError(
            f"PDB ID '{pdbid}' não encontrado no RCSB (HTTP {resp.status_code})."
        )

    outfile = outdir / f"{pdbid}.pdb"
    outfile.write_bytes(resp.content)
    print(f"[prepare] PDB salvo em: {outfile}")
    return outfile


# ---------------------------------------------------------------------------
# Separação receptor / ligante
# ---------------------------------------------------------------------------

def _list_heterogens(pdb_file: Path) -> list:
    """Lista todos os HETATM presentes no PDB (exceto água e metais comuns)."""
    skip = {"HOH", "WAT", "H2O", "ZN", "MG", "CA", "NA", "CL", "K", "FE",
            "MN", "CU", "CO", "NI", "SE", "BR", "IOD", "SO4", "PO4", "EDO"}
    found = {}
    for line in pdb_file.read_text().splitlines():
        if line.startswith("HETATM"):
            resname = line[17:20].strip()
            chain   = line[21].strip()
            resnum  = line[22:26].strip()
            if resname not in skip:
                key = f"{resname}{resnum}{chain}"
                found[key] = {"resname": resname, "chain": chain, "resnum": resnum}
    return list(found.values())


def _separate_receptor_ligand(
    pdb_file: Path,
    ligand_code: str,
    outdir: Path,
) -> tuple:
    """
    Separa receptor e ligante num arquivo PDB.

    Parâmetros
    ----------
    pdb_file     : Path do PDB completo
    ligand_code  : código do ligante (ex: "MRV") — apenas o nome de 3 letras
    outdir       : pasta de saída

    Retorna
    -------
    (receptor_pdb, ligand_pdb) — Paths dos arquivos gerados
    """
    receptor_lines = []
    ligand_lines   = []

    ligand_code = ligand_code.upper().strip()

    for line in pdb_file.read_text().splitlines():
        record = line[:6].strip()

        if record == "HETATM":
            resname = line[17:20].strip()
            if resname == ligand_code:
                ligand_lines.append(line)
            else:
                # Mantém outros HETATM que não são o ligante alvo
                # (ex: cofatores, mas remove água)
                if resname not in {"HOH", "WAT", "H2O"}:
                    receptor_lines.append(line)
        elif record in {"ATOM", "TER"}:
            receptor_lines.append(line)
        elif record in {"REMARK", "HEADER", "TITLE", "SEQRES", "CRYST1"}:
            receptor_lines.append(line)

    if not ligand_lines:
        raise ValueError(
            f"Ligante '{ligand_code}' não encontrado no PDB.\n"
            f"Use 'grep HETATM {pdb_file.name} | awk '{{print $4}}' | sort -u' "
            f"para ver os ligantes disponíveis."
        )

    receptor_pdb = outdir / "receptor_clean.pdb"
    ligand_pdb   = outdir / f"{ligand_code.lower()}_crystal.pdb"

    receptor_pdb.write_text("\n".join(receptor_lines) + "\nEND\n")
    ligand_pdb.write_text("\n".join(ligand_lines) + "\nEND\n")

    print(f"[prepare] Receptor salvo em : {receptor_pdb}")
    print(f"[prepare] Ligante salvo em  : {ligand_pdb} ({len(ligand_lines)} átomos)")

    return receptor_pdb, ligand_pdb


# ---------------------------------------------------------------------------
# Conversão para PDBQT
# ---------------------------------------------------------------------------

def _to_pdbqt(pdb_file: Path, outfile: Path, is_receptor: bool = True):
    """Converte PDB para PDBQT via Open Babel."""
    if not shutil.which("obabel"):
        raise EnvironmentError("obabel não encontrado no PATH.")

    flags = ["-xr"] if is_receptor else []
    cmd = ["obabel", str(pdb_file), "-O", str(outfile)] + flags

    result = subprocess.run(cmd, capture_output=True, text=True)

    if not outfile.exists():
        raise RuntimeError(
            f"Conversão para PDBQT falhou:\n{result.stderr}"
        )
    print(f"[prepare] PDBQT gerado: {outfile}")


# ---------------------------------------------------------------------------
# Função principal
# ---------------------------------------------------------------------------

def prepare_receptor(
    ligand_code: str = None,
    pdbid: str = None,
    pdb_file: str = None,
    outdir: str = ".",
    keep_heterogens: bool = False,
) -> dict:
    """
    Prepara receptor a partir de PDB ID ou arquivo local.
    Se ligand_code for omitido, prepara apenas o receptor sem separar ligante.

    Parâmetros
    ----------
    ligand_code : código de 3 letras do ligante no PDB (ex: "MRV"). Opcional.
    pdbid       : PDB ID para download automático (ex: "4MBS")
    pdb_file    : caminho para arquivo PDB local (alternativa ao pdbid)
    outdir      : pasta do projeto

    Retorna
    -------
    dict com paths gerados
    """
    if not pdbid and not pdb_file:
        raise ValueError("Informe --pdbid ou --pdb-file.")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. Obtém o PDB
    if pdbid:
        pdb_path = _download_pdb(pdbid, outdir)
    else:
        pdb_path = Path(pdb_file)
        if not pdb_path.exists():
            raise FileNotFoundError(f"Arquivo não encontrado: {pdb_file}")
        print(f"[prepare] Usando arquivo local: {pdb_path}")

    # 2. Lista heterogens disponíveis
    heterogens = _list_heterogens(pdb_path)
    if heterogens:
        print(f"[prepare] Ligantes encontrados no PDB:")
        for h in heterogens[:10]:
            print(f"  {h['resname']}{h['resnum']}{h['chain']}")
        if len(heterogens) > 10:
            print(f"  ... e mais {len(heterogens) - 10}")
    else:
        print(f"[prepare] Nenhum ligante encontrado no PDB (apenas receptor).")

    result = {"pdb": pdb_path}

    # 3. Separa receptor e ligante (se ligand_code fornecido)
    if ligand_code:
        receptor_pdb, ligand_pdb = _separate_receptor_ligand(
            pdb_path, ligand_code, outdir
        )
        ligand_ids = [
            f"{h['resname']}{h['resnum']}{h['chain']}"
            for h in heterogens
            if h["resname"] == ligand_code.upper()
        ]
        ligand_id = ligand_ids[0] if ligand_ids else f"{ligand_code}1A"
        result["ligand_pdb"] = ligand_pdb
        result["ligand_id"]  = ligand_id
    else:
        # Sem ligante — prepara receptor removendo apenas água
        print(f"[prepare] Modo sem ligante — preparando receptor completo...")
        receptor_lines = []
        for line in pdb_path.read_text().splitlines():
            record = line[:6].strip()
            if record in {"ATOM", "TER", "REMARK", "HEADER", "TITLE", "SEQRES", "CRYST1"}:
                receptor_lines.append(line)
            elif record == "HETATM":
                resname = line[17:20].strip()
                if resname not in {"HOH", "WAT", "H2O"}:
                    receptor_lines.append(line)
        receptor_pdb = outdir / "receptor_clean.pdb"
        receptor_pdb.write_text("\n".join(receptor_lines) + "\nEND\n")
        print(f"[prepare] Receptor salvo em: {receptor_pdb}")
        ligand_id = None

    # 4. Converte receptor para PDBQT
    receptor_pdbqt = outdir / "receptor.pdbqt"
    _to_pdbqt(receptor_pdb, receptor_pdbqt, is_receptor=True)

    # 5. Instrução de uso
    print(f"\n[prepare] Pronto! Use nas próximas etapas:")
    print(f"  --receptor {receptor_pdbqt}")
    print(f"  --pdb      {pdb_path}")
    if ligand_id:
        print(f"  --autobox-ligand {ligand_id}")
    else:
        print(f"  --blind                     (blind docking)")
        print(f"  --autobox-residues RES1 RES2 (resíduos do sítio ativo)")

    # 6. Atualiza state
    state_file = outdir / "vsdock_state.yaml"
    if state_file.exists():
        import yaml
        state = yaml.safe_load(state_file.read_text())
        state["receptor"] = str(receptor_pdbqt)
        state["pdb"]      = str(pdb_path)
        if ligand_id:
            state["ligand_crystal_id"] = ligand_id
        state_file.write_text(yaml.dump(state))
        print(f"[prepare] vsdock_state.yaml atualizado.")

    result.update({
        "receptor_pdb":   receptor_pdb,
        "receptor_pdbqt": receptor_pdbqt,
    })
    return result
