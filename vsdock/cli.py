"""
vsdock CLI — entry point principal
"""

import argparse
from vsdock import __version__


def cmd_init(args):
    from vsdock.fetch import fetch_ligand_smiles
    import yaml
    from pathlib import Path

    print(f"[vsdock] Inicializando projeto '{args.ligand}' ...")

    dirs = ["ligand", "zinc", "hits", "docking", "plip", "admet", "report"]
    for d in dirs:
        Path(d).mkdir(exist_ok=True)

    smiles = fetch_ligand_smiles(args.ligand)
    state = {
        "ligand_name": args.ligand,
        "ligand_smiles": smiles,
        "receptor": args.target,
        "config": args.config,
    }
    Path("vsdock_state.yaml").write_text(yaml.dump(state))
    print(f"[vsdock] Projeto inicializado. Estado salvo em vsdock_state.yaml")


def cmd_fetch(args):
    from vsdock.fetch import fetch_ligand, fetch_database

    if args.ligand:
        fetch_ligand(args.ligand, outdir="ligand", fmt=args.fmt)

    if args.database:
        fetch_database(
            outdir="zinc",
            source=args.source,
            mw_min=args.mw_min,
            mw_max=args.mw_max,
            logp_min=args.logp_min,
            logp_max=args.logp_max,
            max_mols=args.max_mols,
        )


def cmd_screen(args):
    from vsdock.screen import screen, load_query_from_state
    import glob

    query_smiles = load_query_from_state()

    db_file = args.database
    if not db_file:
        candidates = glob.glob("zinc/*.smi") + glob.glob("zinc/*.csv")
        if not candidates:
            print("[screen] ERRO: nenhum banco encontrado em zinc/. Rode 'vsdock fetch --database' primeiro.")
            raise SystemExit(1)
        db_file = candidates[0]
        print(f"[screen] Banco detectado: {db_file}")

    screen(
        query_smiles=query_smiles,
        database_file=db_file,
        outdir="hits",
        threshold=args.similarity,
        fp_type=args.fingerprint,
        radius=args.radius,
        max_hits=args.max_hits,
    )


def cmd_dock(args):
    print(f"[vsdock] Docking (exhaustiveness={args.exhaustiveness}, top={args.top})")
    print("  -> modulo dock ainda nao implementado")


def cmd_analyze(args):
    from vsdock.pains import filter_pains
    import os

    hits_file = args.hits_file
    if not hits_file:
        candidates = ["hits/hits.csv", "hits/hits_clean.csv"]
        for c in candidates:
            if os.path.exists(c):
                hits_file = c
                break
        if not hits_file:
            print("[analyze] ERRO: nenhum arquivo de hits encontrado. Rode 'vsdock screen' primeiro.")
            raise SystemExit(1)
        print(f"[analyze] Arquivo detectado: {hits_file}")

    filter_pains(
        hits_file=hits_file,
        outdir="hits",
        also_filter_lipinski=not args.no_lipinski,
    )


def cmd_report(args):
    print(f"[vsdock] Gerando relatorio (formato: {args.format})")
    print("  -> modulo report ainda nao implementado")


def main():
    parser = argparse.ArgumentParser(
        prog="vsdock",
        description="Virtual Screening and Docking Pipeline",
    )
    parser.add_argument("--version", action="version", version=f"vsdock {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")
    subparsers.required = True

    # --- init ---
    p_init = subparsers.add_parser("init", help="Inicializa projeto e baixa ligante")
    p_init.add_argument("--target", required=True, help="Arquivo .pdbqt do receptor")
    p_init.add_argument("--ligand", required=True, help="Nome do ligante (PubChem)")
    p_init.add_argument("--config", default="configs/default.yaml")
    p_init.set_defaults(func=cmd_init)

    # --- fetch ---
    p_fetch = subparsers.add_parser("fetch", help="Baixa ligante e/ou banco molecular")
    p_fetch.add_argument("--ligand", default=None, help="Nome do ligante (PubChem)")
    p_fetch.add_argument("--fmt", choices=["sdf", "smiles"], default="sdf")
    p_fetch.add_argument("--database", action="store_true", help="Baixar banco molecular")
    p_fetch.add_argument("--source", choices=["chembl", "zinc"], default="chembl")
    p_fetch.add_argument("--mw-min", type=float, default=150, dest="mw_min")
    p_fetch.add_argument("--mw-max", type=float, default=500, dest="mw_max")
    p_fetch.add_argument("--logp-min", type=float, default=-1, dest="logp_min")
    p_fetch.add_argument("--logp-max", type=float, default=5, dest="logp_max")
    p_fetch.add_argument("--max-mols", type=int, default=500, dest="max_mols")
    p_fetch.set_defaults(func=cmd_fetch)

    # --- screen ---
    p_screen = subparsers.add_parser("screen", help="Triagem por similaridade estrutural")
    p_screen.add_argument("--similarity", type=float, default=0.4)
    p_screen.add_argument("--fingerprint", choices=["morgan", "maccs", "rdkit"], default="morgan")
    p_screen.add_argument("--radius", type=int, default=2)
    p_screen.add_argument("--max-hits", type=int, default=500, dest="max_hits")
    p_screen.add_argument("--database", default=None)
    p_screen.set_defaults(func=cmd_screen)

    # --- dock ---
    p_dock = subparsers.add_parser("dock", help="Docking com AutoDock Vina")
    p_dock.add_argument("--exhaustiveness", type=int, default=8)
    p_dock.add_argument("--top", type=int, default=20)
    p_dock.set_defaults(func=cmd_dock)

    # --- analyze ---
    p_analyze = subparsers.add_parser("analyze", help="Filtro PAINS + Lipinski")
    p_analyze.add_argument("--hits-file", default=None, dest="hits_file",
                           help="CSV de hits (detectado automaticamente se omitido)")
    p_analyze.add_argument("--no-lipinski", action="store_true", dest="no_lipinski",
                           help="Desativa filtro de Lipinski")
    p_analyze.set_defaults(func=cmd_analyze)

    # --- report ---
    p_report = subparsers.add_parser("report", help="Gera relatorio/manuscrito")
    p_report.add_argument("--format", choices=["quarto", "markdown", "html"], default="quarto")
    p_report.set_defaults(func=cmd_report)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
