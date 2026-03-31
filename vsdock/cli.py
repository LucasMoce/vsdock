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

    # Cria estrutura de pastas do projeto
    dirs = ["ligand", "zinc", "hits", "docking", "plip", "admet", "report"]
    for d in dirs:
        Path(d).mkdir(exist_ok=True)

    # Baixa SMILES do ligante e salva em arquivo de estado
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
    print(f"[vsdock] Screening (similaridade >= {args.similarity})")
    print("  → módulo screen ainda não implementado")


def cmd_dock(args):
    print(f"[vsdock] Docking (exhaustiveness={args.exhaustiveness}, top={args.top})")
    print("  → módulo dock ainda não implementado")


def cmd_analyze(args):
    print("[vsdock] Análise PAINS + PLIP + ADMET")
    print("  → módulo analyze ainda não implementado")


def cmd_report(args):
    print(f"[vsdock] Gerando relatório (formato: {args.format})")
    print("  → módulo report ainda não implementado")


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
    p_fetch.add_argument("--source", choices=["chembl", "zinc"], default="chembl",
                         help="Fonte do banco (default: chembl)")
    p_fetch.add_argument("--mw-min", type=float, default=150, dest="mw_min")
    p_fetch.add_argument("--mw-max", type=float, default=500, dest="mw_max")
    p_fetch.add_argument("--logp-min", type=float, default=-1, dest="logp_min")
    p_fetch.add_argument("--logp-max", type=float, default=5, dest="logp_max")
    p_fetch.add_argument("--max-mols", type=int, default=500, dest="max_mols")
    p_fetch.set_defaults(func=cmd_fetch)

    # --- screen ---
    p_screen = subparsers.add_parser("screen", help="Triagem por similaridade")
    p_screen.add_argument("--similarity", type=float, default=0.4)
    p_screen.add_argument("--zinc-filters", default="")
    p_screen.set_defaults(func=cmd_screen)

    # --- dock ---
    p_dock = subparsers.add_parser("dock", help="Docking com AutoDock Vina")
    p_dock.add_argument("--exhaustiveness", type=int, default=8)
    p_dock.add_argument("--top", type=int, default=20)
    p_dock.set_defaults(func=cmd_dock)

    # --- analyze ---
    p_analyze = subparsers.add_parser("analyze", help="PAINS + PLIP + ADMET")
    p_analyze.add_argument("--admet-weights", default="configs/default.yaml")
    p_analyze.set_defaults(func=cmd_analyze)

    # --- report ---
    p_report = subparsers.add_parser("report", help="Gera relatório/manuscrito")
    p_report.add_argument("--format", choices=["quarto", "markdown", "html"], default="quarto")
    p_report.set_defaults(func=cmd_report)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
