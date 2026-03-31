"""
vsdock CLI — entry point principal
"""

import argparse
import sys
from vsdock import __version__


def cmd_init(args):
    print(f"[vsdock] Inicializando projeto...")
    print(f"  Receptor : {args.target}")
    print(f"  Ligante  : {args.ligand}")
    print(f"  Config   : {args.config}")
    print("[vsdock] Pronto. Edite o arquivo de config antes de continuar.")


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
    p_init = subparsers.add_parser("init", help="Inicializa um novo projeto")
    p_init.add_argument("--target", required=True, help="Arquivo .pdbqt do receptor")
    p_init.add_argument("--ligand", required=True, help="Nome ou CID do ligante (PubChem)")
    p_init.add_argument("--config", default="configs/default.yaml", help="Arquivo de configuração")
    p_init.set_defaults(func=cmd_init)

    # --- screen ---
    p_screen = subparsers.add_parser("screen", help="Triagem por similaridade estrutural")
    p_screen.add_argument("--similarity", type=float, default=0.4, help="Threshold de similaridade (Tanimoto)")
    p_screen.add_argument("--zinc-filters", default="", help="Filtros ZINC (ex: 'mw:150-500,logp:-1-5')")
    p_screen.set_defaults(func=cmd_screen)

    # --- dock ---
    p_dock = subparsers.add_parser("dock", help="Docking com AutoDock Vina")
    p_dock.add_argument("--exhaustiveness", type=int, default=8)
    p_dock.add_argument("--top", type=int, default=20, help="Número de hits para dockar")
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
