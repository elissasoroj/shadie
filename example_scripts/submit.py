#!/usr/bin/env python

"""Create replicates and submit to run.
"""


import os
import subprocess
import typer


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
app = typer.Typer(add_completion=False, context_settings=CONTEXT_SETTINGS)


@app.command()
def info(
    organism: str = typer.Argument(..., help="organism name"),
    mode: str = typer.Argument(..., help="run type"),    
    idx: int = typer.Argument(..., help="replicate number"),
    outdir: str = typer.Argument(..., help="directory"),
    slim_binary: str = typer.Argument(..., help="path to slim"),
    ):
    """Start a SliM job on a subprocess."""
    script = get_slim_script(organism, mode)
    script = substitute_outfile(script, outdir, idx)
    typer.secho(f"starting job: {script}", fg=typer.colors.CYAN)
    start_slim_job(script, idx, slim_binary)


def get_slim_script(organism: str, mode: str):
    """Return the script for a selected type."""
    script = f"./target-{organism}-{mode}.slim"
    if not os.path.exists(script):
        raise IOError(f"{script} not found")
    return script


def substitute_outfile(script: str, outdir: str, idx: int):
    """Creates a tmp copy of the file with the outfile modified."""
    with open(script, 'r') as indata:
        lines = indata.readlines()
        for ldx, content in enumerate(lines):
            if content.startswith("\tdefineConstant('OUTPATH',"):
                outfile = content.split("'OUTPATH',")[-1]
                outfile = outfile.rsplit("-", 1)[0]
                outfile = outfile + f"-{idx}.trees"
                lines[ldx] = f"\tdefineConstant('OUTPATH',{outfile}');"

    outname = script.replace(".slim", f"-{idx}.slim")
    outpath = os.path.join(outdir, outname)
    with open(outpath, 'wt') as out:
        out.write("".join(lines))
    return outpath


def start_slim_job(script: str, idx: int, slim_binary: str):
    """Starts a replicate job with a different seed."""
    if not idx:
        seed = "123"
    else:
        seed = f"{idx}{idx}{idx}"
    cmd = [slim_binary, "-s", seed, script]
    return subprocess.run(cmd, check=True)    


if __name__ == "__main__":
    app()
