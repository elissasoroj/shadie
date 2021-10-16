#!/usr/bin/env python

"""Create replicates and submit to run.
"""


import os
import subprocess
import typer


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
app = typer.Typer(add_completion=False, context_settings=CONTEXT_SETTINGS)


@app.command()
def run(
    organism: str = typer.Argument(..., help="organism name"),
    mode: str = typer.Argument(..., help="run type"),    
    idx: int = typer.Argument(..., help="replicate number"),
    outdir: str = typer.Argument(..., help="directory"),
    slim_binary: str = typer.Argument(..., help="path to slim"),
    sim_time: int = typer.Argument(10001, help="length of simulation"),
    ):
    """Start a SliM job on a subprocess."""
    script = get_slim_script(organism, mode)
    script = substitute_outfile(script, outdir, idx, sim_time)
    typer.secho(f"starting job: {script}", fg=typer.colors.CYAN)
    start_slim_job(script, idx, slim_binary)


def get_slim_script(organism: str, mode: str):
    """Return the script for a selected type."""
    script = f"./target-{organism}-{mode}.slim"
    if not os.path.exists(script):
        raise IOError(f"{script} not found")
    return script


def substitute_outfile(script: str, outdir: str, idx: int, sim_time: int):
    """Creates a tmp copy of the file with the outfile modified."""
    name_script = script.replace(".slim", f"-{idx}.slim")
    name_trees = script.replace(".slim", f"-{idx}.trees")

    with open(script, 'r') as indata:
        lines = indata.readlines()
        for ldx, content in enumerate(lines):

            if "defineConstant('OUTPATH'," in content:
                outpath = os.path.join(outdir, name_trees)
                lines[ldx] = f"\tdefineConstant('OUTPATH', '{outpath}');\n"

            if content.startswith("10001 late() {"):
                lines[ldx] = content.replace("10001", str(sim_time))

    outpath = os.path.join(outdir, name_script)
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

    # "bryo-dio", "real", 2, "/tmp", "/usr/local/bin/slim", 10001)
    app()
