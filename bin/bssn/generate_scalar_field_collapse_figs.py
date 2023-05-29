import matplotlib.pyplot as plt
from pathlib import Path
import csv
import numpy as np

SOURCES = Path('build/bssn/scalar_field_collapse/output/')
RUNS_CONFIG = SOURCES.joinpath('runs.config')

def get_runs(src: Path):
    runs = [ l for l in src.read_text().splitlines()]
    return csv.DictReader(runs, skipinitialspace=True)

def make_plot(fig, data):
    x,r,sf,sfM,alpha,cf,logH = data
    ax  = fig.add_subplot(321)
    ax2 = fig.add_subplot(322)
    ax3 = fig.add_subplot(323)
    ax4 = fig.add_subplot(324)
    ax5 = fig.add_subplot(325)

    ax.set_title("Scalar field")
    ax.set_ylabel(r"$\varphi(t,r)$")
    ax.set_xlim(0,20)
    ax.set_ylim(-0.6,0.6)
    ax.plot(r,sf,'k-')
    ax.grid()

    ax2.set_title("Scalar field conjugate momentum")
    ax2.set_ylabel(r"$\Pi(t,r)$")
    ax2.set_xlim(0,20)
    ax2.set_ylim(-1,1)
    ax2.plot(r,sfM,'b-')
    ax2.grid()

    ax3.set_title("Lapse function")
    ax3.set_ylabel(r"$\alpha(t,r)$")
    ax3.set_xlim(0,20)
    ax3.set_ylim(0,1.02)
    ax3.plot(r,alpha,'r-')
    ax3.grid()

    ax4.set_title("Conformal factor")
    ax4.set_xlabel(r"$r$")
    ax4.set_ylabel(r"$W(t,r)$")
    ax4.set_xlim(0,20)
    ax4.set_ylim(0,1.02)
    ax4.plot(r,cf,'g-',label=("$p = 0.043149493$"))
    ax4.grid()

    ax5.set_title("Hamiltonian constraint violation")
    ax5.set_xlabel(r"$r$")
    ax5.set_ylabel(r"$\mathcal{H}(t,r)$")
    ax5.set_xlim(0,20)
    ax5.set_ylim(-16,0)
    ax5.plot(r,logH,'m-')
    ax5.grid()

def plot_results(results: Path, Nx0):
    figs_home = results.joinpath('figs')
    figs_home.mkdir(exist_ok=True)
    for result in sorted(results.rglob(f'out{Nx0}-00*.txt')):
        data = np.loadtxt(result).T
        fig = plt.figure(figsize=(8,6))
        make_plot(fig, data)
        plt.tight_layout()
        plt.savefig(figs_home.joinpath(result.stem+".png"),dpi=150)
        plt.close(fig)

if __name__ == '__main__':
    runs = get_runs(RUNS_CONFIG)
    for count, run in enumerate(runs, start=1):
        run_name = run['name'] if run['name'] else count
        run_results = SOURCES.joinpath(f"run-{run_name}")
        plot_results(run_results, run['Nx0'])
