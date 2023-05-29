# nrpytutorial

Refactor and extend project structure following modern day software architecture
concepts and implementations.

## Quickstart

**Bias warning**: Improves focus is Scalar Field collapse

```
# clone repo
git clone https://github.com/dpopchev/nrpytutorial.git
# work version is python3.8; assuming version control is done by pyenv
pyseed=~/.pyenv/versions/3.8.16/bin/python make development
# populate venv with nrpy modules
make install-nrpy-modules
# populate venv with progress scalar field collapse module dependencies
make install-nrpy-progress-sfcollapse-modules
# generate in progress nrpytutorial for scalar field collapse
make run-nrpy-progress-sfcollapse-bssn
# print help working with project; recommend only use make
make venv
# some more help
make help
# play around with moduels trough CLI
make run-ipython
# play around with notebooks trough Jupyter
make run-jupyter
```

Using bash autocompletion will help tremendously anyone exploring the `Makefile`.

I ensure anyone using the `install-` targets that there is an counter `uninstall-`
that cleans up after itself.

Find in `data/runners/` scripts used by `make` to run executable products.

The default runner takes two arguments, the configuration and executable;
The configuration should be in `csv` format. The script iterates over the rows
and passes the values to the executable preserving the order.

Outputs are stored into `build/.../output/`; Backups are made if runs already were
made.

The figs are generated using the output into `build/.../output/figs`.
