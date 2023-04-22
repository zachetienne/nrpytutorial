# nrpytutorial fork

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
# print help working with project; recommend only use make
make venv
# some more help 
make help
# play around with moduels trough CLI
make run-ipython
# play around with notebooks trough Jupyter
make run-jupyter
```
