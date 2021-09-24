# Contributing

## Local dev setup

Install poetry:

```shell
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/install-poetry.py | python -
# Make sure you have or put ~/.local/bin in PATH
poetry config settings.virtualenvs.in-project true
```

Setup virtualenv, install dependencies:

```shell
poetry install
```

Whenever dependencies change, run `poetry install` again. You can add dependencies with `poetry add <package>`.

In your IDE select `.venv/bin/python` as interpreter. In a shell, you can activate the environment with `. .venv/bin/activate` and deactivate it with `deactivate`. To run in you IDE, select `colabfold.batch` as module to run and the git root as working directory.

To run the tests: 

```shell
pytest
```

Format the code:

```shell
black .
```

## Colab dev setup

**Note: This is work in progress**

```shell
%%bash

#curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/install-poetry.py | python -
#/root/.local/bin/poetry config virtualenvs.create false # Install in the colab env
git clone https://github.com/konstin/Colabfold
(cd Colabfold; pip install .; cd ..)
```

If you also need to patch alphafold:

```shell
%%bash

pip uninstall -y alphafold
git clone https://github.com/konstin/alphafold
(cd alphafold; python setup.py develop; cd ..)
```