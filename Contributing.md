# Contributing

## Local dev setup

Install poetry:

```shell
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/install-poetry.py | python -
# Make sure you have ~/.local/bin in PATH
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

We clone to _directory to avoid python from importing from the directory.

```
%%bash

git clone https://github.com/konstin/Colabfold _colabfold
pip install --use-feature=in-tree-build _colabfold
# Unholy Hack: Use the files from our cloned git repository instead of installed copy
site_packages=$(python -c 'import site; print(site.getsitepackages()[0])')
rm -r ${site_packages}/colabfold
ln -s $(pwd)/_colabfold/colabfold ${site_packages}/colabfold
```

If you also need to patch alphafold:

```
%%bash

pip uninstall -y alphafold
git clone https://github.com/konstin/alphafold _alphafold
pip install -e ./_alphafold
```

After that, restart the runtime.

When you changed a file in the `colabfold` package, you need to reload the modules you were using with `importlib.reload()`, e.g.

```python
import colabfold.batch
import importlib

importlib.reload(colabfold.batch)
```
