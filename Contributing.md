# Contributing to colabfold and our alphafold fork

Install poetry (once per machine). Please consult the [poetry docs](https://python-poetry.org/docs/), they are well written:

```shell
curl -sSL https://install.python-poetry.org | python3 -
poetry config virtualenvs.in-project true
```

Clone and install the dependencies:

```shell
git clone https://github.com/sokrypton/ColabFold
cd ColabFold
poetry install -E alphafold
```

Activate the environment (everytime you want to run some python or install something):

```shell
source .venv/bin/activate
```

Install jax; You need to repeat after every `poetry install`/`poetry lock`/`poetry update` unfortunately:

```shell
pip install -q "jax[cuda]>=0.3.8,<0.4" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
```

If you also want to modify our alphafold fork

```shell
git clone https://github.com/steineggerlab/alphafold
pip install -e alphafold
```

## Edit a dependency

Edit the corresponding `pyproject.toml`, then run `poetry lock --no-update` in the directory of the `pyproject.toml`.

## Edit colabfold

You can run the tests with

```
pytest tests
```

## Edit alphafold

 * switch to the alphafold folder
 * Make edits to alphafold
 * With the `pip install -e` install, you can directly test them in colabfold
 * Set the last digit of `version=` in setup.py one higher, e.g. to `2.1.1234`
 * git commit and push as usual
 * `git tag v2.1.1234 -m v2.1.1234`, make sure it's the correct number (if you don't have permission, ask one of the team to push the tag)
 * `git push --tags`
 * In colabfold: update the number of `alphafold-colabfold = { version = "` in pyproject.toml, commit and push

# Colab dev setup

While it's generally easier to edit locally, you can also test and develop directly in google colab.

We clone to _directory to avoid python from importing from the directory.

```
%%bash

pip install -U pip
git clone https://github.com/sokrypton/Colabfold _colabfold
pip install _colabfold
# Unholy Hack: Use the files from our cloned git repository instead of installed copy
site_packages=$(python -c 'import site; print(site.getsitepackages()[0])')
rm -r ${site_packages}/colabfold
ln -s $(pwd)/_colabfold/colabfold ${site_packages}/colabfold
```

If you also need to patch alphafold:

```
%%bash

pip uninstall -y alphafold
git clone https://github.com/sokrypton/alphafold _alphafold
pip install -e ./_alphafold
```

After that, restart the runtime.

When you changed a file in the `colabfold` package, you need to reload the modules you were using with `importlib.reload()`, e.g.

```python
import colabfold.batch
import importlib

importlib.reload(colabfold.batch)
```

