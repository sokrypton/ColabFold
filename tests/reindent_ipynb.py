import json

for notebook in ["AlphaFold2.ipynb", "batch/AlphaFold2_batch.ipynb"]:
    data = json.loads(open(notebook).read())
    open(notebook, "w").write(json.dumps(data, indent=2))
