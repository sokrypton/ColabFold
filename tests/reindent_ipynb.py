import json
from pathlib import Path

if __name__ == "__main__":
    for notebook in Path(".").rglob("*.ipynb"):
        print(notebook)
        data = json.loads(open(notebook).read())
        open(notebook, "w").write(json.dumps(data, indent=2, ensure_ascii=False))
