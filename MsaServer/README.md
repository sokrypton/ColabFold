# Setting up your local ColabFold API server

Here you will find two examples of how to setup your own API server on a Linux (or macOS for testing) machine.

## `setup-and-start-local.sh`

The `setup-and-start-local.sh` script will execute most of the steps to get a server running for you.
It will do the following steps:
* check that all required software is installed (`curl`, `aria2c`, `rsync`, `aws`)
* download pinned **MMseqs2** and **mmseqs-server** binaries for your platform (Linux x86\_64/arm64, macOS universal)
* download the databases (UniRef30 and ColabFoldDB, this might take some time)
* download the API server and compile its binary
* start the API server

The script can be called repeatedly to (re)start the server. It avoids unnecessary work and only re-downloads components when the pinned commit changed.

### CPU/GPU and platform detection

* Uncomment the `export GPU=1` line to enable GPU mode (Linux only).
* The script adds the parameters `--paths.colabfold.gpu.gpu 1 --paths.colabfold.gpu.server 1`. See `config.json` for more details.

### Choosing a PDB rsync mirror

At the top of the script you can set the PDB mirror to use (RCSB, PDBe or PDBj).
Uncomment the pair you want. The script exits if no mirror is selected.

### Configuration

Edit `config.json` as needed. Common tweaks:

* `server.address` — change the bind address/port (we recommend putting `nginx` in front for gzip/SSL).
* `local.workers` — number of local job workers.
* Optional GPU block under `paths.colabfold.gpu` lets you pin device IDs per DB when you run multi-GPU.
* A `server.ratelimit` example is included and can be enabled.

### Run

```
./setup-and-start-local.sh
```

If `DEBUG_MINI_DB=1` is set, the server starts with templates disabled and a tiny DB for quick tests.

## Setup a systemd service
To better manage the ColabFold API server, we recommend to setup a systemd service. It will automatically restart on failure and lets you use `journalctl`/`systemctl`.

1. First run `setup-and-start-local.sh` once to get the folder structure and binaries.
2. Adjust the `systemd-example-mmseqs-server.service` example and point it to your paths:
3. Enable and start `./restart-systemd.sh`

## Forcing databases to stay resident in system memory

The ColabFold MSA API server will only achieve response time of few seconds if the search database are held fully within system memory. We use vmtouch (https://github.com/hoytech/vmtouch) to keep the precomputed database index file within system memory. In CPU mode, this is the most expensive part of the MSA API server, as the two default databases (UniRef30+ColabFoldDB) require currently 768GB-1024GB RAM to stay resident in RAM and have enough RAM spare for worker processes.

After installing `vmtouch`, you can execute the following command to make sure that the search databases are not evicted from the system cache:

```
cd databases
sudo vmtouch -f -w -t -l -d -m 1000G *.idx
```

## Using a custom API server

You can pass the server URL to `colabfold_batch` via `--host-url`.
In notebooks, add `host_url=https://yourserver.example.org` to the `run()` call in the *Run Prediction* cell.
