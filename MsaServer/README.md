# Setting up your local ColabFold API server

Here you will find two examples of how to setup your own API server on a Linux machine.

## `setup-and-start-local.sh`

The `setup-and-start-local.sh` script will execute most of the steps to get a server running for you.
It will do the following steps:
* check that all required software is installed (go, git, aria2c, curl)
* download a specific MMseqs2 version that we tested to work with ColabFold
* download the databases (UniRef30 and ColabFoldDB, this might take some time)
* download the API server and compile its binary
* start the API server

The script can be called repeatedly to start the server. It will avoid doing any unnecessary setup work.

### Tweaking `config.json`
You can tweak the provided example `config.json` file. Two values you will likely want to change are the:
* `server.address` field to specify a custom port. We recommend putting a `nginx` server infront of the ColabFold API server to deal with gzip compression, SSL etc.
* `local.workers` field to specify how many job workers are allowed to run in parallel.

## Setup a systemd service
To better manage the ColabFold API server, we recommend to setup a systemd service. It will automatically restart the server in case of failure and allow managing the server with common operating system tools (like `journalctl`, `systemctl`, etc.).

You can first execute the `setup-and-start-local.sh` script to get a working server. Then tweak the `systemd-example-mmseqs-server.service` file and adjust the directories to the MSA server binary/directories.

Afterwards copy the service file to `/etc/systemd/system/mmseqs-server.service` (this might vary by distribution).

Then call:
```
sudo systemctl daemon-reload
sudo systemctl start mmseqs-server.service
```

The `restart-systemd.sh` script contains an example how to stop the server, clear the job cache and start it again.

## Forcing databases to stay resident in system memory

The ColabFold MSA API server will only achieve response time of few seconds if the search database are held fully within system memory. We use vmtouch (https://github.com/hoytech/vmtouch) to keep the precomputed database index file within system memory. This is the most expensive part of the MSA API server, as the two default databases (UniRef30+ColabFoldDB), require currently 768GB-1024GB RAM to stay resident in RAM and have enough RAM spare for worker processes. If you are only running batch searches or are using the command line tool with our API server, system requirements are much much lower.

After installing `vmtouch`, you can execute the following command to make sure that the search databases are not evicted from the system cache:

```
cd databases
sudo vmtouch -f -w -t -l -d -m 1000G *.idx
```

This assumes that precomputed database index was created without splits. Check that there are no `uniref30_2103_db.idx.{0,1,...}` or `colabfold_envdb_202108_db.idx.{0,1,...}` files in the databases folder. If these files are there, you should recreate the precomputed database indices with the following command:

```
cd databases
rm uniref30_2103_db.idx* colabfold_envdb_202108_db.idx*
mmseqs createindex uniref30_2103_db tmp --remove-tmp-files 1 --split 1
mmseqs createindex colabfold_envdb_202108_db tmp --remove-tmp-files 1 --split 1
```

## Using a custom API server

You can now pass the server URL to `colabfold_batch`'s `--host-url` parameter. If you want to use the notebook with a custom API server add a `host_url=https://yourserver.example.org` parameter to the `run()` call in the *Run Prediction* cell.

Templates are still requested from our server (if the `--templates` flag is used). We will improve templates in a future release.
