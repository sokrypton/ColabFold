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

## Using a custom API server

You can now pass the server URL to `colabfold_batch`'s `--host-url` parameter. If you want to use the notebook with a custom API server add a `host_url=https://yourserver.example.org` parameter to the `run()` call in the *Run Prediction* cell.

Templates are still requested from our server (if the `--templates` flag is used). We will improve templates in a future release.
