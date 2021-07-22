# ColabFold template server

ColabFold submits a list of templates to this server, which returns a subset of the PDB70 containing only these hits and the respctive mmCIF files.

## Installation

You need Go (1.12 or later) to compile the server.

```
go build -o template-server
```

Create a `config.json` file and enter the path to the PDB70 database in the `database` field:
```
{
    "verbose": true,
    "server" : {
        "address"    : "127.0.0.1:3002",
        "cors"       : true
    },
    "paths" : {
        "database"    : "/path-to/pdb70"
    }
}
```

The template server requires in addition to the usual PDB70 database files (`pdb70_{a3m,hhm,cs219}.ff{data,index}`) another FFindex database `pdb70_cif.ff{data,index}` containing all mmCIF files. We will provide a script to generate this database at some point.
