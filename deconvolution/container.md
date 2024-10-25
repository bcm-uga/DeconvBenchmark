# Using an Apptainer container

In order to have a fixed environment, one can use an Apptainer container, an HPC friendly alternative to Docker containers.

[Full documentation](https://docs.sylabs.io/guides/main/user-guide/index.html)

## Building a container

Building a container can be done through the use of a **definition** file (.def) defining the operating system of the container, the commands to be run at build time, the ones to execute at run time, etc...
In our case (see [container2.def](container2.def) file) the OS is Ubuntu/Jammy onto which we install **R** and all the required packages.  
Once the definition file is readyn the **container** can be built into a .sif file as simply as:

```shell
sudo apptainer build container2.sif container2.def
```

In the cluster you can run
```shell
apptainer build --remote container_test.sif container2.def
```

Due to access restrictions to GitHub API it can be tricky to install a lot of R packages via *devtools::install_github* during the build, causing it to fail.  
In that case it might be preferable to:

- build an editable container without installing R packages from GitHub (by commenting the corresponding line of [container2.def](container2.def)),
```shell
sudo apptainer build --sandbox container2 container2.def
```

- *connect* to the container and install the GitHub hosted packages one by one (via R),
```shell
sudo apptainer shell --writable container2
# inside the container
R -e 'devtools::install_github(...)'
# once done
exit
```

- build the final container from the editable one.
```shell
sudo apptainer build container2.sif container2
```

### Updating the container

It might be interesting to use an editable container to quickly build a new .sif container after having installed new dependencies or packages into the editable one, without the need to rebuild everything from scratch.

## Running experiments inside the container

Commands can now be executed inside the container via:

```shell
apptainer exec container2.sif some command
```

If we want to execute several commands it is easier to create an executable shell script to be run in the container:

```shell
# insert some commands
cat > run.sh
#!/bin/bash
cd path/to/folder
Rscript script.R

# make it executable
chmod +x run.sh

# execute run.sh inside container2.sif
apptainer exec --bind /bettik:/bettik container2.sif ./run.sh
```

## Using the container on the cluster
The .sif file can be copied to the cluster and then used in a submission script (or in an interactive session).