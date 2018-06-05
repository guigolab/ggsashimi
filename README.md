# ggsashimi

[![Build Status](https://travis-ci.org/guigolab/ggsashimi.svg?branch=master)](https://travis-ci.org/guigolab/ggsashimi)

Command-line tool for the visualization of splicing events across multiple samples

## Installation

The `ggsashimi` script can be directly downloaded from this repository:
```
wget https://raw.githubusercontent.com/guigolab/ggsashimi/master/sashimi-plot.py
```
Change the execution permissions:
```
chmod u+x sashimi-plot.py
```
Provided all dependencies are already installed (see below), you can directly execute the script:
```
./sashimi-plot.py --help
```
To download the entire repository, which includes the dockerfile and example files:
```
git clone https://github.com/guigolab/ggsashimi.git
```

To avoid dependencies issues, the script is also available through a docker image.

### Download docker image

A public `ggsashimi` Docker image is available in the [Docker Hub](https://hub.docker.com/r/guigolab/ggsashimi/) and can be downloaded as follows:
```
docker pull guigolab/ggsashimi
```
Alternatively, we provide the Dockerfile if you want to build your local docker image.


### Build docker image
After downloading the repository, move inside the repository folder:
```
cd ggsashimi
```
To build the docker image run the following command:
```
docker build -f docker/Dockerfile -t guigolab/ggsashimi .
```
This can take several minutes. 


### Use docker image
Once the image is downloaded or built, to execute ggsashimi with docker:
```
docker run guigolab/ggsashimi --help
```
Because the image is used in a docker container which has its own file system, to use the program with local files, a host data volume needs to be mounted.

As an example, you can run this command from the main repository folder:
```
docker run -w $PWD -v $PWD:$PWD guigolab/ggsashimi -b examples/input_bams.tsv -c chr10:27040584-27048100
```
The '-w' option sets the working directory inside the container to the current directory.
The '-v' option mounts the current working directory and all child folders inside the container to the same path (host_path:container_path).
If your files are in another folder, for example the annotation file is stored in a different folder then the one containing the bam file, you can mount extra folders like this:
```
f="$DIR/annotation.gtf"
docker run -w $PWD -v $PWD:$PWD -v $DIR:$DIR guigolab/ggsashimi -b examples/input_bams.tsv -c chr10:27040584-27048100 -g $f
```
You can even mount a single file:
```
docker run -w $PWD -v $PWD:$PWD -v $f:$f guigolab/ggsashimi -b examples/input_bams.tsv -c chr10:27040584-27048100 -g $f
```

## Dependencies

In order to run `ggsashimi` the following software components and packages are required:

- python (2.7 or 3)
- samtools (>=1.3)
- R (>=3.3)
  - ggplot2 (>=2.2.1)
  - data.table (>=1.10.4)
  - gridExtra (>=2.2.1)   

Additional required R packages `grid` and `gtable` should be automatically installed when installing R and `ggplot2`, respectively.

## Usage
Execute the script with `--help` option for a complete list of options.  
An example of usage can be found at `examples/example_run.sh`
