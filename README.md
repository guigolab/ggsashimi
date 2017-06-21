# ggsashimi
Command-line tool for the visualization of splicing events across multiple samples

## Installation
The ggsashimi script can be directly downloaded from this repository:
```
wget https://raw.githubusercontent.com/guigolab/ggsashimi/master/sashimi-plot.py
```
Change the execution permissions:
```
chmod u+x sashimi-plot.py
```
Provided all dependencies are already installed, you can directly execute the script:
```
./sashimi-plot.py --help
```
To download the entire repository, which includes the dockerfile and example files:
```
git clone https://github.com/guigolab/ggsashimi.git
```
To avoid dependecies issues, the script is also available through a docker image.

## Docker image

A public `ggsashimi` Docker image is available in the [Docker Hub](https://hub.docker.com/r/guigolab/ggsashimi/) and can be downloaded as follows:

```
docker pull guigolab/ggsashimi
```
To execute ggsashimi with docker:
```
docker run guigolab/ggsashimi --help
```
Alternatively, we provide the Dockerfile if you want to build your local docker image.

## Build docker image

After downloading the repository, move inside the repository folder:
```
git clone https://github.com/guigolab/ggsashimi.git
cd ggsashimi
```
To build the docker image run the following command:
```
docker build -f docker/Dockerfile -t guigolab/ggsashimi .
```
This can take several minutes. After the image is built, ggsashimi can be executed as follows:
```
docker run guigolab/ggsashimi --help
```
