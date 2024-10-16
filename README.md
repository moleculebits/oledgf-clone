### Disclaimer
This repository is a clone of the institutional GitLab repository where the OLEDgf project is being developed. It is intended for demonstration purposes ONLY and should not be forked. The first full release to GitHub is expected soon.
### oledgf (OLED green's function)
Is a high-performance library for simulating and fitting the energy emission behavior of light-emitting material stacks. This library is under active development by T. Marcato and M. Nouman at ETH Zurich.

### oledgf Installation

## Creating Development environment

To ensure a standardized environment for development the project has to be built within a Docker container.

In order to do so it requires:

- Docker
- An X11 server to forward the application UI 

### Install Docker

- Download and install Docker Desktop ([Windows](https://docs.docker.com/desktop/install/windows-install/) or [Mac](https://docs.docker.com/desktop/install/mac-install/)).

- Open Docker Desktop to start the Docker daemon.

- (On Windows) If the Docker engine cannot run make sure that you have enabled Hardware Virtualization (see [this guide](https://docs.docker.com/desktop/troubleshoot/topics/#virtualization))

### Install X11 Server

#### Windows

- Download and install [XServer](https://sourceforge.net/projects/vcxsrv/)

#### Mac

- Download and install [XQuartz](https://www.xquartz.org/)
- In XQuartz settings, under Security tab, check "Allow connections from network clients"
- Restart XQuartz
- In a terminal on the host, run `xhost +localhost`

### Create Dev Container Image

First clone the main repo and enter root dir of project
```
git clone https://github.com/moleculebits/oledgf
cd oledgf
```

Open terminal and run
```
docker build -t test:0.0.1 -f Dockerfile .
```
to build the docker image `test` that will contain the dev environment.

### Run the container

In the terminal

```
docker run -it --rm --name=test --mount type=bind,source=${PWD},target=/src test:0.0.1 bash
```
If you are running docker from Git Bash on Windows you should run 

```
winpty docker run -it --rm --name=test --mount type=bind,source=${PWD},target=/src test:0.0.1 bash
```

You will be now in the bash terminal of the container (based on Ubuntu).

Navigate to the `src` directory which is mapped to the project root.

```
cd ./src
```

### Build the project

Build the project

```
cmake -S . -B ./build && cd ./build
```

and compile it

```
make
```

# Build Details

The build system is based on the [cmake_template](https://github.com/cpp-best-practices/cmake_template) by Jason Turner.

For now it includes

- High level warnings and Warnings as Errors (`OFF` by default)
- [clang-format](https://clang.llvm.org/docs/ClangFormat.html)
- [CPM](https://github.com/cpm-cmake/CPM.cmake) for external dependencies
- Sanitizers
- Static check with CPPcheck

## External dependencies

The project relies on the following dependencies

- [fmt](https://github.com/fmtlib/fmt): For formatting
- [Eigen](https://gitlab.com/libeigen/eigen): For Matrices and linear algebra
- [Matplot++](https://github.com/alandefreitas/matplotplusplus): For plotting (gnuplot backend)
