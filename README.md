# Scaling Package Queries to a Billion Tuples via Hierarchical Partitioning and Customized Optimization

## Abstract

A package query returns a package — a multiset of tuples — that maximizes or minimizes a linear objective function subject to linear constraints, thereby enabling in-database decision support. Prior work has established the equivalence of package queries to Integer Linear Programs (ILPs) and developed the Sketch Refine algorithm for package query processing. While this algorithm was an important first step toward supporting prescriptive analytics scalably inside a relational database, it struggles when the data size grows beyond a few hundred million tuples or when the constraints become very tight. In this paper, we present Progressive Shading, a novel algorithm for processing package queries that can scale efficiently to billions of tuples and gracefully handle tight constraints. Progressive Shading solves a sequence of optimization problems over a hierarchy of relations, each resulting from an ever-finer partitioning of the original tuples into homogeneous groups until the original relation is obtained. This strategy avoids the premature discarding of high-quality tuples that can occur with Sketch Refine. Our novel partitioning scheme, Dynamic Low Variance, can handle very large relations with multiple attributes and can dynamically adapt to both concentrated and spread-out sets of attribute values, provably outperforming traditional partitioning schemes such as kd-tree. We further optimize our system by replacing our off-the-shelf optimization software with customized ILP and LP solvers, called Dual Reducer and Parallel Dual Simplex respectively, that are highly accurate and orders of magnitude faster.

## Installation and Setup

### Install Python packages

This code needs Python 3.8 or higher.
```bash
pip3 install -r requirements.txt
```

### Setup Postgres server and hardware information

This code needs Postgres v14.7 or higher. Configure the file *config.txt* with the appropriate information 
 
Parameters for Postgres: hostname, port, username, database, password, schema
  
Parameters for Hardware information:
  
* physical_core: Number of cores in CPU
  
* main_memory_size: Main memory allowed in GB
  

### Data generation

```bash
cd resource
python3 ssds.py
python3 tpch.py
```

### Setup C++ and Gurobi

#### Install CMake

Download CMake

```bash
wget https://github.com/Kitware/CMake/releases/download/v3.23.0-rc1/cmake-3.23.0-rc1-linux-x86_64.tar.gz
tar -xvf cmake-3.23.0-rc1-linux-x86_64.tar.gz
```

Add CMake path

```bash
vim .bashrc
insert into .bashrc the line "export PATH=$<CMakeFolder>/cmake-3.23.0-rc1-linux-x86_64/:$PATH"
source .bashrc
```

#### Install Conan

```bash
pip3 install conan
cd PackageQuery/build
conan install .
```

#### Install Gurobi

This code needs Gurobiv9.5.2 or higher. 

1. Request an academic license from [Gurobi](https://www.gurobi.com/academia/academic-program-and-licenses/).
2. Follow the instructions to download and unzip Gurobi as well as set up the license.
3. Edit the file *cmake/FindGUROBI.cmake*:
   * Change the line
   ```bash
   set(GUROBI_HOME "/home/alm818/downloads/gurobi950/linux64")
   ```
     to the Gurobi's home in your system
   * Change the line
   ```bash
   find_library(GUROBI_LIBRARY
    NAMES gurobi gurobi95
    HINTS ${GUROBI_DIR} ${GUROBI_HOME}
    PATH_SUFFIXES lib)
   ```
     to the correct Gurobi's version in your system

## Result Reproduction

Build files

```bash
  cd build
  cmake -DCMAKE_BUILD_TYPE=Release .. 
```

Compile files and run

```bash
  cd build
  cmake --build . --config Release
  bin/main
```

This will generate 3 files in *plots* folder: A3.csv, A4_tpch.csv and A4_ssds.csv

You can draw the related graph in the paper by running Experiment.ipynb in the *plots* folder

