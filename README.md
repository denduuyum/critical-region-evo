# critical-region-evo

## How to install

First, we need to obtain networkit source code and apply our patch which added a distance parameter to the centrality calculation. So that it can find a distance based betweenness values.
It is best to do it in python virtual environment.

### Requirement

These requirements are taken from networkit repository. They are required to build networkit from source with our patch.

  1. A modern C++ compiler, e.g.: g++ (>= 6.1), clang++ (>= 3.9) or MSVC (>= 14.13)
     On ubuntu, `build-essential` package suffices.
	 ```shell
	 sudo apt install build-essential
	 ```
  2. Python3 (3.6 or higher) with development libraries. 
     On ubuntu, it can be installed with the following command.
      ```shell
	  sudo apt install python3-dev
	  ```
  3. [Pip][https://pypi.python.org/pypi/pip]
  4. Cmake 3.6 or higher.
  5. Build system: [Make][https://www.gnu.org/software/make/]
  6. Cython version 0.29 or higher (e.g. `pip3 install cython`)
  7. Python3 venv module.
     On ubuntu, it can be install with the following command.
	 ```shell
	 sudo apt install python3-venv
	 ```
  
### Applying patch and building networkit

First let's create a virtual environment for this task.
```shell
python3 -m venv cr_evo
source cr_evo/bin/activate
```
Then, clone CR-evo repository and begin building networkit with the patch.

```shell
git clone https://github.com/denduuyum/critical-region-evo.git
cd critical-region-evo
git submodule update --init --recursive
patch -d ./networkit/ -p1 < networkit.patch
cd networkit
python3 setup.py build_ext -j2
pip3 install -e .
cd ..
```

## How to run

Once networkit library is installed, we can run `cr-evo.py` with the following command.

```shell

```
