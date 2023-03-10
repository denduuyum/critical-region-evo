# critical-region-evo

## How to run

First, we need to obtain networkit source code and apply our patch which added a distance parameter to the centrality calculation. So that it can find a distance based betweenness values.
It is best to do it in python virtual environment.

### Requirement

These requirements are taken from networkit repository. They are required to build networkit from source with our patch.

  1. A modern C++ compiler, e.g.: g++ (>= 6.1), clang++ (>= 3.9) or MSVC (>= 14.13)
  2. Python3 (3.6 or higher) with development libraries. 
     On ubuntu, it can be installed with the following command.
      ```shell
	  sudo apt install python3-dev
	  ```
  3. [Pip][https://pypi.python.org/pypi/pip]
  4. Cmake 3.6 or higher.
  5. Build system: [Make][https://www.gnu.org/software/make/]
  6. Cython version 0.29 or higher (e.g. `pip3 install cython`)
  
### Applying patch and building networkit

First let's create a virtual environment for this task.
```shell
python3 -m venv cr_evo
. cr_evo/bin/activate
```
Then, clone netwokit repository and apply the patch.

```shell
git clone https://github.com/networkit/networkit.git
cd networkit
git checkout 154216e  9ff008476898ca54f1d06105adf25f66e
```
