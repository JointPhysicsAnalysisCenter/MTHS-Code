# MTHS Code

Repository contains the code for exercise sessions of Modern Techniques in Hadron Spectroscopy school.

### Clone the Code Repository 

```bash
cd ~/Documents
git clone https://github.com/JointPhysicsAnalysisCenter/MTHS-code MTHS-code
cd MTHS-code
```

## Python installation via Conda

1. ****Install `cond**a`**

Make sure you have `conda` installed. The easiest way for homebrew users is:
```bash
brew install --cask anaconda
```
For Linux, follow the instructions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).

Create a new conda environment and install the required packages:
```bash
conda create --name lattice-qcd python=3.xxx
conda activate lattice-qcd
conda install numpy scipy matplotlib iminuit jupyterlab
```

2. **Start Jupyter Lab and Open the Notebook**

Launch Jupyter Lab from the terminal:
```bash
jupyter lab
```
In your web browser, navigate to the Jupyter Lab interface. Open the notebook file inside the folder you just downloaded.

Open the `.ipynb` and Follow the Instructions in the NotebookNote that the first two cells just load packages and fixed variables, as described in the notebook.

### Jupyter Notebook Cheat Sheet

For basic commands and shortcuts in Jupyter Notebook, refer to the [Jupyter Notebook Cheat Sheet](https://jupyter-notebook.readthedocs.io/en/stable/notebook.html).



## Julia and Pluto installation

1. **Download and install Julia** from the [official website](https://julialang.org/downloads/).

2. Open Julia and **install `Pluto.jl`**:

```julia
] add Pluto
```

3. **Start Pluto and Open the Notebook**

Launch Pluto from the Julia REPL:
```julia
using Pluto
Pluto.run()
```
In your web browser, navigate to the Pluto interface.
Open the notebook files.


# README

## Practical Exercise: Correlated Fits to Monte Carlo Sampled Two-Point Correlators

In this practical exercise, we perform correlated fits to Monte Carlo sampled two-point correlators, calculated from Wilson lines for a $L=32$ lattice with periodic boundary conditions.

### Notes and Tips

- For Julia notes, visit [JuliaNotes.jl](https://m3g.github.io/JuliaNotes.jl/stable/).
- For `Pluto.jl` tips, check out the [Pluto.jl GitHub repository](https://github.com/fonsp/Pluto.jl?ref=juliafordatascience.com).


For more detailed information and additional resources, refer to the code repository.

## C++ and openQCD installation




