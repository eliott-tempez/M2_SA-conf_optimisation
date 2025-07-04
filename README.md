This is the git for the details of my master's project around improving SA-conf. For the official SA-conf github, click [here](https://github.com/eliott-tempez/SA-conf)

## SA-conf Optimisaiton
**Eliott TEMPEZ - Université Paris Cité**

**M2 Bioinformatique**

### Long project - in collaboration with Leslie REGAD
Protein flexibility is often implied in binding with different partners and is essential for protein function. The analysis of structural variability through available redundant structures of a target, called multiple target conformations (MTC) is one way to explore protein flexibility. *SA-conf* is a tool dedicated to capturing and linking the amino acid and local structure variability of proteins, and analyzing the target structural variability space. It was developped by Leslie REGAD *et al.* in 2017 **(1)**.

The aim of this project is to implement new features in SA-conf, so that it allows for a more in-depth analysis. 


### Installation
* First, clone the repository from GitHub (only the main branch):

```
git clone --branch main https://github.com/eliott-tempez/M2_SA-conf_optimisation.git
cd M2_SA-conf_optimisation
```

* Create and install the Conda environment from the env.yaml file:

```
conda env create -f env.yaml
```

* Activate the environment:

```
conda activate SA-conf
```

* You can now use the program. The guidelines are available in the saconf-tutorial pdf.


### New additions (full report : `report.pdf`)
The new additions include:
- Switching from Python 2 to Python 3
- Distinguishing unmodelled residues from absent residues bescause of DNA deletions in the graphs
- Taking uncertain coordinates into account via the RSRZ and adding it to the main graphs
- Taking the B-factor into account and adding it to the *neq* graph
- Clustering the proteins according to their structural letters sequence

The examples for the differences between the old program and the new are in the `/example_outputs` folder.


### Case study
The program and its new additions were tested on the Bcl-Xl protein. You can find the whole report in the `report.pdf` file. The input data is in the `/data` folder, and the results in the `/results` folder.


### References
**(1)** Regad, L.; Chéron, J.-B.; Triki, D.; Senac, C.; Flatters, D.; Camproux, A.-C. Exploring the Potential of a Structural Alphabet-Based Tool for Mining Multiple Target Conformations and Target Flexibility Insight. *PLOS ONE* 2017, 12 (8), e0182972. https://doi.org/10.1371/journal.pone.0182972.
