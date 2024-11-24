## SA-conf Optimisaiton
**Eliott TEMPEZ - Université Paris Cité**

**M2 Bioinformatique**

### Long project - in collaboration with Leslie REGAD
Protein flexibility is often implied in binding with different partners and is essential for protein function. The analysis of structural variability through available redundant structures of a target, called multiple target conformations (MTC) is one way to explore protein flexibility. *SA-conf* is a tool dedicated to capturing and linking the amino acid and local structure variability of proteins, and analyzing the target structural variability space. It was developped by Leslie REGAD *et al.* in 2017 **(1)**.

The aim of this project is to implement new features in SA-conf, so that it allows for a more in-depth analysis. 

<small>**(1)** Regad, L.; Chéron, J.-B.; Triki, D.; Senac, C.; Flatters, D.; Camproux, A.-C. Exploring the Potential of a Structural Alphabet-Based Tool for Mining Multiple Target Conformations and Target Flexibility Insight. *PLOS ONE* 2017, 12 (8), e0182972. https://doi.org/10.1371/journal.pone.0182972.</small>


### Installation
* First, clone the repository from GitHub:

```
git clone https://github.com/eliott-tempez/M2_SA-conf_optimisation.git
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
