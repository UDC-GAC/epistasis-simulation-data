# Epistasis simulation data

This repository serves three purposes:

1. Document the simulation design process used in our work __*[Citation missing,
   paper under revision]*__.
2. Publish our developed data sets, allowing other researchers to use the same
   data.
3. Ensure that our results are reproducible by any other scientist.


## Simulation design

The objective of our work is to compare the performance of several epistasis
detection methods, both in terms of detection power and false positives. To do
that, a collection of data sets was created with the intention of providing a
wide variety of population characteristics resembling real human traits, on
which the different methods will be evaluated.

All our simulation was carried out using [GAMETES][][1], the most used epistasis
simulation package in recent literature. Our simulation includes both epistasis
with marginal effects and with no marginal effects. Simulation in GAMETES can be
summarized in three steps:

1. Calculate a penetrance table describing the epistatic interaction to be
   included in the data.
2. Simulate population samples of the interacting SNPs following its penetrance
   table.
3. Simulate population samples of the rest of the non-interacting SNPs.

Our criteria was to create third and fourth order interactions, with MAFs of the
interacting SNPs of 0.10, 0.25 and 0.40, resulting in a trait with
heritabilities of 0.10, 0.25, 0.50 and 0.80, and prevalences above 1E-06.

### Penetrance tables

GAMETES can create penetrance tables describing interactions following no
epistasis model and showing no marginal effects. Therefore, to generate
penetrance tables following a specific epistasis model, [Toxo][][2] was used.
Toxo is a MATLAB library capable of calculating penetrance tables of any
bivariate epistasis model under certain conditions. In this study we decided to
use the additive and threshold models proposed my Machini *et al.* in [3] for
third and fourth order epistasis.

The following MATLAB code snippet shows how the different penetrance tables with
marginal effects were obtained using Toxo. The directory
[epistasis_models/](epistasis_models/) contains the used Marchini's models in
CSV format, as required by Toxo.

```matlab
addpath('<path to Toxo>/src/');
m = toxo.Model('epistasis_models/<model>');
pt = m.find_max_prevalence([<MAFs>], <heritability>);
pt.write('<path to penetrance table output>',
         toxo.PTable.format_gametes,
         <MAFs>);
```

Following this procedure, the following penetrance tables were obtained:

| Model     | Order | MAF  | h²   | P(D)     | File                |
|-----------|-------|------|------|----------|---------------------|
| Additive  | 3     | 0.10 | 0.10 | 0.000012 | [Link][additive1]   |
| Additive  | 3     | 0.10 | 0.25 | 0.000004 | [Link][additive2]   |
| Additive  | 3     | 0.10 | 0.50 | 0.000002 | [Link][additive3]   |
| Additive  | 3     | 0.10 | 0.80 | 0.000001 | [Link][additive4]   |
| Additive  | 3     | 0.25 | 0.10 | 0.005370 | [Link][additive5]   |
| Additive  | 3     | 0.25 | 0.25 | 0.001153 | [Link][additive6]   |
| Additive  | 3     | 0.25 | 0.50 | 0.000504 | [Link][additive7]   |
| Additive  | 3     | 0.25 | 0.80 | 0.000306 | [Link][additive8]   |
| Additive  | 3     | 0.40 | 0.10 | 0.254558 | [Link][additive9]   |
| Additive  | 3     | 0.40 | 0.25 | 0.022186 | [Link][additive10]  |
| Additive  | 3     | 0.40 | 0.50 | 0.008545 | [Link][additive11]  |
| Additive  | 3     | 0.40 | 0.80 | 0.005091 | [Link][additive12]  |
| Additive  | 4     | 0.10 | 0.10 | < 1E-6   | Unavailable         |
| Additive  | 4     | 0.10 | 0.25 | < 1E-6   | Unavailable         |
| Additive  | 4     | 0.10 | 0.50 | < 1E-6   | Unavailable         |
| Additive  | 4     | 0.10 | 0.80 | < 1E-6   | Unavailable         |
| Additive  | 4     | 0.25 | 0.10 | 0.000234 | [Link][additive13]  |
| Additive  | 4     | 0.25 | 0.25 | 0.000068 | [Link][additive14]  |
| Additive  | 4     | 0.25 | 0.50 | 0.000031 | [Link][additive15]  |
| Additive  | 4     | 0.25 | 0.80 | 0.000019 | [Link][additive16]  |
| Additive  | 4     | 0.40 | 0.10 | 0.036282 | [Link][additive17]  |
| Additive  | 4     | 0.40 | 0.25 | 0.003383 | [Link][additive18]  |
| Additive  | 4     | 0.40 | 0.50 | 0.001374 | [Link][additive19]  |
| Additive  | 4     | 0.40 | 0.80 | 0.000822 | [Link][additive20]  |
| Threshold | 3     | 0.10 | 0.10 | 0.064602 | [Link][threshold1]  |
| Threshold | 3     | 0.10 | 0.25 | 0.025561 | [Link][threshold2]  |
| Threshold | 3     | 0.10 | 0.50 | 0.013270 | [Link][threshold3]  |
| Threshold | 3     | 0.10 | 0.80 | 0.008417 | [Link][threshold4]  |
| Threshold | 3     | 0.25 | 0.10 | 0.477516 | [Link][threshold5]  |
| Threshold | 3     | 0.25 | 0.25 | 0.267707 | [Link][threshold6]  |
| Threshold | 3     | 0.25 | 0.50 | 0.154539 | [Link][threshold7]  |
| Threshold | 3     | 0.25 | 0.80 | 0.102529 | [Link][threshold8]  |
| Threshold | 3     | 0.40 | 0.10 | 0.780354 | [Link][threshold9]  |
| Threshold | 3     | 0.40 | 0.25 | 0.586967 | [Link][threshold10] |
| Threshold | 3     | 0.40 | 0.50 | 0.415395 | [Link][threshold11] |
| Threshold | 3     | 0.40 | 0.80 | 0.307526 | [Link][threshold12] |
| Threshold | 4     | 0.10 | 0.10 | 0.012563 | [Link][threshold13] |
| Threshold | 4     | 0.10 | 0.25 | 0.005140 | [Link][threshold14] |
| Threshold | 4     | 0.10 | 0.50 | 0.002590 | [Link][threshold15] |
| Threshold | 4     | 0.10 | 0.80 | 0.001623 | [Link][threshold16] |
| Threshold | 4     | 0.25 | 0.10 | 0.275518 | [Link][threshold17] |
| Threshold | 4     | 0.25 | 0.25 | 0.132034 | [Link][threshold18] |
| Threshold | 4     | 0.25 | 0.50 | 0.070683 | [Link][threshold19] |
| Threshold | 4     | 0.25 | 0.80 | 0.041819 | [Link][threshold20] |
| Threshold | 4     | 0.40 | 0.10 | 0.668428 | [Link][threshold21] |
| Threshold | 4     | 0.40 | 0.25 | 0.446405 | [Link][threshold22] |
| Threshold | 4     | 0.40 | 0.50 | 0.287337 | [Link][threshold23] |
| Threshold | 4     | 0.40 | 0.80 | 0.201273 | [Link][threshold24] |

Four penetrance tables from the fourth order additive model were discarded due
to resulting in very low prevalence values, since prevalences that low are
unrealistic for real human populations.

Penetrance tables with no epistatic model can be obtained directly from GAMETES,
by running the following shell command:

```sh
java -jar <path to GAMETES> \
    -M "-h <heritability> <MAFs> -o <output model file>" \
    -q 1 \
    -p 100 \
    -t 100000 \
    -r <seed>
```

Unlike Toxo, GAMETES implements a stochastic algorithm and thus penetrance
tables can vary depending on the seed used by the pseudorandom number generator.
Therefore, we also picked seeds at random to guarantee that results are
reproducible:

| Order | MAF  | h²   | P(D)   | GAMETES seed | Link                |
|-------|------|------|--------|--------------|---------------------|
| 3     | 0.10 | 0.10 | -      |              | Unavailable         |
| 3     | 0.10 | 0.25 | -      |              | Unavailable         |
| 3     | 0.10 | 0.50 | -      |              | Unavailable         |
| 3     | 0.10 | 0.80 | -      |              | Unavailable         |
| 3     | 0.25 | 0.10 | 0.5860 | -1643481676  | [Link][nme1]        |
| 3     | 0.25 | 0.25 | 0.4923 | 1764873474   | [Link][nme2]        |
| 3     | 0.25 | 0.50 | 0.4223 | 1893932570   | [Link][nme3]        |
| 3     | 0.25 | 0.80 | -      |              | Unavailable         |
| 3     | 0.40 | 0.10 | 0.5163 | -1568883956  | [Link][nme4]        |
| 3     | 0.40 | 0.25 | 0.5644 | 2089582692   | [Link][nme5]        |
| 3     | 0.40 | 0.50 | 0.5019 | 1343608856   | [Link][nme6]        |
| 3     | 0.40 | 0.80 | 0.4970 | 1329415446   | [Link][nme7]        |
| 4     | 0.10 | 0.10 | -      |              | Unavailable         |
| 4     | 0.10 | 0.25 | -      |              | Unavailable         |
| 4     | 0.10 | 0.50 | -      |              | Unavailable         |
| 4     | 0.10 | 0.80 | -      |              | Unavailable         |
| 4     | 0.25 | 0.10 | 0.4201 | -600481346   | [Link][nme8]        |
| 4     | 0.25 | 0.25 | 0.5910 | -965873914   | [Link][nme9]        |
| 4     | 0.25 | 0.50 | -      |              | Unavailable         |
| 4     | 0.25 | 0.80 | -      |              | Unavailable         |
| 4     | 0.40 | 0.10 | 0.4356 | 913749367    | [Link][nme10]       |
| 4     | 0.40 | 0.25 | 0.4720 | 203584226    | [Link][nme11]       |
| 4     | 0.40 | 0.50 | -      |              | Unavailable         |
| 4     | 0.40 | 0.80 | -      |              | Unavailable         |

Many of the penetrance tables could not be obtained. The larger the interaction,
the lower the MAFs and the higher the heritability are, the more unlikely that
GAMETES is to converge on a solution.

### Data generation

Using previous penetrance tables, 100 data sets were simulated for each
penetrance table containing 500 SNPs (including the interacting SNPs) of 2000
individuals (1000 cases and 1000 controls), with MAFs of non-interacting loci
uniformly sampled from the interval [0.05, 0.5]. This can be achieved by running
in a shell, for each penetrance table:

```sh
java -jar <path to GAMETES> \
    -i <path to penetrance table> \
    -D "-n <lower MAF>
        -x <upper MAF>
        -a <number of SNPs>
        -s <number of cases>
        -w <number of controls>
        -r <number of datasets>
        -o <output folder>"
```

The following table lists all epistatic data sets used during the evaluation:

| Model      | Marginal Effects | Order | MAF  | h²   | Data set               |
|------------|------------------|-------|------|------|------------------------|
| Additive   | Yes              | 3     | 0.10 | 0.10 | [Folder][dataset1]     |
| Additive   | Yes              | 3     | 0.10 | 0.25 | [Folder][dataset2]     |
| Additive   | Yes              | 3     | 0.10 | 0.50 | [Folder][dataset3]     |
| Additive   | Yes              | 3     | 0.10 | 0.80 | [Folder][dataset4]     |
| Additive   | Yes              | 3     | 0.25 | 0.10 | [Folder][dataset5]     |
| Additive   | Yes              | 3     | 0.25 | 0.25 | [Folder][dataset6]     |
| Additive   | Yes              | 3     | 0.25 | 0.50 | [Folder][dataset7]     |
| Additive   | Yes              | 3     | 0.25 | 0.80 | [Folder][dataset8]     |
| Additive   | Yes              | 3     | 0.40 | 0.10 | [Folder][dataset9]     |
| Additive   | Yes              | 3     | 0.40 | 0.25 | [Folder][dataset10]    |
| Additive   | Yes              | 3     | 0.40 | 0.50 | [Folder][dataset11]    |
| Additive   | Yes              | 3     | 0.40 | 0.80 | [Folder][dataset12]    |
| Additive   | Yes              | 4     | 0.25 | 0.10 | [Folder][dataset13]    |
| Additive   | Yes              | 4     | 0.25 | 0.25 | [Folder][dataset14]    |
| Additive   | Yes              | 4     | 0.25 | 0.50 | [Folder][dataset15]    |
| Additive   | Yes              | 4     | 0.25 | 0.80 | [Folder][dataset16]    |
| Additive   | Yes              | 4     | 0.40 | 0.10 | [Folder][dataset17]    |
| Additive   | Yes              | 4     | 0.40 | 0.25 | [Folder][dataset18]    |
| Additive   | Yes              | 4     | 0.40 | 0.50 | [Folder][dataset19]    |
| Additive   | Yes              | 4     | 0.40 | 0.80 | [Folder][dataset20]    |
| Threshold  | Yes              | 3     | 0.10 | 0.10 | [Folder][dataset21]    |
| Threshold  | Yes              | 3     | 0.10 | 0.25 | [Folder][dataset22]    |
| Threshold  | Yes              | 3     | 0.10 | 0.50 | [Folder][dataset23]    |
| Threshold  | Yes              | 3     | 0.10 | 0.80 | [Folder][dataset24]    |
| Threshold  | Yes              | 3     | 0.25 | 0.10 | [Folder][dataset25]    |
| Threshold  | Yes              | 3     | 0.25 | 0.25 | [Folder][dataset26]    |
| Threshold  | Yes              | 3     | 0.25 | 0.50 | [Folder][dataset27]    |
| Threshold  | Yes              | 3     | 0.25 | 0.80 | [Folder][dataset28]    |
| Threshold  | Yes              | 3     | 0.40 | 0.10 | [Folder][dataset29]    |
| Threshold  | Yes              | 3     | 0.40 | 0.25 | [Folder][dataset30]    |
| Threshold  | Yes              | 3     | 0.40 | 0.50 | [Folder][dataset31]    |
| Threshold  | Yes              | 3     | 0.40 | 0.80 | [Folder][dataset32]    |
| Threshold  | Yes              | 4     | 0.10 | 0.10 | [Folder][dataset33]    |
| Threshold  | Yes              | 4     | 0.10 | 0.25 | [Folder][dataset34]    |
| Threshold  | Yes              | 4     | 0.10 | 0.50 | [Folder][dataset35]    |
| Threshold  | Yes              | 4     | 0.10 | 0.80 | [Folder][dataset36]    |
| Threshold  | Yes              | 4     | 0.25 | 0.10 | [Folder][dataset37]    |
| Threshold  | Yes              | 4     | 0.25 | 0.25 | [Folder][dataset38]    |
| Threshold  | Yes              | 4     | 0.25 | 0.50 | [Folder][dataset39]    |
| Threshold  | Yes              | 4     | 0.25 | 0.80 | [Folder][dataset40]    |
| Threshold  | Yes              | 4     | 0.40 | 0.10 | [Folder][dataset41]    |
| Threshold  | Yes              | 4     | 0.40 | 0.25 | [Folder][dataset42]    |
| Threshold  | Yes              | 4     | 0.40 | 0.50 | [Folder][dataset43]    |
| Threshold  | Yes              | 4     | 0.40 | 0.80 | [Folder][dataset44]    |
| None       | No               | 3     | 0.25 | 0.10 | [Folder][dataset45]    |
| None       | No               | 3     | 0.25 | 0.25 | [Folder][dataset46]    |
| None       | No               | 3     | 0.25 | 0.50 | [Folder][dataset47]    |
| None       | No               | 3     | 0.40 | 0.10 | [Folder][dataset48]    |
| None       | No               | 3     | 0.40 | 0.25 | [Folder][dataset49]    |
| None       | No               | 3     | 0.40 | 0.50 | [Folder][dataset50]    |
| None       | No               | 3     | 0.40 | 0.80 | [Folder][dataset51]    |
| None       | No               | 4     | 0.25 | 0.10 | [Folder][dataset52]    |
| None       | No               | 4     | 0.25 | 0.25 | [Folder][dataset53]    |
| None       | No               | 4     | 0.40 | 0.10 | [Folder][dataset54]    |
| None       | No               | 4     | 0.40 | 0.25 | [Folder][dataset55]    |


## Downloading

Performing a `git clone` is the preferred method of downloading this repository,
although, due to its large size, it might take a while to complete. If you are
interested in only downloading a specific folder, you may want to use
subversion instead:

    svn checkout https://github.com/UDC-GAC/epistasis-simulation-data/trunk/<folder>

For example, if you are only interested in epistasis data with no marginal
effects, you can download only those datasets by running `svn checkout
https://github.com/UDC-GAC/epistasis-simulation-data/trunk/datasets/epistasis/no_model`.


## References

[1] Urbanowicz, R.J., Kiralis, J., Sinnott-Armstrong, N.A. et al. GAMETES: a
fast, direct algorithm for generating pure, strict, epistatic models with random
architectures. BioData Mining 5, 16 (2012).
https://doi.org/10.1186/1756-0381-5-16

[2] Ponte-Fernández, C., González-Domínguez, J., Carvajal-Rodríguez, A. et al.
Toxo: a library for calculating penetrance tables of high-order epistasis
models. BMC Bioinformatics 21, 138 (2020).
https://doi.org/10.1186/s12859-020-3456-3

[3] Marchini, J., Donnelly, P. & Cardon, L. Genome-wide strategies for detecting
multiple loci that influence complex diseases. Nat Genet 37, 413–417 (2005).
https://doi.org/10.1038/ng1537


<!-------------------------------- Links -------------------------------->
[GAMETES]: https://sourceforge.net/projects/gametes/
[Toxo]: https://github.com/UDC-GAC/toxo
[additive1]: penetrance_tables/additive_model/3_0.10_0.10.txt
[additive2]: penetrance_tables/additive_model/3_0.10_0.25.txt
[additive3]: penetrance_tables/additive_model/3_0.10_0.50.txt
[additive4]: penetrance_tables/additive_model/3_0.10_0.80.txt
[additive5]: penetrance_tables/additive_model/3_0.25_0.10.txt
[additive6]: penetrance_tables/additive_model/3_0.25_0.25.txt
[additive7]: penetrance_tables/additive_model/3_0.25_0.50.txt
[additive8]: penetrance_tables/additive_model/3_0.25_0.80.txt
[additive9]: penetrance_tables/additive_model/3_0.40_0.10.txt
[additive10]: penetrance_tables/additive_model/3_0.40_0.25.txt
[additive11]: penetrance_tables/additive_model/3_0.40_0.50.txt
[additive12]: penetrance_tables/additive_model/3_0.40_0.80.txt
[additive13]: penetrance_tables/additive_model/4_0.25_0.10.txt
[additive14]: penetrance_tables/additive_model/4_0.25_0.25.txt
[additive15]: penetrance_tables/additive_model/4_0.25_0.50.txt
[additive16]: penetrance_tables/additive_model/4_0.25_0.80.txt
[additive17]: penetrance_tables/additive_model/4_0.40_0.10.txt
[additive18]: penetrance_tables/additive_model/4_0.40_0.25.txt
[additive19]: penetrance_tables/additive_model/4_0.40_0.50.txt
[additive20]: penetrance_tables/additive_model/4_0.40_0.80.txt
[threshold1]: penetrance_tables/threshold_model/3_0.10_0.10.txt
[threshold2]: penetrance_tables/threshold_model/3_0.10_0.25.txt
[threshold3]: penetrance_tables/threshold_model/3_0.10_0.50.txt
[threshold4]: penetrance_tables/threshold_model/3_0.10_0.80.txt
[threshold5]: penetrance_tables/threshold_model/3_0.25_0.10.txt
[threshold6]: penetrance_tables/threshold_model/3_0.25_0.25.txt
[threshold7]: penetrance_tables/threshold_model/3_0.25_0.50.txt
[threshold8]: penetrance_tables/threshold_model/3_0.25_0.80.txt
[threshold9]: penetrance_tables/threshold_model/3_0.40_0.10.txt
[threshold10]: penetrance_tables/threshold_model/3_0.40_0.25.txt
[threshold11]: penetrance_tables/threshold_model/3_0.40_0.50.txt
[threshold12]: penetrance_tables/threshold_model/3_0.40_0.80.txt
[threshold13]: penetrance_tables/threshold_model/4_0.10_0.10.txt
[threshold14]: penetrance_tables/threshold_model/4_0.10_0.25.txt
[threshold15]: penetrance_tables/threshold_model/4_0.10_0.50.txt
[threshold16]: penetrance_tables/threshold_model/4_0.10_0.80.txt
[threshold17]: penetrance_tables/threshold_model/4_0.25_0.10.txt
[threshold18]: penetrance_tables/threshold_model/4_0.25_0.25.txt
[threshold19]: penetrance_tables/threshold_model/4_0.25_0.50.txt
[threshold20]: penetrance_tables/threshold_model/4_0.25_0.80.txt
[threshold21]: penetrance_tables/threshold_model/4_0.40_0.10.txt
[threshold22]: penetrance_tables/threshold_model/4_0.40_0.25.txt
[threshold23]: penetrance_tables/threshold_model/4_0.40_0.50.txt
[threshold24]: penetrance_tables/threshold_model/4_0.40_0.80.txt
[nme1]: penetrance_tables/no_model/3_0.25_0.10.txt
[nme2]: penetrance_tables/no_model/3_0.25_0.25.txt
[nme3]: penetrance_tables/no_model/3_0.25_0.50.txt
[nme4]: penetrance_tables/no_model/3_0.40_0.10.txt
[nme5]: penetrance_tables/no_model/3_0.40_0.25.txt
[nme6]: penetrance_tables/no_model/3_0.40_0.50.txt
[nme7]: penetrance_tables/no_model/3_0.40_0.80.txt
[nme8]: penetrance_tables/no_model/4_0.25_0.10.txt
[nme9]: penetrance_tables/no_model/4_0.25_0.25.txt
[nme10]: penetrance_tables/no_model/4_0.40_0.10.txt
[nme11]: penetrance_tables/no_model/4_0.40_0.25.txt
[dataset1]: datasets/epistasis/additive_model/3/0.10/0.10/
[dataset2]: datasets/epistasis/additive_model/3/0.10/0.25/
[dataset3]: datasets/epistasis/additive_model/3/0.10/0.50/
[dataset4]: datasets/epistasis/additive_model/3/0.10/0.80/
[dataset5]: datasets/epistasis/additive_model/3/0.25/0.10/
[dataset6]: datasets/epistasis/additive_model/3/0.25/0.25/
[dataset7]: datasets/epistasis/additive_model/3/0.25/0.50/
[dataset8]: datasets/epistasis/additive_model/3/0.25/0.80/
[dataset9]: datasets/epistasis/additive_model/3/0.40/0.10/
[dataset10]: datasets/epistasis/additive_model/3/0.40/0.25/
[dataset11]: datasets/epistasis/additive_model/3/0.40/0.50/
[dataset12]: datasets/epistasis/additive_model/3/0.40/0.80/
[dataset13]: datasets/epistasis/additive_model/4/0.25/0.10/
[dataset14]: datasets/epistasis/additive_model/4/0.25/0.25/
[dataset15]: datasets/epistasis/additive_model/4/0.25/0.50/
[dataset16]: datasets/epistasis/additive_model/4/0.25/0.80/
[dataset17]: datasets/epistasis/additive_model/4/0.40/0.10/
[dataset18]: datasets/epistasis/additive_model/4/0.40/0.25/
[dataset19]: datasets/epistasis/additive_model/4/0.40/0.50/
[dataset20]: datasets/epistasis/additive_model/4/0.40/0.80/
[dataset21]: datasets/epistasis/threshold_model/3/0.10/0.10/
[dataset22]: datasets/epistasis/threshold_model/3/0.10/0.25/
[dataset23]: datasets/epistasis/threshold_model/3/0.10/0.50/
[dataset24]: datasets/epistasis/threshold_model/3/0.10/0.80/
[dataset25]: datasets/epistasis/threshold_model/3/0.25/0.10/
[dataset26]: datasets/epistasis/threshold_model/3/0.25/0.25/
[dataset27]: datasets/epistasis/threshold_model/3/0.25/0.50/
[dataset28]: datasets/epistasis/threshold_model/3/0.25/0.80/
[dataset29]: datasets/epistasis/threshold_model/3/0.40/0.10/
[dataset30]: datasets/epistasis/threshold_model/3/0.40/0.25/
[dataset31]: datasets/epistasis/threshold_model/3/0.40/0.50/
[dataset32]: datasets/epistasis/threshold_model/3/0.40/0.80/
[dataset33]: datasets/epistasis/threshold_model/4/0.10/0.10/
[dataset34]: datasets/epistasis/threshold_model/4/0.10/0.25/
[dataset35]: datasets/epistasis/threshold_model/4/0.10/0.50/
[dataset36]: datasets/epistasis/threshold_model/4/0.10/0.80/
[dataset37]: datasets/epistasis/threshold_model/4/0.25/0.10/
[dataset38]: datasets/epistasis/threshold_model/4/0.25/0.25/
[dataset39]: datasets/epistasis/threshold_model/4/0.25/0.50/
[dataset40]: datasets/epistasis/threshold_model/4/0.25/0.80/
[dataset41]: datasets/epistasis/threshold_model/4/0.40/0.10/
[dataset42]: datasets/epistasis/threshold_model/4/0.40/0.25/
[dataset43]: datasets/epistasis/threshold_model/4/0.40/0.50/
[dataset44]: datasets/epistasis/threshold_model/4/0.40/0.80/
[dataset45]: datasets/epistasis/no_model/3/0.25/0.10/
[dataset46]: datasets/epistasis/no_model/3/0.25/0.25/
[dataset47]: datasets/epistasis/no_model/3/0.25/0.50/
[dataset48]: datasets/epistasis/no_model/3/0.40/0.10/
[dataset49]: datasets/epistasis/no_model/3/0.40/0.25/
[dataset50]: datasets/epistasis/no_model/3/0.40/0.50/
[dataset51]: datasets/epistasis/no_model/3/0.40/0.80/
[dataset52]: datasets/epistasis/no_model/4/0.25/0.10/
[dataset53]: datasets/epistasis/no_model/4/0.25/0.25/
[dataset54]: datasets/epistasis/no_model/4/0.40/0.10/
[dataset55]: datasets/epistasis/no_model/4/0.40/0.25/
<!----------------------------------------------------------------------->
