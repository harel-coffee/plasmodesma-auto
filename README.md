# plasmodesma
A project to automatically process a set of 1D and 2D NMR experiments

*This program is associated to the publication*

# Automatic differential analysis of NMR experiments in complex samples

*Laure Margueritte,
Petar Markov,
Lionel Chiron,
Jean-Philippe Starck,
Catherine Vonthron-Sénécheau,
Mélanie Bourjot,
and Marc-André Delsuc*

*Magn. Reson. Chem.*, (2018) **80** (5), 1387. http://doi.org/10.1002/mrc.4683

---

## Abstract

Liquid state NMR is commonly used for the analysis of mixtures containing unknown molecules,
and this capacity has been used in many analytical approaches:
metabolomics, identification of active compounds in natural extracts, characterization of species.
Such studies require the acquisition of many diverse NMR measurements on series of samples.

The Plasmodesma program allows the autonomous, unsupervised processing of a large corpus of 1D, 2D and DOSY experiments from a series of samples acquired in different conditions.
The program relies on the SPIKE library ( https://bitbucket.org/delsuc/spike ) for processing and analysis.
It concentrate on the detection and identification step of unknown compounds, as is usually required in the identification of an active molecule in a natural extracts starts with its detection and then its characterization of an unknown compond or eventually a family of related species.

It provides all the signal processing steps, as well as peak-picking and bucketing of 1D and 2D spectra.
The main step allowing a strong reduction of the size of the data is the bucketing of the area under the spectra.
Plasmodesma is a novel tool which should be useful to decipher complex mixtures, particularly in the discovery of biologically active natural products from plants extracts, but can also in drug discovery or metabolomics studies.

---

*it was further developped for the publication*
## Automatised pharmacophoric deconvolution of plant extracts application to cinchona bark crude extract
*Margueritte L., Duciel L., Bourjot M., Vonthron-Sénécheau C., Delsuc M-A.*

Faraday Discussions (2019) **218**, 441-458 
http://doi.org/10.1039/c8fd00242h

---

## Files
The project contains several files

**Plasmodesma.py** The program itself.

run it with
```
python Plasmodesma.py -d folder_containing_all_data_Sets
```
or
```
python Plasmodesma.py -d folder_containing_all_data_Sets -n 4
```
if you want to parrallelized on 4 processors


**BucketUtilities.py**  a set of utilities used for the analysis step

**spike_plugins** plugins for SPIKE library

## Analysis examples
- **Analysis.ipynb** the *ARTE* sample: Detection of varying amount of artemisin spiked in a plant extract.
	- 5 fractions were separated from the extract and spiked with various amount of artemisinine
	- artemisinine signals are detected in 2D COSY pairs of spectra by bucket comparison
	- linear regression on spiked concentration shows the quality of the analysis obtained on COSY, HSQC and HMBC experiments.

- **Analysis2.ipynb** the *SMARTE* sample: The same plant extract with artemisinine added before fractionation, and the biological activity of each sample was measured.
	- linear regression on activity obtained on DOSY, HSQC and combining DOSY and HSQC.


## DATA
The data presented in the publication, and used in the Analysis.ipynb file are not in this repository for the moment.
They can be downloaded from our laboratory web site, at the following address :
https://pydio.igbmc.fr/pydio/public/8b07c3c6a58c381


---
*This Program is provided under the [http://www.cecill.info/licences.en.html](Licence CeCILL 2.1)* 

*This program HAS NOT been tested intensively, it is believed to do what it is supposed to do, However, you are welcome to check it on your own data.*

    Authors : Petar Markov, Laure Margueritte, Marc-André Delsuc (madelsuc@unistra.fr)
    Version : 1.0   Date : June 2017 (corresponding to Plasmodesma_v6_2)
    Version : 1.1   Date : September 2017  (corresponding to Plasmodesma_v6_3)
    Version : 1.2   Date : October 2017  (corresponding to manuscript revision)
    Version : 1.3   Date : February 2018  (corresponding to Plasmodesma_v6_4)
