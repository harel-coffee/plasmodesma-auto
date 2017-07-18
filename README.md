# plasmodesma
A project to automatically process a set of 1D and 2D NMR experiments

*This program is associated to the Manuscript*

# Automatic differential analysis of NMR experiments in complex samples

*Laure Margueritte,
Petar Markov,
Lionel Chiron,
Jean-Philippe Starck,
Catherine Vonthron-Sénécheau,
Mélanie Bourjot,
and Marc-André Delsuc*

*to be submitted*

----
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

## Files
The project contains several files

The project contains several files

**Plasmodesma.py** The program itself.

run it with
```
python Plasmodesma.py -d folder_containing_all_data_Sets
```

**BucketUtilities.py**  a set of utilities used for the analysis step

**Analysis.ipynb** an example of analysis of a set of experiments. Detection of varying amount of artemisin in a plant extract is used as an example

**spike_plugins** plugins for SPIKE library



