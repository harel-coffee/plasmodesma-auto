# Documentation

*This is a very preliminary documentation, aimed at installing and parametring the process...*

## Principle of operation
The principle of the prgm is that the all the experiments are processed automatically from default parameters.
Most processing methods are hard-wired in the program, using sensible algorithms and parameters.
A set of parameters can be modified by the user.

## Installation

Use `git_hub.com/delsuc/Plasmodesma`

it requires Spike and its dependences (numpy scipy) and uses python 3.9 or higher.

## Running the program
### Command syntax

Typically:

        python Plasmodesma -D MyProject

where `MyProject` contains the NMR experiments

or        

        python Plasmodesma --help


options:
```
  -h, --help            show this help message and exit
  -D DIREC, --Data DIREC
                        DIRECTORY_with_NMR_experiments, default is "."
  -N NPROC              number of processors to use, default=1
  -n, --dry             list parameters and do not run
  -T, --template        Generate default config files templates (parameters.json RunConfig.json)
```

### Data organisation
The NMR data-sets have to be organised in a specific manner.

The principle is the following:

- say you have 4 samples with varying conditions:

        sample1 sample2 sample3

each in a Topspin experiment

- you run on each a set of 1D and 2D experiments, each in a given `EXPNO`, eg:
    - 1: 1D
    - 10: COSY
    - 20: HSQC
    - 30: DOSY
    - 100: final 1D to compare with 1 to check if there is sample degradation during the experiment

You create a Project folder to hold everything, say `MyProject` with the following structure
*(the  files ` RunConfig.json, parameters.json` are described in the next section)*

```
MyProject/
    RunConfig.json
    parameters.json
    sample1/
        1
        10
        20
        30
        100
    sample2/
        1
        10
        20
        30
        100
    sample3/
        1
        10
        20
        30
        100
```
After running the program with:

        python Plasmodesma.py -D MyProject

you will get an additional set of files in `MyProject`
```
MyProject/
    report.csv
    analysis.csv
    Config.dump
    Results/
        sample1/
            1D/
                1.pdf
                1.png
                1_pp.pdf
                1_pp.png
                1_peaklist.csv
                1_bucketlist.csv
                100.pdf
                100.png
                ... etc
            2D/
                cosy_10.pdf
                cosy_10.png
                cosy_10_pp.pdf
                cosy_10_pp.png
                cosy_10_peaklist.csv
                cosy_10_bucketlist.csv
                hsqc_20.pdf
                hsqc_20.png
                ... etc
        sample2/
            ... etc
```

- `Results/` contains the the result of the processing, in 1D and 2D, graphics spectra (with ( `*_pp.pdf` )  and without peak picking) as well as peaklist and bucketlist in `.csv` format.
- `Config.dump` is a json dump of the configuration used for processing
- `report.csv` contains a summary of the experiments, 
- `analysis.csv` details the result of the processing for each experiment (number of detected peak, buckelist statistics, etc...)

## Parametrisation of the processing
The parameters used for the processing can be modified by the user, there are set-up in two different files.
The files are located alongside the sample folders containing the dataset, and are thus specific to this set of experiments.

All entries in both files are fully optionnal and the program uses default values when missing.

### general parameters -- `RunConfig.json`
The `RunConfig.json` file to modified the values.
Documented list of possible entries - given here with the default values.

```
{
    'NPROC' : 1,            # The default number of processors for calculation, if  value >1 will activate multiprocessing mode
                            # for best results keep it below your actual number of cores ! (MKL and hyperthreading !).
    'BC_ALGO' : 'Spline',   # baseline correction algo, either 'None', 'Coord', 'Spline' or 'Iterative'   
    'BC_ITER' : 5,          # Used by 'Iterative' baseline Correction; It is advisable to use a larger number for iterating, e.g. 5
    'BC_CHUNKSZ' : 1000,    # chunk size used by 'Iterative' baseline Correction,
    'BC_NPOINTS' : 8,       # number of pivot points used by automatic 'Spline' baseline Correction
    'BC_COORDS' : [],       # coordinates of pivot points in ppm used by 'Coord' baseline Correction
    'TMS' : True,           # if true, TMS (or any 0 ppm reference) is supposed to be present and used for ppm calibration
    'LB_1H' : 1.0,          # exponential linebroadening in Hz used for 1D 1H 
    'LB_13C' : 3.0,         # exponential linebroadening in Hz used for 1D 13C
    'LB_19F' : 4.0,         # exponential linebroadening in Hz used for 1D 19F
    'MODUL_19F' : False,    # 19F are processed in modulus
    'SANERANK' : 20,        # used for denoising of 2D experiments, sane is an improved version of urQRd
                            # typically 10-50 form homo2D; 5-15 for HSQC, setting to 0 deactivates denoising
                            # takes time !  and time is proportional to SANERANK (hint more is not better !)
    'DOSY_LAZY' : False,    # if True, will not reprocess DOSY experiment if an already processed file is on the disk
    'PALMA_ITER' : 20000,   # used for processing of DOSY
    'BCK_1H_1D' : 0.01,     # bucket size for 1D 1H
    'BCK_1H_2D' : 0.03,     # bucket size for 2D 1H
    'BCK_1H_LIMITS' : [0.5, 9.5],   # limits of zone to  bucket and display in 1H
    'BCK_13C_LIMITS' : [-10, 150],  # limits of zone to  bucket and display in 13C
    'BCK_13C_1D' : 0.03,    # bucket size for 1D 13C
    'BCK_13C_2D' : 1.0,     # bucket size for 2D 13C
    'BCK_19F_LIMITS' : [-220, -40],  # limits of zone to  bucket and display in 19F
    'BCK_19F_1D' : 0.1,    # bucket size for 1D 19F
    'BCK_19F_2D' : 1.0,     # bucket size for 2D 19F
    'BCK_DOSY' : 0.1,       # bucket size for vertical axis of DOSY experiments
    'BCK_PP' : False,        # if True computes number of peaks per bucket (different from global peak-picking)
    'BCK_SK' : False,       # if True computes skewness and kurtosis over each bucket
    'TITLE': False,         # if true, the title file will be parsed for standard values (see documentation in Bruker_Report.py)
    'PNG': True,            # Figures of computed spectra are stored as PNG files
    'PDF': False,            # Figures of computed spectra are stored as PDF files
    'addpar': [],           # additional parameters for report.csv : eg ['D2', 'D12', 'P31']
    'add2Dpar': [],
    'addDOSYpar': []
}
```

In `Plasmopera`, you've got in addition:
```
    'GB_19F' : 32.0,        # gaussian linebroadening in Hz used for 1D 19F
    ...
    'safe_NS_scaling': 256, # All OPERA_safe experiments are scaled to this NS
```
Be carefull that the `json` format does not accept comments (which should all be removed if you use the example above), and that the syntax is strict (balacing `{}` `[]` and commas).

A template is created for you with the command

        python Plasmodesma.py -T -D Target

### experiment specific parameters -- `parameters.json`
Additional processing parameters can be assigned per experiments.
These parameters are given in a file called `parameters.json` located at the place than `RunConfig.json`
It has the following organisation (given here as an example):

```
{
"sample1/1": {
        "remark": "1D experiment",
        "ph0": 33.2,
        "ph1": 70
    },
"sample1/20": {
        "remark": "HSQC experiment",
        "SANERANK" : 20
    },
"sample2/10": {
        "remark": "COSY experiment"
    },
"sample3/100": {
        "remark": "1D experiment",
        'LB_1H' : 0.3,
        'ppm_offset' : 0.12,
        "ph0": -13
        "ph1": 0
    }
}
```
The only parameters recognized so far are:

- `ph0`, `ph1` : an additional phase correction applied to the one stored into dataset
- `ppm_offet` : and offset in ppm to the experiment in the acquisition dimension
- all global parameters found in `RunConfig.json` which are overiden for this experiment only

In this file, all entries and parameters are optionnal, and all entries are Ok, the "remark" entry is available to write any usefull comment.

A template for both files is created for you with the command

        python Plasmodesma.py -T -D Target


