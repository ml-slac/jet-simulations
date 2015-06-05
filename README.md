# `jet-images` --  a framework for jet image event generation and processing.

## Dependencies

* `numpy`, `matplotlib`, and `scikit-image` for python.
* 

## Event generation.

Event generation is handled in the `./event-gen/` folder. If you go in that folder and type `make`, everything should build.
This generates the low level event generation script, which can be invoked as `./event-gen`. However, to speed things up, there is a multicore wrapper around this generation process in `generateEvents.py`. Calling `python generateEvents.py --help` yields:

```
usage: generateEvents.py [-h] [--outfile OUTFILE] [--nevents NEVENTS]
                         [--ncpu NCPU] [--process PROCESS] [--pixels PIXELS]
                         [--range RANGE] [--pileup PILEUP]
                         [--pt_hat_min PT_HAT_MIN] [--pt_hat_max PT_HAT_MAX]
                         [--bosonmass BOSONMASS]

optional arguments:
  -h, --help            show this help message and exit
  --outfile OUTFILE
  --nevents NEVENTS
  --ncpu NCPU
  --process PROCESS     Can be one of ZprimeTottbar, WprimeToWZ_lept,
                        WprimeToWZ_had, or QCD
  --pixels PIXELS
  --range RANGE
  --pileup PILEUP
  --pt_hat_min PT_HAT_MIN
  --pt_hat_max PT_HAT_MAX
  --bosonmass BOSONMASS
```







```
python jetconverter.py --signal=wprime --save=newfile.npy --plot=jet_image *.root
```


