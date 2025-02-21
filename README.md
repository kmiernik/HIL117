# HIL117

[![Build Status](https://github.com/k.a.miernik@gmail.com/HIL117.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/k.a.miernik@gmail.com/HIL117.jl/actions/workflows/CI.yml?query=branch%3Amain)

This project is intended for analysis of data collected in the Heavy Ion Laboratory
at University of Warsaw experiments using EAGLE, NEDA and DIAMANT detectors.

The projected in tailored for a specific experiment (HIL117 as you might 
guess). Adjusting it to another one might require some tinkering. Please 
feel free to contact the author if you wish to do so. 

Most of the settings are done through a human-readable config file (TOML).
Base config is needed to start scanning the data, some information is acquired 
automatically during a pre-scan phase (time shifts, energy calibration).

## Nomenclature
Whenever word ''hit'' is used it means individual channel (detector) information
containing energy, time and possibly other data (e.g. pid).
''Hits'' are grouped into ''events'' - physically related groups of hits.

## Workflow

## Channels setup
Each channel has its number (label) in config file. The label is calculated
as (board + 16*channel + 1). Notice that DIAMANT crystal/board are recalculated
to match that scheme. Mapping the diamant board number is done in `diadat.jl`, 
see `board_translation` dict.  

There are number of configuration fields for each label e.g.
```
[label.9]
	type = "Ge" - type of detector
	ring = "143" - position in spectrometer
	phi = "0" - also, but currently not used
	dt = 213.851 - time shift, calculated automatically
	cal = [0.5178388739719042, 0.16172422866441732] - calibration see below
	valid = true - if set to false, the channel will be ignored
	comment = " " - just human readable comment
```

Calibration may have different meaning for different detectors. For HPGe it
is a energy calibration, which is calculated automatically (the values are
ignored unless pre-scan is omitted). For NEDA first number gives low limit 
for gamma events, second for neutrons. For DIAMANT the values are used for
calculating PID. The relevant functions for DIAMANT/NEDA are in `dev.jl`
(neda_lims, dia_banana_fit), and in `scripts/` for a graphic version. The
''calibration'' for NEDA and DIAMANT needs to be set, as it is not recalculated
during the pre-scan.

### Calibration
Ge calibration uses two selected lines, indicated in this section
```
[calibration]
	E1 = 511.0
	E2 = 1256.69
```
First line should be the maximum in spectrum (511.0 keV is quite obvious), 
and easy to locate. Second one should be significant enough to be found for 
each detector, its location is calculated based on the first one.

### Particles identification
PID is governed by the following section. This code uses only the on-board
PID for NEDA and DIAMANT, limits of given type of particles are self-describing
```
[pid]
	g_low = 0.71
	g_high = 0.95
	n_low = 0.5
	n_high = 0.69
	a_low = 0.42
	a_high = 0.56
	p_low = 0.56
	p_high = 0.67
```
As mentioned before NEDA has also low limit in channel number for gamma/neutron
set in calibration field, while diamant uses calibration to calculate PID so
alpha particles are around value of 0.5 and independent of energy.

## Event building
Events window may very depend on the detector type. For each one it should
be set in the following section.
```
[event]
	ge_low = 30.0
	bgo_low = 2.0
	neda_low = 70.0
	dia_low = 100.0
	beam_period = 74.0
	ge_dt = 148.0
	bgo_dt = 54.0
	neda_g_dt = 17.0
	neda_n_dt = 37.0
	dia_dt = 54.0
	t_delay = 20.0
```
The *_low fields are for removal of low-channel noise before any further
 processing. All channels, except of NEDA, are delayed by the `t_delay` value, 
 so they can match beam spill detected by the NEDA, and the natural spread (width)
 if events is somehow included - notice that the automatically calculated shift 
 is set to the maximum of time difference. The timing can be verified in the 
 `tloc` spectra (time-location).

## Scanning
Once the setup is done, the main function for scanning - `scan_run`. 
- can be called. It will first run a pre-scan will to find shifts and 
calibration.  Default output is a 
[HDF5](https://www.hdfgroup.org/solutions/hdf5/) file containing various spectra. 
Their parameters can be adjusted by the following section
```
[spectra]
	dE = 1.0            (energy bin in keV or ch)
	Emax = 4096.0       (max. energy in signles)
	E2max = 2560.0      (max. energy in gamma-gamma)
	dt = 1.0            (time bin in ns)
	tmax = 255.0        (max. time difference)
	dpid = 0.01         (PID step)
	pidmax = 1.27       (max. PID value)
	Mmax = 10           (max. multiplicity in gamma-multiplicity spectrum)
```

The spectra are defined in the `prepare_spectra_file()` function in `scan.jl`. 
If you need to add or modify something follow the examples of existing ones.
The `spectra` variable is of type NamedTuple.  Spectra are updated in 
the `update_spectra!()` function (`eventbuilder.jl`), the syntax is basically:
```
spectra.name[i, j] += 1
```
or something similar depending on the case.

## EDF mode
EDF - Event Data Format - is a simple and basic binary format for storing 
data on event by event basis. Each event starts with a header.
```
mutable struct EDFHeader
    Mraw::UInt8
    M::UInt8
    ge::UInt8
    gammas::UInt8
    neutrons::UInt8
    protons::UInt8
    alphas::UInt8
end
```

and then is followed by `M` hits (`Mraw` indicates number of hits before 
removal of BGO suppresed HPGe, neutron scattering etc., number of written hits 
`M` is always lower or equal to `Mraw`).  Field give number
of HPGe, neda-gamma, neda-neutrons, diamant-protons, and diamant-alphas hits.
```
struct EDFHit
    loc::UInt8
    type::UInt8
    E::UInt16
    t::UInt16
end
```
field `loc` gives hit label as defined in the config file. `type` indicates
particle type (gamma-1, neutron-2, proton-3, alpha-4, see `structs.jl` 
enum type `ParticleType`).

If you are on the dark side of the power called CERN ROOT, you might consider
creating your precious trees, branches, and leafs out of EDF files, or even
instead of, which might require some
 [C++ wrapper](https://github.com/JuliaInterop/CxxWrap.jl?tab=readme-ov-file),
 or using scheme suggested by the Julian library for reading ROOT data -
 [UnROOT](https://juliahep.github.io/UnROOT.jl/stable/exampleusage/).