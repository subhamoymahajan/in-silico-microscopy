# Tutorial 11: Generating in silico images with fluorophores undergoing photophysical processes

## Format for photophysical constants filed (.pp)

line 1: {Number of fluorophore type} {Timestep of photophysical process}

Timestep of photophysical process can be different from timestep of coordinates. The time axis of photophysical processes is scaled so that timestep of photophysical process is equal to timestep of the coordinates. Replace with the value, do not use curly brackets.

Line one is followed by `{number of fluorophores type}` blocks of information associated with each fluorophore type. The lines associated with each box is shown below.

line x: [ Primary emission wavelength ]

Keep the square bracket. Keep the space after [ and before ]. 

line x+1: {Number of photophysical states} {Number of photophysical states that produce emissions}

If states are S0, S1, S2, T1, T2, {number of photophysical states} is 5. If transitions S1->S0, S2->S0, T1->S0 produce emissions, {number of photophysical states that produce emissions} is 1 and is denoted by S0f.  
If transitions S1->S0 and S2->S1 produce emissions then {number of photophysical states that produce emissions} is 2 and is denoted by S0f and S1f.

line x+2: {list all states} {list of emission states}.

If states are S0, S1, S2, T1, T2, and transitions S1->S0 and S2->S1 produce emissions then the line is `S0 S1 S2 T1 T2 S0f S1f`.

line x+3 to x+n1: {state i} {state j} {rate constant}

The unit of {rate constant} is the inverse unit of {timestep of photophysical process}.
If {state i} to {state j} produces emission replace {state j} with {state j}f. For example, If S1 -> S0 produces emission and transition occurs at 1E+6 second inverse, then the line would be `S1 S0f 1E+6`,


line x+n1+1 to x+n2: {statei} {statej} -1 {emission wavelenght in nm}

{state i} is an emissive state, and {state j} is the corresponding non-emissive state. Example, `S0f S0 -1 670`.

 
Example of a 3 state triplet photophysics.
```bash
1 2E-5
[ 670 ]
3   1
S0  S1  T1 S0f
S0  S1  2E+7
S1  T1  3.333333333333E+6
T1  S0  2E+5
S1  S0f 2E+7
S0f S0  -1   670
```


## Differences in parameters.dat for using photophysical processes
1. Add `pp_file =  [.pp file]`
2. Change `lam[i] = {primary wavelenght}` to `lam[i] = {primary wavelenght} ; {emissive wavelengths including primary}`
3. Add `pos_prec = [precision of .gro positions]`. 
4. The line `max_Nf = {maximum number of emissions associated with primary wavelength 1} ; {maximum number of emissions associated with primary wavelength 2}`. This line should be automatically generated afeter step (2) shown below.
5. Make sure the color mixing scheme is `mix_type = nomix`.

 
## 1. Extract GRO files
Download `Struct.tar.gz` from University of Alberta Dataverse DOI: https://doi.org/10.7939/DVN/F3JKZH, version 3.0
```bash
mkdir -p GROs
tar -xzf Struct.tar.gz -C GROs
```

## 2. Generate Specimen (.spm) files 
```bash
mkdir -p SPMs
siliscopy gen_spm --data spm_data.dat --paramfile parameters.dat --output SPMs/out  
```

## 3. Generate Image data
```bash
siliscopy gen_mono_pp --data imggen.dat --multiprocess
```

### 4. Generating Image/Video 

Follow any of the previous tutorials (Tutorial 1-8). 

For example,
```bash
siliscopy plot --file pp_img --paramfile parameters.da --method color2dt --calc specific --output img --type tiff8
```


