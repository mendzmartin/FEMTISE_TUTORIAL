Inside `~/quantum_harmonic_oscillator_2D/010_script/` folder you have two script

+ `input_generator.sh`: to this script we generate input file inside `~/quantum_harmonic_oscillator_2D/020_input/` folder (see `input.dat` file).
+ `run_generator.sh`: to this script we generate run Julia code inside `~/quantum_harmonic_oscillator_2D/030_src/` folder (see `run.jl` file).


Is mandatory run `input_generator.sh` script as following:
```bash
    @prompt$: ./input_generator.sh
```
if we want to run analysis's notebook inside `~/quantum_harmonic_oscillator_2D/050_analysis/` folder because `input.dat` file has specific path files wich depends on what computer we are running the script.