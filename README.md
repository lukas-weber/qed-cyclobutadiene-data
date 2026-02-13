# Undoing the Jahn-Teller effect in molecules using magnetic cavity QED

This is the data repository for our paper [Undoing the Jahn-Teller effect in molecules using magnetic cavity QED](TODO).

The data is in JSON format in the `data` directory.

To regenerate the plots, you need Julia. Running

```bash
julia --project -e "using Pkg; Pkg.instantiate()"
julia --project scripts/make_plots.jl
```

should install all further dependencies and write plots to the `plots` directory.
