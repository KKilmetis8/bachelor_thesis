## CompOSE Solver

CompOSE is an online library of various tabulated EoS.\
CompOSE integrates with LORENE. In our case, some -minimal- data-cleanup\
is required since the pipeline is fully automated.

Here I've put the solvers for the TOV eqs for each EoS I studied. \
The structure is interchangeable, except for BSK20. BSK20 is *not* \
included in CompOSE. The authors who derived it, have provided their \
own tables in SIMBAD.

`plot_diffeos.py` contains some plotting tools which rely on the solvers.

The lesson here is that the french are great and their gifts plentiful.

---

Getting into the nitty gritty of it, the BSK20 Table has already derived the *pressure* and the *density* while the CompOSE EoS provide us with:
* Pressure by baryon number density 
$ P/n_b $
* Scaled internal energy per bayron 
$ e/(n_bm_n) - 1 $
* Baryon density
$ n_b $

Combining these quantities we can derive *pressure* and *density*.