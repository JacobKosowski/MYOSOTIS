import params_clean as params
#Constants

Mbolsun = 4.77

Draine3p1KappaV=8.551E+03
Draine4KappaV=8.492E+03
Draine5p5KappaV=7.313E+03

if (params.Rv == 3.1): 
    DraineKappaV=Draine3p1KappaV
elif (params.Rv == 4.0): 
    DraineKappaV=Draine4KappaV
elif (params.Rv == 5.5): 
    DraineKappaV=Draine5p5KappaV
else: print('For Dmodel, R_V should be 3.1 or 4.0 or 5.5. If you need other Rv values please choose Fmodel')