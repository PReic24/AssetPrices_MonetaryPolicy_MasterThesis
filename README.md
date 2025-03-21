

# AssetPrices_MonetaryPolicy_MasterThesis
This repository contains the code needed to replicate the findings, data and figures of the Master Thesis "On the Relation between Asset Prices and the Stance of Monetary Policy". 

The file structure and code components are derived from the replication package of 
"*Pflueger and Rinaldi (2022), Why does the Fed move markets so much? A model of monetary policy and time-varying risk aversion, In: The Journal of Financial Economics*" published on Caroline Pflueger's GitHub page under:
> https://github.com/cpflueger/ProgrammingPackage_public/tree/main/ModelFOMCResponses
>

Files in this repository contain the names "Pflueger" or "Reichert" in the header-description of each file, stating the authorship of the files creation.

## Structure of the Code

**Clean_Simulation.m** 
This file is the main script of the code, based largely on *minimum_working_example_NK.m* of Pflueger & Rinaldi (2022). Further information on the functionality of the replication package by Pflueger & Rinaldi (2022) can be found in the file *Minimum_Working_Example_Manual_NK.pdf*. By design, *Clean_Simulation.m* loops over the set of *\gamma_x* and *\gamma_pi* parameters, solving the model and generating asset moments. The different simulation strategies mentioned in Section 4.2 and 4.3 are run by the use of the dummy variable *simulation_run*.
- asset moments for Base Calibration and Strategy A are printed in the console vs. the empirical data, when setting *simulation_run* == 0 (base case) or == 1 (Strategy A). To replicate asset moments of Strategies B.1 to C.3 in table form, similar to Table 3, 4 and 5 of the thesis, enter the corresponding parameter pair under *elseif simulation_run == 1* in line 89. 
- By setting *simulate_irfs* == 1, impulse response functions with 5000 iterations are computed. **Caution**: This slows down the code considerably, and did not change asset moments. If set to 0, only structural impulse response functions are computed, speeding up the simulation times. 

**simulateMoments_wthout_FOMC.m** 

This function is based on the function *simulateMoments* in the script *asset_p.m* and was adjusted to only include components that are necessary for this thesis. Therefore, components that were included by Pflueger & Rinaldi (2022)  to generate results mentioned in Section 4 of Pflueger & Rinaldi (2022), are not included in this function. The exclusion of said components was tested not to distort the moments generated.

**Plotting Impulse Response functions** 

The following functions are employed to plot impulse response functions for a given calibration strategy discussed in Section 4.2:

 - **plot_strategy_a_irfs.m**
 - **plot_strategy_b_irfs.m**
 - **plot_strategy_c_irfs.m**

These functions plot the impulse response functions for variables:

 - Output Gap
 - Nominal Interest Rate
 - Equity Returns
 - 10-Year Nominal Bond Yields

**Caution**: The code of the functions depends on the current format of Calibration Strategies in the script *Clean_Simulation.m*. Changes applied to this format, might result in the failure to generate figures.

**Plotting Asset Moment Surfaces**
For plotting the asset moment surfaces,  the function **plot_asset_moment_surfaces.m**, which uses the MatLab Function *surf(x, y, z)* has been employed. To simulate the asset moments on the the grid following Strategy D, set *simulation_run* == 5. After simulating, the asset moment surfaces are plotted automatically.

## Required, unchanged scripts and functions from Pflueger & Rinaldi (2020)

The following scripts and functions are indispensable for the correct functioning of the code. Other than that, these scripts and functions are almost always identical to the originals provided in the replication package of Pflueger & Rinaldi (2020).

 - asset_p.m:
 - eqinterpnND2.m
 - GaussLegendre.m
 - macro_dyn.m
 - num_set.m
 - senshat.m
 - Uncon_var.m 

## Legacy scripts and functions from Pflueger & Rinaldi (2020)

 - minimum_working_example_NK.m: Original main-file of the replication package. Generates Tables 3 and 4 and Figure 3 of Pflueger & Rinaldi (2022)
 - minimum_working_example_NK.pdf: Details the functionality of the replication package
 - plot_StructuralIRF_STMP.m: Legacy plotting function to replicate Figure 3 in Pflueger & Rinaldi (2022). Can be used to test, that *Clean_Simulation.m* is able to replicate results of Pflueger & Rinaldi (2022), when setting *simulation_run*  == 0.

