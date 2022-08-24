# Model_MH

If you find some bug on the code, have questions or suggestions feel free to send an email!

Model_MH_v4 changes - 24/08/2022
- All calculations implemented in NumPy, improving the speed of the code
- The code will generate all plots and data on c_H and %wt H units
- Now the configurational entropy is ploted
- All code was rewrite to be easily understood, using the same logic as explained on the paper: https://doi.org/10.1016/j.ijhydene.2022.07.179
- Small corrections on user input parameters treatment were made
- In order to take less time, the determination of chemical equilibrium between the phases changed to be done using NumPy now
- The step on Hydrogen composition of each phase for the calculations can be set for each phase separately by changing the values of the variables: alpha.cH_step, beta.cH_step and delta.cH_step, by default the three are the same and equal to 0.00025 or the value given by the user when code is runned
