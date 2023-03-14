1. Make a code that calculates alpha and beta and plot them

	The code we worked on in class to accomplish this contains an array Vm that holds
	voltages that start at -120 and end at 50, parsed into 100 evenly spaced voltages
	between using linspace. The Hodgkin Huxley alpha and beta values for m,n, and h were
	intialized as an array of 100 zeros to hold the values as the code parses through each
	voltage. The same was done for the noble model alpha and betas. 
	
	A for loop was used to iterate through the 100 steps we created using our voltage array.
	The hodgkin huxley voltage was intialized as the selected first voltage (-120) plus 56
	because in the HH model 0 voltage was assigned to the resting membrane potential which
	was actually -56 mV. The alpha and betas for each HH gateare then calculated using the given
	equations at each voltage in the array, starting with the first voltage on the intial iteration.
	The Noble model initial voltage is then set to -120 mV and the same calculations are done for
	alpha and beta using the noble model equations and voltage. Once this loop has run for every
	gate's alpha and beta, a plot is constructed to display them.

Did the conclusions we made in class still hold?

	Our conclusions from class were that there may be a slight change in the kinetics between
	the two models for sodium, but the main difference used to get a plateau in the action potential is the
	conductance. We also concluded that potassium was changed to be much slower in order to create the plateau. 
	This appears to be true for the m gate, where the plots for alpha and beta are
	similar for both models. The main difference in kinetics is that the sodium m gate closes more
	drastically at low voltages than the noble model. This probably does not make that much of a difference.
	The activation is about the same and the kinetics are probably not driving the plateau. The h (inactivation)
	gate going open to close is about the same. So we can conclude that the kinetics of the h gate are not
	playing much of a role in the plateau either. This leads us to believe that conductance is in fact
	what is driving the plateau of the action potential in the Nobel model. However, for the n (potassium)
	gates, the scale of alpha and beta are a lot smaller which means the kinetics are a lot different from
	the HH model and the kinetics in this case may be driving potassium's part in creating the plateau.

2. Add the L-type Ca2+ channel from Luo and Rudy's 1994 model to the Noble pacemaker model using the given equations
	
	The noble pacemaker model (with Luo and Rudy's added channel) was coded using the solve_ivp differntial
	equation solver function. The variables t (time) and X (V, m, h, n, d, f) were used as inputs. Constants such as
	F,R,T, valence (z), and constants used in the equations (VF/RT) were intialized. Other constants such
	as conductance, nernst potential, permeability (for calcium), concentration of calcium, and capacitance were
	all initialized. As discussed in class there is a factor for gk1 used to make the code work properly.	The given kmca vlaue was also initialized.
	In order to code in current pulses, some values were initialized. They were the start time of the stimulation cycle,
	the period of the stimulation (frequency with which stim was applied), deltat (the stimulation duration of application)
	and Iappamp which is the amplitude of current being applied. If statements were used so that when time is in the stimulation time
	but smaller than the duration of the stimulation the Iapp would be the desired amplitude, otherwise it would be set to 0.
	The sodium channel equations were coded (conductance, current, aplha and beta, the change in m and h over time). The same was done for
	the relevant potassium equations which included two conductances (gk1 and gk2), the potassium current, alpha and beta
	and the change in n over time.
	
	Then, the L-type calcium equations were added. The first was d_infinity, then tau_d, alpha and beta, and the change in d
	over time. Then the f gate was added using the equations for f_infinity, tau_f, alpha and beta, and the change in f.
	over time. The f calcium gate was also added using the given equation. Calcium current was also coded as given. 

	The equations for leak current and membrane potential were also included. The variable ROC was established to be
	the change in voltage, m, h, n, d, and f. The command was then written to return ROC.

	Outside of the solve_ivp function tspan was intialized to give the time span and time step for the differentiation.
	The variable results was made to hold all of the outputs from solve_ivp. Initial values for V, m, h, n, d, and f were selected. I chose
	-80, .1, .8, .01, .9, .5.

	The results where then plotted. One plot has membrane potential. The second plot has m and h gates as well as the n gates plotted. The and f gates
	were plotted in the 3rd plot.



Change something on the pacing or conductance to see a different behavior and explain why you are seeing what you see.

	I altered the pacing of the potassium gate's opening and closing rates (alpha and beta) back to their original hodgkin huxley leading coefficients to see what would happen to the behavior of the gates as well as the overall action potential. One notable change to the action potential is that it still maintains a plateau, however, it is quite small in comparison to the orginally coded Luo-Rudy model coded. There is a notable sharp spike in the beginning of the action potential that was not so prominent in the unchanged model. This is most likely because changing the coefficients back allows the potassium gating to happen much more quickly (HH model is quick compared to the LR model) which allows the potassium current to repolarize the cell right as it spikes instead of a more delayed potassium response (due to the LR coefficients). The reason a slight platea/bump is still seen is because the presence of the calcium is still keeping the voltage somewhat depolarized, its effect is simply dampened by the hyperpolarizing potassium current. This can be seen in the plots of the gates opening/closing which shows the n gate making a slow rise to .75 and a slow descent to close over the time period in the LR model. In the HH model the n gate quickly rises to 1 and nearly as quickly falls back to closed in the first .25 second of the 2 second time period. The d and f gates which are changed due to the differing voltage profile produced by the kinetic changes in the potassium gate, they turn on and off to the same degree (same proportion between 0 and 1), but much more swiftly due to the changed voltage. The same can be said for the sodium gates which open and close much more quickly. This explains why the plateau isn't as sustained because both calcium and sodium contribute to the sustained depolarization that we see as the plateau/bump. Such a small change in dynamics of a single channel can impact the entire voltage profile through slight change in voltages that change the kinetics of all the channel gates.

	