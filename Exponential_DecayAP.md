The code generated for this exercise uses the python ODE solver with a custom function to model
the equation that governs the voltage of a cell that has free flowing sodium and potassium.
Considerations are psi (membrane voltage), the nernst voltage of each ion, the capacitance of
the cell (assumed to be 1e-6), and the conductance of sodium and potassium. These conductances
govern the flow or "current" of the two ions. Inward current will cause the voltage to depolarize
(increase) while outward current will cause the voltage to hyperpolarize (decrease). The Nernst
potential of sodium (which is part of its driving force (which is the gap between the current voltage
and the Nernst voltage)) is 55 mV. The Nernst potential of Potassium is -80 mV. We begin 
the simulation at an intial membrane potential psi of -70 mV. The conductances of sodium and potassium
also begin at 0. To simulate the gating seen in action potentials, if statements based on time
periods in the simulation control the values of each ion's conductance. The first time period
has the conductance of potassium remain at 0 while Na raises to 4e-3. This causes an inward flow 
of current that depolarizes the voltage to over 40 mV. Then at the next time period, the conductance
of potassium raises to 1e-3 while the conductance of sodium reduces to .1e-3 which causes the voltage
to hyperpolarize back towards rest which was around the nernst of potassium by design. The equation
dpsi/dt =(-gna*(psi-NernstNa)-(gk*(psi-Nernstk)))/C is entered and the function is asked to
return dpsi/dt at each timestep. I create 1000 timesteps between 0 and .01 to solve. I then use
odeint built in function calling my model, intial membrane voltage condition, and my time vector t
to run the code. At the end I plot the voltage versus time to visualize this solution. 

I played around with the code to introduce a repetitive spiking action. After the first spike 
was finished a made sure the voltage was back at rest using another time segmented if statement
that set dpsi/dt to a set value of -80 mV. I then use the rest of my time (now 1000 segments
between 0 and .1) to redo the segments that made the first spike to generate a second consecutive spike.
Another interesting way to introduce repetitive spiking would be to mimic something known as the LIF
(leaky integrate and fire) model. In this model once a certain voltage threshold is reached (an
if statement based on the current value of psi rather than time), the membrane voltage (dpsidt) would
be set back to rest instantaneously, and then the code would go back to functioning as before to hit
the peak of the next spike. If accomplished with a for or a while loop, the cycle could continue
over and over again over any selected period of time. This may be easier with voltage gated channels
rather than time dependent conductance changes like we have here.
