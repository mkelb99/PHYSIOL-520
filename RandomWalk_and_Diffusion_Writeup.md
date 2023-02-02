I began this assignment using the random walk and diffusion code that we built and
discussed during the class period. We discovered through exploration of set variables
that depending on things like the selected time and space steps, the profiles may look different,
but there is a set of variables that make the two profiles line up near identically, as they 
both represent the same phenomena in different ways. When N was halved from the value that
matched the two profiles(dx increases 2 fold), N_frames (timesteps) for the random walk had 
to decrease 4 fold in order for the profiles to match up again. This helps us determine the
relationship between the time and space steps and the diffusion coefficient D in the diffusion 
equation. The units of D are length^2 (cm^2) divided by unit time. This is why in order for the
random walk to match the diffusion equation again, any change in dx must be balanced by a 1/(change in dx)^2 
change in timestep (to get the ratio of dx^2/time to equal 1, this must occur).