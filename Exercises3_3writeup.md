The framework of this exercise was the example we went over in class on Thursday which had the libraries imported as well 
as the my_Gaussian function as defined in the textbook. Values of D and M were selected in this code and the value for t
was commented out as this problem wants multiple values for t. The first thing I did to complete part B was select random, somewhat far
apart values of t to explore the width of the gaussian as t changes. These values were stored as a numpy array called times in order to easily access
them later. X was varied between -.35 to .35 in 500 evenly spaced increments as it had been in the prior class example. This is accomplished 
using the linspace command from the numpy library. In order to easily test multiple values of t in a function that also uses multiple values
of x, a for loop was created to iterate through and create the different plots together. The syntax for t in times allows t to take on each
value in the array times (one for each iteration of the loop). Then, I call the command my_Gaussian and name the output of that function Y
as it is the Y-axis value for the Gaussian plot. On the first iteration, it passes the argument X into the function, so the Gaussian result
Y is an array with each value of Y outputted for the corresponding X value. These Y Values are then plotted against the X (distance) values
in the first subplot which is called using the plt command subplot (1,2,1 means it is 1 row and 2 columns of plots and this plot is mapping 
to the plot in the first column). I then label the axes. This occurs for each value of t one at a time in each loop. For instance t has
a value of 600 on the first iteration, which is used as the variable t in the my_Gaussian function. This is why the result is a plot with
3 gaussians of varying heights and widths, because each plot corresponds to a different t. Since I had 3 values of t, there are 3 gaussians
plotted. All lines that are part of the for loop come after line 49 and are indented slightly to indicate they are part of the loop. The
next unindented section are independent lines of code.

I used these plots to find the inflection points for the exercise. Inflection point was described as the change from the concave down 
to concave up parts of the graph. For me I called this the midpoint between the vertexes of the concave up and concave down sections. I 
guessed the x-value of these points visually as instructed in the exercise. These points were to find the width and since the gaussian is 
symmetrical at x=0, each plot has 2 inflection points that are equally spaced from 0. To find the width I simply multiplied the x-value
of one inflection point by 2, as seen in line 59. I called the values for all 3 plots "inflection" and then once again saved them as an
array, this time so I could easily plot these values against my array of times. Subplot was called to put the next plot in the second
column plotting space, and I plotted the inflection points against my 3 values of time. The axes were labeled and plt.show was used to 
display the plots in the spyder environment. What can be seen in this second plot is that larger values of time result in wider Gaussian 
plots/profiles (inflection's value increases with increasing values of time).