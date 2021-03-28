# savitzky_golay_filter
The function applies Savitzky-Glolay filter to noisy timeseries data at non-uniform spacing using the weights as well.

Applies a Savitzky-Golay filter to y with non-uniform spacing
as defined in x

This is based on 
#https://dsp.stackexchange.com/questions/1676/savitzky-golay-smoothing-filter-for-not-equally-spaced-data
#Similar function is available for R as well at
#https://rdrr.io/github/ranghetti/sen2rts/src/R/w_savgol.R
#rangetti's R script also uses weights in addition to non-uniform spacing
#So I followed his R script to include the weights parameter in my function
The borders are interpolated like scipy.signal.savgol_filter would do

Parameters
x : array_like
	List of floats representing the x values of the data
y : array_like
	List of floats representing the y values. Must have same length
	as x
weight: array_like 
	List of floats between 0 and 1 representing y values. Must have 
	same length as x and y
window : int (odd)
	Window length of datapoints. Must be odd and smaller than x
polynom : int
	The order of polynom used. Must be smaller than the window size

Returns
np.array of float
	The smoothed y values
