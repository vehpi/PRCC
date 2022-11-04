# PRCC
Function for calculating Partial (Rank) Correlation Coefficients


@author: Vehpi Yildirim, v.yildirim@amsterdamumc.nl, https://github.com/vehpi

partial_corr(A): 
Returns the sample linear partial or rank correlation coefficients and associated p-values between pairs of variables in A, controlling 
for the remaining variables in A.

partial_corr(A,B): 
Returns the sample linear partial or rank correlation coefficients and associated p-values between pairs of variables in A and B, controlling 
for the remaining variables in A.

partial_corr(A,B,C): 
Returns the sample linear partial or rank correlation coefficients and associated p-values between pairs of variables in A and B, controlling 
for the the variables in C.

Parameters
----------
A : array-like, shape (n, pa)

B : array-like, shape (n, pb) (optional)

C : array-like, shape (n, pc) (optional)

namesA : aray-like, shape(1,pa) names of the variables in columns of A (optional)

namesB : aray-like, shape(1,pb) names of the variables in columns of B (optional)

method : string value='Pearson' or 'Spearman' default='Spearman'
	if method='Pearson': linear partial correlation coefficients and assocated p-values are returned.
	if method='Spearman': linear partial rank correlaltion coefficients and associated p-values are returned. 

Returns
-------
Cor_coefs : pandas data frame, table of correlation coefffiecents

P_vals : pandas data frame, table of p-values
