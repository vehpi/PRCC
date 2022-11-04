#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 09:36:05 2022

@author: Vehpi Yildirim, v.yildirim@amsterdamumc.nl, https://github.com/vehpi

partial_corr(A)
Returns the sample linear partial or rank correlation coefficients and associated p-values between pairs of variables in A, controlling 
for the remaining variables in A.
partial_corr(A,B)
Returns the sample linear partial or rank correlation coefficients and associated p-values between pairs of variables in A and B, controlling 
for the remaining variables in A.
partial_corr(A,B,C)
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

"""
def partial_corr(A,B=None,C=None,method='Spearman',namesA=None,namesB=None,significance_th=0.05,labelA='A',labelB='B',box_plt=False):
	
	def box_plots(A,B,namesA,namesB,labelA,labelB,box_pl):
		fig1=plt.figure(np.random.randint(low=0,high=1000000))
		ax1=fig1.add_subplot(111)
		if box_pl=='log':
			ax1.boxplot(np.log10(A),vert=False,labels=namesA)
			trans='log transformed'
		elif box_pl=='standardized':
			ax1.boxplot(A/np.mean(A,axis=0),vert=False,labels=namesA)
			trans='standardized'
		elif box_pl=='normal':
			ax1.boxplot(A,vert=False,labels=namesA)
			trans=''
		ax1.set_xlabel(f'{trans} value')
		fig1.suptitle(f'Box plots of {trans} {labelA}')
		
		if B is not None:
			fig2=plt.figure(np.random.randint(low=0,high=1000000))
			ax2=fig2.add_subplot(111)
			
			if box_pl=='log':
				ax2.boxplot(np.log10(B),vert=False,labels=namesB)
				trans='log transformed'
			elif box_pl=='standardized':
				ax2.boxplot(B/np.mean(B,axis=0),vert=False,labels=namesB)
				trans='standardized'
			elif box_pl=='normal':
				ax2.boxplot(B,vert=False,labels=namesB)
				trans=''
			ax2.set_xlabel(f'{trans} value')
			fig2.suptitle(f'Box plots of {trans} {labelB}')
			
	if namesA:
		rownames=namesA
	else:
		rownames=[f'r{i}' for i in range(A.shape[1])]
	
	if B is not None:
		if namesB:
			colnames=namesB
		else:
			colnames=[f'c{i}' for i in range(B.shape[1])]
	else:
		if namesA:
			colnames=namesA
		else:
			colnames=[f'c{i}' for i in range(A.shape[1])]
			
	if C is not None and B is not None:
		if box_plt:
			box_plots(A=A,B=B,namesA=namesA,namesB=namesB,labelA=labelA,labelB=labelB,box_pl=box_plt)
			
		if np.shape(A)[0]==np.shape(B)[0]==np.shape(C)[0]:
			if method=='Spearman':
				A=rankdata(A,axis=0)
				B=rankdata(B,axis=0)
				C=rankdata(C,axis=0)
			A=np.asarray(A)
			B=np.asarray(B)
			C=np.asarray(C)
			n = np.shape(A)[0]
			pa = A.shape[1]
			pb = B.shape[1]
			pc = C.shape[1]
			Cor_coefs = np.zeros((pa, pb), dtype=float)
			P_vals = np.zeros((pa, pb), dtype=float)
			
			C=np.append(np.ones([n,1]),C,axis=1)
			for i in range(pa):
				for j in range(pb):
					beta_a = linalg.lstsq(C, A[:, i])[0]
					beta_b = linalg.lstsq(C, B[:, j])[0]
					
					res_a = A[:, i] - C.dot( beta_a)
					res_b = B[:, j] - C.dot(beta_b)
					
					res_obj=stats.pearsonr(res_a, res_b)
					corr=res_obj[0]
					#pval=res_obj[1] 
					# p-value directly calculated by stats.pearsonr, which choses df=N-2
					
					df = len(A) - 2 - (np.shape(C)[1])
					t =abs( corr * np.sqrt(df/(1 - corr**2)))
					pval=2*tdist.sf(t, df) # p-value corrected df for parameters varied, which uses df=N-2-pc, where p number of columns of A
					
					Cor_coefs[i,j] = corr
					P_vals[i,j] = pval
					
			last_d=np.nan*np.zeros(np.shape(Cor_coefs))
			last_d[P_vals<significance_th]=Cor_coefs[P_vals<significance_th]
			last_table=pd.DataFrame(data=last_d,index=rownames,columns=colnames)
			last_table.columns.name='Sig. r-values'
			
			Cor_coefs=pd.DataFrame(data=Cor_coefs,index=rownames,columns=colnames)
			P_vals=pd.DataFrame(data=P_vals,index=rownames,columns=colnames)
			Cor_coefs.columns.name='r-values'
			P_vals.columns.name='P-values'
			
			return Cor_coefs, P_vals, last_table
		else:
			print('Row numbers of A, B and C do not match')
	
	elif B is not None:
		if box_plt:
			box_plots(A=A,B=B,namesA=namesA,namesB=namesB,labelA=labelA,labelB=labelB,box_pl=box_plt)
		if method=='Spearman':
			A=rankdata(A,axis=0)
			B=rankdata(B,axis=0)
		A=np.asarray(A)
		B=np.asarray(B)
		if np.shape(A)[0]==np.shape(B)[0]:
			n = np.shape(A)[0]
			pa = A.shape[1]
			pb = B.shape[1]
			Cor_coefs = np.zeros((pa, pb), dtype=float)
			P_vals = np.zeros((pa, pb), dtype=float)
			
			A_extd=np.append(np.ones([n,1]),A,axis=1)
			
			for i in range(pa):
				idx = np.ones(1+pa, dtype=bool)
				idx[i+1] = False
				for j in range(pb):
					beta_a = linalg.lstsq(A_extd[:,idx], A[:, i])[0]
					beta_b = linalg.lstsq(A_extd[:,idx], B[:, j])[0]
					
					res_a = A[:, i] - A_extd[:,idx].dot( beta_a)
					res_b = B[:, j] - A_extd[:,idx].dot(beta_b)
					res_obj=stats.pearsonr(res_a, res_b)
					corr=res_obj[0]
					#pval=res_obj[1] 
					# p-value directly calculated by stats.pearsonr, which choses df=N-2
					
					df = len(A) - 2 - (np.shape(A)[1]-1)
					t =abs( corr * np.sqrt(df/(1 - corr**2)))
					pval=2*tdist.sf(t, df) # p-value corrected df for parameters varied, which uses df=N-2-(pa-1), where p number of columns of A
					
					Cor_coefs[i,j] = corr
					P_vals[i,j] = pval
			
			last_d=np.nan*np.zeros(np.shape(Cor_coefs))
			last_d[P_vals<significance_th]=Cor_coefs[P_vals<significance_th]
			last_table=pd.DataFrame(data=last_d,index=rownames,columns=colnames)
			last_table.columns.name='Sig. r-values'
			
			Cor_coefs=pd.DataFrame(data=Cor_coefs,index=rownames,columns=colnames)
			P_vals=pd.DataFrame(data=P_vals,index=rownames,columns=colnames)
			Cor_coefs.columns.name='r-values'
			P_vals.columns.name='P-values'
			
			return Cor_coefs, P_vals,last_table
		else:
			print('Row numbers of A and B does not match')
	else:
		if box_plt:
			box_plots(A=A,B=None,namesA=namesA,namesB=namesB,labelA=labelA,labelB=labelB,box_pl=box_plt)
		if method=='Spearman':
			A=rankdata(A,axis=0)
		A=np.asarray(A)
		[n,pa] = A.shape
		Cor_coefs = np.zeros((pa, pa), dtype=float)
		P_vals = np.zeros((pa, pa), dtype=float)
		A_extd=np.append(np.ones([n,1]),A,axis=1)
		for i in range(pa):
			Cor_coefs[i, i] = 1
			P_vals[i, i] = 0
			for j in range(i+1,pa):
				idx = np.ones(1+pa, dtype=bool)
				idx[i+1] = False
				idx[j+1] = False
				
				beta_i = linalg.lstsq(A_extd[:, idx], A[:, j])[0]
				beta_j = linalg.lstsq(A_extd[:, idx], A[:, i])[0]
	
				res_j = A[:, j] - A_extd[:, idx].dot( beta_i)
				res_i = A[:, i] - A_extd[:, idx].dot(beta_j)
				
				res_obj=stats.pearsonr(res_i, res_j)
				corr=res_obj[0]
				#pval=res_obj[1] 
				# p-value directly calculated by stats.pearsonr, which chooses df=N-2
				
				df = len(A) - 2 - (np.shape(A)[1]-2)
				t =abs( corr * np.sqrt(df/(1 - corr**2)))
				pval=2*tdist.sf(t, df) # p-value corrected df for parameters varied, which uses df=N-2-(p-2), where p number of columns of A
				
				Cor_coefs[i,j] = Cor_coefs[j,i] = corr
				P_vals[i,j] = P_vals[j,i] = pval
				
		
		last_d=np.nan*np.zeros(np.shape(Cor_coefs))
		last_d[P_vals<significance_th]=Cor_coefs[P_vals<significance_th]
		last_table=pd.DataFrame(data=last_d,index=rownames,columns=colnames)
		last_table.columns.name='Sig. r-values'
		
		Cor_coefs=pd.DataFrame(data=Cor_coefs,index=rownames,columns=colnames)
		P_vals=pd.DataFrame(data=P_vals,index=rownames,columns=colnames)
		Cor_coefs.columns.name='r-values'
		P_vals.columns.name='P-values'
		
		return Cor_coefs, P_vals, last_table