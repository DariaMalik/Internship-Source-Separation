% B = rescale(A,l,u,'InputMin',inmin,'InputMax',inmax) uses the formula
% 
% l + [(A-inmin)./(inmax-inmin)].*(u-l)
% 
% to scale the elements of an array A.
% 
% If l and u are not specified, then rescale uses the default values 0 and 1, respectively.
% 
% If the 'InputMin' name-value pair is not specified, then rescale sets its value to the default min(A(:)).
% 
% If the 'InputMax' name-value pair is not specified, then rescale sets its value to the default max(A(:)).