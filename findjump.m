function [trans,found,test] = findjump(dU_c,thresh,init)
% [trans,found,test] = findjump(dU_c,thresh,init)
% Finds a statistical jump in a vector data. THe statistical threshold is
% mean(dU_c)+thresh*std(dU_c).
% INPUT:
% - dU_c : data vector
% - thresh : number of standard deviations the statistical step must exceed
% -init : positive interger, first element of the vector on which the search 
% can be conducted. If set to one, the statistical test won't be significative.
% 
% OUTPUT:
% - trans : index of the significative steps
% - found : logical. True if a significative step has been found
% - test : value or the threshold for significative test. 

idx = init;
found = 0;
trans = NaN;
while (found<1 && idx<numel(dU_c)-1)
   stddiff = std(dU_c(1:idx));
   test = mean(dU_c(1:idx)) + thresh*stddiff;
   idx=idx+1;
   test1 = dU_c(idx)>test;
       if test1
           trans = find(dU_c>test);
           found = 1;
       end
end
   
end