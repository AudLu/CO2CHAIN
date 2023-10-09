function [bedload_sort,discharge_sort]=sort_sameQ(bedload,discharge)
% [bedload_sort,discharge_sort]=sort_sameQ(bedload,discharge)
% sorts an input dataset into files of rising discharge.


[discharge_sort,discharge_index] = sort(discharge);
bedload_sort = bedload(discharge_index);
