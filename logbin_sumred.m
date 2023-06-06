function [x_med,y_med] = logbin_sumred(xval,yval,option,bintype,nobins)
%[x_med,y_med] = logbin_sumred(xval,yval,option,bintype,nobins)
%reduced version of logbin_sum
% Inputs: xval: column vector of independent variable
%         yval: column vector of dependent variable
%         option:  0 default
%                  1 bins without content are removed
%                  2 zero values of the x-vector are removed
%                  3 both 1 and 2
%         bintype: if bintype = 1 linearly distributed bins are used
%                  otherwise logarithmically distributed bins are used
%                  (default)
%         nobins:  number of bins
%                  if no number is given, the squareroot of the number of
%                  values is used as default.
%
% Outputs: bin median for independent and dependant variable


% establish size

[m,n]=size(xval);
[u,v]=size(yval);

if n~=1 && v~=1 && m~=u
	error('input files are not the same size')
end

% set default values if not all arguments are given
if nargin < 4
    bintype = 0;
end
if nargin < 3
    option = 0;
end

	
[yval,xval] = sort_sameQ(yval,xval);

% remove zero values and nan in y vector

if option == 2 || option == 3
    dummyx = zeros(m,n);
    dummyy = zeros(u,v);

    k = 0;
    for i=1:m
        if yval(i) ~= 0 && xval(i) ~= 0 && isnan(yval(i)) == 0
            k = k+1;
            dummyy(k,1) = yval(i);
            dummyx(k,1) = xval(i);
        end
    end
    yval = dummyy(1:k,1);
    xval = dummyx(1:k,1);

    clear dummyx
    clear dummyy
else
    k = m;
end
if nargin < 5
    % set default to squareroot of the number of values.
    nobins = round(sqrt(k));
end

% set bin boundaries
%   binbound = 10.^(linspace(log10(min(xval)),log10(max(xval)),nobins+1));
if bintype ~= 1
      binbound = 10.^(linspace(log10(min(xval)),log10(max(xval)),nobins+1));
% if max(xval)/min(xval) > 20 , 
% binbound = logspace(1,7,nobins+1) ;
% else
%      binbound = logspace(2.5,3.5,nobins+1) ;
% end
 else
    binbound = linspace(min(xval),max(xval),nobins+1);
end

binbound = binbound';

% initialise parameters

counterbin = zeros(nobins,1);		% counts number of data points in this bin
currentbin = 1;		% counts number of bins for this discharge
zerocount = 0;


x_med = zeros(nobins,1);
y_med = zeros(nobins,1);


maxind = zeros(nobins,1);
minind = zeros(nobins,1);


for i=1:nobins
    smallerval = find(xval<binbound(i+1));
    largerval = find(xval>=binbound(i));
    if isempty(smallerval), smallerval=-1;end
    if isempty(largerval), largerval=0;end
       
    if i~=nobins
        maxind(i) = max(smallerval);
    else
        maxind(i) = k;
    end
    if i~=1
        minind(i) = min(largerval);
    else
        minind(i) = 1;
    end
    if maxind(i) >= minind(i) && minind(i)>0
        counterbin(i) = maxind(i)-minind(i)+1;
       
        x_med(i,1) = median(xval(minind(i):maxind(i)));
        y_med(i,1) = median(yval(minind(i):maxind(i)));
    
        perc=maxind(i)-minind(i);
         
    else
        counterbin(i) = 0;
        zerocount = zerocount+1;
        x_med(currentbin) = median([binbound(currentbin),binbound(currentbin+1)]);
    end
end

% remove bins without entries

if (option == 1 || option == 3) && zerocount > 0
    
    counterbin2 = zeros(nobins-zerocount,1);
    
    x_med2 = zeros(nobins-zerocount,1);
    y_med2 = zeros(nobins-zerocount,1);
 
        minind2 = zeros(nobins-zerocount,1);
    maxind2 = zeros(nobins-zerocount,1);
    j = 0;
    
    for i=1:nobins
        if counterbin(i) ~= 0
            j = j+1;
            counterbin2(j) = counterbin(i);
           
            x_med2(j,1) = x_med(i,1);
            y_med2(j,1) = y_med(i,1);
           
            minind2(j,1) = minind(i,1);
            maxind2(j,1) = maxind(i,1);
            
        end
    end
    counterbin = counterbin2;
    
    x_med = x_med2;
    y_med = y_med2;
    
    minind = minind2;
    maxind = maxind2;
end