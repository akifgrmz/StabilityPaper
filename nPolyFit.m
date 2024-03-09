function yFitted=nPolyFit(x,y,n)
% this function fits a dataset to an n degree polynomial and returns the
% fitted y values
%n: poly degree
%x:function input
%y:function output

p=polyfit(x,y,n);

yFitted=zeros(1,length(x));
    
    for iP=1:n+1
        
        yFitted=p(iP).*x.^(n+1-iP)+yFitted;
        
    end
end