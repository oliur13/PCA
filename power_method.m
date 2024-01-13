%Inputs matrix and guess eigenvector, returns largest eigenvalue and
%resulting eigenvector
%Engr 592 Winter 2018
%Original code developed by rayryeng at Stackoverflow
%http://stackoverflow.com/questions/29198277/power-method-in-matlab

function [m,x]=power_method(A,x)
m=0;
n=length(x);
tol=1e-10;
while(1)
     mold = m;
     x=A*x;
     m=max(x);
     x=x/m;
     if abs(m-mold) < tol
         x = x/(norm(x,2));
         break;
     end
end
end