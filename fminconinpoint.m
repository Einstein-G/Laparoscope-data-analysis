function [x,fval,exitflag,output,lambda,grad,hessian] = fminconinpoint(x0,lb,ub)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = optimoptions('fmincon');
%% Modify options setting
options = optimoptions(options,'Display', 'off','TolFun',1E-1);
options = optimoptions(options,'PlotFcns', {  @optimplotx @optimplotfval });
options = optimoptions(options,'Algorithm', 'interior-point');
[x,fval,exitflag,output,lambda,grad,hessian] = ...
fmincon(@weber,x0,[],[],[],[],lb,ub,[],options);
display(x);display(fval);