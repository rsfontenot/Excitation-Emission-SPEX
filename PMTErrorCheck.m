function []=PMTErrorCheck(x,obj)
%This program peforms checks on the PMT response. 0 means debug mode, i.e.
%error is generated no matter the response. 1 means error means that VA
%will allow the program to proceed. obj is the object you wish to check.
if x==1
    tStart=tic;
while obj.BytesAvailable<2
    tStop=toc(tStart);
    if tStop>10
        error('Stuck in Infinite loop')
    end
end
sol=fread(obj,obj.BytesAvailable)';
if char(sol)=='VA'
elseif char(sol)=='BC'
    error('Bad Command')
elseif char(sol)=='BA'
    error('Bad Argument')
else
    error('Really bad argument')
end
else
    tStart=tic;
while obj.BytesAvailable<2
    tStop=toc(tStart);
    if tStop>30
        error('Stuck in Infinite loop')
    end
end
sol=fread(obj,obj.BytesAvailable)';
if char(sol)=='VA'
    error('Valid Response')
elseif char(sol)=='BC'
    error('Bad Command')
elseif char(sol)=='BA'
    error('Bad Argument')
else
    error('Really bad argument')
end
end