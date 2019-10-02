%Notch filter
function out = notch60(sig, rate)

% sig: array with signal or matrix of signals to be filtered
% rate: sampling rate 

tstep = 1/rate;
Fc = 60*tstep;
Bandwidth=4;

d = exp(-2*pi*(Bandwidth/2)*tstep);
b = (1 + d*d)*cos(2*pi*Fc);
a0 = 1;
a1 = -b;
a2 = d*d;
a = (1 + d*d)/2;
b0 = 1;
b1 = -2*cos(2*pi*Fc);
b2 = 1;

a=[a0,a1,a2];
b=[b0,b1,b2];

%aux1 = filter(b,a,sig,[],2);
%aux2 = filter(b,a,fliplr(aux1),[],2);
%out = fliplr(aux2);

for i=1:size(sig,1)
    sig(i,:)=filtfilt(b,a,sig(i,:));
end

out=sig;

return