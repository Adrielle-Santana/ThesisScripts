%This script loads and plots the EN coefficients chosen for the groups
%of AEP against the physical and psychophysical expected responses

% Adrielle C. Santana
%19-08-2019
clear all
close all

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%Moving average filter parameters
windowSize = 100; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
gain=15;
sel=[1,3,4,5,6,7,9,10];

%groups
g=7; %7 é o melhor

%Uncomment the condition you want to analyze

%load(sprintf('Todos/SIG-matrices-all-VOTpass-sqrt-%d-groups.mat',g))
%load(sprintf('Todos/SIG-matrices-all-VOTact-sqrt-%d-groups.mat',g))
%load(sprintf('Todos/SIG-matrices-all-FormPass-sqrt-%d-groups.mat',g))
load(sprintf('Todos/SIG-matrices-all-FormAct-sqrt-%d-groups.mat',g))

phy=mean(coefs_phy(sel,:));
psy=mean(coefs_psy(sel,:));

%For one subject in sel, use this
%phy=coefs_phy(sel,:);
%psy=coefs_psy(sel,:);

sigSize=length(phy)/2;

t=(0:sigSize-1)/5000;

ymin=-10; 
ymax=10;

subplot(2,2,1)
plot(t,phy(1,1:sigSize))
hold on
title('Phy-Left')
ylabel('uV')
xlabel('Time(s)')
axis([0 inf ymin ymax])
y_phy_l = filter(b,a,phy(1,1:sigSize));
plot(t,gain*y_phy_l,'r')
hold off
grid

subplot(2,2,2)
plot(t,phy(1,(sigSize+1):end))
hold on
title('Phy-Right')
ylabel('uV')
xlabel('Time(s)')
axis([0 inf ymin ymax])
y_phy_r = filter(b,a,phy(1,(sigSize+1):end));
plot(t,gain*y_phy_r,'r')
hold off
grid

subplot(2,2,3)
plot(t,psy(1,1:sigSize))
hold on
title('Psy-Left')
ylabel('uV')
xlabel('Time(s)')
axis([0 inf ymin ymax])
y_psy_l = filter(b,a,psy(1,1:sigSize));
plot(t,gain*y_psy_l,'r')
hold off
grid

subplot(2,2,4)
plot(t,psy(1,(sigSize+1):end))
hold on
title('Psy-Right')
ylabel('uV')
xlabel('Time(s)')
axis([0 inf ymin ymax])
y_psy_r = filter(b,a,psy(1,(sigSize+1):end));
plot(t,gain*y_psy_r,'r')
hold off
grid

c=corrcoef(y_psy_l,y_psy_r);
cPsL_PsR=c(1,2);

c=corrcoef(y_psy_l,y_phy_l);
cPsL_PhL=c(1,2);

c=corrcoef(y_psy_l,y_phy_r);
cPsL_PhR=c(1,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=corrcoef(y_phy_l,y_phy_r);
cPhL_PhR=c(1,2);

c=corrcoef(y_phy_l,y_psy_r);
cPhL_PsR=c(1,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=corrcoef(y_phy_r,y_psy_r);
cPhL_PhR=c(1,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Mean AIC
phyLAIC=mean(AIC_coef.phy_left)
phyRAIC=mean(AIC_coef.phy_right)
psyLAIC=mean(AIC_coef.psy_left)
psyRAIC=mean(AIC_coef.psy_right)
AICall=mean([phyLAIC,phyRAIC,psyLAIC,psyRAIC])
