%This script loads and plots the EN coefficients chosen for the groups
%of AEP against the physical and psychophysical expected responses

% Adrielle C. Santana
%16-09-2019
clear all
close all

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%Moving average filter parameters
windowSize = 50; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
gain=20;
sel=[1,3,4,5,6,7,9,10];

%Uncomment the condition you want to analyze
%Feito com sequência 0,0.05,0.5,0.95 e 100% e número de estímulos
%Interpretabilidade ruim. Plotar para reportar que fiz, na Tese.

load('SIG-matrices-all-VOTpass-psyold.mat')
%load('SIG-matrices-all-VOTact-psyold.mat')
%load('SIG-matrices-all-FormPass-psyold.mat')
%load('SIG-matrices-all-FormAct-psyold.mat')

phy=coef_phy;
psy=coef_psy;

%For one subject in sel, use this
%phy=coefs_phy(sel,:);
%psy=coefs_psy(sel,:);

sigSize=length(phy)/2;

col=0.4
t=(0:sigSize-1)/5000;

ymin=-0.05; 
ymax=0.05;

subplot(2,2,1)
plot(t,phy(1:sigSize),'k')
hold on
title('Phy-Left')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf ymin ymax])
y_phy_l = filter(b,a,phy(1:sigSize));
plot(t,gain*y_phy_l,'Color', [0.9 col col])
hold off
grid

subplot(2,2,2)
plot(t,phy((sigSize+1):end),'k')
hold on
title('Phy-Right')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf ymin ymax])
y_phy_r = filter(b,a,phy((sigSize+1):end));
plot(t,gain*y_phy_r,'Color', [0.9 col col])
hold off
grid

subplot(2,2,3)
plot(t,psy(1:sigSize),'k')
hold on
title('Psy-Left')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf ymin ymax])
y_psy_l = filter(b,a,psy(1:sigSize));
plot(t,gain*y_psy_l,'Color', [0.9 col col])
hold off
grid

subplot(2,2,4)
plot(t,psy((sigSize+1):end),'k')
hold on
title('Psy-Right')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf ymin ymax])
y_psy_r = filter(b,a,psy((sigSize+1):end));
plot(t,gain*y_psy_r,'Color', [0.9 col col])
hold off
grid

orient(figure(1),'landscape')
print(figure(1),'Formact.pdf','-dpdf')