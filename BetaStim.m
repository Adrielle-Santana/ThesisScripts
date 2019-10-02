% Computa estatística para comparar os coeficientes EN para o beta de
% diferentes estímulos
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
%Moving average filter parameters
windowSize = 25; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
gain=20;

col=0.4;

ad={'VOTact.mat', 'VOTpass.mat', 'FormAct.mat', 'FormPass.mat'};

pos=4;

close all 

ax=0.06;

load(sprintf('psy1/SIG-matrices-all-%s',ad{pos}))
left=coef_psy(1:1750);
right=coef_psy(1751:end);

load(sprintf('psy2/SIG-matrices-all-%s',ad{pos}))
left2=coef_psy(1:1750);
right2=coef_psy(1751:end);

load(sprintf('psy3/SIG-matrices-all-%s',ad{pos}))
left3=coef_psy(1:1750);
right3=coef_psy(1751:end);

load(sprintf('psy4/SIG-matrices-all-%s',ad{pos}))
left4=coef_psy(1:1750);
right4=coef_psy(1751:end);

load(sprintf('psy5/SIG-matrices-all-%s',ad{pos}))
left5=coef_psy(1:1750);
right5=coef_psy(1751:end);

t=(0:(length(left)-1))/5000;
 
%plot(t,left,t,left2,'r')

figure(1)
subplot(2,3,1)
plot(t,left,'k')
hold on
title('Psy-Left')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf -ax ax])
f_left = filtfilt(b,a,left);
plot(t,gain*f_left,'Color', [0.9 col col])
hold off
grid

subplot(2,3,2)
plot(t,left2,'k')
hold on
title('Psy-Left2')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf -ax ax])
f_left2 = filtfilt(b,a,left2);
plot(t,gain*f_left2,'Color', [0.9 col col])
hold off
grid

subplot(2,3,3)
plot(t,left3,'k')
hold on
title('Psy-Left3')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf -ax ax])
f_left3 = filtfilt(b,a,left3);
plot(t,gain*f_left3,'Color', [0.9 col col])
hold off
grid

subplot(2,3,4)
plot(t,left4,'k')
hold on
title('Psy-Left4')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf -ax ax])
f_left4 = filtfilt(b,a,left4);
plot(t,gain*f_left4,'Color', [0.9 col col])
hold off
grid

subplot(2,3,5)
plot(t,left5,'k')
hold on
title('Psy-Left5')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf -ax ax])
f_left5 = filtfilt(b,a,left5);
plot(t,gain*f_left5,'Color', [0.9 col col])
hold off
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
subplot(2,3,1)
plot(t,right,'k')
hold on
title('Psy-Right')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf -ax ax])
f_right = filtfilt(b,a,left);
plot(t,gain*f_right,'Color', [0.9 col col])
hold off
grid

subplot(2,3,2)
plot(t,right2,'k')
hold on
title('Psy-Right2')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf -ax ax])
f_right2 = filtfilt(b,a,right2);
plot(t,gain*f_right2,'Color', [0.9 col col])
hold off
grid

subplot(2,3,3)
plot(t,right3,'k')
hold on
title('Psy-Right3')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf -ax ax])
f_right3 = filtfilt(b,a,right3);
plot(t,gain*f_right3,'Color', [0.9 col col])
hold off
grid

subplot(2,3,4)
plot(t,right4,'k')
hold on
title('Psy-Right4')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf -ax ax])
f_right4 = filtfilt(b,a,right4);
plot(t,gain*f_right4,'Color', [0.9 col col])
hold off
grid

subplot(2,3,5)
plot(t,right5,'k')
hold on
title('Psy-Right5')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf -ax ax])
f_right5 = filtfilt(b,a,right5);
plot(t,gain*f_right5,'Color', [0.9 col col])
hold off
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigSize=length(coef_phy)/2;

t=(0:sigSize-1)/5000;

figure(3)

subplot(1,2,1)
plot(t,coef_phy(1:sigSize),'k')
hold on
title('Phy-Left')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf -60 60])
y_phy_l = filtfilt(b,a,coef_phy(1:sigSize));
plot(t,gain*y_phy_l,'Color', [0.9 col col])
hold off
grid

subplot(1,2,2)
plot(t,coef_phy((sigSize+1):end),'k')
hold on
title('Phy-Right')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf -60 60])
y_phy_r = filtfilt(b,a,coef_phy((sigSize+1):end));
plot(t,gain*y_phy_r,'Color', [0.9 col col])
hold off
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
claro=mean([left,left2,left4,left5,right,right2,right4,right5],2);
amb=mean([left3,right3],2);

figure(4)

ax=0.05

subplot(1,2,1)
plot(t,claro,'k')
hold on
title('Claros')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf -ax ax])
y_phy_l = filtfilt(b,a,claro);
plot(t,gain*y_phy_l,'Color', [0.9 col col])
hold off
grid

subplot(1,2,2)
plot(t,amb,'k')
hold on
title('Amb')
ylabel('Elastic Net Coefficient')
xlabel('Time(s)')
axis([0 inf -ax ax])
y_phy_r = filtfilt(b,a,amb);
plot(t,gain*y_phy_r,'Color', [0.9 col col])
hold off
grid