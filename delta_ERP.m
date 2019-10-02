%Reproduce regression between delta_ERP and slope of psychometric curve made
%by Bidelman 2017
close all
clear all

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%Subject
s=[1,2,3,4,5,6,7,8,9,10,11];
nb_subj=length(s);

bt(1:2,1:nb_subj)=0;

med=[1,2,4,5];

mN1=900;

aux=zeros(2,nb_subj);
aux2=zeros(2,nb_subj);
aux3=zeros(2,nb_subj);
aux4=zeros(2,nb_subj);
aux5=zeros(2,nb_subj);
aux6=zeros(2,nb_subj);
aux7=zeros(2,nb_subj);
aux8=zeros(2,nb_subj);

box1=zeros(2,nb_subj);
box2=zeros(2,nb_subj);
box3=zeros(2,nb_subj);
box4=zeros(2,nb_subj);

p1=zeros(nb_subj,5);
picos=zeros(nb_subj,5);
picos_val=zeros(nb_subj,5);
vales=zeros(nb_subj,5);
vales_val=zeros(nb_subj,5);
diffsp1=zeros(nb_subj,5);
diffs=zeros(nb_subj,5);
deltas=zeros(nb_subj,1);
left=struct('p1',p1,'picos',picos,'vales',vales,'picos_val',picos_val,'vales_val',vales_val,'diffsp1',diffsp1,'diffs',diffs,'deltas',deltas);
right=struct('p1',p1,'picos',picos,'vales',vales,'picos_val',picos_val,'vales_val',vales_val,'diffsp1',diffsp1,'diffs',diffs,'deltas',deltas);

ampFormpass=struct('left',left,'right',right);
ampFormact=struct('left',left,'right',right);
ampVOTact=struct('left',left,'right',right);
ampVOTpass=struct('left',left,'right',right);

for subj=s

        load(sprintf('../Sujeito%d/VOT/Identification/beta.mat',subj))
        bt(1,subj)=beta;
        
        load(sprintf('../Sujeito%d/Formantes/Identification/beta.mat',subj))
        bt(2,subj)=beta;

       for stim=1:5
        
        load(sprintf('../Sujeito%d/VOT/Passivo/ChanStimICA%d',subj,stim))
        
        ampVOTpass.right.vales(subj,stim)=(find(-Xica_DWT_rec(2,mN1:1500)==min(-Xica_DWT_rec(2,mN1:1500))))+(mN1-1);
        ampVOTpass.right.p1(subj,stim)=(find(-Xica_DWT_rec(2,550:1050)==max(-Xica_DWT_rec(2,550:1050))))+549;
        ampVOTpass.right.picos(subj,stim)=(find(-Xica_DWT_rec(2,1100:2000)==max(-Xica_DWT_rec(2,1100:2000))))+1099;
        ampVOTpass.right.diffsp1(subj,stim)= -Xica_DWT_rec(2,ampVOTpass.right.p1(subj,stim))-(-Xica_DWT_rec(2,ampVOTpass.right.vales(subj,stim)));
        ampVOTpass.right.diffs(subj,stim)= -Xica_DWT_rec(2,ampVOTpass.right.picos(subj,stim))-(-Xica_DWT_rec(2,ampVOTpass.right.vales(subj,stim)));
        ampVOTpass.right.picos_val(subj,stim)=-Xica_DWT_rec(2,ampVOTpass.right.picos(subj,stim));
        ampVOTpass.right.vales_val(subj,stim)=-Xica_DWT_rec(2,ampVOTpass.right.vales(subj,stim));
        
        ampVOTpass.left.vales(subj,stim)=(find(-Xica_DWT_rec(17,mN1:1500)==min(-Xica_DWT_rec(17,mN1:1500))))+(mN1-1);
        ampVOTpass.left.p1(subj,stim)=(find(-Xica_DWT_rec(17,550:1050)==max(-Xica_DWT_rec(17,550:1050))))+549;
        ampVOTpass.left.picos(subj,stim)=(find(-Xica_DWT_rec(17,1100:2000)==max(-Xica_DWT_rec(17,1100:2000))))+1099;
        ampVOTpass.left.diffsp1(subj,stim)= -Xica_DWT_rec(17,ampVOTpass.left.p1(subj,stim))-(-Xica_DWT_rec(17,ampVOTpass.left.vales(subj,stim)));
        ampVOTpass.left.diffs(subj,stim)= -Xica_DWT_rec(17,ampVOTpass.left.picos(subj,stim))-(-Xica_DWT_rec(17,ampVOTpass.left.vales(subj,stim)));
        ampVOTpass.left.picos_val(subj,stim)=-Xica_DWT_rec(17,ampVOTpass.left.picos(subj,stim));
        ampVOTpass.left.vales_val(subj,stim)=-Xica_DWT_rec(17,ampVOTpass.left.vales(subj,stim));
        
       
        load(sprintf('../Sujeito%d/VOT/Ativo/ChanStimICA%d',subj,stim))
        
        ampVOTact.right.vales(subj,stim)=(find(-Xica_DWT_rec(2,mN1:1500)==min(-Xica_DWT_rec(2,mN1:1500))))+(mN1-1);
        ampVOTact.right.p1(subj,stim)=(find(-Xica_DWT_rec(2,550:1050)==max(-Xica_DWT_rec(2,550:1050))))+549;
        ampVOTact.right.picos(subj,stim)=(find(-Xica_DWT_rec(2,1100:2000)==max(-Xica_DWT_rec(2,1100:2000))))+1099;
        ampVOTact.right.diffsp1(subj,stim)= -Xica_DWT_rec(2,ampVOTact.right.p1(subj,stim))-(-Xica_DWT_rec(2,ampVOTact.right.vales(subj,stim)));
        ampVOTact.right.diffs(subj,stim)= -Xica_DWT_rec(2,ampVOTact.right.picos(subj,stim))-(-Xica_DWT_rec(2,ampVOTact.right.vales(subj,stim)));
        ampVOTact.right.picos_val(subj,stim)=-Xica_DWT_rec(2,ampVOTact.right.picos(subj,stim));
        ampVOTact.right.vales_val(subj,stim)=-Xica_DWT_rec(2,ampVOTact.right.vales(subj,stim));
        
        ampVOTact.left.vales(subj,stim)=(find(-Xica_DWT_rec(17,mN1:1500)==min(-Xica_DWT_rec(17,mN1:1500))))+(mN1-1);
        ampVOTact.left.p1(subj,stim)=(find(-Xica_DWT_rec(17,550:1050)==max(-Xica_DWT_rec(17,550:1050))))+549;
        ampVOTact.left.picos(subj,stim)=(find(-Xica_DWT_rec(17,1100:2000)==max(-Xica_DWT_rec(17,1100:2000))))+1099;
        ampVOTact.left.diffsp1(subj,stim)= -Xica_DWT_rec(17,ampVOTact.left.p1(subj,stim))-(-Xica_DWT_rec(17,ampVOTact.left.vales(subj,stim)));
        ampVOTact.left.diffs(subj,stim)= -Xica_DWT_rec(17,ampVOTact.left.picos(subj,stim))-(-Xica_DWT_rec(17,ampVOTact.left.vales(subj,stim)));
        ampVOTact.left.picos_val(subj,stim)=-Xica_DWT_rec(17,ampVOTact.left.picos(subj,stim));
        ampVOTact.left.vales_val(subj,stim)=-Xica_DWT_rec(17,ampVOTact.left.vales(subj,stim));
    
        
        load(sprintf('../Sujeito%d/Formantes/Passivo/ChanStimICA%d',subj,stim))
       
        ampFormpass.right.vales(subj,stim)=(find(-Xica_DWT_rec(2,853:1500)==min(-Xica_DWT_rec(2,853:1500))))+(853-1);
        ampFormpass.right.p1(subj,stim)=(find(-Xica_DWT_rec(2,550:1050)==max(-Xica_DWT_rec(2,550:1050))))+549;
        ampFormpass.right.picos(subj,stim)=(find(-Xica_DWT_rec(2,1100:2000)==max(-Xica_DWT_rec(2,1100:2000))))+1099;
        ampFormpass.right.diffsp1(subj,stim)= -Xica_DWT_rec(2,ampFormpass.right.p1(subj,stim))-(-Xica_DWT_rec(2,ampFormpass.right.vales(subj,stim)));
        ampFormpass.right.diffs(subj,stim)= -Xica_DWT_rec(2,ampFormpass.right.picos(subj,stim))-(-Xica_DWT_rec(2,ampFormpass.right.vales(subj,stim)));
        ampFormpass.right.picos_val(subj,stim)=-Xica_DWT_rec(2,ampFormpass.right.picos(subj,stim));
        ampFormpass.right.vales_val(subj,stim)=-Xica_DWT_rec(2,ampFormpass.right.vales(subj,stim));
        
        ampFormpass.left.vales(subj,stim)=(find(-Xica_DWT_rec(17,853:1500)==min(-Xica_DWT_rec(17,853:1500))))+(853-1);
        ampFormpass.left.p1(subj,stim)=(find(-Xica_DWT_rec(17,550:1050)==max(-Xica_DWT_rec(17,550:1050))))+549;
        ampFormpass.left.picos(subj,stim)=(find(-Xica_DWT_rec(17,1100:2000)==max(-Xica_DWT_rec(17,1100:2000))))+1099;
        ampFormpass.left.diffsp1(subj,stim)= -Xica_DWT_rec(17,ampFormpass.left.p1(subj,stim))-(-Xica_DWT_rec(17,ampFormpass.left.vales(subj,stim)));
        ampFormpass.left.diffs(subj,stim)= -Xica_DWT_rec(17,ampFormpass.left.picos(subj,stim))-(-Xica_DWT_rec(17,ampFormpass.left.vales(subj,stim)));
        ampFormpass.left.picos_val(subj,stim)=-Xica_DWT_rec(17,ampFormpass.left.picos(subj,stim));
        ampFormpass.left.vales_val(subj,stim)=-Xica_DWT_rec(17,ampFormpass.left.vales(subj,stim));

              
        load(sprintf('../Sujeito%d/Formantes/Ativo/ChanStimICA%d',subj,stim))
      
        ampFormact.right.vales(subj,stim)=(find(-Xica_DWT_rec(2,mN1:1500)==min(-Xica_DWT_rec(2,mN1:1500))))+(mN1-1);
        ampFormact.right.p1(subj,stim)=(find(-Xica_DWT_rec(2,550:1050)==max(-Xica_DWT_rec(2,550:1050))))+549;
        ampFormact.right.picos(subj,stim)=(find(-Xica_DWT_rec(2,1100:2000)==max(-Xica_DWT_rec(2,1100:2000))))+1099;
        ampFormact.right.diffsp1(subj,stim)= -Xica_DWT_rec(2,ampFormact.right.p1(subj,stim))-(-Xica_DWT_rec(2,ampFormact.right.vales(subj,stim)));
        ampFormact.right.diffs(subj,stim)= -Xica_DWT_rec(2,ampFormact.right.picos(subj,stim))-(-Xica_DWT_rec(2,ampFormact.right.vales(subj,stim)));
        ampFormact.right.picos_val(subj,stim)=-Xica_DWT_rec(2,ampFormact.right.picos(subj,stim));
        ampFormact.right.vales_val(subj,stim)=-Xica_DWT_rec(2,ampFormact.right.vales(subj,stim));
        
        ampFormact.left.vales(subj,stim)=(find(-Xica_DWT_rec(17,mN1:1500)==min(-Xica_DWT_rec(17,mN1:1500))))+(mN1-1);
        ampFormact.left.p1(subj,stim)=(find(-Xica_DWT_rec(17,550:1050)==max(-Xica_DWT_rec(17,550:1050))))+549;
        ampFormact.left.picos(subj,stim)=(find(-Xica_DWT_rec(17,1100:2000)==max(-Xica_DWT_rec(17,1100:2000))))+1099;
        ampFormact.left.diffsp1(subj,stim)= -Xica_DWT_rec(17,ampFormact.left.p1(subj,stim))-(-Xica_DWT_rec(17,ampFormact.left.vales(subj,stim)));
        ampFormact.left.diffs(subj,stim)= -Xica_DWT_rec(17,ampFormact.left.picos(subj,stim))-(-Xica_DWT_rec(17,ampFormact.left.vales(subj,stim)));
        ampFormact.left.picos_val(subj,stim)=-Xica_DWT_rec(17,ampFormact.left.picos(subj,stim));
        ampFormact.left.vales_val(subj,stim)=-Xica_DWT_rec(17,ampFormact.left.vales(subj,stim));
        
        end
        
        ampVOTpass.right.deltas(subj)=mean(ampVOTpass.right.diffs(subj,med))-ampVOTpass.right.diffs(subj,3);
        ampVOTpass.left.deltas(subj)=mean(ampVOTpass.left.diffs(subj,med))-ampVOTpass.left.diffs(subj,3);
        aux(1,subj)=mean(ampVOTpass.left.diffs(subj,med));
        aux2(1,subj)=ampVOTpass.left.diffs(subj,3);
        aux(2,subj)=mean(ampVOTpass.right.diffs(subj,med));
        aux2(2,subj)=ampVOTpass.right.diffs(subj,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        box1(1,subj)=mean(ampVOTpass.left.diffs(subj,1:5));
        box1(2,subj)=mean(ampVOTpass.right.diffs(subj,1:5));
        
        ampVOTact.right.deltas(subj)=mean(ampVOTact.right.diffs(subj,med))-ampVOTact.right.diffs(subj,3);
        ampVOTact.left.deltas(subj)=mean(ampVOTact.left.diffs(subj,med))-ampVOTact.left.diffs(subj,3);
        aux3(1,subj)=mean(ampVOTact.left.diffs(subj,med));
        aux4(1,subj)=ampVOTact.left.diffs(subj,3);
        aux3(2,subj)=mean(ampVOTact.right.diffs(subj,med));
        aux4(2,subj)=ampVOTact.right.diffs(subj,3);
        
        box2(1,subj)=mean(ampVOTact.left.diffs(subj,1:5));
        box2(2,subj)=mean(ampVOTact.right.diffs(subj,1:5));
        
        ampFormpass.right.deltas(subj)=mean(ampFormpass.right.diffs(subj,med))-ampFormpass.right.diffs(subj,3);
        ampFormpass.left.deltas(subj)=mean(ampFormpass.left.diffs(subj,med))-ampFormpass.left.diffs(subj,3);
        aux5(1,subj)=mean(ampFormpass.left.diffs(subj,med));
        aux6(1,subj)=ampFormpass.left.diffs(subj,3);
        aux5(2,subj)=mean(ampFormpass.right.diffs(subj,med));
        aux6(2,subj)=ampFormpass.right.diffs(subj,3);
        
        box3(1,subj)=mean(ampFormpass.left.diffs(subj,1:5));
        box3(2,subj)=mean(ampFormpass.right.diffs(subj,1:5));
        
        ampFormact.right.deltas(subj)=mean(ampFormact.right.diffs(subj,med))-ampFormact.right.diffs(subj,3);
        ampFormact.left.deltas(subj)=mean(ampFormact.left.diffs(subj,med))-ampFormact.left.diffs(subj,3);
        aux7(1,subj)=mean(ampFormact.left.diffs(subj,med));
        aux8(1,subj)=ampFormact.left.diffs(subj,3);
        aux7(2,subj)=mean(ampFormact.right.diffs(subj,med));
        aux8(2,subj)=ampFormact.right.diffs(subj,3);
        
        box4(1,subj)=mean(ampFormact.left.diffs(subj,1:5));
        box4(2,subj)=mean(ampFormact.right.diffs(subj,1:5));
end

%Normalize values of deltas
ampFormact.left.deltas=ampFormact.left.deltas/rssq(ampFormact.left.deltas);
ampFormpass.left.deltas=ampFormpass.left.deltas/rssq(ampFormpass.left.deltas);
ampVOTact.left.deltas=ampVOTact.left.deltas/rssq(ampVOTact.left.deltas);
ampVOTpass.left.deltas=ampVOTpass.left.deltas/rssq(ampVOTpass.left.deltas);

ampFormact.right.deltas=ampFormact.right.deltas/rssq(ampFormact.right.deltas);
ampFormpass.right.deltas=ampFormpass.right.deltas/rssq(ampFormpass.right.deltas);
ampVOTact.right.deltas=ampVOTact.right.deltas/rssq(ampVOTact.right.deltas);
ampVOTpass.right.deltas=ampVOTpass.right.deltas/rssq(ampVOTpass.right.deltas);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
b=polyfit(ampVOTpass.left.deltas,bt(1,:)',1);
subplot(1,2,1)
plot(ampVOTpass.left.deltas,bt(1,:),'o')
hold on
plot(ampVOTpass.left.deltas,polyval(b,ampVOTpass.left.deltas),'red')
hold off
xlabel('Deltas (uV)')
ylabel('Psychometric Curve Slope')
resid=bt(1,:)'-polyval(b,ampVOTpass.left.deltas);
SSresid = sum(resid.^2);
SStotal = (length(bt(1,:))-1) * var(bt(1,:));
R2 = 1 - (SSresid/SStotal);% * (length(beta(1,:))-1)/(length(beta(1,:))-length(b)))
title(sprintf('VOT Passive Left R2: %.4f', R2))

format long
b=polyfit(ampVOTpass.right.deltas,bt(1,:)',1);
subplot(1,2,2)
plot(ampVOTpass.right.deltas,bt(1,:),'o')
hold on
plot(ampVOTpass.right.deltas,polyval(b,ampVOTpass.right.deltas),'red')
hold off
xlabel('Deltas (uV)')
ylabel('Psychometric Curve Slope')
resid=bt(1,:)'-polyval(b,ampVOTpass.right.deltas);
SSresid = sum(resid.^2);
SStotal = (length(bt(1,:))-1) * var(bt(1,:));
R2 = 1 - (SSresid/SStotal);% * (length(beta(1,:))-1)/(length(beta(1,:))-length(b)))
title(sprintf('VOT Passive Right R2: %.4f', R2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
format long
b=polyfit(ampVOTact.left.deltas,bt(1,:)',1);
subplot(1,2,1)
plot(ampVOTact.left.deltas,bt(1,:),'o')
hold on
plot(ampVOTact.left.deltas,polyval(b,ampVOTact.left.deltas),'red')
hold off
xlabel('Deltas (uV)')
ylabel('Psychometric Curve Slope')
resid=bt(1,:)'-polyval(b,ampVOTact.left.deltas);
SSresid = sum(resid.^2);
SStotal = (length(bt(1,:))-1) * var(bt(1,:));
R2 = 1 - (SSresid/SStotal);% * (length(beta(1,:))-1)/(length(beta(1,:))-length(b)))
title(sprintf('VOT Active Left R2: %.4f', R2))


format long
b=polyfit(ampVOTact.right.deltas,bt(1,:)',1);
subplot(1,2,2)
plot(ampVOTact.right.deltas,bt(1,:),'o')
hold on
plot(ampVOTact.right.deltas,polyval(b,ampVOTact.right.deltas),'red')
hold off
xlabel('Deltas (uV)')
ylabel('Psychometric Curve Slope')
resid=bt(1,:)'-polyval(b,ampVOTact.right.deltas);
SSresid = sum(resid.^2);
SStotal = (length(bt(1,:))-1) * var(bt(1,:));
R2 = 1 - (SSresid/SStotal);% * (length(beta(1,:))-1)/(length(beta(1,:))-length(b)))
title(sprintf('VOT Active Right R2: %.4f', R2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
format long
b=polyfit(ampFormpass.left.deltas,bt(2,:)',1);
subplot(1,2,1)
plot(ampFormpass.left.deltas,bt(2,:),'o')
hold on
plot(ampFormpass.left.deltas,polyval(b,ampFormpass.left.deltas),'red')
hold off
xlabel('Deltas (uV)')
ylabel('Psychometric Curve Slope')
resid=bt(2,:)'-polyval(b,ampFormpass.left.deltas);
SSresid = sum(resid.^2);
SStotal = (length(bt(2,:))-1) * var(bt(2,:));
R2 = 1 - (SSresid/SStotal);% * (length(beta(1,:))-1)/(length(beta(1,:))-length(b)))
title(sprintf('Formants Passive Left R2: %.4f', R2))

format long
b=polyfit(ampFormpass.right.deltas,bt(2,:)',1);
subplot(1,2,2)
plot(ampFormpass.right.deltas,bt(2,:),'o')
hold on
plot(ampFormpass.right.deltas,polyval(b,ampFormpass.right.deltas),'red')
hold off
xlabel('Deltas (uV)')
ylabel('Psychometric Curve Slope')
resid=bt(2,:)'-polyval(b,ampFormpass.right.deltas);
SSresid = sum(resid.^2);
SStotal = (length(bt(2,:))-1) * var(bt(2,:));
R2 = 1 - (SSresid/SStotal);% * (length(beta(1,:))-1)/(length(beta(1,:))-length(b)))
title(sprintf('Formants Passive Right R2: %.4f', R2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4)
format long
b=polyfit(ampFormact.left.deltas,bt(2,:)',1);
subplot(1,2,1)
plot(ampFormact.left.deltas,bt(2,:),'o')
hold on
plot(ampFormact.left.deltas,polyval(b,ampFormact.left.deltas),'red')
hold off
xlabel('Deltas (uV)')
ylabel('Psychometric Curve Slope')
resid=bt(2,:)'-polyval(b,ampFormact.left.deltas);
SSresid = sum(resid.^2);
SStotal = (length(bt(2,:))-1) * var(bt(2,:));
R2 = 1 - (SSresid/SStotal);% * (length(beta(1,:))-1)/(length(beta(1,:))-length(b)))
title(sprintf('Formants Actie Left R2: %.4f', R2))

format long
b=polyfit(ampFormact.right.deltas,bt(2,:)',1);
subplot(1,2,2)
plot(ampFormact.right.deltas,bt(2,:),'o')
hold on
plot(ampFormact.right.deltas,polyval(b,ampFormact.right.deltas),'red')
hold off
xlabel('Deltas (uV)')
ylabel('Psychometric Curve Slope')
resid=bt(2,:)'-polyval(b,ampFormact.right.deltas);
SSresid = sum(resid.^2);
SStotal = (length(bt(2,:))-1) * var(bt(2,:));
R2 = 1 - (SSresid/SStotal);% * (length(beta(1,:))-1)/(length(beta(1,:))-length(b)))
title(sprintf('Formants Active Right R2: %.4f', R2))

close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)
boxplot([aux(1,:)',aux2(1,:)', aux3(1,:)',aux4(1,:)', aux5(1,:)', aux6(1,:)', aux7(1,:)', aux8(1,:)'],'Labels',{'V 1/5 p','V 3 p','V 1/5 a','V 3 a', 'F 1/5 p','F 3 p','F 1/5 a','F 3 a'})
ylabel('Amplitude N1-P2 (uV)')
title('Mean stimuli 1,2,4,5 vs 3, p=passive, a=active, V=VOT, F=Formants Left')

figure(6)
boxplot([aux(2,:)',aux2(2,:)', aux3(2,:)',aux4(2,:)', aux5(2,:)', aux6(2,:)', aux7(2,:)', aux8(2,:)'],'Labels',{'V 1/5 p','V 3 p','V 1/5 a','V 3 a', 'F 1/5 p','F 3 p','F 1/5 a','F 3 a'})
ylabel('Amplitude N1-P2 (uV)')
title('Mean stimuli 1,2,4,5 vs 3, p=passive, a=active, V=VOT, F=Formants Right')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Média por sujeito, lado, continuum e condição 
col={'r','b','k','green','c','yellow','green','k','magenta','c','b'};
col2 = [0 0 0; 0.2 0.2 0.2; 0.4 0.4 0.4; 0.6 0.6 0.6; 0.8 0.8 0.8];
Leg=cell(5,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nome='VOT/Ativo';
subj=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7)
set(gcf, 'Position',  [200, 200, 1300, 1000])
for i=1:5
    load(sprintf(strcat('../Sujeito%d/',nome,'/ChanStimICA%d'),subj,i))
    
    subplot(1,2,1)
    plot(-Xica_DWT_rec(17,:),col{i})
    title(sprintf('Left Subject %d',subj))
    grid
    hold on
    
    subplot(1,2,2)
    plot(-Xica_DWT_rec(2,:),col{i})
    title(sprintf('Right Subject %d',subj))
    Leg{i}=sprintf('%d',i);
    grid
    hold on
end
hold off
legend(Leg,'Location','northeast')

% Média das diferenças N1-P2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nome = ampFormact.left;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8)
set(gcf, 'Position',  [200, 200, 1000, 1000])
for j=1:nb_subj
    plot(1:5,nome.diffs(j,:),strcat(col{j},'.'),'markersize',40)
    text(1-0.1,nome.diffs(j,1),num2str(j))
    hold on
    plot(1:5,nome.diffs(j,:),strcat(col{j},'-'))
end
plot(1:5,mean(nome.diffs),'magenta*')
plot(1:5,mean(nome.diffs),'magenta-')
title('Diffs N1-P2')
axis([0.5 5.5 -inf inf])
hold off
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Média das diferenças P1-N1

figure(9)
set(gcf, 'Position',  [200, 200, 1000, 1000])
for i=1:5
    for j=1:nb_subj
    plot(i,nome.diffsp1(j,i),strcat(col{i},'o'))
    hold on
    end
end
plot(1:5,mean(nome.diffsp1),'magenta*')
plot(1:5,mean(nome.diffsp1),'magenta')
title('Diffs P1-N1')
axis([0.5 5.5 -inf inf])
hold off
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
enderecos={'VOT/Passivo','VOT/Ativo','Formantes/Passivo','Formantes/Ativo'};
auxDiffs(1:4,1:2,1:5,1:5)=0; %endereco, lado, linhas(picosP2, vales, diffN1P2, picosP1, diffP1N1), stims
%Média geral por lado, condição e continuum 
figure(10)
col2 = [0 0 0; 0 0 0; 0.4 0.4 0.4; 0.8 0.8 0.8; 0.8 0.8 0.8];
sty2={'-','--','.','--','-'};

for ed=1:4

XicaMat(1:2,1:nb_subj,1:4096)=0; %lado,sujeito,sinal
set(gcf, 'Position',  [100, 100, 1300, 800])
for i=1:5
    for j=1:nb_subj
    load(sprintf('../Sujeito%d/%s/ChanStimICA%d',j,enderecos{ed},i))
    XicaMat(1,j,:)=-Xica_DWT_rec(17,:);
    XicaMat(2,j,:)=-Xica_DWT_rec(2,:);
    end
    
    subplot(1,2,1)
    Z=squeeze(XicaMat(1,:,:));
    aux=mean(Z);
    if (ed==4)
    if i == 3
    h=plot(1:10:4096,aux(1:10:end));%col{i}
    set(h,'linewidth',2,'Color',col2(i,:), 'Marker',sty2{i},'LineStyle','none')   
    else
        h=plot(aux);%col{i}
        set(h,'linewidth',2,'Color',col2(i,:), 'LineStyle',sty2{i})
    end
    grid
    title('Left')
    xlabel('Time (s)')
    ylabel('EEG Amplitude (uV)')
    hold on
    end
    auxDiffs(ed,1,1,i)=find(aux(1,1100:2000)==max(aux(1,1100:2000)))+1099;
    auxDiffs(ed,1,2,i)=find(aux(1,mN1:1500)==min(aux(1,mN1:1500)))+(mN1-1);
    auxDiffs(ed,1,3,i)= aux(1,auxDiffs(ed,1,1,i))-aux(1,auxDiffs(ed,1,2,i));
    auxDiffs(ed,1,4,i)=(find(aux(1,550:1050)==max(aux(1,550:1050))))+549;
    auxDiffs(ed,1,5,i)= aux(1,auxDiffs(ed,1,4,i))-aux(1,auxDiffs(ed,1,2,i));
    auxDiffs(ed,1,1,i)=aux(find(aux(1,1100:2000)==max(aux(1,1100:2000)))+1099);
    auxDiffs(ed,1,2,i)=aux(find(aux(1,mN1:1500)==min(aux(1,mN1:1500)))+(mN1-1));
    auxDiffs(ed,1,4,i)=aux(find(aux(1,550:1050)==max(aux(1,550:1050)))+549);
    
    subplot(1,2,2)
    Z=squeeze(XicaMat(2,:,:));
    aux=mean(Z);
    if (ed==4)
    plot(aux,col{i})
    grid
    title('Right')
    xlabel('Time (s)')
    ylabel('EEG Amplitude (uV)')
    Leg{i}=sprintf('%d',i);
    hold on
    end
    auxDiffs(ed,2,1,i)=find(aux(1,1100:2000)==max(aux(1,1100:2000)))+1099;
    auxDiffs(ed,2,2,i)=find(aux(1,mN1:1500)==min(aux(1,mN1:1500)))+(mN1-1);
    auxDiffs(ed,2,3,i)= aux(1,auxDiffs(ed,2,1,i))-aux(1,auxDiffs(ed,2,2,i));
    auxDiffs(ed,2,4,i)=find(aux(1,550:1050)==max(aux(1,550:1050)))+549;
    auxDiffs(ed,2,5,i)= aux(1,auxDiffs(ed,2,4,i))-aux(1,auxDiffs(ed,2,2,i));
    auxDiffs(ed,2,1,i)=aux(find(aux(1,1100:2000)==max(aux(1,1100:2000)))+1099);
    auxDiffs(ed,2,2,i)=aux(find(aux(1,mN1:1500)==min(aux(1,mN1:1500)))+(mN1-1));
    auxDiffs(ed,2,4,i)=aux(find(aux(1,550:1050)==max(aux(1,550:1050)))+549);
end
end
hold off
legend(Leg,'Location','northeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nomes={'VOTpass','VOTact','FormPass','FormAct'};
figure(11)
set(gcf, 'Position',  [200, 200, 1000, 1000])
for ed=1:4
%Plota N1-P2 das curvas médias da figura anterior
aux10=squeeze(auxDiffs(ed,1,3,:));
aux20=squeeze(auxDiffs(ed,2,3,:));
for i=1:5
    subplot(1,2,1)
    plot(i,aux10(i),strcat(col{i},'.'),'markersize',30)
    title('Left')
    grid
    hold on
    
    subplot(1,2,2)
    plot(i,aux20(i),strcat(col{i},'.'),'markersize',30)
    title('Right')
    Leg{i}=sprintf('%d',i);
    grid
    hold on
end
subplot(1,2,1)

plot(1:5,aux10,'magenta')
text(5.05,aux10(5),nomes{ed})
title('Diffs N1-P2 all')

subplot(1,2,2)
plot(1:5,aux20,'magenta')
text(5.05,aux20(5),nomes{ed})
title('Diffs N1-P2 all')
axis([0.5 5.5 -inf inf])
end
hold off
legend(Leg,'Location','best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(12)
%Plota P1 médio de acordo com estímulo
set(gcf, 'Position',  [200, 200, 1000, 1000])
for ed=1:4
aux10=abs(squeeze(auxDiffs(ed,1,4,:)));
aux20=abs(squeeze(auxDiffs(ed,2,4,:)));
for i=1:5
    subplot(1,2,1)
    plot(i,aux10(i),strcat(col{i},'.'),'markersize',30)
    title('Left')
    grid
    hold on
    
    subplot(1,2,2)
    plot(i,aux20(i),strcat(col{i},'.'),'markersize',30)
    title('Right')
    Leg{i}=sprintf('%d',i);
    grid
    hold on
end
subplot(1,2,1)
plot(1:5,aux10,'magenta')
text(5.05,aux10(5),nomes{ed})
title('P1 all')
subplot(1,2,2)
plot(1:5,aux20,'magenta')
text(5.05,aux20(5),nomes{ed})
title('P1 all')
end
hold off
legend(Leg,'Location','best')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(13)
%Plota P1-N1 médio de acordo com estímulo
set(gcf, 'Position',  [200, 200, 1000, 1000])
for ed=1:4
aux10=abs(squeeze(auxDiffs(ed,1,5,:)));
aux20=abs(squeeze(auxDiffs(ed,2,5,:)));
for i=1:5
    subplot(1,2,1)
    plot(i,aux10(i),strcat(col{i},'.'),'markersize',30)
    title('Left')
    grid
    hold on
    
    subplot(1,2,2)
    plot(i,aux20(i),strcat(col{i},'.'),'markersize',30)
    title('Right')
    Leg{i}=sprintf('%d',i);
    grid
    hold on
end
subplot(1,2,1)
plot(1:5,aux10,'magenta')
text(5.05,aux10(5),nomes{ed})
title('P1-N1 all')
subplot(1,2,2)
plot(1:5,aux20,'magenta')
text(5.05,aux20(5),nomes{ed})
title('P1-N1 all')
end
hold off
legend(Leg,'Location','best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(14)
%Plota N1 médio de acordo com estímulo
set(gcf, 'Position',  [200, 200, 1000, 1000])
for ed=1:4
aux10=abs(squeeze(auxDiffs(ed,1,2,:)));
aux20=abs(squeeze(auxDiffs(ed,2,2,:)));
for i=1:5
    subplot(1,2,1)
    plot(i,aux10(i),strcat(col{i},'.'),'markersize',30)
    title('Left')
    grid
    hold on
    
    subplot(1,2,2)
    plot(i,aux20(i),strcat(col{i},'.'),'markersize',30)
    title('Right')
    Leg{i}=sprintf('%d',i);
    grid
    hold on
end
subplot(1,2,1)
plot(1:5,aux10,'magenta')
text(5.05,aux10(5),nomes{ed})
title('N1 all')
subplot(1,2,2)
plot(1:5,aux20,'magenta')
text(5.05,aux20(5),nomes{ed})
title('N1 all')
end
hold off
legend(Leg,'Location','best')       
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(15)
%Plota P2 médio de acordo com estímulo
set(gcf, 'Position',  [200, 200, 1000, 1000])
for ed=1:4
aux10=squeeze(auxDiffs(ed,1,1,:));
aux20=squeeze(auxDiffs(ed,2,1,:));
for i=1:5
    subplot(1,2,1)
    plot(i,aux10(i),strcat(col{i},'.'),'markersize',30)
    title('Left')
    grid
    hold on
    
    subplot(1,2,2)
    plot(i,aux20(i),strcat(col{i},'.'),'markersize',30)
    title('Right')
    Leg{i}=sprintf('%d',i);
    grid
    hold on
    
end
subplot(1,2,1)
plot(1:5,aux10,'magenta')
text(5.05,aux10(5),nomes{ed})
title('P2 all')
subplot(1,2,2)
plot(1:5,aux20,'magenta')
text(5.05,aux20(5),nomes{ed})
title('P2 all')
end
hold off
legend(Leg,'Location','best')       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Adaptação para código do Hani
nome=ampVOTact;

figure(16)

set(gcf, 'Position',  [200, 200, 1000, 1000])

YL=nome.left.picos_val;
YR=nome.right.picos_val;
X=1:5;
subplot(1,2,1)
plotmeanstd(X,YL)
title('P2 Left')
grid
subplot(1,2,2)
plotmeanstd(X,YR)
title('P2 Right')
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(17)

boxplot([box1(1,:)',box2(1,:)', box3(1,:)',box4(1,:)'],'Labels',{'V pass','V act','F pass','F act'})
ylabel('Amplitude N1-P2 (uV)')
title('Mean stimuli 1,2,3,4,5 V=VOT, F=Formants Left')

figure(18)

boxplot([box1(2,:)',box2(2,:)', box3(2,:)',box4(2,:)'],'Labels',{'V pass','V act','F pass','F act'})
ylabel('Amplitude N1-P2 (uV)')
title('Mean stimuli 1,2,3,4,5 V=VOT, F=Formants Right')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 box1(1,:)=mean(ampVOTpass.left.picos_val,2)';
 box1(2,:)=mean(ampVOTpass.right.picos_val,2)';
    
 box2(1,:)=mean(ampVOTact.left.picos_val,2)';
 box2(2,:)=mean(ampVOTact.right.picos_val,2)';
    
 box3(1,:)=mean(ampFormpass.left.picos_val,2)';
 box3(2,:)=mean(ampFormpass.right.picos_val,2)';
    
 box4(1,:)=mean(ampFormact.left.picos_val,2)';
 box4(2,:)=mean(ampFormact.right.picos_val,2)';


figure(19)

boxplot([box1(1,:)',box2(1,:)', box3(1,:)',box4(1,:)'],'Labels',{'V pass','V act','F pass','F act'})
ylabel('Amplitude P2 (uV)')
title('Mean stimuli 1,2,3,4,5 V=VOT, F=Formants Left')

figure(20)

boxplot([box1(2,:)',box2(2,:)', box3(2,:)',box4(2,:)'],'Labels',{'V pass','V act','F pass','F act'})
ylabel('Amplitude P2 (uV)')
title('Mean stimuli 1,2,3,4,5 V=VOT, F=Formants Right')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 box1(1,:)=mean(ampVOTpass.left.vales_val,2)';
 box1(2,:)=mean(ampVOTpass.right.vales_val,2)';
    
 box2(1,:)=mean(ampVOTact.left.vales_val,2)';
 box2(2,:)=mean(ampVOTact.right.vales_val,2)';
    
 box3(1,:)=mean(ampFormpass.left.vales_val,2)';
 box3(2,:)=mean(ampFormpass.right.vales_val,2)';
    
 box4(1,:)=mean(ampFormact.left.vales_val,2)';
 box4(2,:)=mean(ampFormact.right.vales_val,2)';

figure(21)

boxplot([box1(1,:)',box2(1,:)', box3(1,:)',box4(1,:)'],'Labels',{'V pass','V act','F pass','F act'})
ylabel('Amplitude N1 (uV)')
title('Mean stimuli 1,2,3,4,5 V=VOT, F=Formants Left')

figure(22)

boxplot([box1(2,:)',box2(2,:)', box3(2,:)',box4(2,:)'],'Labels',{'V pass','V act','F pass','F act'})
ylabel('Amplitude N1 (uV)')
title('Mean stimuli 1,2,3,4,5 V=VOT, F=Formants Right')