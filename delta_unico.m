%Repruce regression between delta_ERP and slope of psychometric curve made
%by Bidelman 2017 using the mean of the frontocentral electrodes
s=[1,2,3,4,5,6,7,8,9,10];
nb_subj=length(s);
med=[1,2,4,5];

bt(1:2,1:nb_subj)=0;

mN1=900;

picos=zeros(nb_subj,5);
vales=zeros(nb_subj,5);
diffs=zeros(nb_subj,5);
deltas=zeros(nb_subj,1);

ampFormpass=struct('picos',picos,'vales',vales,'diffs',diffs,'deltas',deltas);
ampFormact=struct('picos',picos,'vales',vales,'diffs',diffs,'deltas',deltas);
ampVOTact=struct('picos',picos,'vales',vales,'diffs',diffs,'deltas',deltas);
ampVOTpass=struct('picos',picos,'vales',vales,'diffs',diffs,'deltas',deltas);

for subj=1:nb_subj
        load(sprintf('../Sujeito%d/VOT/Identification/beta.mat',subj))
        bt(1,subj)=beta;
        
        load(sprintf('../Sujeito%d/Formantes/Identification/beta.mat',subj))
        bt(2,subj)=beta;
        
        for stim=1:5
        
        load(sprintf('../Sujeito%d/VOT/Passivo/ChanStimICA%d',subj,stim))
        sig=-mean(Xica_DWT_rec([2,4,6,7,8,11,14,16,17],:));
        ampVOTpass.vales(subj,stim)=(find(sig(1,mN1:1500)==min(sig(1,mN1:1500))))+(mN1-1);
        ampVOTpass.picos(subj,stim)=(find(sig(1,1000:1800)==max(sig(1,1000:1800))))+999;
        ampVOTpass.diffs(subj,stim)= sig(1,ampVOTpass.picos(subj,stim))-(sig(1,ampVOTpass.vales(subj,stim)));
        
        load(sprintf('../Sujeito%d/VOT/Ativo/ChanStimICA%d',subj,stim))
        sig=-mean(Xica_DWT_rec([2,4,6,7,8,11,14,16,17],:));
        ampVOTact.vales(subj,stim)=(find(sig(1,mN1:1500)==min(sig(1,mN1:1500))))+(mN1-1);
        ampVOTact.picos(subj,stim)=(find(sig(1,1000:1800)==max(sig(1,1000:1800))))+999;
        ampVOTact.diffs(subj,stim)= sig(1,ampVOTact.picos(subj,stim))-(sig(1,ampVOTact.vales(subj,stim)));
        
        load(sprintf('../Sujeito%d/Formantes/Passivo/ChanStimICA%d',subj,stim))
        sig=-mean(Xica_DWT_rec([2,4,6,7,8,11,14,16,17],:));
        ampFormpass.vales(subj,stim)=(find(sig(1,mN1:1500)==min(sig(1,mN1:1500))))+(mN1-1);
        ampFormpass.picos(subj,stim)=(find(sig(1,1000:1800)==max(sig(1,1000:1800))))+999;
        ampFormpass.diffs(subj,stim)= sig(1,ampFormpass.picos(subj,stim))-(sig(1,ampFormpass.vales(subj,stim)));
        
        load(sprintf('../Sujeito%d/Formantes/Ativo/ChanStimICA%d',subj,stim))
        sig=-mean(Xica_DWT_rec([2,4,6,7,8,11,14,16,17],:)); 
        ampFormact.vales(subj,stim)=(find(sig(1,mN1:1500)==min(sig(1,mN1:1500))))+(mN1-1);
        ampFormact.picos(subj,stim)=(find(sig(1,1000:1800)==max(sig(1,1000:1800))))+999;
        ampFormact.diffs(subj,stim)= sig(1,ampFormact.picos(subj,stim))-(sig(1,ampFormact.vales(subj,stim)));
        
        end
        
        ampVOTpass.deltas(subj)=mean(ampVOTpass.diffs(subj,med))-ampVOTpass.diffs(subj,3);
        
        ampVOTact.deltas(subj)=mean(ampVOTact.diffs(subj,med))-ampVOTact.diffs(subj,3);
        
        ampFormact.deltas(subj)=mean(ampFormact.diffs(subj,med))-ampFormact.diffs(subj,3);
        
        ampFormpass.deltas(subj)=mean(ampFormpass.diffs(subj,med))-ampFormpass.diffs(subj,3);
        
end

ampFormact.deltas=ampFormact.deltas/rssq(ampFormact.deltas);
ampFormpass.deltas=ampFormpass.deltas/rssq(ampFormpass.deltas);
ampVOTact.deltas=ampVOTact.deltas/rssq(ampVOTact.deltas);
ampVOTpass.deltas=ampVOTpass.deltas/rssq(ampVOTpass.deltas);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1)
b=polyfit(ampVOTpass.deltas,bt(1,:)',1);
plot(ampVOTpass.deltas,bt(1,:),'o')
hold on
plot(ampVOTpass.deltas,polyval(b,ampVOTpass.deltas),'red')
hold off
xlabel('Deltas (uV)')
ylabel('Psychometric Curve Slope')
resid=bt(1,:)'-polyval(b,ampVOTpass.deltas);
SSresid = sum(resid.^2);
SStotal = (length(bt(1,:))-1) * var(bt(1,:));
R2 = 1 - (SSresid/SStotal);% * (length(beta(1,:))-1)/(length(beta(1,:))-length(b)))
title(sprintf('VOT Passive R2: %.4f', R2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,2)
b=polyfit(ampVOTact.deltas,bt(1,:)',1);
plot(ampVOTact.deltas,bt(1,:),'o')
hold on
plot(ampVOTact.deltas,polyval(b,ampVOTact.deltas),'red')
hold off
xlabel('Deltas (uV)')
ylabel('Psychometric Curve Slope')
resid=bt(1,:)'-polyval(b,ampVOTact.deltas);
SSresid = sum(resid.^2);
SStotal = (length(bt(1,:))-1) * var(bt(1,:));
R2 = 1 - (SSresid/SStotal);% * (length(beta(1,:))-1)/(length(beta(1,:))-length(b)))
title(sprintf('VOT Active R2: %.4f', R2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,3)
format long
b=polyfit(ampFormpass.deltas,bt(2,:)',1);
plot(ampFormpass.deltas,bt(2,:),'o')
hold on
plot(ampFormpass.deltas,polyval(b,ampFormpass.deltas),'red')
hold off
xlabel('Deltas (uV)')
ylabel('Psychometric Curve Slope')
resid=bt(1,:)'-polyval(b,ampFormpass.deltas);
SSresid = sum(resid.^2);
SStotal = (length(bt(2,:))-1) * var(bt(2,:));
R2 = 1 - (SSresid/SStotal);% * (length(beta(1,:))-1)/(length(beta(1,:))-length(b)))
title(sprintf('Formants Passive R2: %.4f', R2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,4)
format long
b=polyfit(ampFormact.deltas,bt(2,:)',1);
plot(ampFormact.deltas,bt(2,:),'o')
hold on
plot(ampFormact.deltas,polyval(b,ampFormact.deltas),'red')
hold off
xlabel('Deltas (uV)')
ylabel('Psychometric Curve Slope')
resid=bt(1,:)'-polyval(b,ampFormact.deltas);
SSresid = sum(resid.^2);
SStotal = (length(bt(2,:))-1) * var(bt(2,:));
R2 = 1 - (SSresid/SStotal);% * (length(beta(1,:))-1)/(length(beta(1,:))-length(b)))
title(sprintf('Formants Active R2: %.4f', R2))

