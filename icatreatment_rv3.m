%This script executes the ICA treatment at the mean signals and save the
%resulting matrix with channels after the treatment to be used in the
%averaring. This is performed twice, being once for passive and another 
%for active data. Saves are made in the corresponding folders. This version
%uses the function ica develped by professor Hani

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Here a matrix with independent components (IC's) is saved to be denoised
%with DWT in denoising.r before reconstruction for all subjects and stages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Passive stage data treatment
clear all
close all

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

local={'VOT/Passivo/','VOT/Ativo/','Formantes/Passivo/','Formantes/Ativo/'};

for (sujeito=1:11)
    for (ad=1:4)
        for (stim=1:5)
            address=strcat(sprintf('../Sujeito%d/',sujeito),local(ad));
            load(strcat(address{1},sprintf('ChanStim%d.mat',stim)));
            [N1,M] = size(chan_stim);
            N = N1-1;
            ch_ref = 12;
            X=chan_stim([1:(ch_ref-1) (ch_ref+1):N1],:)-repmat(chan_stim(ch_ref,:),N,1);
            X0 = X - repmat(mean(X,2),1,M);
            X0norm = X0./repmat(std(X')',1,M)*3;

            %Compute ICA
            [W,w0,Xica]=ica(X0norm,[],[],200,0.01);
            save(strcat(address{1},sprintf('Xica%d',stim)),'Xica','W','w0','X')
        end
    end
end

sprintf('Rode agora o script denoising.r e depois volte para continuar esse aqui')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HERE CALL THE R SCRIPT DENOISING.R BEFORE CONTINUE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for (sujeito=1:11)
    for (ad=1:4)
        for (stim=1:5)
            address=strcat(sprintf('../Sujeito%d/',sujeito),local(ad));
            load(strcat(address{1},sprintf('Xica_DWT%d.mat',stim)));
            
            %Reconstruction of DWT denoised signal
            Xica_DWT_rec = W\(Xica_DWT - repmat(w0,1,M));
            Xica_DWT_rec=Xica_DWT_rec.*repmat(std(X')',1,M)/3;
            
            %Baseline correction
            for (ch=1:size(Xica_DWT_rec,1))
                baseline=mean(Xica_DWT_rec(ch,1:(0.1*rate)));
                Xica_DWT_rec(ch,:)=Xica_DWT_rec(ch,:)-baseline;
            end

            %Save treated signals
            nome=strcat(address{1},sprintf('ChanStimICA%d',stim));
            save(nome,'Xica_DWT_rec')
        end
    end
end

%Eliminates bad signals
idx = zeros(1,N);
spike = [max(abs(Xica),[],2) std(Xica,[],2) max(abs(Xica),[],2)./std(Xica,[],2)];
idx = find(spike(:,3) > 20);
Xica_signal = Xica;
Xica_signal([idx],:) = 0;

%Reconstruction of only Xica signal
Xica_rec = W\(Xica_signal - repmat(w0,1,M));

Xica_rec=Xica_rec.*repmat(std(X')',1,M)/3;
%Reconstructed ICA signals baseline correction
for (ch=1:size(Xica_rec,1))
    baseline=mean(Xica_rec(ch,1:(0.1*rate)));
    Xica_rec(ch,:)=Xica_rec(ch,:)-baseline;
end

%Filtering Reconstructed signals
fc = 156.25; %Hz
[b,a] =  butter(6,fc/(rate/2));
aux1 = filter(b,a,Xica_rec,[],2);
aux2 = filter(b,a,fliplr(aux1),[],2);
Xica_rec_filt = fliplr(aux2);

%Original signals baseline correction
for (ch=1:size(X,1))
    baseline=mean(X(ch,1:(0.1*rate)));
    X(ch,:)=X(ch,:)-baseline;
end
%Filtering original signals
aux1 = filter(b,a,X,[],2);
aux2 = filter(b,a,fliplr(aux1),[],2);
X_filt=fliplr(aux2);

t_size=M/rate;
t=-0.1:1/rate:((t_size-0.1)-1/rate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot results of the last subject
chan={'T5','Tp10','T6','Fz','Oz','F7','Fp2','F4','C4','T4','F3','T3','C3','Fp1','Pz','F8','Tp9'};

%Original signals
figure(1)
for i=1:N
  subplot(3,6,i)
  plot(t,X(i,:))
  axis([-0.1 0.4 -inf inf]);
  xlabel('Time (s)')
  ylabel('Magnitude (uV)')
  title(sprintf('OrgChannel %s',chan{i}))
  grid 
end

%ICA signals baseline correction
for (ch=1:size(Xica,1))
    baseline=mean(Xica(ch,1:(0.1*rate)));
    Xica(ch,:)=Xica(ch,:)-baseline;
end
%ICA components 
figure(2)
for i=1:N
  subplot(3,6,i)
  plot(t,Xica(i,:))
  axis([-0.1 0.4 -inf inf]);
  xlabel('Time (s)')
  ylabel('Magnitude (uV)')
  title(sprintf('IC %s',chan{i}))
  grid 
end

%Reconstructed signals
snrIca(1:N)=0;
figure(3)
for i=1:N
  subplot(3,6,i)
  plot(t,X(i,:),'r-')
  hold on
  plot(t,Xica_rec(i,:),'b-','linewidth',2)
  axis([-0.1 0.4 min(min(Xica_rec)) max(max(Xica_rec))]);
  xlabel('Time (s)')
  ylabel('Magnitude (uV)')
  title(sprintf('Rec ICA %s',chan{i}))
  grid
  hold off
  v=snr(Xica_rec(i,:), rate);
  snrIca(i)=v;
  Leg=sprintf('SNR:%g',v);
  text(0,(max(max(Xica_rec))-0.5),Leg)
  if (i==N)
     text(0.7,(max(max(Xica_rec))-1), '--- Original signal','Color','r') 
     text(0.7,(max(max(Xica_rec))-2), '--- Treated signal','Color','b')
     text(0.7,(max(max(Xica_rec))-3),sprintf('MeanSNR: %g',mean(snrIca)))
  end
end

%Reconstructed signals filtered
snrXica(1:N)=0;
figure(4)
for i=1:N
  subplot(3,6,i)
  plot(t,X(i,:),'r-')
  hold on
  plot(t,Xica_rec_filt(i,:),'b-','linewidth',2)
  axis([-0.1 0.4 min(min(Xica_rec_filt)) max(max(Xica_rec_filt))]);
  xlabel('Time (s)')
  ylabel('Magnitude (uV)')
  title(sprintf('Rec ICA Filt. %s',chan{i}))
  grid
  hold off
  v=snr(Xica_rec_filt(i,:), rate);
  snrXica(i)=v;
  Leg=sprintf('SNR:%g',v);
  text(0,(max(max(Xica_rec_filt))-0.5),Leg)
  if (i==N)
     text(0.7,(max(max(Xica_rec_filt))-0.5), '--- Original signal','Color','r') 
     text(0.7,(max(max(Xica_rec_filt))-1.5), '--- Treated signal','Color','b')
     text(0.7,(max(max(Xica_rec_filt))-3),sprintf('MeanSNR: %g',mean(snrXica)))
  end
end

%Reconstructed DWT not filtered signals
snrXDWT(1:N)=0;
figure(5)
for i=1:N
  subplot(3,6,i)
  plot(t,X(i,:),'r-')
  hold on
  plot(t,Xica_DWT_rec(i,:),'b-','linewidth',2)
  axis([-0.1 0.4 min(min(Xica_DWT_rec)) max(max(Xica_DWT_rec))]);
  xlabel('Time (s)')
  ylabel('Magnitude (uV)')
  title(sprintf('Rec ICA DWT %s',chan{i}))
  grid
  hold off
  v=snr(Xica_DWT_rec(i,:), rate);
  snrXDWT(i)=v;
  Leg=sprintf('SNR:%g',v);
  text(0,(max(max(Xica_DWT_rec))-0.5),Leg)
  if (i==N)
     text(0.7,(max(max(Xica_DWT_rec))-1), '--- Original signal','Color','r') 
     text(0.7,(max(max(Xica_DWT_rec))-2), '--- Treated signal','Color','b')
     text(0.7,(max(max(Xica_DWT_rec))-3),sprintf('MeanSNR: %g',mean(snrXDWT)))
  end
end

%Signals just filtered
snrXfilt(1:N)=0;
figure(6)
for i=1:N
  subplot(3,6,i)
  plot(t,X(i,:),'r-')
  hold on
  plot(t,X_filt(i,:),'b-','linewidth',2)
  axis([-0.1 0.4 min(min(X_filt)) max(max(X_filt))]);
  xlabel('Time (s)')
  ylabel('Magnitude (uV)')
  title(sprintf('Orig Filt. %s',chan{i}))
  grid
  hold off
  v=snr(X_filt(i,:), rate);
  snrXfilt(i)=v;
  Leg=sprintf('SNR:%g',v);
  text(0,(max(max(X_filt))-0.5),Leg)
  if (i==N)
     text(0.7,(max(max(X_filt))-1), '--- Original signal','Color','r') 
     text(0.7,(max(max(X_filt))-2), '--- Treated signal','Color','b')
     text(0.7,(max(max(X_filt))-3),sprintf('MeanSNR: %g',mean(snrXfilt)))
  end
end

%Difference between reconstructed DWT signals and just filtered signals
figure(7)
for i=1:N
  subplot(3,6,i)
  plot(t,Xica_DWT_rec(i,:)-X_filt(i,:))
  axis([-0.1 0.4 min(min(X_filt)) max(max(X_filt))]);
  xlabel('Time (s)')
  ylabel('Magnitude (uV)')
  title(sprintf('Diff. DWT-X filt %s',chan{i}))
  grid 
end

%Difference between reconstructed DWT signals and just Xica_filt
figure(8)
for i=1:N
  subplot(3,6,i)
  plot(t,Xica_DWT_rec(i,:)-Xica_rec_filt(i,:))
  axis([-0.1 0.4 min(min(Xica_rec_filt)) max(max(Xica_rec_filt))]);
  xlabel('Time (s)')
  ylabel('Magnitude (uV)')
  title(sprintf('Diff. DWT-Xica filt. %s',chan{i}))
  grid 
end

%Difference between reconstructed DWT signals and original signals
figure(9)
for i=1:N
  subplot(3,6,i)
  plot(t,Xica_DWT_rec(i,:)-X(i,:))
  axis([-0.1 0.4 min(min(X)) max(max(X))]);
  xlabel('Time (s)')
  ylabel('Magnitude (uV)')
  title(sprintf('Diff. DWT-X %s',chan{i}))
  grid 
end

%Difference between reconstructed ICA signals filtered and original signals
figure(10)
for i=1:N
  subplot(3,6,i)
  plot(t,Xica_rec_filt(i,:)-X(i,:))
  axis([-0.1 0.4 min(min(X)) max(max(X))]);
  xlabel('Time (s)')
  ylabel('Magnitude (uV)')
  title(sprintf('Diff. XicaFilt-X %s',chan{i}))
  grid 
end

%Difference between signals filtered and original signals
figure(11)
for i=1:N
  subplot(3,6,i)
  plot(t,X_filt(i,:)-X(i,:))
  axis([-0.1 0.4 min(min(X)) max(max(X))]);
  xlabel('Time (s)')
  ylabel('Magnitude (uV)')
  title(sprintf('Diff. X Filt-X %s',chan{i}))
  grid 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
