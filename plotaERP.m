%Plota potenciais evocados de todos
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

clear all
close all

s=[1,3,4,5,6,7,9,10,11];
nb_subj=size(s,2);
fs=5000;
subj=4; %se for plotar um sujeito específico, coloque sua posição no vetor s aqui

enderecos={'VOT/Passivo','VOT/Ativo','Formantes/Passivo','Formantes/Ativo'};
enderecos2={'VOT','VOT','Formantes','Formantes'};
nome={'Figuras/VOTpass.pdf','Figuras/VOTativ.pdf','Figuras/Formpass.pdf','Figuras/Formativ.pdf'};

%Informação para inicializar vetores
load(sprintf('../Sujeito%d/%s/ChanStimICA%d.mat',1,enderecos{1},1))
Xica_DWT_rec=Xica_DWT_rec(:,501:end);
L=size(Xica_DWT_rec,2);
t=(1:L)/fs;

col2 = [0 0 0; 0.2 0.2 0.2; 0.4 0.4 0.4; 0.6 0.6 0.6; 0.6 0.6 0.6];
sty2={'-','--','.','--','-'};
stim1={'da','da?','da-ta?','ta?','ta'};
stim2={'pa','pa?','pa-pe?','pe?','pe'};
Leg=cell(1,5);

figure(1)

for ad=1:4
    close all
XicaMat(1:2,1:nb_subj,1:L)=0; %lado,sujeito,sinal
set(gcf, 'Position',  [100, 100, 1200, 700])
for i=1:5
    count=0;
    for j=s
    count=count+1;
    [medL, medR, sdL, sdR]=normaliza(enderecos2{ad},j);
    load(sprintf('../Sujeito%d/%s/ChanStimICA%d.mat',j,enderecos{ad},i))
    Xica_DWT_rec=Xica_DWT_rec(:,501:end);
    XicaMat(1,count,:)=-((Xica_DWT_rec(17,:)-medR)/sdR);
    XicaMat(2,count,:)=-((Xica_DWT_rec(2,:)-medR)/sdR);
    end
    
    subplot(1,2,1)
    Z=squeeze(XicaMat(1,:,:));
    
    aux=mean(Z);
    %aux=Z(subj,:);
    if i == 3
        for k=1:size(aux,2)
            if (mod(k,15)==0)
                h=plot(t(1,k),aux(1,k));
                set(h,'linewidth',2,'Color',col2(i,:), 'Marker',sty2{i},'LineStyle','none')
            end
        end
    else
        h=plot(t,aux);
        set(h,'linewidth',2,'Color',col2(i,:), 'LineStyle',sty2{i})
    end
    grid
    title('Grand Average Left')
    xlabel('Time (s)')
    ylabel('EEG Amplitude (uV)')
    hold on
    axis([-inf inf -2.6 3.3])
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    subplot(1,2,2)
    Z=squeeze(XicaMat(2,:,:));
    aux=mean(Z);
    %aux=Z(subj,:);
    if i == 3
        aux2=aux;
        h=plot(t(1,1),aux(1,1));
        set(h,'linewidth',2,'Color',col2(i,:), 'Marker',sty2{i},'LineStyle','none')
    else
        h=plot(t,aux);
        set(h,'linewidth',2,'Color',col2(i,:), 'LineStyle',sty2{i})
    end
    grid
    title('Grand Average Right')
    xlabel('Time (s)')
    ylabel('EEG Amplitude (uV)')
    
    hold on
    axis([-inf inf -2.6 3.3])
    if ad<3
        Leg{1,i}=sprintf(stim1{i});
    else
        Leg{1,i}=sprintf(stim2{i});
    end
end

for k=1:size(aux2,2)
            if (mod(k,15)==0)
                h=plot(t(1,k),aux2(1,k));
                set(h,'linewidth',2,'Color',col2(3,:), 'Marker',sty2{3},'LineStyle','none')
            end
end
hold off
legend(Leg,'Location','northeast')

orient(figure(1),'landscape')
print(figure(1),nome{ad},'-dpdf')
end
