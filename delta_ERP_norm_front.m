%Reproduce regression between delta_ERP and slope of psychometric curve made
%by Bidelman 2017
close all
clear all

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

%Subjects
s=[1,3,5,6,7,8,9,10,11];
sel=[1,2,3,4,5,6,7,8,9,10,11];
nb_subj=length(s);

bt(1:2,1:nb_subj)=0;

med=[1,2,4,5];

mN1=353;

fs=5000;

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
picos_lat=zeros(nb_subj,5);
vales=zeros(nb_subj,5);
vales_val=zeros(nb_subj,5);
vales_lat=zeros(nb_subj,5);
p1_val=zeros(nb_subj,5);
N2=zeros(nb_subj,5);
N2_val=zeros(nb_subj,5);
diffsp1=zeros(nb_subj,5);
diffs=zeros(nb_subj,5);
deltas=zeros(nb_subj,1);
deltasp2=zeros(nb_subj,1);
meio_val=zeros(nb_subj,5);
meio_lat=zeros(nb_subj,5);
left=struct('meio_lat',meio_lat,'meio_val',meio_val,'N2',N2,'N2_val',N2_val,'p1',p1,'picos',picos,'vales',vales,'picos_val',picos_val,'vales_val',vales_val,'diffsp1',diffsp1,'diffs',diffs,'deltas',deltas,'deltasp2',deltasp2,'picos_lat',picos_lat,'vales_lat',vales_lat);
right=struct('meio_lat',meio_lat,'meio_val',meio_val,'N2',N2,'N2_val',N2_val,'p1',p1,'picos',picos,'p1_val',p1_val,'vales',vales,'picos_val',picos_val,'vales_val',vales_val,'diffsp1',diffsp1,'diffs',diffs,'deltas',deltas,'deltasp2',deltasp2,'picos_lat',picos_lat,'vales_lat',vales_lat);

ampFormpass=struct('left',left,'right',right);
ampFormact=struct('left',left,'right',right);
ampVOTact=struct('left',left,'right',right);
ampVOTpass=struct('left',left,'right',right);

count=0;
for subj=s
count=count+1;
        load(sprintf('../../Sujeito%d/VOT/Identification/beta.mat',subj))
        bt(1,count)=beta;
        
        load(sprintf('../../Sujeito%d/Formantes/Identification/beta.mat',subj))
        bt(2,count)=beta;

        [medL, medR, sdL, sdR]=normaliza('VOT',subj);
        
        for stim=1:5
        
        load(sprintf('../../Sujeito%d/VOT/PassivoFront/ChanStimICA%d',subj,stim))
        Xica_DWT_rec=Xica_DWT_rec(:,501:end);
        Left=mean([Xica_DWT_rec(14,:);Xica_DWT_rec(11,:);Xica_DWT_rec(5,:)]);
        Right=mean([Xica_DWT_rec(16,:);Xica_DWT_rec(8,:);Xica_DWT_rec(7,:)]);
        Right=(Right-medR)/sdR;
        Left=(Left-medL)/sdL;
        ampVOTpass.right.vales(count,stim)=(find(-Right(mN1:1000)==min(-Right(mN1:1000))))+(mN1-1);
        ampVOTpass.right.p1(count,stim)=(find(-Right(50:550)==max(-Right(50:550))))+49;
        ampVOTpass.right.N2(count,stim)=(find(-Right(1300:1900)==min(-Right(1300:1900))))+1299;
        ampVOTpass.right.N2_val(count,stim)=-Right(ampVOTpass.right.N2(count,stim));
        ampVOTpass.right.picos(count,stim)=(find(-Right(500:1500)==max(-Right(500:1500))))+499;
        ampVOTpass.right.diffsp1(count,stim)= -Right(ampVOTpass.right.p1(count,stim))-(-Right(ampVOTpass.right.vales(count,stim)));
        ampVOTpass.right.diffs(count,stim)= -Right(ampVOTpass.right.picos(count,stim))-(-Right(ampVOTpass.right.vales(count,stim)));
        ampVOTpass.right.picos_val(count,stim)=-Right(ampVOTpass.right.picos(count,stim));
        ampVOTpass.right.vales_val(count,stim)=-Right(ampVOTpass.right.vales(count,stim));
        ampVOTpass.right.p1_val(count,stim)=-Right(ampVOTpass.right.p1(count,stim));
        ampVOTpass.right.picos_lat(count,stim)=ampVOTpass.right.picos(count,stim)/fs;
        ampVOTpass.right.vales_lat(count,stim)=ampVOTpass.right.vales(count,stim)/fs;
        ampVOTpass.right.meio_lat(count,stim)=(ampVOTpass.right.vales(count,stim)+round((ampVOTpass.right.picos(count,stim)-ampVOTpass.right.vales(count,stim))/2))/fs;
        ampVOTpass.right.meio_val(count,stim)=-Right((ampVOTpass.right.vales(count,stim)+round((ampVOTpass.right.picos(count,stim)-ampVOTpass.right.vales(count,stim))/2)));
        
        ampVOTpass.left.vales(count,stim)=(find(-Left(mN1:1000)==min(-Left(mN1:1000))))+(mN1-1);
        ampVOTpass.left.p1(count,stim)=(find(-Left(50:550)==max(-Left(50:550))))+49;
        ampVOTpass.left.N2(count,stim)=(find(-Left(1300:1900)==min(-Left(1300:1900))))+1299;
        ampVOTpass.left.N2_val(count,stim)=-Left(ampVOTpass.right.N2(count,stim));
        ampVOTpass.left.picos(count,stim)=(find(-Left(500:1500)==max(-Left(500:1500))))+499;
        ampVOTpass.left.diffsp1(count,stim)= -Left(ampVOTpass.left.p1(count,stim))-(-Left(ampVOTpass.left.vales(count,stim)));
        ampVOTpass.left.diffs(count,stim)= -Left(ampVOTpass.left.picos(count,stim))-(-Left(ampVOTpass.left.vales(count,stim)));
        ampVOTpass.left.picos_val(count,stim)=-Left(ampVOTpass.left.picos(count,stim));
        ampVOTpass.left.vales_val(count,stim)=-Left(ampVOTpass.left.vales(count,stim));
        ampVOTpass.left.p1_val(count,stim)=-Left(ampVOTpass.left.p1(count,stim));
        ampVOTpass.left.picos_lat(count,stim)=ampVOTpass.left.picos(count,stim)/fs;
        ampVOTpass.left.vales_lat(count,stim)=ampVOTpass.left.vales(count,stim)/fs;
        ampVOTpass.left.meio_lat(count,stim)=(ampVOTpass.left.vales(count,stim)+round((ampVOTpass.left.picos(count,stim)-ampVOTpass.left.vales(count,stim))/2))/fs;
        ampVOTpass.left.meio_val(count,stim)=-Left((ampVOTpass.left.vales(count,stim)+round((ampVOTpass.left.picos(count,stim)-ampVOTpass.left.vales(count,stim))/2)));
        
        end
        
        for stim=1:5
        
        load(sprintf('../../Sujeito%d/VOT/AtivoFront/ChanStimICA%d',subj,stim))
        Xica_DWT_rec=Xica_DWT_rec(:,501:end);
        Left=mean([Xica_DWT_rec(14,:);Xica_DWT_rec(11,:);Xica_DWT_rec(5,:)]);
        Right=mean([Xica_DWT_rec(16,:);Xica_DWT_rec(8,:);Xica_DWT_rec(7,:)]);
        Right=(Right-medR)/sdR;
        Left(:)=(Left(:)-medL)/sdL;
        ampVOTact.right.vales(count,stim)=(find(-Right(mN1:1000)==min(-Right(mN1:1000))))+(mN1-1);
        ampVOTact.right.p1(count,stim)=(find(-Right(50:550)==max(-Right(50:550))))+49;
        ampVOTact.right.N2(count,stim)=(find(-Right(1300:1900)==min(-Right(1300:1900))))+1299;
        ampVOTact.right.N2_val(count,stim)=-Right(ampVOTact.right.N2(count,stim));
        ampVOTact.right.picos(count,stim)=(find(-Right(500:1500)==max(-Right(500:1500))))+499;
        ampVOTact.right.diffsp1(count,stim)= -Right(ampVOTact.right.p1(count,stim))-(-Right(ampVOTact.right.vales(count,stim)));
        ampVOTact.right.diffs(count,stim)= -Right(ampVOTact.right.picos(count,stim))-(-Right(ampVOTact.right.vales(count,stim)));
        ampVOTact.right.picos_val(count,stim)=-Right(ampVOTact.right.picos(count,stim));
        ampVOTact.right.vales_val(count,stim)=-Right(ampVOTact.right.vales(count,stim));
        ampVOTact.right.p1_val(count,stim)=-Right(ampVOTact.right.p1(count,stim));
        ampVOTact.right.picos_lat(count,stim)=ampVOTact.right.picos(count,stim)/fs;
        ampVOTact.right.vales_lat(count,stim)=ampVOTact.right.vales(count,stim)/fs;
        ampVOTact.right.meio_lat(count,stim)=(ampVOTact.right.vales(count,stim)+round((ampVOTact.right.picos(count,stim)-ampVOTact.right.vales(count,stim))/2))/fs;
        ampVOTact.right.meio_val(count,stim)=-Right((ampVOTact.right.vales(count,stim)+round((ampVOTact.right.picos(count,stim)-ampVOTact.right.vales(count,stim))/2)));
                
        ampVOTact.left.vales(count,stim)=(find(-Left(mN1:1000)==min(-Left(mN1:1000))))+(mN1-1);
        ampVOTact.left.p1(count,stim)=(find(-Left(50:550)==max(-Left(50:550))))+49;
        ampVOTact.left.N2(count,stim)=(find(-Left(1300:1900)==min(-Left(1300:1900))))+1299;
        ampVOTact.left.N2_val(count,stim)=-Left(ampVOTact.right.N2(count,stim));
        ampVOTact.left.picos(count,stim)=(find(-Left(500:1500)==max(-Left(500:1500))))+499;
        ampVOTact.left.diffsp1(count,stim)= -Left(ampVOTact.left.p1(count,stim))-(-Left(ampVOTact.left.vales(count,stim)));
        ampVOTact.left.diffs(count,stim)= -Left(ampVOTact.left.picos(count,stim))-(-Left(ampVOTact.left.vales(count,stim)));
        ampVOTact.left.picos_val(count,stim)=-Left(ampVOTact.left.picos(count,stim));
        ampVOTact.left.vales_val(count,stim)=-Left(ampVOTact.left.vales(count,stim));
        ampVOTact.left.p1_val(count,stim)=-Left(ampVOTact.left.p1(count,stim));
        ampVOTact.left.picos_lat(count,stim)=ampVOTact.left.picos(count,stim)/fs;
        ampVOTact.left.vales_lat(count,stim)=ampVOTact.left.vales(count,stim)/fs;
        ampVOTact.left.meio_lat(count,stim)=(ampVOTact.left.vales(count,stim)+round((ampVOTact.left.picos(count,stim)-ampVOTact.left.vales(count,stim))/2))/fs;
        ampVOTact.left.meio_val(count,stim)=-Left((ampVOTact.left.vales(count,stim)+round((ampVOTact.left.picos(count,stim)-ampVOTact.left.vales(count,stim))/2)));
        
        end
        
        [medL, medR, sdL, sdR]=normaliza('Formantes',subj);
        
        for stim=1:5
        
        load(sprintf('../../Sujeito%d/Formantes/PassivoFront/ChanStimICA%d',subj,stim))
        Xica_DWT_rec=Xica_DWT_rec(:,501:end);
        Left=mean([Xica_DWT_rec(14,:);Xica_DWT_rec(11,:);Xica_DWT_rec(5,:)]);
        Right=mean([Xica_DWT_rec(16,:);Xica_DWT_rec(8,:);Xica_DWT_rec(7,:)]);
        Right=(Right-medR)/sdR;
        Left(:)=(Left(:)-medL)/sdL;
        ampFormpass.right.vales(count,stim)=(find(-Right(353:1000)==min(-Right(353:1000))))+352;
        ampFormpass.right.p1(count,stim)=(find(-Right(50:550)==max(-Right(50:550))))+49;
        ampFormpass.right.N2(count,stim)=(find(-Right(1300:1900)==min(-Right(1300:1900))))+1299;
        ampFormpass.right.N2_val(count,stim)=-Right(ampFormpass.right.N2(count,stim));
        ampFormpass.right.picos(count,stim)=(find(-Right(500:1500)==max(-Right(500:1500))))+499;
        ampFormpass.right.diffsp1(count,stim)= -Right(ampFormpass.right.p1(count,stim))-(-Right(ampFormpass.right.vales(count,stim)));
        ampFormpass.right.diffs(count,stim)= -Right(ampFormpass.right.picos(count,stim))-(-Right(ampFormpass.right.vales(count,stim)));
        ampFormpass.right.picos_val(count,stim)=-Right(ampFormpass.right.picos(count,stim));
        ampFormpass.right.vales_val(count,stim)=-Right(ampFormpass.right.vales(count,stim));
        ampFormpass.right.p1_val(count,stim)=-Right(ampFormpass.right.p1(count,stim));
        ampFormpass.right.picos_lat(count,stim)=ampFormpass.right.picos(count,stim)/fs;
        ampFormpass.right.vales_lat(count,stim)=ampFormpass.right.vales(count,stim)/fs;
        ampFormpass.right.meio_lat(count,stim)=(ampFormpass.right.vales(count,stim)+round((ampFormpass.right.picos(count,stim)-ampFormpass.right.vales(count,stim))/2))/fs;
        ampFormpass.right.meio_val(count,stim)=-Right((ampFormpass.right.vales(count,stim)+round((ampFormpass.right.picos(count,stim)-ampFormpass.right.vales(count,stim))/2)));
                
        ampFormpass.left.vales(count,stim)=(find(-Left(253:1000)==min(-Left(253:1000))))+352;
        ampFormpass.left.p1(count,stim)=(find(-Left(50:550)==max(-Left(50:550))))+49;
        ampFormpass.left.N2(count,stim)=(find(-Left(1300:1900)==min(-Left(1300:1900))))+1299;
        ampFormpass.left.N2_val(count,stim)=-Left(ampFormpass.right.N2(count,stim));
        ampFormpass.left.picos(count,stim)=(find(-Left(500:1500)==max(-Left(500:1500))))+499;
        ampFormpass.left.diffsp1(count,stim)= -Left(ampFormpass.left.p1(count,stim))-(-Left(ampFormpass.left.vales(count,stim)));
        ampFormpass.left.diffs(count,stim)= -Left(ampFormpass.left.picos(count,stim))-(-Left(ampFormpass.left.vales(count,stim)));
        ampFormpass.left.picos_val(count,stim)=-Left(ampFormpass.left.picos(count,stim));
        ampFormpass.left.vales_val(count,stim)=-Left(ampFormpass.left.vales(count,stim));
        ampFormpass.left.p1_val(count,stim)=-Left(ampFormpass.left.p1(count,stim));
        ampFormpass.left.picos_lat(count,stim)=ampFormpass.left.picos(count,stim)/fs;
        ampFormpass.left.vales_lat(count,stim)=ampFormpass.left.vales(count,stim)/fs;
        ampFormpass.left.meio_lat(count,stim)=(ampFormpass.left.vales(count,stim)+round((ampFormpass.left.picos(count,stim)-ampFormpass.left.vales(count,stim))/2))/fs;
        ampFormpass.left.meio_val(count,stim)=-Left((ampFormpass.left.vales(count,stim)+round((ampFormpass.left.picos(count,stim)-ampFormpass.left.vales(count,stim))/2)));
        
        end

        for stim=1:5
        
        load(sprintf('../../Sujeito%d/Formantes/AtivoFront/ChanStimICA%d',subj,stim))
        Xica_DWT_rec=Xica_DWT_rec(:,501:end);
        Left=mean([Xica_DWT_rec(14,:);Xica_DWT_rec(11,:);Xica_DWT_rec(5,:)]);
        Right=mean([Xica_DWT_rec(16,:);Xica_DWT_rec(8,:);Xica_DWT_rec(7,:)]);
        Right=(Right-medR)/sdR;
        Left(:)=(Left(:)-medL)/sdL;
        ampFormact.right.vales(count,stim)=(find(-Right(mN1:1000)==min(-Right(mN1:1000))))+(mN1-1);
        ampFormact.right.p1(count,stim)=(find(-Right(50:550)==max(-Right(50:550))))+49;
        ampFormact.right.N2(count,stim)=(find(-Right(1300:1900)==min(-Right(1300:1900))))+1299;
        ampFormact.right.N2_val(count,stim)=-Right(ampFormact.right.N2(count,stim));
        ampFormact.right.picos(count,stim)=(find(-Right(500:1500)==max(-Right(500:1500))))+499;
        ampFormact.right.diffsp1(count,stim)= -Right(ampFormact.right.p1(count,stim))-(-Right(ampFormact.right.vales(count,stim)));
        ampFormact.right.diffs(count,stim)= -Right(ampFormact.right.picos(count,stim))-(-Right(ampFormact.right.vales(count,stim)));
        ampFormact.right.picos_val(count,stim)=-Right(ampFormact.right.picos(count,stim));
        ampFormact.right.vales_val(count,stim)=-Right(ampFormact.right.vales(count,stim));
        ampFormact.right.p1_val(count,stim)=-Right(ampFormact.right.p1(count,stim));
        ampFormact.right.picos_lat(count,stim)=ampFormact.right.picos(count,stim)/fs;
        ampFormact.right.vales_lat(count,stim)=ampFormact.right.vales(count,stim)/fs;
        ampFormact.right.meio_lat(count,stim)=(ampFormact.right.vales(count,stim)+round((ampFormact.right.picos(count,stim)-ampFormact.right.vales(count,stim))/2))/fs;
        ampFormact.right.meio_val(count,stim)=-Right((ampFormact.right.vales(count,stim)+round((ampFormact.right.picos(count,stim)-ampFormact.right.vales(count,stim))/2)));
                
        ampFormact.left.vales(count,stim)=(find(-Left(mN1:1000)==min(-Left(mN1:1000))))+(mN1-1);
        ampFormact.left.p1(count,stim)=(find(-Left(50:550)==max(-Left(50:550))))+49;
        ampFormact.left.N2(count,stim)=(find(-Left(1300:1900)==min(-Left(1300:1900))))+1299;
        ampFormact.left.N2_val(count,stim)=-Left(ampFormact.right.N2(count,stim));
        ampFormact.left.picos(count,stim)=(find(-Left(500:1500)==max(-Left(500:1500))))+499;
        ampFormact.left.diffsp1(count,stim)= -Left(ampFormact.left.p1(count,stim))-(-Left(ampFormact.left.vales(count,stim)));
        ampFormact.left.diffs(count,stim)= -Left(ampFormact.left.picos(count,stim))-(-Left(ampFormact.left.vales(count,stim)));
        ampFormact.left.picos_val(count,stim)=-Left(ampFormact.left.picos(count,stim));
        ampFormact.left.vales_val(count,stim)=-Left(ampFormact.left.vales(count,stim));
        ampFormact.left.p1_val(count,stim)=-Left(ampFormact.left.p1(count,stim));
        ampFormact.left.picos_lat(count,stim)=ampFormact.left.picos(count,stim)/fs;
        ampFormact.left.vales_lat(count,stim)=ampFormact.left.vales(count,stim)/fs;
        ampFormact.left.meio_lat(count,stim)=(ampFormact.left.vales(count,stim)+round((ampFormact.left.picos(count,stim)-ampFormact.left.vales(count,stim))/2))/fs;
        ampFormact.left.meio_val(count,stim)=-Left((ampFormact.left.vales(count,stim)+round((ampFormact.left.picos(count,stim)-ampFormact.left.vales(count,stim))/2)));
            
        
        end
        
        %Compute deltas of N1-P2 mean(stim1~5)-stim3
        ampVOTpass.right.deltas(count)=mean(ampVOTpass.right.diffs(count,med))-ampVOTpass.right.diffs(count,3);
        ampVOTpass.left.deltas(count)=mean(ampVOTpass.left.diffs(count,med))-ampVOTpass.left.diffs(count,3);
        ampVOTpass.right.deltasp2(count)=mean(ampVOTpass.right.picos_val(count,med))-ampVOTpass.right.picos_val(count,3);
        ampVOTpass.left.deltasp2(count)=mean(ampVOTpass.left.picos_val(count,med))-ampVOTpass.left.picos_val(count,3);
        
        aux(1,count)=mean(ampVOTpass.left.diffs(count,med));
        aux2(1,count)=ampVOTpass.left.diffs(count,3);
        aux(2,count)=mean(ampVOTpass.right.diffs(count,med));
        aux2(2,count)=ampVOTpass.right.diffs(count,3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        box1(1,count)=mean(ampVOTpass.left.diffs(count,1:5));
        box1(2,count)=mean(ampVOTpass.right.diffs(count,1:5));
        
        ampVOTact.right.deltas(count)=mean(ampVOTact.right.diffs(count,med))-ampVOTact.right.diffs(count,3);
        ampVOTact.left.deltas(count)=mean(ampVOTact.left.diffs(count,med))-ampVOTact.left.diffs(count,3);
        ampVOTact.right.deltasp2(count)=mean(ampVOTact.right.picos_val(count,med))-ampVOTact.right.picos_val(count,3);
        ampVOTact.left.deltasp2(count)=mean(ampVOTact.left.picos_val(count,med))-ampVOTact.left.picos_val(count,3);
        aux3(1,count)=mean(ampVOTact.left.diffs(count,med));
        aux4(1,count)=ampVOTact.left.diffs(count,3);
        aux3(2,count)=mean(ampVOTact.right.diffs(count,med));
        aux4(2,count)=ampVOTact.right.diffs(count,3);
        
        box2(1,count)=mean(ampVOTact.left.diffs(count,1:5));
        box2(2,count)=mean(ampVOTact.right.diffs(count,1:5));
        
        ampFormpass.right.deltas(count)=mean(ampFormpass.right.diffs(count,med))-ampFormpass.right.diffs(count,3);
        ampFormpass.left.deltas(count)=mean(ampFormpass.left.diffs(count,med))-ampFormpass.left.diffs(count,3);
        ampFormpass.right.deltasp2(count)=mean(ampFormpass.right.picos_val(count,med))-ampFormpass.right.picos_val(count,3);
        ampFormpass.left.deltasp2(count)=mean(ampFormpass.left.picos_val(count,med))-ampFormpass.left.picos_val(count,3);
        aux5(1,count)=mean(ampFormpass.left.diffs(count,med));
        aux6(1,count)=ampFormpass.left.diffs(count,3);
        aux5(2,count)=mean(ampFormpass.right.diffs(count,med));
        aux6(2,count)=ampFormpass.right.diffs(count,3);
        
        box3(1,count)=mean(ampFormpass.left.diffs(count,1:5));
        box3(2,count)=mean(ampFormpass.right.diffs(count,1:5));
        
        ampFormact.right.deltas(count)=mean(ampFormact.right.diffs(count,med))-ampFormact.right.diffs(count,3);
        ampFormact.left.deltas(count)=mean(ampFormact.left.diffs(count,med))-ampFormact.left.diffs(count,3);
        ampFormact.right.deltasp2(count)=mean(ampFormact.right.picos_val(count,med))-ampFormact.right.picos_val(count,3);
        ampFormact.left.deltasp2(count)=mean(ampFormact.left.picos_val(count,med))-ampFormact.left.picos_val(count,3);
        aux7(1,count)=mean(ampFormact.left.diffs(count,med));
        aux8(1,count)=ampFormact.left.diffs(count,3);
        aux7(2,count)=mean(ampFormact.right.diffs(count,med));
        aux8(2,count)=ampFormact.right.diffs(count,3);
        
        box4(1,count)=mean(ampFormact.left.diffs(count,1:5));
        box4(2,count)=mean(ampFormact.right.diffs(count,1:5));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
b=polyfit(ampVOTpass.left.deltas,bt(1,:)',1);
subplot(1,2,1)
plot(ampVOTpass.left.deltas,bt(1,:),'o')
hold on
plot(ampVOTpass.left.deltas,polyval(b,ampVOTpass.left.deltas),'red')
hold off
[R2,p]=corrcoef(polyval(b,ampVOTpass.left.deltas),bt(1,:));
xlabel('Deltas (uV)')
ylabel('Psychometric Curve Slope')
title(sprintf('VOT Passive Left R2: %.4f p:%.4f', R2(2,1),p(2,1)))

format long
b=polyfit(ampVOTpass.right.deltas,bt(1,:)',1);
subplot(1,2,2)
plot(ampVOTpass.right.deltas,bt(1,:),'o')
hold on
plot(ampVOTpass.right.deltas,polyval(b,ampVOTpass.right.deltas),'red')
hold off
xlabel('Deltas (uV)')
ylabel('Psychometric Curve Slope')
[R2,p]=corrcoef(polyval(b,ampVOTpass.right.deltas),bt(1,:));
title(sprintf('VOT Passive Right p:%.4f', R2(2,1),p(2,1)))

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
[R2,p]=corrcoef(polyval(b,ampVOTact.left.deltas),bt(1,:));
title(sprintf('VOT Active Left R2: %.4f p:%.4f', R2(2,1),p(2,1)))


format long
b=polyfit(ampVOTact.right.deltas,bt(1,:)',1);
subplot(1,2,2)
plot(ampVOTact.right.deltas,bt(1,:),'o')
hold on
plot(ampVOTact.right.deltas,polyval(b,ampVOTact.right.deltas),'red')
hold off
xlabel('Deltas (uV)')
ylabel('Psychometric Curve Slope')
[R2,p]=corrcoef(polyval(b,ampVOTact.right.deltas),bt(1,:));
title(sprintf('VOT Active Right R2: %.4f p:%.4f', R2(2,1),p(2,1)))

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
[R2,p]=corrcoef(polyval(b,ampFormpass.left.deltas),bt(2,:));
title(sprintf('Formants Passive Left R2: %.4f p:%.4f', R2(2,1),p(2,1)))

format long
b=polyfit(ampFormpass.right.deltas,bt(2,:)',1);
subplot(1,2,2)
plot(ampFormpass.right.deltas,bt(2,:),'o')
hold on
plot(ampFormpass.right.deltas,polyval(b,ampFormpass.right.deltas),'red')
hold off
xlabel('Deltas (uV)')
ylabel('Psychometric Curve Slope')
[R2,p]=corrcoef(polyval(b,ampFormpass.right.deltas),bt(2,:));
title(sprintf('Formants Passive Right R2: %.4f p:%.4f', R2(2,1),p(2,1)))

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
[R2,p]=corrcoef(polyval(b,ampFormact.left.deltas),bt(2,:));
title(sprintf('Formants Active Left R2: %.4f p:%.4f', R2(2,1),p(2,1)))

format long
b=polyfit(ampFormact.right.deltas,bt(2,:)',1);
subplot(1,2,2)
plot(ampFormact.right.deltas,bt(2,:),'o')
hold on
plot(ampFormact.right.deltas,polyval(b,ampFormact.right.deltas),'red')
hold off
xlabel('Deltas (uV)')
ylabel('Psychometric Curve Slope')
[R2,p]=corrcoef(polyval(b,ampFormact.right.deltas),bt(2,:));
title(sprintf('Formants Active Right R2: %.4f p:%.4f', R2(2,1),p(2,1)))

%close all

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
col={'r','b','k','magenta','green','yellow','green','k','magenta','c','b'};
Leg=cell(5,1);

nome='VOT';
nome2='/AtivoFront';
subj=8;
L=size(Xica_DWT_rec,2);
t=(1:L)/fs;

[medL, medR, sdL, sdR]=normaliza(nome,subj);
figure(7)
set(gcf, 'Position',  [200, 200, 1300, 1000])
for i=1:5
    load(sprintf(strcat('../../Sujeito%d/',nome,nome2,'/ChanStimICA%d.mat'),subj,i))
    Xica_DWT_rec=Xica_DWT_rec(:,501:end);
    Right(:)=(Right(:)-medR)/sdR;
    Left(:)=(Left(:)-medL)/sdL;

    subplot(1,2,1)
    plot(t,-Left(:),col{i})
    title('Left')
    grid
    hold on
    
    subplot(1,2,2)
    plot(t,-Right(:),col{i})
    title('Right')
    Leg{i}=sprintf('%d',i);
    grid
    hold on
end
hold off
legend(Leg,'Location','northeast')

% Média das diferenças N1-P2
nome = ampFormact.left.diffs;
figure(8)
set(gcf, 'Position',  [200, 200, 1000, 1000])
for j=1:nb_subj
    plot(1:5,nome(j,:),strcat(col{j},'.'),'markersize',40)
    text(1-0.1,nome(j,1),num2str(sel(j)))
    hold on
    plot(1:5,nome(j,:),strcat(col{j},'-'))
end
plot(1:5,mean(nome),'magenta*')
plot(1:5,mean(nome),'magenta-')
title('Diffs N1-P2')
axis([0.5 5.5 -inf inf])
hold off
grid

% Média das diferenças P1-N1
nome=ampVOTact.right.diffsp1;
figure(9)
set(gcf, 'Position',  [200, 200, 1000, 1000])
for i=1:5
    for j=1:nb_subj
    plot(i,nome(j,i),strcat(col{i},'o'))
    hold on
    end
end
plot(1:5,mean(nome),'magenta*')
plot(1:5,mean(nome),'magenta')
title('Diffs P1-N1')
axis([0.5 5.5 -inf inf])
hold off
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
enderecos={'VOT/PassivoFront','VOT/AtivoFront','Formantes/PassivoFront','Formantes/AtivoFront'};
enderecos2={'VOT','VOT','Formantes','Formantes'};

auxDiffs(1:4,1:2,1:10,1:5)=0; %endereco, lado, linhas(picosP2, vales, 
%diffN1P2, picosP1, diffP1N1, N2_val, pico_lat, vale_lat, meio_lat, meio_val), stims
%Média geral por lado, condição e continuum 
figure(10)

for ed=1:4

XicaMat(1:2,1:nb_subj,1:L)=0; %lado,sujeito,sinal
set(gcf, 'Position',  [100, 100, 1200, 700])
for i=1:5
    for j=s
    [medL, medR, sdL, sdR]=normaliza(enderecos2{ed},j);
    load(sprintf('../../Sujeito%d/%s/ChanStimICA%d.mat',j,enderecos{ed},i))
    Xica_DWT_rec=Xica_DWT_rec(:,501:end);
    Left=mean([Xica_DWT_rec(14,:);Xica_DWT_rec(11,:);Xica_DWT_rec(5,:)]);
    Right=mean([Xica_DWT_rec(16,:);Xica_DWT_rec(8,:);Xica_DWT_rec(7,:)]);
    XicaMat(1,j,:)=-((Left-medR)/sdR);
    XicaMat(2,j,:)=-((Right-medR)/sdR);
    end
    
    subplot(1,2,1)
    Z=squeeze(XicaMat(1,:,:));
    aux=mean(Z);
    if (ed==4)
    plot(t,aux,col{i})
    grid
    title('Left')
    hold on
    end
    auxDiffs(ed,1,1,i)=find(aux(1,500:1500)==max(aux(1,500:1500)))+499;
    auxDiffs(ed,1,2,i)=find(aux(1,mN1:1000)==min(aux(1,mN1:1000)))+(mN1-1);
    auxDiffs(ed,1,7,i)=auxDiffs(ed,1,1,i)/fs;
    auxDiffs(ed,1,8,i)=auxDiffs(ed,1,2,i)/fs;
    auxDiffs(ed,1,3,i)= aux(1,auxDiffs(ed,1,1,i))-aux(1,auxDiffs(ed,1,2,i));
    auxDiffs(ed,1,9,i)= (auxDiffs(ed,1,2,i)+round((auxDiffs(ed,1,1,i)-auxDiffs(ed,1,2,i))/2))/fs;
    auxDiffs(ed,1,10,i)= aux((auxDiffs(ed,1,2,i)+round((auxDiffs(ed,1,1,i)-auxDiffs(ed,1,2,i))/2)));
    auxDiffs(ed,1,4,i)=(find(aux(1,50:550)==max(aux(1,50:550))))+49;
    auxDiffs(ed,1,5,i)= aux(1,auxDiffs(ed,1,4,i))-aux(1,auxDiffs(ed,1,2,i));
    auxDiffs(ed,1,6,i)=(find(aux(1,1300:1900)==min(aux(1,1300:1900))))+1299;
    auxDiffs(ed,1,1,i)=aux(find(aux(1,500:1500)==max(aux(1,500:1500)))+499);
    auxDiffs(ed,1,2,i)=aux(find(aux(1,mN1:1000)==min(aux(1,mN1:1000)))+(mN1-1));
    auxDiffs(ed,1,4,i)=aux(find(aux(1,50:550)==max(aux(1,50:550)))+49);
    auxDiffs(ed,1,6,i)=aux(find(aux(1,1300:1900)==max(aux(1,1300:1900)))+1299);
    
    subplot(1,2,2)
    Z=squeeze(XicaMat(2,:,:));
    aux=mean(Z);
    if (ed==4)
    plot(t,aux,col{i})
    grid
    title('Right')
    Leg{i}=sprintf('%d',i);
    hold on
    end
    auxDiffs(ed,2,1,i)=find(aux(1,500:1500)==max(aux(1,500:1500)))+499;
    auxDiffs(ed,2,2,i)=find(aux(1,mN1:1000)==min(aux(1,mN1:1000)))+(mN1-1);
    auxDiffs(ed,2,7,i)=auxDiffs(ed,2,1,i)/fs;
    auxDiffs(ed,2,8,i)=auxDiffs(ed,2,2,i)/fs;
    auxDiffs(ed,2,9,i)= (auxDiffs(ed,2,2,i)+round((auxDiffs(ed,2,1,i)-auxDiffs(ed,2,2,i))/2))/fs;
    auxDiffs(ed,2,10,i)= aux((auxDiffs(ed,2,2,i)+round((auxDiffs(ed,2,1,i)-auxDiffs(ed,2,2,i))/2)));
    auxDiffs(ed,2,3,i)= aux(1,auxDiffs(ed,2,1,i))-aux(1,auxDiffs(ed,2,2,i));
    auxDiffs(ed,2,4,i)=find(aux(1,50:550)==max(aux(1,50:550)))+49;
    auxDiffs(ed,2,5,i)= aux(1,auxDiffs(ed,2,4,i))-aux(1,auxDiffs(ed,2,2,i));
    auxDiffs(ed,2,1,i)=aux(find(aux(1,500:1500)==max(aux(1,500:1500)))+499);
    auxDiffs(ed,2,2,i)=aux(find(aux(1,mN1:1000)==min(aux(1,mN1:1000)))+(mN1-1));
    auxDiffs(ed,2,4,i)=aux(find(aux(1,50:550)==max(aux(1,50:550)))+49);
    auxDiffs(ed,2,6,i)=aux(find(aux(1,1300:1900)==max(aux(1,1300:1900)))+1299);
end
end
hold off
legend(Leg,'Location','northeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plota N1-P2 das curvas médias da figura anterior
nomes={'VOTpass','VOTact','FormPass','FormAct'};
cor={'magenta','blue','red','black'};
simb={'diamond','square','*','pentagram','hexagram'};
figure(11)
set(gcf, 'Position',  [200, 200, 1000, 1000])
for ed=1:4
aux10=squeeze(auxDiffs(ed,1,3,:));
aux20=squeeze(auxDiffs(ed,2,3,:));
for i=1:5
    subplot(1,2,1)
    plot(i,aux10(i),strcat(col{i},simb{i}),'markersize',10,'MarkerFaceColor',col{i})
    title('Left')
    grid
    hold on
    
    subplot(1,2,2)
    plot(i,aux20(i),strcat(col{i},simb{i}),'markersize',10,'MarkerFaceColor',col{i})
    title('Right')
    Leg{i}=sprintf('Stim. %d',i);
    grid
    hold on
end
subplot(1,2,1)
plot(1:5,aux10,cor{ed})
text(5.05,aux10(5),nomes{ed},'Color',cor{ed})
title('Left')
ylabel('Mean N1-P2 amplitude (uV)')
xlabel('Stimulus')
axis([0.5 5.5 1.1 4.2])

subplot(1,2,2)
plot(1:5,aux20,cor{ed})
text(5.05,aux20(5),nomes{ed},'Color',cor{ed})
title('Right')
ylabel('Mean N1-P2 amplitude (uV)')
xlabel('Stimulus')
axis([0.5 5.5 1.1 4.2])
end
hold off
subplot(1,2,1)
grid
subplot(1,2,2)
grid
legend(Leg,'Location','best')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
plot(1:5,aux10,cor{ed})
text(5.05,aux10(5),nomes{ed})
title('P1 all')
subplot(1,2,2)
plot(1:5,aux20,cor{ed})
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
plot(1:5,aux10,cor{ed})
text(5.05,aux10(5),nomes{ed})
title('P1-N1 all')
subplot(1,2,2)
plot(1:5,aux20,cor{ed})
text(5.05,aux20(5),nomes{ed})
title('P1-N1 all')
end
hold off

legend(Leg,'Location','best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
plot(1:5,aux10,cor{ed})
text(5.05,aux10(5),nomes{ed})
axis([0.5 5.5 0 1.8])
title('N1 all')
subplot(1,2,2)
plot(1:5,aux20,cor{ed})
text(5.05,aux20(5),nomes{ed})
axis([0.5 5.5 0 1.8])
title('N1 all')
end
hold off
subplot(1,2,1)
grid
subplot(1,2,2)
grid
legend(Leg,'Location','best')       
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(15)
%Plota P2 médio de acordo com estímulo
set(gcf, 'Position',  [200, 200, 1000, 1000])
cor={'magenta','blue','red','black'};
simb={'diamond','square','*','pentagram','hexagram'};
for ed=1:4
aux10=squeeze(auxDiffs(ed,1,1,:));
aux20=squeeze(auxDiffs(ed,2,1,:));
for i=1:5
    subplot(1,2,1)
    plot(i,aux10(i),strcat(col{i},simb{i}),'markersize',10,'MarkerFaceColor',col{i})
    title('Left')
    hold on
    
    subplot(1,2,2)
    plot(i,aux20(i),strcat(col{i},simb{i}),'markersize',10,'MarkerFaceColor',col{i})
    title('Right')
    Leg{i}=sprintf('%d',i);
    hold on
    
end
subplot(1,2,1)
grid
plot(1:5,aux10,cor{ed})
axis([0.5 5.5 0.7 2.4])
text(5.05,aux10(5),nomes{ed})
title('P2 all')

subplot(1,2,2)
grid
plot(1:5,aux20,cor{ed})
axis([0.5 5.5 0.7 2.4])
text(5.05,aux20(5),nomes{ed})
title('P2 all')
end
hold off
hold off
subplot(1,2,1)
grid
subplot(1,2,2)
grid
legend(Leg,'Location','best')       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(16)
%Plota N2 médio de acordo com estímulo
set(gcf, 'Position',  [200, 200, 1000, 1000])
for ed=1:4
aux10=abs(squeeze(auxDiffs(ed,1,6,:)));
aux20=abs(squeeze(auxDiffs(ed,2,6,:)));
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
plot(1:5,aux10,cor{ed})
text(5.05,aux10(5),nomes{ed},'Color',cor{ed})
title('N2 all')
subplot(1,2,2)
plot(1:5,aux20,cor{ed})
text(5.05,aux20(5),nomes{ed},'Color',cor{ed})
title('N2 all')
end
hold off
legend(Leg,'Location','best')  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Primeiro média dos estímulos e aqui o box une os sujeitos
%Comentado está um método em que tem-se primeiro a média dos sujeitos e o
%box une os estímulos.
figure(17)
subplot(1,2,1)
%boxplot([squeeze(auxDiffs(1,1,3,:)),squeeze(auxDiffs(2,1,3,:)), squeeze(auxDiffs(3,1,3,:)),squeeze(auxDiffs(4,1,3,:))],'Labels',{'V pass','V act','F pass','F act'})
boxplot([box1(1,:)',box2(1,:)', box3(1,:)',box4(1,:)'],'Labels',{'V pass','V act','F pass','F act'})
ylabel('Amplitude N1-P2 (uV)')
title('Mean stimuli 1,2,3,4,5 V=VOT, F=Formants LEFT')
grid

%figure(18)
subplot(1,2,2)
boxplot([box1(2,:)',box2(2,:)', box3(2,:)',box4(2,:)'],'Labels',{'V pass','V act','F pass','F act'})
ylabel('Amplitude N1-P2 (uV)')
title('Mean stimuli 1,2,3,4,5 V=VOT, F=Formants RIGHT')
grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 box1(1,:)=mean(ampVOTpass.left.picos_lat,2)';
 box1(2,:)=mean(ampVOTpass.right.picos_lat,2)';
   
 box2(1,:)=mean(ampVOTact.left.picos_lat,2)';
 box2(2,:)=mean(ampVOTact.right.picos_lat,2)';
    
 box3(1,:)=mean(ampFormpass.left.picos_lat,2)';
 box3(2,:)=mean(ampFormpass.right.picos_lat,2)';
    
 box4(1,:)=mean(ampFormact.left.picos_lat,2)';
 box4(2,:)=mean(ampFormact.right.picos_lat,2)';


figure(19)

boxplot([box1(1,:)',box2(1,:)', box3(1,:)',box4(1,:)'],'Labels',{'V pass','V act','F pass','F act'})
ylabel('Latency P2 (s)')
title('Mean stimuli 1,2,3,4,5 V=VOT, F=Formants LEFT')
grid

figure(20)

boxplot([box1(2,:)',box2(2,:)', box3(2,:)',box4(2,:)'],'Labels',{'V pass','V act','F pass','F act'})
ylabel('Latency P2 (s)')
title('Mean stimuli 1,2,3,4,5 V=VOT, F=Formants RIGHT')
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 box1(1,:)=mean(ampVOTpass.left.vales_lat,2)';
 box1(2,:)=mean(ampVOTpass.right.vales_lat,2)';
    
 box2(1,:)=mean(ampVOTact.left.vales_lat,2)';
 box2(2,:)=mean(ampVOTact.right.vales_lat,2)';
    
 box3(1,:)=mean(ampFormpass.left.vales_lat,2)';
 box3(2,:)=mean(ampFormpass.right.vales_lat,2)';
    
 box4(1,:)=mean(ampFormact.left.vales_lat,2)';
 box4(2,:)=mean(ampFormact.right.vales_lat,2)';

figure(21)

boxplot([box1(1,:)',box2(1,:)', box3(1,:)',box4(1,:)'],'Labels',{'V pass','V act','F pass','F act'})
ylabel('Latency N1 (s)')
title('Mean stimuli 1,2,3,4,5 V=VOT, F=Formants LEFT')
grid

figure(22)

boxplot([box1(2,:)',box2(2,:)', box3(2,:)',box4(2,:)'],'Labels',{'V pass','V act','F pass','F act'})
ylabel('Latency N1 (s)')
title('Mean stimuli 1,2,3,4,5 V=VOT, F=Formants Right')
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Relação com beta com left and right juntos
figure(23)
set(gcf, 'Position',  [100, 100, 400, 900])
nbeta=[bt(2,:) bt(2,:)]; %bt(1,:) para VOT e bt(2,:) para formantes
for s=1:5
nome=[ampFormact.left.vales_val(:,s);ampFormact.right.vales_val(:,s)];
subplot(5,1,s)
format long
b=polyfit(nome,nbeta',1);
plot(nome,nbeta,'o')
hold on
plot(nome,polyval(b,nome),'red')
hold off
[R2,p]=corrcoef(nome,nbeta');
xlabel('P2 (uV)')
ylabel('Beta')
title(sprintf('VOT active R2: %.4f p-value:%.4f', R2(2,1),p(2,1)))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(24)
%Plota latência P2 média de acordo com estímulo
set(gcf, 'Position',  [200, 200, 1000, 1000])
cor={'magenta','blue','red','black'};
simb={'diamond','square','*','pentagram','hexagram'};
n={ampVOTpass, ampVOTact, ampFormpass, ampFormact,};
for ed=1:4
aux10=abs(squeeze(auxDiffs(ed,1,7,:)));
aux20=abs(squeeze(auxDiffs(ed,2,7,:)));
%aux10=mean(n{ed}.left.picos_lat);
%aux20=mean(n{ed}.right.picos_lat);
for i=1:5
    subplot(1,2,1)
    plot(i,aux10(i),strcat(col{i},simb{i}),'markersize',10,'MarkerFaceColor',col{i})
    title('Left')
    grid
    hold on
    
    subplot(1,2,2)
    plot(i,aux20(i),strcat(col{i},simb{i}),'markersize',10,'MarkerFaceColor',col{i})
    title('Right')
    Leg{i}=sprintf('Stim. %d',i);
    grid
    hold on
end
subplot(1,2,1)
plot(1:5,aux10,cor{ed})
axis([0.5 5.5 0.14 0.235])
text(5.05,aux10(5),nomes{ed},'Color',cor{ed})
ylabel('Mean P2 Latency (s)')
xlabel('Stimulus')
title('Left')
subplot(1,2,2)
plot(1:5,aux20,cor{ed})
axis([0.5 5.5 0.14 0.235])
text(5.05,aux20(5),nomes{ed},'Color',cor{ed})
ylabel('Mean P2 Latency (s)')
xlabel('Stimulus')
title('Right')
end
hold off
subplot(1,2,1)
grid
subplot(1,2,2)
grid
legend(Leg,'Location','best')  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(25)
%Plota latência N1 média de acordo com estímulo
set(gcf, 'Position',  [200, 200, 1000, 1000])
cor={'magenta','blue','red','black'};
for ed=1:4
aux10=abs(squeeze(auxDiffs(ed,1,8,:)));
aux20=abs(squeeze(auxDiffs(ed,2,8,:)));
for i=1:5
    subplot(1,2,1)
    plot(i,aux10(i),strcat(col{i},simb{i}),'markersize',10,'MarkerFaceColor',col{i})
    title('Left')
    grid
    hold on
    
    subplot(1,2,2)
    plot(i,aux20(i),strcat(col{i},simb{i}),'markersize',10,'MarkerFaceColor',col{i})
    title('Right')
    Leg{i}=sprintf('Stim. %d',i);
    grid
    hold on
end
subplot(1,2,1)
plot(1:5,aux10,cor{ed})
text(5.05,aux10(5),nomes{ed},'Color',cor{ed})
ylabel('Mean N1 Latency (s)')
xlabel('Stimulus')
title('Left')
axis([0.5 5.5 0.07 0.16])
subplot(1,2,2)
plot(1:5,aux20,cor{ed})
text(5.05,aux20(5),nomes{ed},'Color',cor{ed})
ylabel('Mean N1 Latency (s)')
xlabel('Stimulus')
title('Right')
axis([0.5 5.5 0.07 0.16])
end
hold off
subplot(1,2,1)
grid
subplot(1,2,2)
grid
legend(Leg,'Location','best')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Relação com beta com left and right juntos

figure(26)
set(gcf, 'Position',  [100, 100, 400, 900])
nomeL=ampFormact.left.vales_val;
nomeR=ampFormact.right.vales_val;
a=2;

nbeta=[bt(a,:) bt(a,:) bt(a,:) bt(a,:)]; %bt(1,:) para VOT e bt(2,:) para formantes
claros=[nomeL(:,1);nomeR(:,1);nomeL(:,5);nomeR(:,5)];
subplot(3,1,1)
format long
b=polyfit(claros,nbeta',1);
plot(claros,nbeta','o')
hold on
plot(claros,polyval(b,claros),'red')
hold off
[R2,p]=corrcoef(claros,nbeta');
xlabel('Stim1+Stim5 (uV)')
ylabel('Beta')
title(sprintf('VOT active R2: %.4f p-value:%.4f', R2(2,1),p(2,1)))

ambiguos=[nomeL(:,2);nomeR(:,2);nomeL(:,4);nomeR(:,4)];
subplot(3,1,2)
format long
b=polyfit(ambiguos,nbeta',1);
plot(ambiguos,nbeta','o')
hold on
plot(ambiguos,polyval(b,ambiguos),'red')
hold off
[R2,p]=corrcoef(ambiguos,nbeta');
xlabel('Stim2+Stim4 (uV)')
ylabel('Beta')
title(sprintf('VOT active R2: %.4f p-value:%.4f', R2(2,1),p(2,1)))

nbeta=[bt(a,:) bt(a,:)];
amb=[nomeL(:,3);nomeR(:,3)];
subplot(3,1,3)
format long
b=polyfit(amb,nbeta',1);
plot(amb,nbeta,'o')
hold on
plot(amb,polyval(b,amb),'red')
hold off
[R2,p]=corrcoef(amb,nbeta');
xlabel('Stim3 (uV)')
ylabel('Beta')
title(sprintf('VOT active R2: %.4f p-value:%.4f', R2(2,1),p(2,1)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Relação com beta com left and right juntos
figure(27)
set(gcf, 'Position',  [100, 100, 400, 900])
nomeL=ampFormpass.left.vales_val;
nomeR=ampFormpass.right.vales_val;
a=2; %1 para VOT e 2 para Form

nbeta=[bt(a,:) bt(a,:)]; %bt(1,:) para VOT e bt(2,:) para formantes
claros=[mean([nomeL(:,1)';nomeL(:,2)';nomeL(:,4)';nomeL(:,5)']),mean([nomeR(:,1)';nomeR(:,2)';nomeR(:,4)';nomeR(:,5)'])];
subplot(2,1,1)
format long
b=polyfit(claros,nbeta,1);
plot(claros,nbeta,'o')
hold on
plot(claros,polyval(b,claros),'red')
hold off
[R2,p]=corrcoef(claros,nbeta);
xlabel('Stim1~Stim5 (uV)')
ylabel('Beta')
title(sprintf('Form passive R2: %.4f p-value:%.4f', R2(2,1),p(2,1)))

nbeta=[bt(a,:) bt(a,:)];
amb=[nomeL(:,3);nomeR(:,3)];
subplot(2,1,2)
format long
b=polyfit(amb,nbeta',1);
plot(amb,nbeta,'o')
hold on
plot(amb,polyval(b,amb),'red')
hold off
[R2,p]=corrcoef(amb,nbeta');
xlabel('Stim3 (uV)')
ylabel('Beta')
title(sprintf('Form passive R2: %.4f p-value:%.4f', R2(2,1),p(2,1)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%orient(figure(15),'landscape')
%print(figure(15),'P2.pdf','-dpdf')