%Creates mat files for test of RHD data with old scripts from LPNC
%Adrielle de Carvalho Santana
%Rodar com passivo.rhd
%22/05/2019

clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('idx_passive.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num=sort(abs(unique(idx)));

read_Intan_RHD2000_file2('passivo_190529_110950.rhd')

amplifier=[amplifier_data(1:11,:); zeros(1,size(amplifier_data,2)); amplifier_data(12:15,:); zeros(1,size(amplifier_data,2)); amplifier_data(16:17,:)];

left0=amplifier(13,:)-amplifier(19,:);
right0=amplifier(13,:) - amplifier(3,:);
stim0=amplifier(1,:);

rate=frequency_parameters.amplifier_sample_rate;

left=left0(1,[1:2145167 2146732:end]);
right=right0(1,[1:2145167 2146732:end]);
stim=stim0(1,[1:2145167 2146732:end]);
amplifier=amplifier(:,[1:2145167 2146732:end]);

left=notch60(left,rate);
right=notch60(right,rate);


% Find stimuli
%floor=find(stim>30);

% Find stimuli
%count=1;
%samp=[floor(1)];

%for i=2:length(floor)
 %   if ((floor(i)-floor(i-1))>300)
 %      samp=[samp floor(i)];
 %      count=count+1;
  %  end
%end

%plot(stim)
%hold on
%stem(samp, stim(samp), 'red')

%samp= [samp(1,1:1605) 9251947 9257522 9263558 9269129 9275168 9281205 9286328 9291886 9297923 9303037 9308487 9314646 9320679 9326262 9332290 9338308 9343435 9349042 9354135 9359688 9366197 9372213 9377796 9382876 9475717 9481361 9487396 9493432 samp(1,1662:end)];
load('sampPass.mat')

idx=[idx(1,1:1629) idx(1,1645:end)];

t=-0.1:1/rate:(0.6-1/rate);
colors=colormap(lines);

Leg=cell(length(num),1);

aux=0;
n=size(amplifier,1)-1;

data_size=4096; %múltiplo de 2 para usar na DWT
t_size=4096/rate; %em termos de tempo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
canal=left;
top=find(abs(canal)>80); %Oz channel
bad=[];
for j=2:length(top)
    if (top(j)-top(j-1)>600)
        x=(top(j)-2500):(top(j)+510);
        for k=1:length(samp)
            g=find(x==samp(k));
            if g>0
                bad=[bad x(g)];
            end
        end
    end
end

for i=1:length(samp)
    if (abs(canal(samp(i)))>38)
        bad=[bad samp(i)];
    end
end
bad=[bad 9087552 9098698 9359688 9343435];
%In yellow are trials that stay and in red the ones which are removed
close all
plot(canal)
hold on
stem(samp, canal(samp), 'yellow')
stem(bad, canal(bad), 'red')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l1=zeros(1,5);
c=0;
for i=1:length(num)
    if mod(i,2)==1
        c=c+1;
        st=samp(find(idx==num(i) | idx==-num(i)));
                       
        for k=1:length(bad)
            st=st(find(st~=bad(k)));
        end
        l1(c)=length(st);
    end
end
            
for i=1:length(num)
    if mod(i,2)==1
            aux=aux+1;
            chan_stim(1:n,1:(t_size*rate))=0;
            st=samp(find(idx==num(i) | idx==-num(i)));
                       
            for k=1:length(bad)
                st=st(find(st~=bad(k)));
            end
            
            if (length(st)>min(l1))
                dif=length(st)-min(l1);
                st=st(randperm(length(st)));
                st=st(1:(end-dif));
            end
               
            reps=length(st);
            channels(1:n,1:(t_size*rate),1:reps)=0;
            
            for c=1:n
                for p=1:length(st)
                    channels(c,:,p)=amplifier((c+1),(st(p)-(0.1*rate)):(st(p)+((t_size-0.1)*rate)-1));
                end
                chan_stim(c,:)=mean(channels(c,:,:),3);
            end
            chan_stim=notch60(chan_stim,rate);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    name=sprintf('Passivo/ChanStim%g',aux);
    save(name,'chan_stim','rate')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l1=zeros(1,10);
c=0;
for i=1:length(num)
        c=c+1;
        s=samp(find(idx==num(i)));
                       
        for k=1:length(bad)
            s=s(find(s~=bad(k)));
        end
        l1(c)=length(s);
 end

for i=1:length(num)
    
    s=samp(find(idx==num(i)));
    for k=1:length(bad)
          s=s(find(s~=bad(k)));
    end
    
    if (length(s)>min(l1))
           dif=length(s)-min(l1);
           s=s(randperm(length(s)));
           s=s(1:(end-dif));
    end

    reps=length(s);
    
    samples(1:2,1:(0.7*rate),1:reps)=0;
    
    for k=1:length(s)
        samples(1,:,k)=left((s(k)-(0.1*rate)):(s(k)+(0.6*rate)-1));
        samples(2,:,k)=right((s(k)-(0.1*rate)):(s(k)+(0.6*rate)-1));
    
        base=mean(samples(1,1:(0.1*rate),k));
        samples2(1,:,k)=left(s(k):(s(k)+(t_size*rate)-1))-base;
        base=mean(samples(2,1:(0.1*rate),k));
        samples2(2,:,k)=right(s(k):(s(k)+(t_size*rate)-1))-base;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    name=sprintf('Passivo/Stim%g',i);
    name2=sprintf('Passivo/Stim2_%g',i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save(name,'samples')
    save(name2,'samples2')
    
    
    left_base=mean(mean(samples(1,1:(0.1*rate),:),3));
    right_base=mean(mean(samples(2,1:(0.1*rate),:),3));
    
    figure(1)
    subplot(2,5,i)
    plot(t,mean(samples(1,:,:),3)-left_base)
    title(sprintf('left %d',num(i)))
    grid  
    
    figure(2)
    subplot(2,5,i)
    plot(t,mean(samples(2,:,:),3)-right_base)
    title(sprintf('right %d',num(i)))
    grid 
    
    figure(3)
    plot(t,mean(samples(1,:,:),3),'color',[colors(i,:)])
    title('Left')
    Leg{i}=sprintf('%d',num(i));
    hold on

end
figure(3)
hold off
grid
legend(Leg,'Location','northeast')

%Plota as médias da parte positiva com a negativa e filtra em 80Hz

load('low.mat');

count=0;
for j=1:length(num)
    
    if mod(j,2)==1
        count=count+1;
        load(sprintf('Passivo/Stim%d.mat',j))
        baseline=mean(mean(samples(1,1:(0.1*rate),:),3));
        s1_left=mean(samples(1,:,:),3)-baseline;
        
        baseline=mean(mean(samples(2,1:(0.1*rate),:),3));
        s1_right=mean(samples(2,:,:),3)-baseline;
        
        load(sprintf('Passivo/Stim%d.mat',j+1))
        baseline=mean(mean(samples(1,1:(0.1*rate),:),3));
        s2_left=mean(samples(1,:,:),3)-baseline;
        
        baseline=mean(mean(samples(2,1:(0.1*rate),:),3));
        s2_right=mean(samples(2,:,:),3)-baseline;
        
        m_left=(s1_left+s2_left)/2;
        
        m_right=(s1_right+s2_right)/2;
        
        m_filt_left=dtrend(filtfilt(L.sosMatrix, L.ScaleValues, m_left));
        m_filt_right=dtrend(filtfilt(L.sosMatrix, L.ScaleValues, m_right));
        
        figure(4)
        subplot(2,3,count)
        plot(t,m_filt_left)
        title(sprintf('left %d',num(j)))
        grid 
        
        figure(5)
        subplot(2,3,count)
        plot(t,m_filt_right)
        title(sprintf('right %d',num(j)))
        grid 
        
        figure(6)
        plot(t,m_left, 'color',[colors(count,:)])
        title('Left sem filtro')
        Leg2{count}=sprintf('%d',num(j));
        hold on
        
        figure(7)
        plot(t,m_right, 'color',[colors(count,:)])
        title('Right sem filtro')
        hold on
        
        figure(8)
        plot(t,m_filt_left, 'color',[colors(count,:)])
        title('Left com filtro')
        hold on
        
        figure(9)
        plot(t,m_filt_right, 'color',[colors(count,:)])
        title('Right com filtro')
        hold on
    
    end
end

figure(6)
hold off
grid
legend(Leg2,'Location','northeast')
figure(7)
hold off
grid
legend(Leg2,'Location','northeast')
figure(8)
hold off
grid
legend(Leg2,'Location','northeast')
figure(9)
hold off
grid
legend(Leg2,'Location','northeast')

