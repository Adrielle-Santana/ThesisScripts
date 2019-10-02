%Plot stimuli vs response time for all subjects and the mean betweet them
%for the VOT and Formant conditions (active task)
clear all

for ad=1:2
    x(1:11,1:5)=0;
    y(1:11,1:5)=0;
    for i=1:11
        address={sprintf('../Sujeito%d/VOT/',i), sprintf('../Sujeito%d/Formantes/',i)};

        if (i==2 && ad==2)
            load(strcat(address{ad},'rt1.mat'))
            auxrt=rt;
            load(strcat(address{ad},'idx_active1.mat'))
            auxidx=idx;
            load(strcat(address{ad},'rt2.mat'))
            rt=[auxrt rt];
            load(strcat(address{ad},'idx_active2.mat'))
            idx=[auxidx idx];
            
        elseif (i==6 && ad==2)
            load(strcat(address{ad},'rt.mat'))
            auxrt=rt;
            load(strcat(address{ad},'idx_active.mat'))
            auxidx=idx;
            load(strcat(address{ad},'rt50.mat'))
            rt=[auxrt rt];
            load(strcat(address{ad},'idx_active50.mat'))
            idx=[auxidx idx];
            
        elseif(i==8 && ad==2)
            load(strcat(address{ad},'rt1.mat'))
            auxrt=rt;
            load(strcat(address{ad},'idx_active1.mat'))
            auxidx=idx;
            load(strcat(address{ad},'rt2.mat'))
            rt=[auxrt rt];
            load(strcat(address{ad},'idx_active2.mat'))
            idx=[auxidx idx];
            
        else
            load(strcat(address{ad},'rt.mat'))
            load(strcat(address{ad},'idx_active.mat'))
        end
        
        if (i==10 && ad==2)
            x(i,:)=[1 sort(unique(abs(idx)))];
        else
            x(i,:)=sort(unique(abs(idx)));
        end
        
        for k=1:5
            st=rt(find(idx==x(i,k) | idx==-x(i,k)));
            if (ad==1)
                y(i,k)=mean(st)-0.22;
            else
                y(i,k)=mean(st)-0.191;
            end
        end
        
        stem(x(i,:),y(i,:),':b*')
        ylabel('Time (s)')
        xlabel('Stimulus continuum number')
        title(sprintf('Response Time Subject%d', i))
        hold on
        plot(x(i,:),y(i,:),'--red')
        hold off
                
        saveas(gcf,strcat(address{ad},'resp_time.png'))
        
    end
    
        stem(mean(x),mean(y),':b*')
        ylabel('Time (s)')
        xlabel('Stimulus continuum number')
        title('Mean response time all subjects')
        hold on
        plot(mean(x),mean(y),'--red')
        hold off
        
        if (ad==1)
            saveas(gcf,'resp_time_VOT.png')
        else
            saveas(gcf,'resp_time_Form.png')
        end
end