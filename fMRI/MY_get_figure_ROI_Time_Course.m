function F = MY_get_figure_ROI_Time_Course(Time_course,temporal_resolution,legend_txt,title_txt)
%% Time_course : cell{Time_course} ; Multiple Time_course datasets
% F: Figure

F = figure;
suptitle(title_txt);
set(gca,'FontWeight','bold');
set(F,'position',[100 100 900 900/5]);

color = {[138 043 226];[002 152 182];[135 199 000];[248 199 000];[005 139 085];[000 000 255]};
rand = [1 2 4 8];
%% rp
for count = 1:numel(Time_course)
    serial = Time_course{count};
    for kk = 1:4
        stimulus = rand(kk);
        eval(['TimeCourse_mean = transpose(serial.stimu_',num2str(stimulus),'sec);']);
        number = size(TimeCourse_mean,2);
        Time_series = mean(TimeCourse_mean,2);
        initial_time_point = 2/temporal_resolution;
        upper = Time_series+std(TimeCourse_mean,0,2)/sqrt(number);
        lower = flipud(Time_series-std(TimeCourse_mean,0,2)/sqrt(number));

        subplot(numel(Time_course),numel(rand),(count-1)*numel(rand)+kk);
        axis_index = [1:numel(Time_series)]*temporal_resolution;
        f = plot(axis_index,Time_series,'color',color{count}/255,'LineWidth',1);
        hold on;
        plot(axis_index,axis_index*0,'k--','LineWidth',0.5);
        x = [1:numel(Time_series),fliplr(1:numel(Time_series))]';
        y = [upper;lower];
        patch(x*temporal_resolution,y,color{count}/255,'FaceAlpha',.2,'EdgeColor','none');  
        hold on;
        
        t_control = TimeCourse_mean(initial_time_point,:);
        clear sig;
        for timepoint = initial_time_point:numel(Time_series)%(2+stimulus+13)/temporal_resolution
            [h,~,~,~] = ttest2(t_control,TimeCourse_mean(timepoint,:));
            sig(timepoint) = h;
        end
        sig_x = find(sig==1)-1;
        sig_y = upper(sig==1);
%         text(sig_x*temporal_resolution,sig_y,'*');

        if kk == 1
            ylabel('BOLD response (%)');
        end
        if count == numel(Time_course)
            xlabel('Time (s)');
        end
        set(gca,'FontSize',8,'FontWeight','bold');
        set(gca, 'LineWidth',1.5);
        
        y_min =-0.5;
        y_max = 2.5;
        %% stimulus
        stimulus_RAM = repmat([zeros(1,2) ones(1,stimulus) zeros(1,13)],[1/temporal_resolution,1]);
        stimulus_T = reshape(stimulus_RAM,[1,numel(stimulus_RAM)]);
        stimulus_deT = stimulus_T - [stimulus_T(2:end) stimulus_T(1)];
        patch_x = [find(stimulus_deT==-1);find(stimulus_deT==1);find(stimulus_deT==1);find(stimulus_deT==-1)];
        patch_y = [y_min y_min y_max y_max]';
        patch(patch_x*temporal_resolution,patch_y,'r','FaceAlpha',.2,'EdgeColor','none'); 
        hold off;
        axis([0*temporal_resolution numel(Time_series)*temporal_resolution y_min y_max])
        if kk == 4
            h = legend(f,legend_txt{count});
            pos = 0.20 + 0.12*(count-1);
            set(h,'Position',[pos 0.03 0.070 0.05]);
            set(h,'Box','off');
        end
        box off;
    end
end
text(10,-1.8,'(Mean¡ÀSEM)','FontWeight','bold')
end