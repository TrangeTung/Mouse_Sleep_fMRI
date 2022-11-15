function F = MY_get_figure_Task_related_Displacement(rp,temporal_resolution,y_axis_min_max,legend_txt)
%% rp : cell{rp} ; Multiple Headmotion datasets
% F: Figure

F = figure;
suptitle('Task related Displacement');
set(gca,'FontWeight','bold');
set(F,'position',[0 0 900 900]);

color = {[138 043 226];[002 152 182];[135 199 000];[248 199 000];[005 139 085];[000 000 255]};
rand = [1 2 4 8];
%% rp
for count = 1:numel(rp)
    rp_displacement = rp{count};
    for kk = 1:4
        stimulus = rand(kk);
        eval(['TimeCourse_mean = transpose(rp_displacement.stimu_',num2str(stimulus),'sec);']);
        number = size(TimeCourse_mean,2);
        Time_series = mean(TimeCourse_mean,2);
        upper = Time_series+std(TimeCourse_mean,0,2)/sqrt(number);
        lower = flipud(Time_series-std(TimeCourse_mean,0,2)/sqrt(number));
    
        subplot(numel(rp),numel(rand),(count-1)*numel(rand)+kk);
        f = plot([1:numel(Time_series)]*temporal_resolution,Time_series,'color',color{count}/255,'LineWidth',1);
        x = [1:numel(Time_series),fliplr(1:numel(Time_series))]';
        y = [upper;lower];
        patch(x*temporal_resolution,y,color{count}/255,'FaceAlpha',.2,'EdgeColor','none');  
        hold on;
        
        initial_time_point = 2/temporal_resolution;
        t_control = TimeCourse_mean(initial_time_point,:);
        clear sig;
        for timepoint = initial_time_point:(2+stimulus+5)/temporal_resolution
            [h,~,~,~] = ttest2(t_control,TimeCourse_mean(timepoint,:));
            sig(timepoint) = h;
        end
        sig_x = find(sig==1)-1;
        sig_y = upper(sig==1);
        text(sig_x*temporal_resolution,sig_y,'*');

        if kk == 1
            ylabel('Displacement (A.U.)');
        end
        if count == numel(rp)
            xlabel('Time (s)');
        end
        set(gca,'FontSize',8,'FontWeight','bold');
        set(gca, 'LineWidth',1.5);
        y_min = y_axis_min_max(count,1);
        y_max = y_axis_min_max(count,2);
        %% stimulus
        patch_x = [2+[0 stimulus stimulus 0]]/temporal_resolution';
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
text(10,0.255,'(Mean¡ÀSEM)','FontWeight','bold')
end