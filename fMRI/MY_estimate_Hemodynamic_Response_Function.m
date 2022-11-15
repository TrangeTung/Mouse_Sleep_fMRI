function [F,HRF] = MY_estimate_Hemodynamic_Response_Function(Time_series,temporal_resolution) %#ok<INUSL>
%% Time_series : the input struct Array
% HRF : the output
color = {[138 043 226];[002 152 182];[135 199 000];[248 199 000];[005 139 085];[000 000 255]};
%% Signal extraction
stimulus = [];
for stimu = [1 2 4 8]
    eval(['RAM = Time_series.stimu_',num2str(stimu),'sec;'])
    RAM(isnan(RAM))=0;
    Time_series_RAM = mean(RAM,1);
    initial_time_point = 2/temporal_resolution;
%     Time_series_RAM = Time_series_RAM - Time_series_RAM(initial_time_point);
    eval(['dS_',num2str(stimu,'%02d'),'=Time_series_RAM;'])
    stimulus = [stimulus zeros(1,2) ones(1,stimu) zeros(1,13)];
end
stimulus_RAM = repmat(stimulus,[1/temporal_resolution,1]);
stimulus_T = reshape(stimulus_RAM,[1,numel(stimulus_RAM)]);
BOLD_T = [dS_01,dS_02,dS_04,dS_08];
stimulus_T(end+(numel(BOLD_T)-numel(stimulus_T))) = 0;
%% Gamma Fitting (FSL)
Signal_Gamma = dS_01;
EX = sum(Signal_Gamma(:))/numel(Signal_Gamma);
DX = sum(Signal_Gamma(:).^2)/numel(Signal_Gamma)-EX^2;
param(1) = EX*EX/DX;    param(2) = EX/DX;  param(3) = max(dS_01(:)); 
lb = param*0.5;
ub = param*100;

options = optimset('LargeScale','on',...  
           'Algorithm', 'trust-region-reflective',...
     'TolFun',10^(-100),'TolX',10^(-10),...
  'MaxFunEvals',1000000,'MaxIter',1000000);

variable{1} = 1:numel(dS_01);
variable{2} = stimulus_T;
%
fun_FSL = @(param,variable)conv(param(3)*gampdf(variable{1},param(1),param(2)),variable{2}*10000);
actual_BOLD = [BOLD_T*10000 zeros(1,(numel(dS_01)-1))];
[param0,~,~,~] = lsqcurvefit(fun_FSL,param,variable,actual_BOLD,lb,ub,options);
RAM = [param0-lb,ub-param0];
if find(RAM<=0.001)
   error('Fitting Bug in Gamma_FSL!'); 
end
HRF_FSL = param0(3)*gampdf(variable{1},param0(1),param0(2));
simu_FSL = conv(HRF_FSL,stimulus_T);
simu_FSL(numel(BOLD_T)+1:end) = [];
R2_FSL = corr(BOLD_T',simu_FSL')^2;
HRF_FSL = HRF_FSL/sum(HRF_FSL(:));
HRF.FSL = HRF_FSL;
HRF.FSL_parameter = param0;

%% Gamma Fitting (SPM)
param(1:3) = param0; param(5:6) = param0(2:3)*1;param(4) = param0(1)*16/6;
lb = param*0.2;
ub = param*10;
fun_SPM = @(param,variable)conv([param(3)*gampdf(variable{1},param(1),param(2))-...
    param(6)*gampdf(variable{1},param(4),param(5))],variable{2}*10000);
actual_BOLD = [BOLD_T*10000 zeros(1,(numel(dS_01)-1))];
[param0,~,~,~] = lsqcurvefit(fun_SPM,param,variable,actual_BOLD,lb,ub,options);
RAM = [param0-lb,ub-param0];
if find(RAM<=0.001)
   error('Fitting Bug in Gamma_SPM!'); 
end
HRF_SPM = param0(3)*gampdf(variable{1},param0(1),param0(2))-param0(6)*gampdf(variable{1},param0(4),param0(5));
simu_SPM = conv(HRF_SPM,stimulus_T);
simu_SPM(numel(BOLD_T)+1:end) = [];
R2_SPM = corr(BOLD_T',simu_SPM')^2;
HRF_SPM = HRF_SPM/sum(HRF_SPM(:));
HRF.SPM = HRF_SPM;
HRF.SPM_parameter = param0;

%% Stimulus Patch
y_min = -0.5; y_max = 2.5;
stimulus_dev = stimulus_T - [stimulus_T(2:end) stimulus_T(1)];
stimulus_up = find(stimulus_dev==1);
stimulus_down = find(stimulus_dev==-1);
patch_x = [stimulus_down;stimulus_down;stimulus_up;stimulus_up]*temporal_resolution;
patch_y = [repmat(y_min,1,numel(stimulus_down));repmat(y_max,1,numel(stimulus_up));
           repmat(y_max,1,numel(stimulus_up));repmat(y_min,1,numel(stimulus_down));];

F = figure;
set(F,'position',[200 200 900 200]*1.5);
%% SubFigure 01 HRF
subplot(2,6,[1,2,7,8]);
plot(variable{1}*temporal_resolution,HRF_FSL,'LineWidth',2,'color',color{1}/255);hold on;
plot(variable{1}*temporal_resolution,HRF_SPM,'LineWidth',2,'color',color{2}/255);hold off;
xlim([0*temporal_resolution numel(HRF_SPM)*temporal_resolution]);
h = legend({'FSL';'SPM'},'Location','best');
set(h,'Box','off');
title('\fontsize{13}Hemodynamic Response Function')
set(gca,'FontSize',8,'FontWeight','bold');
xlabel('Time (s)','Fontsize',10)
set(gca, 'LineWidth',1.5);
box off;
%% SubFigure 02 FSL_Fitting
subplot(2,6,[3,4,5,6]);
plot([1:numel(BOLD_T)]*temporal_resolution,BOLD_T,'k','LineWidth',1);hold on;
plot([1:numel(BOLD_T)]*temporal_resolution,simu_FSL,'color',color{1}/255,'LineWidth',2);
patch(patch_x,patch_y,'r','FaceAlpha',.2,'EdgeColor','none');hold off;
ylim([y_min y_max]);
text(75,2,['R^2 = ',num2str(R2_FSL,'%4.2f')],'FontWeight','Bold')
title('\fontsize{13}Fitting between actual and FSL predicted BOLD responses')
set(gca,'FontSize',8,'FontWeight','bold');
set(gca, 'LineWidth',1.5);
box off;
%% SubFigure 03 SPM_Fitting
subplot(2,6,[9,10,11,12]);
plot([1:numel(BOLD_T)]*temporal_resolution,BOLD_T,'k','LineWidth',1);hold on;
plot([1:numel(BOLD_T)]*temporal_resolution,simu_SPM,'color',color{2}/255,'LineWidth',2);
patch(patch_x,patch_y,'r','FaceAlpha',.2,'EdgeColor','none');hold off;
ylim([y_min y_max]);
text(75,2,['R^2 = ',num2str(R2_SPM,'%4.2f')],'FontWeight','Bold')
title('\fontsize{13}Fitting between actual and SPM predicted BOLD responses')
set(gca,'FontSize',8,'FontWeight','bold');
set(gca, 'LineWidth',1.5);
box off;

end