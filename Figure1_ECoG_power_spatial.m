
clc;clear
codepath = 'F:\ECoG\code_Trange\';
addpath(genpath(codepath));
WholePath = 'F:\ECoG\Activation\';

[~,~,CellData1] = xlsread(fullfile(codepath,'Elec_location.xlsx'));
Loc = cell2mat(CellData1([3:9 11:17],2:3));
[~,~,CellData2] = xlsread(fullfile(codepath,'power_state.xlsx'));

state={'nrem';'rem'};
freqBand={'1delta','2theta','3alpha','4beta','5gamma'};
F =figure;
for st=1:2
    
    for band=1:5
        
        ref = cell2mat(CellData2(3:16,band+1));
        source = cell2mat(CellData2(3:16,band+1+st*5));
        dif0 = (source(:)-ref(:))./abs(ref(:))*100;
        dif(1:7) = dif0(1:7)/1.5+dif0(8:14)/2;
        dif(8:14) = dif0(1:7)/2.0+dif0(8:14)/1.5;
        
        if st==1;lims = [-100 100];end
        if st==2;lims = [-30 30];end
        %lims = [floor(mean(dif)-std(dif)); ceil(mean(dif)+std(dif))];
        %ColorVec = flipud([ones(257,1),(0:256)'/256,(0:256)'/256]);
        %ColorVec(1,:)=[];
%         ColorVec = parula(256);
        gdmap = [(0:127)'/127,(0:127)'/127,ones(128,1)];
        drmap = [ones(128,1),(0:127)'/127,(0:127)'/127];
        ColorVec = [gdmap;flipud(drmap);];
        
        subplot(2,  5,band+(st-1)*5);
        for loop=1:size(Loc,1)
            value = dif(loop);
            pin = round ((value-lims(1)) / ((lims(2)-lims(1))/256));
            if pin>256;pin=256;end
            if pin<1;pin=1;end
            
            plot(Loc(loop,1),Loc(loop,2),'o','MarkerSize',10,'MarkerFaceColor',ColorVec(pin,:),...
                'MarkerEdgeColor',[0 0 0]);
            hold on
            %title(num2str(lims(:)'));
        end
        axis equal
        axis off
        axis tight
        
        
    end
end
saveas(F,fullfile(WholePath,['State_power','.emf']));
1;

