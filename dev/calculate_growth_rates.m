%calculate growth rates for donors, recipients and conjugated cells in the
%experiments. First 'cellClassesCC.m' scripts needs to by run to get the
%cell classes table! 06-06-2022.

expname = 'EXP-22-BY4448';

load(['/hdd2/RecBCD2/codedev/Analysis/' expname '/therun/prod/expInfoObj.mat'])
celltable = readtable(['/hdd2/RecBCD2/codedev/data/' expname '/cell_data.csv']);
spottable = readtable(['/hdd2/RecBCD2/codedev/data/' expname '/spot_data.csv']);

% cellrange = celltable.id(celltable.posnum==2 & celltable.venus & celltable.cherry &~celltable.dead);

% %try to screen cells based on the spot detection. 
% cellrange = [];
% for posnum = 1%:length(expInfoObj.positions)
%     
%     % store cell id, cell images, ourline, frame indexes,
%     mCellsMatFile = expInfoObj.getMCellMatPath(expInfoObj.positions{posnum});
%     mCells = Cell.MCell.loadMCells(mCellsMatFile);
%     subspottable = spottable(spottable.posnum == posnum,:);
%     subcelltable = celltable(celltable.posnum == posnum,:);
%     
%     for celliterator = 1:length(mCells)
%         tmp = subspottable(subspottable.id == mCells(celliterator).id,:);
%         
%         if isempty(tmp), continue, end
%         
%         if max(unique(tmp.frame)) - min(unique(tmp.frame)) > mCells(celliterator).lifeTime / 2
%             
%             if subcelltable.venus(subcelltable.id==mCells(celliterator).id)
%                 cellrange = [cellrange celliterator];
%             end
%                 
%         end
%         
%     end
%     
% end
% %%

res = {};

for posnum = 1:length(expInfoObj.positions)
    
    % store cell id, cell images, ourline, frame indexes,
    mCellsMatFile = expInfoObj.getMCellMatPath(expInfoObj.positions{posnum});
    mCells = Cell.MCell.loadMCells(mCellsMatFile);
    
    %1-donors, 2-recipients, 3-conjugates
    
    out=zeros(length(mCells),3);
    % posnum = 1;
    subcelltable = celltable(celltable.posnum==posnum,:);
    % figure
    % hold on
    parfor celliterator = 1:length(mCells)
        celli = mCells(celliterator);
        
        if isempty(celli.descendants), continue, end
        
        if celli.isBadCell~=0, continue, end
        
        a = celli.areas;
        t = celli.birthFrame:celli.lastFrame;
        gind = celli.badSegmentations==0;
        a = a(gind);
        t = t(gind);
        
        if length(t) < 10, continue, end
        
        [ft, gof] = fit(t',a','exp1');
        
        if gof.rsquare < 0.5, continue, end
            
        db_tm = log(2)/ft.b;
        
        %ft.b, since t=1 minute then the growth rate is 1xft.b; ft.b. Doubling
        %time is ln(2)/grRate
        
        cell_in_table = find(subcelltable.id==celli.id);
        
        if subcelltable.dead(cell_in_table), continue, end
        
        if subcelltable.venus(cell_in_table) & subcelltable.cherry(cell_in_table)
            out(celliterator,:) = [0 0 db_tm];
        elseif ~subcelltable.venus(cell_in_table)
            out(celliterator,:) = [db_tm 0 0];
        elseif subcelltable.venus(cell_in_table)
            out(celliterator,:) = [0 db_tm 0];
        end
    end
    
    res{posnum} = out;
end
%
res_cat = vertcat(res{:});


figure
subplot(1,2,1)
hold on
histogram( res_cat(res_cat(:,3)>0,3),20:1:100,'normalization','pdf','EdgeColor','k','DisplayStyle','stairs','LineWidth',2)
histogram( res_cat(res_cat(:,1)>0,1),20:1:100,'normalization','pdf','EdgeColor','m','DisplayStyle','stairs','LineWidth',2)
histogram( res_cat(res_cat(:,2)>0,2),20:1:100,'normalization','pdf','EdgeColor','g','DisplayStyle','stairs','LineWidth',2)
xlabel('doubling time (min)')
ylabel('pdf')

legend({'trans','donor','recipinet'})

subplot(1,2,2)
boxplotdata =  [res_cat(res_cat(:,3)>0,3);...
                res_cat(res_cat(:,1)>0,1);...
                res_cat(res_cat(:,2)>0,2)];
            
boxplotlabels = [repmat('trans',[length(res_cat(res_cat(:,3)>0,3)),1]);...
                 repmat('donor',[length(res_cat(res_cat(:,1)>0,1)),1]);...
                 repmat('recip',[length(res_cat(res_cat(:,2)>0,2)),1])];


boxplot(boxplotdata,boxplotlabels,'Notch','on')
ylim([20 100])
ylabel('doubling time (min)')
