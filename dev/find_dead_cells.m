expInfoObj_path = '/hdd2/RecBCD2/codedev/Analysis/EXP-22-BY4448/therun/prod/expInfoObj.mat';
load(expInfoObj_path)

posName = 'Pos11'
mCellsMatFile = expInfoObj.getMCellMatPath(posName);
mCells = Cell.MCell.loadMCells(mCellsMatFile);
    
cell_indexes_dead = find_dying_cells(mCells)


function res = find_dying_cells(mCells)
 %selection of stressed/dying cells. The cells is classified as
 %stressed/dying if it stops growing (when exp1 fit is poor and its size
 %can be better fit by a parabola), or when it is 3 times larger than a
 %typical cell at division.
 %
 % input: mCells structure from ImAnalysis pipeline
 %
 % output: res - vector where each index corresponds to mCells.id field and is
 %           1 if the cell is stressed/dead, 0 otherwise

cell_indexes = {};

med_lifetime = median([mCells.lifeTime]);
med_size = median(cellfun(@(x) x(end), {mCells.areas}));
size_diff_thresh = -med_size/3; %big drop in size

for ci = 1:length(mCells)
    
    this_cell = mCells(ci);
    
    %cells shouldnt have descendants
    if ~isempty(this_cell.descendants)
        continue
    end
        
    %cells should live longer than typical
    if this_cell.lifeTime < med_lifetime
        continue
    end
    
    a = this_cell.areas;
    a = medfilt1(a,5);
    
    %if a cell is two times larger than typical is selected
    if sum(a > med_size*2) > 3
         cell_indexes{ci} = this_cell.id;
         continue
    end

    t = 1:length(a);
    [~, outlierInd] = rmoutliers(a);
    a = a(~outlierInd);
    t = t(~outlierInd);
    
    if length(a) <= 5
        continue
    end
    
    [ft, gof] = fit(t', a', 'exp1');
    
    %cells with poor exp fit are selected
    if gof.rsquare < 0.5
        [ft2, gof2] = fit(t', a', 'poly2');        
        if  ft2.p1 < 0 &&  min(diff(a)) > size_diff_thresh
            cell_indexes{ci} = this_cell.id;
%             nexttile
%             hold on
%             plot(t,a,'ok')
%             plot(t,ft2(t),'b')
%          
        end
    end
    %
end

cell_indexes = horzcat(cell_indexes{:});
res = zeros(length(mCells),1);
for tc = mCells'  
    if ismember(tc.id,cell_indexes)
        res(tc.id) = 1;
    end
end

end
