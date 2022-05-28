expInfoObj_path = '/hdd2/RecBCD2/codedev/Analysis/EXP-22-BY4448/therun/prod/expInfoObj.mat';
load(expInfoObj_path)
fluoChanName = expInfoObj.fluoChannelNames{1};

pj = 4;

posName = expInfoObj.positions{pj};
sourceDir = expInfoObj.getPathName('source');
outputDir = expInfoObj.getPathName('output');
mCellsMatFile = expInfoObj.getMCellMatPath(posName);
mCells = Cell.MCell.loadMCells(mCellsMatFile);
mFrames = Cell.MCell.mFrames;
fluoIndices = expInfoObj.getIndices(posName, fluoChanName);
phaseChanName = expInfoObj.getChannelNames('phase');
phaseRange = expInfoObj.getRange(posName, phaseChanName);
Ts = expInfoObj.getAllT(posName, fluoChanName);

pathFluo = fullfile(sourceDir,posName,fluoChanName);
listFluo = dir(fullfile(pathFluo,'*.tif*'));
listFluo = {listFluo.name};
pathSeg = fullfile(outputDir,posName,'SegmentedPhase');
listSeg = dir(fullfile(pathSeg,'*.tif*'));
%listSeg = {listSeg(fluoIndices).name};
listSeg = {listSeg(:).name};
    
% cell_indexes_venus = select_fluo_cells(mCells,2,2,'stdev',0); %change this function to accept expInfoObj?
% cell_indexes_dead = find_dying_cells(mCells);
cell_indexes_dead = [mCells([mCells.lifeTime] < 7).id]
%parfor vi = 1 : length(expInfoObj.imRange{pj, find(strcmp(expInfoObj.getChannelNames,fluoChanName))})

for vi = 1 : length(expInfoObj.imRange{pj, 1})
%     parfor vi = 1 : length(expInfoObj.imRange{pj, 1})

    seg  = imread(fullfile(pathSeg, listSeg{vi}));
    this_frame = vi;

    seg_selected_cherry = zeros(size(seg));
    seg_selected_venus  = zeros(size(seg));
    seg_selected_dead   = zeros(size(seg));

    all_blobs = [];
    for ci = 1 : length(mCells)
        tt = mCells(ci).birthFrame : mCells(ci).lastFrame;
        if ismember(vi,tt)
            thisblob = mCells(ci).blobLabels(tt==vi);
            if thisblob == 0 %sometimes cells can disappear
                continue,
            end
            all_blobs = [all_blobs; [thisblob mCells(ci).id]];
        end
    end

    for blob = 1:size(all_blobs,1)
        this_cell = all_blobs(blob,2); %what if there are 2 cells in one blob?
        %select cells with venus intensity
        if ismember(this_cell,cell_indexes_dead) 
            seg_selected_dead(seg==all_blobs(blob,1)) = 1;
        end
    end

    %find out which cells belong to which blob in the seg image. 
    rgb = cat(3,repmat((seg>0).*0.5,[1,1,3]));
    rgb(:,:,1) = rgb(:,:,1) + (seg_selected_cherry>0).*0.5;
    rgb(:,:,2) = rgb(:,:,2) + (seg_selected_venus>0).*0.5;
    
    sel_cells = imerode(seg_selected_dead,strel('disk',4));
    
    rgb(:,:,1) = rgb(:,:,1)-sel_cells;
    rgb(:,:,2) = rgb(:,:,2)-sel_cells;
    
    imwrite(rgb, fullfile('stack',[num2str(vi) '.tiff']))

end

function cell_inds = select_fluo_cells(mCells, fluo_chan_index, min_fluo_time, algorithm, plot_it)
    %  selecting cells that contain fluorescnce in a one channel given that
    %  two populations of cells are mixed in one mCells structure. Selection
    %  is done by fitting 2-gaussian model onto fluorescent data. Requires
    %  fitgmdist function from statistic toolbox
    %
    %
    %  input: 
    %       mCells - Cell.MCell object
    %   
    %       fluo_chan_index - index of fluorescent channel to analyse as in
    %           mChells.fluoIntensities field
    %
    %       min_fluo_time - minimum number of frames at which cell
    %           displayed fluorescence signal above threshold
    %       
    %        plot_it - bool, if distributions should be plotted.
    %
    %    output:
    %       indexes of cells with fluoresncence above threshold as in
    %       mCells(i).id

    assert(~isempty([mCells.fluoIntensities]), 'fluo intensities are empty')
    
    assert(any(strcmp({'gaussians','stdev'},algorithm)), 'wrong algorithm specified, can be "gaussians", or "stdev"')
    
    if nargin<5
        plot_it = 0;
    end

    dout = cell(length(mCells),1);

    for ni = 1:length(mCells)
        t = mCells(ni).birthFrame:mCells(ni).lastFrame;
        f = mCells(ni).fluoIntensities(fluo_chan_index,:);
        good_ind = ~isnan(f);
        dout{ni} = [t(good_ind)' f(good_ind)' repmat(ni,[sum(good_ind) 1])];
    end

    dout = vertcat(vertcat(dout{:}));
    t = dout(:,1);
    f = dout(:,2);
    
    uni_t = unique(t);
    thresholds = zeros(1,length(uni_t));
    selected_cells = cell(length(uni_t),1);

    switch algorithm
        %maybe mix those approaches in case gaussian fails?
        
        case 'gaussians'
        %1st approach - use gaussian mixture model - needts ditgmidst function
        
        for ti = 1:length(uni_t)
            this_t = uni_t(ti);
            this_f = f(t==this_t);
            this_f = rmoutliers(this_f);
            GMModel = fitgmdist(this_f,2,'Replicates',9);
            mus = sort(GMModel.mu);
            xspace = floor(mus(1)):ceil(mus(2));
            yspace = GMModel.pdf(xspace');
            [~, min_ind] = min(yspace);
            thresholds(ti) = xspace(min_ind);
            selected_cells{ti} = dout(t==uni_t(ti) & f > thresholds(ti), 3); %cell index is 3
        end
    
        case 'stdev'
        %2nd approach - split data in 2 and minimize standard deviation
        
        for ti = 1:length(uni_t)
            this_t = uni_t(ti);
            this_f = f(t==this_t);

            this_f = rmoutliers(this_f);
            this_f_sorted = sort(this_f);        
            min_val = inf;
            min_ind = nan;
            %this loop can be optimized
            for nf = 2:length(this_f_sorted)-1
                a = std(this_f_sorted(1:nf));
                b = std(this_f_sorted(nf+1:end));
                if a+b < min_val
                    min_val = a+b;
                    min_ind = nf+1;
                end
            end
            thresholds(ti) = mean(this_f_sorted(min_ind-1):this_f_sorted(min_ind+1));
            selected_cells{ti} = dout(t==uni_t(ti) & f > thresholds(ti), 3); %cell index is 3
        end
    end
    
    preselected_cells = vertcat(selected_cells{:});
    cells_above_threshold = unique(preselected_cells);
    
    cell_inds = [];
    for cindex = cells_above_threshold'
        if sum(preselected_cells==cindex) > min_fluo_time
            cell_inds = [cell_inds mCells(cindex).id];
        end
    end
    
    if plot_it
       figure
       hold on
       plot(t+randn([length(t) 1])./2,f,'.k','MarkerSize',2), 
       xlim([0 max(t)])
       plot(uni_t, thresholds,'r','LineWidth',2)
       title('Fluorescnce signal and threshold')
       xlabel('Frame number')
       ylabel('Fluorescence intensity (au)')
    end

end


function [cell_indexes] = find_dying_cells(mCells)
% Selection of cells that are dying. It does not work very well yet,
% but seems to have low number of false positives, but still misses
% many cells that stop growing

%try incorporate the concept of good detection frames here?

cell_indexes = {};

size_diff_thresh = -500; %big drop in size
med_lifetime = median([mCells.lifeTime]);
med_size = median(cellfun(@(x) x(end), {mCells.areas}));
% figure

for ci = 1:length(mCells)
    
    this_cell = mCells(ci);
    if ~isempty(this_cell.descendants)
        continue
    end
    
%     this_cell.id
    
    if this_cell.lifeTime < med_lifetime
        continue
    end
    
    a = this_cell.areas;
    a = medfilt1(a,5);
    
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
end

