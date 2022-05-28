
codePath = '/hdd2/RecBCD2/codedev/ImAnalysis/';
addpath(codePath);
ImAnalysis_setup(); 

expInfoObj_base= '/hdd2/RecBCD2/codedev/Analysis/EXPHANDLE/therun/prod/expInfoObj.mat';


experiments = {'EXP-22-BY4442','EXP-22-BY4448'};

for expnum = 1:length(experiments)
    
disp(' ')   
disp(experiments{expnum})

expInfoObj_path = strrep(expInfoObj_base, 'EXPHANDLE', experiments{expnum});

load(expInfoObj_path)

numpos = length(expInfoObj.positions);

channel_name_mCherry = expInfoObj.fluoChannelNames{contains(expInfoObj.fluoChannelNames,'594')};
channel_name_SSB = expInfoObj.fluoChannelNames{contains(expInfoObj.fluoChannelNames,'514') | contains(expInfoObj.fluoChannelNames,'venus')};
SSB_channel_index = find(contains(expInfoObj.fluoChannelNames,'514') | contains(expInfoObj.fluoChannelNames,'venus'));

save_dir_base = fullfile('/hdd2/RecBCD2/codedev/Analysis/output', experiments{expnum});
mkdir(save_dir_base)

out = {};

for pj = 1 : numpos
%for pj = 7
    
    %
    %initialize experiment parameters
    %
    
    posName = expInfoObj.positions{pj};
    sourceDir = expInfoObj.getPathName('source');
    outputDir = expInfoObj.getPathName('output');
    mCellsMatFile = expInfoObj.getMCellMatPath(posName);
    mCells = Cell.MCell.loadMCells(mCellsMatFile);
    mFrames = Cell.MCell.mFrames;
%     fluoIndices = expInfoObj.getIndices(posName, fluoChanName);
%     phaseChanName = expInfoObj.getChannelNames('phase');
%     phaseRange = expInfoObj.getRange(posName, phaseChanName);
%     Ts = expInfoObj.getAllT(posName, fluoChanName);
    
%     pathFluo = fullfile(sourceDir,posName,fluoChanName);
%     listFluo = dir(fullfile(pathFluo,'*.tif*'));
%     listFluo = {listFluo.name};

    pathSeg = fullfile(outputDir,posName,'SegmentedPhase');
    listSeg = dir(fullfile(pathSeg,'*.tif*'));
%     listSeg = {listSeg(fluoIndices).name};
    listSeg = {listSeg(:).name};
    
    %
    %classification part
    %
    
    cell_indexes_venus = select_fluo_cells(mCells, SSB_channel_index,2,'stdev',1); %change this function to accept expInfoObj?
    % get gcf, save, close
    saveas(gcf,fullfile(save_dir_base, [posName '_venusThreshold.png']))
    close(gcf)

    cell_indexes_cherry = find_cells_with_spots(expInfoObj,posName,channel_name_mCherry); %improved version
    % get gcf, save, close
    saveas(gcf,fullfile(save_dir_base, [posName '_mCherryDetection.png']))
    close(gcf)
    
    dead_indexes_cells = find_dying_cells(mCells);
    
    out{pj} = {[mCells.id], cell_indexes_venus, cell_indexes_cherry, dead_indexes_cells, [mCells.birthFrame]};
    
    %
    %plotting part
    %
    
    plotit = 0;
    if pj == 1
        plotit = 1;
        stack_save_dir = fullfile(save_dir_base,['stack_' posName]);
        mkdir(stack_save_dir)
    end
    
    if ~plotit
        continue
    end
    
	%parfor vi = 1 : length(expInfoObj.imRange{pj, find(strcmp(expInfoObj.getChannelNames,fluoChanName))})
    parfor (vi = 1 : length(expInfoObj.imRange{pj, 1}), 5)
%     parfor vi = 1 : length(expInfoObj.imRange{pj, 1})
        
        %fluo = imread(fullfile(pathFluo, listFluo{vi}));
        seg  = imread(fullfile(pathSeg, listSeg{vi}));
        %this_frame = fluoIndices(vi);
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
            
%             if mCells(this_cell).isBadCell ~= 0
%                 continue
%             end
            
            %select cells with spots
            if ismember(this_cell,cell_indexes_cherry) 
                seg_selected_cherry(seg==all_blobs(blob,1)) = 1;
            end
            
             %select cells with venus intensity
            if ismember(this_cell,cell_indexes_venus) 
                seg_selected_venus(seg==all_blobs(blob,1)) = 1;
            end
            
            %select dead cells
            if ismember(this_cell,dead_indexes_cells) 
                seg_selected_dead(seg==all_blobs(blob,1)) = 1;
            end
                        
        end
                
            
        rgb = cat(3,repmat((seg>0).*0.5,[1,1,3]));
        rgb(:,:,1) = rgb(:,:,1) + (seg_selected_cherry>0).*0.5;
        rgb(:,:,2) = rgb(:,:,2) + (seg_selected_venus>0).*0.5;

        sel_cells = imerode(seg_selected_dead,strel('disk',4));

        rgb(:,:,1) = rgb(:,:,1)-sel_cells;
        rgb(:,:,2) = rgb(:,:,2)-sel_cells;

        %write image
        imwrite(rgb, fullfile(stack_save_dir, [num2str(vi) '.tiff']))
            
%         T = Ts(:,:,phaseRange(vi));
% 
%         transFluoIm = imwarp(fluo,affine2d(T),'outputView',imref2d(size(seg)));
%         transFluoIm = transFluoIm - mean(transFluoIm(seg==0));
% 
%         %imshowpair(edge(seg>0) + edge(seg_selected>0), imadjust(transFluoIm))
%         %imshowpair(edge(seg>0),seg_selected)
%         rgb = cat(3,repmat(imadjust(transFluoIm),[1,1,3]));
%         rgb(:,:,2) = min(rgb(:,:,2) + uint16(edge(seg_selected_venus>0)).*((2^16)-1), ((2^16)-1));
%         rgb(:,:,1) = min(rgb(:,:,1) + uint16(edge(seg_selected_cherry>0)).*((2^16)-1), ((2^16)-1));
%         imshow(rgb,[])
%         %imwrite(uint16(rgb), fullfile('stack',[num2str(vi) '.png']))

    end
    
end

T = plot_stats(out);
writetable(T,fullfile('/crex/proj/uppstore2018129/elflab/Projects/CRISPR_conjugation/codedev/data',[experiments{expnum},'.csv']))

% histogram([mCells(intersect(cell_indexes_venus,cell_indexes_cherry)).birthFrame])
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%        functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function cell_inds = select_spots_cells(mCells,channel_name,min_spot_time)
    
    %documentation. Its important here to keep the min_spot_time as
    % '1' because mCherry is imaged so infrequently, cells dont have much
    % time to accumulate frames with spots.
    %OLD VERSION!
    
    assert(~isempty([mCells.particles]), 'fluo intensities are empty')

    cell_inds = [];
    for ci = 1:length(mCells)
        spot_indexes = find(strcmp({mCells(ci).particles.channelName},channel_name));
        spot_struct = mCells(ci).particles(spot_indexes);
        frames_with_foci = ~cellfun(@isempty,{spot_struct.id});
        if length(unique([spot_struct(frames_with_foci).cellDetectionFrames])) > min_spot_time
            cell_inds = [cell_inds mCells(ci).id];
        end
    end
       
end

function cell_indexes = find_cells_with_spots(expInfoObj,pos_name,fluoChanName)
% Function to detect ParB foci for CRISPR conjugation project. The spots are
% first detected usign radial symmetry algorightm and then are filtered by
% applying 2d gaussian fit and discarding foci with poor fit, or with broad
% fit. 
%
%   cell_indexes = findCellsWithSpots(expInfoObj,'Pos0','fluo594')
%
%   Input:
%       expInfoObj - struct, expInfoObj structure from ImAnalysis pipeline
%       pos_index  - string, name of the position as in expInfoObj
%       chan_index - string, name of fluorescenc channel, as in expInfoObj
%
%   Output:
%       cell_indexes - indexes of cells with spots as in mCell.id

pos_index = find(strcmp(expInfoObj.getPositionList,pos_name));
if isempty(pos_index)
    error(['Wrong position: ' pos_name])
end

chan_index = find(strcmp(expInfoObj.channels,fluoChanName));
if isempty(chan_index)
    error(['Wrong channel: ' fluoChanName])
end

posName = expInfoObj.positions{pos_index};
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
listSeg = {listSeg(fluoIndices).name};

%radial symmetry parameters
radii = [2 3];
alpha =  3;
beta = 0;
gaussKernRadFactor = 1;
radius = 5;
thresh = 1.5;
subPixel = 0;

% dotDetectionAlgorithm = radialSymmetry
% radii = 2 3 4
% alpha = 3
% beta = 0
% radThreshold = 5
% radSymThreshold = 1.6
% gaussKernRadFactor = 1.5
% cellInclusionMargin = 2

%for gaussian fitting initial guesses: this is radius of searched spots
gauss_kren_size = 5;
fttype = fittype( @(a,b1,s1,b2,s2,x,y) a*exp(-(((x-b1).^2)/(2*s1^2) + ((y-b2).^2)/(2*s2^2))) ,'numindep', 2);
[X,Y] = meshgrid(1:1+2*gauss_kren_size);

cell_indexes = {};

for vi = 1:length(expInfoObj.imRange{pos_index, chan_index})
    
    fluo = imread(fullfile(pathFluo, listFluo{vi}));
    seg  = imread(fullfile(pathSeg, listSeg{vi}));
    this_frame = fluoIndices(vi);
        
    T = Ts(:,:,phaseRange(vi));
    
    transFluoIm = imwarp(fluo,affine2d(T),'outputView',imref2d(size(seg)));
        
    S = BuildTrajectories.fastradial(transFluoIm, radii, alpha, beta, gaussKernRadFactor);
    [yy,xx] = BuildTrajectories.nonmaxsuppts(S, radius, thresh, subPixel); %gives position of the spots
        
    %looking for spots that have good fit based on 'r^2' and are quite
    %tight - that is with small std.dev.
    keep_spots = {};
    sigmasX = {};
    sigmasY = {};
        
    parfor sp = 1:length(xx)
        
        this_crop = transFluoIm(yy(sp)-gauss_kren_size : yy(sp)+gauss_kren_size ,...
            xx(sp)-gauss_kren_size : xx(sp)+gauss_kren_size);
        
        %fit gaussian to spot
        [ft,gof] = fit([X(:) Y(:)], double(this_crop(:)), fttype,...
            'StartPoint', [max(this_crop(:)),...
            gauss_kren_size,...
            gauss_kren_size,...
            gauss_kren_size,...
            gauss_kren_size]);
        
        keep_spots{sp} = gof.rsquare;
        sigmasX{sp} = ft.s1;
        sigmasY{sp} = ft.s2;
    end
    
    keep_spots = [keep_spots{:}];
    sigmasX = [sigmasX{:}];
    sigmasY = [sigmasY{:}];
    sel_condition = keep_spots > 0.5 & sigmasX < 1+2*gauss_kren_size & sigmasY < 1+2*gauss_kren_size;
    
    %swith plotting on/off
    plotit = 0;
    if vi == 5
        plotit = 1;
    end
    
    if plotit
        figure
        imshow(transFluoIm,prctile(double(transFluoIm(:)),[0, 99]))
        hold on
        plot(xx,yy,'ro')
    end
    
    xx = xx(sel_condition);
    yy = yy(sel_condition);
    
    if plotit
        plot(xx,yy,'go')
        drawnow
    end

    %disp('hey')

    %connect blobs to cells - save blob id, cell id, and number of detected spots
    all_blobs = [];
    for this_cell = mCells'
        
        tt = this_cell.birthFrame : this_cell.lastFrame;
        
        if ismember(this_frame,tt)
          
            thisblob = this_cell.blobLabels(tt==this_frame);
            
            if thisblob == 0 %sometimes cells can disappear
                continue,
            end
            
            all_blobs = [all_blobs; [thisblob this_cell.id 0]];
            
        end
        
    end
    
    T = array2table(all_blobs, 'VariableNames', {'blobId','cellId','spotsNum'});
    
    %put spots inside cells - use region props to speed up
    spot_mat = zeros(size(seg));
   
    for eachspot = 1:length(yy)
   
        spot_mat(yy(eachspot),xx(eachspot)) = eachspot;
   
    end
    
    seg_fat = imdilate(seg, strel('disk',2));
    propspot = regionprops(seg_fat,spot_mat,'PixelValues');
    
    %cell_spot_list = {};
%     testim = zeros(size(seg)); %for plotting cells with spots
    
    for bj = 1:length(propspot)
        
        tmp = unique(propspot(bj).PixelValues);
        
        if sum(tmp>0) > 0
        
            %cell_spot_list{bj} = [xx(tmp(tmp>0)) yy(tmp(tmp>0))]; %xy position of spots
            T.spotsNum(T.blobId == bj) = sum(tmp>0);
%             testim = testmat + (seg == bj);
      
        end
        
    end
    
%     imshowpair(transFluoIm.*20,edge(testim))
%     disp('breakpoint')
    
    cell_indexes{vi} =  T.cellId(T.spotsNum>0);
    
end

cell_indexes = unique(vertcat(cell_indexes{:}))';

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
            thresholds(ti) = this_f_sorted(min_ind);
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
end


function T = plot_stats(out)
% a quick wrapper on the out cell data structure to print some stats about
% the experiment

% cell_indexes_venus = select_fluo_cells(mCells,2,2,'stdev',0);
% cell_indexes_cherry = select_spots_cells(mCells,'fluor_594_cherry',1);
% dead_indexes_dead = find_dying_cells(mCells);

num_all = 0;
num_recipient = 0;
num_donor = 0;
num_conj = 0;
num_unselected = 0;
num_dead_venus = 0;
num_dead_cherry = 0;
    
res = {};

for ii = 1:length(out)
    %1 all cells
    %2 venus
    %3 cherry
    %4 dead
    %5 birthframe
    if isempty(out{ii})
        continue
    end
    
    res_tmp = cell(length(out{ii}{1}),5);
    
    all_cells = [out{ii}{1}];
    cell_indexes_venus = out{ii}{2};
    cell_indexes_cherry = out{ii}{3};
    dead_indexes_dead = out{ii}{4};
    
    unselected_cells = setxor(all_cells,unique([cell_indexes_venus cell_indexes_cherry]));
    only_cherry_cells = setxor(cell_indexes_cherry , intersect(cell_indexes_cherry, cell_indexes_venus));
    dead_venus = intersect(cell_indexes_venus,dead_indexes_dead);
    dead_cherry = intersect(only_cherry_cells ,dead_indexes_dead);
    
    num_all = num_all + length(all_cells);
    num_recipient = num_recipient + length(cell_indexes_venus);
    num_donor = num_donor + length(cell_indexes_cherry);
    num_conj = num_conj + length(intersect(cell_indexes_venus,cell_indexes_cherry));
    num_unselected = num_unselected + length(unselected_cells);
    num_dead_venus = num_dead_venus + length(dead_venus);
    num_dead_cherry = num_dead_cherry + length(dead_cherry);
    
    for eachcell = all_cells
    	res_tmp{eachcell,1} = ii;
        res_tmp{eachcell,2} = eachcell;
        res_tmp{eachcell,3} = ismember(eachcell,cell_indexes_venus);
        res_tmp{eachcell,4} = ismember(eachcell,cell_indexes_cherry);
        res_tmp{eachcell,5} = ismember(eachcell, dead_indexes_dead);
    end
    
    res{ii} = res_tmp;
    
end
res = vertcat(res{:});

%make into table and save the table
T = cell2table(res, 'VariableNames',{'posname', 'id', 'SSB', 'mCherry', 'stressed'});

fprintf('\n========= Cell stats ===========\n')
fprintf('Num. of cells: %d\n', num_all)
fprintf('Num. of recipient cells (v+/c+-): %d\n', num_recipient)
fprintf('Num. of donor cells (v-/c+): %d\n', num_donor)
fprintf('Num. of conjugants cells (v+/c+): %d\n', num_conj)
fprintf('Num. of unclassified cells (v-/c+): %d\n', num_unselected)
fprintf('Num. of dead/SOS recipients cells: %d\n', num_dead_venus)
fprintf('Num. of dead donors %d\n', num_dead_cherry)
end