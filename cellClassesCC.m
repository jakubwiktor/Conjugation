%uppmax code
% codePath = '/home/bobdylan/elflab/Projects/CRISPR_conjugation/codedev/ImAnalysis/';
% addpath(codePath);
% ImAnalysis_setup(); 
% expInfoObj_base= '/home/bobdylan/elflab/Projects/CRISPR_conjugation/Analysis2/EXPHANDLE/therun/prod/expInfoObj.mat';
% experiments = {'EXP-22-BJ7093','EXP-22-BY4412','EXP-22-BY4413','EXP-22-BY4414','EXP-22-BY4415','EXP-22-BY4419','EXP-22-BY4436','EXP-22-BY4440','EXP-22-BY4441','EXP-22-BY4442','EXP-22-BY4444','EXP-22-BY4447','EXP-22-BY4448'};

%local code
codePath = '/hdd2/RecBCD2/codedev/ImAnalysis/';
addpath(codePath);
ImAnalysis_setup(); 
expInfoObj_base= '/hdd2/RecBCD2/codedev/Analysis/EXPHANDLE/therun/prod/expInfoObj.mat';
experiments = {'EXP-22-BY4442','EXP-22-BY4448'};

% for expnum = 1:length(experiments)
for expnum = 2
    
expInfoObj_path = strrep(expInfoObj_base, 'EXPHANDLE', experiments{expnum});

load(expInfoObj_path)

numpos = length(expInfoObj.positions);

channel_name_mCherry = expInfoObj.fluoChannelNames{contains(expInfoObj.fluoChannelNames,'594')};
channel_name_SSB = expInfoObj.fluoChannelNames{contains(expInfoObj.fluoChannelNames,'514') | contains(expInfoObj.fluoChannelNames,'venus')};
SSB_channel_index = find(contains(expInfoObj.fluoChannelNames,'514') | contains(expInfoObj.fluoChannelNames,'venus'));

% save_dir_base = fullfile('/crex/proj/uppstore2018129/elflab/Projects/CRISPR_conjugation/codedev/data', experiments{expnum});
save_dir_base = fullfile('/hdd2/RecBCD2/codedev/Analysis/output', experiments{expnum});
mkdir(save_dir_base)

out = {};
res = {};
% for pj = 1 : numpos
for pj = 2 %: numpos
    
    %
    %initialize experiment parameters
    %
    
    posName = expInfoObj.positions{pj};
    sourceDir = expInfoObj.getPathName('source');
    outputDir = expInfoObj.getPathName('output');
    mCellsMatFile = expInfoObj.getMCellMatPath(posName);
    mCells = Cell.MCell.loadMCells(mCellsMatFile);
    mFrames = Cell.MCell.mFrames;
   
    %
    %classification part
    %
    
    figure
    cell_indexes_venus = select_fluo_cells(mCells,SSB_channel_index,'stdev',1); %change this function to accept expInfoObj?
    
    %save figure
    saveas(gcf,fullfile(save_dir_base, [posName '_venusThreshold.png']))
    close(gcf)
    
    figure
    [cell_indexes_cherry, Tspot] = find_cells_with_spots(expInfoObj,posName,channel_name_mCherry); %improved version
    
    %save figure
    saveas(gcf,fullfile(save_dir_base, [posName '_mCherryDetection.png']))
    close(gcf)
        
    dead_indexes_cells = find_dying_cells(mCells);
    
    res{pj} = [repmat(pj,[length(mCells),1]) [mCells.id]' cell_indexes_venus cell_indexes_cherry dead_indexes_cells];
    
    %
    %plotting part
    %
    if pj==1
        stack_save_dir = fullfile(save_dir_base,['stack_' posName]);
        mkdir(stack_save_dir)
        make_stack(expInfoObj, posName, res{pj}, stack_save_dir)
    end
    
end

%spell out the results
disp('')
disp('Results for:')   
disp(experiments{expnum})
T = array2table(vertcat(res{:}), 'VariableNames',{'posnum', 'id', 'venus', 'cherry', 'dead'});
fprintf('\n========= Cell stats ===========\n')
fprintf('Num. of cells: %d\n', height(T))
fprintf('Num. of recipient cells (v+/c-): %d\n', sum(T.venus))
fprintf('Num. of donor cells (v-/c+): %d\n', sum(~T.venus & T.cherry))
fprintf('Num. of conjugants cells (v+/c+): %d\n', sum(T.venus & T.cherry))
fprintf('Num. of unclassified cells (v-/c+): %d\n', sum(~T.venus & ~T.cherry))
fprintf('Num. of dead/SOS recipients cells: %d\n', sum(T.dead & T.venus))
fprintf('Num. of dead donors %d\n', sum(T.dead & (~T.venus & T.cherry)))

writetable(Tspot,fullfile(save_dir_base,'spot_data.csv'))
writetable(T,fullfile(save_dir_base,'cell_data.csv'))

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%        functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function make_stack(expInfoObj, posName, cell_classes, stack_save_dir)
% Construct RGB image for each segmentation image and color cells depending
% on the selection criteria. Green - 'venus high' cells, red - 'mcherry'
% cells, and yellow - 'transconjugated' cells. Cells with just an outline
% and dark inside are classified as dead/stressed. 
%
% inupt: expInfoObj - expInfo structure from ImAnalysis pipeline
%        posName - string, name of the position as in expInfoObj.positions
%        cell_classes - n-by-5 matrix where the columns are: 1-position
%        number, 2-cellId (as in mCells.id), 3-to-5 - bool cell clases
%        specifying if given cell is classified as venus, cherry, or dead.
%
% output: none, the function saves the images in spacified folder.

outputDir = expInfoObj.getPathName('output');
mCellsMatFile = expInfoObj.getMCellMatPath(posName);
mCells = Cell.MCell.loadMCells(mCellsMatFile);

pj = find(ismember(expInfoObj.positions,posName));
pathSeg = fullfile(outputDir,posName,'SegmentedPhase');
listSeg = dir(fullfile(pathSeg,'*.tif*'));
listSeg = {listSeg(:).name};

T = array2table(cell_classes, 'VariableNames',{'posnum', 'id', 'venus', 'cherry', 'dead'});

parfor (vi = 1 : length(expInfoObj.imRange{pj, 1}), 5)
% for vi = 1 : length(expInfoObj.imRange{pj, 1})

    seg  = imread(fullfile(pathSeg, listSeg{vi}));
   
    %vectorize cell selection - 20x faster than loops
    props = regionprops(seg,'PixelList');
    
    %calculate links between blobs and cells and classes
    all_blobs = zeros(length(mCells)-1,5);
    for ci = 1 : length(mCells)
        tt = mCells(ci).birthFrame : mCells(ci).lastFrame;
        if ismember(vi,tt)
            thisblob = mCells(ci).blobLabels(tt==vi);
            
            if thisblob == 0 %sometimes cells can disappear
                continue,
            end
            %cell id, blobNumber, cell clasees
            all_blobs(ci,:) = [mCells(ci).id thisblob T.venus(ci) T.cherry(ci) T.dead(ci)];
            
        end
    end
    all_blobs(all(all_blobs == 0,2),:) = []; %remove empty rows
    
    seg_selected_venus  = zeros(size(seg));
    tmp = vertcat(props(all_blobs(all_blobs(:,3)==1,2)).PixelList);
    seg_selected_venus(sub2ind(size(seg),tmp(:,2),tmp(:,1))) = 1;
    
    seg_selected_cherry = zeros(size(seg));
    tmp = vertcat(props(all_blobs(all_blobs(:,4)==1,2)).PixelList);
    seg_selected_cherry(sub2ind(size(seg),tmp(:,2),tmp(:,1))) = 1;
    
    seg_selected_dead   = zeros(size(seg));
    tmp = vertcat(props(all_blobs(all_blobs(:,5)==1,2)).PixelList);
    seg_selected_dead(sub2ind(size(seg),tmp(:,2),tmp(:,1))) = 1;
   
    rgb = cat(3,repmat((seg>0).*0.5,[1,1,3]));
    rgb(:,:,1) = rgb(:,:,1) + (seg_selected_cherry>0).*0.5;
    rgb(:,:,2) = rgb(:,:,2) + (seg_selected_venus>0).*0.5;
    
    sel_cells = imerode(seg_selected_dead,strel('disk',4));
    
    rgb(:,:,1) = rgb(:,:,1)-sel_cells;
    rgb(:,:,2) = rgb(:,:,2)-sel_cells;
    
    %imshow(rgb)
    %drawnow

    imwrite(rgb, fullfile(stack_save_dir, [num2str(vi) '.tiff']))
    
end
end


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

function [res, Tspot] = find_cells_with_spots(expInfoObj,pos_name,fluoChanName)
% Function to detect ParB foci for CRISPR conjugation project. The spots are
% first detected usign radial symmetry algorightm and then are filtered by
% applying 2d gaussian fit and discarding foci with poor fit, or with broad
% fit. 
%
%   res = findCellsWithSpots(expInfoObj,'Pos0','fluo594')
%
%   Input:
%       expInfoObj - struct, expInfoObj structure from ImAnalysis pipeline
%       pos_index  - string, name of the position as in expInfoObj
%       chan_index - string, name of fluorescenc channel, as in expInfoObj
%
%   Output:
%       res - vector where each index corresponds to mCells.id field and is
%           1 if the cell had at least a single detected spot during its
%           lifetime, and 0 when not.
%       Tspot - table with frame number, cellId, and x and y position of
%          detected spots

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
    
    %use regionprops to put spots inside cells fast
    
    %this makes 2d image of with pixel value related to the spot index
    spot_mat = zeros(size(seg));
    for eachspot = 1:length(yy)
        spot_mat(yy(eachspot),xx(eachspot)) = eachspot;
    end
    
    seg_fat = imdilate(seg, strel('disk',2));
    propspot = regionprops(seg_fat,spot_mat,'PixelValues');
    
    %cell_spot_list = {};
	%testim = zeros(size(seg)); %for plotting cells with spots
    
    spot_pos_mat = []; %also save all spot detections together with the frame number and cellId
    for bj = 1:length(propspot)
        tmp = unique(propspot(bj).PixelValues);
        if sum(tmp>0) > 0
            T.spotsNum(T.blobId == bj) = sum(tmp>0);
            if ~isempty(T.spotsNum(T.blobId == bj))
                cellid = T.cellId(find(T.blobId == bj));
                
                %store positions of spots
                for this_cellid = cellid'
                    spot_pos_mat = [spot_pos_mat; [repmat(this_frame,[sum(tmp>0) 1]),...
                                                   repmat(this_cellid,[sum(tmp>0) 1]),...
                                                   xx(tmp(2:end)),...
                                                   yy(tmp(2:end))]];
                end
            end
%             testim = testmat + (seg == bj);
        end
    end
    
%     imshowpair(transFluoIm.*20,edge(testim))
%     disp('breakpoint')

    spot_positions{vi} = spot_pos_mat;
    cell_indexes{vi} =  T.cellId(T.spotsNum>0);
    
end

%reconstruct the output matrix to match n-by-2
cell_indexes = unique(vertcat(cell_indexes{:}));

%prepare the output matrix
res = zeros(length(mCells),1);
for tc = mCells'  
    if ismember(tc.id,cell_indexes)
        res(tc.id) = 1;
    end
end

%prepare table of detected spots
spot_positions = vertcat(spot_positions{:});
Tspot = array2table(spot_positions, 'VariableNames', {'frame','cellId','x','y'});
Tspot = sortrows(Tspot,["frame","cellId"]);
end

function [res] = select_fluo_cells(mCells, fluo_chan_index, algorithm, plot_it)
    %  selecting cells that contain fluorescnce in a one channel given that
    %  two populations of cells are mixed in one mCells structure.
    %  3 different approaches can be selected - minimizing standar
    %  deviation between 2 populations, fitting 2-gaussian model, or by
    %  clustering the data (the last approach consolidates fluorescence and
    %  time and is suitable for conditions when one population dissapears
    %  during the experiments - e.g. when cells are pushed out of the
    %  growth chamber). Cells is fluorescent when at least half of its
    %  frames are in high fluorescnce population.
    %
    %  input: 
    %       mCells - mCells object from ImAnalysis pipeline
    %   
    %       fluo_chan_index - index of fluorescent channel to analyse as in
    %           mChells.fluoIntensities field
    %       
    %        plot_it - bool, if distributions should be plotted.
    %        algorithm - 'stdev,'gaussian','clustering'
    %
    %    output:
    %       res = vector where indexes correspond to mCells.id field and
    %       the value is 1 if cell is classified as high-fluorescent, and
    %       0 otherwise

    assert(~isempty([mCells.fluoIntensities]), 'fluo intensities are empty')
    
    assert(any(strcmp({'gaussians','stdev','clustering'},algorithm)), 'wrong algorithm specified, can be "gaussians", or "stdev"')
    
    if nargin<4
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
    cell_classes = zeros(size(dout,1),1);
    
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
            cell_classes(t==uni_t(ti) & f > thresholds(ti)) = 1;
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
            cell_classes(t==uni_t(ti) & f > thresholds(ti)) = 1;
        end

        case 'clustering'
        %3rd approach - cluster data - works if one population
        %dissapears during the exp., when cells are not mixed while
        %loading
        fitdata = [];
        for tc = mCells'
            f = tc.fluoIntensities(fluo_chan_index,:);
            gi = ~isnan(f);
            f = f(gi);
            t = tc.birthFrame:tc.lastFrame;
            t = t(gi);
            [~, idxs] = rmoutliers(f);
            fitdata = [fitdata; [t(~idxs)' f(~idxs)']];
        end

        %correct for bleaching
        ft = fit(fitdata(:,1),fitdata(:,2),'exp2');
        x = unique(dout(:,1));
        fitvals = ft(x);
        corrvals = fitvals./fitvals(1);
        fluo_corrected = dout(:,2);
        for tj = x'
            fluo_corrected(dout(:,1)==tj) = fluo_corrected(dout(:,1)==tj) ./ corrvals(x==tj);
        end

        %cluster the data into 2 populations
        cell_classes = clusterdata([dout(:,1),fluo_corrected],'Linkage','ward','SaveMemory','on','Maxclust',2);

    end   
    
    [~,hig_fluo_ind] = max([mean(dout(cell_classes==1,2)), mean(dout(cell_classes==2,2))]);
    high_fluo_pop = dout(cell_classes==hig_fluo_ind,3);

    %prepare the output matrix
    res = zeros(length(mCells),1);
    for tc = mCells'  
        fluo_frames = tc.fluoIntensities(fluo_chan_index,:);
        fluo_lifetime = sum(~isnan(fluo_frames));
        if sum(high_fluo_pop == tc.id) >= 0.5*fluo_lifetime
            res(tc.id) = 1;
        end
    end

    if plot_it
        gscatter(dout(:,1)+randn(length(dout(:,2)),1)./2,dout(:,2),cell_classes,['k','r'],['.','.'],[2,2])
        drawnow
    end

end

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
            nexttile
            hold on
            plot(t,a,'ok')
            plot(t,ft2(t),'b')
            xlabel('frame')
            ylabel('size (px^2)')
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
