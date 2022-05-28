%this version of the code was used to plot moveis of the recipient and
%donor strains on 19th May 2022.

codePath = '/hdd2/RecBCD2/codedev/ImAnalysis/';
addpath(codePath);
ImAnalysis_setup(); 

expInfoObj_path = '/hdd2/RecBCD2/codedev/Analysis/EXP-22-BY4448/therun/prod/expInfoObj.mat';
load(expInfoObj_path)

numpos = length(expInfoObj.positions);

fluoChanName = expInfoObj.fluoChannelNames{2};

disp(fluoChanName)

% cherryT = load('/crex/proj/uppstore2018129/elflab/Projects/CRISPR_conjugation/codedev/matlab/cherry_to_venus_trans_20220513.mat');

for pj = 1 %: numpos
    
    %initialize experiment parameters
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
    listSeg = {listSeg(fluoIndices).name};

    
    cell_indexes_venus = select_fluo_cells(mCells,2,2,'stdev',0);
    
    % tic
    % cell_indexes_venus = select_fluo_cells(mCells,2,2,'gaussians',1);
    % toc

    cell_indexes_cherry = select_spots_cells(mCells,'fluor_594_cherry',2);
    
    [dead_indexes_dead] = find_dying_cells(mCells);

    for vi = 100 : length(expInfoObj.imRange{pj, find(strcmp(expInfoObj.getChannelNames,fluoChanName))})
        
        fluo = imread(fullfile(pathFluo, listFluo{vi}));
        seg  = imread(fullfile(pathSeg, listSeg{vi}));
        this_frame = fluoIndices(vi);
        
        seg_selected_cherry = zeros(size(seg));
        seg_selected_venus = zeros(size(seg));
        
        for blob = mFrames.getBlobsInFrame(this_frame)
            this_cell = mFrames.getCellsInBlob(this_frame,blob).id; %what if there are 2 cells in one blob?
            
%             if mCells(this_cell).isBadCell ~= 0
%                 continue
%             end
            
            %select cells with spots
            if ismember(this_cell,cell_indexes_cherry) 
                seg_selected_cherry(seg==blob) = 1;
            end
            
             %select cells with venus intensity
            if ismember(this_cell,cell_indexes_venus) 
                seg_selected_venus(seg==blob) = 1;
            end
                        
        end
        
        %find out which cells belong to which blob in the seg image. 
        
        T = Ts(:,:,phaseRange(vi));
        
        transFluoIm = imwarp(fluo,affine2d(T),'outputView',imref2d(size(seg)));
        transFluoIm = transFluoIm - mean(transFluoIm(seg==0));
        

        rgb = cat(3,repmat(imadjust(transFluoIm),[1,1,3]));
        rgb(:,:,2) = min(rgb(:,:,2) + uint16(edge(seg_selected_venus>0)).*((2^16)-1), ((2^16)-1));
        rgb(:,:,1) = min(rgb(:,:,1) + uint16(edge(seg_selected_cherry>0)).*((2^16)-1), ((2^16)-1));
        imshow(rgb,[])
        
%         imwrite(uint16(rgb), fullfile('stack',[num2str(vi) '.png']))
    end
    
end

% histogram([mCells(intersect(cell_indexes_venus,cell_indexes_cherry)).birthFrame])




function cell_inds = select_spots_cells(mCells,channel_name,min_spot_time)
    
    %documentation. Its important here to keep the min_spot_time as
    % '1' because mCherry is imaged so infrequently, cells dont have much
    % time to accumulate frames with spots.
    
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

% previous version
%     %adjust by exponential model
%     fit_model = fit(t,f,'exp2');
%     corr_values = fit_model(t);
%     corr_values = corr_values ./ fit_model(1);
%     f_adj = f./corr_values;

%     %find peak by fitting 2 gaussian model to fluo values
%     [hist_val,edges] =  histcounts(f_adj);
%     hist_x = edges(1:end-1) + mean(diff(edges))/2;
%     f_gauss = fit(hist_x',hist_val', 'gauss2');
    
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
            this_f = f(t==this_t);[cell_indexes] = find_dying_cells(mCells)

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
    % documentation
    
    out = [];
    for ci = 1  : length(mCells)

        if mCells(ci).lifeTime < 5
            continue
        end

        a = mCells(ci).areas;
        t = 1:length(a);
        [~, outlierInd] = rmoutliers(a);

        a = a(~outlierInd);
        t = t(~outlierInd);

        if length(a) < 3
            continue
        end

        [~, gof] = fit(t', a', 'exp1');


        if gof.rsquare < 0.5 && ...
                mCells(ci).lifeTime > 40 && ...
                length(a) > 4 &&...
                max(abs(diff(a))) < 0.5 * (max(a) - min(a))

            qfun = fittype( @(a, b, c, x) a*x.^2 + b*x + c ); %sad fasce shape

            [ft2, gof2] = fit(t', a', qfun);

            if gof2.rsquare > 0.6 && ft2.a < 0
                out = [out mCells(ci).id];
                
%                 nexttile
%                 hold on
%                 plot(t,a,'ok')
%                 plot(t,ft(t),'b')
%                 plot(t,ft2(t),'r')
%                 title([num2str(ci) ' ' num2str(gof2.rsquare)])
%                 counter = counter + 1;

            end
        end
    end
end