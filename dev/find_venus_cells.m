
listpos = dir(fullfile('/hdd2/RecBCD2/codedev/Analysis/EXP-22-BY4447/therun/prod/','Pos*'))

fluo_chan_index = 2;
figure

for ii = listpos'
    ii.name
    load(strrep('/hdd2/RecBCD2/codedev/Analysis/EXP-22-BY4447/therun/prod/CHANGE/trackedCells.mat','CHANGE',ii.name))


    nexttile
    [out1] = select_fluo_cells(mCells, fluo_chan_index, 'stdev', 1);
    nexttile
    [out2] = select_fluo_cells(mCells, fluo_chan_index, 'gaussians', 1);
    nexttile 
    [out3] = select_fluo_cells(mCells, fluo_chan_index, 'clustering', 1);
    break
end

function [res] = select_fluo_cells(mCells, fluo_chan_index, algorithm, plot_it)
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
    %        algorithm - 'stdev,'gaussian','clustering'
    %
    %    output:
    %       indexes of cells with fluoresncence above threshold as in
    %       mCells(i).id

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
        gscatter(dout(:,1),dout(:,2),cell_classes,['k','r'],['.','.'],[2,2])
        drawnow
    end

end

