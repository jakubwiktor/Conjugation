load('/hdd2/RecBCD2/codedev/Analysis/EXP-22-BY4448/therun/prod/Pos11/trackedCells.mat')
fluo_chan_index = 2;
[~, out] = select_fluo_cells(mCells, fluo_chan_index, 2, 'stdev', 0);

function [cell_inds, out] = select_fluo_cells(mCells, fluo_chan_index, min_fluo_time, algorithm, plot_it)
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
    
    out = [t f];
    
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
