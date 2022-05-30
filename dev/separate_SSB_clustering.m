load('/hdd2/RecBCD2/codedev/Analysis/EXP-22-BY4442/therun/prod/Pos12/trackedCells.mat');

fluo_chan_idx = 2;

select_fluo_cells(mCells,2,1)

function res = select_fluo_cells(mCells, fluo_chan_idx, plotit)
    %  Separating cells in the mCells structure into two population - high
    %  fluorescence, and low fluorescence. For the channel specified by
    %  fluo_chan_ids as in mCells(n).fluoIntensities. The cell is qualified
    %  as fluorescnt if during at least half of the frames it was 
    %  classified as high fluorescent. 
    %
    %  input: 
    %       mCells - Cell.MCell object
    %   
    %       fluo_chan_index - index of fluorescent channel to analyse as in
    %           mChells.fluoIntensities field
    %
    %    output:
    %       n by 2 matrix where first column is cell id (as in mCells.id)
    %       and the 2n is a class - 1 - fluoresncent, 0 - not fluorescent
    
    assert(~isempty([mCells.fluoIntensities]), 'fluo intensities are empty')
        
    if nargin<3
        plot_it = 0;
    end
    
    dout = [];
    fitdata = [];

    for tc = mCells'

        f = tc.fluoIntensities(fluo_chan_idx,:);

        gi = ~isnan(f);

        f = f(gi);
        t = tc.birthFrame:tc.lastFrame;
        t = t(gi);

        dout = [dout; [t' f' repmat(tc.id,[length(t), 1])]];

        if tc.isBadCell==0 
            [~, idxs] = rmoutliers(f);
            fitdata = [fitdata; [t(~idxs)' f(~idxs)']];
        end

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
    T = clusterdata([dout(:,1),fluo_corrected],'Linkage','ward','SaveMemory','off','Maxclust',2);

    [~,hig_fluo_ind] = max([mean(dout(T==1,2)), mean(dout(T==2,2))]);
    high_fluo_pop = dout(T==hig_fluo_ind,3);

    %prepare the output matrix
    res = zeros(length(mCells),1);
    for tc = mCells'  
        fluo_frames = tc.fluoIntensities(fluo_chan_idx,:);
        fluo_lifetime = sum(~isnan(fluo_frames));
        if sum(high_fluo_pop == tc.id) >= 0.5*fluo_lifetime
            res(tc.id) = 1;
        end
    end

    if plotit
        gscatter(dout(:,1),dout(:,2),T,['k','r'],['.','.'],[2,2])
        drawnow
    end
end