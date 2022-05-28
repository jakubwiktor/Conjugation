root = '/hdd2/RecBCD2/codedev/Analysis/EXP-22-BY4442/therun/prod';

poslist = dir(fullfile(root,'Pos*'));
fluo_chan_idx = 2;

% poslist = poslist(2);
for p = poslist'
        
    load(fullfile(p.folder,p.name,'trackedCells.mat'))
    
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
   
    %try to correct for bleaching
    ft = fit(fitdata(:,1),fitdata(:,2),'exp2');
    x = unique(dout(:,1));
    fitvals = ft(x);
    corrvals = fitvals./fitvals(1);
    fluo_corrected = dout(:,2);
    for tj = x'
        fluo_corrected(dout(:,1)==tj) = fluo_corrected(dout(:,1)==tj) ./ corrvals(x==tj);
    end

%     plot(dout(:,1),fluo_corrected,'.k')
%     nexttile

    T = clusterdata([dout(:,1),fluo_corrected],'Linkage','ward','SaveMemory','off','Maxclust',2);

    [~,hig_fluo_ind] = max([mean(dout(T==1,2)), mean(dout(T==2,2))]);
    high_fluo_pop = dout(T==hig_fluo_ind,3);
    
    res = zeros(length(mCells),1);
    for tc = mCells'  
        fluo_frames = tc.fluoIntensities(fluo_chan_idx,:);
        fluo_lifetime = sum(~isnan(fluo_frames));
        if sum(high_fluo_pop == tc.id) > 0.5*fluo_lifetime
            res(tc.id) = 1;
        else
            res(tc.id) = 0;
        end
    end
    
end
    
        
%     gscatter(dout(:,1),dout(:,2),T,['k','r'],['.','.'],[2,2])
%     drawnow

end
