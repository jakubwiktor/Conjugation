%calculate doubling times (growth rates) for cell selected on the first
%fraome - the ones that were growing next to the donors, and the ones that
%were growing NOT in contact with any donor on the frame 1. The assumption
%is that there is a lot of conjugation initially, and then later maybe
%less.


codePath = '/hdd2/RecBCD2/codedev/ImAnalysis/';
addpath(codePath);
ImAnalysis_setup();
expInfoObj_base= '/hdd2/RecBCD2/codedev/Analysis/EXPHANDLE/therun/prod/expInfoObj.mat';
experiments = {'EXP-22-BY4440','EXP-22-BY4442','EXP-22-BY4444','EXP-22-BY4448'};

res_all = {};

for expnum = 2%1:length(experiments)
% for expnum = 1:2
    
    expInfoObj_path = strrep(expInfoObj_base, 'EXPHANDLE', experiments{expnum});
    
    load(expInfoObj_path)
    
    numpos = length(expInfoObj.positions);
    
    channel_name_mCherry = expInfoObj.fluoChannelNames{contains(expInfoObj.fluoChannelNames,'594')};
    channel_name_SSB = expInfoObj.fluoChannelNames{contains(expInfoObj.fluoChannelNames,'514') | contains(expInfoObj.fluoChannelNames,'venus')};
    SSB_channel_index = find(contains(expInfoObj.fluoChannelNames,'514') | contains(expInfoObj.fluoChannelNames,'venus'));
    
    % save_dir_base = fullfile('/crex/proj/uppstore2018129/elflab/Projects/CRISPR_conjugation/codedev/data', experiments{expnum});
    % save_dir_base = fullfile('/hdd2/RecBCD2/codedev/Analysis/output', experiments{expnum});
    % mkdir(save_dir_base)
    T = readtable(fullfile('/hdd2/RecBCD2/codedev/data',experiments{expnum},'cell_data.csv'));

    res_neigh = {};
    res_noneigh = {};
    
    for pj = 1 : numpos
        
        posName = expInfoObj.positions{pj};
        sourceDir = expInfoObj.getPathName('source');
        outputDir = expInfoObj.getPathName('output');
        mCellsMatFile = expInfoObj.getMCellMatPath(posName);
        mCells = Cell.MCell.loadMCells(mCellsMatFile);
        mFrames = Cell.MCell.mFrames;
        
        outputDir = expInfoObj.getPathName('output');
        mCellsMatFile = expInfoObj.getMCellMatPath(posName);
        mCells = Cell.MCell.loadMCells(mCellsMatFile);
        
        pathSeg = fullfile(outputDir,posName,'SegmentedPhase');
        listSeg = dir(fullfile(pathSeg,'*.tif*'));
        listSeg = {listSeg(:).name};
        
        pathVenus = fullfile(sourceDir,posName,'fluor_YFP_venus');
        listFluo = dir(fullfile(pathVenus,'*.tif*'));
        listFluo = {listFluo(:).name};
        Ts = expInfoObj.getAllT(posName, 'fluor_YFP_venus');
        fluo = imread(fullfile(pathVenus, listFluo{1}));
        Tarr = Ts(:,:,1);
    
        Tt = T(T.posnum==pj,:);
        
        cell_contacts = {}; %store the contacts between cells

        %         parfor (vi = 1 : length(expInfoObj.imRange{pj, 1}), 5)
        for vi = 1% : length(expInfoObj.imRange{pj, 1})
            
            cell_contacts{vi} = nan(length(mCells),1);
            
            seg  = imread(fullfile(pathSeg, listSeg{vi}));
                  
            fluo = imread(fullfile(pathVenus, listFluo{vi}));
            transFluoIm = imwarp(fluo,affine2d(Tarr),'outputView',imref2d(size(seg)));

            %vectorize cell selection - 20x faster than loops
            props = regionprops(seg,'PixelList');
            
            %calculate links between blobs and cells in this frame
            all_blobs = zeros(length(mCells)-1,2);
            for ci = 1 : length(mCells)
                tt = mCells(ci).birthFrame : mCells(ci).lastFrame;
                if ismember(vi,tt)
                    thisblob = mCells(ci).blobLabels(tt==vi);
                    
                    if thisblob == 0 %sometimes cells can disappear
                        continue,
                    end
                    
                    %cell id, blobNumber
                    all_blobs(ci,:) = [mCells(ci).id thisblob];
                    
                end
            end
            all_blobs(all(all_blobs == 0,2),:) = []; %remove empty rows
          
            se1 = strel('disk',5);
            se2 = strel('disk',9);
            
            props = regionprops(seg,'Centroid','BoundingBox');
            fatseg =  imdilate(seg,se1);
            
%             imshow(seg>0,[])
%             imshowpair(imadjust(uint16(transFluoIm)),edge(seg>0))

%             hold on
            overlap_threshold = 200;
            
            %this loops computes connectin with cells in mCells struc.
            %Sometime a cell can be a donor, but not being tracket, then
            %the loop omits it. It may introduce unexpected results to the
            %analysis. Next step is to try: plot only the contacts between
            %donor and recipients.
            
            %Check how many cells are in contact with donors, and then how
            %many of those became conjugated - is there a difference in
            %induced crispr??
            
            for tb = 1:1:size(all_blobs,1)
                                
                % 'imcrop' seems to be handling condition when bbox is outside the limits of the size of the image
                bbox = props(all_blobs(tb,2)).BoundingBox;
                bbox = bbox + [-10 -10 20 20];
                cell_im = imcrop(seg,bbox)==all_blobs(tb,2);
                cell_im = imdilate(cell_im,se2);
                
                neigh_im = imcrop(fatseg,bbox);
                neigh_vals = neigh_im(cell_im);

                neigh_ids = unique(neigh_vals(neigh_vals~=0 & neigh_vals~=all_blobs(tb,2))); %is this slowing it down? Do unique before to shrink the matrix?
                
                if isempty(neigh_ids), continue, end
                                
                neigh_ids = [neigh_ids zeros(length(neigh_ids),1)];
                for ni = 1:size(neigh_ids,1)
                    overlap_sum = sum(neigh_vals == neigh_ids(ni,1));
                    neigh_ids(ni,2) = overlap_sum;
                end

                %indexes of good cells, and blobs, maybe put this on one matrix not to confuse it later?
%                 neigh_cell_ids = all_blobs(ismember(all_blobs(:,2),neigh_ids(:,1)),1);
                good_blobs = all_blobs(ismember(all_blobs(:,2),neigh_ids(:,1)),2);
                                
                this_center = props(all_blobs(tb,2)).Centroid;
                centers = vertcat(props(good_blobs).Centroid);
                
                this_cell = all_blobs(tb,1);
                [~,~,contained_blobs] = intersect(all_blobs(:,2),neigh_ids(:,1));
                neigh_ids = neigh_ids(contained_blobs,:);
                contacting_cells = all_blobs(ismember(all_blobs(:,2),neigh_ids(:,1)),1);

                if ~Tt.venus(this_cell), continue, end
                                        
                cells_over_overlap_thresh = contacting_cells(neigh_ids(:,2) > overlap_threshold);
                if any(~Tt.venus(cells_over_overlap_thresh))
                   cell_contacts{vi}(this_cell) = 1;
%                    plot(this_center(1), this_center(2), '.m','MarkerSize',16)                
                elseif all(Tt.venus(contacting_cells))
                   cell_contacts{vi}(this_cell) = 0;
%                    plot(this_center(1), this_center(2), 'xm','MarkerSize',16)
                end
            end
            
        end
        aaa = find(cell_contacts{vi}==1)';
        break
       %here are recipients that were close to donors (defined by the
       %overlap threshold, but better would be to compare the angles and
       %select the ones that have a difference in angles close to 0?
       %        figure
       tmp_mat = [];
       selected_recipients = find(cell_contacts{1}==1);
       for reci_iterator = 1:length(selected_recipients)
           this_recipient = mCells(selected_recipients(reci_iterator));
           if this_recipient.isBadCell ~= 0, continue, end
           aa = this_recipient.areas;
           if length(aa) < 10, continue, end
           tt = this_recipient.birthFrame:this_recipient.lastFrame;
           [ft, gof] = fit(tt',aa','exp1');
           if gof.rsquare < 0.5, continue, end
           tmp_mat = [tmp_mat ft.b];
       end
       res_neigh{pj} = tmp_mat;
       
       tmp_mat = [];
       selected_recipients = find(cell_contacts{1}==0);
       for reci_iterator = 1:length(selected_recipients)
           this_recipient = mCells(selected_recipients(reci_iterator));
           if this_recipient.isBadCell ~= 0, continue, end
           aa = this_recipient.areas;
           if length(aa) < 10, continue, end
           tt = this_recipient.birthFrame:this_recipient.lastFrame;
           [ft, gof] = fit(tt',aa','exp1');
           if gof.rsquare < 0.5, continue, end
           tmp_mat = [tmp_mat ft.b];
       end
       res_noneigh{pj} = tmp_mat;
    end
    
    gr_neigth = horzcat(res_neigh{:});
    gr_noneigth = horzcat(res_noneigh{:});
    res_all{expnum}{1} = gr_neigth;
    res_all{expnum}{2} = gr_noneigth;
%    
%     nexttile
%     subplot(1,2,1)
%     hold on
%     histogram(log(2)./gr_neigth, 20:5:100,'normalization','pdf','EdgeColor','r','DisplayStyle','stairs','LineWidth',2)
%     histogram(log(2)./gr_noneigth, 20:5:100,'normalization','pdf','EdgeColor','k','DisplayStyle','stairs','LineWidth',2)
%     legend('donor contact','no donor contact')
%     xlabel('doubling time (min)')
%     ylabel('pfd')
%     
%     subplot(1,2,2)
%     out = [log(2)./gr_neigth'; log(2)./gr_noneigth'];
%     lbls = [repmat('neigh   ',[length(gr_neigth),1]);repmat('no neigh',[length(gr_noneigth),1])];
%     boxplot(out,lbls,'notch','on')
%     ylim([20 100])
%     ylabel('doubling time (min)')

end

%%
out = [];
lbls = {};
for expnum = 1:length(experiments)
    exp_prefix = strsplit(experiments{expnum},'-');
    exp_prefix = exp_prefix{end};
    
    out =  [out; [res_all{expnum}{1}';  res_all{expnum}{2}']];
    lbls = [lbls; [repmat({fullfile(exp_prefix,'neigh')}, size(res_all{expnum}{1}')) ;...
                   repmat({fullfile(exp_prefix,'noneigh')}, size(res_all{expnum}{2}'))]];
end
lbls = categorical(lbls);
figure
boxplot(log(2)./out,lbls,'notch','on')
ylim([20 100])
set(gca,'XTickLabelRotation',30)
grid on
grid minor
