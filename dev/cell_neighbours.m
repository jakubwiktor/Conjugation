%local code
codePath = '/hdd2/RecBCD2/codedev/ImAnalysis/';
addpath(codePath);
ImAnalysis_setup();
expInfoObj_base= '/hdd2/RecBCD2/codedev/Analysis/EXPHANDLE/therun/prod/expInfoObj.mat';
experiments = {'EXP-22-BY4440','EXP-22-BY4442','EXP-22-BY4444','EXP-22-BY4448'};

for expnum = 1:length(experiments)
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
    out = {};
    res = {};
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
        
        Tt = T(T.posnum==pj,:);
        
        cell_contacts = {}; %store the contacts between cells

        %         parfor (vi = 1 : length(expInfoObj.imRange{pj, 1}), 5)
        parfor vi = 1 : length(expInfoObj.imRange{pj, 1})
            
            cell_contacts{vi} = nan(length(mCells),1);
            
            seg  = imread(fullfile(pathSeg, listSeg{vi}));
            
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
%             hold on
            overlap_threshold = 50;
            
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
                    
                %remove cells with small contact threshold
                good_overlap_index = zeros(length(neigh_ids),1);
                for ni = 1:length(neigh_ids)
                    overlap_sum = sum(neigh_vals == neigh_ids(ni));
                    if overlap_sum < overlap_threshold, continue, end
                    good_overlap_index(ni) = 1;
                end
                neigh_ids(~good_overlap_index) = [];
                
                %indexes of good cells, and blobs, maybe put this on one matrix not to confuse it later?
                neigh_cell_ids = all_blobs(ismember(all_blobs(:,2),neigh_ids),1);
                good_blobs = all_blobs(ismember(all_blobs(:,2),neigh_ids),2);
                                
                this_center = props(all_blobs(tb,2)).Centroid;
                centers = vertcat(props(good_blobs).Centroid);
                
                this_cell = all_blobs(tb,1);
                contacting_cells = all_blobs(ismember(all_blobs(:,2),neigh_ids),1);
                
%                 if ~(Tt.venus(this_cell) & ~Tt.cherry(this_cell)), continue, end
%                 if ~all(Tt.cherry(ismember(Tt.id,contact_cells))==0), continue, end
               
                if any(Tt.cherry(contacting_cells) & ~Tt.venus(contacting_cells)) %only true donors, not conjugants
                   cell_contacts{vi}(this_cell) = sum(Tt.cherry(contacting_cells) & ~Tt.venus(contacting_cells));
                end

%                 if Tt.venus(this_cell) && Tt.cherry(this_cell)
%                     plot(this_center(1), this_center(2), '.m','MarkerSize',16)
%                 elseif Tt.venus(this_cell) && ~Tt.cherry(this_cell)
%                     plot(this_center(1), this_center(2), '.g','MarkerSize',16)
%                 elseif ~Tt.venus(this_cell) && Tt.cherry(this_cell)
%                     plot(this_center(1), this_center(2), '.r','MarkerSize',16)
%                 end
                
%                 for cc = 1:size(centers,1)
%                     
%                     if Tt.venus(Tt.id==neigh_cell_ids(cc)) && Tt.cherry(Tt.id==neigh_cell_ids(cc))
%                         plot([this_center(1) centers(cc,1)], [this_center(2) centers(cc,2)],'m','LineWidth',2)
%                     elseif Tt.venus(Tt.id==neigh_cell_ids(cc)) && ~Tt.cherry(Tt.id==neigh_cell_ids(cc))
%                         plot([this_center(1) centers(cc,1)], [this_center(2) centers(cc,2)],'g','LineWidth',2)
%                     elseif Tt.cherry(Tt.id==neigh_cell_ids(cc)) && ~Tt.venus(Tt.id==neigh_cell_ids(cc))
%                         plot([this_center(1) centers(cc,1)], [this_center(2) centers(cc,2)],'r','LineWidth',2)
%                     end
%                     
%                 end
            end
            
        end
        
        res = nanmean(horzcat(cell_contacts{:}),2);
        %what are the cell fractions? sum(~Tt.cherry(find(isnan(res)))) / sum(isnan(res))
        % isnan(res) - is the cells that never had a nighobour that was +cherry
        % ~isnan(res) - cells that had at least one cherry neighbour
        donors_with_contact = ~isnan(res) & Tt.venus;
        conj = sum(Tt.cherry(donors_with_contact)) / sum(donors_with_contact); %those are cells that had cherry neighbour and got conjugated
        
        fprintf('%s,%0.4f\n',experiments{expnum},conj)
        
    end
end

