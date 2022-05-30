cell_classes = res{2};
load( '/hdd2/RecBCD2/codedev/Analysis/EXP-22-BY4448/therun/prod/expInfoObj.mat')
stack_save_dir = '/hdd2/RecBCD2/codedev/stack';
posName = expInfoObj.positions{2};
make_stack(expInfoObj, posName, cell_classes)

function make_stack(expInfoObj, posName, cell_classes)
% Construct RGB image for each segmentation image and color cells depending
% on the selection criteria. Green - 'venus high' cells, red - 'mcherry'
% cells, and yellow - 'transconjugated' cells. Cells with just an outline
% and dark inside are classified as dead/stressed. 
%
% inupt: expInfoObj - expInfo structure from ImAnalysis pipeline
%        posName - string, name of the position as in expInfoObj.positions
%        
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
   
    %vectorize cell selection - fast way to reconstruct image
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
    
%     imshow(rgb)
%     drawnow
    
    %write image
    imwrite(rgb, fullfile( '/hdd2/RecBCD2/codedev/stack', [num2str(vi) '.tiff']))
    
end
end