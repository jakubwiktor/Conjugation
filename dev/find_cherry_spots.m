
load('/hdd2/RecBCD2/codedev/Analysis/EXP-22-BY4448/therun/prod/expInfoObj.mat')

out = findCellsWithSpots(expInfoObj,'Pos12','fluor_594_cherry')

function cell_indexes = findCellsWithSpots(expInfoObj,pos_name,fluoChanName)
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

%for gaussian fitting initial guesses: this is radius of searched spots
gauss_kren_size = radius;
fttype = fittype( @(a,b1,s1,b2,s2,x,y) a*exp(-(((x-b1).^2)/(2*s1^2) + ((y-b2).^2)/(2*s2^2))) ,'numindep', 2);
[X,Y] = meshgrid(1:1+2*gauss_kren_size);

cell_indexes = {};

for vi = 10%:length(expInfoObj.imRange{pos_index, 2})
    
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
    
    figure
    imshow(transFluoIm,prctile(double(transFluoIm(:)),[0, 99]))
    hold on
    plot(xx,yy,'rx')
    
    xx = xx(sel_condition);
    yy = yy(sel_condition);
    
    plot(xx,yy,'gx')
    disp('hey')

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

cell_indexes = unique(vertcat(cell_indexes{:}));

end
