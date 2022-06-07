
load('/hdd2/RecBCD2/codedev/Analysis/EXP-22-BY4442/therun/prod/expInfoObj.mat')

% store cell id, cell images, ourline, frame indexes, 

%cells slected from 'first_frame_contacts.m' for this experiment and position
%EXP-22-BY4442, pos1
% cellrange = [4,8,15,25,35,40,41,46,49,51,53,74,83,123,141,199,222,285,309,327,335,343,366,372,430,438,461,474,483,489,546,547,620,653,710,735,737,823,1020,1163,1283,1336,1350,1959,2331,2408,2830,3009,3076,3492,3499,3501];
% posname = 'Pos11';
% pos2
cellrange = [1,2,3,8,9,11,23,24,31,35,41,43,47,52,56,89,90,134,138,162,182,183,185,203,239,292,297,314,323,324,326,327,374,385,412,439,454,471,494,513,515,538,539,606,607,669,894,1022,1024,1595,2055,2084,2085,2373,2375,2418,2494,2497,2574,2881,3012,3017,3154,3155];
posname = 'Pos12';
%EXP-22-BY4440, pos1
% cellrange = [26,30,37,50,55,79,95,169,196,197,256,263,272,276,298,314,334,349,350,373,374,391,466,468,539,566,596,597,638,695,709,725,727,742,744,845,863,867,919,940,942,951,1034,1087,1088,1148,1167,1168,1169,1278,1281,1299,1330,1333,1335,1358,1428,1446,1447,1448,1473,1489,1548,1623,1639,1754,1798,1815,1823,1825,1881,2141,2177,2411,2629,2661,2937,3071,3145,3264,3268,3269,3273,3374,3377,3387,3473];
% posname = 'Pos11';
%EXP-22-BY4444, pos1
% cellrange = [4,5,6,9,16,19,22,25,29,30,31,45,48,63,64,65,66,72,75,76,77,78,81,83,89,110,112,113,149,180,216,217,218,250,275,278,287,312,317,320,335,338,360,362,373,375,376,388,396,399,418,419,422,427,428,429,447,491,503,590,600,604,615,617,618,648,654,699,707,708,741,776,792,808,832,873,957,979,1045,1046,1120,1146,1265,1287,1304,1326,1737,1738,1789,2240,2289,2351,2357,2359,2360,2528,2576,2659,2660,2700,2709,2736,2908,2909,2910,3176,3179,3429,3430,3432,3660];
% posname = 'Pos11';
%EXP-22-BY4448 pos1
% cellrange = [6,11,13,14,23,34,40,58,62,68,69,79,83,85,128,158,164,203,210,211,252,261,264,279,295,298,325,331,366,413,437,438,515,516,533,534,542,565,595,632,675,705,731,771,780,803,865,963,977,1048,1049,1138,1139,1248,1280,1302,1326,1348,1435,1462,1489,1510,1589,1639,1755,1829,1830,2089,2093,2134,2136,2140,2214,2341,2403,2811,2815];
% posname = 'Pos11';

mCellsMatFile = expInfoObj.getMCellMatPath(posname);
mCells = Cell.MCell.loadMCells(mCellsMatFile);

cellrange = intersect(cellrange,find([mCells.isBadCell]==0 & [mCells.lifeTime]>20 & [mCells.lifeTime]<70));
% cellrange = find([mCells.isBadCell] == 0, 25);

% celltable = readtable('/hdd2/RecBCD2/codedev/data/EXP-22-BY4442/cell_data.csv');
% cellrange = celltable.id(celltable.posnum==1 & celltable.venus & celltable.cherry &~celltable.dead);
% cellrange = cellrange(1:50)'; %14,27,42,45,58,67

% cellrange = intersect(cellrange,find([mCells.isBadCell]==0));


resCherry = find_cells_with_spots(expInfoObj,posname,'fluor_594_cherry',cellrange);
resYFP = find_cells_with_spots(expInfoObj,posname,'fluor_YFP_venus',cellrange);

% reconstruct 2-color images
for ci = 1 : length(resYFP)
    
    %i know YFP is imaged more frequently, so I loop through every frame
    %from min(YFP):max(YFP) and once it is also in mCehrry, I add it.
    thisYFP = resYFP{ci};
    thisCherry = resCherry{ci};
    
    ysize = cellfun(@(x) size(x,1), thisYFP(:,1));
    xsize = cellfun(@(x) size(x,2), thisYFP(:,1));
    montageYFP = zeros(max(ysize),sum(xsize));
    montageCherry = zeros(max(ysize),sum(xsize));
    xstart = 1;
    
    bbs = {};
    
    for ij = 1:length(thisYFP(:,1))
        imi = thisYFP{ij,1};
        ystart = max(1, floor((max(ysize)-size(imi,1))/2));
        montageYFP(ystart : size(imi,1) + ystart - 1 , xstart : xstart + size(imi,2)-1) = imi;
        bbs{ij} = [thisYFP{ij,2}(:,1) + xstart, thisYFP{ij,2}(:,2) + ystart];
        
        if ismember(thisYFP{ij,3},[thisCherry{:,3}])
            imc = thisCherry{ismember([thisCherry{:,3}],thisYFP{ij,3}),1};
            montageCherry(ystart : size(imc,1) + ystart - 1 , xstart : xstart + size(imc,2)-1) = imc;
        end
        xstart = xstart+size(imi,2);
        
    end
    
    imyfp     = imadjust(uint16(montageYFP),stretchlim(uint16(montageYFP),[0.01 0.999]));
    imycherry = imadjust(uint16(montageCherry),stretchlim(uint16(montageCherry),[0.01 0.9999]));
    
    figure
    imshow([imyfp;imycherry])
    hold on
    for ij = 1:length(thisYFP(:,1))
        plot(bbs{ij}(:,1), bbs{ij}(:,2),'w')
        plot(bbs{ij}(:,1), bbs{ij}(:,2)+size(imyfp,1),'w')
    end
    pause(0.1)

end


function res = find_cells_with_spots(expInfoObj,pos_name,fluoChanName,cellrange)
% make cell montages

%define constant for radian <-> degrees conversion
radconst = (pi/180);
res = {};

pos_index = find(strcmp(expInfoObj.getPositionList,pos_name));
if isempty(pos_index)
    error(['Wrong position: ' pos_name])
end

chan_index = find(strcmp(expInfoObj.channels,fluoChanName));
if isempty(chan_index)
    error(['Wrong channel: ' fluoChanName])
end

if size(cellrange,1) > 1
    error('Wrong shape of cell range, shuld be 1-by-n, where n is cell indexes')
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


%load the fluo stack
imstack = {};
for vi = 1:length(expInfoObj.imRange{pos_index, chan_index})
    fluo = imread(fullfile(pathFluo, listFluo{vi}));
    seg  = imread(fullfile(pathSeg, listSeg{vi}));
    this_frame = fluoIndices(vi);
    T = Ts(:,:,phaseRange(vi));
    transFluoIm = imwarp(fluo,affine2d(T),'outputView',imref2d(size(seg)));
    imstack{vi} = transFluoIm;
end

res = {};
for ci = 1:length(cellrange)
    
    thiscell = mCells(cellrange(ci));

%     if thiscell.isBadCell ~= 0, continue, end
    
    cell_life = thiscell.birthFrame : thiscell.lastFrame;
    
    [~,life_index, fluo_index] =  intersect(cell_life,fluoIndices);
    
    bboxes = thiscell.boundingBoxes;
    outlines = thiscell.boundaries;
    
    
    for bbi = 1:length(fluo_index)
        
        %
        %find out the orientation of the cell and crop the image
        %
        
        if isempty(outlines{life_index(bbi)})
            res{ci}{bbi,1} = 0;
            res{ci}{bbi,2} = [0 0];
            res{ci}{bbi,3} = life_index(bbi);
            res{ci}{bbi,4} = cellrange(ci);
            continue
        end
        
        x = outlines{life_index(bbi)}(:,1);
        y = outlines{life_index(bbi)}(:,2);
        trans_vals = [mean(x) mean(y)];
        
        x = x - trans_vals(1);
        y = y - trans_vals(2);
        
        tmp = zeros(180,1);
        
        for theta = 1:180
            xtmp = x*cosd(theta) - y*sind(theta);
            tmp(theta) = max(xtmp)-min(xtmp);
        end
        
        [~,thetam] = min(tmp);
        
        %this parts deals with the rotationi of the entire image and cell
        %boundaries. It crops a new bounding box on the rotated call
        %boundaries.
        bbox = double(bboxes(:,life_index(bbi)));
        
        im = imstack{fluo_index(bbi)};
        thisim = imrotate(im, -thetam, 'loose');

        mbx = mean(outlines{life_index(bbi)}(:,1)+bbox(1)) - size(im,2)/2;
        mby = mean(outlines{life_index(bbi)}(:,2)+bbox(2)) - size(im,1)/2;

        x0 = (mbx*cosd(thetam) - mby*sind(thetam)) + size(thisim,2)/2;
        y0 = (mbx*sind(thetam) + mby*cosd(thetam)) + size(thisim,1)/2;
        x1 = [x*cosd(thetam) - y*sind(thetam)] + x0;
        y1 = [x*sind(thetam) + y*cosd(thetam)] + y0;

%         figure
%         imshow(thisim,prctile(thisim(:),[0 99]));
%         hold on
%         plot(x1,y1,'g')
        
        padvalue = 11;
        cropbox = [round(min(x1)-padvalue),...
                   round(min(y1)-padvalue),...
                   round(max(x1)-min(x1)) + padvalue*2,...
                   round(max(y1)-min(y1)) + padvalue*2];
               
        imcropped = imcrop(thisim, cropbox);
%         imcropped = thisim(round(min(y1)-padvalue) : round(max(y1)+padvalue) ,...
%                            round(min(x1)-padvalue) : round(max(x1)+padvalue)); 

        %group images in a cell - cell index, then rotated image, adjusted
        %boundaries, frame, cell index
        res{ci}{bbi,1} = imcropped;
        res{ci}{bbi,2} = [x1-min(x1)+padvalue y1-min(y1)+padvalue];
        res{ci}{bbi,3} = life_index(bbi);
        res{ci}{bbi,4} = cellrange(ci);

    end

end


end
