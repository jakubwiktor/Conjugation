%pipeline wrapper

anal_dir = '/hdd2/RecBCD2/codedev/Analysis/EXP-22-BY4444';
exp_dir = '/hdd2/RecBCD2/EXP-22-BY4443';
exp_name = {'therun'};
parameter_file_path = fullfile('/hdd2/RecBCD2/codedev','ParameterFile.txt');

codePath ='/hdd2/RecBCD2/codedev/ImAnalysis/';
addpath(codePath);
ImAnalysis_setup(); 

for exp = 1
    tstart = tic;
    vars.this_exp = exp_name{exp};
%     vars.phase_name = 'phase';
%     vars.fluo_names = {'fluor_594_cherry','fluor_YFP_venus'};    
    vars.add_fluo_chans = {'fluor_YFP_venus'};
    vars.trans_mat_file = '';
%     vars.barcodeROI = [];
%     vars.ROI = '';
    vars.param_file_base = parameter_file_path;
    vars.code_dir = codePath;
    vars.sourceDir = fullfile(exp_dir,exp_name{exp});
    vars.outputDir = fullfile(anal_dir,exp_name{exp});    
    vars.imposeTranformation = 0; %PAY ATTENTION TO THIS!
    analyseExperiment(vars)
    tseg_end = datestr(seconds(toc(tstart)),'HH:MM:SS');
    disp(['experiment analysed in whooping: ' tseg_end])
end

function analyseExperiment(vars)
%wrapper to analyse the whole experiment
%1-Create directories and parameter file
%2-Set initial conditions
%3-Preprocess the files
%4-Segment with NestedUnet
%5-Modify parameter file and put cells in channels
%6-Track cells
%7-Adjust transformation matrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Preprocessing part
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%make analysis folder
checkMakeDir(vars.outputDir)

% 1 prepare parameter file - this part should be reworked to input pairs of
% fields with values

mod_params = {'doPreprocessing',1,...
              'doSegmentation',2,...
              'doBlobProcessing',0,...
              'doDotDetection',0,...
              'doCellMeasurements',0,...
              'rangePositions','',...
              'parNumWorkers',16};
 
param_file = fullfile(vars.outputDir,['ParameterFile_' vars.this_exp '.txt']);
modParamFile(vars.param_file_base,mod_params(1:2:end),mod_params(2:2:end),param_file) %initial parameter file

%make prod folder and preprocess the images
checkMakeDir(fullfile(vars.outputDir,'prod'))

processCells(vars.sourceDir, fullfile(vars.outputDir,'prod'), param_file, vars.code_dir)

%load expInfoObj and pass parameters structure from now on
load(fullfile(vars.outputDir,'prod','expInfoObj.mat'),'expInfoObj');
params = expInfoObj.parameters;


%-----------------------------
%
%     SEGMENTATION PART
%
%-----------------------------

seg_start = tic;

%get list of images - assume that its caalled PreprocessedPhase
all_images = dir(fullfile(vars.outputDir,'**/PreprocessedPhase/*.tif*')); %recursively, 
all_images_list = fullfile({all_images.folder},{all_images.name},':');

image_packet_size = 64;
chunks = 0 : image_packet_size : length(all_images_list); %divide imagaes in packets of 64

if chunks(end) ~= length(all_images_list)
    chunks = [chunks length(all_images_list)];
end

%make folders, assume that we switch PreprocessedPhase to SegmentedPhase
preproc_folder_names = unique({all_images.folder});
seg_folder_names = strrep(preproc_folder_names,'PreprocessedPhase','SegmentedPhase');
for each_folder = seg_folder_names
    checkMakeDir(char(each_folder))
end

% parfor (nij = 1:length(chunks)-1, 6)      %rackham
for nij = 1:length(chunks)-1                %if snowy then just single core and GPU
    
    %run python - write a script to pass the arguments to the data loader
    this_images_concat = [all_images_list{chunks(nij)+1 : chunks(nij+1)}];
    this_images_concat = this_images_concat(1:end-2);
    numcores = 0;

    system_command = ['python3 ',...
                      '/hdd2/RecBCD2/codedev/python/omniposeModel/main.py ',...
                      this_images_concat,...
                      ' PreprocessedPhase',...
                      ' SegmentedPhase ',...
                      num2str(numcores)];
                  
    system(system_command)
end
seg_end = datestr(seconds(toc(seg_start)),'HH:MM:SS');
disp(['Segmentation finished in: ' seg_end])

%-----------------------------
%
%  SEGMENTATION POSTPROCESSING
%  release touching blobs from omnipose, tghis is a reasonably fast way to do
%  it, other solution would require rewriting the blob computation code
%
%-----------------------------
seg_postprocessing_start = tic;
disp('Cleaning up segmentation images')

all_images = dir(fullfile(vars.outputDir,'**/SegmentedPhase/*.tif*')); %recursively, 
parfor imi = 1:length(all_images)
    thisim_name = fullfile(all_images(imi).folder, all_images(imi).name);
    thisim = imread(thisim_name);
    im2 = imclearborder(thisim);
    seedmat = zeros(size(im2));
    ulabels = unique(im2);
    for ul = 2:length(ulabels) %first is background
       thiscellim = im2 == ulabels(ul);
       seedmat = seedmat + imerode(thiscellim,strel('disk',1));
    end
    distmat = double(im2>0);
    distmat(seedmat>0) = -Inf;
    w = watershed(distmat);
    im2(w==0)=0;
    im2 = bwlabel(im2>0,4);
    imwrite(uint16(im2),thisim_name,'Compression','deflate');
end
seg_postprocessing_end = datestr(seconds(toc(seg_postprocessing_start)),'HH:MM:SS');
disp(['Segmentation images cleaned in: ' seg_postprocessing_end])


%-----------------------------
%
%          TRACKING
%
%-----------------------------

params.doPreprocessing = 2;
params.doSegmentation = 2;
params.doBlobProcessing = 1;
params.doCellMeasurements = 0;
params.doDotDetection = 0;
params.shortTrackCutOff = 5;
params.centerOfMassMovementCutOff = 35;
params.parforPositions = 1;
params.parNumWorkers = 12;
processCells(vars.sourceDir, fullfile(vars.outputDir,'prod'), params, vars.code_dir)

%-----------------------------
%
%  Dot detections - but no dot tracking
%
%-----------------------------

% params.doBlobProcessing = 2;
% params.doDotDetection = 1;
% params.dotDetectionAlg = 'wavelet';
% params.dotTrackingAlg = '';
% params.cellInclusionMargin = 1;
% params.noiseThreshold = 3;
% params.minSpotArea = 4;
% params.maxAxesRatio = 1.5;
% params.parforPositions = 0;
% params.parNumWorkers = 16;
% processCells(vars.sourceDir, fullfile(vars.outputDir,'prod'), params, vars.code_dir)

%-----------------------------
%
%  FLUORESCENCE INTEGRATION
%
%-----------------------------


%add fluorescence measurements
disp('starting fluo calculations')
tFluo = tic;
for c = 1:length(vars.add_fluo_chans)
addFluorescence(fullfile(vars.outputDir,'prod'),vars.add_fluo_chans{c}, 0)
end
tFluo_end = datestr(seconds(toc(tFluo)),'HH:MM:SS');
disp(['fluo added  in ' tFluo_end])

end

%{
%adjust tranformation matrices
if vars.imposeTranformation
imposeTranformationMatrix(fullfile(vars.outputDir,'prod'), 'fluo594', 'CFP')
% imposeTranformationMatrix(fullfile(vars.outputDir,'prod'), 'fluo635', 'cy5')
end

%copy entire ImAnalysis foder and pipeling_bolts to the analysis forlder
dirs = {'/home/skynet/code/ImAnalysis','/home/skynet/code/pipeline_bolts'};
copyPipeline(fullfile(vars.outputDir,'code_version'),dirs)

end
%}



%--------------------------------------------------------------------------
%all the other functions 
%--------------------------------------------------------------------------



function modParamFile(param_file, mod_fields, mod_vals, save_name)
%open parameter file, modify it and save. 
%mod_fields - a cell of strings with fields to modify
%mod_vals - a cell of cells with strings to set at fields. Can be multiple 
%           strings in one cell.

fid = fopen(param_file);
pfile = [];
fOut = fopen(save_name, 'w');
while (~feof(fid))
   ln = fgetl(fid);
   splt_ln = strsplit(ln);
   flag = find(strcmp(mod_fields, splt_ln{1}));
   if ~isempty(flag)
       %a bit of spaggetti line but matlab is stubborn with strings
       ln = char(join([splt_ln{1} '=' string(mod_vals{flag})]));
   end   
%    pfile = horzcat(pfile,[ln '\n']);
   fprintf(fOut,'%s\n',ln);
end
% fOut = fopen(save_name, 'w');
% fprintf(fOut,'%s',pfile);
fclose(fOut);
fclose(fid);
end

function checkMakeDir(dir_path)
%checks if there is folder, or makes it if not
    if exist(dir_path,'dir') ~= 7
        mkdir(dir_path)
    end
end

function nGrowthChannels = linkCellsChannels(exp_dir)
%function to put blobs in corresponding channels after NestedUnet
%segmnentation. 
%exp_dir - is the 'prod' analysis folder.
%SegmentedChannels folder from nested unet segmentation has to be there
%
%updated on 02-06-2021 - changed to median filtering of the channel
%locations to deal with the problem where edges of the channels were
%expanded and mutliple channels were connected.

% exp_dir = fullfile(exp_dir, 'prod');
expInfoObj_placeholder = load(fullfile(exp_dir, 'expInfoObj.mat'));
expInfoObj = expInfoObj_placeholder.expInfoObj;
pos_name = expInfoObj.positions;
%chan_im = 'SegmentedChannels/img_000000000_KubNet.tiff';

num_frames = length(expInfoObj.imRange{strcmp(expInfoObj.channels,'aphase')});

chan_mat = zeros(1, 1, length(pos_name));

yoffset = 10; %specify where the channels masks are starting in the y axis.
dilateFactor = 1; %how much channels will bloat
for p = 1 : length(pos_name)
    this_ims = dir(fullfile(exp_dir,char(pos_name(p)),'SegmentedChannels','*tiff'));
    this_im = imread(fullfile(this_ims(1).folder, this_ims(1).name));
   
    %find peaks using median filtering. Should be resistant to noise
    chanMask = [0 median(this_im>0) 0]; %push 0 into the mask to isolate edge channels
    pksDiff = diff(chanMask);
    chanLocs = [find(pksDiff==1)' find(pksDiff==-1)'];
    chanWidth = median(chanLocs(:,2)-chanLocs(:,1));    
    bbox = [chanLocs(:,1)-dilateFactor,...
                repmat(yoffset,size(chanLocs,1),1),...
                chanLocs(:,2)-chanLocs(:,1)+dilateFactor*2,...
                repmat(size(this_im,1) - yoffset, size(chanLocs,1),1)]; %dilateFactor may be an error if the channel location is at 1!!
    
    %%previous method using region props
    %props = regionprops(this_im, 'BoundingBox');
    %bbox = reshape([props.BoundingBox],4,length(props))';
    chan_mat(1:size(bbox,1), 1:size(bbox,2), p) = bbox;
    %expInfoObj.setChannelLocations(char(pos_l(p)), repmat(bbox, [1,1,num_frames]))
end

%2-update parameters to match the number of channels to pass the artument
%later to the expInfObj.

nGrowthChannels = size(chan_mat,1);
expInfoObj.setParameters('nGrowthChannels', size(chan_mat,1))

for p = 1 :  length(pos_name)
    expInfoObj.setChannelLocations(char(pos_name(p)), repmat(chan_mat(:,:,p), [1,1,num_frames]))
end

expInfoObj.saveExpInfo

end

function imposeTranformationMatrix(exp_dir, fluo_chan_name, filter_cube_name)
%function to correct the pixel shift in Nikon 6. Currently its tested and
%works for the mCherry images shot with CFP cube. Can also correct mcherry
%imgaes show with venus cube. Input is the 
%exp_dir - analysis prod directory,
%fluo_chan_name - is supposed to be mCherry channel name,
%filter_cube_name - points to the transformation matrix

eio = load(fullfile(exp_dir,'expInfoObj.mat'));
expInfoObj = eio.expInfoObj;
% this_fluo = find(strcmp(expInfoObj.channels,'fluo594'));
this_fluo = find(strcmp(expInfoObj.channels,fluo_chan_name));

if strcmp(filter_cube_name,'CFP')
%     tf = load('/home/skynet/code/pipeline_bolts/tform_mcherry_imaged_on_CFP.mat');
    tf = load('/home/skynet/code/pipeline_bolts/tform_mcherry_imaged_on_CFP_20201201.mat');
elseif strcmp(filter_cube_name,'venus')
    tf = load('/home/skynet/code/pipeline_bolts/tform_mcherry_imaged_on_venus.mat');
elseif strcmp(filter_cube_name,'cy5')
    tf = load('/home/skynet/code/pipeline_bolts/tform_SYTO60_imaged_on_CFP.mat');
end

tform = tf.tform;

for p = 1:length(expInfoObj.transfMatrices)
    this_t = expInfoObj.transfMatrices{p}{this_fluo};
    for fr = 1:size(this_t,3)
        this_t(1:2,1:2,fr) = tform.T(1:2,1:2);
        this_t(3,1:2,fr) = this_t(3,1:2,fr) + tform.T(3,1:2);
    end
    expInfoObj.setAllT(expInfoObj.positions{p}, expInfoObj.channels{this_fluo}, this_t);
end
expInfoObj.saveExpInfo()
end

function addFluorescence(exp_dir,fluo_chan_name, bg_window)
%add fluorescence to the tracked cells. 
%exp_dir - the 'prod' directory
%fluo_chan_name- name of fluorescence channel
%bg_window - if not empty will be used to subract background calculation

eIO = load(fullfile(exp_dir,'expInfoObj.mat'));
expInfoObj = eIO.expInfoObj;
for pi = 1:length(expInfoObj.positions)
    pos_name = expInfoObj.positions{pi};
    %using sigma parameter from Spartak addition
    if ~isempty(bg_window)
        computeCellFluorescence(expInfoObj, pos_name, fluo_chan_name, bg_window);
    else
        computeCellFluorescence(expInfoObj, pos_name, fluo_chan_name);
    end
end
end

function copyPipeline(anal_dir,dirs)
%function to copy image processing files into analysis directory
for d = dirs
    ending = strsplit(d{1},'/');
    ending = ending{end};
    status = copyfile(d{1}, fullfile(anal_dir,ending));
    if status == 0
        disp('the analysis folder couldnt be copied')
    end
end
end
