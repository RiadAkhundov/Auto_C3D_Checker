%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Automatic C3D Checker                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Riad Akhundov <r.akhundov@griffith.edu.au> and/or <Riad.Akhundov@uon.edu.au>

%Main script for... 

%%%Requirements: 
%1) MATLAB 2019b or newer (made with MATLAB version 2023a)
%2) Deep Learning Toolbox (introduced in version 2018a)
%3) Parallel Computing Toolbox (also called Distributed Computing Toolbox)

%Version: v0.23.08.14

%%%ToDo:
%!!! include naming of fp channels in xml !!!
%!!! make user guide !!!
%!!! make GUI !!!

%%%ToDo Later:
% Find MotionDirection through foot orientation through majority of trial
% Combine with EMG class to excitation script 
% Add Manual ReNaming script



clc; clearvars; close all;
tic
%% 0) Initial Setup
warning('off', 'MATLAB:xlswrite:AddSheet')
warning('off', 'MATLAB:MKDIR:DirectoryExists')
disp('%%% Script started %%%');disp('%')

%Change Current Folder to m-file location (if "Add to Path" used instead of "Change Folder")
if ~isdeployed
  cd(fileparts(which(mfilename)));
end

%Create diary (command window logger)
currentDateTime = datestr(now,'dd-mm-yyyy-HH-MM');
diaryName = ['log_', currentDateTime]; 
diary(diaryName); 
diary on

%User inputs - Paths to base folder, template .xml, etc.
baseFolderPath = [pwd, '\Sample Data\Base Folder']; %Base folder containing subjects subfolders with .c3d files
xmlTemplate = [pwd, '\templatesXML\acquisition_example.xml']; %.xml containing acquisition information
overwriteEMGNames = true;
runEMGClass = true;
%!!! add user inputable paths with GUI

%Add paths to btk & functions
addpath(genpath([pwd, '\btk & functions'])); %Path to btk and other functions used in this script
addpath(genpath([pwd, '\templatesXML'])); %Path to template .xml used in this script

%Get lab orientation, subjects, instumented legs, important markers, foot markers, and EMG info from .xml
prefXMLRead.Str2Num = 'never';
tree = xml_read(xmlTemplate, prefXMLRead);
% -> add user inputable path to template .xml (+same with GUI)

% labOrientation = tree.Laboratory.CoordinateSystemOrientation; %Not used yet
subjectInstrumented(:,1) = split(tree.Subjects.SubjectCodes);
subjectInstrumented(:,2) = split(tree.Subjects.InstrumentedLeg);

% importantMarkersStatic = split(tree.MarkersProtocol.MarkersSetStaticTrials); %Not used yet
importantMarkersDynamic = split(tree.MarkersProtocol.MarkersSetDynamicTrials);

footMarkersRight = split(tree.MarkersProtocol.RightFootMarkers);
footMarkersLeft = split(tree.MarkersProtocol.LeftFootMarkers);

emgNamesOriginal = split(tree.EMGs.OriginalChannels);
emgNamesAdjusted = split(tree.EMGs.RenamedChannels);
emgDelay = str2double(tree.EMGs.EMGDelay); %0.2 sec is the MOtoNMS default value for electromechanical + hardware delay

thresholdFP = str2double(tree.ScriptSettings.thresholdFP);
analysisWindow = str2double(tree.ScriptSettings.analysisWindow);
padNum = str2double(tree.ScriptSettings.paddingFrames); %How many frames to pad .c3d with

excelTemplateFilePath = [pwd, '\template.xlsx']; %Excel results prep
excelFilePath = [baseFolderPath, '\Results.xlsx'];
copyfile(excelTemplateFilePath, excelFilePath);
excelEmptyArray = cell(1000,100);

if runEMGClass; load('combo_xception.mat'); disp('%% Loaded neural network %%'); end


%% 1) Loop Through Subjects & Trials 
subjectFolders = dir(baseFolderPath); subjectFolders(1:2) = []; %Get only the subject subfolders in baseFolderPath
subjectFolders = subjectFolders([subjectFolders.isdir]); %Delete any non-folders

for s = 1:length(subjectFolders)
    disp(['%% ', subjectFolders(s).name, ' %%']);
    c3dChosenFilesPath = [subjectFolders(s).folder, '\', subjectFolders(s).name, '\InputData\'];
    calibrationFilesPath = [c3dChosenFilesPath, subjectFolders(s).name, '\1Calibration\'];
    executionFilesPath = [c3dChosenFilesPath, subjectFolders(s).name, '\2Execution\'];
    mkdir(c3dChosenFilesPath); mkdir(calibrationFilesPath); mkdir(executionFilesPath);

    idxSubject = contains(subjectInstrumented(:,1),subjectFolders(s).name); %Find subject idx in subjectInstrumented list
    instrumentedSide = subjectInstrumented{idxSubject,2};
    countCalibration = 0;
    countExecution = 0;
    countUnusable = 0;

    %Check if instrumentedSide is defined for the subject
    if  strcmp(instrumentedSide,'R')
        instrumentedLeg = 'Right';      
    elseif strcmp(instrumentedSide,'L')
        instrumentedLeg = 'Left';        
    else
        error(['%%% No EMG side defined for ', subjectFolders(s).name, ' ! %%%'])
    end

    if runEMGClass
        dirOutput_Figures = [subjectFolders(s).folder, '\', subjectFolders(s).name, '\EMG Figures']; %Make Figure folder at script location
        mkdir(dirOutput_Figures);
        dirOutput_Figures_Good = [dirOutput_Figures, '\1_Good']; mkdir(dirOutput_Figures_Good);
        dirOutput_Figures_Noisy = [dirOutput_Figures, '\2_Noisy']; mkdir(dirOutput_Figures_Noisy);
        dirOutput_Figures_Bad = [dirOutput_Figures, '\3_Bad']; mkdir(dirOutput_Figures_Bad);
        dirOutput_Figures_Missing = [dirOutput_Figures, '\4_Missing']; mkdir(dirOutput_Figures_Missing);
    end
    
    currentTrials = dir([subjectFolders(s).folder, '\', subjectFolders(s).name, '\*.c3d']); %Get subject trials
    classifications = cell(length(currentTrials), length(emgNamesOriginal)); %Preallocating classifications cell array (one of the few instances where preallocation is actually useful/required) 

    for t = 1:length(currentTrials)       
        h = btkReadAcquisition([currentTrials(t).folder, '\', currentTrials(t).name]);
        disp(['% ', currentTrials(t).name, ' %']);
%         h = btkReadAcquisition('E:\Melbourne ACLR Dataset\_Sorted Gait Data\__Outlier examples\FP1 Artefact - CTRL04 Forward hop 95.c3d');
              
        emgNamesOriginal = split(tree.EMGs.OriginalChannels); %Calling it again, but it's needed if overwriteEMGNames is enabled, don't see a better solution

        originalFirstFrame = btkGetFirstFrame(h);
        btkSetFirstFrame(h, 1); %Set first frame of .c3d to be 1 to reduce complexity and avoid bugs
        lastFrame = btkGetLastFrame(h);
        
        %Marker data
        [markers, markersInfo, ~] = btkGetMarkers(h);
        markersRate = markersInfo.frequency;  
        requiredExtraFrames = emgDelay*markersRate+analysisWindow+1; 

        %Forceplate data
        [forceplates, ~] = btkGetForcePlatforms(h); %Read original FP data

        %EMG data
        [analogs, analogsInfo] = btkGetAnalogs(h);
        analogRate = analogsInfo.frequency; 
    
        %Adjust EMG names in .c3d 
        if overwriteEMGNames    
            for e = 1:length(emgNamesOriginal)
                if isfield(analogs,(emgNamesOriginal{e})) == 1
                    [analogs, ~] = btkSetAnalogLabel(h, emgNamesOriginal{e}, emgNamesAdjusted{e});
                end
            end
            emgNamesOriginal = emgNamesAdjusted;
        end

        fp2markerRate = btkGetAnalogSampleNumberPerFrame(h);
        numFP = length(forceplates);
        paddingRequiredAtStart = zeros(1,numFP);
        paddingRequiredAtEnd = zeros(1,numFP);       
        startPadded = 'No';
        endPadded = 'No';
        fpReNameCheck = 0;

        
        %% 2) Find Foot On FP
        for fp = 1:numFP
            res = 4+(fp-1)*6; %Results cell array idx helper
            fpStr = num2str(fp);
            fpFieldNames = fieldnames(forceplates(fp).channels);
            runTrialPaddingCheck = false;
            
            %Check/rename FP channels in .c3d
            [h, fpReNameCheck(fp), currentFPData] = fpChannelCheck(h, fp, fpStr, fpFieldNames, forceplates);
 
            %Change FP channel labels if they are not correctly named in .c3d
            if fpReNameCheck(fp)
                [forceplates, ~] = btkGetForcePlatforms(h);
                fpFieldNames = fieldnames(forceplates(fp).channels);
            end
            
            %Filter/adjust data before calculating tStrike_fpFrames and tOff                      
            filtFPData = movmean(currentFPData, 10); %Moving average filter | moving window of 10 works well
            filtFPData = filtFPData + prctile(abs(filtFPData),1); %Remove zero offset | prctile works better than min
              
            %%% figure; 
            %%% hold on;
            %%% plot(abs(currentFPData));
            %%% plot(abs(filtFPData));
            %%% plot(abs(filtFPData - prctile(abs(filtFPData),1)));

            %Check FP for footstrike and if footstrike is present get the first frame
            tStrike_fpFrames{fp} = (find(abs(filtFPData)>thresholdFP,1)); %Footstrike on FP in FP frame rate
            
            %Find first foot off frame on FP after foot strike (won't find multiple footstrikes per FP!)          
            tOff_fpFrames{fp} = (find(abs(filtFPData(tStrike_fpFrames{fp}:end,1))<thresholdFP,1)+tStrike_fpFrames{fp});

            timeOnFP{fp} = ceil((tOff_fpFrames{fp} - tStrike_fpFrames{fp})/fp2markerRate)/markersRate;    
            
            if ~isempty(tStrike_fpFrames{fp}) && isempty(timeOnFP{fp}) || ~isempty(tStrike_fpFrames{fp}) && timeOnFP{fp} > 10/markersRate  
                tStrike{fp} = ceil(tStrike_fpFrames{fp}/fp2markerRate); %Footstrike on FP in marker frames              
                fpOrigin{fp} = (sum(forceplates(fp).corners,2)./4)'; %Calculate correct FP origin               
                             
                framesBefore_tStrike{fp} = tStrike{fp}-1; %Frames before footstrike
                
                %Frames after footoff (if not on FP at the end of trial)
                if isempty(tOff_fpFrames{fp})
                    tOff{fp} = [];
                    framesAfter_tOff{fp} = 'On FP';
                    framesOnFP{fp} = lastFrame-tStrike{fp}-1;
                else
                    tOff{fp} = ceil(tOff_fpFrames{fp}/fp2markerRate); %Footoff from FP in marker frames
                    framesAfter_tOff{fp} = lastFrame-ceil(tOff_fpFrames{fp}/fp2markerRate)-1;
                end
                
                %Find foot on FP
                for m = 1:length(footMarkersRight) 
                    %Calculate distance of foot markers at tStrike
                    currentMarkerDist_R(m,1) = norm(fpOrigin{fp} - markers.(footMarkersRight{m})(tStrike{fp},:));
                    currentMarkerDist_L(m,1) = norm(fpOrigin{fp} - markers.(footMarkersLeft{m})(tStrike{fp},:));  
                    
                    %Calculate if foot marker is on FP
                    fullFootOnFP_R(m,1) = inpolygon(markers.(footMarkersRight{m})(tStrike{fp},1),...
                        markers.(footMarkersRight{m})(tStrike{fp},2),forceplates(fp).corners(1,:)',forceplates(fp).corners(2,:)');
                    fullFootOnFP_L(m,1) = inpolygon(markers.(footMarkersLeft{m})(tStrike{fp},1),...
                        markers.(footMarkersLeft{m})(tStrike{fp},2),forceplates(fp).corners(1,:)',forceplates(fp).corners(2,:)');
                end
                
                markerDist_R{fp} = mean(currentMarkerDist_R); %Average distance of footMarkersRight from origin at tStrike          
                markerDist_L{fp} = mean(currentMarkerDist_L); 
                
                fullFootOnFP_R = all(fullFootOnFP_R); %Is full foot on FP 
                fullFootOnFP_L = all(fullFootOnFP_L);
                
                %% 3) Compile FP Results
                %Leg on FP & full foot on FP
                if abs(markerDist_L{fp}-markerDist_R{fp}) < 200 
                    results{t,res} = 'Both';
                    if fullFootOnFP_R && fullFootOnFP_L
                        results{t,res+5} = 'Yes';
                    else
                        results{t,res+5} = 'No';
                    end
                elseif markerDist_L{fp} > markerDist_R{fp}
                    results{t,res} = 'Right';
                    if fullFootOnFP_R
                        results{t,res+5} = 'Yes';                        
                        runTrialPaddingCheck = true;          
                    else
                        results{t,res+5} = 'No';
                    end
                elseif markerDist_L{fp} < markerDist_R{fp} 
                    results{t,res} = 'Left';
                    if fullFootOnFP_L
                        results{t,res+5} = 'Yes';                        
                        runTrialPaddingCheck = true;
                    else
                        results{t,res+5} = 'No';
                    end
                else
                    results{t,res} = 'Unknown';
                end
                                
                %Trial padding checks
                if runTrialPaddingCheck
                    %If instrumented leg fully hits this fp, check if trial needs padding
                    if strcmp(results{t,res}, instrumentedLeg) 
                        %Note that instrumented leg fully hits this fp (for cases where instrumented leg fully hits multiple fp during trial)
                        paddingRequiredAtStart(fp) = 1;  
                        paddingRequiredAtEnd(fp) = 1;

                        %Start padding
                        if framesBefore_tStrike{fp} < requiredExtraFrames 
                            if framesBefore_tStrike{fp} < 3
                                paddingRequiredAtStart(fp) = 3; %No padding possible, instrumented leg on FP at trial start
                                startPadded = 'Error';
                            else
                                paddingRequiredAtStart(fp) = 2; %FP needs start padding
                            end
                        end

                        %End padding
                        if framesAfter_tOff{fp} < requiredExtraFrames
                            if framesAfter_tOff{fp} < (requiredExtraFrames - analysisWindow)
                                paddingRequiredAtEnd(fp) = 3; %Not enough EMG data for padding
                                endPadded = 'Error';
                            else
                                paddingRequiredAtEnd(fp) = 2; %FP needs end padding
                            end
                            
                        %Special case for balance/hop trials    
                        elseif strcmp(framesAfter_tOff{fp},'On FP')                                                           
                            if framesOnFP{fp} < (requiredExtraFrames*5)                               
                                paddingRequiredAtEnd(fp) = 3; %Not enough EMG data for padding
                                endPadded = 'Error'; 
                            end 
                        end %End padding 
                    end
                end %Padding check

                results{t,res+1} = tStrike{fp}; %tStrike frame
                results{t,res+2} = tOff{fp}; %tOff frame
                results{t,res+3} = framesBefore_tStrike{fp}; %Frames before tStrike
                results{t,res+4} = framesAfter_tOff{fp}; %Frames after tOff
                
            elseif ~isempty(tStrike_fpFrames{fp}) && timeOnFP{fp} <= 10/markersRate
                results{t,res} = 'Artefact';               
            else
                results{t,res} = 'None';
            end %tStrike 
        end %numFP
        
        results{t,1} = currentTrials(t).name(1:end-4); %Trial name
        results{t,2} = instrumentedLeg; %Instrumented leg
        
        %Check if instrumented leg hits any fp
        if any(paddingRequiredAtStart)
            results{t,3} = 'Yes';            
        else
            results{t,3} = 'No';
        end        
        
        %Find first FP not requiring start and/or end padding
        chosenFP{t,1} = find(paddingRequiredAtStart == 1 & paddingRequiredAtEnd == 1, 1);
                
        if isempty(chosenFP{t,1})            
            if any(paddingRequiredAtStart == 2) %If start padding is required
                chosenFP{t,1} = find(paddingRequiredAtStart == 2, 1);
                [h, startPadded] = paddingAtStart(h,importantMarkersDynamic, fp2markerRate, padNum, tStrike{chosenFP{t,1}});
                
                res = 4+(chosenFP{t,1}-1)*6; %Results cell array idx helper for chosenFP
                results{t,res+1} = results{t,res+1} + padNum; %Account for padding in tStrike frame
                results{t,res+2} = results{t,res+2} + padNum;
                results{t,res+3} = results{t,res+3} + padNum;
                
                if strcmp(startPadded, 'Error'); chosenFP{t,1} = []; end %Don't use trial if paddingAtStart produced an error 
            end
            
            if (isempty(chosenFP{t,1}) && any(paddingRequiredAtEnd == 2)) ||... %No need to pad end for different FP if start already padded
                    ~isempty(find(paddingRequiredAtStart == 2 & paddingRequiredAtEnd == 2, 1)) %If both start and end padding for same FP is required
                chosenFP{t,1} = find(paddingRequiredAtEnd == 2, 1);
                lastFrame_padded = lastFrame + padNum;
                btkSetFrameNumber(h, lastFrame_padded); %Set new last frame, thus padding the end of the trial
                endPadded = 'Yes'; 
                
                res = 4+(chosenFP{t,1}-1)*6; %Results cell array idx helper for chosenFP
                results{t,res+4} = results{t,res+4} + padNum; %Account for padding in framesAfter_tOff
            end 
        else
            startPadded = 'No'; %If FP not requiring start and/or end padding exists
            endPadded = 'No';
        end
        
        chosenFP{t,2} = startPadded;
        chosenFP{t,3} = endPadded;

        %Trial usability for FP
        if any(chosenFP{t,1}) 
            if strcmp(chosenFP{t,2},'No') && strcmp(chosenFP{t,3},'No')
                chosenFP{t,4} = 'Calibration'; %Full foot on fp, instrumented leg, no padding required
                countCalibration = countCalibration + 1;
                btkWriteAcquisition(h, [calibrationFilesPath, currentTrials(t).name]); %Write calibration quality .c3d file 
            elseif strcmp(chosenFP{t,2},'Error') || strcmp(chosenFP{t,3},'Error')
                chosenFP{t,4} = 'Unusable'; %Error in padding
                countUnusable = countUnusable +1;
            else
                chosenFP{t,4} = 'Execution'; %Full foot on fp, instrumented leg, but padding required
                countExecution = countExecution +1;
                btkWriteAcquisition(h, [executionFilesPath, currentTrials(t).name]); %Write execution quality .c3d file
            end
        else
            chosenFP{t,4} = 'Unusable'; %Not instrumented leg
            countUnusable = countUnusable +1;
        end
                

        %% 4) Make EMG Figures For Calibration/Execution Trials & Classify EMG             
        if runEMGClass && strcmp(chosenFP{t,4},'Calibration') || runEMGClass && strcmp(chosenFP{t,4},'Execution')           
            for i = 1:length(emgNamesOriginal)
                if isfield(analogs, (emgNamesOriginal{i}))
                    if any(analogs.(emgNamesOriginal{i}))
                        %Data for images
                        currentEMG = analogs.(emgNamesOriginal{i}) - mean(analogs.(emgNamesOriginal{i})); %Remove zero offset
                        currentEMG = currentEMG ./ max(abs(currentEMG)); %Self-Normalize                
        
                        %EMG envelope
                        BPFiltEMG = ZeroLagButtFiltfilt((1/analogRate), [30,300], 2, 'bp', currentEMG);
                        RectBPFiltEMG = abs(BPFiltEMG);
                        LinEnvEMG = ZeroLagButtFiltfilt((1/analogRate), 18, 2, 'lp', RectBPFiltEMG); 
                                
                        %Power spectral density
                        n = length(currentEMG); %Number of samples
                        fft_range = (0:(n-1)/2)/(n/2)*(analogRate/2); %Frequency range of real signal 
                    
                        y = fft(currentEMG); %FFT calculation
                        fft_power = abs(y(1:length(fft_range))).^2/n;
                    
                        % fft_power_25 = obw(fft_power,fft_range,[],25);
                        fft_power_50 = obw(fft_power,fft_range,[],50);
                        fft_power_75 = obw(fft_power,fft_range,[],75);
                        fft_power_99 = obw(fft_power,fft_range); %Estimate 99% occupied bandwidth of the power spectral density
    
    
                        % Plots
                        %3x Combo Plot      
                        f = figure('Color', [1 1 1], 'Visible', 'off'); 
                        %%% f = figure('Color', [1 1 1]);
                        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2.99 2.99]);
                    
                        %P3           
                        subplot(4,1,4);
                        hold on; 
                        plot(fft_range,fft_power,'LineWidth',1,'Color',[0 0.4470 0.7410]);
                        x1 = xline(fft_power_50,'-','LineWidth',2,'Color',[0 0 0]);
                        x2 = xline(fft_power_75,'-','LineWidth',2,'Color',[0 0 0]);
                        x3 = xline(fft_power_99,'-','LineWidth',2,'Color',[0 0 0]);
                        axis tight;
                        xlim([0 analogRate/2]);
                        set(gca,'xtick',[],'ytick',[]);
                        set(gca, 'Box', 'off');
                        set(gca,'YColor','none','XColor','none')
                        set(gca,'LooseInset',get(gca,'TightInset'));
                        set(gca,'Position',[0 0 1 0.25]);
                    
                        %P1
                        subplot(4,1,[1 2]);
                        hold on; 
                        plot(0:1/analogRate:(length(currentEMG)-1)/analogRate, currentEMG,'LineWidth',1,'Color',[0 0.4470 0.7410]); 
                        plot(0:1/analogRate:(length(LinEnvEMG)-1)/analogRate,LinEnvEMG,'LineWidth',2,'Color',[0.8500 0.3250 0.0980]); 
                        axis tight;
                        ylim([-1 1]);
                        set(gca,'xtick',[],'ytick',[]);
                        set(gca, 'Box', 'off');
                        set(gca,'YColor','none','XColor','none')
                        set(gca,'LooseInset',get(gca,'TightInset'));
                        set(gca,'Position',[0 0.5 1 0.5]);
                    
                        %P2
                        subplot(4,1,3);
                        plot(0:1/analogRate:(length(currentEMG)-1)/analogRate, currentEMG,'LineWidth',1,'Color',[0 0.4470 0.7410]);
                        xlim([0 3]);
                        ylim([-1 1]);
                        set(gca,'xtick',[],'ytick',[]);
                        set(gca, 'Box', 'off');
                        set(gca,'YColor','none','XColor','none')
                        set(gca,'LooseInset',get(gca,'TightInset'));
                        set(gca,'Position',[0 0.25 1 0.25]);
        
        
                        % Save image
                        print ([dirOutput_Figures, '\', currentTrials(t).name(1:end-4), '_', emgNamesOriginal{i}, '.jpg'], '-djpeg', '-r100');
            
                        close(f)
                    else                        
                        %Plot missing data
                        f = figure('Color', [1 1 1], 'Visible', 'off'); 
                        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2.99 2.99]);
        
                        plot(0:1/analogRate:(length(analogs.(emgNamesOriginal{i}))-1)/analogRate, analogs.(emgNamesOriginal{i})); 
                        axis tight;
        
                        %Save image
                        print ([dirOutput_Figures_Missing, '\', currentTrials(t).name(1:end-4), '_', emgNamesOriginal{i}, '.jpg'], '-djpeg', '-r100');

                        close(f)
                    end %If EMG non-zero
                end %If EMG exists
            end %EMGs

            %Create image datastore for better memory management
            imds = imageDatastore(dirOutput_Figures, 'IncludeSubfolders', false);
            
            %Classify images
            [labels, ~] = classify(combo_xception, imds); %The classification line
            % [labels, percentages] = classify(combo_xception, imds);
            %%% percentages_round = floor(percentages*10000)/100;      
            
            % Saving Results by Sorting Images
            for la = 1:length(labels)
                if labels(la,1) == '1'
                    movefile(imds.Files{la}, dirOutput_Figures_Good);                  
                elseif labels(la,1) == '2'
                    movefile(imds.Files{la}, dirOutput_Figures_Noisy);                
                elseif labels(la,1) == '3'
                    movefile(imds.Files{la}, dirOutput_Figures_Bad);        
                end %Class
            end %Images
    
            %Save label cell array
            labels = double(labels');
            for em = 1:length(emgNamesOriginal)    
                classifications{t,em} = labels(1,em);
            end
    
            %Update trial usability with EMG results
            numUsableEMG = sum(labels<3)/length(labels)*100;     
            if numUsableEMG < 75 && strcmp(chosenFP{t,4},'Calibration')  
                chosenFP{t,4} = 'Execution';
            end
        end %runEMGClass
        
        btkCloseAcquisition(h); %Close trial   
    end %Trials


    %% 5) Export Results To Excel    
    %Create headers for each FP
    for fp_res = 1:numFP
        fpResultsHeader{1,fp_res} = {['Stance on FP', num2str(fp_res)],'Footstrike Frame', 'Footoff Frame', 'Frames Before Footstrike',...
            'Frames After Footoff', 'Full Foot On FP?'};
    end
    
    %Add all headers to results cell array   
    if runEMGClass
        results = [[{'Trials','Instrumented Leg', 'Instrumented Leg Hit FP?'}, [fpResultsHeader{:}], emgNamesOriginal',...
            {'Chosen FP','Start Padded?', 'End Padded?', 'Trial Usability'}]; [results, classifications, chosenFP]];
    else
        results = [[{'Trials','Instrumented Leg', 'Instrumented Leg Hit FP?'}, [fpResultsHeader{:}],...
            {'Chosen FP','Start Padded?', 'End Padded?', 'Trial Usability'}]; [results, chosenFP]];
    end

    %Tally trial usability
    numCalibration(s,1) = countCalibration;
    numExecution(s,1) = countExecution;
    numUnusable(s,1) = countUnusable;
    
    countCalibration = {[num2str(round(countCalibration/length(currentTrials)*100)), '% (', num2str(countCalibration), ') Calibration']};
    countExecution = {[num2str(round(countExecution/length(currentTrials)*100)), '% (', num2str(countExecution), ') Execution']};
    countUnusable = {[num2str(round(countUnusable/length(currentTrials)*100)), '% (', num2str(countUnusable), ') Unusable']};
    
    %Update results with trial usability tally
    trialUsability = cell(4,size(results,2));
    trialUsability(end-2:end,end) = [countCalibration; countExecution; countUnusable];
    results = [results; trialUsability];

    %Excel magic
    newSheetName = subjectFolders(s).name;

    if s == 1
        %Open Activex server for Excel
        excelServer = actxserver('Excel.Application');
        % Make the application invisible
        set(excelServer, 'Visible', 0);
        % Make excel not display alerts
        set(excelServer,'DisplayAlerts',0);
    end
    
    %Open the worksheet.
    excelFile = excelServer.Workbooks.Open(excelFilePath);

    excelServer.ActiveWorkbook.ActiveSheet.Copy([],  excelServer.ActiveWorkbook.ActiveSheet);
    excelServer.ActiveWorkbook.ActiveSheet.Name = newSheetName;

    if s == length(subjectFolders)
        excelFile.Sheets.Item(1).Delete;
        excelFile.Save;
        excelFile.Close;
        excelServer.Quit;
        delete(excelServer);   
    else
        excelFile.Save;
        excelFile.Close;      
    end    
    
    %Save results in excel
    xlswrite(excelFilePath,excelEmptyArray,newSheetName,'A2:CU1002'); %"Depricated", but takes half the time compared to writecell
    xlswrite(excelFilePath,results,newSheetName);


    %% 6) Move Usable Participant Data Into Unified InputData Folder
    movefile(c3dChosenFilesPath,fileparts(baseFolderPath));
    
    clear results chosenFP classifications imds labels
end %Subjects

movefile(excelFilePath,[fileparts(baseFolderPath),'\InputData\']);


%% 7) Close and Move Diary
userview = memory; mem5 = userview.MemUsedMATLAB/1024^3;
t1 = toc;
disp(['%% Automatic C3D Checker ran successfully in ' num2str(floor(t1/60)) ' minutes and ' num2str(rem(t1,60))...
    ' seconds, using ' num2str(mem5) 'GB of RAM' ' %%'])

diary off
try
    movefile(diaryName, [fileparts(baseFolderPath), '\', diaryName, '.txt']); %Creates diary in folder selected at the beginning
catch err
    disp('%%%  Error while moving diary. %%%');
    disp(err.message)
end


%% Functions
function [h, fpReNameCheck] = fpChannelRename(h,fpFieldNames,fpChannelNames)
    %Rename FP channels if necessary
    if ~all(strcmp(fpFieldNames,fpChannelNames))       
        for ch = 1:length(fpFieldNames)
            try
                btkSetAnalogLabel(h, fpFieldNames{ch}, fpChannelNames{ch});
            catch ERR
                disp(['%%%  Error in fpChannelRename: ', ERR.message, ' %%%']);
                continue;  %Jump to next iteration in case of error
            end
        end %Channels  
        fpReNameCheck = true;
    else
        fpReNameCheck = false;
    end 
end


function [h, fpReNameCheck, currentFPData] = fpChannelCheck(h, fp, fpStr, fpFieldNames, forceplates)
    switch forceplates(fp).type                
        case 1
            fpChannelNames = {['Fx', fpStr]; ['Fy', fpStr]; ['Fz', fpStr]; ['Px', fpStr]; ['Py', fpStr]; ['Pz', fpStr]}; %Expected names
            
            [h, fpReNameCheck] = fpChannelRename(h,fpFieldNames,fpChannelNames);
            
            currentFPData = forceplates(fp).channels.(fpFieldNames{3}); %Fz is always the third channel for FP type 1, 2, and 4 
            
        case {2, 4} %FP type 2 and 4 have the same channel names
            fpChannelNames = {['Fx', fpStr]; ['Fy', fpStr]; ['Fz', fpStr]; ['Mx', fpStr]; ['My', fpStr]; ['Mz', fpStr]};
            
            [h, fpReNameCheck] = fpChannelRename(h,fpFieldNames,fpChannelNames);
            
            currentFPData = forceplates(fp).channels.(fpFieldNames{3});
            
        case 3
            fpChannelNames = {['F', fpStr, 'x12']; ['F', fpStr, 'x34']; ['F', fpStr, 'y14']; ['F', fpStr, 'y23'];...
                ['F', fpStr, 'z1']; ['F', fpStr, 'z2']; ['F', fpStr, 'z3']; ['F', fpStr, 'z4']};
            
            [h, fpReNameCheck] = fpChannelRename(h,fpFieldNames,fpChannelNames);           
            
            currentFPData = (forceplates(fp).channels.(fpFieldNames{5})+forceplates(fp).channels.(fpFieldNames{6})+...
                forceplates(fp).channels.(fpFieldNames{7})+forceplates(fp).channels.(fpFieldNames{8})); %For FP type 3: Fz=Fz1+Fz2+Fz3+Fz4
            
        otherwise                   
            error(['%%% FP type is undefined for ', subjectFolders(s).name, ': ', currentTrials(t).name, '! %%%'])                
    end
end


function [h, startPadded] = paddingAtStart(h,importantMarkersDynamic, fp2markerRate, padNum, tStrike)
    %% First Good Frame
    [points, ~] = btkGetPoints(h);

    %Make marker matrix
    marker_matrix = points.(importantMarkersDynamic{1}); %The crucial markers have to be present from the start of the analysis window
    for m = 2:length(importantMarkersDynamic)   
        marker_matrix(:,end+1:end+3) = points.(importantMarkersDynamic{m});
    end

    firstRealFrame = find(~any(marker_matrix == 0, 2), 1); %Find first frame with all markers present

    if firstRealFrame > tStrike
        disp('%%%  Error: Important markers only present after tStrike, cannot use this trial. %%%');
        startPadded = 'Error';
        return;
    elseif isempty(firstRealFrame)
        disp('%%%  Error: Important markers not present in trial, cannot use this trial. %%%');
        startPadded = 'Error';
        return;
    else
        btkSetFirstFrame(h, firstRealFrame); %Set first frame of .c3d to firstRealFrame 
        btkSetFirstFrame(h, 1); %Label first frame as 1 for consistency
        startPadded = 'Yes';
    end

    lastFrame = btkGetLastFrame(h); %Get updated last frame
    lastFrame_padded = lastFrame + padNum;

    btkSetFrameNumber(h, lastFrame_padded); %Set new last frame padding the end of the trial

    %Get updated data from .c3d (data was cropped by btkSetFirstFrame and btkSetFrameNumber)
    [points, ~] = btkGetPoints(h);
    [analogs, ~] = btkGetAnalogs(h);
    
    fields_points = fieldnames(points);
    fields_analogs = fieldnames(analogs);

    old_lentgh_points = lastFrame;
    old_lentgh_analogs = lastFrame*fp2markerRate;

    %Pad marker coordinates with data from firstRealFrame
    for i=1:length(fields_points) 
        %Switch padding from end to start
        points_new.(fields_points{i}) = [points.(fields_points{i})(old_lentgh_points+1:end,:) ; points.(fields_points{i})(1:old_lentgh_points,:)];   
        points_new.(fields_points{i})(1:padNum,:) = ones(padNum,1)*points_new.(fields_points{i})(padNum+1,:); %Replace padded values with first real values  
        btkSetPointValues(h, i, points_new.(fields_points{i}));
    end
    %Pad FP and EMG with zeroes
    for i=1:length(fields_analogs)        
        analogs_new.(fields_analogs{i}) = [analogs.(fields_analogs{i})(old_lentgh_analogs+1:end,:) ; analogs.(fields_analogs{i})(1:old_lentgh_analogs,:)];   
        btkSetAnalogValues(h, i, analogs_new.(fields_analogs{i}));
    end
end
