%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Auto Elaboration XML Maker                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Riad Akhundov <r.akhundov@griffith.edu.au> and/or <Riad.Akhundov@uon.edu.au>

%Bonus script for automatic elaboration.xml creation for each subject & session. 
%Part of the Auto_C3D_Checker toolbox.

%%%Requirements: 
%1) MATLAB 2019b or newer (made with MATLAB version 2023a)
%2) Output from the Auto_C3D_Checker for the same dataset

%Version: v0.23.08.25

%%%ToDo:
% *) Share a sunset


clc; clearvars; close all;
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
warning('off', 'MATLAB:MKDIR:DirectoryExists')
%% Manual Inputs
disp('%%% Script started %%%');disp('%')

%Paths
inputDataPath = [pwd, '\..\Sample Data\InputData'];
elaboratedDataPath = [pwd, '\..\Sample Data\ElaboratedData'];
resultsExcelPath = dir([inputDataPath, '\Results', '*.xlsx']);
resultsExcelPath = [resultsExcelPath.folder, '\', resultsExcelPath.name];
templateXMLPath = [pwd, '\..\templatesXML\autoElaborationTemplate_example_'];
autoC3DxmlTemplate = [pwd, '\..\templatesXML\autoC3Dsetup_example.xml'];

%Filter cutoff frequencies
fcutMarkers = 18;
fcutFP = 18;


%% 2) Make acquisition.xml From The Above Manual Inputs & Auto C3D Checker Results Excel
%Subject ID & emgSide (L or R) taken from autoC3Dsetup
prefXMLRead.Str2Num = 'never';
prefXmlWrite.StructItem = false;
prefXmlWrite.CellItem   = false;

autoC3DxmlTree = xml_read(autoC3DxmlTemplate);
autoC3DxmlMarkers = autoC3DxmlTree.MarkersProtocol.MarkersSetDynamicTrials;
autoC3DxmlEMGDelay = autoC3DxmlTree.EMGs.EMGDelay;
autoC3DxmlWindow = autoC3DxmlTree.ScriptSettings.analysisWindow;

subjects = split(autoC3DxmlTree.Subjects.SubjectCodes);
subjectEMG = split(autoC3DxmlTree.Subjects.InstrumentedLeg);
nSubjects = length(subjects);

for s=1:nSubjects
    disp(['%% ', subjects{s}, ' %%']);
    currentDirInput = [elaboratedDataPath, '\', subjects{s}, '\'];
    
    %Find all sessions in this participant folder
    currentDirInput_Sub = dir(currentDirInput);
    currentDirInput_Sub(ismember({currentDirInput_Sub.name}, {'.', '..'})) = [];
    currentDirInput_Sub = currentDirInput_Sub([currentDirInput_Sub.isdir]);
    nSubjectSessions = length(currentDirInput_Sub);    

    %Read in results excel sheet for this participant
    resultsExcelTable = readtable(resultsExcelPath,'Sheet', subjects{s});
    resultsExcelTable(end-3:end,:) = [];
    resultsExcelColumns = fieldnames(resultsExcelTable);
 
    for i = 1:nSubjectSessions
        subjSession = currentDirInput_Sub(i).name;
        currentDirInput_Session = [currentDirInput, subjSession, '\sessionData\'];
        currentDirOutput_Session =  [currentDirInput, subjSession, '\dynamicElaborations\'];
        mkdir(currentDirOutput_Session);

        currentXMLTree = xml_read([templateXMLPath, subjectEMG{s}, '.xml'], prefXMLRead);

        load([currentDirInput_Session, 'trialsName.mat']);

        trialsDynamic = trialsName(~contains(upper(trialsName),'STATIC'));

        %Set FolderName
        currentXMLTree.FolderName=['.\InputData\', subjects{s}, '\', subjSession];

        %Set trial names
        trialsDynamicNames = strjoin(trialsDynamic);
        currentXMLTree.Trials = trialsDynamicNames;

        %Set filtering parameters and analysis windows
        for fl = 1:length(trialsDynamic)
            currentXMLTree.Filtering.Trial(1,fl).Name = trialsDynamic{1,fl};
            currentXMLTree.Filtering.Trial(1,fl).Fcut.Markers = fcutMarkers;
            currentXMLTree.Filtering.Trial(1,fl).Fcut.Forces = fcutFP;

            currentXMLTree.WindowSelectionProcedure.Manual.TrialWindow(1,fl).TrialName = trialsDynamic{1,fl};

            %Find chosen FP for trial and get start stop frames for it
            startFrameColumn = resultsExcelColumns{find(strcmp(resultsExcelColumns, ['StanceOnFP', num2str(resultsExcelTable.ChosenFP(strcmp(resultsExcelTable.Trials,...
                trialsDynamic{1,fl})))]))+1,1}; %Find the columnname of startFrame for chosen FP

            endFrameColumn = resultsExcelColumns{find(strcmp(resultsExcelColumns, ['StanceOnFP', num2str(resultsExcelTable.ChosenFP(strcmp(resultsExcelTable.Trials,...
                trialsDynamic{1,fl})))]))+2,1}; %Find the columnname of endFrame for chosen FP

            %Adjust start/stop by analysis window offset
            currentStartFrame = resultsExcelTable.(startFrameColumn)(strcmp(resultsExcelTable.Trials, trialsDynamic{1,fl})) - autoC3DxmlWindow;
            currentEndFrame = resultsExcelTable.(endFrameColumn)(strcmp(resultsExcelTable.Trials, trialsDynamic{1,fl})) + autoC3DxmlWindow;

            %Set analysis window
            currentXMLTree.WindowSelectionProcedure.Manual.TrialWindow(1,fl).StartFrame = currentStartFrame;
            currentXMLTree.WindowSelectionProcedure.Manual.TrialWindow(1,fl).EndFrame = currentEndFrame;  
        end

        %Set Markers
        currentXMLTree.Markers = autoC3DxmlMarkers;
        
        %Set EMGOffset
        currentXMLTree.EMGOffset = autoC3DxmlEMGDelay;
        
        %Write xml
        xml_write([currentDirOutput_Session 'elaboration.xml'],currentXMLTree,'elaboration',prefXmlWrite);
        clear trialsDynamic

        disp(['% elaboration.xml for ', subjSession, ' done %']);
    end %Sessions
end %Subjects

disp('%'); disp('%%% Script finished successfully %%%');
