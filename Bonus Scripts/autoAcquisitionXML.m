%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Auto Acquisition XML Maker                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Riad Akhundov <r.akhundov@griffith.edu.au> and/or <Riad.Akhundov@uon.edu.au>

%Bonus script for automatic acquisition.xml creation for each subject & session. 
%Part of the Auto_C3D_Checker toolbox.

%%%Requirements: 
%1) MATLAB 2019b or newer (made with MATLAB version 2023a)
%2) Output from the Auto_C3D_Checker for the same dataset

%Version: v0.23.08.25

%%%ToDo:
% *) Dance under the stars


clc; clearvars; close all;
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
%% Manual Inputs
disp('%%% Script started %%%');disp('%')

%Paths
inputDataPath = [pwd, '\..\Sample Data\InputData'];
resultsExcelPath = dir([inputDataPath, '\Results', '*.xlsx']);
resultsExcelPath = [resultsExcelPath.folder, '\', resultsExcelPath.name];
templateXMLPath = [pwd, '\..\templatesXML\autoAcquisitionTemplate_example.xml'];
autoC3DxmlTemplate = [pwd, '\..\templatesXML\autoC3Dsetup_example.xml'];

%Height in m and mass in kg (needs to be in identical order to subjectEmg)
heightMass=[...
1.685, 62.5; %ACLR11
1.685, 75.0; %ACLR28
1.682, 48.5; %CTRL09
    ]; 


%% 2) Make acquisition.xml From The Above Manual Inputs & Auto C3D Checker Results Excel
%Subject ID & emgSide (L or R) taken from autoC3Dsetup
prefXMLRead.Str2Num = 'never';
treeAutoC3D = xml_read(autoC3DxmlTemplate, prefXMLRead);
subjects = split(treeAutoC3D.Subjects.SubjectCodes);
subjectEMG = split(treeAutoC3D.Subjects.InstrumentedLeg);
nSubjects = length(subjects);

templateXMLTree = xml_read(templateXMLPath, prefXMLRead);

for s=1:nSubjects
    disp(['%% ', subjects{s}, ' %%']);
    currentDirInput = [inputDataPath, '\', subjects{s}, '\'];
    
    %Find all sessions in this participant folder
    currentDirInput_Sub = dir(currentDirInput);
    currentDirInput_Sub(ismember({currentDirInput_Sub.name}, {'.', '..', 'EMG Figures'})) = [];
    currentDirInput_Sub = currentDirInput_Sub([currentDirInput_Sub.isdir]);
    nSubjectSessions = length(currentDirInput_Sub);

    %Read in results excel sheet for this participant
    resultsExcelTable = readtable(resultsExcelPath,'Sheet', subjects{s});
    resultsExcelTable(end-3:end,:) = [];

    for i = 1:nSubjectSessions
        subjSession = currentDirInput_Sub(i).name;
        currentDirInput_Session = [currentDirInput, subjSession, '\'];

        %Adjust anthropometric values
        currentXMLTree = templateXMLTree;
        currentXMLTree.Subject.FirstName = subjects{s};
        currentXMLTree.Subject.Code = subjects{s};
        currentXMLTree.Subject.Weight = heightMass(s,2);
        currentXMLTree.Subject.Height = heightMass(s,1);
        if strcmp(subjectEMG{s}, 'L')
            currentXMLTree.EMGs.Protocol.InstrumentedLeg = 'Left';
        elseif strcmp(subjectEMG{s}, 'R')
            currentXMLTree.EMGs.Protocol.InstrumentedLeg = 'Right';
        else
            error(['ERROR: Specify EMG leg for ' subjects{s}]);
        end

        %Find all trials in session
        trialNames = dir([currentDirInput_Session, '*.c3d']);
        trialNames = {trialNames.name}';
        
        trialsDynamic = trialNames(~contains(upper(trialNames),'STATIC'));
        trialStatic = trialNames(contains(upper(trialNames),'STATIC'));
        
        if i ==1 && isempty(trialStatic) %Need a static trial at least in each 1Calibration folder (I prefer to put one in all session folders)
            error(['ERROR: No static trial for ', currentDirInput_Session]);    
        end

        for j = 1:length(trialsDynamic)+1
            if j == 1 %Static trial is always first for consistency   
                %Repetition and type
                trialType = erase(trialStatic{1}, '.c3d');

                currentXMLTree.Trials.Trial(j).Type = trialType; 
                currentXMLTree.Trials.Trial(j).RepetitionNumber = [];

                %MotionDirection                
                currentXMLTree.Trials.Trial(j).MotionDirection = char(resultsExcelTable.MotionDirection(strcmp(resultsExcelTable.Trials,...
                    trialStatic{1}(1:end-4))));

                %StanceOnFP
                for fp = 1:length(currentXMLTree.Trials.Trial(j).StancesOnForcePlatforms.StanceOnFP)
                    currentXMLTree.Trials.Trial(j).StancesOnForcePlatforms.StanceOnFP(fp).Forceplatform = fp;
                    
                    currentStance = char(resultsExcelTable.(['StanceOnFP', num2str(fp)])(strcmp(resultsExcelTable.Trials,...
                        trialStatic{1}(1:end-4))));
                    if isempty(currentStance)
                        currentXMLTree.Trials.Trial(j).StancesOnForcePlatforms.StanceOnFP(fp).Leg = 'None';
                    else
                        currentXMLTree.Trials.Trial(j).StancesOnForcePlatforms.StanceOnFP(fp).Leg = currentStance;
                    end                  
                end %fp

            else
                %Repetition and type
                trialRepetition = regexp(trialsDynamic{j-1}(end-6:end),'\d*','Match');
                trialRepetition = trialRepetition{1};
    
                trialType = erase(trialsDynamic{j-1}, [trialRepetition,'.c3d']);

                currentXMLTree.Trials.Trial(j).Type = trialType; 
                currentXMLTree.Trials.Trial(j).RepetitionNumber = trialRepetition;

                %MotionDirection               
                currentXMLTree.Trials.Trial(j).MotionDirection = char(resultsExcelTable.MotionDirection(strcmp(resultsExcelTable.Trials,...
                    trialsDynamic{j-1}(1:end-4))));

                %StanceOnFP
                for fp = 1:length(currentXMLTree.Trials.Trial(1).StancesOnForcePlatforms.StanceOnFP)
                    currentXMLTree.Trials.Trial(j).StancesOnForcePlatforms.StanceOnFP(fp).Forceplatform = fp;

                    currentStance = char(resultsExcelTable.(['StanceOnFP', num2str(fp)])(strcmp(resultsExcelTable.Trials,...
                        trialsDynamic{j-1}(1:end-4))));
                    if isempty(currentStance)
                        currentXMLTree.Trials.Trial(j).StancesOnForcePlatforms.StanceOnFP(fp).Leg = 'None';
                    else
                        currentXMLTree.Trials.Trial(j).StancesOnForcePlatforms.StanceOnFP(fp).Leg = currentStance;
                    end  
                end %fp
            end %Static/Dynamic
        end %Trials

        %Write xml
        prefXmlWrite.StructItem = false;
        prefXmlWrite.CellItem   = false;
        xml_write([currentDirInput_Session, 'acquisition.xml'], currentXMLTree, 'acquisition', prefXmlWrite);
        clear trialsDynamic trialStatic

        disp(['% acquisition.xml for ', subjSession, ' done %']);
    end %Sessions      
end %Subjects

disp('%'); disp('%%% Script finished successfully %%%');
