# Auto_C3D_Checker
Automatic .c3d evaluation and modification pipeline for neuromusculoskeletal modeling.

# Read the User Guide.docx to start using this tool.

## What is this?
This is a tool for automatically evaluating the quality of a biomechanical trial (in .c3d format) involving the use of marker-driven motion capture, forceplates (from here on referred to as FP), and surface EMG.
This tool evaluates the .c3d in terms of its appropriateness for neuromusculoskeletal (from here on referred to as NMS) modelling. Specifically, this tool automatically determines:

1.	Whether Markers crucial for tracking are present throughout the trial
  - If yes, what is the motion direction throughout the trial

2.	If the instrumented leg hits a FP during the trial 
  - If yes, which FP is hit and does the foot/both feet hit this FP fully or only partially
  - If FP is hit, at what frame do foot strike and foot off happen
  - If all above are found, the required trial length (foot strike to foot off, plus extra frames for analysis window and EMG delay) is calculated 

3.	If the recorded trial is long enough for the required trial length
  - If not, is there sufficient data to account for EMG delay while padding the trial
  - If yes does start/end/both need padding
  - If required, padding is then performed

4.	Whether EMG signal quality is sufficient for EMG-assisted NMS modelling
  - If not, how many of the total EMG in the trial are of insufficient quality

Based on these .c3d quality criteria, trials are modified (i.e., start padding) and divided into three groups:

*	Calibration – Highest quality trials (where the instrumented leg fully hits a FP and trials have required length) without a need for any data modifications and with at least 75% of EMG being of sufficient quality. These trials can be used to calibrate the NMS model.

*	Execution – Trials (where the instrumented leg fully hits a FP) with a need for data modifications (i.e., padding) or with insufficient EMG quality. These trials can be used to run the NMS model, but not to calibrate it.

*	Unusable – Trials where the instrumented leg does not hit the FP (or not fully), or where not enough data is present to account for EMG delay for padding. 

The checker outputs a “Results.xlsx” containing the results of all the abovementioned checks for all participants. Further, it creates an “EMG Figures” folder for each participant, where the classified EMG can be found. Furthermore, it creates an “InputData” folder and places the modified Calibration/Execution .c3d trials into appropriate folders.
