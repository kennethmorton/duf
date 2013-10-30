%%
clear all;
close all;
clc;

eleFile = 'C:\Users\pete\Documents\MATLAB\toolboxes\duf\score\exampleData\APH-L61_03-26_07_Run3_EDU4_20cm_S10_Asphalt_CE_SP_HPF-A.ele';
ascFile = 'C:\Users\pete\Documents\MATLAB\toolboxes\duf\score\exampleData\tlane61.asc';

alarmSet = eleRead(dufAlarmSet,eleFile);
objectSet = ascRead(dufObjectSet,ascFile);

plot(objectSet);
hold on;
plot(alarmSet);
%%
excursion = dufExcursion('alarmSet',alarmSet,'objectSet',objectSet);
excursion = excursion.score;
plot(excursion.scoring.nfa,excursion.scoring.pd)

%%
eleFile2 = 'C:\Users\pete\Documents\MATLAB\toolboxes\duf\score\exampleData\Lane_54_N0247_HDCV_URVIS-3_20120511_radar_1_prescreen.ele';
ascFile2 = 'C:\Users\pete\Documents\MATLAB\toolboxes\duf\score\exampleData\tlane54_may2012.asc';

alarmSet2 = eleRead(dufAlarmSet,eleFile2);
objectSet2 = ascRead(dufObjectSet,ascFile2);

plot(objectSet2);
hold on;
plot(alarmSet2);

excursion2 = dufExcursion('alarmSet',alarmSet2,'objectSet',objectSet2);
%%
eSet = dufExcursionSet('excursionArray',cat(1,excursion,excursion2));