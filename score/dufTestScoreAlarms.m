%%
clear all;
close all;
clc;

eleFile = fullfile(dufRoot,'score','exampleData','APH-L61_03-26_07_Run3_EDU4_20cm_S10_Asphalt_CE_SP_HPF-A.ele');
ascFile = fullfile(dufRoot,'score','exampleData','tlane61.asc');

alarmSet = eleRead(dufAlarmSet,eleFile);
objectSet = ascRead(dufObjectSet,ascFile);

plot(objectSet);
hold on;
plot(alarmSet);
%%
