%% Resample acceleration time history with variable time step size

%% Statement of the problem
% * The linear dynamic response history analysis (DRHA) algorithms that are
% included in this package work only for acceleration time histories
% defined in terms of time step with constant size. This may not be the
% case in acceleration time histories with nonuniform time step (i.e.
% variable time step size). In this case, the time history needs to be
% resampled, so that an equivalent time history is defined with constant
% time step size, suitable for use in the various functions of this
% package. Here an example is provided for converting an acceleration time
% history with nonuniform time step into an equivalent acceleration time
% history with uniform (constant size) time step.
%
%% Load earthquake data
% Initial acceleration time history
D=load('elcentro_truncated.dat');
dt=D(2,1)-D(1,1);
t=D(:,1);
xgtt=D(:,2);
%% Define resampling parameters
% Set the desired time step of the new acceleration time history that is
% produced after resampling.
dt_new=0.02;
%%
% Upsampling factor
p=1;
%%
% Downsampling factor
q=1;
%% Resampling procedure
% Sample rate
fs = 1/dt_new;
%%
% Resample the original acceleration time history into a new one with
% constant time step size equal to dt_new
[xgtt_new,t_new] = resample(xgtt,t,fs,p,q);
%% Plot original and new acceleration time history
%
FigHandle=figure('Name','Acceleration time histories','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 500, 300]);
plot(t,xgtt,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[0 0 0],'markeredgecolor','k')
hold on
plot(t_new,xgtt_new,'LineWidth',1.,'Marker','.',...
    'MarkerSize',1,'Color',[1 0 0],'markeredgecolor','r')
grid on
xlim([0,30])
ylim([-3.5,3.5])
xlabel('Time (sec)','FontSize',10);
ylabel('A (m/s^2)','FontSize',10);
title(['Acceleration time histories'],...
    'FontSize',10)
legend('Original time history','Resampled time history')
%%
% Check uniformity of the time step of the original acceleration time
% history
if any(diff(diff(t))>1e-14)
    disp('The original acceleration time history is nonuniform')
else
    disp('The original acceleration time history is uniform')
end
%%
% Check uniformity of the time step of the resampled acceleration time
% history
if any(diff(diff(t_new))>1e-14)
    disp('The resampled acceleration time history is nonuniform')
else
    disp('The resampled acceleration time history is uniform')
end

%% Copyright
%
% Copyright (c) 2015-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%
