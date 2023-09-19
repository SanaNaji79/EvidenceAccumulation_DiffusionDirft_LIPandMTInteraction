%% README
%evidence accumulation
%diffusion drift model (DDM) for evidence accumulation
%investigating time interval influences, bias, time limit, and start point
%effect on final decision, reaction time and error rate.
% in the end of part one we implement the Race diffusion model, an extension of drift diffusion model
%part2: Implementing the model proposed by Shadlen and Newsome (2001) for the interaction between area MT and LIP
%%
clc
clear all
%% part1 question 1
B = 1 ;
sigma = 1 ;
T =1 ;
dt = 0.1 ;
start = 0 ;
p = 1 ;
all_choice = [] ;
last_move = [] ;
figure ;
for i = 1:20
    subplot(4 , 5 , i) ;
    [choice , x] = simple_model(B , sigma , dt , T , start , p) ;
    all_choice = [all_choice , choice] ;
    last_move = [last_move , x(end)] ;
    title(['Choice experiment', num2str(i)]) ;
end
%% part 1 question 2
figure ;
histogram(last_move , 5) ;
title('last step histogram') ;
figure ;
histogram(all_choice) ;
title('choice histogram') ;
%% part 1 question 2
B = 1 ;
sigma = 1 ;
T =10 ;
dt = 0.1 ;
start = 0 ;
p = 1 ;
B1 = [1 10 0 0.1 -1 -2] ;
all_choice = [] ;
last_move = [] ;
figure ;
for i = 1:6
    subplot(2 , 3 , i) ;
    [choice , x] = simple_model(B1(i) , sigma , dt , T , start , p) ;
    all_choice = [all_choice , choice] ;
    last_move = [last_move , x(end)] ;
    title(['Choice experiment B = ', num2str(B1(i))]) ;
end
%% part 1 question 3
T1 = (1:1:40) ;
B = 0.1 ;
sigma = 1 ;
dt = 0.1 ;
start = 0 ;
p = 0 ;
error1 = [] ;
trial = 50 ;
for T = T1
    all_choice = [] ;
    last_move = [] ;
    for i = 1:trial
        [choice , x] = simple_model(B , sigma , dt , T , start , p) ;
        all_choice = [all_choice , choice] ;
        last_move = [last_move , x(end)] ;
    end
    fault = length(find(all_choice == -1)) ;
    e1 = fault/trial ;
    error1 = [error1 , e1] ;
end
plot(T1 , error1) ;
ylabel('Error rate') ;
xlabel('Time interval (s)') ;
title('Error rate through time') ;
xlim([min(T1) max(T1)]) ;
%% part 1 question 4
T = 10 ;
B = 0.1 ;
sigma = 1 ;
dt = 0.1 ;
start = 0 ;
p = 1 ;
all_choice = [] ;
last_move = [] ;
path = [] ;
figure ;
for i = 1:10
    [choice , x] = simple_model(B , sigma , dt , T , start , p) ;
    hold on ;
    all_choice = [all_choice , choice] ;
    last_move = [last_move , x(end)] ;
    title(['Choice experiment', num2str(i)]) ;
    path = [path ; x] ;
end
plot((0:dt:T) , mean(path , 1) , 'LineWidth' , 3) ;
title('trajectories') ;
xlabel('Time (s)') ;
%%
plot((0:dt:T) , var(path , 1)) ;
title('trajectories variance') ;
xlabel('Time (s)') ;
%% part 1 question 5
for i = 1:10
    subplot(2 , 5 , i)
    [choice , x] = simple_model(B , sigma , dt , T , 2*i/10 , p) ;
    xlabel('Time (s)') ;
    title(['starting point=' num2str(2*i/10)]) ;
end
%% part 1 question 6
teta1 = 10 ;
teta2 = -10 ;
[choice , x , react] = two_choice_trial(B , sigma , dt , start  , teta1 , teta2 , p) ;
%% part1 question 7
B = 0.2 ;
sigma = 1 ;
dt = 0.1 ;
start = 0 ;
p = 0 ;
teta1 = 10 ;
teta2 = -10 ;
error = [] ;
RT = [] ;
RT1 = [] ;
RT2 = [] ;
for i = 1:1000
    [choice , x , react] = two_choice_trial(B , sigma , dt , start  , teta1 , teta2 , p) ;
    RT = [RT , react] ;
    error = [error , choice] ;
    if choice == 1
        RT1 = [RT1 , react] ;
    else
        RT2 = [RT2 , react] ;
    end
end
%%
scatter(RT*dt , error , 25 , 'filled') ;
xlabel('Reaction time(s)') ;
ylabel('choice') ;
ylim([-1.5 1.5]) ;
%%
scatter(RT*dt , ones(1 , 100) , 25 , error) ;
xlabel('Reaction time(s)') ;
ylabel('choice') ;
%%legend('true' , 'false') ;
%%ylim([-1.5 1.5]) ;
%%
figure 
histogram(RT1 , 40) ;
hold on ;
histogram(RT2 , 40) ;
legend('true' , 'false') ;
title('reaction time histogram') ;
%%
figure 
histogram(RT1 , 40 , 'Normalization' , 'pdf') ;
hold on ;
histogram(RT2 , 40 , 'Normalization' , 'pdf') ;
legend('true' , 'false') ;
title('reaction time histogram') ;
%% part 1 question 8&9
B1 = 0.2 ;
sigma1 = 1 ;
B2 = 0.3 ;
sigma2 = 1.1 ;
dt = 0.1 ;
start = 0 ;
p = 1 ;
T = 30 ;
teta1 = 5 ;
teta2 = 7 ;
for i = 1:6
    subplot(2 , 3 , i) ;
    [choice , x1 , x2 , react] = race_trial(B1 , sigma1 , B2 , sigma2 , dt , T , start  , teta1 , teta2 , p) ;
end
%% part 2 question 1

LIP_weights = [1 ; -1] ;
LIP_threshold = 0.9 ;
M = 10 ;
for i = 1:20
    MT_p_values = [0.5 + 0.025*i ; 0.5 - 0.025*i] ;
    [rt , ex_act , in_act , LIP_act , T] = lip_activity(MT_p_values,LIP_weights,LIP_threshold,M) ;
    subplot(4 , 5 , i) ;
    plotRaster([logical(LIP_act) ; logical(in_act) ; logical(ex_act)]  , T) ;
    xlim([0 rt]) ;
    xlabel('Time(s)') ;
end
%% part 2 question 2
T1 = 1000 ;
percent1 = 0.5 ;
percent2 = 0.5 ;
M = 20 ;
rt_all = [] ;
win_all = [] ;
for i = 1:15
    stimulus = generate(percent1 , percent2 , T1) ;
    [dir1 , dir2] = convert(stimulus) ;
    [rt , ex_act , in_act , LIP_act1 , LIP_act2 , T , win] = lip2_activity(stimulus,LIP_weights,LIP_threshold,M) ;
    subplot(3 , 5 , i) ;
    plotRaster([dir1(1:length(LIP_act1)) ; dir2(1:length(LIP_act1)) ; logical(LIP_act1) ; logical(LIP_act2) ; logical(in_act) ; logical(ex_act)]  , T) ;
    xlabel('Time(s)') ;
    xlim([0 rt]) ;
end
%%
for i = 1:100
    stimulus = generate(percent1 , percent2 , T1) ;
    [dir1 , dir2] = convert(stimulus) ;
    [rt , ex_act , in_act , LIP_act1 , LIP_act2 , T , win] = lip2_activity(stimulus,LIP_weights,LIP_threshold,M) ;
    rt_all = [rt_all , rt] ;
    win_all = [win_all , win] ;
end
figure ; 
histogram(rt_all , 20) ;
xlabel('Reaction Time(s)') ;
correct = length(find(win_all == 1))/length(win_all) ;
%%
correct_all = [] ;
for t = 0:20
    percent1 = t/20 ;
    percent2 = 1 - percent1 ;
    win_all = [] ;
    rt_all = [] ;
    for i = 1:200
        stimulus = generate(percent1 , percent2 , T1) ;
        [dir1 , dir2] = convert(stimulus) ;
        [rt , ex_act , in_act , LIP_act1 , LIP_act2 , T , win] = lip2_activity(stimulus,LIP_weights,LIP_threshold,M) ;
        rt_all = [rt_all , rt] ;
        win_all = [win_all , win] ;
    end
    subplot(3 , 7 , t+1) ;
    h1 = rt_all(find(win_all == 1)) ;
    h2 = rt_all(find(win_all == 2)) ;
    histogram(h1 , 15) ;
    hold on ;
    histogram(h2 , 15) ;
    xlabel('Reaction Time(s)') ;
    title([num2str(t*5) '% right frames']) ;
    correct = length(find(win_all == 1))/length(win_all) ;
    correct_all = [correct_all , correct] ;
end
figure ;
plot(0:5:100 , correct_all) ;
xlabel('percentage of right frames') ;
ylabel('Correct label') ;
%% functions declaration
function [choice , x] = simple_model(B , sigma , dt , T , start , p)
T1 = (0:dt:T) ;
dW = normrnd(0 , sqrt(dt) , [1 , length(T1)-1]) ;
dx = B*dt + sigma*dW ;
x = cumsum([start dx]) ;
if x(end)>0
    choice = 1 ;
else
    choice = -1 ;
end
if p
    plot(T1 , x) ;
    xlabel('Time (s)') ;
end
end
function [choice , x , react] = two_choice_trial(B , sigma , dt , start  , teta1 , teta2 , p)
react = 0 ;
x = [start] ;
t = 0 ;
while react == 0
    dW = normrnd(0 , sqrt(dt)) ;
    dx = B*dt + sigma*dW ;
    x = [x , x(end)+dx] ;
    t = t+ dt ;
    if x(end)>teta1
        react = t ;
        choice = 1 ;
    elseif x(end)<teta2
        react = t ;
        choice = -1 ;
    end
end
if p 
    plot((0:dt:t) , x) ;
    hold on ;
    yline(teta1 , 'LineWidth' , 2) ;
    yline(teta2 , 'LineWidth' , 2) ;
    xlabel('Time (s)') ;
    ylim([teta2-1 teta1+1]) ;
    hold off ;
end
end
function [choice , x1 , x2 , react] = race_trial(B1 , sigma1 , B2 , sigma2 , dt , T , start  , teta1 , teta2 , p)
react = 0 ;
x1 = [start] ;
x2 = [start] ;
t = 0 ;
while (react == 0)&&(t<T)
    dW1 = normrnd(0 , sqrt(dt)) ;
    dx1 = B1*dt + sigma1*dW1 ;
    x1 = [x1 , x1(end)+dx1] ;
    dW2 = normrnd(0 , sqrt(dt)) ;
    dx2 = B2*dt + sigma2*dW2 ;
    x2 = [x2 , x2(end)+dx2] ;
    t = t+ dt ;
    if x1(end)>teta1
        react = t ;
        choice = 1 ;
    elseif x2(end)>teta2
        react = t ;
        choice = -1 ;
    end
end
if react == 0
    if x1(end)>x2(end)
        choice = 1 ;
        react = T ;
        t = T ;
    else
        choice = -1 ;
        react = T ;
        t = T ;
    end
end
if p 
    plot(linspace(0 , t , length(x1)) , x1) ;
    hold on ;
    plot(linspace(0 , t , length(x1)) , x2) ;
    yline(teta1 , 'LineWidth' , 2) ;
    yline(teta2 , 'LineWidth' , 2) ;
    xlabel('Time (s)') ;
    hold off ;
    xlim([0 t]) ;
end
end

function [rt , ex_act , in_act , LIP_act , T] = lip_activity(MT_p_values,LIP_weights,LIP_threshold,M)
% Parameters:
% MT_p_values - a vector with 2 elements, firing probabilities for the
% excitatory and inhibitory neurons, resp.
% LIP_weights - a length 2 vector of weighting factors for the evidence
% from the excitatory (positive) and
% inhibitory (negative) neurons
% LIP_threshold - the LIP firing rate that represents the choice threshold criterion
% use fixed time scale of 1 ms
dt = 0.001;
N = [0 0]'; % plus is first, minus is second
rate = 0.0;
LIP_event_times = [];
ex_act = [] ;
in_act = [] ;
LIP_act = [] ;
N_LIP = 0 ; 
T = [] ;
t = 0 ;
while rate < LIP_threshold
    t = t + dt;
    T = [T , t] ;
    dN = rand(2 , 1) < MT_p_values;
    ex_act = [ex_act , dN(1)] ;
    in_act = [in_act , dN(2)] ;
    N = N + dN;
    p_LIP = sum(N.*LIP_weights);
    LIP_event = rand < p_LIP;
    LIP_act = [LIP_act , LIP_event] ;
    if LIP_event
        LIP_event_times = [LIP_event_times t];
        N_LIP = N_LIP + 1 ;
    end
    % check LIP mean rate for last M spikes
    if N_LIP>M
        rate = M/(t - LIP_event_times(N_LIP - M));
        rate = rate/1000 ;
    end
end
rt = t;
end
function [] = plotRaster(spikeMat, tVec)
hold all;
for trialCount = 1:size(spikeMat,1)
    spikePos = tVec(spikeMat(trialCount, :));
    for spikeCount = 1:length(spikePos)
        plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-0.2 trialCount+0.2], 'Color',[(trialCount == 6) , trialCount == 5 ,(trialCount == 3) + (trialCount == 4)] );
%         plot([spikePos(spikeCount) spikePos(spikeCount)], ...
%             [trialCount-0.2 trialCount+0.2], 'Color',[(trialCount == 2) , trialCount == 1 ,(trialCount == 4) + (trialCount == 3)] );
    end
end
ylim([0 size(spikeMat, 1)+1]);
end
function [stimulus] = generate(percent1 , percent2 , T)
stimulus = zeros(T , 1) ;
for i = 1:T
    a = rand ;
    if a<percent1
        stimulus(i) = 1 ;
    elseif a<percent1+percent2
        stimulus(i) = -1 ;
    end
end
end
function [rt , ex_act , in_act , LIP_act1 , LIP_act2 , T , win] = lip2_activity(stimulus,LIP_weights,LIP_threshold,M)
dt = 0.001;
N = [0 0]'; 
rate1 = 0.0;
rate2 = 0.0;
LIP1_event_times = [];
LIP2_event_times = [];
ex_act = [] ;
in_act = [] ;
LIP_act1 = [] ;
LIP_act2 = [] ;
N_LIP1 = 0 ; 
N_LIP2 = 0 ;
T = [] ;
t = 0 ;
while (rate1 < LIP_threshold)&&(rate2 < LIP_threshold)
    t = t + dt;
    T = [T , t] ;
    if stimulus(round(t/dt))==1
        dN = rand(2 , 1)<[0.9 ; 0.1] ;
    elseif stimulus(round(t/dt))==-1
        dN = rand(2 , 1)<[0.1 ; 0.9] ;
    else
        a = length(find(stimulus(1:t/dt)/(t/dt)==1)) ;
        b = length(find(stimulus(1:t*dt)/(t/dt)==-1)) ;
        dN = rand(2 , 1)<[a ; b] ;
    end
    ex_act = [ex_act , dN(1)] ;
    in_act = [in_act , dN(2)] ;
    N = N + dN;
    p_LIP1 = sum(N.*LIP_weights);
    p_LIP2 = -sum(N.*LIP_weights);
    LIP1_event = rand < p_LIP1;
    LIP2_event = rand < p_LIP2;
    LIP_act1 = [LIP_act1 , LIP1_event] ;
    LIP_act2 = [LIP_act2 , LIP2_event] ;
    if LIP1_event
        LIP1_event_times = [LIP1_event_times t];
        N_LIP1 = N_LIP1 + 1 ;
    end
    if LIP2_event
        LIP2_event_times = [LIP2_event_times t];
        N_LIP2 = N_LIP2 + 1 ;
    end
    % check LIP mean rate for last M spikes
    if N_LIP1>M
        rate1 = M/(t - LIP1_event_times(N_LIP1 - M));
        rate1 = rate1/1000 ;
    end
    if N_LIP2>M
        rate2 = M/(t - LIP2_event_times(N_LIP2 - M));
        rate2 = rate2/1000 ;
    end
end
rt = t;
if (rate1 >= LIP_threshold)
    win = 1 ;
else 
    win = 2 ;
end
end
function [output1 , output2] = convert(stimulus)
T = length(stimulus) ;
output1 = zeros(T) ;
output2 = zeros(T) ;
output1(find(stimulus==1)) = 1 ;
output2(find(stimulus==-1)) = 1 ;
output1 = logical(output1) ;
output2 = logical(output2) ;
end