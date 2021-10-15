clear all
close all

%% define the data paths
cwd = pwd;
addpath(genpath('~/Dropbox/Utils/fieldtrip/external/dmlt/external/gpstuff/misc/'));
pathToData = '~/Documents/Data/Max EEG/Github/data/';
data_path_pre=[pathToData 'Classifier_results/Pre_training/'];
data_path_post=[pathToData 'Classifier_results/Post_training/'];

%% define all subjects; 1-20= group 2; 21-40=group 1 
Subjects_all={'101','103', '106', '108','113','114','115','118','119','120','123','125','127','131','134','136','137','139','141','145','102','104','105', '107', '109','110','111','116','117', '121','122','124','130','135','138','140','142', '143','144','146'};

%% define settings
burn_in=30; % define minimum of time points that have to be used for calculating evidence accumulation (here 150ms)
start_index=1; 

smoothed=1; 
smooth_average=0; 
smoothin_kernel=10;


Matrix_slope=[];
Matrix_intercept=[];
Matrix_unsigned_slope=[];
Matrix_unsigned_intercept=[];
Matrix_stimulus_representation=[]; 
Matrix_choice=[];
Matrix_choice2=[];
Matrix_change=[];
Matrix_conf=[];
Matrix_conf2=[];
Matrix_acc=[];
Matrix_acc2=[];
Matrix_direction=[];
Matrix_post=[];
Matrix_RT=[];
Matrix_RT2=[];
Matrix_subject=[];
Matrix_session=[];
Matrix_group=[]; 
Matrix_decoding=[];
Matrix_difficulty=[];
Matrix_frequency_tagging=[];

cd(data_path_pre)

for sj_all=1:length(Subjects_all)
    
    load([data_path_pre 'classifier_results' Subjects_all{sj_all} '_pre3.mat'])
    
    %% smooth decoding accuracy to find a reliable timepoint of the highest decoding accuracy
    [smoothed_accuracy,smoothin_kernel]  = smoothdata(AUC_gen);
    [end_point end_index]=max(smoothed_accuracy);
    mean_dec=mean(smoothed_accuracy(1:end_index)); % define overall decodability for this session 
  
    if end_index<burn_in %if heighest decodability is too early (too little data for calculating evidence processing) find a later time-point of optimal decodability
        [end_point end_index]=max(smoothed_accuracy(burn_in:end));
        end_index=end_index+burn_in
    end
    if end_index==171
    end_index=170
    end
    
    %% Calcualte slope and intercept for each trial
    intercept=[];
    slope=[];
    outcome_glm=prob_outcome_glm;

    for trial=1:length(outcome_glm(1,:))
        fit_trial=fitlm([start_index:length(outcome_glm(1: end_index,trial))],outcome_glm(start_index: end_index,trial));
        intercept(trial)=fit_trial.Coefficients.Estimate(1);
        slope(trial)=fit_trial.Coefficients.Estimate(2);
        
    end
    
    %% zscore slope and intercept and flip them for leftwards decisiosn to get a stimulus indepent measure of evidence processing
    zscored_slope=zscore(slope);
    zscored_intercept=zscore(intercept);
    unsigned_slope=zscored_slope;
    unsigned_slope(behT.Motion_direction(good_trial)==1)=-unsigned_slope(behT.Motion_direction(good_trial)==1);
    unsigned_intercept=zscored_intercept;
    unsigned_intercept(behT.Motion_direction(good_trial)==1)=-unsigned_intercept(behT.Motion_direction(good_trial)==1);

    %% save neural variables for the hierarchical model
    Matrix_slope=[Matrix_slope; zscored_slope'];
    Matrix_unsigned_slope=[Matrix_unsigned_slope; unsigned_slope'];
    Matrix_unsigned_intercept=[Matrix_unsigned_intercept; unsigned_intercept'];
    Matrix_intercept=[Matrix_intercept;  zscored_intercept'];
    Matrix_decoding=[Matrix_decoding;repmat(mean_dec, length(zscored_intercept), 1)];
    
    %% change coding of behavioral varibales for usage in the hierarchical regression
    behT.Accuracy_initial( behT.Accuracy_initial==0)=-1;
    behT.Confidence_initial(behT.Confidence_initial==1)=-1;  
    behT.Confidence_initial(behT.Confidence_initial==2)=1;  
    behT.Motion_direction(behT.Motion_direction==1)=-1;
    behT.Motion_direction(behT.Motion_direction==2)=1;
    behT.Confidence_final(behT.Confidence_final==1)=0;
    behT.Confidence_final(behT.Confidence_final==2)=1;
    
    %% save the behavioral data as predictors for the hierarchical regression
    if sj_all<21 
        Matrix_group=[Matrix_group; repmat(1, length(zscored_intercept), 1)];
    else
              Matrix_group=[Matrix_group; repmat(-1, length(zscored_intercept), 1)];  
    end
    Matrix_subject=[Matrix_subject; repmat(sj_all, length(zscored_intercept), 1)];
    Matrix_session=[Matrix_session; repmat(-1, length(zscored_intercept), 1)];
    Matrix_choice=[Matrix_choice; ClassifierLabels];
    Matrix_acc=[Matrix_acc; behT.Accuracy_initial(good_trial)];
    Matrix_conf=[Matrix_conf; behT.Confidence_initial(good_trial)];
    Matrix_direction=[Matrix_direction; behT.Motion_direction(good_trial)];
    Matrix_post=[Matrix_post; behT.Coherence_Post(good_trial)];
    Matrix_RT2=[Matrix_RT2;  behT.RT_final_Type1_decision(good_trial)];
    Matrix_RT=[Matrix_RT;  behT.RT_Initial_Type1_decision(good_trial)];
    Matrix_difficulty=[Matrix_difficulty;  behT.Coherence_pre_strength(good_trial)];
    Matrix_frequency_tagging=[Matrix_frequency_tagging;behT.Frequency_tagging(good_trial)]
    
    Matrix_conf2=[Matrix_conf2; behT.Confidence_final(good_trial)];
    Matrix_acc2=[Matrix_acc2; behT.Accuracy_final(good_trial)];
    Matrix_choice2=[Matrix_choice2;  behT.Final_Type1_decision(good_trial)];
    
    %% assign overall metacognitive ability, confirmation bias and neural evidence integration
    config.beh_path = [pathToData 'Behavior_EEG_Session/InitialAssessment/'];
    behT2=load([config.beh_path 'Behavior_metaD_pre_subject' Subjects_all{sj_all} '.mat'])
    metacognitive_ability_pre(sj_all)=behT2.behT.meta_d./behT2.behT.d_prime;
    beh_post_integration_pre(sj_all)=behT2.behT.d_prime_2-behT2.behT.d_prime;
    staircase_variability_pre(sj_all)=std(behT.Coherence_pre_strength);
    difficulty_pre(sj_all)=mean(behT.Coherence_pre_strength);

    
    
    %% save summary statistics for each condition (this was already calculate for each participant in the pre-processing
    left_sj_pre(sj_all, :)=left;
    right_sj_pre(sj_all, :)=right; 
    
    out_confirm_sj_pre(sj_all, :)=out_confirm;
    out_disconfirm_sj_pre(sj_all, :)=out_dicconfirm;

    
    highest_decodability_time_pre(sj_all)=end_index;
    highest_decodability_pre(sj_all)=end_point;
    
    %% number of trials per participant per condition (will be used for calculating weighted means)

    amount_confirm_pre(sj_all)=length(find(behT.Accuracy_initial(good_trial)==1));
    amount_disconfirm_pre(sj_all)=length(find(behT.Accuracy_initial(good_trial)==-1));

    
    clear prob_outcome_glm slope intercept stim_rep alternative_intercept alternative_slope end_points_window
end


for sj_all=1:length(Subjects_all)
    
    load([data_path_post 'classifier_results' Subjects_all{sj_all} '_post3.mat'])
    
    %% smooth decoding accuracy to find a reliable timepoint of the highest decoding accuracy
    smoothed_accuracy = smoothdata(AUC_gen);
    [end_point end_index]=max(smoothed_accuracy);
    mean_dec=mean(smoothed_accuracy(1:end_index));% define overall decodability for this session       
    if end_index<burn_in
        [end_point end_index]=max(smoothed_accuracy(burn_in:end));  %if heighest decodability is too early (too little data for calculating evidence processing) find a later time-point of optimal decodability
        end_index=end_index+burn_in
    end
    if end_index==171
    end_index=170
    end
    
    
    %% Calcualte slope and intercept for each trial   
    intercept=[];
    slope=[];
    outcome_glm=prob_outcome_glm;
    for trial=1:length(outcome_glm(1,:))
        fit_trial=fitlm([start_index:length(outcome_glm(1: end_index,trial))],outcome_glm(start_index: end_index,trial));
        intercept(trial)=fit_trial.Coefficients.Estimate(1);
        slope(trial)=fit_trial.Coefficients.Estimate(2);        
    end
    
    %% zscore slope and intercept and flip them for leftwards decisiosn to get a stimulus indepent measure of evidence processing
    zscored_slope=zscore(slope);
    zscored_intercept=zscore(intercept);
    unsigned_slope=zscored_slope;
    unsigned_slope(behT.Motion_direction(good_trial)==1)=-unsigned_slope(behT.Motion_direction(good_trial)==1);
    unsigned_intercept=zscored_intercept;
    unsigned_intercept(behT.Motion_direction(good_trial)==1)=-unsigned_intercept(behT.Motion_direction(good_trial)==1);
   
    %% save neural variables for the hierarchical model
    Matrix_slope=[Matrix_slope; zscored_slope'];
    Matrix_unsigned_slope=[Matrix_unsigned_slope; unsigned_slope'];
    Matrix_unsigned_intercept=[Matrix_unsigned_intercept; unsigned_intercept'];
    Matrix_intercept=[Matrix_intercept;  zscored_intercept'];
    Matrix_decoding=[Matrix_decoding;repmat(mean_dec, length(zscored_intercept), 1)];
       
    %% change coding of behavioral varibales for usage in the hierarchical regression
    behT.Accuracy_initial( behT.Accuracy_initial==0)=-1;
    behT.Confidence_initial(behT.Confidence_initial==1)=-1;  
    behT.Confidence_initial(behT.Confidence_initial==2)=1;  
    behT.Motion_direction(behT.Motion_direction==1)=-1;
    behT.Motion_direction(behT.Motion_direction==2)=1;
    behT.Confidence_final(behT.Confidence_final==1)=0;
    behT.Confidence_final(behT.Confidence_final==2)=1;
    
    %% save the behavioral data as predictors for the hierarchical regression
    if sj_all<21
        Matrix_group=[Matrix_group; repmat(1, length(zscored_intercept), 1)];
    else
        Matrix_group=[Matrix_group; repmat(-1, length(zscored_intercept), 1)];
    end
    
    Matrix_subject=[Matrix_subject; repmat(sj_all, length(zscored_intercept), 1)];
    Matrix_session=[Matrix_session; repmat(1, length(zscored_intercept), 1)];

    Matrix_choice=[Matrix_choice; ClassifierLabels];
    Matrix_acc=[Matrix_acc; behT.Accuracy_initial(good_trial)];
    Matrix_conf=[Matrix_conf; behT.Confidence_initial(good_trial)];
    Matrix_direction=[Matrix_direction; behT.Motion_direction(good_trial)];
    Matrix_post=[Matrix_post; behT.Coherence_Post(good_trial)];
    Matrix_RT2=[Matrix_RT2;  behT.RT_final_Type1_decision(good_trial)];
    Matrix_RT=[Matrix_RT;  behT.RT_Initial_Type1_decision(good_trial)];
    Matrix_conf2=[Matrix_conf2; behT.Confidence_final(good_trial)];
    Matrix_acc2=[Matrix_acc2; behT.Accuracy_final(good_trial)];
    Matrix_choice2=[Matrix_choice2;  behT.Final_Type1_decision(good_trial)];
    Matrix_difficulty=[Matrix_difficulty;  behT.Coherence_pre_strength(good_trial)];
    Matrix_frequency_tagging=[Matrix_frequency_tagging;behT.Frequency_tagging(good_trial)];


    config.beh_path = [pathToData 'Behavior_EEG_Session/FinalAssessment/'];
    behT2=load([config.beh_path 'Behavior_metaD_post_subject' Subjects_all{sj_all} '.mat'])
    metacognitive_ability_post(sj_all)=behT2.behT.meta_d./behT2.behT.d_prime;
    beh_post_integration_post(sj_all)=behT2.behT.d_prime_2-behT2.behT.d_prime;

    staircase_variability_post(sj_all)=std(behT.Coherence_pre_strength);
    difficulty_post(sj_all)=mean(behT.Coherence_pre_strength);

    
    
    %% save summary statistics for each condition (this was already calculate for each participant in the pre-processing
    left_sj_post(sj_all, :)=left;
    right_sj_post(sj_all, :)=right; 
    
    out_confirm_sj_post(sj_all, :)=out_confirm;
    out_disconfirm_sj_post(sj_all, :)=out_dicconfirm;
        
    highest_decodability_time_post(sj_all)=end_index;
    highest_decodability_post(sj_all)=end_point;
    
    
    %% number of trials per participant per condition (will be used for calculating weighted means)
    amount_confirm_post(sj_all)=length(find(behT.Accuracy_initial(good_trial)==1));
    amount_disconfirm_post(sj_all)=length(find(behT.Accuracy_initial(good_trial)==-1));

    clear prob_outcome_glm slope intercept stim_rep alternative_intercept alternative_slope end_points_window
end

%% conduct hierarchical regression to test the change in neural processing between sessions
input_table_Confidence = table(Matrix_unsigned_slope,Matrix_unsigned_intercept, Matrix_acc, Matrix_conf,Matrix_post,Matrix_session,  Matrix_subject, Matrix_group,Matrix_decoding,Matrix_difficulty, Matrix_frequency_tagging);
input_table_Confidence.Properties.VariableNames = {'unsigned_slope','unsigned_intercept','acc','confidence','post','session','Subject', 'group','decodability','difficulty','frequency_tag'};
fit_slope=fitglme(input_table_Confidence,'unsigned_slope ~ 1 +confidence*acc*session*group+difficulty*post+decodability+unsigned_intercept+(1 | Subject)','Distribution','Normal')


%% calculate the weighted group averages of the decoding timeline for figures 3B-D
time_highest_decodability=median([highest_decodability_time_pre highest_decodability_time_post]);

% calculate for each participant the amount of trials in each condition
amount_confirm_pre_use=amount_confirm_pre./sum(amount_confirm_pre);
amount_disconfirm_pre_use=amount_disconfirm_pre./sum(amount_disconfirm_pre);
amount_confirm_post_use=amount_confirm_post./sum(amount_confirm_post);
amount_disconfirm_post_use=amount_disconfirm_post./sum(amount_disconfirm_post);
amount_confirm_all=amount_confirm_post+amount_confirm_pre./sum(amount_confirm_post+amount_confirm_pre);
amount_disconfirm_all=amount_disconfirm_post+amount_disconfirm_pre./sum(amount_disconfirm_post+amount_disconfirm_pre);

for k=1:size(left_sj_pre,2)

% generate a group average value for each time point
out_confirm_sj_pre_use=out_confirm_sj_pre;
out_disconfirm_sj_pre_use=out_disconfirm_sj_pre;

out_confirm_sj_post_use=out_confirm_sj_post;
out_disconfirm_sj_post_use=out_disconfirm_sj_post;    

average_confirm_pre(k)=wmean(out_confirm_sj_pre_use(:, k), amount_confirm_pre_use');
average_disconfirm_pre(k)=wmean(out_disconfirm_sj_pre_use(:, k), amount_disconfirm_pre_use');

average_confirm_post(k)=wmean(out_confirm_sj_post_use(:, k), amount_confirm_post_use');
average_disconfirm_post(k)=wmean(out_disconfirm_sj_post_use(:, k), amount_disconfirm_post_use');

average_confirm_diff(k)=average_confirm_post(k)-average_confirm_pre(k)
average_disconfirm_diff(k)=average_disconfirm_post(k)-average_disconfirm_pre(k);
end

cd(cwd);

%% Figure 3B
figure(1)
hold on
fit5=fitlm([1:time_highest_decodability], average_disconfirm_pre(1:time_highest_decodability))
ypred = predict(fit5,[1:time_highest_decodability]')
change=plot([1:time_highest_decodability], average_disconfirm_pre(1:time_highest_decodability),'Color', [51/255 153/255 255/255],'LineWidth',2)
plot([1:time_highest_decodability], ypred, '-k','LineWidth',1.3)
fit6=fitlm([1:time_highest_decodability], average_confirm_pre(1:time_highest_decodability))
ypred = predict(fit6,[1:time_highest_decodability]')
no=plot([1:time_highest_decodability], average_confirm_pre(1:time_highest_decodability),'Color', [255/255 153/255 51/255],'LineWidth',2)
plot([1:time_highest_decodability], ypred, '-k','LineWidth',1.3)
ylim([-.04 .07])
xlim([1 100])
pbaspect([.75 1 1])
legend([no change], {'confirming', 'disconfirming'}, 'Location', 'NorthWest','box','off')
title('Before training')
ylabel('representation motion direction')
xlabel('post-decision period (ms)')
set(gca, 'FontSize', 16,'FontName','Arial','FontWeight','bold','box','off', 'YTick',[-.15 -.1 -.05 0 .05 .1 .15], 'XTick',[20 40 60 80 100],'XTickLabel',{'100','200', '300','400','500'})

%% Figure 3C
figure(2)
hold on
fit5=fitlm([1:time_highest_decodability], average_disconfirm_post(1:time_highest_decodability))
ypred = predict(fit5,[1:time_highest_decodability]')
change=plot([1:time_highest_decodability], average_disconfirm_post(1:time_highest_decodability),'Color', [51/255 153/255 255/255],'LineWidth',2)
plot([1:time_highest_decodability], ypred, '-k','LineWidth',1.3)
fit6=fitlm([1:time_highest_decodability], average_confirm_post(1:time_highest_decodability))
ypred = predict(fit6,[1:time_highest_decodability]')
no=plot([1:time_highest_decodability], average_confirm_post(1:time_highest_decodability),'Color', [255/255 153/255 51/255],'LineWidth',2)
plot([1:time_highest_decodability], ypred, '-k','LineWidth',1.3)
ylim([-.04 .07])
xlim([1 100])
pbaspect([.75 1 1])
legend([no change], {'confirming', 'disconfirming'}, 'Location', 'NorthWest','box','off')
title('After training')
ylabel('representation motion direction')
xlabel('post-decision period (ms)')
set(gca, 'FontSize', 16,'FontName','Arial','FontWeight','bold','box','off', 'YTick',[-.15 -.1 -.05 0 .05 .1 .15], 'XTick',[20 40 60 80 100],'XTickLabel',{'100','200', '300','400','500'})

%% Figure 3D
figure(3)
hold on
fit5=fitlm([1:time_highest_decodability], average_disconfirm_diff(1:time_highest_decodability))
ypred = predict(fit5,[1:time_highest_decodability]')
change=plot([1:time_highest_decodability], average_disconfirm_diff(1:time_highest_decodability),'Color', [51/255 153/255 255/255],'LineWidth',2)
plot([1:time_highest_decodability], ypred, '-k','LineWidth',1.3)
fit6=fitlm([1:time_highest_decodability], average_confirm_diff(1:time_highest_decodability))
ypred = predict(fit6,[1:time_highest_decodability]')
no=plot([1:time_highest_decodability], average_confirm_diff(1:time_highest_decodability),'Color', [255/255 153/255 51/255],'LineWidth',2)
plot([1:time_highest_decodability], ypred, '-k','LineWidth',1.3)
ylim([-.04 .07])
xlim([1 100])
pbaspect([.75 1 1])
legend([no change], {'confirming', 'disconfirming'}, 'Location', 'NorthWest','box','off')
title('Training effect')
ylabel('representation motion direction')
xlabel('post-decision period (ms)')
set(gca, 'FontSize', 16,'FontName','Arial','FontWeight','bold','box','off', 'YTick',[-.15 -.1 -.05 0 .05 .1 .15], 'XTick',[20 40 60 80 100],'XTickLabel',{'100','200', '300','400','500'})





%% Figure 3E
for sj_all=1:length(Subjects_all)
    Mbefore(sj_all)=mean(Matrix_unsigned_slope(Matrix_subject==sj_all& Matrix_session==-1));
    Mafter(sj_all)=mean(Matrix_unsigned_slope(Matrix_subject==sj_all& Matrix_session==1));
end

processing_before=mean(Mbefore);
processing_after=mean(Mafter);

err_before=std(Mbefore)/sqrt(length(Mbefore));
err_after=std(Mafter)/sqrt(length(Mafter));

figure(4)
hold on
bar([1], [processing_before  ], 'BarWidth', .3)
bar([2], [ processing_after ], 'BarWidth', .3)
plot(repmat(1, 1, length(Mbefore)), Mbefore,'o','MarkerSize',3, 'MarkerEdgeColor',[.6 .6 .6],'MarkerFaceColor',[.6 .6 .6])
plot(repmat(2, 1, length(Mafter)), Mafter,'o','MarkerSize',3, 'MarkerEdgeColor',[.6 .6 .6],'MarkerFaceColor',[.6 .6 .6])
errorbar([1,2], [processing_before processing_after], [err_before err_after], '.k' ,'LineWidth', 1.7)
pbaspect([.7 1 1])
xlim([.5 2.5])
ylim([-.132 .19])
ylabel('Neural evidence integration')
set(gca, 'FontSize', 16,'FontName','Arial','FontWeight','bold','box','off', 'XTick',[1:2], 'XTickLabel',{'Before training','After training'})
fix_xticklabels(gca,2,{'FontSize',16,'FontName','Arial','FontWeight','bold'});




%% relate individual differences in metacognitive training effects to individual differences in neural training effects
% include random effects per participant for sesssion to get a neural
% measure of the change in post-decision evidence processing
fit_slope_ind=fitglme(input_table_Confidence,'unsigned_slope ~ 1 +confidence*acc*session*group+difficulty*post+decodability+unsigned_intercept+(1 | Subject)+(session | Subject)','Distribution','Normal')
[b,Bnames,stats]=randomEffects(fit_slope_ind)
change_neuro=b((length(b)-(length(Subjects_all)-1)*2)-1:2:(length(b))-1)

% define metacognitive training effects
change_metacognition=metacognitive_ability_post-metacognitive_ability_pre; 


% predict neural changes with training effects in metacognition
fitlm([zscore(change_metacognition)'], zscore(change_neuro)', 'RobustOpts','on')


%% Figure 3F
z_scored_neural_change = zscore(change_neuro);
figure(5)
hold on
scatter(change_metacognition, zscore(change_neuro), 'o', 'filled')
lsline
 xlabel('Training effect metacognition')
ylabel('Training effect neural integration')
set(findall(gca, 'Type', 'Line'),'LineWidth',3)
set(gca, 'FontSize', 16,'FontName','Arial','FontWeight','bold','box','off')
xlim([-.6 1.5])
plot(change_metacognition(1:20), z_scored_neural_change(1:20),'o','MarkerSize',5, 'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[.4 .4 .4])
plot(change_metacognition(21:end), z_scored_neural_change(21:end),'o','MarkerSize',5, 'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[.4 .4 .4]) % different color scheme [1 .5 .3]

  