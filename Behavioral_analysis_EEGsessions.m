
close all
clear all
%% add paths
pathToData = '~/Documents/Data/Max EEG/Github/data/';

% add path for the mediation toolbox
addpath(genpath('~/Dropbox/Utils/fmri/CanlabCore'))
addpath(genpath('~/Dropbox/Utils/fmri/MediationToolbox'))

%% setup config
config = [];

%% load the subjects for each group
Subjects_control_group={'101','103', '106', '108','113','114','115','118','119','120','123','125', '127', '131','134', '136', '137','139','141','145'}; % group 1
Subjects_experimental_group={'102','104','105', '107', '109','110','111','116','117', '121','122','124','130','135','138','140','142', '143','144','146'}; % group 2

%% load data from initial assessment 


%% group 1
for subj=1:length(Subjects_control_group)
    
config.ID =Subjects_control_group{subj}
config.beh_path = [pathToData 'Behavior_EEG_Session/InitialAssessment/'];
load([config.beh_path 'Behavior_metaD_pre_subject' config.ID '.mat'])


%% define variables
meta_d_control_group_pre(subj)=behT.meta_d;
d_prime_control_group_pre(subj)=behT.d_prime;
M_ratio_control_group_pre(subj)=meta_d_control_group_pre(subj)/d_prime_control_group_pre(subj);
profit_new_evidence_control_group_pre(subj)=behT.d_prime_2-behT.d_prime;
change_error_control_group_pre(subj)=behT.change_error;
change_correct_control_group_pre(subj)=behT.change_correct;

staircase_variability_control_pre(subj)=std(behT.Coherence_pre_strength);
difficulty_control_pre(subj)=mean(behT.Coherence_pre_strength);


end

%% group 2
for subj=1:length(Subjects_experimental_group)
    
config.ID =Subjects_experimental_group{subj}
config.beh_path = [pathToData 'Behavior_EEG_Session/InitialAssessment/'];
load([config.beh_path 'Behavior_metaD_pre_subject' config.ID '.mat']) % load behavioral data 


%% define variables

meta_d_exp_group_pre(subj)=behT.meta_d;
d_prime_exp_group_pre(subj)=behT.d_prime;
M_ratio_exp_group_pre(subj)=meta_d_exp_group_pre(subj)/d_prime_exp_group_pre(subj);
profit_new_evidence_exp_group_pre(subj)=behT.d_prime_2-behT.d_prime;
change_error_exp_group_pre(subj)=behT.change_error;
change_correct_exp_group_pre(subj)=behT.change_correct;

staircase_variability_exp_pre(subj)=std(behT.Coherence_pre_strength);
difficulty_exp_pre(subj)=mean(behT.Coherence_pre_strength);

end

%% Load data from the second assessment
%% group 1
for subj=1:length(Subjects_control_group)
    
config.ID =Subjects_control_group{subj}
config.beh_path = [pathToData 'Behavior_EEG_Session/FinalAssessment/'];
load([config.beh_path 'Behavior_metaD_post_subject' config.ID '.mat'])

%% define variables 
meta_d_control_group_post(subj)=behT.meta_d;
d_prime_control_group_post(subj)=behT.d_prime;
M_ratio_control_group_post(subj)=meta_d_control_group_post(subj)/d_prime_control_group_post(subj);
profit_new_evidence_control_group_post(subj)=behT.d_prime_2-behT.d_prime;
change_error_control_group_post(subj)=behT.change_error;
change_correct_control_group_post(subj)=behT.change_correct;

staircase_variability_control_post(subj)=std(behT.Coherence_pre_strength);
difficulty_control_post(subj)=mean(behT.Coherence_pre_strength);

end

%% group 2
for subj=1:length(Subjects_experimental_group)
    
config.ID =Subjects_experimental_group{subj}
config.beh_path = [pathToData 'Behavior_EEG_Session/FinalAssessment/'];
load([config.beh_path 'Behavior_metaD_post_subject' config.ID '.mat'])

%% define variables 
meta_d_exp_group_post(subj)=behT.meta_d;
d_prime_exp_group_post(subj)=behT.d_prime;
M_ratio_exp_group_post(subj)=meta_d_exp_group_post(subj)/d_prime_exp_group_post(subj);

profit_new_evidence_exp_group_post(subj)=behT.d_prime_2-behT.d_prime;
change_error_exp_group_post(subj)=behT.change_error;
change_correct_exp_group_post(subj)=behT.change_correct;

staircase_variability_exp_post(subj)=std(behT.Coherence_pre_strength);
difficulty_exp_post(subj)=mean(behT.Coherence_pre_strength);
end

%% combine data from both groups for ANOVA's
meta_d_pre=[M_ratio_exp_group_pre M_ratio_control_group_pre];
meta_d_post=[M_ratio_exp_group_post M_ratio_control_group_post];
confirmation_pre=[profit_new_evidence_exp_group_pre profit_new_evidence_control_group_pre];
confirmation_post=[profit_new_evidence_exp_group_post profit_new_evidence_control_group_post];
change_error_pre=[change_error_exp_group_pre change_error_control_group_pre];
change_error_post=[change_error_exp_group_post change_error_control_group_post];
change_correct_pre=[change_correct_exp_group_pre change_correct_control_group_pre];
change_correct_post=[change_correct_exp_group_post change_correct_control_group_post];
var_stair_pre=[staircase_variability_exp_pre staircase_variability_control_pre];
var_stair_post=[staircase_variability_exp_post staircase_variability_control_post];
mean_stair_pre=[difficulty_exp_pre difficulty_control_pre];
mean_stair_post=[difficulty_exp_post difficulty_control_post];
performance_pre=[d_prime_exp_group_pre d_prime_control_group_pre];
performance_post=[d_prime_exp_group_post d_prime_control_group_post];
change_var_stair=[var_stair_pre-var_stair_post]';
change_mean_stair=[mean_stair_pre-mean_stair_post]';
change_performance=[performance_pre-performance_post]';


%% define dependent variables for ANOVA's
Y=[meta_d_pre' meta_d_post'];
Y2=[confirmation_pre' confirmation_post'];
Y3=[change_error_pre' change_error_post'];
Y4=[change_correct_pre' change_correct_post'];
Y5=[performance_pre' performance_post'];

%% define independent variables for  ANOVA's
session=[repmat(-1,1,40), repmat(1,1,40)];
Group=[repmat('E',length(Subjects_experimental_group),1); repmat('C',length(Subjects_control_group),1)];
subject=[repmat(1:40,1,2)];
time=[1, 2];

%% Conduct Analysis for Figures 2A-D

distance_bars=1.65;

%% ANOVA for meta-d/d'
t = table(Group,Y(:,1),Y(:,2),change_var_stair,change_mean_stair,'VariableNames',{'Group','t1','t2','variability','mean'});
rm = fitrm(t,'t1-t2 ~ Group','WithinDesign',time);
ranovatbl = ranova(rm)
% ANOVA for meta-d/d' controlling for stimulus characteristics
rm = fitrm(t,'t1-t2 ~ Group*variability*mean','WithinDesign',time);
ranovatbl = ranova(rm)

%% Figure 2A
figure(1)
hold on
pre=bar([1],[mean(meta_d_pre)],.35)
post=bar([distance_bars],[mean(meta_d_post)],.35)
errorbar([1],[mean(meta_d_pre)],[std(meta_d_pre)/sqrt(length(meta_d_pre))],'.','Color', [0 0 0], 'MarkerSize',2,'MarkerFaceColor',   [0 0 0],'LineWidth',1.5)
errorbar([distance_bars],[mean(meta_d_post)],[std(meta_d_post)/sqrt(length(meta_d_post))],'.','Color', [0 0 0], 'MarkerSize',2,'MarkerFaceColor',   [0 0 0],'LineWidth',1.5)
xlim([0.5 distance_bars+.5])
ylim([0 1.15])
ylabel('Metacognitive efficiency')
pbaspect([.65 1 1])
set(gca, 'FontSize', 14,'FontName','Arial','FontWeight','bold','box','off', 'XTick',[1, distance_bars], 'XTickLabel',{'Before training','After training'})
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5)
fix_xticklabels(gca,2,{'FontSize',14,'FontName','Arial','FontWeight','bold'});
% export_fig('D:\Train Confirmation Bias\Figure\Meta_d',  '-pdf','-nocrop', '-painters', '-transparent', [gca])


%% ANOVA for post-decision evidence integration
t2 = table(Group,Y2(:,1),Y2(:,2),change_var_stair,change_mean_stair,change_performance,'VariableNames',{'Group','t1','t2','variability','mean','performance'});
rm = fitrm(t2,'t1-t2 ~ Group','WithinDesign',time);
ranovatbl = ranova(rm)
% ANOVA for post-decision evidence integration controlling for stimulus characteristics
rm = fitrm(t2,'t1-t2 ~ Group*mean*performance','WithinDesign',time);
ranovatbl = ranova(rm)

%% Figure 2B
figure(2)
hold on
pre=bar([1],[mean(confirmation_pre)],.35)
post=bar([distance_bars],[mean(confirmation_post)],.35)
errorbar([1],[mean(confirmation_pre)],[std(confirmation_pre)/sqrt(length(confirmation_pre))],'.','Color', [0 0 0], 'MarkerSize',2,'MarkerFaceColor',   [0 0 0],'LineWidth',1.5)
errorbar([distance_bars],[mean(confirmation_post)],[std(confirmation_post)/sqrt(length(confirmation_post))],'.','Color', [0 0 0], 'MarkerSize',2,'MarkerFaceColor',   [0 0 0],'LineWidth',1.5)
xlim([0.5 distance_bars+.5])
ylim([0 .7])
ylabel('Post-decision integration')
pbaspect([.65 1 1])
set(gca, 'FontSize', 14,'FontName','Arial','FontWeight','bold','box','off', 'XTick',[1, distance_bars], 'XTickLabel',{'Before training','After training'})
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5)
fix_xticklabels(gca,2,{'FontSize',14,'FontName','Arial','FontWeight','bold'});
% export_fig('D:\Train Confirmation Bias\Figure\Post_integration',  '-pdf','-nocrop', '-painters', '-transparent', [gca])


%% ANOVA for revisions of errors
t3 = table(Group,Y3(:,1),Y3(:,2),'VariableNames',{'Group','t1','t2'});
rm = fitrm(t3,'t1-t2 ~ Group','WithinDesign',time);
ranovatbl = ranova(rm)

%% Figure 2C lower panel 
figure(3)
hold on
pre=bar([1],[mean(change_error_pre)],.35)
post=bar([distance_bars],[mean(change_error_post)],.35)
errorbar([1],[mean(change_error_pre)],[std(change_error_pre)/sqrt(length(change_error_pre))],'.','Color', [0 0 0], 'MarkerSize',2,'MarkerFaceColor',   [0 0 0],'LineWidth',1.5)
errorbar([distance_bars],[mean(change_error_post)],[std(change_error_post)/sqrt(length(change_error_post))],'.','Color', [0 0 0], 'MarkerSize',2,'MarkerFaceColor',   [0 0 0],'LineWidth',1.5)
xlim([0.5 distance_bars+.5])
ylim([0 .65])
ylabel('Revised mistakes')
pbaspect([.65 1 1])
set(gca, 'FontSize', 14,'FontName','Arial','FontWeight','bold','box','off', 'XTick',[1, distance_bars], 'XTickLabel',{'Before training','After training'})
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5)
fix_xticklabels(gca,2,{'FontSize',14,'FontName','Arial','FontWeight','bold'});
% export_fig('D:\Train Confirmation Bias\Figure\ReviseError',  '-pdf','-nocrop', '-painters', '-transparent', [gca])


%% ANOVA for revisions of correct
t4 = table(Group,Y4(:,1),Y4(:,2),'VariableNames',{'Group','t1','t2'});
rm = fitrm(t4,'t1-t2 ~ Group','WithinDesign',time);
ranovatbl = ranova(rm)

%% Figure 2C upper panel 
figure(4)
hold on
pre=bar([1],[mean(change_correct_pre)],.35)
post=bar([distance_bars],[mean(change_correct_post)],.35)
errorbar([1],[mean(change_correct_pre)],[std(change_correct_pre)/sqrt(length(change_correct_pre))],'.','Color', [0 0 0], 'MarkerSize',2,'MarkerFaceColor',   [0 0 0],'LineWidth',1.5)
errorbar([distance_bars],[mean(change_correct_post)],[std(change_correct_post)/sqrt(length(change_correct_post))],'.','Color', [0 0 0], 'MarkerSize',2,'MarkerFaceColor',   [0 0 0],'LineWidth',1.5)
xlim([0.5 distance_bars+.5])
ylim([0 .65])
ylabel('Revise correct')
pbaspect([.65 1 1])
set(gca, 'FontSize', 14,'FontName','Arial','FontWeight','bold','box','off', 'XTick',[1, distance_bars], 'XTickLabel',{'Before training','After training'})
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5)
fix_xticklabels(gca,2,{'FontSize',14,'FontName','Arial','FontWeight','bold'});
% export_fig('D:\Train Confirmation Bias\Figure\ReviseCorrect',  '-pdf','-nocrop', '-painters', '-transparent', [gca])


%% ANOVA for initial performance
t5 = table(Group,Y5(:,1),Y5(:,2),'VariableNames',{'Group','t1','t2'});
rm = fitrm(t5,'t1-t2 ~ Group','WithinDesign',time);
ranovatbl = ranova(rm)

%% Figure 2D
figure(5)
hold on
pre=bar([1],[mean(performance_pre)],.35)
post=bar([distance_bars],[mean(performance_post)],.35)
errorbar([1],[mean(performance_pre)],[std(performance_pre)/sqrt(length(performance_pre))],'.','Color', [0 0 0], 'MarkerSize',2,'MarkerFaceColor',   [0 0 0],'LineWidth',1.5)
errorbar([distance_bars],[mean(performance_post)],[std(performance_post)/sqrt(length(performance_post))],'.','Color', [0 0 0], 'MarkerSize',2,'MarkerFaceColor',   [0 0 0],'LineWidth',1.5)
xlim([0.5 distance_bars+.5])
ylim([0 1.4])
ylabel('D-prime (initial decision)')
pbaspect([.65 1 1])
set(gca, 'FontSize', 14,'FontName','Arial','FontWeight','bold','box','off', 'XTick',[1, distance_bars], 'XTickLabel',{'Before training','After training'})
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5)
fix_xticklabels(gca,2,{'FontSize',14,'FontName','Arial','FontWeight','bold'});
% export_fig('D:\Train Confirmation Bias\Figure\D_prime',  '-pdf','-nocrop', '-painters', '-transparent', [gca])



%% calculate individual differences in the training effect and correlate these
change_M_ratio=[meta_d_post-meta_d_pre]';
change_confirmation_bias=[confirmation_post-confirmation_pre]';

fitlm([zscore(change_M_ratio)], zscore(change_confirmation_bias), 'RobustOpts', 'on')
fitlm([zscore(change_M_ratio) zscore(change_var_stair) zscore(change_mean_stair) zscore(change_performance)], zscore(change_confirmation_bias), 'RobustOpts', 'on')

fitlm([zscore(change_var_stair) zscore(change_mean_stair) zscore(change_performance)], zscore(change_confirmation_bias), 'RobustOpts', 'on')
fitlm([zscore(change_mean_stair)], zscore(change_M_ratio), 'RobustOpts', 'on')

%% Figure 2E
figure(6)
hold on
scatter(change_M_ratio, change_confirmation_bias, 'o', 'filled')
lsline
 xlabel('Training effeect metacognition')
ylabel('Training effect post-decision integration')
set(findall(gca, 'Type', 'Line'),'LineWidth',3)
set(gca, 'FontSize', 16,'FontName','Arial','FontWeight','bold','box','off')
plot([0 0], [-.6 1.5], 'k:')
plot([-.6 1.5], [0 0], 'k:')
xlim([-.6 1.5])
ylim([-.6 1.5])

% export_fig('D:\Train Confirmation Bias\Figure\Correlation',  '-pdf','-nocrop', '-painters', '-transparent', [gca])



%% Mediation analysis

training_mediation=[repmat([0],1,40), repmat([1],1,40)]; % independent variable = training session
meta_mediaton=[meta_d_pre, meta_d_post]-repmat(meta_d_pre, 1, 2);  % mediator= change in metacognition
confirm_mediaton=[confirmation_pre, confirmation_post]-repmat(confirmation_pre, 1, 2); % dependent variable= change in confirmatoin bias

% covariates for mediation
var_stair_mediaton=[var_stair_pre, var_stair_post]-repmat(var_stair_pre, 1, 2); 
mean_stair_mediaton=[mean_stair_pre, mean_stair_post]-repmat(mean_stair_pre, 1, 2); 
performance_mediation=[performance_pre, performance_post]-repmat(performance_pre, 1, 2); 

[paths, stats] = mediation(training_mediation', confirm_mediaton', meta_mediaton','covs',[var_stair_mediaton', mean_stair_mediaton', performance_mediation'] ,'robust','verbose','boot', 'bootsamples', 10000)

[paths, stats] = mediation(training_mediation', meta_mediaton',confirm_mediaton','covs',[var_stair_mediaton', mean_stair_mediaton', performance_mediation'] ,'robust','verbose','boot', 'bootsamples', 10000)

