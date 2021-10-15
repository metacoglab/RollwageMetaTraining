%% Load data 
pathToData = '~/Documents/Data/Max EEG/Github/data/'
load([pathToData 'Climate/Climate_metacognition2.mat']) % data from perceptual decision-making and climate change knowledge

%% Data measuring climate change attitudes
% load each file separately (as the script otherwise crashes)
% change the variable "QuestionKey" to a "text" variable in the dropdown
% menu, hit import selection
uiopen('Climate_belief_control_initial.csv')
uiopen('Climate_belief_control_final.csv')
uiopen('Climate_belief_exp_initial.csv')
uiopen('Climate_belief_exp_final.csv')


%% assign data of group 2 to the variables
data=Climatebeliefcontrolinitial;
data2=Climatebeliefcontrolfinal;

a=meta_table2.subject_EEG(1:20); % identify the IDs of interest

Responses_Participant=[];
group=[]
%% define number of questions and how they are coded
Number_of_items=4;
answer_coding=[1,1,-1,-1];


for i=1:length(a)
%% get answers of specific participant
b=find(data.ParticipantPublicID==a(i)); 
questions=data.QuestionKey(b);
response=data.Response(b);

%% Loop over answers for this participant to get their initial answers 
for question=1:4
index_question=strfind(questions, ['ClimateBeliefInit' num2str(question) '-quantised']);
ind_question_use=[];
for entries=1:length(index_question)
    if index_question{entries}==1
       ind_question_use=[ind_question_use, entries];
    end
end
decision_initial(question)=response(ind_question_use(1));

end
%% Loop over answers for this participant to get their final answers 

b=find(data2.ParticipantPublicID==a(i));
questions=data2.QuestionKey(b);
response=data2.Response(b);

for question=1:4
index_question=strfind(questions, ['ClimateBeliefInit' num2str(question) '-quantised']);
ind_question_use=[];
for entries=1:length(index_question)
    if index_question{entries}==1
       ind_question_use=[ind_question_use, entries];
    end
end
decision_final(question)=response(ind_question_use(1));
end
beliefs_initial(i,:)=decision_initial;
beliefs_final(i,:)=decision_final;

%% calulate belief shift for this participant
belief_shift=decision_final-decision_initial; %directed belief shift
belief_shift_sub(i,:)=belief_shift;
initial_attitude(i)=answer_coding*decision_initial'; % higher values indicate more severity
final_attitude(i)=answer_coding*decision_final'; % higher values indicate more severity
directed_belief_shift(i)=answer_coding*belief_shift'; % higher values indicate that people think more severe about climate change after training
absolute_belief_update(i)=sum(abs(belief_shift)); % absolute belief shift
group(i)=0;
end

%% assign data of group 1 to the variables

data3=Climatebeliefexpinitial;
data4=Climatebeliefexpfinal;

z=meta_table2.subject_EEG(21:40); % identify the IDs of interest

for f=1:length(z)
k=i+f % have an index for both samples 
%% get answers of specific participant
b=find(data3.ParticipantPublicID==z(f)); 
questions=data3.QuestionKey(b);
response=data3.Response(b);

%% Loop over answers for this participant to get their initial answers 
for question=1:4
index_question=strfind(questions, ['ClimateBeliefInit' num2str(question) '-quantised']);
ind_question_use=[];
for entries=1:length(index_question)
    if index_question{entries}==1
       ind_question_use=[ind_question_use, entries];
    end
end
decision_initial(question)=response(ind_question_use(1));
end

b=find(data4.ParticipantPublicID==z(f));
questions=data4.QuestionKey(b);
response=data4.Response(b);

%% Loop over answers for this participant to get their final answers 
for question=1:4
index_question=strfind(questions, ['ClimateBeliefInit' num2str(question) '-quantised']);
ind_question_use=[];
for entries=1:length(index_question)
    if index_question{entries}==1
       ind_question_use=[ind_question_use, entries];
    end
end
decision_final(question)=response(ind_question_use(1));
end
beliefs_initial(k,:)=decision_initial;
beliefs_final(k,:)=decision_final;

%% calulate belief shift for this participant

belief_shift=decision_final-decision_initial; % directed belief shift
belief_shift_sub(k,:)=belief_shift;
initial_attitude(k)=answer_coding*decision_initial';
final_attitude(k)=answer_coding*decision_final';
directed_belief_shift(k)=answer_coding*belief_shift'; % higher values indicate that people think more severe about climate change after information
absolute_belief_update(k)=sum(abs(belief_shift)); % absolute belief shift
group(k)=1;
end

metacognitive_training=meta_table2.meta_post-meta_table2.meta_pre;
metacognition_post=meta_table2.meta_post; 
metacognition_pre=meta_table2.meta_pre; 

Update_belief_abs=log(absolute_belief_update+1);
climate_knowledge=meta_table2.knowledge_climate;
climate_confidence=meta_table2.confidence_climate;

directed_feedback=meta_table2.feedback; % higher values indicate that the feedback indicate more severity than expected


fitlm([zscore(metacognitive_training)],zscore(Update_belief_abs),'RobustOpts','on')
fitlm([zscore(metacognitive_training), climate_knowledge, climate_confidence],zscore(Update_belief_abs),'RobustOpts','on')
fitlm([metacognition_pre, zscore(metacognition_post), climate_knowledge, climate_confidence],zscore(Update_belief_abs),'RobustOpts','on')
fitlm([metacognition_pre, climate_knowledge, climate_confidence],zscore(Update_belief_abs),'RobustOpts','on')

fitlm([zscore(directed_feedback), zscore(metacognitive_training), zscore(directed_feedback).*zscore(metacognitive_training),climate_knowledge,climate_confidence],zscore(directed_belief_shift),'RobustOpts','on')
fitlm([zscore(directed_feedback), zscore(directed_feedback).*zscore(metacognitive_training),climate_knowledge,climate_confidence],zscore(directed_belief_shift),'RobustOpts','on')
fitlm([zscore(directed_feedback), zscore(metacognitive_training),climate_knowledge,climate_confidence],zscore(directed_belief_shift),'RobustOpts','on')


%% Figure 4A

figure(1)
hold on
scatter(metacognitive_training, Update_belief_abs, 'o', 'filled')
lsline
xlabel('Change metacognition')
ylabel('Change post-decision integration')
set(findall(gca, 'Type', 'Line'),'LineWidth',3)
set(gca, 'FontSize', 16,'FontName','Arial','FontWeight','bold','box','off')
% plot([-.3 1], [-.3 1],'k-')
plot([0 0], [-.6 2.6], 'k:')
plot([-.6 2.6], [0 0], 'k:')
xlim([-.53 1.5])
ylim([-.2 2.6])

