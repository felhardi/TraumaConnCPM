# demographic information 

trauma %>% select(sex) %>% count(sex)
trauma %>% select(age) %>% summarize(mean = mean(age), sd = sd(age), min = min(age), max = max(age))
trauma %>% select(edu) %>% count(edu)
trauma %>% select(total_5) %>% summarize(mean(total_5), sd(total_5))

trauma %>% select(sex, study) %>% filter(study=="SL") %>% count(sex)
trauma %>% select(age, study) %>% filter(study=="SL") %>% summarize(mean(age), sd(age), min = min(age), max = max(age))
trauma %>% select(edu, study) %>% filter(study=="SL") %>% count(edu)
trauma %>% select(total_5, study) %>% filter(study=="SL") %>% summarize(mean(total_5), sd(total_5))

trauma_raceth = read.csv("./trauma_raceth.csv")

trauma_raceth %>% count(race)
trauma_raceth %>% count(ethnicity)
trauma_raceth %>% filter(sub<200) %>% count(race)
trauma_raceth %>% filter(sub<200) %>% count(ethnicity)

cortmem_demo = read_csv("./CortMem_Intake_November 25, 2024_12.13.csv") %>% select(subj_id, age, sex, race, hispanic_latino, education)
colnames(cortmem_demo) = c("sub", "age", "sex", "race", "ethnicity", "education")
cortmem_demo = cortmem_demo[-c(1:2),]
cortmem_demo$sub = as.numeric(cortmem_demo$sub)
cortmem_sub = read.csv("./cortmem_sub.csv")
cortmem_demo = left_join(cortmem_sub, cortmem_demo)

cortmem_demo %>% summarise(mean(as.numeric(age)), sd(as.numeric(age)), min = min(as.numeric(age)), max = max(as.numeric(age)))
cortmem_demo %>% count(sex) # 1 male, 2 female
cortmem_demo %>% count(education) 
# 10 = Some college (no degree)
# 12 = Bachelor's degree
# 13 = Some graduate or professional studies (completed Bachelor's degree but not graduate degree)
# 14 = Completed Master's degree or equivalent or higher graduate degree
cortmem_demo %>% count(race) # 1=white, 2=black, 4=Asian
cortmem_demo %>% count(ethnicity) # 1=yes, 2=no
cortmem_demo$study = "cortmem"


########## stat sample comparison ###########

samp_age = rbind(cortmem_demo %>% select(sub, study, age), trauma %>% select(sub, study, age))
samp_age$samp = ifelse(samp_age$study=="cortmem", "cort", "full")
samp_age = rbind(samp_age, cbind(trauma %>% select(sub, study, age), samp="secpt") %>% filter(study=="SL"))

summary(aov(age ~ samp, samp_age))
TukeyHSD(aov(age ~ samp, samp_age))

sex = matrix(
  c(102, 68,  54, 38,  11, 16),
  nrow = 2, ncol = 3) 
chisq.test(sex) 

edu = matrix(
  c(1,1,18,39,2,42,27,40,  1,0,12,20,1,19,17,22,  0,0,0,4,0,10,5,8),
  nrow = 8, ncol = 3) 
chisq.test(edu, simulate.p.value = T) 

race = matrix(
  c(6,41,23,89,11,  3,34,12,36,7, 0,4,2,19,2),
  nrow = 5, ncol = 3) 
chisq.test(race, simulate.p.value = T)

eth = matrix(
  c(23, 147,  10, 82,  6, 21),
  nrow = 2, ncol=3)
chisq.test(eth, simulate.p.value = T)


#### validation sample 

# external sample
demo = read.csv("~/Projects/TraumaConn/cpm_shapes/QuestionnaireDataCom-Demosdd_DATA_LABELS_2025-06-03_1153.csv")
demo$ucla_a_id = demo$Subject.ID
demo = left_join(shapes_fin, demo)
demo$sex = ifelse(is.na(demo$Please.indicate.your.assigned.sex.at.birth.)==T, demo$Please.indicate.your.assigned.sex.at.birth..1, demo$Please.indicate.your.assigned.sex.at.birth.)
demo$Hispanic = ifelse(demo$Your.ethnic.group.or.race..choice.Hispanic.or.Latino. == "Checked", 1, 0)
demo$white = ifelse(demo$Your.ethnic.group.or.race..choice.Non.Hispanic.White.or.Caucasian. == "Checked", 1, 0)
demo$black = ifelse(demo$Your.ethnic.group.or.race..choice.Black.or.African.American. == "Checked", 1, 0)
demo$asian = ifelse(demo$Your.ethnic.group.or.race..choice.Asian. == "Checked", 1,0)
demo$natam = ifelse(demo$Your.ethnic.group.or.race..choice.Native.American.== "Checked", 1,0)
demo$natapi = ifelse(demo$Your.ethnic.group.or.race..choice.Native.Hawaiian.or.Pacific.Islander.== "Checked", 1,0)
demo$other = ifelse(demo$Your.ethnic.group.or.race..choice.Other. == "Checked", 1, 0)

demo %>% count(sex)
demo %>% count(What.is.the.highest.grade..or.year..of.regular.school.you.have.completed.)
demo %>% count(Hispanic)

demo$race = demo$natam + demo$asian + demo$black + demo$white + demo$natapi + demo$other
demo$multi = ifelse(demo$race > 1, 1, 0)
demo %>% filter(multi==0) %>% count(natam)
demo %>% filter(multi==0) %>% count(asian)
demo %>% filter(multi==0) %>% count(black)
demo %>% filter(multi==0) %>% count(white)
demo %>% filter(multi==0) %>% count(natapi)
demo %>% filter(multi==0) %>% count(other)

# internal sample
trauma %>% select(sex, study) %>% filter(!study=="SL") %>% count(sex)
trauma %>% select(age, study) %>% filter(!study=="SL") %>% summarize(mean(age), sd(age), min = min(age), max = max(age))
trauma %>% select(edu, study) %>% filter(!study=="SL") %>% count(edu)
trauma_raceth %>% filter(sub>200) %>% count(race)
trauma_raceth %>% filter(sub>200) %>% count(ethnicity)
