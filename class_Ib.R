#TASK: making folders
subfolders <- c("raw_data", "clean_data", "scripts", "results", "tests", "plots", "workspaces")
for (name in subfolders) {
  dir.create(name)
}

#TASK: Patient info management
library(readr)
patient_info <- read_csv("raw_data/patient_info.csv")
View(patient_info)
summary(patient_info)
str(patient_info)

patient_info$gender <- as.factor(patient_info$gender)
patient_info$diagnosis <- as.factor(patient_info$diagnosis)
patient_info$diagnosis <- factor(patient_info$diagnosis, levels = c("Normal", "Cancer"))
patient_info$smoker <- as.factor(patient_info$smoker)

patient_info$smoker_bin <- ifelse(patient_info$smoker == "Yes", 1, 0)

View(patient_info)
summary(patient_info)
str(patient_info)

write.csv(patient_info, "./clean_data/patient_info_clean.csv")
save.image(file = "./workspaces/Ayisha_Class_Ib_Assignment.RData")
