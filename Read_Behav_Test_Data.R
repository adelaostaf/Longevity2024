# Analyse Longevity Project data

library("tidyverse")
library("broom")
library("ggplot2")
library("gghalves")
library("raincloudplots")

# First let's get all the ingredients together to flexibly read in the data

# First, let's set the path to the data; you may have to change this path to the directory on your local machine;
# If you get "Error in setwd(curfldr) : cannot change working directory" then you didn't set the path correctly 
# ... remember that depending on whether you run this on Windows, Linux, Mac, you may have to use different conventions for setting the path; the below is for a windows machine;
# ... if in doubt 
datpth = "T:\\longevity_2024\\data\\beh_data" 
grp = c("short_group", "long_group")
# Select which group you want to analyse
datsrc = c("RatingData", "SyncData", "TestData")

alldat=data.frame()
for (g in 1:length(grp)) {
  curfldr = paste(datpth, "\\", grp[g], "\\", datsrc[3], sep = "") # The argument is needed at the end to avoid having white spaces
  setwd(curfldr)
  subs <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
  for (i in 1:length(subs)) {
    cursubfldr = paste(curfldr, "\\", subs[i], sep = "")
    setwd(cursubfldr)
    blcks <- list.files(pattern="block",path = ".", full.names = FALSE, recursive = FALSE)
    for (ii in 1:length(blcks)) {
      tmpdat <- read_csv(blcks[ii])
      tmpdat <- data.frame(tmpdat, PiD=rep(subs[i], times=length(tmpdat))) 
      tmpdat$testrespMat_.2 <- as.character(tmpdat$testrespMat_.2) # use $ signs to addres specific data frames
      tmpdat$Grp = c(g) # add group info
      alldat <- bind_rows(alldat, tmpdat)
      
    }
    print(i)
  }
}

# Get correct (hit) trials and store in new variable
alldat$Hit <- ifelse(alldat$testrespMat_15 == alldat$testrespMat_14, 1, 0)
alldat <- rename(alldat, "PNum" = `testrespMat_.1`, "PhaseCon" = `testrespMat_.8`)
alldat <- unite(alldat, combined, PhaseCon, sep = "-")
alldat <- rename(alldat, "Cond" = combined)

# Get rid of unnecessary variables
dat_long <- alldat %>% select("PiD", "PNum", "Cond", "Hit", "Grp")

# combine to get accuracy for each participant, group, and condition
dat_shrt <- dat_long %>% 
  group_by(PNum, Cond, Grp) %>% 
  summarise("Acc" = mean(Hit))

# get overall memory performance to check whether participants are above chance level
mn_mem_perfor <- dat_long %>% 
  group_by(PNum, Grp) %>%
  summarise("Acc" = mean(Hit))

# Calculate chance level using a binomial test
# Number of trials
n <- 32
num_successes <-  13# Insert the number of 1s you observed in your experiment
# Probability of success (chance)
p <- 0.25
# Perform binomial test
binom_test_result <- binom.test(num_successes, n, p = p, alternative = "greater")

guessing_thrshld = num_successes / n

# Plot raincloud plots Accuracies per group; left immediate, right long delay group
p1 <- ggplot(mn_mem_perfor, aes(x = Grp, y = Acc, fill = Grp, colour = Grp, group = Grp)) +
  geom_half_violin(data = dat_shrt, position = position_nudge(x = -0.35, y = 0), adjust = 1) +
  geom_point(position = position_jitter(width = 0.15), size = 5) +
  geom_hline(yintercept = guessing_thrshld, linetype = "dashed", color = "red") +
#  geom_boxplot(aes(x = CondNum - 0.25, y = Acc, fill = CondNum, colour = CondNum), outlier.shape = NA, alpha = 0.3, width = 0.1) +
#  geom_line(data = sum_dat, aes(x = CondNum - 0.25, y = mnAcc, group=grp$grp, colour = grp$grp), linetype=2) +
  ylab("Accuracy") + xlab("Condition") + 
  ggtitle("Showing average memory performance for long and short delay group - red line indicates guess level cut off")
print(p1)



