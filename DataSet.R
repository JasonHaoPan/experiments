

adenocarcinoma_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/adenocarcinoma.csv", header = TRUE, sep = ",")
adenocarcinoma.cluster = adenocarcinoma_data[,1]+1
adenocarcinoma_data <- adenocarcinoma_data[,-1]


Brain_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/Brain.csv", header = TRUE, sep = ",")
Brain.cluster = Brain_data[,1]+1
Brain_data <- Brain_data[,-1]


Brain_Tumor1_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/Brain_Tumor1.csv", header = TRUE, sep = ",")
Brain_Tumor1.cluster = Brain_Tumor1_data[,1]+1
Brain_Tumor1_data <- Brain_Tumor1_data[,-1]


Breast2classes_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/Breast2classes.csv", header = TRUE, sep = ",")
Breast2classes.cluster = Breast2classes_data[,1]+1
Breast2classes_data <- Breast2classes_data[,-1]


Breast3classes_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/Breast3classes.csv", header = TRUE, sep = ",")
Breast3classes.cluster = Breast3classes_data[,1]+1
Breast3classes_data <- Breast3classes_data[,-1]



DLBCL_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/DLBCL.csv", header = TRUE, sep = ",")
DLBCL.cluster = DLBCL_data[,1]+1
DLBCL_data <- DLBCL_data[,-1]



Leukemia1_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/Leukemia1.csv", header = TRUE, sep = ",")
Leukemia1.cluster = Leukemia1_data[,1]+1
Leukemia1_data <- Leukemia1_data[,-1]


Lung_cancer_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/Lung_cancer.csv", header = TRUE, sep = ",")
Lung_cancer.cluster = Lung_cancer_data[,1]+1
Lung_cancer_data <- Lung_cancer_data[,-1]


Lymphoma_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/Lymphoma.csv", header = TRUE, sep = ",")
Lymphoma.cluster = Lymphoma_data[,1]+1
Lymphoma_data <- Lymphoma_data[,-1]


NCI_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/NCI.csv", header = TRUE, sep = ",")
NCI.cluster = NCI_data[,1]+1
NCI_data <- NCI_data[,-1]


Prostate_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/Prostate.csv", header = TRUE, sep = ",")
Prostate.cluster = Prostate_data[,1]+1
Prostate_data <- Prostate_data[,-1]


Prostate_Tumor_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/Prostate_Tumor.csv", header = TRUE, sep = ",")
Prostate_Tumor.cluster = Prostate_Tumor_data[,1]+1
Prostate_Tumor_data <- Prostate_Tumor_data[,-1]


SRBCT_txt_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/SRBCT_txt.csv", header = TRUE, sep = ",")
SRBCT_txt.cluster = SRBCT_txt_data[,1]+1
SRBCT_txt_data <- SRBCT_txt_data[,-1]


Tumors9_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/Tumors9.csv", header = TRUE, sep = ",")
Tumors9.cluster = Tumors9_data[,1]+1
Tumors9_data <- Tumors9_data[,-1]




CNS_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/CNS.csv", header = TRUE, sep = ",")
CNS.cluster = CNS_data[,7130]+1
CNS_data <- CNS_data[,-7130]


Colon_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/Colon.csv", header = TRUE, sep = ",")
Colon.cluster = as.numeric(Colon_data[,2001])
Colon_data <- Colon_data[,-2001]



Leukemia_3c_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/Leukemia_3c.csv", header = TRUE, sep = ",")
Leukemia_3c.cluster = as.numeric(Leukemia_3c_data[,7130])
Leukemia_3c_data <- Leukemia_3c_data[,-7130]


Leukemia_4c_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/Leukemia_4c.csv", header = TRUE, sep = ",")
Leukemia_4c.cluster = as.numeric(Leukemia_4c_data[,7130])
Leukemia_4c_data <- Leukemia_4c_data[,-7130]


LeukemiaTXT_data <- read.csv("C:/Users/user/Desktop/Yixin/SZU/1_Dataset/LeukemiaTXT.csv", header = TRUE, sep = ",")
LeukemiaTXT.cluster = as.numeric(LeukemiaTXT_data[,1])
LeukemiaTXT_data <- LeukemiaTXT_data[,-1]
