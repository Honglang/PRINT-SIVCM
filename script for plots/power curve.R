library(ggplot2)

nsim=500
pv1<-rep(0,nsim)
pv2<-rep(0,nsim)
bv1<-rep(0,nsim)
bv2<-rep(0,nsim)
bv3<-rep(0,nsim)
bv4<-rep(0,nsim)
bv5<-rep(0,nsim)
nrandom=1
for(ii in c(1:1000)){
  if(file.exists(paste('C:/Users/cuiyi/Desktop/results/result_',ii,'.RData',sep = ''))){
    load(paste('C:/Users/cuiyi/Desktop/results/result_',ii,'.RData',sep = ''))
    bv4[nrandom]=beta_RKHS
    bv5[nrandom]=beta_WC
    pv1[nrandom]=p1
    pv2[nrandom]=p2
    nrandom=nrandom+1
    if(nrandom>nsim){
      break
    }
  }
}
sum(pv1[abs(bv4-0.4)<0.5]<0.05)/sum(abs(bv4-0.4)<0.5)
sum(pv2[abs(bv5-0.4)<0.5]<0.05)/sum(abs(bv5-0.4)<0.5)

library(ggplot2)
library(gridExtra)

x1=rep(c(0,0.1,0.2,0.3),4)
y1=c(0.082,0.452,0.88,0.988,0.07,0.634,0.99,1,0.064,0.356,0.696,0.836,0.07,0.522,0.874,0.94)
group=rep(c("PMRC(n=1000)","PMRC(n=2000)","PCAL(n=1000)","PCAL(n=2000)"),each=4)
data1=data.frame(x1=x1,y1=y1,group=group)

x2=rep(c(0,0.1,0.2,0.3),4)
y2=c(0.074,0.33,0.74,0.946,0.068,0.49,0.94,1,0.064,0.314,0.638,0.798,0.07,0.462,0.808,0.924)
data2=data.frame(x2=x2,y2=y2,group=group)

par(mfrow = c(1,2))

p1<-ggplot(data1, aes(x = x1, y = y1,fill=group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Power for Set 1", x = expression(x^T * beta), y = "Reject Rate") +
  theme_minimal()
p2<-ggplot(data2, aes(x = x2, y = y2,fill=group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Power for Set 2", x = expression(x^T * beta), y = "Reject Rate") +
  theme_minimal()

gridExtra::grid.arrange(p1, p2, ncol = 2)

#power for project 2

x1=rep(c(0,0.1,0.2,0.3),4)
y1=c(0.082,0.452,0.88,0.988,0.07,0.634,0.99,1,0.064,0.356,0.696,0.836,0.07,0.522,0.874,0.94)
maingroup=rep(rep(c("Sample Size: 1000","Sample Size: 2000"),each=4),2)
subgroup=rep(c("PMRC","PCAL"),each=8)
data1=data.frame(x1=x1,y1=y1,maingroup=maingroup,subgroup=subgroup)

x2=rep(c(0,0.1,0.2,0.3),4)
y2=c(0.074,0.33,0.74,0.946,0.068,0.49,0.94,1,0.064,0.314,0.638,0.798,0.07,0.462,0.808,0.924)
data2=data.frame(x2=x2,y2=y2,maingroup=maingroup,subgroup=subgroup)

p1<-ggplot(data1, aes(x = x1, y = y1, fill = subgroup)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, color = "red", linetype = 1, size = 0.2) + 
  facet_grid(~ maingroup, scales = "free", space = "free") + 
  labs(title = "Scenario 1", x = expression(x^T*beta), y = "Power", fill = "Methods") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"))

p2<-ggplot(data2, aes(x = x2, y = y2, fill = subgroup)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, color = "red", linetype = 1, size = 0.2) + 
  facet_grid(~ maingroup, scales = "free", space = "free") + 
  labs(title = "Scenario 2", x = expression(x^T*beta), y = "Power", fill = "Methods") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"))
gridExtra::grid.arrange(p1, p2, nrow = 2)

x3=rep(c(0,0.1,0.2,0.3),4)
y3=c(0.052,0.32,0.7,0.94,0.054,0.458,0.922,0.998,0.064,0.216,0.394,0.804,0.08,0.292,0.716,0.942)
data3=data.frame(x3=x3,y3=y3,maingroup=maingroup,subgroup=subgroup)

p3<-ggplot(data3, aes(x = x3, y = y3, fill = subgroup)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, color = "red", linetype = 1, size = 0.2) + 
  facet_grid(~ maingroup, scales = "free", space = "free") + 
  labs(title = "Scenario 3", x = expression(x^T*beta), y = "Power", fill = "Methods") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"))

gridExtra::grid.arrange(p1, p2,p3, nrow = 3)
############################################
#project 2 power curve with gamma=1

x4=rep(c(0,0.1,0.2,0.3),4)
y4=c(0.07,0.436, 0.856 , 0.992,0.056,0.618 , 0.976 , 1,0.068,0.456 , 0.874, 0.992,0.06,0.682, 0.988, 1)
maingroup=rep(rep(c("Sample Size: 1000","Sample Size: 2000"),each=4),2)
subgroup=rep(c("PMRC","PCAL"),each=8)
data4=data.frame(x4=x4,y4=y4,maingroup=maingroup,subgroup=subgroup)

x5=rep(c(0,0.1,0.2,0.3),4)
y5=c(0.07,0.36, 0.758, 0.934,0.074,0.494, 0.95, 0.996,0.068,0.406,0.806, 0.978,0.064,0.594, 0.966, 1)
data5=data.frame(x5=x5,y5=y5,maingroup=maingroup,subgroup=subgroup)

p4<-ggplot(data4, aes(x = x4, y = y4, fill = subgroup)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, color = "red", linetype = 1, size = 0.2) + 
  facet_grid(~ maingroup, scales = "free", space = "free") + 
  labs(title = "Scenario 1", x = expression(x^T*beta), y = "Power", fill = "Methods") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"))

p5<-ggplot(data5, aes(x = x5, y = y5, fill = subgroup)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, color = "red", linetype = 1, size = 0.2) + 
  facet_grid(~ maingroup, scales = "free", space = "free") + 
  labs(title = "Scenario 2", x = expression(x^T*beta), y = "Power", fill = "Methods") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"))
gridExtra::grid.arrange(p4, p5, nrow = 2)

x7=rep(c(0,0.1,0.2,0.3),4)
y7=c(0.072,0.34, 0.752, 0.930,0.05,0.482, 0.93, 0.998,0.07,0.268, 0.506, 0.770,0.082,0.336, 0.706, 0.926)
data7=data.frame(x7=x7,y7=y7,maingroup=maingroup,subgroup=subgroup)

p7<-ggplot(data7, aes(x = x7, y = y7, fill = subgroup)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, color = "red", linetype = 1, size = 0.2) + 
  facet_grid(~ maingroup, scales = "free", space = "free") + 
  labs(title = "Scenario 3", x = expression(x^T*beta), y = "Power", fill = "Methods") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"))

gridExtra::grid.arrange(p4, p5,p7, nrow = 3)
#############################################

x=rep(c(0,0.1,0.2,0.3),8)
y=c(0.082,0.452,0.88,0.988,0.07,0.634,0.99,1,0.064,0.356,0.696,0.836,0.07,0.522,0.874,0.94,0.074,0.33,0.74,0.946,0.068,0.49,0.94,1,0.064,0.314,0.638,0.798,0.07,0.462,0.808,0.924)
maingroup=rep(rep(c("Sample Size n=1000","Sample Size n=2000"),each=4),4)
subgroup=c(rep(c("PMRC(Set 1)","PCAL(Set 1)"),each=8),rep(c("PMRC(Set 2)","PCAL(Set 2)"),each=8))
data=data.frame(x=x,y=y,maingroup=maingroup,subgroup=subgroup)

ggplot(data, aes(x = x, y = y, fill = subgroup)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(~ maingroup, scales = "free", space = "free") + 
  labs(title = "Power Curves", x = expression(x^T * beta), y = "Reject Rate") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"))  # Make category labels bold

#power for project 3

x1=rep(c(0,0.1,0.2,0.3),4)
y1=c(0.054,0.356, 0.818 , 0.982,0.064,0.544,0.978,1,0.058,0.396,0.848,0.988,0.06,0.588,0.986,1)
maingroup=rep(rep(c("Sample Size: 1000","Sample Size: 2000"),each=4),2)
subgroup=rep(c("MRC","CAL"),each=8)
data1=data.frame(x1=x1,y1=y1,maingroup=maingroup,subgroup=subgroup)

x2=rep(c(0,0.1,0.2,0.3),4)
y2=c(0.042,,0.284,0.7,0.912,0.066,0.44,0.942,1,0.044,0.326,0.774,0.962,0.06,0.506,0.97,1)
data2=data.frame(x2=x2,y2=y2,maingroup=maingroup,subgroup=subgroup)

p1<-ggplot(data1, aes(x = x1, y = y1, fill = subgroup)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, color = "red", linetype = 1, size = 0.2) + 
  facet_grid(~ maingroup, scales = "free", space = "free") + 
  labs(title = "Scenario 1", x = expression(x^T*beta), y = "Power", fill = "Methods") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"))  # Make category labels bold

p2<-ggplot(data2, aes(x = x2, y = y2, fill = subgroup)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, color = "red", linetype = 1, size = 0.2) + 
  facet_grid(~ maingroup, scales = "free", space = "free") + 
  labs(title = "Scenario 2", x = expression(x^T*beta), y = "Power", fill = "Methods") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"))  # Make category labels bold

gridExtra::grid.arrange(p1, p2, nrow = 2)


x3=rep(c(0,0.1,0.2,0.3),4)
y3=c(0.056,0.29,0.678,0.892,0.058,0.416,0.91,0.996,0.06,0.224,0.454,0.756,0.052,0.302,0.648,0.904)
data3=data.frame(x3=x3,y3=y3,maingroup=maingroup,subgroup=subgroup)

p3<-ggplot(data3, aes(x = x3, y = y3, fill = subgroup)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.05, color = "red", linetype = 1, size = 0.2) + 
  facet_grid(~ maingroup, scales = "free", space = "free") + 
  labs(title = "Scenario 3", x = expression(x^T*beta), y = "Power", fill = "Methods") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"))

gridExtra::grid.arrange(p1, p2,p3, nrow = 3)
