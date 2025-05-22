library(sigmoid)
library(ggplot2)
library(dplyr)
x<-seq(-4,4,by=0.01)
y1=sigmoid(x)
y2=sigmoid(x/0.1)
y3=sigmoid(x/0.01)

# Create 3 data sets
df1 <- data.frame(x = x, y = y1, group = "1")
df2 <- data.frame(x = x, y = y2, group = "0.1")
df3 <- data.frame(x = x, y = y3, group = "0.01")

# Combine into one data frame
df_all <- bind_rows(df1, df2, df3)

ggplot(df_all, aes(x = x, y = y, color = group)) +
  geom_line() +             # lines connecting points
  theme_minimal() +
  scale_color_manual(
    values = c("1" = "red", "0.1" = "blue", "0.01" = "green"),
    labels = c(
      expression(sigma[n] == 1),
      expression(sigma[n] == 0.1),
      expression(sigma[n] == 0.01)
    )
  ) +
  labs(
    title = "Sigmoid Function",
    x = "u",
    y = expression(s[n](u)),
    color = NULL
  )
