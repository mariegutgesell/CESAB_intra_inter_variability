##phylopic silhouettes

#Pkgs
library(tidyverse)
library(RColorBrewer)
library(beepr)
library(rphylopic)
library(rsvg)
library(grid)

##### Code #####
#Create data frame for img coords
df <- data.frame(x = c(2, 4), y = c(10, 20))

##Salvelinus fontinalis
ggplot(df) +
  geom_phylopic(aes(x = x, y = y, name = "Salvelinus fontinalis"), fill = "black",color="black", size = 5) +
  coord_cartesian(xlim = c(1,3), ylim = c(5, 30)) +
  theme_void()

##Catostomus commersonii
ggplot(df) +
  geom_phylopic(aes(x = x, y = y, name = "Catostomus commersonii"), fill = "black",color="black", size = 5) +
  coord_cartesian(xlim = c(1,3), ylim = c(5, 30)) +
  theme_void()

##Gambusia holbrooki
ggplot(df) +
  geom_phylopic(aes(x = x, y = y, name = "Gambusia holbrooki"), fill = "black",color="black", size = 5) +
  coord_cartesian(xlim = c(1,3), ylim = c(5, 30)) +
  theme_void()

##Leiopotherapon unicolor
ggplot(df) +
  geom_phylopic(aes(x = x, y = y, name = "Leiopotherapon unicolor"), fill = "black",color="black", size = 5) +
  coord_cartesian(xlim = c(1,3), ylim = c(5, 30)) +
  theme_void()
