#### world map plots
[use R to plot world map](https://bitbucket.org/MAVERICLab/global-rna-virus-evolution-2021/src/master/Figure4/mapper.R)
	
	#!/usr/bin/R
	library("maps")
	library("tidyverse")
	#library("ggtree")
	library("ggplot2")
	library("ggimage")
	
	world = map_data("world")
	world_map = ggplot(world, fill = T, aes(long, lat, group = group)) +
  	geom_polygon(fill = "#ABABAD") +
  	coord_cartesian(xlim = c(-170, 175), ylim = c(-65, 80)) +
  	#  coord_cartesian(xlim = c(-140, 80),ylim = c(-65, 50)) +
  	labs(x = "Longitude", y = 'Latitude') +
  	theme(
    	legend.position = "none",
    	plot.title = element_text(face = "bold", size = 14),
    	axis.title = element_text(face = 'bold', size = 12)
  	)
	plot(world_map)
