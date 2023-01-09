library("tidyverse")
library("patchwork")

file_name <- "output.csv"

the_data <- read_delim(file = file_name,delim = ";")



plot_size0 <- ggplot(data=the_data,
       mapping=aes(x=generation,y=mean_size0)) +
  geom_line()

plot_size1 <- ggplot(data=the_data
                     ,mapping=aes(x=generation,y=mean_size1)) +
  geom_line()

plot_size0 / plot_size1

ggsave("some_graph_file.pdf")