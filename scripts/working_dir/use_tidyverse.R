library(dplyr)

## see dplyr in firefox bookmarks for more
## see C:\Users\e.bordron\Desktop\CGH-scoring\M2_internship_Bergonie\docs\docs_I_found\dplyr_cheatsheet.pdf

x = data(mtcars)
glimpse(mtcars)

#--------- select = manipulate columns
## two ways of using a tidyverse function
select(mtcars, qsec)
mtcars %>% select(qsec)

## select a range of columns
mtcars %>% select(qsec:gear)
mtcars %>% select(2:4)
mtcars %>% select(!(mpg:am))

## use select helpers. see https://tidyselect.r-lib.org/reference/language.html
mtcars %>% select(starts_with("c"))
mtcars %>% select(starts_with("c") | ends_with("m"))


#--------- filter = manipulate rows
mtcars %>% filter(mpg<17)
