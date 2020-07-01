d1 <- load(url("https://github.com/kcampbel/neocovid-app/raw/master/app/data/forApp.Rda"))
write.csv(predictions,'./data/ctl/predictions.csv')
