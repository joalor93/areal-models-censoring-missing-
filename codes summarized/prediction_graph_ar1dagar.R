
######################################################
########## Lines to install the auxiliarcpp package##
####################################################


#setwd("C:/Users/ALEJANDRO/OneDrive/Documentos/Alejo/postdoc20242/censored_network_dagar/datasets/Beijing/codes sumarized")


library(mvtnorm)
library(invgamma)
library(tmvtnorm)
library(data.table)
library(auxiliarcpp)
library(tidyverse)

setwd("C:/Users/ALEJANDRO/OneDrive/Documentos/Alejo/postdoc20242/censored_network_dagar/datasets/Beijing/codes sumarized/")
data = fread(paste0("Estaciones_unidas_rushhour.csv"))
adj_matrix = readRDS("adj_matrix.rds")

###################################################
##### Preparing the dataset #######################
##################################################

#### complete adjacency matrix (CAR model)
adj_matcom = adj_matrix
adj_matcom = (adj_matrix + t(adj_matrix))


datapre = data

regions = c("Aotizhongxin", "Changping", "Dingling", "Dongsi",
            "Guanyuan", "Gucheng", "Huairou", "Nongzhanguan",
            "Shunyi", "Tiantan","Wanliu", "Wanshouxigong")

datapre2 = datapre%>%group_by("station")%>%
  mutate (aggtemp = if_else(is.na(TEMP),mean(TEMP,na.rm = TRUE),TEMP),
          aggrain = if_else(is.na(aggrain),mean(aggrain,na.rm = TRUE),aggrain),
          aggsumrain = if_else(is.na(aggsumrain),mean(aggsumrain,na.rm = TRUE),aggsumrain),
          aggpres = if_else(is.na(PRES),mean(PRES,na.rm = TRUE),PRES),
          aggwsmp = if_else(is.na(WSPM),mean(aggwsmp,na.rm = TRUE),WSPM))%>%ungroup()


######observed data #####
datapre2obs = datapre2%>%filter(data>= "2016-11-25" & data<="2017-02-25")


##### complete data ########
datapre2 = datapre2%>%filter(data>= "2016-11-25" & data<="2017-02-28")
datapre2 = datapre2%>%mutate(indipred = if_else(data>= "2017-02-25" & data<="2017-02-28",1,0))




xtot = cbind (1,datapre2$aggtemp - mean(datapre2$aggtemp),datapre2$aggwsmp - mean(datapre2$aggwsmp), (datapre2$aggpres - mean(datapre2$aggpres))/100)

#####################################################


###################################################
##### Including predictions on the dataset #######
#################################################

resar1dagar = readRDS("cad1ar1dagar_Beijing.rds")



pred_res_ar1_dagar = readRDS("results_pred_ar1dagar.rds")


datapreinterest = datapre2%>%select(data,station,CO)
datapreinterest = datapreinterest%>%filter(data>= "2016-11-25" & data<="2017-02-28")

datapreinterest = datapreinterest%>%mutate(indipred = if_else(data>= "2017-02-25" & data<="2017-02-28",1,0))
datapreinterest = datapreinterest%>%group_by("station")%>%mutate (log_co = log(CO))

datapreinterest = datapreinterest%>%mutate(indicens3 = as.numeric(is.na(log_co) & indipred ==0))

datapreinterest$log_co[is.na(datapreinterest$log_co) & datapreinterest$indipred ==1] = mean(datapreinterest$log_co, na.rm = TRUE)


####### common predictions #####
pred_common = function(x,beta){
  return(x%*%beta)
}

predicommon = apply(resar1dagar$betaF,1,pred_common,x = xtot)
ypred2 =  apply(predicommon,1,mean)
prelim_pred = apply(predicommon,1,quantile,probs = c(0.025,0.5,0.975))
prelim_pred = t(prelim_pred)


datapreinterest = datapreinterest%>%mutate(ypred = ypred2, ypred_25 = prelim_pred[,1], ypred_975 = prelim_pred[,3])


pred_dagar = apply(pred_res_ar1_dagar,1,quantile,probs = c(0.025,0.5,0.975))
pred_dagar = t(pred_dagar)

datapreinterest$ypred[datapreinterest$indipred==1] = pred_dagar[,2]
ycens = apply(resar1dagar$ycensF,2,mean)
datapreinterest$ypred[datapreinterest$indicens3==1 ] = ycens

datapreinterest$ypred_25[datapreinterest$indipred==1] = pred_dagar[,1]
datapreinterest$ypred_975[datapreinterest$indipred==1] = pred_dagar[,3]

datapreinterestfiltered = datapreinterest%>%filter(data>="2017-02-01")

datapreinterestfiltered = datapreinterestfiltered%>%select(data,station, log_co,ypred_25,ypred,ypred_975)



#Force English locale for month/day labels####

Sys.setlocale("LC_TIME", "C")   #"C" = English POSIX locale, works on all OS



#### auxiliar dataset to use ggplot ########
df_all <- datapreinterestfiltered %>%
  filter(as.POSIXct(data) >= as.POSIXct("2017-02-18 00:00:00")) %>%
  mutate(
    lower = pmin(ypred_25, ypred_975),
    upper = pmax(ypred_25, ypred_975),
    lower_show = if_else(as.Date(data) >= as.Date("2017-02-25"), lower, NA_real_),
    upper_show = if_else(as.Date(data) >= as.Date("2017-02-25"), upper,  NA_real_)
  )

# Long format for lines
lines_long_all <- df_all %>%
  select(station, data, log_co, ypred) %>%
  pivot_longer(c(log_co, ypred), names_to = "Series", values_to = "Value") %>%
  mutate(Series = recode(Series,
                         log_co = "Observed Log(CO)",
                         ypred  = "Predicted Log(CO)"))

#Find where confidence bands start
pred_start <- df_all %>%
  filter(!is.na(lower_show) | !is.na(upper_show)) %>%
  summarise(start_date = min(data, na.rm = TRUE)) %>%
  pull(start_date)

cat("Prediction starts at:", pred_start, "\n")

################################################
######### PREDICTION GRAPH ####################
##############################################

p <- ggplot() +
  # confidence band
  geom_ribbon(
    data = df_all,
    aes(x = data, ymin = lower_show, ymax = upper_show),
    fill = "gray60", alpha = 0.45, na.rm = TRUE, show.legend = FALSE
  ) +
  # predicted line (gray)
  geom_line(
    data = filter(lines_long_all, Series == "Predicted Log(CO)"),
    aes(x = data, y = Value, colour = Series),
    linewidth = 1.2
  ) +
  # observed line (black)
  geom_line(
    data = filter(lines_long_all, Series == "Observed Log(CO)"),
    aes(x = data, y = Value, colour = Series),
    linewidth = 1.2
  ) +
  geom_point(
    data = filter(lines_long_all, Series == "Observed Log(CO)"),
    aes(x = data, y = Value, colour = Series),
    size = 2
  ) +
  # vertical dotted line where confidence band starts
  geom_vline(
    xintercept = pred_start,
    linetype = "dotted", colour = "gray30", linewidth = 0.8,
    show.legend = FALSE
  ) +
  # grayscale legend
  scale_colour_manual(
    name = "",
    values = c("Observed Log(CO)" = "black",
               "Predicted Log(CO)" = "gray30")
  ) +
  facet_wrap(~station, ncol = 4, scales = "free_y") +
  labs(
    title = "",
    x = "Date",
    y = "Log(CO)"
  ) +
  scale_x_datetime(date_breaks = "3 days", date_labels = "%b %d") +
  theme_minimal(base_size = 18) +  # increased global base font
  theme(
    plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    legend.key = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 18)
  ) +
  coord_cartesian(ylim = c(2, 9))


p


#ggsave(
 # filename = "pred_beijing.pdf",
  #plot = p,
  #device = cairo_pdf,
  #width = 12, height = 10, units = "in",
  #dpi = 600,                # applies when raster layers exist
  #fallback_resolution = 600 # ensures any raster part (geom_ribbon alpha) is 600 dpi
#)


