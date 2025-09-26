library(ncdf4)
library(plspm)
library(ggplot2)
library(tidyr)
# 打开 NetCDF 文件
gwr_nc <- nc_open("E:/ISIMIP/1CWatM/mean/year126_GDE_mean.nc")
gwr_nc <- ncvar_get(gwr_nc, "qr")
gwr_mean <- apply(gwr_nc,3,mean,na.rm=TRUE)

qs_nc <- nc_open("E:/ISIMIP/bios/qs/mean/year126_GDE_mean.nc")
qs_nc <- ncvar_get(qs_nc, "qs")
qs_mean <- apply(qs_nc,3,mean,na.rm=TRUE)

snd_nc <- nc_open("E:/ISIMIP/bios/snd/mean/year126_GDE_mean.nc")
snd_nc <- ncvar_get(snd_nc, "snd")
snd_mean <- apply(snd_nc,3,mean,na.rm=TRUE)

soilmoist_nc <- nc_open("E:/ISIMIP/bios/soilmoist/mean/year126_GDE_mean.nc")
soilmoist_nc <- ncvar_get(soilmoist_nc, "soilmoist")
soilmoist_mean <- apply(soilmoist_nc,3,mean,na.rm=TRUE)

tmax_nc <- nc_open("E:/ISIMIP/bios/tasmax/mean/year126_GDE_mean.nc")
tmax_nc <- ncvar_get(tmax_nc, "tasmax")
tmax_mean <- apply(tmax_nc,3,mean,na.rm=TRUE)

tmin_nc <- nc_open("E:/ISIMIP/bios/tasmin/mean/year126_GDE_mean.nc")
tmin_nc <- ncvar_get(tmin_nc, "tasmin")
tmin_mean <- apply(tmin_nc,3,mean,na.rm=TRUE)

rlds_nc <- nc_open("E:/ISIMIP/bios/rlds/mean/year126_GDE_mean.nc")
rlds_nc <- ncvar_get(rlds_nc, "rlds")
rlds_mean <- apply(rlds_nc,3,mean,na.rm=TRUE)

evap_nc <- nc_open("E:/ISIMIP/bios/evap/mean/year126_GDE_mean.nc")
evap_nc <- ncvar_get(evap_nc, "evap")
evap_mean <- apply(evap_nc,3,mean,na.rm=TRUE)

lai_nc <- nc_open("E:/ISIMIP/bios/lai/mean/year126_GDE_mean.nc")
lai_nc <- ncvar_get(lai_nc, "lai")
lai_mean <- apply(lai_nc,3,mean,na.rm=TRUE)

pr_nc <- nc_open("E:/ISIMIP/bios/pr/mean/year126_GDE_mean.nc")
pr_nc <- ncvar_get(pr_nc, "pr")
pr_mean <- apply(pr_nc,3,mean,na.rm=TRUE)

ppl_nc <- nc_open("E:/ISIMIP/bios/population/ppl_GDE_ssp126.nc")
ppl_nc <- ncvar_get(ppl_nc, "ppl")
ppl_mean <- apply(ppl_nc,3,mean,na.rm=TRUE)

gdp_nc <- nc_open("E:/ISIMIP/bios/population/gdp_GDE_ssp126.nc")
gdp_nc <- ncvar_get(gdp_nc, "gdp")
gdp_mean <- apply(gdp_nc,3,mean,na.rm=TRUE)
#co2
co2 <- c(399.95, 403.12, 405.75, 408.59, 411.42, 414.23, 417.04, 419.81, 422.5, 425.12, 427.67, 430.17, 432.6, 434.97, 437.29, 439.56, 441.78, 443.93, 445.99, 447.97, 449.87, 451.68, 453.43, 455.09, 456.68, 458.2, 459.65, 461.02, 462.31, 463.54, 464.68, 465.75, 466.75, 467.68, 468.53, 469.31, 470.02, 470.66, 471.25, 471.78, 472.25, 472.66, 473.02, 473.32, 473.56, 473.75, 473.88, 473.96, 474.0, 474.0, 473.96, 473.87, 473.75, 473.58, 473.36, 473.11, 472.81, 472.46, 472.04, 471.56, 471.02, 470.41, 469.75, 469.02, 468.24, 467.39, 466.48, 465.54, 464.56, 463.56, 462.53, 461.47, 460.38, 459.26, 458.1, 456.92, 455.71, 454.5, 453.32, 452.16, 451.02, 449.91, 448.81, 447.73, 446.67, 445.62)
#
fbnf_nc <- nc_open("E:/ISIMIP/bios/fbnf/mean/year126_GDE_mean.nc")
fbnf_nc <- ncvar_get(fbnf_nc, "fbnf")
fbnf_mean <- apply(fbnf_nc,3,mean,na.rm=TRUE)

# 生成模拟数据（标准化数据）
data <- data.frame(
  gwr       = gwr_mean,
  qs        = qs_mean,
  snd       = snd_mean,
  soilmoist = soilmoist_mean,
  tmax      = tmax_mean,
  tmin      = tmin_mean,
  rlds      = rlds_mean,
  evap      = evap_mean,
  lai       = lai_mean,
  ppl       = ppl_mean,
  gdp       = gdp_mean,
  pr        = pr_mean,
  co2       = co2,
  fbnf      = fbnf_mean
)
# 路径矩阵（下三角矩阵）
path_matrix <- rbind(
  c(0, 0, 0, 0, 0),  # HA
  c(1, 0, 0, 0, 0),  # CT
  c(1, 1, 0, 0, 0),  # HY
  c(1, 1, 1, 0, 0),  # VG
  c(1, 1, 1, 1, 0)   # Fbnf 受 HY、CT、HA、VG 影响
)
colnames(path_matrix) <- rownames(path_matrix) <- c("HA", "CT", "HY", "VG", "Fbnf")
# 义测量模型（block 指定哪些观测变量对应哪个潜变量）
blocks <- list(
  c("ppl", "gdp", "co2"),                   # HA
  c("tmax", "tmin", "rlds", "evap", "pr"), # CT
  c("gwr", "qs", "soilmoist", "snd"),       # HY
  c("lai"),                           # VG
  c("fbnf")                           # Fbnf
)

modes <- c("A", "A", "A", "A", "A")
# 运行 PLS-SEM
pls_model <- plspm(data, path_matrix, blocks, modes = modes)
# 输出结果
print(summary(pls_model))   # 模型摘要

