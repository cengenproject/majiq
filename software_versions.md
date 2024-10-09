
Awk:
```
GNU Awk 4.2.1, API: 2.0 (GNU MPFR 3.1.6-p2, GNU MP 6.2.0)
```

gtf2gff3.pl version 0.1


MAJIQ version 2.4.dev102+g2cae1507. Full Conda environment:
```
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                  2_kmp_llvm    conda-forge
boltons                   23.0.0          py311h06a4308_0
brotlipy                  0.7.0           py311h5eee18b_1002
bzip2                     1.0.8                h7b6447c_0
c-ares                    1.19.0               h5eee18b_0
ca-certificates           2023.7.22            hbcca054_0    conda-forge
certifi                   2023.7.22          pyhd8ed1ab_0    conda-forge
cffi                      1.15.1          py311h5eee18b_3
charset-normalizer        2.0.4              pyhd3eb1b0_0
conda                     23.5.2          py311h38be061_0    conda-forge
conda-content-trust       0.1.3           py311h06a4308_0
conda-libmamba-solver     23.7.0             pyhd8ed1ab_0    conda-forge
conda-package-handling    2.1.0           py311h06a4308_0
conda-package-streaming   0.8.0           py311h06a4308_0
cryptography              39.0.1          py311h9ce1e76_2
fmt                       9.1.0                hdb19cb5_0
icu                       58.2                 he6710b0_3
idna                      3.4             py311h06a4308_0
jsonpatch                 1.32               pyhd3eb1b0_0
jsonpointer               2.1                pyhd3eb1b0_0
krb5                      1.20.1               h143b758_1
ld_impl_linux-64          2.38                 h1181459_1
libarchive                3.6.2                h6ac8c49_2
libcurl                   8.2.1                h251f7ec_0
libedit                   3.1.20221030         h5eee18b_0
libev                     4.33                 h7f8727e_1
libffi                    3.4.4                h6a678d5_0
libgcc-ng                 13.2.0               h807b86a_0    conda-forge
libmamba                  1.5.0                h658169a_0    conda-forge
libmambapy                1.5.0           py311h527f279_0    conda-forge
libnghttp2                1.52.0               h2d74bed_1
libsolv                   0.7.24               he621ea3_0
libssh2                   1.10.0               hdbd6064_2
libstdcxx-ng              13.2.0               h7e041cc_0    conda-forge
libuuid                   1.41.5               h5eee18b_0
libxml2                   2.10.3               hcbfbd50_0
llvm-openmp               12.0.1               h4bd325d_1    conda-forge
lz4-c                     1.9.4                h6a678d5_0
mamba                     1.5.0           py311h3072747_0    conda-forge
ncurses                   6.4                  h6a678d5_0
openssl                   3.1.2                hd590300_0    conda-forge
packaging                 23.0            py311h06a4308_0
pcre2                     10.42                hebb0a14_0
pip                       23.1.2          py311h06a4308_0
pluggy                    1.0.0           py311h06a4308_1
pybind11-abi              4                    hd3eb1b0_1
pycosat                   0.6.4           py311h5eee18b_0
pycparser                 2.21               pyhd3eb1b0_0
pyopenssl                 23.0.0          py311h06a4308_0
pysocks                   1.7.1           py311h06a4308_0
python                    3.11.4               h955ad1f_0
python_abi                3.11                    2_cp311    conda-forge
readline                  8.2                  h5eee18b_0
reproc                    14.2.4               h295c915_1
reproc-cpp                14.2.4               h295c915_1
requests                  2.29.0          py311h06a4308_0
ruamel.yaml               0.17.21         py311h5eee18b_0
setuptools                67.8.0          py311h06a4308_0
six                       1.16.0             pyhd3eb1b0_1
sqlite                    3.41.2               h5eee18b_0
tk                        8.6.12               h1ccaba5_0
toolz                     0.12.0          py311h06a4308_0
tqdm                      4.65.0          py311h92b7b1e_0
tzdata                    2023c                h04d1e81_0
urllib3                   1.26.16         py311h06a4308_0
wheel                     0.38.4          py311h06a4308_0
xz                        5.4.2                h5eee18b_0
yaml-cpp                  0.7.0                h295c915_1
zlib                      1.2.13               h5eee18b_0
zstandard                 0.19.0          py311h5eee18b_0
zstd                      1.5.5                hc292b87_0
```

MAJIQ ran with additional HTSlib version 1.12.



R analysis (dpsi.R script), sessionInfo:
```
R version 4.4.0 (2024-04-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 22631)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] wbData_0.9.2    lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4    
 [6] purrr_1.0.2     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1  
[11] tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] eulerr_7.0.2        rappdirs_0.3.3      utf8_1.2.4          generics_0.1.3     
 [5] cengenDataSC_0.1.0  polylabelr_0.2.0    lattice_0.22-6      stringi_1.8.3      
 [9] hms_1.1.3           magrittr_2.0.3      RColorBrewer_1.1-3  grid_4.4.0         
[13] timechange_0.3.0    Matrix_1.7-0        wormDatasets_0.2.0  ggrepel_0.9.5      
[17] writexl_1.5.0       mgcv_1.9-1          fansi_1.0.6         scales_1.3.0       
[21] RApiSerialize_0.1.2 cli_3.6.2           rlang_1.1.3         crayon_1.5.2       
[25] polyclip_1.10-6     splines_4.4.0       bit64_4.0.5         munsell_0.5.1      
[29] withr_3.0.0         tools_4.4.0         qs_0.26.1           parallel_4.4.0     
[33] tzdb_0.4.0          colorspace_2.1-0    vctrs_0.6.5         printMat_0.1.0     
[37] R6_2.5.1            matrixStats_1.3.0   lifecycle_1.0.4     stringfish_0.16.0  
[41] bit_4.0.5           vroom_1.6.5         pkgconfig_2.0.3     RcppParallel_5.1.7 
[45] pillar_1.9.0        gtable_0.3.5        glue_1.7.0          Rcpp_1.0.12        
[49] tidyselect_1.2.1    rstudioapi_0.16.0   farver_2.1.1        nlme_3.1-164       
[53] ggsci_3.0.3         labeling_0.4.3      pheatmap_1.0.12     compiler_4.4.0 
```
