# Source code for ''Wasserstein distributionally robust IMM algorithm"

## Files Description
* img: Output images (.fig file)
* src: All the functions code (.m file)
* main_unpurt.m: The script to make a simulation without disturbance
* cvct_u.m: The script to make a simulation under good situation
* cvct_o.m: The script to make a simulation under medium situation
* cvct_n.m: The script to make a simulation under bad situation
* cvct_u.mat: The saved output workspace produced by cvct_u.m
* cvct_o.mat: The saved output workspace produced by cvct_o.m
* cvct_n.mat: The saved output workspace produced by cvct_n.m
* huatu.m: The script to draw all the pictures

## Download
You can download it directly or just clone it by using Git
```
git clone https://github.com/KudoSZhang/WIMM.git
```

## Usage
Firstly you should make sure that src folder is added to your working path in matlab. Then run cvct_u.m, cvct_o.m and cvct_n.m. Then save the output workspace as cvct_u.mat, cvct_o.mat and cvct_n.mat. Finally, run huatu.m to draw all the pictures.
