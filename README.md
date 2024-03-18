
# Fun  
![image](https://zenodo.org/badge/641802007.svg)  
**"Fun"** is the algorithm to identify **replication fountains**. The current version is still preliminary and needs further refinement.  

![image](https://github.com/zzdzr/Fun/blob/master/image/Fun.png)

# Contents
- [Description](#description)
- [Requirements](#requirements)
- [Usage](#usage)
- [Output](#output)
- [Version](#version)
  
# Description
- A high-quality fountain possesses two features: a remarkable **signal-over-noise (SoN)** ratio between the fountain center and its flanking regions and continuous interactions extending along the vertical direction against the diagonal.
- Three major steps are included in the “Fun” pipeline, including **summit identification**, **elongation evaluation** and **end determination**, and further validation with quality control.
- The .hic or .mcool files can be obtained from https://doi.org/10.17632/4254mwwbhk.1.
  
# Requirements
  Please ensure that you have already installed the following packages:
  ```
  Python(3.8+) packages: numpy, pandas, scipy, click, logging, numba, cooler, bioframe...
    
  ```

# Usage
## Data processing for Repli-HiC
The sequencing library of Repli-HiC requires some preliminary data processing:
- Firstly, the first 9 bases for R1 of the sequencing library need to be removed.
  `cutadapt -u 9 -o trimmed_R1.fq input_R1.fq.gz`

- Then, the **trimLinker** tool should be used to eliminate potential linker sequences within the paired reads.
  `trimLinker -m 1 -k 2 -l 16 -o output -n name -A ACGCGATATCTTATC -B AGTCAGATAAGATAT trimmed_R1.fq input_R2.fq.gz`
  > trimLinker can be obtained from ChIA-PET2(https://github.com/GuipengLi/ChIA-PET2).

- Afterwards, the processed sequencing library can be integrated with 4DN Hi-C pipelines for further analysis.
  > The 4DN Hi-C data processing pipeline includes alignment, filtering, and matrix aggregation and normalization steps. https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline

## Fountain identification
Generally, there are three major steps in identification of replication fountains. The current version requires step-by-step execution (a more comprehensive version will be updated subsequently, please give me a moment to refine ^_^).
- **Calculate the signal-over-noise (SoN)**.
This calculation module is separated and can be used to calculate the SoN at a given resolution independently. In the present version, you have the ability to define the **search extent** for the sampling box, specify the **width of the sampling box**, and constrain the **offset**. For example, if you want to calculate SoN at 10kb resolution, with search extent at 500kb, padding width at 20kb and offset at 50kb, You can use the following example code.
   ```
   Fun calculate-son-score input.mcool::resolutions/10000 --out_dir /output_dir --coverage_ratio 0 --chromsize_path ChromInfo.txt 
   --ext_length 500000 --padding_width 2 --offset 50000 --integrate True --use_mean True
   ```

- **Identify potential summits of fountains**.
In current version, we attempt to find summits based on an algorithm from cooltools. This calculation module is based on the results of the previous SoN calculation. Therefore, before executing this module, please ensure that the SoN track can be correctly outputted.
   ```
  Fun generate-summits input.mcool::resolutions/10000 --track SoN_10000_merged.bedgraph --out_dir /output_dir
   ```

- **Identify fountains**.
Before identifying the fountains, please remove the summits that fall into low-quality genomic regions. In this module, you can perform algorithm within a Hi-C matrix of given normalization method and resolution. You can specify the **width of the sampling box, the length of the offset**, and also set the **step size for the sliding layer (--extension_pixels)**, **threshold for the p-value** (--p_value) and **fold change** (--signal_noise_background).
   ```
   Fun find-fountains input.mcool::resolutions/10000 --ext_length 500000 --half_width 2 --norm VC_SQRT --region_path Summits_10000_merged.bed
   --extension_pixels 10 100 5 --offset 50000 --interval_length 50000 --coverage_ratio 0 --p_value 0.05 --signal_noise_background 1.3 --output /output_dir/fountains_10kb
   ```
# Output
### Result Files:

  - ***_1.3.tab**: This file contains identified fountains after quality control. In current version, it contains 14 columns:


    
|chrom |start|end  |name |score |strand  |perc_res_list |max_extension |signal_noise_upstream|signal_noise_downstream|signal_noise_average_background| p_values |
|----|-----|-----|----|------|-----|------|------|------|------|----|---|
|chr1|4820000 |4830000|.|1.99 |.|nan, ... , 0.45, 0.45|150|4.60|4.07|4.32|1.1e-08|  

#### Column Explanation:
| Parameter | Description |
| --- | --- |
| **chrom, start, end** | Genomic coordinates of summits. |
| **score** | Signal-over-noise (SoN) score for summits. |
| **perc_res_list** | List containing the proportion of pixels where the central signal region is dominant over the background region at a given distance (The 'nan' values are due to the offset). |
| **max_extension** | The extension length of identified fountains. |
| **signal_noise_upstream** | The ratio of the signal in the central signal region to the signal in the upstream background sampling region. |
| **signal_noise_downstream** | The ratio of the signal in the central signal region to the signal in the downstream background sampling region. |
| **signal_noise_average_background** | The ratio of the signal in the central signal region to the average signal in the upstream and downstream background regions. |
| **p_values** | The p-value obtained from the Kolmogorov-Smirnov (K-S) test performed using the signal from the central region and the average signal from the background region. |

<br>

  - ***_1.3.bedpe**: This file is converted based on the results of the .tab file


| chr1 | x1 | x2 | chr2 | y1 | y2 |
|----|----|----|----|----|----|
|chr1|4675000|4685000|chr1|4975000|4985000|

<br>
  
### Visualization:
- The current version does not yet offer visualization method, but you can generate the following display through matplotlib.  


![image](https://github.com/zzdzr/Fun/blob/master/image/Fountains.png)
  
# Version
  - Fun v1.0.1
