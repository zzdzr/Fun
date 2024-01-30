
# FUN
**"FUN"** is the algorithm to identify **replication fountains**. The current version is still preliminary and needs further refinement.

# Contents
- [Dependencies](#dependencies)
- [Description](#description)
- [Usage](#usage)

# Description
- A high-quality fountain possesses two features: a remarkable **signal-over-noise (SoN)** ratio between the fountain center and its flanking regions and continuous interactions extending along the vertical direction against the diagonal.
- Three major steps are included in the “Fun” pipeline, including **summit identification**, **elongation evaluation** and **end determination**, and further validation with quality control.


# Usage
## Data processing for Repli-HiC
The sequencing library of Repli-HiC requires some preliminary data processing:
- Firstly, the first 9 bases for R1 of the sequencing library need to be removed.
  `cutadapt -u 9 -o trimmed_R1.fq input_R1.fq.gz`

- Then, the **trimLinker** tool should be used to eliminate potential linker sequences within the paired reads.
  `trimLinker -m 1 -k 2 -l 16 -o output -n name -A ACGCGATATCTTATC -B AGTCAGATAAGATAT trimmed_R1.fq input_R2.fq.gz`
  > trimLinker can be obtained from ChIA-PET2.

- Afterwards, the processed sequencing library can be integrated with 4DN Hi-C pipelines for further analysis.
  > The 4DN Hi-C data processing pipeline includes alignment, filtering, and matrix aggregation and normalization steps. https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline

## Fountain identification
Generally, there are three major steps in identification of replication fountains. The current version requires step-by-step execution (a more comprehensive version will be updated subsequently, please give me a moment to refine ^_^).
- **Calculate the signal-over-noise (SoN)**.
This calculation module is separated and can be used to calculate the SoN at a given resolution independently. In the present version, you have the ability to define the **search extent** for the sampling box, specify the **width of the sampling box**, and constrain the **offset magnitude**. For example, if you want to calculate SoN at 10kb resolution, with search extent at 500kb, padding width at 20kb and offset at 50kb, You can use the following example code.
   ```
   Replihic calculate-son-score input.mcool::resolutions/10000 --out_dir /output_dir/ --coverage_ratio 0
   --ext_length 500000 --padding_width 2 --offset 50000 --merge True --use_mean True
   ```

- **Identify potential summits of fountains**.
In current version, we attempt to find summits based on an algorithm from cooltools. This calculation module is based on the results of the previous SoN calculation. Therefore, before executing this module, please ensure that the SoN track can be correctly outputted.
   ```
  Replihic generate-summits input.mcool::resolutions/10000 --track SoN_10000_merged.bedgraph --out_dir /output_dir/
   ```

- **Identify fountains**.
Before finally identifying the fountains, please remove the summits that fall into low-quality genomic regions. In this module, you can set the search step for the sampling box based on the summits that have been identified, within a Hi-C matrix of given normalization method and resolution. You can specify the **width of the sampling box, the length of the offset**, and also set the **step size for the layer (--extension_pixels)**, **threshold for the p-value** (--p_value) and the **fold change** (--signal_noise_background).
   ```
   Replihic find-fountains input.mcool::resolutions/10000 --ext_length 500000 --half_width 2 --norm VC_SQRT --region_path Summits_10000_merged.bed
   --extension_pixels 10 100 5 --offset 50000 --interval_length 50000 --coverage_ratio 0 --p_value 0.05 --signal_noise_background 1.1 1.2 1.3 1.4 1.5 --output /output_dir/result_10kb
   ```
# Version
  FUN v1.0
