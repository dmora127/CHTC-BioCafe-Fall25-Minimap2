# Mapping Reads on CHTC with Minimap2

This tutorial will walk you through a complete long-read sequencing read mapping exercise using [Minimap2](). You'll learn how to:

* Map your reads to a reference genome using Minimap2
* Breakdown massive bioinformatics workflows into many independent smaller tasks
* Submit hundreds to thousands of jobs with a few simple commands
* Use the Open Science Data Federation (OSDF) to manage file transfer during job submission

All of these steps are distributed across hundreds (or thousands!) of jobs using the HTCondor workload manager and Apptainer containers to run your software reliably and reproducibly at scale. The tutorial is built around realistic genomics use cases and emphasizes performance, reproducibility, and portability. You'll work with real data and see how high-throughput computing (HTC) can accelerate your genomics workflows.

>[!NOTE]
>If you're brand new to running jobs on CHTC, we recommend completing the HTCondor ["Hello World"](https://chtc.cs.wisc.edu/uw-research-computing/htcondor-job-submission) exercise before diving into this tutorial.

**Let’s get started!**

Jump to...
<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->
- [Tutorial Setup](#tutorial-setup)
   * [Assumptions](#assumptions)
   * [Materials](#materials)
- [Setting up our software environment - Optional](#setting-up-our-software-environment---optional)
- [Mapping Whole Genome Sequencing Reads to a Reference Genome](#mapping-whole-genome-sequencing-reads-to-a-reference-genome)
   * [Indexing the Reference Genome - Generating `G_californianus_813.fasta.mmi`](#indexing-the-reference-genome---generating-g_californianus_813fastammi)
   * [Splitting Our Sequencing Reads](#splitting-our-sequencing-reads)
   * [Running Minimap to Map Reads to the Reference Genome](#running-minimap-to-map-reads-to-the-reference-genome)
- [Next Steps](#next-steps)

<!-- TOC end -->

## Tutorial Setup

### Assumptions

This tutorial assumes that you:

* Have basic command-line experience (e.g., navigating directories, using bash, editing text files).
* Have a working CHTC account and can log into an Access Point (e.g., ap2002.chtc.wisc.edu).
* Are familiar with HTCondor job submission, including writing simple .sub files and tracking job status with condor_q.
* Understand the general workflow of long-read sequencing analysis: basecalling → mapping → variant calling.

>[!TIP]
>You do not need to be a genomics expert to follow this tutorial. The commands and scripts are designed to be beginner-friendly and self-contained, while still reflecting real-world research workflows.

### Materials

To obtain a copy of the files used in this tutorial, you can

* Clone the repository, with 
  
  ```
  git clone https://github.com/dmora127/CHTC-BioCafe-Fall25-Minimap2.git
  ```

  or the equivalent for your device

* Set up the `/home` and `/staging` directories:

  ```
  mkdir -p ~/genomics_tutorial/
  cp -r ~/CHTC-BioCafe-Fall25-Minimap2/genomics_tutorial/home/* ~/genomics_tutorial/
  
  mkdir -p /staging/<NetID>/genomics_tutorial/
  cp -r ~/CHTC-BioCafe-Fall25-Minimap2/genomics_tutorial/staging/* /staging/<NetID>/genomics_tutorial/
  ```

* Download the toy dataset and container from the `/staging/groups/chtc_staff/` directory: 
  
    ```
    cp /staging/groups/chtc_staff/sd_0001_sub.fastq.gz ~/genomics_tutorial/inputs/
    cp /staging/groups/chtc_staff/G_californianus_813.fasta ~/genomics_tutorial/inputs/
    cp /staging/groups/chtc_staff/minimap2_08OCT2025_v1.sif /staging/<NetID>/genomics_tutorial/software/
    ```
    **or, if you want to copy all three files in parallel (faster! but less readable):**
    ```
    (cp /staging/groups/chtc_staff/sd_0001_sub.fastq.gz ~/genomics_tutorial/inputs/ & cp /staging/groups/chtc_staff/G_californianus_813.fasta ~/genomics_tutorial/inputs/ & cp /staging/groups/chtc_staff/minimap2_08OCT2025_v1.sif /staging/$USER/genomics_tutorial/software/ & wait)
    ```

<!--TODO: Add to pelican readable origin-->

## Setting up our software environment - Optional
Before we can begin mapping our reads, we need to setup our software environment to run Minimap2. We are going to setup our environment using an Apptainer container. 

1. First, let's login to our CHTC Account

    ```
    ssh <NetID>@<ap2002/ap2001>.chtc.wisc.edu
    ```

2. We now need to write up a definition file for singularity to build our Minimap2 container. Copy and paste this block of
text to a new file titled `minimap2.def`. You can open up a text editor, such as `vim` or `nano` using a command like: `vim minimap2.def`.

    ```
    Bootstrap: docker
    From: continuumio/miniconda3:latest
    
    %post
        conda install -c bioconda -c conda-forge minimap2 samtools bedtools -y
    ```

    This definition file uses the latest Anaconda (formerly Continuumio) `Miniconda3` base image from Docker and conda installs `minimap2`, `samtools`, and `bedtools` from the `bioconda` and `conda-forge` channels.

>[!TIP]
> A pre-built minimap2 container is available in the `/staging/groups/chtc_staff/` directory. You can copy it to your `/staging/<NetID>/genomics_tutorial/software/` directory using the command: `cp /staging/groups/chtc_staff/minimap2_08OCT2025_v1.sif /staging/<NetID>/genomics_tutorial/software/`. If copying the pre-built container, you can skip steps 3-6 below.

3. Next we need to write our interactive apptainer build job submission script. This submit file allows you to _build_ a container image interactively on a dedicated build node. Copy and paste this block of text to a new file titled `minimap2_build.sub`. You can open up a text editor, such as `vim` or `nano` using a command like: `vim minimap2_build.sub`. **Make sure to replace `<image.def>` with the path to your `minimap2.def` file you created in step 2.**

    ```
    # build.sub
    # For building an Apptainer container
    
    universe = vanilla
    log = build.log
    
    # If you have additional files in your /home directory that are required for your container, add them to the transfer_input_files line as a comma-separated list.
    #transfer_input_files = <image.def>
    
    requirements = (HasCHTCStaging == true)
    
    +IsBuildJob = true
    request_cpus = 4
    request_memory = 16GB
    request_disk = 16GB
    
    queue
    ```
>[!CAUTION]
> You should **not** submit a standard job to build your container. Building containers requires special permissions and resources that are only available on dedicated build nodes. If you submit a standard job, it will likely fail due to insufficient permissions or resources.
> Similarly, you should **never** build containers on an Access Point (AP). APs are a shared resource and are not configured to support container builds. Building containers on an AP can lead to performance issues for other users and may violate the terms of service of CHTC. Repeated attempts to build containers on an AP may result in your account being disabled for misuse.

4. Submit your interactive apptainer build job to CHTC by running the following command:
    ```
    condor_submit -i minimap2_build.sub
   ```
   
    The `-i` flag tells HTCondor to run the job interactively, allowing you to see the build process in real-time. This is useful for debugging and ensuring that the container builds correctly.
   
5. Once your prompt changes to the build node, build your apptainer container on the Access Point (AP) by running the following command:
    ```
    apptainer build minimap2_<08OCT2025>_<v1>.sif minimap2.def
   ```
   
6. Move your finalized container image `minimap2_08OCT2025_v1.sif` to your `/staging` directory
    
    ```
   mv minimap2_08OCT2025_v1.sif /staging/<NetID>/
   ```
    
> [!TIP]
> You should always use unique filenames for all your files in the `/staging/` directory. The `/staging/` directory is used by the OSDF. **Files transferred with the OSDF are aggressively cached**, so if you use the same filename for different files, you may inadvertently use a cached version of a file you did not intend to use. **This can cause your jobs to go on hold repeatedly or deliver incorrect/inaccurate results.** A good convention is to include the date and a version number in your filenames. For example, `minimap2_<08OCT2025>_<v1>.sif`. When you update your container, increment the version number (e.g., `v2`, `v3`, etc.) and the date to always keep a new name record for the file.
    
## Mapping Whole Genome Sequencing Reads to a Reference Genome

### Indexing the Reference Genome - Generating `G_californianus_813.fasta.mmi`

Before we can map our reads to the reference genome, we need to index the reference genome using `minimap2`. This indexing step only needs to be done once per reference genome and creates a file with the `.mmi` extension that `minimap2` uses for mapping.

1.  Create `executables/minimap2_index.sh` using either `vim` or `nano`
    ```
    #!/bin/bash
    minimap2 -x map-ont -d "${1}.mmi" "$1"
    ```
    
> [!NOTE]
> This script takes one argument from the submit file (see step 2): the name of the reference genome FASTA file. It uses `minimap2` to create an index file with the `.mmi` extension.
   
2. Create `minimap2_index.sub` using either `vim` or `nano`. Make sure to replace `<input_ref_genome_fasta_file_name>` with the name of your reference genome FASTA file (e.g., `G_californianus_813.fasta`). Replace `<NetID>` with **your** actual NetID.
    ```
    container_image        = osdf:///chtc/staging/<NetID>/genomics_tutorial/software/minimap2_08OCT2025_v1.sif
    
    executable		       = ./executables/minimap2_index.sh
    arguments              = <input_ref_genome_fasta_file_name>
    
    transfer_input_files   = ./inputs/<input_ref_genome_fasta_file_name>
    
    transfer_output_files  = <input_ref_genome_fasta_file_name>.mmi
    transfer_output_remaps = "<input_ref_genome_fasta_file_name>.mmi = /staging/<NetID>/genomics_tutorial/inputs/<input_ref_genome_fasta_file_name>.mmi"
    
    output                 = ./logs/$(Cluster)_$(Process)_minimap2_indexing_step1.out
    error                  = ./logs/$(Cluster)_$(Process)_minimap2_indexing_step1.err
    log                    = ./logs/$(Cluster)_minimap2_indexing_step1.log
    
    request_cpus           = 4
    request_disk           = 10 GB
    request_memory         = 24 GB
    
    queue 1
    ```
    
> [!IMPORTANT]
> Replace `<input_ref_genome_fasta_file_name>` with the actual name of your reference genome FASTA file (e.g., `G_californianus_813.fasta`), not the path to the file. The path to the file is specified in the `transfer_input_files` attribute. Our executable script will look for the FASTA file in the top-level working directory of the job, which is where HTCondor places the transferred input files.

3. Submit your `minimap2_index.sub` job to the CHTC Pool.
    ```
    condor_submit minimap2_index.sub
    ```
   
> [!WARNING]  
> Index will take a few minutes to complete, **do not proceed until your indexing job is completed.** You can check the status of your job using `condor_watch_q` or by checking the log files specified in your submit file. Once the job is complete, you should see a new file with the `.mmi` extension in your `/staging/<NetID>/genomics_tutorial/inputs/` directory.

### Splitting Our Sequencing Reads

To get ready for our mapping step, we need to prepare our freshly sequenced reads. We will split our reads into smaller chunks to take advantage of CHTC's high-throughput computing capabilities. This will allow us to map our reads in parallel, significantly speeding up the overall process.

1. Split your reads into smaller chunks using `split`. You can adjust the `-p` parameter to change the number of chunks. Here, we are splitting our reads into 10 chunks.

    ```
    cd ~/genomics_tutorial/inputs/
    split -l 4000 <condor_wgs_reads>.fastq subset_
    ```
    
    This command splits the `cali_condor_wgs_reads.fastq` file into smaller files, each containing 4000 lines (which corresponds to 1000 reads, since each read in a FASTQ file is represented by 4 lines). It prepends the prefix `subset_` to each split output file, this will help us in the next step when listing our input files for HTCondor. The output files will be named `subset_aa_cali_condor_wgs_reads.fastq`, `subset_ab_cali_condor_wgs_reads.fastq`, `subset_ac_cali_condor_wgs_reads.fastq`, etc. 

> [!IMPORTANT]
> Since our read subsets are relatively small and not repeatedly used, we are **not** saving them to the `/staging/` directory. Instead, we are placing them directly in our `/home/<NetID>/genomics_tutorial/inputs/` directory. This is because the `/staging` directory is optimized for larger files that benefit from caching and repeated access. Small files like these read subsets do not benefit from the caching mechanisms and can be transferred more efficiently by placing them directly in our home directory. **Repeated writing to `/staging` can also lead to significant performance issues across all users, so we should avoid using it for small, temporary files.**

2. Generate a list of the split fastq files. Save it as `listofReads.txt` in your project base directory in your `/home/<NetID>` path. 

    ```
   cd ~/genomics_tutorial/inputs/
   ls subset_* > ~/genomics_tutorial/listofReads.txt
   ```
   
> [!TIP]
> This generates a list of all the split FASTQ files and saves it to `listofReads.txt` which will represent our list of jobs for HTCondor to execute. Each line in this file represents a separate job that HTCondor will execute. This file will be used in our HTCondor submit file to iterate over each split FASTQ file for mapping. The concept of a **_list of jobs_** is central to how HTCondor and HTC works. 
   
### Running Minimap to Map Reads to the Reference Genome
 
Now that we have our reference genome indexed and our reads split into smaller chunks, we can proceed to map our reads to the reference genome using `minimap2`. We will submit a cluster of jobs to the CHTC Pool, where each job will map one of the split FASTQ files to the reference genome in parallel.

1.  Create `executables/minimap2_mapping.sh` using either `vim` or `nano`

    ```
    #!/bin/bash
    # Use minimap2 to map the reads to the reference genome
    ./minimap2 -ax map-ont "$1" "$2" > "mapped_${2}_reads_to_genome.sam"
    
    # Use samtools to sort our mapped reads BAM, required for downstream analysis
    samtools sort "mapped_${2}_reads_to_genome.sam" -o "mapped_${2}_reads_to_genome_sam_sorted.bam"
    ```
    
>  [!NOTE]  
> Notice that this script takes two arguments from the submit file (see step 2): the name of the reference genome index file ($1) and the name of the split FASTQ file ($2). It uses `minimap2` to map the reads in the FASTQ file to the reference genome and outputs a SAM file. It then uses `samtools` to sort the SAM file and outputs a sorted BAM file, which is required for downstream analysis.

2. Create `minimap2_mapping.sub` using either `vim` or `nano`. Replace `<ref_genome_mmi_file_name>` with the name of your reference genome index file (e.g., `G_californianus_813.fasta.mmi`) in the `arguments` and `transfer_input_files` lines. Replace `<NetID>` with **your** actual NetID.
    ```
    container_image        = osdf:///chtc/staging/<NetID>/genomics_tutorial/software/minimap2_08OCT2025_v1.sif
    
    executable		       = ./executables/minimap2_mapping.sh
    arguments              = <ref_genome_mmi_file_name> $(READ_SUBSET)
   
    transfer_input_files   = osdf:///chtc/staging/<NetID>/genomics_tutorial/inputs/<ref_genome_mmi_file_name>, inputs/$(READ_SUBSET)
    
    transfer_output_files  = mapped_$(READ_SUBSET)_reads_to_genome_sam_sorted.bam
    transfer_output_remaps = "mapped_$(READ_SUBSET)_reads_to_genome_sam_sorted.bam = outputs/mapped_$(READ_SUBSET)_reads_to_genome_sam_sorted.bam"
    
    output                 = ./logs/$(Cluster)_$(Process)_mapping_$(READ_SUBSET)_step2.out
    error                  = ./logs/$(Cluster)_$(Process)_mapping_$(READ_SUBSET)_step2.err
    log                    = ./logs/$(Cluster)_mapping_step2.log
    
    request_cpus           = 2
    request_disk           = 5 GB
    request_memory         = 10 GB
    
    queue READ_SUBSET from listofReads.txt
    ```
    
     In this step, we **are not** transferring our outputs to the `/staging` directory. The mapped/sorted BAM files are intermediate temporary files in our analysis and do not benefit from the aggressive caching of the OSDF. By default, HTCondor will transfer outputs to the directory where we submitted our job from. Since we want to transfer the sorted mapped BAMs to a specific directory, we can use the `transfer_output_remaps` attribute on our submission script. The syntax of this attribute is:
   
     ```transfer_output_remaps = "<file_on_execution_point>=<desired_path_to_file_on_access_point>``` 
    
3. Submit your cluster of minimap2 jobs to the CHTC Pool.
   
   ```
   condor_submit minimap2_mapping.sub
   ```

## Next Steps

Now that you've completed the read mapping tutorial on CHTC, you're ready to adapt these workflows for your own data and research questions. Here are some suggestions for what you can do next:

🧬 Apply the Workflow to Your Own Data
* Replace the tutorial datasets with your own files and reference genome.
* Modify the mapping submit files to fit your data size, read type (e.g., ONT vs. PacBio), and resource needs.

🧰 Customize or Extend the Workflow
* Incorporate quality control steps (e.g., filtering or read statistics) using FastQC.
* Use other mappers or variant callers, such as ngmlr, pbsv, or cuteSV.
* Add downstream tools for annotation, comparison, or visualization (e.g., IGV, bedtools, SURVIVOR).

📦 Create Your Own Containers
* Extend the Apptainer containers used here with additional tools, reference data, or dependencies.
* For help with this, see our [Containers Guide](https://chtc.cs.wisc.edu/uw-research-computing/software-overview-htc#main).

🚀 Run Larger Analyses
* Submit thousands of mapping jobs across the CHTC pool.
* Explore data staging best practices using the OSDF for large-scale genomics workflows.
* Consider using workflow managers (e.g., [DAGman](https://chtc.cs.wisc.edu/uw-research-computing/htc/dagman-workflows#main)) with HTCondor.

🧑‍💻 Get Help or Collaborate
* Reach out to [chtc@cs.wisc.edu](mailto:chtc@cs.wisc.edu) for one-on-one help with scaling your research.
* Attend office hours or training sessions—see the [CHTC Get Help Page](https://chtc.cs.wisc.edu/uw-research-computing/get-help.html) for details.

### Software

In this tutorial, we created several *starter* apptainer containers, including tools like: Dorado, SAMtools, Minimap, and Sniffles2. These containers can serve as a *jumping-off* for you if you need to install additional software for your workflows. 

Our recommendation for most users is to use "Apptainer" containers for deploying their software.
For instructions on how to build an Apptainer container, see our guide [Using Apptainer/Singularity Containers](https://chtc.cs.wisc.edu/uw-research-computing/apptainer-htc#main).
If you are familiar with Docker, or want to learn how to use Docker, see our guide [Using Docker Containers](https://chtc.cs.wisc.edu/uw-research-computing/docker-jobs#main).

This information can also be found in our guide [Using Software on CHTC](https://chtc.cs.wisc.edu/uw-research-computing/software-overview-htc#main).

### Data

The ecosystem for moving data to, from, and within the HTC system can be complex, especially if trying to work with large data (> gigabytes).
For guides on how data movement works on the HTC system, see our [Data Staging and Transfer to Jobs](https://chtc.cs.wisc.edu/uw-research-computing/htc-job-file-transfer#main) guides.

### GPUs

CHTC has GPU nodes available for common use, like the ones used in this tutorial. If you would like to learn more about our GPU capacity, please visit our [GPU Guide on the CHTC Documentation Portal](https://chtc.cs.wisc.edu/uw-research-computing/gpu-jobs#main).

## Getting Help

The CHTC Research Computing Facilitators are here to help researchers using CHTC for their research. We provide a broad swath of research facilitation services, including:

* **Web guides**: [CHTC - HTC Cluster Guides](https://chtc.cs.wisc.edu/uw-research-computing/htc/guides.html) - instructions and how-tos for using the CHTC and OSDF.
* **Email support**: get help within 1-2 business days by emailing [chtc@cs.wisc.edu](mailto:chtc@cs.wisc.edu).
* **Virtual office hours**: live discussions with facilitators - see the [Email, Office Hours, and 1-1 Meetings](https://chtc.cs.wisc.edu/uw-research-computing/get-help.html) page for current schedule.
* **One-on-one meetings**: dedicated meetings to help new users, groups get started on the system; email [chtc@cs.wisc.edu](mailto:chtc@cs.wisc.edu) to request a meeting.

This information, and more, is provided in our [Get Help](https://chtc.cs.wisc.edu/uw-research-computing/get-help.html) page.