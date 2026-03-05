# VIPER - Viral Identification Pipeline for Emergency Response

<img src="images/VIPER_ico.png" width="420" style="display: block; margin: 0 auto;"/>

[Leia em português](README-pt.md)

VIPER is a user-friendly toolkit for assembling and identifying viral genomes (SARS-CoV-2, Dengue, and Influenza) from Illumina sequencing data. Created by the [CeVIVAS](https://cevivas.butantan.gov.br/) bioinformatics team at the Butantan Institute, VIPER combines an intuitive Windows GUI with robust command-line pipelines.

## Index

- [VIPER - Viral Identification Pipeline for Emergency Response](#viper---viral-identification-pipeline-for-emergency-response)
  - [Index](#index)
  - [Highlights](#highlights)
  - [Supported pipelines](#supported-pipelines)
  - [Installation](#installation)
    - [Windows (GUI installer)](#windows-gui-installer)
    - [Linux (CLI, no GUI)](#linux-cli-no-gui)
  - [Usage](#usage)
    - [Windows GUI](#windows-gui)
    - [Linux command line](#linux-command-line)
  - [Pipeline overviews](#pipeline-overviews)
  - [Future implementations](#future-implementations)
  - [Copyright and licence](#copyright-and-licence)

## Highlights

- Multi-virus pipelines with curated references and lineage assignment.
- Windows GUI for non-technical users; Linux command-line for servers and pipelines.
- Resource-aware execution: tested in low computational resources setups.
- Open source and customizable for lab or surveillance workflows.

## Supported pipelines

- SARS-CoV-2 assembly (Illumina)
- Dengue assembly (Illumina)
- Influenza assembly (Illumina)

## Installation

### Windows (GUI installer)

- Download the Windows installer (`.exe`) from the project [Releases](../../releases). The installer bundles the GUI and sets up the WSL-based backend.
- Requirements: Windows 10 2004+ (build 19041+) or Windows 11, admin permissions, and ~10 GB free disk space. An internet connection is required to download dependencies during setup.
- Antivirus tools can flag installers; if that happens, whitelist the VIPER installer.

**Step-by-step with screenshots**

1) Choose optional tasks (desktop shortcut).  
<img src="images/VIPER_Install_step_1.png" width="640"/>

2) Confirm installation.  
<img src="images/VIPER_Install_step_2.png" width="640"/>

3) Copying files.  
<img src="images/VIPER_Install_step_3.png" width="640"/>

4) Core/WSL bootstrap.  
<img src="images/VIPER_Install_setp_4.png" width="640"/>

5) Allow WSL enablement/restart if prompted.  
<img src="images/VIPER_Install_step_5.png" width="420"/>

6) Windows enabling WSL.  
<img src="images/VIPER_Install_step_6.png" width="900"/>

7) Bootstrap inside WSL (micromamba + environment).  
<img src="images/VIPER_Install_step_7.png" width="900"/>

8) Packages being installed.  
<img src="images/VIPER_Install_step_8.png" width="900"/>

9) Finalizing links and PATH inside WSL.  
<img src="images/VIPER_Install_step_9.png" width="900"/>

10) Finish the installer.  
<img src="images/VIPER_Install_step_10.png" width="640"/>

After installation, launch VIPER from the Start menu or run `VIPER.exe` from the installation folder. If you need to uninstall, use Windows Apps & Features and, if requested, run `wsl --unregister VIPER-Core` from an elevated PowerShell.

### Linux (CLI, no GUI)

On Linux, VIPER runs via command line only.

**Prerequisites**

- Bash, git, curl, and `sudo` (for optional symlinks).
- Micromamba or Conda (micromamba recommended).
- Python 3.8 (handled by the environment YAML).

**1) Get the code**

```bash
git clone https://github.com/alex-ranieri/viper.git
cd viper
```

**2) Install micromamba (if you do not already have micromamba/conda)**  
Example for bash:

```bash
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
./bin/micromamba shell init -s bash -p ~/.local/share/micromamba
source ~/.bashrc
```

**3) Create the VIPER environment**

The YAML pins all dependencies. It sets the prefix to `~/.local/share/mamba/envs/VIPERGenomeAssembler`.

```bash
micromamba create -y -f Pipelines/VIPERGenomeAssembler.yaml
micromamba activate VIPERGenomeAssembler
# If activation by name fails, activate by path:
# micromamba activate ~/.local/share/mamba/envs/VIPERGenomeAssembler
```

**4) Make the pipeline scripts available**

Option A (recommended, creates `/usr/local/bin/VIPER_*.sh` wrappers):

```bash
sudo bash update.d/post-update.sh
```

Option B (no sudo): add the modules folder to PATH and expose `$PIPELINE`:

```bash
cd Pipelines/Modules
chmod +x configurePATH.sh
./configurePATH.sh
# Reload your shell or source ~/.bashrc to pick up PATH/PIPELINE.
```

## Usage

### Windows GUI

Launch `VIPER.exe`. The GUI has a basic and an advanced view:

<img src="images/VIPER_GUI_basic.png" width="640"/>  
<img src="images/VIPER_GUI_Advanced.png" width="640"/>

In the Advanced tab, you can change the maximum RAM VIPER may use: enter the value in the required format, click Apply, then use Update to confirm it applied.

VIPER refreshes viral databases before running; you can also trigger database updates manually when needed.

Choose how many threads VIPER will use for assembly before starting a run.

Run pipeline steps:

⚠️ Ensure your FASTQ files follow the Illumina naming convention used in BaseSpace: [Illumina FASTQ naming guide](https://support.illumina.com/help/BaseSpace_Sequence_Hub_OLH_009008_2/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm)

1) Select the folder containing your `*.fastq.gz` pairs.  
<img src="images/VIPER_GUI_Run_pipeline_step_1.png" width="720"/>  

1.1) Choose the viral pipeline.  
<img src="images/VIPER_GUI_Run_pipeline_step_1.1.png" width="720"/>

2) After selecting folder and pipeline, review the summary and continue.  
<img src="images/VIPER_GUI_Run_pipeline_step_2.png" width="720"/>

3) Pipeline executed successfully; confirmation message displayed.  
<img src="images/VIPER_GUI_Run_pipeline_step_3.png" width="720"/>

4) Choose whether to open a simplified dashboard to view results.  
<img src="images/VIPER_GUI_Run_pipeline_step_4.png" width="720"/>

5) Dashboard view; export results to Excel if needed.  
<img src="images/VIPER_GUI_Run_pipeline_step_5.png" width="720"/>

6) Generated outputs in the pipeline folder (.fasta genomes plus per-sample subfolders with additional files).  
<img src="images/VIPER_GUI_Run_pipeline_step_6.png" width="720"/>

Manual database update (optional):  
<img src="images/VIPER_GUI_Update_database.png" width="640"/>

### Linux command line

Activate the environment first:

```bash
micromamba activate VIPERGenomeAssembler
```

If you ran `update.d/post-update.sh`, use the wrappers:

```bash
# {threads} = threads per sample, {samples} = parallel samples
VIPER_CoV.sh {threads} {samples}       # SARS-CoV-2
VIPER_DENV.sh {threads} {samples}      # Dengue
VIPER_FLU.sh {threads} {samples}       # Influenza
```

If you used `configurePATH.sh` without sudo, run directly via `$PIPELINE`:

```bash
bash "$PIPELINE/SARS-CoV-2/Exec_assembly_pipeline_Illumina_v8_bowtie2_ref_iVar_CeVIVAS.sh" {threads} {samples}
bash "$PIPELINE/DENV/Exec_assembly_pipeline_Illumina_v5_bwa_mem_ref_iVar.sh" {threads} {samples}
bash "$PIPELINE/Influenza/Exec_assembly_pipeline_Illumina_v3_Vapor_SPAdes.sh" {threads} {samples}
```

Results are written next to the input FASTQ directory. Test datasets are available under `Test_samples/`.

## Pipeline overviews

**SARS-CoV-2 assembly overview**  
<img src="images/SARS-CoV-2 iSNVs.png" width="700" style="display: block; margin: 20px auto;"/>

**Dengue assembly overview**  
<img src="images/DENV_ASSEMBLY.png" width="700" style="display: block; margin: 20px auto;"/>

**Influenza assembly overview**  
<img src="images/Influenza assembly pipeline.png" width="700" style="display: block; margin: 20px auto;"/>

## Future implementations

- Pipelines for additional viruses (e.g., Chikungunya, RSV).
- Automatic detection of the virus type being analyzed.
- Phylogeny module with interactive tree visualization via Auspice.
- Reworked assembly modules for iSNV detection.
- Allow users to supply custom reference sequences for assembly.
- Pipeline optimization for workflow managers (e.g., Snakemake, Nextflow).

## Copyright and licence

VIPER was created by the [CeVIVAS](https://cevivas.butantan.gov.br/) bioinformatics team at the Butantan Institute:

- Alex Ranieri J. Lima
- Gabriela Ribeiro
- Vinicius Carius De Souza
- Isabela Carvalho Brcko
- Igor Santana Ribeiro
- James Siqueira Pereira
- Vincent Louis Viala

Supervised by:

- Maria Carolina Quartim Barbosa Elias Sabbaga
- Sandra Coccuzzo Sampaio

VIPER is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

VIPER is distributed in the hope that it will be useful, but **without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose**. See the GNU General Public License for more details: http://www.gnu.org/licenses/.
