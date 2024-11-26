# Calibrated Differential RNA Editing Scanner (CADRES)
An Analytical Pipeline for Precise Identification of Differential RNA Editing Events Across Varied Biological Conditions
# Cardes Project

## Overview

This repository contains scripts for running two main components of the Cardes project: `cardes_boost.sh` and `cardes_DVR.sh`. These scripts are designed to perform specific tasks in the Cardes workflow.

## Getting Started

### Prerequisites

Before you begin, ensure you have the following installed on your system:

- Bash (usually pre-installed on Unix-based systems)
- Any additional dependencies required by the scripts (e.g., Python, R, etc.)

### Installation

Clone this repository to your local machine using Git:

```
git clone https://github.com/yourusername/cardes.git

cd cardes
```
Running the Scripts
Step 1: 

Run cardes_boost.sh
The first step is to run the cardes_boost.sh script. This script performs the initial setup and processing for the Cardes project.
```
./cardes_boost.sh
```
Step 2: Run cardes_DVR.sh
After the cardes_boost.sh script has completed successfully, run the cardes_DVR.sh script. This script continues the processing based on the output of the first script.
```
./cardes_DVR.sh
```
Dependencies
If there are any dependencies required by these scripts, list them here with their installation instructions.

Additional Information
Ensure you have the necessary permissions to execute the scripts.
The scripts assume the working directory is set to the root of the cloned repository.
