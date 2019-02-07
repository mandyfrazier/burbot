Downloading Data
================

Tried to copy sequences to Seagate backup, but could not because it's
Read Only. Revisit this after break. The sequence data is on the Mac in
the lab and backed up via the Genome Center (they guarantee storage for
a few months or so). Will backup to harddrives in January. Checked that
all the data downloaded and there are the correct number of files in
each directory, so it appears that it downloaded correctly.

Re-Formating Hard Drive
-----------------------

The reason that we couldn't back up the data is because the hard drive
is formatting incorrectly for Mac OS (even though it works with some
files... mysteries of the universe). Solution: 1. Copy all files on
external hard drive to desktop to temporarily save them. 2. Open "Disk
Utility" on your Mac (search using the toolbar) 3. Select the highest
level of the external hard drive and click on the "Erase" button at the
top. 4. Rename the external hard drive to something useful, select "Mac
OS Extended (Journaled)" under "Format", and "GUID Partition Map" under
"Scheme". Select "Erase".

Sequence QC
===========

Installing MultiQC
------------------

Must have Anaconda and Python to download and use MultiQC.

1.  Download Anaconda for Python 2.7:
    <https://www.anaconda.com/download/#macos>
2.  Followed the installation wizard.  
3.  Closed Terminal so that it would take effect.
    <https://conda.io/docs/user-guide/install/macos.html>
4.  Downloaded MultiQC through the terminal following the instructions
    (<http://multiqc.info/docs/>) Tried to download using the "conda"
    method, but that didn't work. Instead, I used the "pip" method,
    which worked.

Installing FastQC
-----------------

Intsructions from the website, with my comments:


    chmod 755 fastqc
        
    #..but once you have done that you can run it directly:
        
    ./fastqc
        
    #..or place a link in /usr/local/bin to be able to run the program from any location:  **The following code did not work for me, so it is commented out.
        
    #sudo ln -s /path/to/FastQC/fastqc /usr/local/bin/fastqc  #From the website
        
    #sudo ln -s /Applications/FastQC/fastqc /usr/local/bin/fastqc   #For my computer 
        

Running FastQC works now as long as you're in the fastqc folder:

    ./fastqc /Users/amanda/Desktop/Practice_DE/burbot-Ma01-C_S203_L004_R1_001.fastq.gz /Users/amanda/Desktop/Practice_DE/burbot-Ma01-C_S203_L004_R2_001.fastq.gz -o /Users/amanda/Desktop/Practice_DE/

I couldn't figure out how to get fastqc to run universally, only when
you are in the fastqc directory. To work around this problem I used a
script that specifies the location of the fastqc program and runs it on
every sample in that directory. Since I have four directories (one for
each lane), I put the script in each directory and edited the output
file for each one.

Name of Script: fastq.sh


    #!/bin/bash

    for f in *.gz
        do
            /Users/amanda/Applications/FastQC/fastqc $f -o /Volumes/Burbot/Burbot_DE/8tq5hvxqm/Unaligned/Project_ATMF_L4_H1181P_Frazier/fastqc_output
        done
            

Running MultiQC
---------------

\*\*It seems like MultiQc isn't installed if you just type "multiqc ."
into the Terminal. This is because you must activate an Anaconda
environment to run Multiqc! To activate Anaconda:

    source activate py2.7 

Then, (py2.7) will appear at the beginning of the line in Terminal. You
can then go to the directory with the fastqc files and run multiqc:

    multiqc . 

(The period is telling multiqc to run in the current directory.)