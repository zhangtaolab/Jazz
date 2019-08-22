Jazz 

Non-Histone protein banding site identification

Jazz Dependncies:
samtools http://samtools.sourceforge.net
pysam http://code.google.com/p/pysam/
scipy http://www.scipy.org
numpy http://www.numpy.org

Author: Tao Zhang @ Yangzhou University

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit.
  -d DATAFILE, --data=DATAFILE
                        data file, should be sorted bam format
  -c CONTROLFILE, --control=CONTROLFILE
                        control(input) file, should be sorted bam format
  -n SAMPLENAME, --name=SAMPLENAME
                        NH sample name default=NH_sample
  -b BW, --bandwidth=BW
                        kernel smooth band width, should >20, default==600
  -t THRESHOLD, --threshold=THRESHOLD
                        Hot spots threshold, default=4.0
  -l MINLENGTH, --minlength=MINLENGTH
                        minimum length of hot spots, default=50
  -p PVALUE, --pavlue=PVALUE
                        p-value cutoff for peak identification, default=0.01
  -i INITIAL, --initial=INITIAL
                        Peak's initial length, >5 and <minlength, default=20
  --threads=NTHREADS    threads number or cpu number, default=4
  -w, --wig             whether out put wiggle file, default=False
  -f, --fdr             using FDR instead p-value
  -x EXCLUDECHR, --excludechr=EXCLUDECHR
                        Don't count those DHs, example='-x ChrM,ChrC'
  -g, --gff             whether out put gff file, default=False
  -j JOBTYPE, --jobtype=JOBTYPE
                        job type, such as dh, nhpaired or nhsingle
  -m MAXINSERT, --maxinsert=MAXINSERT
                        when you use paired library, please set the maxinsert
                        size
  --pe                  paired-end reads or single-end reads, default=False
                        (single end)