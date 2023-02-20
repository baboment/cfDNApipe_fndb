# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 09:28:19 2019

@author: zhang
"""

from .StepBase import StepBase
from .cfDNA_utils import commonError
import os
from .Configure import Configure


__metaclass__ = type


class bowtie2(StepBase):
    def __init__(
        self,
        seqInput1=None,
        seqInput2=None,
        ref=None,
        outputdir=None,
        threads=1,
        paired=True,
        other_params={"-q": True, "-N": 1, "--time": True},
        stepNum=None,
        upstream=None,
        verbose=True,
        finaleDB=False,
        **kwargs
    ):
        """
        This function is used for mapping WGS data.
        Note: this function is calling bowtie2.

        bowtie2(seqInput1=None, seqInput2=None, ref=None, outputdir=None, threads=1, paired=True,
                other_params={"-q": True, "-N": 1, "--time": True}, stepNum=None, upstream=None,)
        {P}arameters:
            seqInput1: list, input _1 fastq files.
            seqInput2: list, input _2 fastq files, None for single end.
            ref: bowtie2 reference path.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, how many thread to use.
            paired: True for paired data, False for single end data.
            other_params: dict, other parameters passing to Bismark.
                          "-parameter": True means "-parameter" in command line.
                          "-parameter": 1 means "-parameter 1" in command line.
            stepNum: int or str, step flag for folder name.
            upstream: upstream output results, used for pipeline. This parameter can be True, which means a new pipeline start.
            verbose: bool, True means print all stdout, but will be slow; False means black stdout verbose, much faster.
        """


        super(bowtie2, self).__init__(stepNum, upstream)

        # set sequencing input
        if (upstream is None) or (upstream is True):
            self.setInput("seq1", seqInput1)
            self.setInput("seq2", seqInput2)
        else:
            Configure.configureCheck()
            upstream.checkFilePath()
            if upstream.__class__.__name__ == "inputprocess":
                self.setInput("seq1", upstream.getOutput("fq1"))
                self.setInput("seq2", upstream.getOutput("fq2"))
            elif upstream.__class__.__name__ == "adapterremoval":
                self.setInput("seq1", upstream.getOutput("pair1"))
                self.setInput("seq2", upstream.getOutput("pair2"))
            else:
                raise commonError(
                    "Parameter upstream must from inputprocess or adapterremoval."
                )

            self.checkInputFilePath()

        # set outputdir
        if upstream is None:
            if outputdir is None:
                self.setOutput(
                    "outputdir",
                    os.path.dirname(os.path.abspath(self.getInput("seq1")[0])),
                )
            else:
                self.setOutput("outputdir", outputdir)
        else:
            self.setOutput("outputdir", self.getStepFolderPath())

        # set ref, threads, paired
        if upstream is None:
            self.setParam("ref", ref)
            self.setParam("threads", threads)
            if paired:
                self.setParam("type", "paired")
            else:
                self.setParam("type", "single")

        else:
            self.setParam(
                "ref", os.path.join(Configure.getRefDir(), Configure.getGenome())
            )
            self.setParam("threads", Configure.getThreads())
            self.setParam("type", Configure.getType())

        # check reference for bowtie2
        self.bt2refcheck()

        # paired or single
        if self.getParam("type") == "paired":
            # generate base name
            prefix = []
            for seq1, seq2 in zip(self.getInput("seq1"), self.getInput("seq2")):
                prefix.append(self.getMaxFileNamePrefix(seq1, seq2))
            self.setParam("prefix", prefix)

            self.setParam(
                "outPrefix",
                [
                    os.path.join(self.getOutput("outputdir"), x)
                    for x in self.getParam("prefix")
                ],
            )

            if other_params is None:
                self.setParam("other_params", "")
            else:
                self.setParam("other_params", other_params)

            self.setOutput(
                "bamOutput", [x + ".bam" for x in self.getParam("outPrefix")]
            )
            self.setOutput(
                "unmapped-1", [x + ".unmapped.1.gz" for x in self.getParam("outPrefix")]
            )
            self.setOutput(
                "unmapped-2", [x + ".unmapped.2.gz" for x in self.getParam("outPrefix")]
            )
            self.setOutput(
                "bamOutput", [x + ".bam" for x in self.getParam("outPrefix")]
            )
            self.setOutput(
                "unmapbamOutput", [x + "-unmap.2bam" for x in self.getParam("outPrefix")]
            )

            self.setOutput(
                "sortumbamOutput", [x + "-unmap-sort.bam" for x in self.getParam("outPrefix")]
            )

            self.setOutput(
                "bedout", [x + ".bed" for x in self.getParam("outPrefix")]
            )
            
            if len(self.getInput("seq1")) == len(self.getInput("seq2")):
                multi_run_len = len(self.getInput("seq1"))
            else:
                raise commonError("Paired end Input files are not consistent.")

            all_cmd = []

            perl_script = """perl -ne 'chomp;@f=split " ";if($f[0] ne $f[3]){{next;}}$s=$f[1];$e=$f[5];if($f[8] eq "-"){{$s=$f[4];$e=$f[2];}}if($e>$s){{print "$f[0]\\t$s\\t$e\\t.\\t$f[7]\\t$f[8]\\n";}}' | """

            #perl_script = """perl -ne 'chomp;@f=split " ";if($f[0] ne $f[3]){{next;}}$s=$f[1];$e=$f[5];if($f[8] eq "-"){{$s=$f[4];$e=$f[2];}}if($e>$s){{print "$f[0]\\t$s\\t$e\\t$f[6]\\t$f[7]\\t$f[8]\\n";}}' | """

            for i in range(multi_run_len):
                tmp_cmd = self.cmdCreate(
                    [
                        "bwa mem",
                        "-t",
                        self.getParam("threads"),
                        "/mnt/sas/ref/hg19/v0/Homo_sapiens_assembly19.fasta",
                        self.getInput("seq1")[i],
                        self.getInput("seq2")[i],
                        "|",
                        "samblaster",
                        "|",
                        "samtools view -b",
                        "-@",
                        self.getParam("threads"),
                        "-",
                        ">",
                        self.getOutput("bamOutput")[i],

                    ]

                )
                tmp_cmd2 = self.cmdCreate(
                    [
                        "samtools view -b -f 4",
                        self.getOutput("bamOutput")[i], 
                        ">",
                        self.getOutput("unmapbamOutput")[i],

                    ]
                )
                tmp_cmd3 = self.cmdCreate(
                    [
                        "samtools sort",
                        "-n",
                        self.getOutput("unmapbamOutput")[i],
                        "-o",
                        self.getOutput("sortumbamOutput")[i], 
                    ]
                )
                tmp_cmd4 = self.cmdCreate(
                    [
                        "bedtools bamtofastq",
                        "-i",
                        self.getOutput("sortumbamOutput")[i],  
                        "-fq",
                        self.getOutput("unmapped-1")[i],  
                        "-fq2",
                        self.getOutput("unmapped-2")[i],  
                    ]
                )
                
                tmp_cmd6 = self.cmdCreate(
                    [
                        "samtools view -h -f 3 -F 3852 -q 30",
                        self.getOutput("bamOutput")[i],
                        "|",
                        "bamToBed -bedpe -mate1 -i stdin",
                        "|",
                        perl_script,
                        "sort-bed --max-mem 32G",
                        "-",
                        ">",
                        self.getOutput("bedout")[i],
                    ]
                )


                all_cmd.append(tmp_cmd)
                all_cmd.append(tmp_cmd2)
                all_cmd.append(tmp_cmd3)
                all_cmd.append(tmp_cmd4)
                all_cmd.append(tmp_cmd6)
                



        elif self.getParam("type") == "single":
            self.setParam(
                "prefix",
                [self.getMaxFileNamePrefixV2(x) for x in self.getInput("seq1")],
            )
            self.setParam(
                "outPrefix",
                [
                    os.path.join(self.getOutput("outputdir"), x)
                    for x in self.getParam("prefix")
                ],
            )

            if other_params is None:
                self.setParam("other_params", "")
            else:
                self.setParam("other_params", other_params)

            self.setOutput(
                "bamOutput", [x + ".bam" for x in self.getParam("outPrefix")]
            )

            self.setParam(
                "unmapped", [x + ".unmapped.gz" for x in self.getParam("outPrefix")]
            )

            multi_run_len = len(self.getInput("seq1"))

            all_cmd = []

            for i in range(multi_run_len):
                tmp_cmd = self.cmdCreate(
                    [
                        "bowtie2",
                        "-x",
                        self.getParam("ref"),
                        "-U",
                        self.getInput("seq1")[i],
                        "--un",
                        self.getParam("unmapped")[i],
                        self.getParam("other_params"),
                        "-p",
                        self.getParam("threads"),
                        "|",
                        "samtools view -b -S -@",
                        self.getParam("threads"),
                        "-o",
                        self.getOutput("bamOutput")[i],
                        "-",
                    ]
                )
                all_cmd.append(tmp_cmd)

        else:
            commonError("Wrong data tpye, must be 'single' or 'paired'!")

        finishFlag = self.stepInit(upstream)

        if not finishFlag:
            if verbose:
                self.run(all_cmd)
            else:
                self.multiRun(args=all_cmd, func=None, nCore=1)

        self.stepInfoRec(cmds=[all_cmd], finishFlag=finishFlag)

    # ref check
    def bt2refcheck(self,):
        extension = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
        bt2Ref = [self.getParam("ref") + x for x in extension]
        for filePath in bt2Ref:
            if not os.path.exists(filePath):
                raise commonError("Bowtie2 index file " + filePath + " don not exist!")
