# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 16:45:31 2019

@author: zhang
"""

import os
import time
import urllib.request
from multiprocessing import cpu_count
from .cfDNA_utils import commonError, un_gz, cmdCall, indexCheck
from .Configure import Configure
import glob
import math

__metaclass__ = type


class Configure2:
    __config = {
        "threads": 1,
        "genome": None,
        "refdir": None,
        "outdir": None,
        "data": None,
        "type": "paired",
        "JavaMem": "4G",
        "case": "case",
        "ctrl": "ctrl",
    }

    def __init__(self,):
        """
        threads: int, how many thread to use, default: 1.
        genome: str, which genome you want to use, 'hg19' or 'hg38'.
        refdir: reference folder for aligner (bowtie2 or bismark).
        outdir: overall result folder.
        data: data type, 'WGBS' or 'WGS'.
        type: data type, 'paired' or 'single'.
        "JavaMem": Java memory for every thred, default: 4G.
        case: case name for creating case specific folder.
        ctrl: control name for creating control specific folder.
        """
        raise commonError("Configure2 can not be initialized")

    # get configure names
    @classmethod
    def getConfigs(cls,):
        return cls.__config.keys()

    # get configure through name
    @classmethod
    def getConfig(cls, key):
        return cls.__config[key]

    # set configure through name
    @classmethod
    def setConfig(cls, key, val):
        if key == "threads":
            Configure2.setThreads(val)
        elif key == "genome":
            Configure2.setGenome(val)
        elif key == "outdir":
            Configure2.setOutDir(val)
        elif key == "refdir":
            Configure2.setRefDir(val)
        elif key == "data":
            Configure2.setData(val)
        elif key == "type":
            Configure2.setType(val)
        elif key == "JavaMem":
            Configure2.setJavaMem(val)
        elif key == "case":
            Configure2.setCase(val)
        elif key == "ctrl":
            Configure2.setCtrl(val)
        else:
            cls.__config[key] = val

    # set data
    @classmethod
    def setData(cls, val):
        cls.__config["data"] = val

    # get data
    @classmethod
    def getData(cls):
        return cls.__config["data"]

    # set type
    @classmethod
    def setType(cls, val):
        cls.__config["type"] = val

    # get type
    @classmethod
    def getType(cls):
        return cls.__config["type"]

    # set thread
    @classmethod
    def setThreads(cls, val):
        cls.__config["threads"] = val

    # get thread
    @classmethod
    def getThreads(cls):
        return cls.__config["threads"]

    # set JavaMem
    @classmethod
    def setJavaMem(cls, val):
        cls.__config["JavaMem"] = val

    # get JavaMem
    @classmethod
    def getJavaMem(cls):
        return cls.__config["JavaMem"]

    # set reference path
    @classmethod
    def setRefDir(cls, folderPath):
        Configure2.checkFolderPath(folderPath)
        cls.__config["refdir"] = folderPath

    # get reference path
    @classmethod
    def getRefDir(cls,):
        return cls.__config["refdir"]

    # set overall output directory and sub dir
    @classmethod
    def setOutDir(cls, folderPath):
        Configure2.checkFolderPath(folderPath)
        cls.__config["outdir"] = folderPath
        cls.__config["casedir"] = os.path.join(Configure2.getOutDir(), cls.__config["case"])
        cls.__config["ctrldir"] = os.path.join(Configure2.getOutDir(), cls.__config["ctrl"])

    # get overall output path
    @classmethod
    def getOutDir(cls,):
        return cls.__config["outdir"]

    # create intermediate, final and report folder
    @classmethod
    def pipeFolderInit(cls,):
        Configure2.configureCheck()
        if not os.path.exists(cls.__config["casedir"]):
            os.mkdir(cls.__config["casedir"])
        if not os.path.exists(cls.__config["ctrldir"]):
            os.mkdir(cls.__config["ctrldir"])
        Configure2.checkFolderPath(cls.__config["casedir"])
        Configure2.checkFolderPath(cls.__config["ctrldir"])

    # set genome falg
    @classmethod
    def setGenome(cls, val):
        cls.__config["genome"] = val

    # get genome falg
    @classmethod
    def getGenome(cls):
        return cls.__config["genome"]

    # set case
    @classmethod
    def setCase(cls, val):
        cls.__config["case"] = val

    # get case
    @classmethod
    def getCase(cls):
        return cls.__config["case"]

    # set thread
    @classmethod
    def setCtrl(cls, val):
        cls.__config["ctrl"] = val

    # get thread
    @classmethod
    def getCtrl(cls):
        return cls.__config["ctrl"]

    # get intermediate result path
    @classmethod
    def getTmpPath(cls, foldOrFileName):
        if isinstance(foldOrFileName, list):
            result = []
            for name in foldOrFileName:
                result.append(os.path.join(cls.getTmpDir(), name))
            return result
        else:
            return os.path.join(cls.getTmpDir(), foldOrFileName)

    # check folder legency, existence and accessibility
    @staticmethod
    def checkFolderPath(folderPath):
        if not os.path.isdir(os.path.abspath(folderPath)):
            raise commonError(folderPath + " is not an folder.")
        if not os.path.exists(folderPath):
            raise commonError(folderPath + " is not exist.")
        if not (os.access(folderPath, os.X_OK) and os.access(folderPath, os.W_OK)):
            raise commonError(folderPath + " is not accessible.")
        return True

    # check configure
    @classmethod
    def configureCheck(cls,):
        if Configure2.getType() is None:
            raise commonError("Please set type configure before using.")
        if Configure2.getData() is None:
            raise commonError("Please set data configure before using.")
        if Configure2.getGenome() is None:
            raise commonError("Please set genome configure before using.")
        if Configure2.getRefDir() is None:
            raise commonError("Please set reference configure before using.")
        if Configure2.getConfig("case") is None:
            raise commonError("Please set case configure before using.")
        if Configure2.getConfig("casedir") is None:
            raise commonError("Please set case configure before using.")
        if Configure2.getConfig("ctrl") is None:
            raise commonError("Please set ctrl configure before using.")
        if Configure2.getConfig("ctrldir") is None:
            raise commonError("Please set ctrl configure before using.")

    # check configure
    @classmethod
    def refCheck(cls, build=False):
        Configure2.configureCheck()
        Configure2.genomeRefCheck(build=build)
        Configure2.gitOverAllCheck(build=build)
        if Configure2.getData() == "WGBS":
            Configure2.bismkrefcheck(build)
            print("Background reference check finished!")
        elif Configure2.getData() == "WGS":
            Configure2.bt2refcheck(build)
            print("Background reference check finished!")
        else:
            print("No reference is specified.")

    # bismark ref check
    @classmethod
    def bismkrefcheck(cls, build):
        # check Bismark reference
        CTfiles = [
            os.path.join(Configure2.getRefDir(), "Bisulfite_Genome/CT_conversion/" + x)
            for x in [
                "BS_CT.1.bt2",
                "BS_CT.2.bt2",
                "BS_CT.3.bt2",
                "BS_CT.4.bt2",
                "BS_CT.rev.1.bt2",
                "BS_CT.rev.2.bt2",
                "genome_mfa.CT_conversion.fa",
            ]
        ]
        BAfiles = [
            os.path.join(Configure2.getRefDir(), "Bisulfite_Genome/GA_conversion/" + x)
            for x in [
                "BS_GA.1.bt2",
                "BS_GA.2.bt2",
                "BS_GA.3.bt2",
                "BS_GA.4.bt2",
                "BS_GA.rev.1.bt2",
                "BS_GA.rev.2.bt2",
                "genome_mfa.GA_conversion.fa",
            ]
        ]
        bismkRef = CTfiles + BAfiles
        if not all(map(os.path.exists, bismkRef)):
            print("Bismark index file do not exist or missing some files!")
            if build:
                cmdline = "bismark_genome_preparation " + Configure2.getRefDir()
                # if Configure2.getThreads() > 16:
                #     cmdline = (
                #         "bismark_genome_preparation --parallel "
                #         + str(16)
                #         + " "
                #         + Configure2.getRefDir()
                #     )
                # else:
                #     cmdline = "bismark_genome_preparation " + Configure2.getRefDir()
                print("Start building bismark reference......")
                print("Now, running " + cmdline)
                cmdCall(cmdline)
                print("Finished!")

    # bowtie2 ref check
    @classmethod
    def bt2refcheck(cls, build):
        # bowtie2 ref check
        extension = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
        bt2Ref = [os.path.join(Configure2.getRefDir(), Configure2.getGenome() + x) for x in extension]
        if not all(map(os.path.exists, bt2Ref)):
            print("Bowtie2 index file do not exist or missing some files!")
            if build:
                cmdline = (
                    "bowtie2-build -f --threads "
                    + str(Configure2.getThreads())
                    + " "
                    + Configure2.getConfig("genome.seq")
                    + " "
                    + os.path.join(Configure2.getRefDir(), Configure2.getGenome())
                )
                print("Start building Bowtie2 reference......")
                print("Now, running " + cmdline)
                cmdCall(cmdline)
                print("Finished!")

    # check genome reference
    @classmethod
    def genomeRefCheck(cls, build):
        Configure2.setConfig(
            "genome.seq", os.path.join(Configure2.getRefDir(), Configure2.getGenome() + ".fa"),
        )
        Configure2.setConfig(
            "genome.idx.fai", os.path.join(Configure2.getRefDir(), Configure2.getConfig("genome.seq") + ".fai"),
        )
        Configure2.setConfig(
            "genome.idx.dict", os.path.join(Configure2.getRefDir(), Configure2.getGenome() + ".dict"),
        )
        if not os.path.exists(Configure2.getConfig("genome.seq")):
            print("Reference file " + Configure2.getConfig("genome.seq") + " do not exist!")
            if build:
                url = (
                    "https://hgdownload.soe.ucsc.edu/goldenPath/"
                    + Configure2.getGenome()
                    + "/bigZips/"
                    + Configure2.getGenome()
                    + ".fa.gz"
                )
                print("Download from URL:" + url + "......")
                urllib.request.urlretrieve(
                    url, os.path.join(Configure2.getRefDir(), Configure2.getGenome() + ".fa.gz"),
                )
                print("Uncompressing......")
                un_gz(os.path.join(Configure2.getRefDir(), Configure2.getGenome() + ".fa.gz"))
                print("Finished!")

        if not os.path.exists(Configure2.getConfig("genome.idx.fai")):
            print("Reference file " + Configure2.getConfig("genome.idx.fai") + " do not exist!")
            if build:
                cmdline = "samtools faidx " + Configure2.getConfig("genome.seq")
                print("Start building .fai index for fasta reference......")
                print("Now, running " + cmdline)
                cmdCall(cmdline)
                print("Finished!")

        if not os.path.exists(Configure2.getConfig("genome.idx.dict")):
            print("Reference file " + Configure2.getConfig("genome.idx.dict") + " do not exist!")
            if build:
                cmdline = "gatk CreateSequenceDictionary --REFERENCE " + Configure2.getConfig("genome.seq")
                print("Start building .dict index for fasta reference......")
                print("Now, running " + cmdline)
                cmdCall(cmdline)
                print("Finished!")

    # check github.io file
    @classmethod
    def githubIOFile(cls, configureName, prefix, suffix, gitPath, build):
        fileName = prefix + Configure2.getGenome() + suffix
        fileNameGZ = fileName + ".gz"
        Configure2.setConfig(
            configureName, os.path.join(Configure2.getRefDir(), fileName),
        )
        if not os.path.exists(Configure2.getConfig(configureName)):
            print("Reference file " + Configure2.getConfig(configureName) + " do not exist!")
            if build:
                url = "https://honchkrow.github.io/cfDNAReferences/" + gitPath + "/" + fileNameGZ
                print("Download from URL:" + url + "......")
                urllib.request.urlretrieve(
                    url, os.path.join(Configure2.getRefDir(), fileNameGZ),
                )
                print("Uncompressing......")
                un_gz(os.path.join(Configure2.getRefDir(), fileNameGZ))
                print("Finished!")
                print("Now, waitting for next step......")
                time.sleep(10)

    # check github.io file
    @classmethod
    def gitOverAllCheck(cls, build):
        gitPath = Configure2.getGenome()
        Configure2.githubIOFile(
            configureName="chromSizes", prefix="", suffix=".chrom.sizes", gitPath=gitPath, build=build,
        )
        Configure2.githubIOFile(
            configureName="CpGisland", prefix="cpgIsland_", suffix=".bed", gitPath=gitPath, build=build,
        )
        Configure2.githubIOFile(
            configureName="cytoBand", prefix="cytoBand_", suffix=".txt", gitPath=gitPath, build=build,
        )
        Configure2.githubIOFile(
            configureName="OCF", prefix="OCF_", suffix=".bed", gitPath=gitPath, build=build,
        )
        Configure2.githubIOFile(
            configureName="PlasmaMarker", prefix="plasmaMarkers_", suffix=".txt", gitPath=gitPath, build=build,
        )
        Configure2.githubIOFile(
            configureName="Blacklist", prefix="", suffix="-blacklist.v2.bed", gitPath=gitPath, build=build,
        )
        Configure2.githubIOFile(
            configureName="Gaps", prefix="", suffix=".gaps.bed", gitPath=gitPath, build=build,
        )
        Configure2.githubIOFile(
            configureName="refFlat", prefix="refFlat_", suffix=".txt", gitPath=gitPath, build=build,
        )
        Configure2.githubIOFile(
            configureName="access-mappable", prefix="access-mappable.", suffix=".bed", gitPath=gitPath, build=build,
        )

    # additional function: check virus genome
    @classmethod
    def virusGenomeCheck(cls, folder=None, build=False):
        if folder is None:
            raise commonError('Parameter folder must not be "None"!')

        Configure2.setConfig("snv.folder", folder)

        cf_files = glob.glob(folder + "/*.cf")
        prefix = list(map(lambda x: os.path.basename(x).split(".")[0], cf_files))

        if len(set(prefix)) == 1:
            print("Virus reference " + prefix[0] + " are detected!")
            Configure2.setConfig(
                "virus.ref", os.path.join(folder, prefix[0]),
            )
            return True
        else:
            print("Can not find .cf files, searching NCBI taxonomy files......")

            if not build:
                raise commonError("NCBI taxonomy files need to be downloaded first!")
            else:
                cmdline1 = (
                    "centrifuge-download -o "
                    + os.path.join(folder, "taxonomy")
                    + " -P "
                    + str(math.ceil(Configure2.getThreads() / 4))
                    + " taxonomy"
                )
                cmdline2 = (
                    "centrifuge-download -o "
                    + os.path.join(folder, "library")
                    + " -P "
                    + str(math.ceil(Configure2.getThreads() / 4))
                    + ' -m -d "viral" refseq > '
                    + os.path.join(folder, "seqid2taxid.map")
                )
                cmdline3 = (
                    "cat "
                    + os.path.join(folder, "library/*/*.fna")
                    + " > "
                    + os.path.join(folder, "input-sequences.fna")
                )
                cmdline4 = (
                    "centrifuge-build -p "
                    + str(math.ceil(Configure2.getThreads() / 4))
                    + " --conversion-table "
                    + os.path.join(folder, "seqid2taxid.map")
                    + " --taxonomy-tree "
                    + os.path.join(folder, "taxonomy", "nodes.dmp")
                    + " --name-table "
                    + os.path.join(folder, "taxonomy", "names.dmp")
                    + " "
                    + os.path.join(folder, "input-sequences.fna")
                    + " "
                    + os.path.join(folder, "virus")
                )

                print("********Building Command********")
                print("Step 1:")
                print(cmdline1)
                print("Step 2:")
                print(cmdline2)
                print("Step 3:")
                print(cmdline3)
                print("Step 4:")
                print(cmdline4)
                print("********************************")

                print("Now, downloading NCBI taxonomy files......")
                print(cmdline1)
                cmdCall(cmdline1)

                print("Now, downloading virus genome files......")
                print(cmdline2)
                cmdCall(cmdline2)

                print("Now, merging virus genome files......")
                print(cmdline3)
                cmdCall(cmdline3)

                print("Now, building reference files......")
                print(cmdline4)
                cmdCall(cmdline4)

                print("DONE!")
                Configure2.setConfig(
                    "virus.ref", os.path.join(folder, "virus"),
                )

    # additional function: check SNV reference
    @classmethod
    def snvRefCheck(cls, folder=None, build=False):
        if folder is None:
            raise commonError('Parameter folder must not be "None"!')

        Configure2.setConfig("snv.folder", folder)

        if Configure2.getGenome() == "hg19":
            Configure2.setConfig(
                "snv.ref",
                {
                    "1": os.path.join(folder, "1000G_omni2.5.hg19.sites.vcf"),
                    "2": os.path.join(folder, "1000G_phase1.indels.hg19.sites.vcf"),
                    "3": os.path.join(folder, "1000G_phase1.snps.high_confidence.hg19.sites.vcf"),
                    "4": os.path.join(folder, "dbsnp_138.hg19.vcf"),
                    "5": os.path.join(folder, "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"),
                    "6": os.path.join(folder, "af-only-gnomad.raw.sites.hg19.vcf"),
                    "7": os.path.join(folder, "small_exac_common_3_hg19.SNP_biallelic.vcf"),
                    "8": os.path.join(folder, "somatic-hg19_Mutect2-WGS-panel.vcf"),
                },
            )

            # check for all reference files
            file_exist = list(map(os.path.exists, Configure2.getConfig("snv.ref").values()))

            if all(file_exist):
                print("hg19 SNV reference files are checked!")
            else:
                if not build:
                    raise commonError("SNV reference files are missing, please check reference files!")

                # build required files
                file_exist = list(map(os.path.exists, list(Configure2.getConfig("snv.ref").values())[0:5],))
                if not all(file_exist):
                    raise commonError("SNV reference files are missing, please check reference files!")
                else:  # build index files
                    for tmp_file in list(Configure2.getConfig("snv.ref").values())[0:5]:
                        if not indexCheck(tmp_file, ".idx"):
                            cmd_tmp = "gatk IndexFeatureFile --input " + tmp_file
                            print(cmd_tmp)
                            cmdCall(cmd_tmp)

                Configure2.setConfig(
                    "snv.tmp",
                    {
                        "1": os.path.join(folder, "af-only-gnomad.raw.sites.b37.vcf"),
                        "2": os.path.join(folder, "small_exac_common_3_b37.vcf"),
                        "3": os.path.join(folder, "b37tohg19.chain"),
                        "4": os.path.join(folder, "somatic-b37_Mutect2-WGS-panel-b37.vcf"),
                    },
                )

                if not all(list(map(os.path.exists, Configure2.getConfig("snv.tmp").values()))[0:3]):
                    raise commonError("SNV reference files are missing, please check reference files!")
                else:
                    for tmp_file in list(Configure2.getConfig("snv.tmp").values())[0:2]:
                        if not indexCheck(tmp_file, ".idx"):
                            cmd_tmp = "gatk IndexFeatureFile --input " + tmp_file
                            print(cmd_tmp)
                            cmdCall(cmd_tmp)

                if not os.path.exists(Configure2.getConfig("snv.ref")["6"]):
                    cmdline3 = (
                        "gatk LiftoverVcf "
                        + " --INPUT "
                        + Configure2.getConfig("snv.tmp")["1"]
                        + " --CHAIN "
                        + Configure2.getConfig("snv.tmp")["3"]
                        + " --OUTPUT "
                        + Configure2.getConfig("snv.ref")["6"]
                        + " --REFERENCE_SEQUENCE "
                        + Configure2.getConfig("genome.seq")
                        + " --REJECT "
                        + os.path.join(folder, "reject.variant.2.vcf")
                    )
                    print(cmdline3)
                    cmdCall(cmdline3)

                if not os.path.exists(Configure2.getConfig("snv.ref")["7"]):
                    cmdline1 = (
                        "gatk LiftoverVcf "
                        + " --CHAIN "
                        + Configure2.getConfig("snv.tmp")["3"]
                        + " --INPUT "
                        + Configure2.getConfig("snv.tmp")["2"]
                        + " --OUTPUT "
                        + os.path.join(folder, "small_exac_common_3_hg19.vcf")
                        + " --REFERENCE_SEQUENCE "
                        + Configure2.getConfig("genome.seq")
                        + " --REJECT "
                        + os.path.join(folder, "reject.variant.1.vcf")
                    )
                    print(cmdline1)
                    cmdCall(cmdline1)

                    cmdline2 = (
                        "gatk SelectVariants "
                        + " -R "
                        + Configure2.getConfig("genome.seq")
                        + " -V "
                        + os.path.join(folder, "small_exac_common_3_hg19.vcf")
                        + " --select-type-to-include SNP --restrict-alleles-to BIALLELIC "
                        + " -O "
                        + Configure2.getConfig("snv.ref")["7"]
                    )
                    print(cmdline2)
                    cmdCall(cmdline2)

                # check for files 1~7
                file_exist = list(map(os.path.exists, list(Configure2.getConfig("snv.ref").values())[0:7],))
                if all(file_exist):
                    print("SNV reference for case-control study build finished!")
                else:
                    raise commonError("SNV reference files are missing, please check reference files!")

                if os.path.exists(Configure2.getConfig("snv.ref")["8"]):
                    print("SNV PON reference for single group study is found!")
                else:
                    print("SNV PON reference file somatic-hg19_Mutect2-WGS-panel.vcf is not found!")
                    print("Program try to build it......")

                    if not os.path.exists(Configure2.getConfig("snv.tmp")["4"]):
                        print(
                            "Note: file somatic-b37_Mutect2-WGS-panel-b37.vcf is not found in snv reference folder, please download it from GATK resource."
                        )
                        print("PON file will not be constructed!")
                        Configure2.getConfig("snv.ref")["8"] = None
                    else:
                        print("Note: file somatic-b37_Mutect2-WGS-panel-b37.vcf is found, construct hg19 version.")

                        if not indexCheck(Configure2.getConfig("snv.tmp")["4"], ".idx"):
                            cmd_tmp = "gatk IndexFeatureFile --input " + Configure2.getConfig("snv.tmp")["4"]
                            print(cmd_tmp)
                            cmdCall(cmd_tmp)

                        cmdline4 = (
                            "gatk LiftoverVcf "
                            + " --INPUT "
                            + Configure2.getConfig("snv.tmp")["4"]
                            + " --CHAIN "
                            + Configure2.getConfig("snv.tmp")["3"]
                            + " --OUTPUT "
                            + Configure2.getConfig("snv.ref")["8"]
                            + " --REFERENCE_SEQUENCE "
                            + Configure2.getConfig("genome.seq")
                            + " --REJECT "
                            + os.path.join(folder, "reject.variant.3.vcf")
                        )
                        print(cmdline4)
                        cmdCall(cmdline4)

        elif Configure2.getGenome() == "hg38":
            Configure2.setConfig(
                "snv.ref",
                {
                    "1": os.path.join(folder, "1000G_omni2.5.hg38.vcf"),
                    "2": os.path.join(folder, "1000G_phase1.snps.high_confidence.hg38.vcf"),
                    "3": os.path.join(folder, "dbsnp_146.hg38.vcf"),
                    "4": os.path.join(folder, "hapmap_3.3.hg38.vcf"),
                    "5": os.path.join(folder, "Mills_and_1000G_gold_standard.indels.hg38.vcf"),
                    "6": os.path.join(folder, "af-only-gnomad.hg38.vcf"),
                    "7": os.path.join(folder, "small_exac_common_3_hg38.SNP_biallelic.vcf"),
                    "8": os.path.join(folder, "somatic-hg38_1000g_pon.hg38.vcf"),
                },
            )

            # check for all reference files
            file_exist = list(map(os.path.exists, Configure2.getConfig("snv.ref").values()))

            if all(file_exist):
                print("hg38 SNV reference files are checked!")
            else:
                if not build:
                    raise commonError("SNV reference files are missing, please check reference files!")

                # build required files
                file_exist = list(map(os.path.exists, list(Configure2.getConfig("snv.ref").values())[0:6],))
                if not all(file_exist):
                    raise commonError("SNV reference files are missing, please check reference files!")
                else:  # build index files
                    for tmp_file in list(Configure2.getConfig("snv.ref").values())[0:6]:
                        if not indexCheck(tmp_file, ".idx"):
                            cmd_tmp = "gatk IndexFeatureFile --input " + tmp_file
                            print(cmd_tmp)
                            cmdCall(cmd_tmp)

                Configure2.setConfig(
                    "snv.tmp", {"1": os.path.join(folder, "small_exac_common_3.hg38.vcf")},
                )

                if not all(list(map(os.path.exists, Configure2.getConfig("snv.tmp").values()))):
                    raise commonError("SNV reference files are missing, please check reference files!")
                else:
                    for tmp_file in list(Configure2.getConfig("snv.tmp").values()):
                        if not indexCheck(tmp_file, ".idx"):
                            cmd_tmp = "gatk IndexFeatureFile --input " + tmp_file
                            print(cmd_tmp)
                            cmdCall(cmd_tmp)

                if not os.path.exists(Configure2.getConfig("snv.ref")["7"]):
                    cmdline1 = (
                        "gatk SelectVariants "
                        + " -R "
                        + Configure2.getConfig("genome.seq")
                        + " -V "
                        + Configure2.getConfig("snv.tmp")["1"]
                        + " --select-type-to-include SNP --restrict-alleles-to BIALLELIC "
                        + " -O "
                        + Configure2.getConfig("snv.ref")["7"]
                    )
                    print(cmdline1)
                    cmdCall(cmdline1)

                # check for files 1~7
                file_exist = list(map(os.path.exists, list(Configure2.getConfig("snv.ref").values())[0:7],))
                if all(file_exist):
                    print("SNV reference for case-control study build finished!")
                else:
                    raise commonError("SNV reference files are missing, please check reference files!")

                if os.path.exists(Configure2.getConfig("snv.ref")["8"]):
                    print("SNV PON reference for single group study is found!")
                    if not indexCheck(Configure2.getConfig("snv.ref")["8"], ".idx"):
                        cmd_tmp = "gatk IndexFeatureFile --input " + Configure2.getConfig("snv.ref")["8"]
                        print(cmd_tmp)
                        cmdCall(cmd_tmp)
                else:
                    print("SNV PON reference file somatic-hg38_1000g_pon.hg38.vcf is not found!")
                    Configure2.getConfig("snv.ref")["8"] = None


def switchConfigure(confName=None):
    """
    Switch Configure for case and control, these two situation have different output directory
    parameter confName: one of the Configure name defined in Configure2 (case and ctrl)
    """
    Configure.setData(Configure2.getData())
    Configure.setType(Configure2.getType())
    Configure.setThreads(Configure2.getThreads())
    Configure.setGenome(Configure2.getGenome())
    Configure.setRefDir(Configure2.getRefDir())
    Configure.setJavaMem(Configure2.getJavaMem())
    if confName == Configure2.getCase():
        Configure.setOutDir(Configure2.getConfig("casedir"))
    elif confName == Configure2.getCtrl():
        Configure.setOutDir(Configure2.getConfig("ctrldir"))
    else:
        raise commonError("There is no Configure environment named " + confName + "!")

    Configure.pipeFolderInit()
    Configure.refCheck()


def pipeConfigure2(
    threads=(cpu_count() / 2),
    genome=None,
    refdir=None,
    outdir=None,
    data=None,
    type=None,
    JavaMem=None,
    case=None,
    ctrl=None,
    build=False,
):
    """
    This function is used for setting Configures.
    Note: This function is designed for case control comparison.

    pipeConfigure2(threads=(cpu_count() / 2), genome=None, refdir=None, outdir=None,
                   data=None, type=None, case=None, ctrl=None, build=False,)
    {P}arameters:
        threads: int, how many thread to use, default: 1.
        genome: str, which genome you want to use, 'hg19' or 'hg38'.
        refdir: reference folder for aligner (bowtie2 or bismark).
        outdir: Overall result folder.
        data: data type, 'WGBS' or 'WGS'.
        type: data type, 'paired' or 'single'.
        JavaMem: Java memory for every thred, "10G" like.
        case: case NAME for creating case specific folder.
        ctrl: control NAME for creating control specific folder.
    """
    Configure2.setData(data)
    Configure2.setType(type)
    Configure2.setThreads(threads)
    Configure2.setGenome(genome)
    Configure2.setRefDir(refdir)
    Configure2.setJavaMem(JavaMem)
    Configure2.setCase(case)
    Configure2.setCtrl(ctrl)
    Configure2.setOutDir(outdir)
    Configure2.pipeFolderInit()
    Configure2.refCheck(build=build)
