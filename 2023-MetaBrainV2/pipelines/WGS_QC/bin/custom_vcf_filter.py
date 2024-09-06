#!/opt/conda/envs/eQTLGenPopAssign/bin/python

"""
File:         custom_vcf_filter.py
Created:      2022/10/17
Last Changed: 2024/01/11
Author:       H.J.Westra (Edited by M.Vochteloo)
"""

# Standard imports.
import argparse
import gzip
import sys


# Third party imports.

# Local application imports.

# Metadata
__program__ = "Custom VCF Filter"
__author__ = "Harm-Jan Westra"
__maintainer__ = "Harm-Jan Westra & Martijn Vochteloo"
__email__ = "h.j.westra@umcg.nl"
__license__ = "NA"
__version__ = 5.2
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)

"""
Syntax: 

./custom_vcf_filter.py \
    --input /groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/2022-11-03-FilteredGenotypes/1-bcftools_norm/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_X.recalibrated_variants_norm.vcf.gz \
    --sex /groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/AMP_AD_sexdata.csv \
    --output /groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/2022-11-03-FilteredGenotypes/2-custom_vcffilter/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_X.recalibrated_variants_norm_vcffilter
    
./custom_vcf_filter.py \
    --input /groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/2022-11-03-FilteredGenotypes/1-bcftools_norm/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_22.recalibrated_variants_norm.vcf.gz \
    --sex /groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/AMP_AD_sexdata.csv \
    --output /groups/umcg-biogen/tmp01/input/processeddata/single-cell/AMP-AD/2022-11-03-FilteredGenotypes/2-custom_vcffilter/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_22.recalibrated_variants_norm_vcffilter_test
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.input_path = getattr(arguments, 'input')
        self.output_filename = getattr(arguments, 'output')
        self.sample_filter_set = getattr(arguments, 'samples')
        if isinstance(self.sample_filter_set, list):
            self.sample_filter_set = set(self.sample_filter_set)
        self.sex_path = getattr(arguments, 'sex')
        self.thresh_gq = getattr(arguments, 'genotype_quality')
        self.thresh_ab_lower = getattr(arguments, 'allelic_balance_lower')
        self.thresh_ab_upper = getattr(arguments, 'allelic_balance_upper')
        self.thresh_ib = getattr(arguments, 'inbreeding_coefficient')
        self.check_vqsr_snv = getattr(arguments, 'no_snv_vqsr_check')
        self.thresh_vqsr_snv = getattr(arguments, 'vqsr_snv')
        self.check_vqsr_indel = getattr(arguments, 'no_indel_vqsr_check')
        self.thresh_vqsr_indel = getattr(arguments, 'vqsr_indel')
        self.remove_multiallelic = getattr(arguments, 'remove_multialleic')
        self.remove_non_pass_snv = getattr(arguments, 'remove_non_pass_snv')
        self.remove_non_pass_indel = getattr(arguments, 'remove_non_pass_indel')
        self.remove_low_complexity = getattr(arguments, 'remove_low_complexity')
        self.thresh_maf = getattr(arguments, 'minor_allele_frequency')
        self.thresh_cr = getattr(arguments, 'call_rate')
        self.thresh_hwe = getattr(arguments, 'hardy_weinberg_equilibrium')
        self.thresh_dp = getattr(arguments, 'filtered_depth')
        self.strip_info_col = getattr(arguments, 'keep_info_column')
        self.ignore_homref_stats = getattr(arguments, 'ignore_homref_stats')
        self.replace_poor_quality_genotypes = getattr(arguments, 'replace_poor_quality_genotypes')


        self.debug = False
        self.stopafterlines = 50000

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add other arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("-i",
                            "--input",
                            type=str,
                            required=True,
                            help="The input VCF file.")
        parser.add_argument("-o",
                            "--output",
                            type=str,
                            required=True,
                            help="The output filename prefix.")
        parser.add_argument("-s",
                            "--samples",
                            nargs="*",
                            type=str,
                            required=False,
                            help="The samples to select.")
        parser.add_argument("--sex",
                            type=str,
                            required=False,
                            help="The sample-sex data file.")
        parser.add_argument("-gq",
                            "--genotype_quality",
                            type=float,
                            default=20.,
                            help="The genotype quality threshold. "
                                 "Default: 20.")
        parser.add_argument("-abl",
                            "--allelic_balance_lower",
                            type=float,
                            default=0.2,
                            help="The allelic balance lower-bound threshold. "
                                 "Default: 0.2.")
        parser.add_argument("-abu",
                            "--allelic_balance_upper",
                            type=float,
                            default=0.8,
                            help="The allelic balance upper-bound threshold. "
                                 "Default: 0.8.")
        parser.add_argument("-ib",
                            "--inbreeding_coefficient",
                            type=float,
                            default=-0.3,
                            help="The inbreeding coefficient threshold. "
                                 "Default: -0.3.")
        parser.add_argument("--no_snv_vqsr_check",
                            action='store_false',
                            help="Add this flag to turn off the VQSR check for SNV's."
                                 " Default: check enabled.")
        parser.add_argument("--vqsr_snv",
                            type=float,
                            default=99.8,
                            help="The variant_quality_score_recalibration "
                                 "threshold. Default: 99.8.")
        parser.add_argument("--no_indel_vqsr_check",
                            action='store_false',
                            help="Add this flag to turn off the VQSR check for indel's."
                                 " Default: check enabled.")
        parser.add_argument("--vqsr_indel",
                            type=float,
                            default=99.95,
                            help="The variant_quality_score_recalibration "
                                 "threshold. Default: 99.95.")
        parser.add_argument("--remove_multialleic",
                            action='store_false',
                            help="Add this flag to remove multiallelic "
                                 "variants. Default: False.")
        parser.add_argument("--remove_non_pass_snv",
                            action='store_false',
                            help="Add this flag to remove SNV's without PASS "
                                 "label. Default: False.")
        parser.add_argument("--remove_non_pass_indel",
                            action='store_false',
                            help="Add this flag to remove indel's without PASS "
                                 "label. Default: False.")
        parser.add_argument("--remove_low_complexity",
                            action='store_false',
                            help="Add this flag to remove low complexity [not implemented]"
                                 "variants. Default: False.")
        parser.add_argument("-maf",
                            "--minor_allele_frequency",
                            type=float,
                            default=0.01,
                            help="The minor allele frequency threshold. "
                                 "Default: 0.01.")
        parser.add_argument("-cr",
                            "--call_rate",
                            type=float,
                            default=0.99,
                            help="The call rate threshold. "
                                 "Default: 0.99.")
        parser.add_argument("-hwe",
                            "--hardy_weinberg_equilibrium",
                            type=float,
                            default=1E-6,
                            help="The hardy weinberg equilibrium threshold. "
                                 "Default: 1e-6.")
        parser.add_argument("-dp",
                            "--filtered_depth",
                            type=float,
                            default=10.,
                            help="The filtered depth threshold. "
                                 "Default: 10.")
        parser.add_argument("--keep_info_column",
                            action='store_false',
                            help="Add this flag to keep the INFO column. "
                                 "Default: False.")
        parser.add_argument("--ignore_homref_stats",
                            action='store_true',
                            help="Some VCF files have no DP, AB, or AD information for homozygous reference calls. Exclude such calls from QC tests. "
                                 "Default: True.")
        parser.add_argument("--replace_poor_quality_genotypes",
                            action='store_true',
                            help="Replace poor quality genotypes with ./."
                                 "Default: False.")
        return parser.parse_args()

    def start(self):
        self.print_arguments()

        ismale_dict = None
        if self.sex_path is not None:
            ismale_dict = {}
            fs = open(self.sex_path, 'r')
            for line in fs:
                sample, sex = line.strip("\n").split(",")
                if sex == "M":
                    ismale_dict[sample] = True
                elif sex == "F":
                    ismale_dict[sample] = False
                else:
                    print("Unexpected input in sex file.")
                    exit()

        fh = gzip.open(self.input_path, 'rt')
        fho = gzip.open(self.output_filename + "-filtered.vcf.gz", 'wt')
        fhlog = gzip.open(self.output_filename + "-filtered.log.gz", 'wt')
        fhlog.write("Id\tReason\tPreFilterStats\tPostFilterStats\n")

        sample_mask = []
        sample_male_mask = []

        linectr = 0
        lineswritten = 0
        for line in fh:
            if line.startswith("#CHROM"):
                outln = "##{} tresh_GQ ={};"\
                            "tresh_AB_lower={};"\
                            "tresh_AB_upper={};"\
                            "tresh_IB={};"\
                            "check_VQSR_SNV={};"\
                            "tresh_VQSR_SNV={};"\
                            "check_VQSR_Indel={};"\
                            "tresh_VQSR_Indel={};"\
                            "remove_multiAllelic={};"\
                            "remove_nonPASS_indel={};"\
                            "remove_nonPASS_SNV ={};" \
                            "remove_lowcomplexity={};"\
                            "tresh_MAF ={};" \
                            "tresh_HWE={};"\
                            "tresh_CR={};" \
                            "tresh_DP={};"\
                            "ignore_homref_stats={};"\
                            "replace_poor_quality_genotypes={}\n".format(self.input_path,
                                                self.thresh_gq,
                                                self.thresh_ab_lower,
                                                self.thresh_ab_upper,
                                                self.thresh_ib,
                                                self.check_vqsr_snv,
                                                self.thresh_vqsr_snv,
                                                self.check_vqsr_indel,
                                                self.thresh_vqsr_indel,
                                                self.remove_multiallelic,
                                                self.remove_non_pass_indel,
                                                self.remove_non_pass_snv,
                                                self.remove_low_complexity,
                                                self.thresh_maf,
                                                self.thresh_hwe,
                                                self.thresh_cr,
                                                self.thresh_dp,
                                                self.ignore_homref_stats,
                                                self.replace_poor_quality_genotypes
                                                )
                fho.write(outln)

                # header line with samples
                elems = line.strip().split("\t")
                elems[6] = "FILTER"
                elems[7] = "INFO"
                headerout = "\t".join(elems[0:9])
                n_samples_included = 0
                n_samples_excluded = 0
                for i in range(9, len(elems)):
                    sample = elems[i]
                    if self.sample_filter_set is None or sample in self.sample_filter_set:
                        headerout += "\t" + sample
                        sample_mask.append(True)
                        sample_male_mask.append(None if ismale_dict is None else ismale_dict[sample])
                        n_samples_included += 1
                    else:
                        sample_mask.append(False)
                        n_samples_excluded += 1
                print("{:,} samples selected from VCF header, {:,} excluded".format(n_samples_included, n_samples_excluded))
                fho.write(headerout + "\n")
                lineswritten += 1
            elif line.startswith("#"):
                if line.startswith("##FORMAT") and line[13:15] not in ["GT", "DP", "AD", "AB", "GQ"]:
                    continue
                if line.startswith("##INFO") and self.strip_info_col:
                    continue
                fho.write(line)
                lineswritten += 1
            else:
                # this can now be potentially multithreaded?
                parsed = self.parse_line(
                    line_number=linectr,
                    line=line,
                    sample_mask=sample_mask,
                    sample_male_mask=sample_male_mask
                )  # returns: [line_number, False, logid+"\PASSQC\t"+stats[1], outln]
                fhlog.write(parsed[2])
                if parsed[1]:
                    fho.write(parsed[3])
                    lineswritten += 1
            linectr += 1
            if self.debug and linectr == self.stopafterlines:
                print("DEBUG: stopped after " + str(linectr) + " lines")
                break
            if linectr % 10000 == 0:
                print("{:,} lines parsed, {:,} written".format(linectr, lineswritten), end='\r')
                fho.flush()
                fhlog.flush()
        print("{:,} lines parsed, {:,} written".format(linectr, lineswritten), end='\n')
        print("Done. How about that!")
        fh.close()
        fho.close()
        fhlog.close()

    def parse_line(self, line_number, line, sample_mask, sample_male_mask):
        # line with genotype data
        elems = line.strip().split("\t", 10)

        # CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
        # 0     1       2       3       4       5       6       7       8
        chr = elems[0]
        pos = elems[1]

        id = elems[2]
        logid = chr + ":" + pos + ":" + id
        ref = elems[3]
        alt = elems[4].split(",")

        # check if this is a multiAllelic
        if self.remove_multiallelic and len(alt) > 1:
            return [line_number, False, logid + "\tMultiAllelic\t-\t-\n"]

        # check if this variant is an indel
        is_indel = False
        if len(ref) > 1:
            is_indel = True

        if not is_indel:
            for a in alt:
                if len(a) > 1:
                    is_indel = True

        filter = elems[6]
        # check VQSR first
        if is_indel:
            if filter.startswith("VQSR") and self.check_vqsr_indel:
                if not self.parse_vqsr(filter, is_indel):
                    return [line_number, False, logid + "\tIndelBelowVQSR\t-\t-\n"]
            elif self.remove_non_pass_indel and filter != "PASS":
                return [line_number, False, logid + "\tIndelNonPass\t-\t-\n"]
        else:
            if filter.startswith("VQSR") and self.check_vqsr_snv:
                if not self.parse_vqsr(filter, is_indel):
                    return [line_number, False, logid + "\tSNVBelowVQSR\t-\t-\n"]
            elif self.remove_non_pass_snv and filter != "PASS":
                return [line_number, False, logid + "\tSNVNonPass\t-\t-\n"]

        # parse info string
        infoelems = elems[7].split(";")
        for elem in infoelems:
            vals = elem.split("=")
            if len(vals) > 1:
                key = vals[0]
                value = vals[1]
                if key == "InbreedingCoeff":
                    if value == ".":
                        return [line_number, False, logid + "\tIncorrectInbreedingCoeff:" + value + "\t-\t-\n"]
                    valuef = float(value)
                    if valuef < self.thresh_ib:
                        return [line_number, False, logid + "\tBelowInbreedingCoeff:" + value + "\t-\t-\n"]
            # other filters based on info column can be implemented here...

        # parse genotypes
        elems = line.strip().split("\t")
        format = elems[8].split(":")
        gtcol = -1
        adcol = -1
        dpcol = -1
        gqcol = -1
        plcol = -1
        for c in range(0, len(format)):
            if format[c] == "GT":
                gtcol = c
            elif format[c] == "AD":
                adcol = c
            elif format[c] == "DP":
                dpcol = c
            elif format[c] == "GQ":
                gqcol = c
            elif format[c] == "PL":
                plcol = c

        format_out = ""
        if gtcol > -1:
            format_out = "GT"
        if dpcol > -1:
            format_out += ":DP"
        if adcol > -1:
            format_out += ":AD"
            format_out += ":AB"
        if gqcol > -1:
            format_out += ":GQ"
        if plcol > -1:
            format_out += ":PL"

        # only continue if there are genotypes to parse
        if gtcol < 0:
            return [line_number, False, logid + "\tNoGTCol\t-\t-\n"]

        # get current dosages
        # store split columns temporarily
        sampledataelems = [None] * len(elems)
        dosages_pre_filter = []
        dosage_pre_filter = None
        sampledosages = [-1] * len(elems)
        cctr = 0

        for i in range(9, len(elems)):
            if sample_mask[cctr]:
                # split the sample columns
                sampledata = elems[i].split(":")
                sampledataelems[i] = sampledata
                gt = sampledata[gtcol].split("/")
                if len(gt) < 2:
                    # less than two alleles; perhaps the data is phased
                    gt = sampledata[gtcol].split("|")
                if len(gt) < 2 or len(gt) > 2:
                    # malformatted genotypes
                    # more than two alleles.. skip variant?
                    # fhlog.write(headerInfo[1]+"\tPolyPloid?\t-\n"
                    dosages_pre_filter.append(-1)
                else:
                    if gt[0] == ".":
                        dosages_pre_filter.append(-1)
                    else:
                        dosage_pre_filter = int(gt[0]) + int(gt[1])
                        dosages_pre_filter.append(dosage_pre_filter)
                        sampledosages[i] = dosage_pre_filter
            cctr += 1

        # calculate allele frequency etc before filtering
        stats = self.get_stats(
            dosages=dosages_pre_filter,
            sample_male_mask=sample_male_mask,
            chr=chr
        )
        if not stats[0]:
            return [line_number, False, logid + "\tFailedPrefilterVarStats\t" + stats[1] + "\t-\n"]

        # check whether variant becomes monomorphic after filtering poor calls
        # iterate samples
        dosages_post_filter = []
        sampleinfos = []
        nr_genotypes_replaced = 0
        average_depth = 0
        average_depth_calls = 0
        poor_dp = 0
        poor_gq = 0
        poor_ab_hom_a = 0
        poor_ab_hom_b = 0
        poor_ab_het = 0
        cctr = 0
        for i in range(9, len(elems)):
            sampledata = sampledataelems[i]
            if sampledata is not None:
                dosage_pre_filter = sampledosages[i]
                dosage_post_filter = dosage_pre_filter
                ab = -1
                dp = 0
                # some VCF files have no information specifically for homozygous reference calls. 
                # for those calls it makes sense to ignore the quality statistics
                if self.ignore_homref_stats and dosage_post_filter == 0:
                    dosages_post_filter.append(dosage_post_filter)
                else:   
                    # check genotype qual if the genotype is not missing
                    if dpcol > -1 and dpcol < len(sampledata):  # ugh, VCF really allows all kind of stuff in their genotype format... even if you parse the FORMAT column, there's no telling what is actually provided per sample..
                        # determine read depth
                        if sampledata[dpcol] == ".":
                            dp = 0
                        else:
                            dp = int(sampledata[dpcol])
                            average_depth += dp
                            average_depth_calls += 1
                        if dp < self.thresh_dp:
                            dosage_post_filter = -1
                            poor_dp += 1
                            nr_genotypes_replaced += 1
                    if gqcol > -1 and dosage_post_filter > -1 and gqcol < len(sampledata):
                        gq = float(sampledata[gqcol])
                        if gq < self.thresh_gq:
                            if dosage_post_filter == 0 and gq == 0 and dp >= self.thresh_dp:
                                # bypass potential error in gVCF merging, see: https://github.com/broadinstitute/gatk/issues/5445
                                # best way would to also parse the PL field, if available, and check whether the homRef and het fields are both 0 (unlikely for a homRef call)
                                pass
                            else:
                                dosage_post_filter = -1
                                nr_genotypes_replaced += 1
                                poor_gq += 1
                    # check allelic balance
                    # only need to check this for hets, I guess
                    if dosage_post_filter > -1 and adcol > -1 and adcol < len(sampledata):
                        ad = sampledata[adcol].split(",")
                        ad1 = float(ad[0])
                        ad2 = float(ad[1])
                        if ad1 + ad2 == 0:
                            dosage_post_filter = -1
                            nr_genotypes_replaced += 1
                        else:
                            ab = ad2 / (ad1 + ad2)
                            # data is already corrected for DP, so no need to check again.
                            if dosage_post_filter == 0:
                                if ab > self.thresh_ab_lower:
                                    dosage_post_filter = -1
                                    nr_genotypes_replaced += 1
                                    poor_ab_hom_a += 1
                            if dosage_post_filter == 1:
                                if ab < self.thresh_ab_lower or ab > self.thresh_ab_upper:
                                    dosage_post_filter = -1
                                    nr_genotypes_replaced += 1
                                    poor_ab_het += 1
                            if dosage_post_filter == 2:
                                if ab < self.thresh_ab_upper:
                                    dosage_post_filter = -1
                                    nr_genotypes_replaced += 1
                                    poor_ab_hom_b += 1
                    dosages_post_filter.append(dosage_post_filter)

                sampleinfo = []
                gtout = sampledata[gtcol]

                if self.replace_poor_quality_genotypes:
                    if dosage_post_filter == -1:
                        gtout = "./."
                else:
                    if dosage_pre_filter == -1:
                        gtout = "./."

                sampleinfo.append(gtout)

                # append original information, but note that missing genotypes
                # do not always conform to the FORMAT string
                if dpcol > -1 and dpcol < len(sampledata):
                    sampleinfo.append(sampledata[dpcol])
                if adcol > -1 and adcol < len(sampledata):
                    sampleinfo.append(sampledata[adcol])
                    sampleinfo.append(str(round(ab, 2)))
                if gqcol > -1 and gqcol < len(sampledata):
                    sampleinfo.append(sampledata[gqcol])
                if plcol > -1 and plcol < len(sampledata):
                    sampleinfo.append(sampledata[plcol])

                sampleinfos.append(":".join(sampleinfo))
            cctr += 1

        # average read depth
        if average_depth_calls > 0:
            average_depth /= average_depth_calls
        else:
            average_depth = 0

        # variant level quals
        # calculate allele frequency etc before filtering
        stats_post_filter = self.get_stats(
            dosages=dosages_post_filter,
            sample_male_mask=sample_male_mask,
            chr=chr
        )

        _, _, _, _, _, maf_postfilter = self.get_maf(
            dosages=dosages_post_filter,
            sample_male_mask=sample_male_mask,
            chr=chr
        )

        if maf_postfilter == 0:
            logoutln = "{}\tMonomorphicPostFilter\t-\t{};NrGenosReplaced:{}" \
                       ";PoorDP:{};AvgDP:{:.2f} ({} calls);PoorGQ:{}" \
                       ";PoorABHomA:{};PoorABHomB:{};PoorABHet:{}" \
                       "\n".format(logid,
                                   stats_post_filter[1],
                                   nr_genotypes_replaced,
                                   poor_dp,
                                   average_depth,
                                   average_depth_calls,
                                   poor_gq,
                                   poor_ab_hom_a,
                                   poor_ab_hom_b,
                                   poor_ab_het)
            return [line_number, False, logoutln]

        format_out = ""
        if gtcol > -1:
            format_out = "GT"
        if dpcol > -1:
            format_out += ":DP"
        if adcol > -1:
            format_out += ":AD"
            format_out += ":AB"
        if gqcol > -1:
            format_out += ":GQ"
        if plcol > -1:
            format_out += ":PL"

        elems[8] = format_out

        if self.strip_info_col:
            elems[7] = "."
        outln = "\t".join(elems[0:9])  # construct line header
        # write variant, somehow
        outln += "\t" + "\t".join(sampleinfos) + "\n"

        var_ok = stats_post_filter[0]
        if var_ok:
            return [line_number, True, logid + "\tPASSQC\t" + stats[1] + "\t" + stats_post_filter[1] + "\n", outln]
        else:
            return [line_number, False, logid + "\tFailQCPostFilter\t" + stats[1] + "\t" + stats_post_filter[1] + "\n", outln]

    def parse_vqsr(self, vqsrstring, is_indel):
        if not vqsrstring.startswith("VQSR"):
            return True
        vqsrstring = vqsrstring.replace("VQSRTrancheINDEL", "")
        vqsrstring = vqsrstring.replace("VQSRTrancheSNP", "")
        vqsrstring = vqsrstring.replace("+", "")
        velems = vqsrstring.split("to")
        lower = float(velems[0])
        if is_indel:
            return lower >= self.thresh_vqsr_indel
        else:
            return lower >= self.thresh_vqsr_snv

    def get_stats(self, dosages, sample_male_mask, chr):
        nrcalled, nrdosages, nrhoma, nrhets, nrhomb, maf = self.get_maf(
            dosages=dosages,
            sample_male_mask=sample_male_mask,
            chr=chr
        )
        if nrcalled == 0:
            return [False, "CR=0.000;MAF=0.000;HWE=-1.000e+00"]

        var_ok = True
        callrate = nrcalled / nrdosages
        hwe = -1
        if callrate >= self.thresh_cr and maf >= self.thresh_maf:
            hwe = self.calculate_hwe(
                obs_hets=nrhets,
                obs_hom1=nrhoma,
                obs_hom2=nrhomb
            )  # only calculate HWE when MAF and CR above thresholds
            if hwe < self.thresh_hwe:
                var_ok = False
        else:
            var_ok = False
        statout = "CR={:.3f};MAF={:.3f};HWE={:0.3e}".format(callrate, maf, hwe)

        return [var_ok, statout]

    @staticmethod
    def get_maf(dosages, sample_male_mask, chr):
        nrhoma = 0
        nrhets = 0
        nrhomb = 0
        nrmiss = 0

        for dosage, is_male in zip(dosages, sample_male_mask):
            if chr.lstrip("chr") in ["X", "y"] and (is_male is None or is_male):
                continue

            if dosage == 0:
                nrhoma += 1
            elif dosage == 1:
                nrhets += 1
            elif dosage == 2:
                nrhomb += 1
            else:
                nrmiss += 1

        nrcalled = nrhoma + nrhets + nrhomb
        nrdosages = nrcalled + nrmiss
        maf = 0
        if nrcalled > 0:
            maf = (nrhoma * 2 + nrhets) / (nrcalled * 2)
            if maf > 0.5:
                maf = 1 - maf

        return nrcalled, nrdosages, nrhoma, nrhets, nrhomb, maf

    @staticmethod
    def calculate_hwe(obs_hets, obs_hom1, obs_hom2):
        obs_homc = obs_hom1
        if obs_hom1 < obs_hom2:
            obs_homc = obs_hom2
        obs_homr = obs_hom2
        if obs_hom1 < obs_hom2:
            obs_homr = obs_hom1

        rare_copies = 2 * obs_homr + obs_hets
        l_genotypes = obs_hets + obs_homc + obs_homr

        if l_genotypes == 0:
            return -1

        het_probs = [0] * (rare_copies + 1)
        mid = int(rare_copies) * int(2 * l_genotypes - rare_copies) / int(2 * l_genotypes)
        mid = int(mid)
        if mid % 2 != rare_copies % 2:
            mid += 1
        mid = int(mid)

        curr_hets = mid
        curr_homr = (rare_copies - mid) / 2
        curr_homc = l_genotypes - curr_hets - curr_homr
        het_probs[int(mid)] = 1.0
        sum = het_probs[int(mid)]

        curr_hets = int(mid)
        while curr_hets > 1:
            het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
            sum += het_probs[curr_hets - 2]
            curr_homr += 1
            curr_homc += 1
            curr_hets -= 2

        curr_hets = int(mid)
        curr_homr = (rare_copies - mid) / 2
        curr_homc = l_genotypes - curr_hets - curr_homr
        while curr_hets <= (rare_copies - 2):
            het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
            sum += het_probs[curr_hets + 2]
            curr_homr -= 1
            curr_homc -= 1
            curr_hets += 2

        i = 0
        while i <= rare_copies:
            het_probs[i] /= sum
            i += 1

        p_hwe = 0.0
        i = 0
        while i <= rare_copies:
            if het_probs[i] <= het_probs[obs_hets]:
                p_hwe += het_probs[i]
            i += 1

        if p_hwe > 1:
            p_hwe = 1
        return p_hwe

    def print_arguments(self):
        print("Arguments:")
        print("  > VCF input: {}".format(self.input_path))
        print("  > Output filename: {}".format(self.output_filename))
        print("  > Samples: {}".format("" if self.sample_filter_set is None else "N={:,}\t{:,}".format(len(self.sample_filter_set), self.sample_filter_set)))
        print("  > Sex input: {}".format(self.sex_path))
        print("  > Thresholds:")
        print("    > Genotype quality (GQ): >={}".format(self.thresh_gq))
        print("    > Allelic balance (AB) lower-bound: <{}".format(self.thresh_ab_lower))
        print("    > Allelic balance (AB) upper-bound: >{}".format(self.thresh_ab_upper))
        print("    > Inbreeding coefficient (IB): >={}".format(self.thresh_ib))
        print("    > Minor allele frequency (MAF): >={}".format(self.thresh_maf))
        print("    > Call rate (CR): >={}".format(self.thresh_cr))
        print("    > Hardy-Weinberg equilibrium (HWE): >={}".format(self.thresh_hwe))
        print("    > Filtered depth (DP): >={}".format(self.thresh_dp))
        print("  > Check VQSR")
        print("    > SNV: check {}, <{}".format(self.check_vqsr_snv, self.thresh_vqsr_snv))
        print("    > Indel: check {}, <{}".format(self.check_vqsr_indel, self.thresh_vqsr_indel))
        print("  > Remove multiallelic: {}".format(self.remove_multiallelic))
        print("  > Remove non-PASS SNV: {}".format(self.remove_non_pass_snv))
        print("  > Remove non-PASS Indel: {}".format(self.remove_non_pass_indel))
        print("  > Remove low complexity: {}".format(self.remove_low_complexity))
        print("  > Strip INFO column: {}".format(self.strip_info_col))
        print("  > Ignore homozygous reference stats: {}".format(self.ignore_homref_stats))
        print("  > Replace poor quality genotype calls with missing: {}".format(self.replace_poor_quality_genotypes))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()