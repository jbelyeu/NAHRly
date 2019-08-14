import cyvcf2
import numpy as np

header_tmpl = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="integer copy-number">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="integer genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Float,Description="normalized depth">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=CNV,Description="Copy-number Variant">{contigs}
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"""

variant_tmpl = """{chrom}\t{POS}\t.\tN\t<{svt}>\t{qual}\tPASS\tSVTYPE={svt};END={stop};SVLEN={svlen}\tGT:GQ:CN:DP"""

def get_writer(out_path, sample_names, chroms=None):
    if chroms is not None:
        contigs = "\n" + "\n".join("##contig=<ID=%s>" % c for c in chroms)
    else:
        contigs = ""
    header = header_tmpl + "\t" + "\t".join(sample_names)
    header = header.format(contigs=contigs)
    return cyvcf2.Writer.from_string(out_path, header)

def get_gt(svt, cn):
    if svt == "DEL":
        return ["1/1", "0/1", "0/0"][min(2, cn)]
    if svt == "DUP":
        return ["0/0", "0/1", "1/1"][max(0, min(2, cn - 2))]

def get_gq():
    # TODO:
    return 30

def generate_sample_fields(variant_dict, svt):
    result = []
    for i in range(len(variant_dict["CN"])):
        cn = variant_dict["CN"][i]
        gt = [get_gt(svt, cn)]
        gt.append(str(get_gq()))
        gt.append(str(cn))
        gt.append("%.3f" % variant_dict["DP"][i])
        result.append(":".join(gt))
    return result

def write_variant(wtr, variant_dict):

    variant_dict = variant_dict.copy()

    variant_dict["CN"] = np.asarray(variant_dict["CN"])
    variant_dict["DP"] = np.asarray(variant_dict["DP"])
    variant_dict["svlen"] = variant_dict["stop"] - variant_dict["start"]
    variant_dict["POS"] = variant_dict["start"] + 1

    less = variant_dict["CN"] <= 2
    more = variant_dict["CN"] > 2
    for (svt, sel) in (("DEL", less), ("DUP", more)):
        if svt == "DUP" and sel.sum() == 0: continue
        variant_dict["qual"] = 50 # TODO
        variant_dict["svt"] = svt

        variant = variant_tmpl.format(**variant_dict)
        variant += "\t" + "\t".join(generate_sample_fields(variant_dict, svt))

        vobj = wtr.variant_from_string(variant)
        wtr.write_record(vobj)

if __name__ == "__main__":

    w = get_writer("__test.vcf", ["sampleA", "sampleB"])

    write_variant(w, {"chrom": "chr1", "start": 12354, "stop": 45677, "CN": [0,
        1], "DP": [0.02, 0.55]})
    w.close()

