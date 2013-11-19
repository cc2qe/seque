# awk script to split a vcf into SNVs and indels
# Lifted from: http://www.biostars.org/p/7403/ 2013-11-19
# usage: awk -f split_vcf_snv_indel.awk file.vcf

/^#/    {
    print $0 > "snv.vcf";
    print $0 > "indels.vcf";
    next;
}

/^[^\t]+\t[0-9]+\t[^\t]*\t[atgcATGC]\t[a-zA-Z]\t/   {
    print $0 > "snv.vcf";
    next;
}

{
    print $0 > "indels.vcf";
    next;
}