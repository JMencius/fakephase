from cyvcf2 import VCF, Writer


def output_vcf(invcf: str, outvcf: str, fake_blocks: list, chrom: list) -> None:
    vcf_in = VCF(invcf)
    
    vcf_in.add_format_to_header({
    "ID": "PS",
    "Number": "1",
    "Type": "Integer",
    "Description": "Phase Set identifier"
    })
    
    vcf_out = Writer(outvcf, vcf_in)
    
    for var in vcf_in:
        # first set it to unphased
        var.format["PS"] = '.'
        if var.CHROM in chrom:
            working_index = chrom.index(var.CHROM)
            actions = fake_blocks[working_index]
            if var.POS in actions:
                var.format["PS"] = actions[var.POS][1]
            
                # flip or keep
                left, right, _ = variant.genotypes[0]
                if actions[var.POS][0] == -1:
                    # flip and set to phased
                    variant.genotype[0] = right, left, 1
                else:
                    # keep and set to phased
                    variant.genotype[0] = left, right, 1
        
        # write to new vcf
        vcf_out.write_record(var)
                
        
