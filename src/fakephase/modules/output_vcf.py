import numpy as np
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
        ##var.format["PS"] = '.'
        var.set_format("PS", np.array(["."], dtype = 'S'))
        if var.CHROM in chrom:
            working_index = chrom.index(var.CHROM)
            actions = fake_blocks[working_index]
            if var.POS in actions:
                #var.format["PS"] = actions[var.POS][1]
                var.set_format("PS", np.array([actions[var.POS][1]], dtype = np.int32))

                # flip or keep
                left, right, _ = var.genotypes[0]
                if actions[var.POS][0] == -1:
                    # flip and set to phased
                    new_gt = np.array([[right, left, True]], dtype=np.int32)
                    var.genotypes = new_gt
                else:
                    # keep and set to phased
                    new_gt = np.array([[left, right, True]], dtype=np.int32)
                    var.genotypes = new_gt
        
        # write to new vcf
        vcf_out.write_record(var)
                
        
