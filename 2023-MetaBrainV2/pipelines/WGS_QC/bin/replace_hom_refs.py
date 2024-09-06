#!/opt/conda/envs/eQTLGenPopAssign/bin/python


"""
    usage:
        ./dotrevive.py -i [INPUT VCF FILE] -o [OUTPUT VCF FILE]
"""

# METADATA VARIABLES
__author__ = "Orfeas Gkourlias"
__status__ = "WIP"
__version__ = "0.1"

# IMPORTS
import sys
import argparse
import gzip

# CLASSES
class DotRevive:
    """Class to handle VCF file and retrieve variant & dosgae statistics."""
    def __init__(self, input_handle, output_handle):
        self.input_handle = input_handle
        self.output_handle = output_handle
        
    def walk(self):
        """Walk through the VCF file and retrieve the variants."""
        with gzip.open(self.input_handle, "rt") as vcf_i, gzip.open(self.output_handle, "wt") as vcf_o:
            ctr = 0
            for line in vcf_i:
                if line.startswith("#") == False:
                    ctr += 1
                    vcf_o.write(self.record_extract(line))
                    if ctr % 10000 == 0:
                        print(f"{ctr} Variants Written.", end = "\r")
                        ctr += 1
                else: 
                    vcf_o.write(line)
                    if ctr % 10000 == 0:
                        print(f"{ctr} Variants Written.", end = "\r")
                        ctr += 1

            print(f"A Total of {ctr} Variants Were Written. Done!")
                

    def record_extract(self, variant):
        """Extracts the relevant intersection information from two tabix records.

        Args:
            i_record (list): list containing record entries for input file.
            headdx (dict): header indexes

        Returns:
            string: string to be written to output file.
        """
        # Get a column value for each row if header/column index is present in both files.
        variant = variant.strip()
        variant = variant.split("\t")
        new = []

        for i, entry in enumerate(variant):
            # match i:  # Match the index value.
            
            if i in [0,1,2,3,4,5,6,7]:  # Values on these indexes are not relevant.
                new.append(entry)

            elif i == 8:  # Initialse the format column indicators/headers.
                val_dict = {}

                # Store format indication/header into a dictionary with empty values.
                for ind in entry.split(":"):
                    val_dict[ind] = "."
                
                new.append(entry)


            else:  # Sample rows.
                # Assign values in format order.
                for key, val in zip(
                    val_dict, entry.split(":")
                ):
                    val_dict[key] = val
            
                # Save genotypes in seperate var.
                gt = val_dict["GT"]
                ad = val_dict["AD"]
                sep = gt[1]

                if gt  == f"0{sep}0" and ad == ".":
                    val_dict["GT"] = "./."
                
                new_vals = ":".join(str(dict_val) for dict_val in val_dict.values())
                new.append(new_vals)
    
    # Return a constructed row string for immediate output writing.
        return (
            "\t".join(new) + "\n"
        )
                    

# FUNCTIONS
def arg_parse():

    """
    Function to parse arguments.
    :file: Input VCF file.
    :reference: VCF file to be compared against.
    :input_header: Text file containing tab seperated headers.
    :comp_header: Text file containing tab seperated headers for compare file.
    :chr: Test file named as input chromosome.
    :return: Argparse namespace object.
    """

    # Parser init
    parser = argparse.ArgumentParser(prog="dotrevive")
    parser.add_argument("-i", "--input", type=str)
    parser.add_argument("-o", "--output", type=str)
    return parser.parse_args()


# MAIN
def main(args):
    """Main function"""
    # Get args.
    args = arg_parse()
    # Perform comparison and write output.
    dot_revive = DotRevive(args.input, args.output)
    dot_revive.walk()
    # FINISH
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))