# input genbabk gen_records
# find all proteins for a virus (cds and mat_peptides)
# needs to handle cds, polyprotiens and segmented viruses

from Bio import SeqIO
import os

ref_dir = " "

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref', type=str,
                       help='ref directory')
    parser.add_argument('--loc', type=str,
                       help='csv with location of polyprotiens')
    parser.add_argument('--csv', type=str,
                       help='csv identifying segmented and polyproteins')
    parser.add_argument('--dir', type=str,
                       help='directory of genbank record files')
    return parser.parse_args()


def get_dict(file_path, csv, loc):

    record = SeqIO.read(file_path, "gb")

    d = {}
    d["id"] = file
    d["poly"] = "FALSE"
    d["seg"] = "FALSE"

    d["protein_id"] = []
    d["CDS"] = []

    for feature in record.features:
        # check if segmented
        if feature.type == "source" and 'segment' in feature.qualifiers:
            line.append("seg")
            d["seg"] = "TRUE"

        # check is it's a polyprotein
        if feature.type == 'mat_peptide':
            d["poly"] = "TRUE"

            try:
                loc = feature.location
                prot_accn = feature.qualifiers['protein_id'][0]

                #dict
                d["protein_id"].append((prot_accn,loc))

                #csv file
                line = "{},{},{}\n".format(record.id, prot_accn, loc)
                loc_handle = open(args.loc, 'a')
                loc_handle.write(line)

                #write out the dictionary
                l = ("{},{},{}\n".format(d["id"], d["seg"], d["poly"]))
                csv_handle = open(args.csv, "a")
                csv_handle.write(l)

           except:
                # skip viruses with missing record ID
                print("issue: " + record.id)
                
        if feature.type == 'CDS':
            try:
                prot_accn = feature.qualifiers['protein_id'][0]

            except:
                d[CDS].append(prot_accn)
    return(d)


def main()

    args = parse_args()

    for file in os.listdir(dir):
        print(file)
        file_path = os.path.join(dir,file)

        dict = get_dict(file_path, csv, loc)

        print(dict)
        
        ref_cds_dir_path = os.path.join(ref_dir, file)
        os.mkdir(ref_cds_dir_path)
       
if __name__ == '__main__':
    main()       
