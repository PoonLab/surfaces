import re
import os

def parse_ref_header(ref):
    h, original_rgene = get_ref(ref)
    header_lst = h.split(",")
    ref_loc = header_lst[-1].strip('"')
    ref_loc = ref_loc.split(":")
    start = ref_loc[0]
    stop = ref_loc[1]
    return(start, stop)

def get_ref(ref):
    ref_handle = open(ref)
    for line in ref_handle:
        if line.startswith('>'):
            h = line.strip('>#\n')
        else:
            seq = line.strip('\n').upper()
    return(h, seq)

def parse_csv(file):
    result = []
    handle = open(file)
    for line in handle:
        line = line.strip('\n')
        line = line.split(",")
        result.append(line)
    return(result)

def parse_loc(item):
    loc = item.replace("(+)","").replace('[','').replace(']','')
    cord = loc.split(":")
    start = cord[0]
    stop = cord[1]
    return(start, stop)
  
def main():
    csv_path = "/home/sareh/2020/polyprotein/poly_loc.csv"
    ref_dir = "/home/sareh/2020/polyprotein/ref_cds"

    lines = parse_csv(csv_path)
    for line in lines:
        start,stop = parse_loc(line[2])
        print(line)
        print("start: " + start)
        print("stop: " + stop)

        accn_full = line[0]
        prot_id = line[1]
        accn = re.sub("\.[1-9]","",accn_full)

        ref_path = os.path.join(ref_dir,accn)
        for file in os.listdir(ref_path):
            file_path = os.path.join(ref_path,file)

            ref_start, ref_stop = parse_ref_header(file_path)
            print("ref_start: " + ref_start)
            print("ref_stop: " + ref_stop)

            h, seq = get_ref(file_path)

            cut_start = int(start) - int(ref_start) - 1
            cut_stop = int(stop) - int(ref_start)

            print(cut_start)
            print(cut_stop)

            cut = seq[cut_start:cut_stop]
            print(cut)

            new_h = "{},{},name,num,{}:{}".format(accn_full,prot_id,start,stop)

            outdir = "/home/sareh/2020/polyprotein/ref_cds_cut/" 
            outname = "{}_{}".format(accn_full, prot_id)
            outpath = os.path.join(outdir, outname)
            out_handle = open(outpath,"w")

            out_handle.write(">{}\n{}\n".format(new_h, cut))

            print(outpath)

            print(new_h)


if __name__ == "__main__":
    main()  
  
