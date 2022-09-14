import subprocess
import biolib
from IPython.display import Image

fasta_path = argv[1]

deeptmhmm = biolib.load('DTU/DeepTMHMM')

input = "--fasta {}".format(fasta_path)

deeptmhmm_res = deeptmhmm.cli(args=input)

outfile = "{}_tmhmm/".format(fasta_path)
deeptmhmm_res.save_files(outfile)

plot = "{}.png".format(fasta_path)
plot_path = os.path.join(outfile, plot)
Image(filename=plot_path)
