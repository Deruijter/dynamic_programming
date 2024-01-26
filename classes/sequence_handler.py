# This class handles Sequence data needed for the
# Dynamic Programming assignment for the "Algorithms in Sequence Alignment" course

from entities.sequence import Sequence

class SequenceHandler:

	# This routine parses a fasta file and returns a list of Sequence objects
	# containing the sequences with label and sequence data set
	# Return type: list<Sequence>
	def read_sequences(self,filename):
	    # read sequences from fasta file, and catch error reading file
	    try:
	        lines = open(filename)
	    except IOError:
	        print("ERROR: cannot open or read fasta input file:", fastafile)
	        exit(-1)

	    seqs = []
	    label = None
	    seq_lines = []
	    for line in lines:
	        line = line.strip()      # strip off white space
	        if not line:             # skip empty lines
	            continue
	        if line.startswith(';'): # ignore comment lines
	            continue
	        # check for start of next sequence:
	        if line.startswith('>'): # label line
	            # first, store the previous sequence if we had one:
	            if seq_lines:
	                seqs.append(Sequence(label, ''.join(seq_lines)))
	                seq_lines = []
	            # get the label (name) for the next sequence
	            label = line[1:].strip()
	        else:
	            # collect all lines with sequence information for this sequence:
	            seq_lines.append(line)
	    # take care of the last sequence in the file
	    seqs.append(Sequence(label, ''.join(seq_lines)))

	    for seq in seqs:
	        print(seq)
	    print("")   

	    return(seqs)
