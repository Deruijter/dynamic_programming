#!/usr/bin/python

"""
Dynamic programming excercises.
"""

# for commandline options:
from optparse import OptionParser, OptionGroup

from classes.sequence_handler import SequenceHandler
from classes.alignment_handler import AlignmentHandler

def parse_commandline():
    usage = "%prog <fasta> [options]"
    version = "v1.0.0.0"
    description = \
        "%prog aligns two sequences."
    epilog = \
        "Author: Markus de Ruijter                                                     "\
        "Original by Anton Feenstra -- feenstra@few.vu.nl -- www.few.vu.nl/~feenstra"
    parser = OptionParser(usage=usage, description=description,
                          version="%prog "+version, epilog=epilog)

    # sequence/alignment options:
    parser.add_option("-f", "--fasta",  dest="fasta", metavar="<file>",
                     help="input alignment file (fasta)")
    parser.set_defaults(fasta=None)
    
    parser.add_option("-e", "",  dest="exchange_matrix", 
                     help="Exchange matrix: pam250, blosum62 or identity (%default)")
    parser.set_defaults(exchange_matrix="pam250")
    
    parser.add_option("-l", "",  dest="align_local",  action="store_true",
                     help="align local")
    parser.set_defaults(align_local=False)
    
    parser.add_option("-g", "",  dest="align_global", action="store_true",
                     help="align global")
    parser.set_defaults(align_global=False)
    
    parser.add_option("-s", "",  dest="align_semiglobal", action="store_true",
                     help="align semi-global")
    parser.set_defaults(align_semiglobal=False)
    
    parser.add_option("-p", "",  dest="gap_penalty", type="int",
                     help="Gap penalty (%default)")
    parser.set_defaults(gap_penalty=2)
    
    parser.add_option("-v", "",  dest="verbose", action="store_true",
                     help="Print additional information")
    parser.set_defaults(verbose=False)
    
    #local alignment specific
    parser.add_option("", "--lalt",  dest="local_alignment_limiter_type", type="string",
                     help="Local Alignment Limiter Type (score|length|iterations) (%default)")
    parser.set_defaults(local_alignment_limiter_type="iterations")

    parser.add_option("", "--lalv",  dest="local_alignment_limiter_value", type="int",
                     help="The value to use for the Local Alignment Limiter Type (%default)")
    parser.set_defaults(local_alignment_limiter_value=5)
    
    # get the options:
    (options, args) = parser.parse_args()

    if not options.fasta:
        # check if we have an option left (to be used as input filename):
        if args:
            options.fasta = args.pop()
        else:
            print("Need at least an input file (fasta)")
            print("")
            parser.print_help()
            print("")
            print("ERROR: no input file given")
            exit(-1)

    # check alignment type:
    align_options = [options.align_local, options.align_global, options.align_semiglobal]
    # check if at least one alignment option was true, else choose global
    if align_options.count(True)==0:
        print("No alignment type given, using Global")
        options.align_global=True
    # check if not more than one alignment option was true, else error and exit 
    if align_options.count(True)>1:
        print("ERROR: multiple alignment types chosen")
        exit(-1)

    if options.local_alignment_limiter_value < 0:
        print("ERROR: Local alignment limiter value must be greater than 0")
        exit(-1)

    # check for any leftover command line arguments:
    if len(args):
        warning("ignoring additional arguments "+str(args))
    
    # clean up (recommended):
    del(parser)
    return options

# Get alignment type based on the supplied arguments
def get_alignment_type(options):
    if options.align_global:
        return("global")
    if options.align_semiglobal:
        return("semi-global")
    if options.align_local:
        return("local")

# main function:
def main():
    # get command line options
    options = parse_commandline()

    # Instantiate classes
    sequence_handler = SequenceHandler()
    alignment_handler = AlignmentHandler(get_alignment_type(options)
        , options.exchange_matrix
        , sequence_handler.read_sequences(options.fasta)
        , options.gap_penalty
        , options.verbose
        , options.local_alignment_limiter_type
        , options.local_alignment_limiter_value)

    alignment_handler.align()

if __name__ == "__main__":
    main()
