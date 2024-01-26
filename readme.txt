Purpose:
	The purpose of this application is to align two sequences against each other
	using global, semi-global and local alignment techniques
Usage:
	Store two sequences in a file, in fasta format. Example:
	>seq1
	YICSFADCCF
	>seq2
	FPCKEECA
	(note that only the first 2 sequences are handled at this time)

	Then run the python script and provide the filename as value for the -f argument

	For more information about the script use the -h parameter

Structure:
	The main function is located in align.py
	The Classes folder contains data processors that handle most of the work
	The Entities folder contains custom type objects
	
This was part of my Bioinformatics & Systems biology masters at the University of Amsterdam.
The excercise and initial file was created by Anton Feenstra.
