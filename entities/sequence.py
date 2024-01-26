class Sequence:
    """Stores a sequence object"""
    
    def __init__(self, Label="", Sequence="" ):
        """Initialize a new Sequence object

        Label -- identifier of sequence (text)
        Sequence -- sequence string in single-letter alphabet
        """
        self.Label       = Label
        self.Sequence    = Sequence

    # this makes that you can do 'print sequence' and get nice output:
    def __str__(self):
        """Return string representation of a Sequence object"""
        # newline-delimited values of all the attributes
        return(">%s\n%s" % (self.Label, self.Sequence))
