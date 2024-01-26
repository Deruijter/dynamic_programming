# Original author: Anton Feenstra
# Stripped from original script by: Markus de Ruijter

# Built-in exchange matrices used to determine amino acid similarity
# Matrices are ordered alphabetically
class ExchangeMatrices:

    identity = [
        [ 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,],
        [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,],
    ]

    pam250 = [
        [ 2, 0,-2, 0, 0,-4, 1,-1,-1, 0,-1,-2,-1, 0, 0, 1, 0,-2, 1, 1, 0, 0,-6, 0,-3, 0,],
        [ 0, 2,-4, 3, 2,-5, 0, 1,-2, 0, 1,-3,-2, 2, 0,-1, 1,-1, 0, 0, 0,-2,-5, 0,-3, 2,],
        [-2,-4,12,-5,-5,-4,-3,-3,-2, 0,-5,-6,-5,-4, 0,-3,-5,-4, 0,-2, 0,-2,-8, 0, 0,-5,],
        [ 0, 3,-5, 4, 3,-6, 1, 1,-2, 0, 0,-4,-3, 2, 0,-1, 2,-1, 0, 0, 0,-2,-7, 0,-4, 3,],
        [ 0, 2,-5, 3, 4,-5, 0, 1,-2, 0, 0,-3,-2, 1, 0,-1, 2,-1, 0, 0, 0,-2,-7, 0,-4, 3,],
        [-4,-5,-4,-6,-5, 9,-5,-2, 1, 0,-5, 2, 0,-4, 0,-5,-5,-4,-3,-3, 0,-1, 0, 0, 7,-5,],
        [ 1, 0,-3, 1, 0,-5, 5,-2,-3, 0,-2,-4,-3, 0, 0,-1,-1,-3, 1, 0, 0,-1,-7, 0,-5,-1,],
        [-1, 1,-3, 1, 1,-2,-2, 6,-2, 0, 0,-2,-2, 2, 0, 0, 3, 2,-1,-1, 0,-2,-3, 0, 0, 2,],
        [-1,-2,-2,-2,-2, 1,-3,-2, 5, 0,-2, 2, 2,-2, 0,-2,-2,-2,-1, 0, 0, 4,-5, 0,-1,-2,],
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
        [-1, 1,-5, 0, 0,-5,-2, 0,-2, 0, 5,-3, 0, 1, 0,-1, 1, 3, 0, 0, 0,-2,-3, 0,-4, 0,],
        [-2,-3,-6,-4,-3, 2,-4,-2, 2, 0,-3, 6, 4,-3, 0,-3,-2,-3,-3,-2, 0, 2,-2, 0,-1,-3,],
        [-1,-2,-5,-3,-2, 0,-3,-2, 2, 0, 0, 4, 6,-2, 0,-2,-1, 0,-2,-1, 0, 2,-4, 0,-2,-2,],
        [ 0, 2,-4, 2, 1,-4, 0, 2,-2, 0, 1,-3,-2, 2, 0,-1, 1, 0, 1, 0, 0,-2,-4, 0,-2, 1,],
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
        [ 1,-1,-3,-1,-1,-5,-1, 0,-2, 0,-1,-3,-2,-1, 0, 6, 0, 0, 1, 0, 0,-1,-6, 0,-5, 0,],
        [ 0, 1,-5, 2, 2,-5,-1, 3,-2, 0, 1,-2,-1, 1, 0, 0, 4, 1,-1,-1, 0,-2,-5, 0,-4, 3,],
        [-2,-1,-4,-1,-1,-4,-3, 2,-2, 0, 3,-3, 0, 0, 0, 0, 1, 6, 0,-1, 0,-2, 2, 0,-4, 0,],
        [ 1, 0, 0, 0, 0,-3, 1,-1,-1, 0, 0,-3,-2, 1, 0, 1,-1, 0, 2, 1, 0,-1,-2, 0,-3, 0,],
        [ 1, 0,-2, 0, 0,-3, 0,-1, 0, 0, 0,-2,-1, 0, 0, 0,-1,-1, 1, 3, 0, 0,-5, 0,-3,-1,],
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
        [ 0,-2,-2,-2,-2,-1,-1,-2, 4, 0,-2, 2, 2,-2, 0,-1,-2,-2,-1, 0, 0, 4,-6, 0,-2,-2,],
        [-6,-5,-8,-7,-7, 0,-7,-3,-5, 0,-3,-2,-4,-4, 0,-6,-5, 2,-2,-5, 0,-6,17, 0, 0,-6,],
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
        [-3,-3, 0,-4,-4, 7,-5, 0,-1, 0,-4,-1,-2,-2, 0,-5,-4,-4,-3,-3, 0,-2, 0, 0,10,-4,],
        [ 0, 2,-5, 3, 3,-5,-1, 2,-2, 0, 0,-3,-2, 1, 0, 0, 3, 0, 0,-1, 0,-2,-6, 0,-4, 3,],
    ]

    blosum62 = [
        [ 4,-2, 0,-2,-1,-2, 0,-2,-1, 0,-1,-1,-1,-2, 0,-1,-1,-1, 1, 0, 0, 0,-3, 0,-2,-1,],
        [-2, 4,-3, 4, 1,-3,-1, 0,-3, 0, 0,-4,-3, 3, 0,-2, 0,-1, 0,-1, 0,-3,-4,-1,-3, 1,],
        [ 0,-3, 9,-3,-4,-2,-3,-3,-1, 0,-3,-1,-1,-3, 0,-3,-3,-3,-1,-1, 0,-1,-2,-2,-2,-3,],
        [-2, 4,-3, 6, 2,-3,-1,-1,-3, 0,-1,-4,-3, 1, 0,-1, 0,-2, 0,-1, 0,-3,-4,-1,-3, 1,],
        [-1, 1,-4, 2, 5,-3,-2, 0,-3, 0, 1,-3,-2, 0, 0,-1, 2, 0, 0,-1, 0,-2,-3,-1,-2, 4,],
        [-2,-3,-2,-3,-3, 6,-3,-1, 0, 0,-3, 0, 0,-3, 0,-4,-3,-3,-2,-2, 0,-1, 1,-1, 3,-3,],
        [ 0,-1,-3,-1,-2,-3, 6,-2,-4, 0,-2,-4,-3, 0, 0,-2,-2,-2, 0,-2, 0,-3,-2,-1,-3,-2,],
        [-2, 0,-3,-1, 0,-1,-2, 8,-3, 0,-1,-3,-2, 1, 0,-2, 0, 0,-1,-2, 0,-3,-2,-1, 2, 0,],
        [-1,-3,-1,-3,-3, 0,-4,-3, 4, 0,-3, 2, 1,-3, 0,-3,-3,-3,-2,-1, 0, 3,-3,-1,-1,-3,],
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
        [-1, 0,-3,-1, 1,-3,-2,-1,-3, 0, 5,-2,-1, 0, 0,-1, 1, 2, 0,-1, 0,-2,-3,-1,-2, 1,],
        [-1,-4,-1,-4,-3, 0,-4,-3, 2, 0,-2, 4, 2,-3, 0,-3,-2,-2,-2,-1, 0, 1,-2,-1,-1,-3,],
        [-1,-3,-1,-3,-2, 0,-3,-2, 1, 0,-1, 2, 5,-2, 0,-2, 0,-1,-1,-1, 0, 1,-1,-1,-1,-1,],
        [-2, 3,-3, 1, 0,-3, 0, 1,-3, 0, 0,-3,-2, 6, 0,-2, 0, 0, 1, 0, 0,-3,-4,-1,-2, 0,],
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
        [-1,-2,-3,-1,-1,-4,-2,-2,-3, 0,-1,-3,-2,-2, 0, 7,-1,-2,-1,-1, 0,-2,-4,-2,-3,-1,],
        [-1, 0,-3, 0, 2,-3,-2, 0,-3, 0, 1,-2, 0, 0, 0,-1, 5, 1, 0,-1, 0,-2,-2,-1,-1, 3,],
        [-1,-1,-3,-2, 0,-3,-2, 0,-3, 0, 2,-2,-1, 0, 0,-2, 1, 5,-1,-1, 0,-3,-3,-1,-2, 0,],
        [ 1, 0,-1, 0, 0,-2, 0,-1,-2, 0, 0,-2,-1, 1, 0,-1, 0,-1, 4, 1, 0,-2,-3, 0,-2, 0,],
        [ 0,-1,-1,-1,-1,-2,-2,-2,-1, 0,-1,-1,-1, 0, 0,-1,-1,-1, 1, 5, 0, 0,-2, 0,-2,-1,],
        [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
        [ 0,-3,-1,-3,-2,-1,-3,-3, 3, 0,-2, 1, 1,-3, 0,-2,-2,-3,-2, 0, 0, 4,-3,-1,-1,-2,],
        [-3,-4,-2,-4,-3, 1,-2,-2,-3, 0,-3,-2,-1,-4, 0,-4,-2,-3,-3,-2, 0,-3,11,-2, 2,-3,],
        [ 0,-1,-2,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1,-1, 0,-2,-1,-1, 0, 0, 0,-1,-2,-1,-1,-1,],
        [-2,-3,-2,-3,-2, 3,-3, 2,-1, 0,-2,-1,-1,-2, 0,-3,-1,-2,-2,-2, 0,-1, 2,-1, 7,-2,],
        [-1, 1,-3, 1, 4,-3,-2, 0,-3, 0, 1,-3,-1, 0, 0,-1, 3, 0, 0,-1, 0,-2,-3,-1,-2, 4,],
    ]