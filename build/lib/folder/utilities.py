llista_lletres= ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S",
               "T", "U", "V", "W", "X", "Y", "Z", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"
                , "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"]
llista_lletres2=[]

llista_lletres2= [element+element2 for element in llista_lletres for element2 in llista_lletres]
merged_list= llista_lletres + llista_lletres2
