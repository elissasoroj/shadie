#standard mutation types, with structure:
#(name, dominance, distribution, {following depends on distribution})

from mutations import MutationType
from elements import ElementType

#old strings - depreciated
S_NEUT = '("m100", 0.5, "f", 0.0)'      #neutral mutation
S_SYN = '("m101", 0.5, "f", 0.0)'           #synonymous
S_DEL = '("m102", 0.1, "g", -0.03, 0.2)'    #deleterious
S_BEN = '("m103", 0.8, "e", 0.1)'           #beneficial

#standard genetic elements have the structure:
#(name, mutations, frequency of each mutation)
S_EXON = '("g100", c(m100,m101,m102), c(2,8,0.1))'  #exon
S_INTRON = '("g101", c(m100,m102), c(9,1))'         #intron
S_NONCOD = '("g102", c(m100), 1)'               #non-coding


#saved as MutationType and ElementType objects
NEUT = MutationType(0.5, "f", [0.0])		#neutral mutation
SYN = MutationType(0.5, "f", [0.0])         #synonymous
DEL = MutationType(0.1, "g", [-0.03, 0.2])  #deleterious
BEN = MutationType(0.8, "e", 0.1)          	#beneficial

EXON = ElementType([SYN, DEL, BEN], (2,8,0.1))  #exon
INTRON = ElementType([SYN,DEL], (9,1))         #intron
NONCOD = ElementType(NEUT, 1)              	 #non-coding