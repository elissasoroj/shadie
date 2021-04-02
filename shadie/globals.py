#standard mutation types, with structure:
#(name, dominance, distribution, {following depends on distribution})

from shadie.mutations import MutationType
from shadie.elements import ElementType

#saved as MutationType and ElementType objects
NEUT = MutationType(0.5, "f", 0.0)		#neutral mutation
SYN = MutationType(0.5, "f", 0.0)         #synonymous
DEL = MutationType(0.1, "g", -0.03, 0.2)  #deleterious
BEN = MutationType(0.8, "e", 0.1)          	#beneficial

EXON = ElementType([SYN, DEL, BEN], (2,8,0.1))  #exon
INTRON = ElementType([SYN,DEL], (9,1))         #intron
NONCOD = ElementType(NEUT, 1)              	 #non-coding

#nucleotide defulats

nucEXON = ElementType([SYN, DEL, BEN], (2,8,0.1))

if __name__ == "__main__":

    # generate random chromosome
    print(EXON)