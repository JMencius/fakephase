

class myvariant:
    def __init__(self, chrom, pos, ref, alt, left, right, category):
        self.chrom : str = chrom
        self.pos : int = pos
        self.ref : str = ref
        self.alt : list = alt
        self.left : str = left
        self.right : str = right
        self.category : str = category


    def __str__(self) -> str:
        return (f"CHR:{self.chrom} POS:{self.pos} REF:{self.ref} ALT:{self.alt} LEFT: {self.left} RIGHT: {self.right} category:{self.category}")

