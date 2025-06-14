class RankProfile(object):
    """
        Object to hold information about various (rank <= r) loci of a space of matrices
    """

    def __init__(self, E, take_radicals=True):
        self.matrices = E
        self.polyring = E.base_ring()
        self.take_radicals = take_radicals

        self.has_loci = False

    def compute_loci(self):
        """
            Compute the (rank <= r) loci of this space of matrices. Only take those r that are actually achieved by some matrix in the space
        """

        self.loci = {}
        self.dims = {}
        self.ranks = []
        last_d = None

        for r in range(min(self.matrices.ncols(), self.matrices.nrows())+1):
            I = self.polyring.ideal(self.matrices.minors(r+1))
            if self.take_radicals:
                I = I.radical()
            d = I.dimension()
            if last_d is None or d > last_d:
                self.ranks.append(r)
                self.loci[r] = I
                self.dims[r] = d
                last_d = d

        self.has_loci = True

    def get_next(self, r):
        """
            Return the rank r' in self.ranks such that the (rank <= r) locus actually coincides with (rank <= r').

            This is just the largest r' <= r such that r' actually appears as a rank of some matrix in the space.
        """

        if r < 0:
            raise ValueError("It doesn't make sense to ask about negative rank.")

        result = None
        for rr in self.ranks:
            if rr <= r:
                result = rr
        return result

    def get_locus(self, r):
        """
            Return the ideal of the (rank <= r) locus
        """
        if not self.has_loci:
            self.compute_loci()
        return self.loci[self.get_next(r)]


    def get_dim(self, r):
        """
            Return the dimension of the (rank <= r) locus
        """
        if not self.has_loci:
            self.compute_loci()
        return self.dims[self.get_next(r)]

    def __getitem__(self, r):
        return self.get_locus(r)

    def geometric_rank(self):
        """
            Return the geometric rank of the given tensor (or rather a flattening of a tensor).
        """
        if not self.has_loci:
            self.compute_loci()

        GR = None
        n = self.polyring.ngens()
        for r in self.loci:
            d = self.loci[r].dimension()
            GR_Aj = n - d + r
            if GR is None or GR > GR_Aj:
                GR = GR_Aj
        return GR

    def __str__(self):
        if not self.has_loci:
            self.compute_loci()

        return "\n\n".join(["rank at most %i:\n\t%s\n\t(dim = %i)" % (r, extract_generator_string(self.loci[r]), self.dims[r]) for r in self.ranks[::-1]]) + "\n---------------------------\n" + "geometric rank: %i" % self.geometric_rank()

    def tex_string(self):
        """
            Return a nice LaTeX-ready string, for the purposes of automatically generating some factsheets for any tensors of interest
        """

        ### TODO
        raise NotImplementedError

def extract_generator_string(ideal):
    """
        Extract just the generators from str(ideal), throw away additional text and info about the parent ring.
    """
    s = str(ideal)
    begin = s.find("(")
    end = s.rfind(")")
    return s[begin+1:end]
