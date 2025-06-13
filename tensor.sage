class Tensor(object):
    """
        Class to hold a 3-tensor (very much in bases) and give flattenings as matrices of linear forms
    """

    def __init__(self, data):
        if type(data) == list:  # presume it's a 3-nested list of FIELD entries
            if type(data[0]) == list and type(data[0][0]) == list:
                self.shape = (len(data), len(data[0]), len(data[0][0]))
                self.entries = data
            elif type(data[0]) in (sage.matrix.matrix_rational_dense.Matrix_rational_dense,):
                self.shape = (len(data), data[0].nrows(), data[0].ncols())
                self.entries = [[[data[k][i][j] for j in range(self.shape[2])] for i in range(self.shape[1])] for k in range(self.shape[0])]

        elif type(data) in (sage.matrix.matrix_mpolynomial_dense.Matrix_mpolynomial_dense,):    # it's a matrix of linear forms
            polyring = data.parent().base_ring()
            x = polyring.gens()
            n = len(x)
            self.shape = (n, data.nrows(), data.ncols())
            self.entries = []
            for i in range(n):
                subst = {x[j] : 0 for j in range(n)}
                subst[x[i]] = 1
                newgen = data.substitute(subst).change_ring(FIELD)
                self.entries.append([list(row) for row in newgen.rows()])

        else:
            raise ValueError("Unrecognized input format for Tensor.")

    def __getitem__(self, pos):
        """
            pos should be a container supporting indices 0, 1, 2 (possibly among others)
        """
        return self.entries[pos[0]][pos[1]][pos[2]]

    def flatten(self, A=0, B=1, C=2):
        """
            Return flattening as a matrix of linear forms.

            A, B, C should be a permutation of (1,2,3) -- A will index the variables, B the rows and C the columns
        """

        assert(set([A,B,C]) == set([0,1,2]))
        def perm(i,j,k):
            """
                Auxiliary function to translate between the given permutation
            """
            result = [0] * 3
            result[A] = i
            result[B] = j
            result[C] = k

            return result

        n = self.shape[A]
        h = self.shape[B]
        w = self.shape[C]

        polyring = PolynomialRing(FIELD, self.shape[A], "x")
        x = polyring.gens()
        return Matrix([[sum([self[perm(k,hh,ww)] * x[k] for k in range(n)]) for ww in range(w)] for hh in range(h)])
