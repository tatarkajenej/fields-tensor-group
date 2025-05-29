SAMPLE = 50
FIELD = QQ
FIELD_RANDOM_PARAMS = (20, 10) # num_bound, den_bound

def RND(r, E, random_point_factory=None, points_used=SAMPLE):
    """
        Compute RND_r(E), where r is an integer and E is given as a matrix of linear forms.

        Return it as an abstract submodule object, its elements can be read of as matrices by using .lift() on them.

        random_point_factory should be a function that takes no arguments and returns a point on the rank r locus of E
    """

    h = E.nrows()
    w = E.ncols()
    polyring = E.parent().base_ring()
    x = polyring.gens()
    n = len(x)

    if random_point_factory is None:   # if factory is not given
        # try to take it using a hypersurface
        I = polyring.ideal(E.minors(r+1)).radical()
        # this might very much fail
        if len(I.gens()) > 1:
            raise ValueError("Could create a linear-variable-hypersurface-style factory.")
        f = I.gens()[0]
        random_point_factory = hypersurface_linear_variable_factory(f, E)

    def evalE(*vals):
        """
            Get an element of E with the given coordinates w.r.t. the parametrization by matrix of linear forms
        """
        return E.substitute({x[i] : vals[i] for i in range(n)}).change_ring(FIELD)

    matrix_space = MatrixSpace(FIELD, h, w)
    E_gens = []
    for i in range(n):
        vals = [0] * n
        vals[i] = 1
        E_gens.append(evalE(*vals))
    E_as_submodule = matrix_space.submodule(E_gens)

    def tang_space(e):
        """
            Return the tangent space of sigma_r Segre at e.

            (Presumes e has rank r, will probably throw an error otherwise.)
        """

        s, u, v = e.smith_form()
        uinv = u.inverse()
        vinv = v.inverse()

        # sanity checks
        assert(s == u * e * v)
        for i in range(min(h,w)):
            if i < r:
                assert(s[i,i] != 0)
            else:
                assert(s[i,i] == 0)

        generators = []
        for i in range(r):
            newgen = matrix_space()
            newgen[i,i] = 1
            generators.append(newgen)
            for j in range(i+1,h):
                newgen = matrix_space()
                newgen[j,i] = 1
                generators.append(newgen)
            for j in range(i+1,w):
                newgen = matrix_space()
                newgen[i,j] = 1
                generators.append(newgen)

        return matrix_space.submodule([uinv * gen * vinv for gen in generators])

    assert(SAMPLE > 0)

    result = None
    count = 0
    while count < SAMPLE:
        print("--cutting with sample no.", count)
        e = random_point_factory()
        assert(matrix_space.submodule([e]).is_submodule(E_as_submodule))    # check that the factory didn't cheat and really gave a point in E
        assert(e.rank() <= r)   # ... and that it gave a point of correct rank
        if e.rank() < r:    # try again if its not the rank wanted
            continue

        newcut = E_as_submodule + tang_space(e)
        if result is None:
            result = newcut
        else:
            result = result.intersection(newcut)

        if result is not None and result.dimension() <= E_as_submodule.dimension():
            break
        count += 1

    return result

def rnd_lift_attempt(rnd):
    """
        Take an RND as the new space of matrices E. Beware that this is probably not reasonable unless dim(RND) = dim(E) + 1.
    """
    return Tensor([e.lift() for e in rnd.gens()])

ATTEMPTS = 5
def hypersurface_linear_variable_factory(f, E):
    """
        Create a random_point_factory for a hypersurface cutting the space of matrices E.

        f needs to be from the polynomial ring that entries of E are in. It also needs to be linear with respect to at least one of its variables.
    """

    # identify the linear variable
    l = f.degrees().index(1)    # will throw ValueError if there is no 1 and that is correct
    # print("l is", l)

    polyring = E.parent().base_ring()
    x = polyring.gens()
    n = len(x)

    def factory():
        count = 0
        while count < ATTEMPTS:
            subst = {x[i] : FIELD.random_element(*FIELD_RANDOM_PARAMS) for i in range(n) if i != l}
            fbar = f.substitute(subst)
            if fbar.degree() == 1:
                c = fbar.coefficients()
                # print("just before final subst")
                subst[x[l]] = -c[1] / c[0] # stupid indexing, but it should be absolute / linear coef
                # print("just before in-factory return")
                return E.substitute(subst).change_ring(FIELD)
                break
            count += 1

    return factory

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

