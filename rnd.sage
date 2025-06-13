SAMPLE = 50
FIELD = QQ
FIELD_RANDOM_PARAMS = (20, 10) # num_bound, den_bound
ATTEMPTS = 5

load("factories.sage")

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
        random_point_factory = my_best_factory(r, E)

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
    bad_attempts = 0
    while count < SAMPLE:
        print("--cutting with sample no.", count)
        e = random_point_factory()
        assert(matrix_space.submodule([e]).is_submodule(E_as_submodule))    # check that the factory didn't cheat and really gave a point in E
        assert(e.rank() <= r)   # ... and that it gave a point of correct rank
        if e.rank() < r:    # try again if its not the rank wanted
            bad_attempts += 1
            if bad_attempts > ATTEMPTS:
                raise ValueError("Too many sub-rank results from the factory.")
            continue

        newcut = E_as_submodule + tang_space(e)
        if result is None:
            result = newcut
        else:
            result = result.intersection(newcut)

        if result is not None and result.dimension() <= E_as_submodule.dimension():
            break
        bad_attempts = 0
        count += 1

    return result

load("tensor.sage")

def rnd_lift_attempt(rnd):
    """
        Take an RND as the new space of matrices E. Beware that this is probably not reasonable unless dim(RND) = dim(E) + 1.
    """
    return Tensor([e.lift() for e in rnd.gens()])
