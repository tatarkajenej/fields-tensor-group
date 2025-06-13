"""
Auxiliary file to handle "random point factories" -- functions that randomly generate rational points on varieties.

This is of course hard to handle generally (the variety might not even have rational points), so I just construct some hacked-together solutions for the types of varieties that come up in examples.
"""

FIELD = QQ
FIELD_RANDOM_PARAMS = (20, 10) # num_bound, den_bound
FACTORY_ATTEMPTS = 5

def solve_linear_equation(f, var):
    """
        Solve an equation given by a linear polynomial. The parent ring might have more variables, so we specify it explicitly
    """

    return -f.coefficient({var:0}) / f.coefficient({var:1})

def hypersurface_linear_variable_factory(f, E, return_matrix=True):
    """
        Create a random_point_factory for a hypersurface cutting the space of matrices E.

        f needs to be from the polynomial ring that entries of E are in. It also needs to be linear with respect to at least one of its variables.

        return_matrix: if True, return E.subst(...), if False, return the substitution dictionary
    """

    # identify the linear variable
    l = f.degrees().index(1)    # will throw ValueError if there is no 1 and that is correct

    polyring = E.parent().base_ring()
    x = polyring.gens()
    n = len(x)

    def factory():
        count = 0
        while count < FACTORY_ATTEMPTS:
            subst = {x[i] : FIELD.random_element(*FIELD_RANDOM_PARAMS) for i in range(n) if i != l}
            fbar = f.substitute(subst)
            if fbar.degree() == 1:
                # c = fbar.coefficients()
                subst[x[l]] = solve_linear_equation(fbar, x[l]) # stupid indexing, but it should be absolute / linear coef
                if return_matrix:
                    return E.substitute(subst).change_ring(FIELD)
                else:
                    return subst
                break
            count += 1
        else:
            raise ValueError("Failed to generate a point in the allowed number of attempts.")

    return factory

def trivial_factory(E, return_matrix=True):
    """
        Just take random points in E.
    """

    def factory():
        subst = {E.parent().base_ring().gens()[i] : FIELD.random_element(*FIELD_RANDOM_PARAMS) for i in range(E.parent().base_ring().ngens())}
        if return_matrix:
            return E.substitute(subst)
        else:
            return subst

    return factory

def my_best_factory(r, E, return_matrix=True):
    """
        Attempt at a parent function that can handle as many cases as possible
    """
    # try to take it using a hypersurface
    polyring = E.parent().base_ring()
    x = polyring.gens()
    n = len(x)
    I = polyring.ideal(E.minors(r+1)).radical()
    if I.is_zero():
        return trivial_factory(E, return_matrix)
    elif len(I.gens()) == 1:
        f = I.gens()[0]
        # this might fail if there's no linear variable
        return hypersurface_linear_variable_factory(f, E, return_matrix)
    else:
        variable_gens = []
        other_gens = []
        for f in I.gens():
            if f.is_monomial() and f.degree() == 1:
                variable_gens.append(f)
            else:
                other_gens.append(f)
        if len(other_gens) <= 1:
            if len(other_gens) == 0:
                subfactory = trivial_factory(E, False)
            elif len(other_gens) == 1:
                subfactory = hypersurface_linear_variable_factory(other_gens[0], E, False)
            def new_factory():
                subst = subfactory()
                for m in variable_gens:
                    subst[m] = 0
                if return_matrix:
                    return E.substitute(subst)
                else:
                    return subst
            return new_factory
        else:
            raise ValueError("Couldn't create a factory.")
