load("rnd.sage")

R = PolynomialRing(FIELD, ["x", "y", "z", "w"])
x,y,z,w = R.gens()
Ematrix_multiplication_2x2 = Matrix([
    [x,y,0,0],
    [z,w,0,0],
    [0,0,x,y],
    [0,0,z,w]
])
matrix_mult = Tensor(Ematrix_multiplication_2x2)

R = PolynomialRing(FIELD, 6, "x")
x0, x1, x2, x3, x4, x5 = R.gens()
EDerek6 = Matrix([
    [ 2*x3,  2*x4,  2*x5,     0,     0,     0],
    [  -x1,   -x2,     0,    x4,  2*x5,     0],
    [    0,   -x1,   -x2, -2*x3,   -x4,     0],
    [ 2*x0,     0,     0, -2*x2,     0,  2*x5],
    [    0,  2*x0,     0,    x1,   -x2,   -x4],
    [    0,     0,  2*x0,     0,  2*x1,  2*x3]
])
Derek6 = Tensor(EDerek6)
# There is a rank 4 cubic hypersurface in all three flattenings, and this rank 4 locus is unliftable in all three.


R = PolynomialRing(FIELD, 7, "x")
x0, x1, x2, x3, x4, x5, x6 = R.gens()
EDerek7 = Matrix([
    [ -3*x3, -12*x4, -30*x5, -60*x6,      0,      0,      0],
    [  2*x2,   3*x3,      0, -10*x5, -30*x6,      0,      0],
    [ -2*x1,      0,   3*x3,   4*x4,      0, -12*x6,      0],
    [  3*x0,  -3*x1,  -3*x2,      0,   3*x4,   3*x5,  -3*x6],
    [     0,  12*x0,      0,  -4*x2,  -3*x3,      0,   2*x5],
    [     0,      0,  30*x0,  10*x1,      0,  -3*x3,  -2*x4],
    [     0,      0,      0,  60*x0,  30*x1,  12*x2,   3*x3]
])
Derek7 = Tensor(EDerek7)
# This is liftable in all three flattenings. It can eventually be lifted to the following:
R = PolynomialRing(FIELD, 8, "x")
x0, x1, x2, x3, x4, x5, x6, x7 = R.gens()
E888 = Matrix([
    [ 20*x0,  15*x1,   6*x2,     x3,      0,      0,      0,      0],
    [ 10*x4,  15*x5,      0,      0,   6*x2,     x3,      0,      0],
    [  4*x6,      0,   6*x5,      0,  -6*x1,      0,     x3,      0],
    [  2*x7,      0,      0,    -x5,      0,     x1,     x2,      0],
    [     0,  -3*x7,      0,    -x4,      0,   2*x0,      0,   6*x2],
    [     0,      0,  -3*x7,    -x6,      0,      0,   5*x0, -15*x1],
    [     0,      0,      0,      0, -18*x7,  -6*x6,  15*x4, -90*x5],
    [     0,     x6,    -x4,      0,   2*x0,      0,      0,    -x3]
])
T888 = Tensor(E888)
