package.path = package.path .. ";?/init.lua"

local cmath = require("complex")
local vmath = require("vector")
local mmath = require("matrix")
local qu = require("quantum")

local vector = vmath.vector
local matrix = mmath.matrix
local identity = mmath.identity

----------------------------------------------------------------

local H = vector(1, 0)
local V = vector(0, 1)
local D = (H + V):norm()
local A = (H - V):norm()
local R = (H + V * 1i):norm()
local L = (H - V * 1i):norm()

local sq2 = math.sqrt(2)

----------------------------------------------------------------

local zero, one, plus, minus, plusi, minusi = H, V, D, A, R, L

local zerozero = zero:kron(zero)
local oneone = one:kron(one)
local zeroone = zero:kron(one)
local onezero = one:kron(zero)
local minusplus = minus:kron(plus)
local pluszero = plus:kron(zero)
local plusione = plusi:kron(one)

----------------------------------------------------------------

assert(cmath.abs(qu.sHV:leftvector(H:bra()):dot(H) - 1) < 1e-12)
assert(cmath.abs(qu.sHV:leftvector(R:bra()):dot(R)) < 1e-12)

assert(qu.CNOT:rightvector(minusplus):ratio(minusplus) == 1) -- |-+> is eigenvector for CNOT
assert(qu.CNOT:rightvector(pluszero):ratio((zerozero + oneone) / sq2) == 1) -- entanglement state
assert(qu.CZ:rightvector(oneone):ratio(oneone) == -1) -- -1

local m = H:dot(H:bra()) - V:dot(V:bra())
assert(m:hermit():rightmatrix(m) == identity(2)) -- hermitian operator

local m = -R:dot(D:bra())
local e1, e2 = m:braeigenvalues()
assert(cmath.abs(e1) < 1e-12)
assert(cmath.abs(e2 - (-1 + 1i) / 2) < 1e-12)

local m = -R:dot(D:bra())
assert(not m:iseigenvector(H))
assert(not m:iseigenvector(V))
assert(not m:iseigenvector(D))
assert(m:iseigenvector(A))
assert(m:iseigenvector(R))
assert(not m:iseigenvector(L))

local ph = (H * (-3) + V * 4i):norm()
local sRL = R:dot(R:bra()) - L:dot(L:bra())
local average = (sRL:rightmatrix(sRL)):leftvector(ph:bra()):dot(ph) -- <ph|sRL^2|ph>
assert(cmath.abs(average - 1) < 1e-12)

-- local psi = qu.psi_theta_phi(math.pi / 4, -math.pi / 2) -- ?
-- local H = U(psi, math.pi)
-- print(H * sq2)

assert(cmath.abs(qu.H_op:rightvector(zero):ratio(plus) - 1))
assert(cmath.abs(qu.H_op:rightvector(one):ratio(minus) - 1))
assert(cmath.abs(qu.H_op:rightvector(plus):ratio(zero) - 1))
assert(cmath.abs(qu.H_op:rightvector(minus):ratio(one) - 1))
assert(cmath.abs(qu.H_op:rightvector(plusi):ratio(minusi) - cmath.exp(1i * math.pi / 4)))
assert(cmath.abs(qu.H_op:rightvector(minusi):ratio(plusi) - cmath.exp(1i * math.pi / 4)))

local rho1 = (qu.pH + qu.pV) * 0.5 -- entanglement state
local ph = (H * 2 + R * 3):norm()
assert(cmath.abs(rho1:rightvector(ph):ratio(ph) - 0.5)  < 1e-12) -- 0.5 for any vector => any vector is eigenvector
assert(cmath.abs((qu.pH:rightmatrix(qu.X_phi())):clear():trace()) < 1e-12) -- <X> = Tr[pH*X] = 0
assert(cmath.abs((qu.pH:rightmatrix(qu.Y_phi())):clear():trace()) < 1e-12) -- <Y> = Tr[pH*Y] = 0
assert(cmath.abs((qu.pH:rightmatrix(qu.Z_phi())):clear():trace() - 1) < 1e-12) -- <Z> = Tr[pH*Z] = 1

local state = (zerozero + oneone) / sq2 -- entanglement state
local rho12 = state:dot(state:bra())
local averageZ = (qu.subtrace(rho12, 2):rightmatrix(qu.X_phi())):trace()
assert(cmath.abs(averageZ) < 1e-12)

local id_Z = identity(2):rightkron(qu.Z_phi())
local id_H = identity(2):rightkron(qu.H_op)
local op = id_H:leftmatrix(qu.CZ):leftmatrix(id_Z):leftmatrix(qu.CZ):leftmatrix(id_H)
assert(op:rightvector(zerozero) == zeroone)
assert(op:rightvector(zeroone) == zerozero)
assert(op:rightvector(onezero) == oneone)
assert(op:rightvector(oneone) == onezero)

local F = cmath.abs((zero * 0.5 + one * (0.8 + 0.3i)):bra():dot(plus)) ^ 2
assert(F == 0.89) -- fidelity (<psi|plus>) ^ 2
