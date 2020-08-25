local cmath = require("complex")
local vmath = require("vector")
local mmath = require("matrix")

local qu = {}

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

qu.sHV = H:dot(H:bra()) - V:dot(V:bra())
qu.sDA = D:dot(D:bra()) - A:dot(A:bra())
qu.sRL = R:dot(R:bra()) - L:dot(L:bra())
qu.pH = H:dot(H:bra())
qu.pV = V:dot(V:bra())
qu.pD = D:dot(D:bra())
qu.pA = A:dot(A:bra())
qu.pR = R:dot(R:bra())
qu.pL = L:dot(L:bra())

----------------------------------------------------------------

local zero, one, plus, minus, plusi, minusi = H, V, D, A, R, L

qu.X_phi = function(phi)
	return plus:dot(plus:bra()) + minus:dot(minus:bra()) * ((phi or math.pi) * 1i):exp()
end

qu.Y_phi = function(phi)
	return plusi:dot(plusi:bra()) + minusi:dot(minusi:bra()) * ((phi or math.pi) * 1i):exp()
end

qu.Z_phi = function(phi)
	return zero:dot(zero:bra()) + one:dot(one:bra()) * ((phi or math.pi) * 1i):exp()
end

qu.psi_theta_phi = function(theta, phi)
	return zero * math.cos(theta / 2) + one * (phi * 1i):exp() * math.sin(theta / 2)
end

qu.theta_phi = function(v)
	v = v:norm() * cmath.complexpolar(1, v.x:conj():arg())
	local cost2 = v.x:abs()
	local eiphi = v.y / (1 - cost2 ^ 2) ^ 0.5
	return math.acos(cost2) * 2, eiphi:arg()
end

qu.get_orth = function(psi)
	local theta, phi = qu.theta_phi(psi)
	return qu.psi_theta_phi(math.pi - theta, math.pi + phi)
end

qu.U = function(psi, phi)
	local psiorth = qu.get_orth(psi)
	return psi:dot(psi:bra()) + psiorth:dot(psiorth:bra()) * ((phi or math.pi) * 1i):exp()
end

qu.H_op = mmath.matrix(1, 1, 1, -1) / sq2

----------------------------------------------------------------

qu.NOT = qu.X_phi()

qu.CX = (zero:dot(zero:bra()):rightkron(identity(2)) + one:dot(one:bra()):rightkron(qu.X_phi())):clear()
qu.CY = (zero:dot(zero:bra()):rightkron(identity(2)) + one:dot(one:bra()):rightkron(qu.Y_phi())):clear()
qu.CZ = (zero:dot(zero:bra()):rightkron(identity(2)) + one:dot(one:bra()):rightkron(qu.Z_phi())):clear()
qu.CNOT = qu.CX

qu.CU = function(U_op)
	return zero:dot(zero:bra()):rightkron(identity(2)) + one:dot(one:bra()):rightkron(U_op)
end

local zerozero = zero:kron(zero)
local oneone = one:kron(one)
local zeroone = zero:kron(one)
local onezero = one:kron(zero)
local minusplus = minus:kron(plus)
local pluszero = plus:kron(zero)
local plusione = plusi:kron(one)

qu.SWAP = qu.CNOT:rightmatrix(
	zero:dot(zero:bra()):leftkron(identity(2)) + one:dot(one:bra()):leftkron(qu.X_phi())
):rightmatrix(qu.CNOT):clear()

----------------------------------------------------------------

qu.U_theta_delta = function(theta, delta)
	local ketTheta = H * math.cos(theta) + V * math.sin(theta)
	local ketThetaPi2 = H * math.cos(theta + math.pi / 2) + V * math.sin(theta + math.pi / 2)
	return ketTheta:dot(ketTheta:bra()) + ketThetaPi2:dot(ketThetaPi2:bra()) * (delta * 1i):exp()
end

----------------------------------------------------------------

local basis = {H, V}
local basis2 = {{zerozero, zeroone}, {onezero, oneone}}

local rho_ijmn = function(s, i, j, m, n)
	return s:ratio(basis2[i][m]:dot(basis2[j][n]:bra()))
end

qu.subtrace = function(s, subsystem)
	local e = basis
	local op = matrix(0, 0, 0, 0)
	for i = 1, 2 do
		for j = 1, 2 do
			local sum = 0
			for m = 1, 2 do
				for n = 1, 2 do
					local r = subsystem == 2 and rho_ijmn(s, i, j, m, n) or rho_ijmn(s, m, n, i, j)
					sum = sum + r * e[m]:bra():dot(e[n])
				end
			end
			op = op + e[i]:dot(e[j]:bra()) * sum
		end
	end
	return op
end

----------------------------------------------------------------

return qu
