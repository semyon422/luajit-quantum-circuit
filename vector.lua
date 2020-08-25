local cmath = require("complex")
local abs = cmath.abs
local sqrt = cmath.sqrt
local conj = cmath.conj
local tocomplex = cmath.tocomplex

--------------------------------------------------------------------------------

local mt = {}

local vector = function(...)
	local vector = {...}
	vector.row = false

	setmetatable(vector, mt)

	return vector
end

local isvector = function(a)
	return getmetatable(a) == mt and a
end

--------------------------------------------------------------------------------

local _copy = function(a)
	return vector(unpack(a)):setrow(a.row)
end

local _abs = function(a)
	local sum = 0
	for i = 1, #a do
		sum = sum + abs(a[i]) ^ 2
	end
	return math.sqrt(sum)
end

local _conj = function(a)
	local c = _copy(a)
	for i = 1, #a do
		c[i] = conj(a[i])
	end
	return c
end

local _norm = function(a)
	return a / _abs(a)
end

local _transpose = function(a)
	return _copy(a):setrow(not a.row)
end

local _bra = function(a)
	if not a.row then
		return _conj(_transpose(a))
	end
	return _copy(a)
end

local _ket = function(a)
	if a.row then
		return _conj(_transpose(a))
	end
	return _copy(a)
end

local _dot = function(a, b)
	assert(isvector(b))
	if not a.row and b.row then
		local mmath = require("matrix")
		local matrix = mmath.matrix()
		for i = 1, #a do
			for j = 1, #a do
				matrix[(i - 1) * #a + j] = a[i] * b[j]
			end
		end
		return matrix
	elseif a.row and not b.row then
		local sum = 0
		for i = 1, #a do
			sum = sum + a[i] * b[i]
		end
		return sum
	end
	error("attempt to multiply vectors of the same type")
end

local _kron = function(a, b)
	assert(isvector(b))
	if not a.row and b.row or a.row and not b.row then
		error("attempt to multiply vectors of the different types")
	end
	local c = vector()
	for i = 1, #a do
		for j = 1, #b do
			c[(i - 1) * #a + j] = a[i] * b[j]
		end
	end
	return c
end

local _ratio = function(a, b)
	assert(isvector(b))
	local t = {}
	for i = 1, #a do
		local r = a[i] / b[i]
		if r == r then
			t[#t + 1] = r
		end
	end
	for i = 2, #t do
		if abs(t[i] - t[1]) > 1e-12 then
			return tocomplex(0)
		end
	end
	return t[1]
end

local _setrow = function(a, row)
	a.row = row
	return a
end

local _clear = function(a)
	local c = vector()
	for i = 1, #a do
		local t = tocomplex(a[i])
		local re = math.abs(t.re) < 1e-12 and 0 or t.re
		local im = math.abs(t.im) < 1e-12 and 0 or t.im
		c[i] = re + im * 1i
	end
	assert(a == c)
	return c
end

--------------------------------------------------------------------------------

local vmath = {}
local supervector = {}

vmath.vector = vector
vmath.isvector = isvector

supervector.copy		= _copy
supervector.abs			= _abs
supervector.conj		= _conj
supervector.norm		= _norm
supervector.transpose	= _transpose
supervector.bra			= _bra
supervector.ket			= _ket
supervector.dot			= _dot
supervector.ratio		= _ratio
supervector.setrow		= _setrow
supervector.kron		= _kron
supervector.clear		= _clear

--------------------------------------------------------------------------------

mt.__index = function(_, key)
	return supervector[key]
end

mt.__eq = function(a, b)
	assert(isvector(a))
	assert(isvector(b))
	return _abs(a - b) < 1e-12
end

mt.__unm = function(a)
	local c = _copy(assert(isvector(a)))
	for i = 1, #a do
		c[i] = -a[i]
	end
	return a
end

mt.__add = function(a, b)
	local c = _copy(assert(isvector(a)))
	assert(isvector(b))
	for i = 1, #a do
		c[i] = a[i] + b[i]
	end
	return c
end

mt.__sub = function(a, b)
	local c = _copy(assert(isvector(a)))
	assert(isvector(b))
	for i = 1, #a do
		c[i] = a[i] - b[i]
	end
	return c
end

mt.__mul = function(a, b)
	local c = _copy(assert(isvector(a)))
	b = assert(tocomplex(b), "second argument should be 'complex' or 'number'")
	for i = 1, #a do
		c[i] = a[i] * b
	end
	return c
end

mt.__div = function(a, b)
	local c = _copy(assert(isvector(a)))
	b = assert(tocomplex(b), "second argument should be 'complex' or 'number'")
	for i = 1, #a do
		c[i] = a[i] / b
	end
	return c
end

mt.__mod = function()
	return error("__mod is not implemented")
end

mt.__lt = function()
	return error("__lt is not implemented")
end

mt.__le = function()
	return error("__le is not implemented")
end

mt.__concat = function(a, b)
	return tostring(a) .. tostring(b)
end

mt.__tostring = function(a)
	local d = a.row and ";" or ","
	local t = a[1]
	for i = 2, #a do
		t = t .. d .. a[i]
	end
	return t
end

return vmath
