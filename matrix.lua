local cmath = require("complex")
local floor = math.floor
local conj = cmath.conj
local abs = cmath.abs
local tocomplex = cmath.tocomplex

--------------------------------------------------------------------------------

local mt = {}

local matrix = function(...)
	local matrix = {...}

	setmetatable(matrix, mt)

	return matrix
end

local ismatrix = function(m)
	return getmetatable(m) == mt and m
end

local identity = function(size)
	local m = matrix()
	for i = 1, size do
		for j = 1, size do
			m[(i - 1) * size + j] = i == j and 1 or 0
		end
	end
	return m
end

local assert_int = function(a)
	return assert(a and floor(a) == a, "k should be integer") and a
end

--------------------------------------------------------------------------------

local _copy = function(m)
	local c = matrix()
	for i = 1, #m do
		c[i] = m[i]
	end
	return c
end

local _conj = function(m)
	local c = matrix()
	for i = 1, #m do
		c[i] = conj(m[i])
	end
	return c
end

local _transpose = function(m)
	local c = matrix()
	local size = math.sqrt(#m)
	for i = 1, size do
		for j = 1, size do
			c[(i - 1) * size + j] = m[(j - 1) * size + i]
		end
	end
	return c
end

local _hermit = function(m)
	return _conj(_transpose(m))
end

local _rightvector = function(m, v)
    local vmath = require("vector")
	assert(vmath.isvector(v))
    assert(not v.row, "vector should be a column")
	local size = math.sqrt(#m)
	assert(size == #v, "matrix and vector should be same size")
	local c = vmath.vector():setrow(v.row)
	for i = 1, size do
		local sum = 0
		for j = 1, size do
			sum = sum + m[(i - 1) * size + j] * v[j]
		end
		c[i] = sum
	end
	return c
end

local _leftvector = function(m, v)
    local vmath = require("vector")
	assert(vmath.isvector(v))
    assert(v.row, "vector should be a row")
	local size = math.sqrt(#m)
	assert(size == #v, "matrix and vector should be same size")
	local c = vmath.vector():setrow(v.row)
	for i = 1, size do
		local sum = 0
		for j = 1, size do
			sum = sum + m[(j - 1) * size + i] * v[j]
		end
		c[i] = sum
	end
	return c
end

local _rightmatrix = function(m, n)
	assert(#m == #n, "tables should be same size")
	local size = math.sqrt(#m)
	local c = matrix()
	for i = 1, size do
		for j = 1, size do
			local sum = 0
			for k = 1, size do
				sum = sum + m[(i - 1) * size + k] * n[(k - 1) * size + j]
			end
			c[(i - 1) * size + j] = sum
		end
	end
	return c
end

local _leftmatrix = function(m, n)
	return _rightmatrix(n, m)
end

local _det = function(m)
	assert(#m == 4)
	return m[1] * m[4] - m[2] * m[3]
end

local _keteigenvalues = function(m)
	assert(#m == 4)
	local D = (m[1] + m[4]) ^ 2 - 4 * (m[1] * m[4] - m[2] * m[3])
	return ((m[1] + m[4]) + D ^ 0.5) / 2, ((m[1] + m[4]) - D ^ 0.5) / 2
end

local _braeigenvalues = function(m)
	assert(#m == 4)
	m = _hermit(m)
	local D = (m[1] + m[4]) ^ 2 - 4 * (m[1] * m[4] - m[2] * m[3])
	return ((m[1] + m[4]) + D ^ 0.5) / 2, ((m[1] + m[4]) - D ^ 0.5) / 2
end

local _iseigenvector = function(m, v)
	assert(#m == 4)
    local vmath = require("vector")
	assert(vmath.isvector(v))
	local e1, e2, n1, n2
	if not v.row then
		e1, e2 = _keteigenvalues(m)
		n1 = _rightvector(m - identity(2) * e1, v)
		n2 = _rightvector(m - identity(2) * e2, v)
	else
		e1, e2 = _braeigenvalues(m)
		n1 = _leftvector(m - identity(2) * e1, v)
		n2 = _leftvector(m - identity(2) * e2, v)
	end
	return n1:abs() < 1e-12 or n2:abs() < 1e-12
end

local _rightkron = function(m, n)
	local sizem = math.sqrt(#m)
	local sizen = math.sqrt(#n)
	local size = sizem * sizen
	local c = matrix()
	for i = 1, size do
		for j = 1, size do
			local i_m = math.ceil(i / sizem)
			local j_m = math.ceil(j / sizem)
			local i_n = (i - 1) % sizen + 1
			local j_n = (j - 1) % sizen + 1
			c[(i - 1) * size + j] = m[(i_m - 1) * sizem + j_m] * n[(i_n - 1) * sizen + j_n]
		end
	end
	return c
end

local _leftkron = function(m, n)
	return _rightkron(n, m)
end

local _clear = function(m)
	local c = matrix()
	for i = 1, #m do
		local t = tocomplex(m[i])
		local re = abs(t.re) < 1e-12 and 0 or t.re
		local im = abs(t.im) < 1e-12 and 0 or t.im
		c[i] = re + im * 1i
	end
	assert(m == c)
	return c
end

local _trace = function(m)
	local size = math.sqrt(#m)
	local sum = 0
	for i = 1, size do
		sum = sum + m[(i - 1) * size + i]
	end
	return sum
end

local _ratio = function(a, b)
	assert(ismatrix(b))
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

--------------------------------------------------------------------------------

local mmath = {}
local supermatrix = {}

mmath.matrix = matrix
mmath.ismatrix = ismatrix
mmath.identity = identity

supermatrix.copy			= _copy
supermatrix.conj			= _conj
supermatrix.transpose		= _transpose
supermatrix.hermit  		= _hermit
supermatrix.leftvector		= _leftvector
supermatrix.rightvector		= _rightvector
supermatrix.leftmatrix		= _leftmatrix
supermatrix.rightmatrix		= _rightmatrix
supermatrix.det				= _det
supermatrix.keteigenvalues	= _keteigenvalues
supermatrix.braeigenvalues	= _braeigenvalues
supermatrix.iseigenvector	= _iseigenvector
supermatrix.rightkron		= _rightkron
supermatrix.leftkron		= _leftkron
supermatrix.clear			= _clear
supermatrix.trace			= _trace
supermatrix.ratio			= _ratio
--------------------------------------------------------------------------------

mt.__index = function(_, key)
	return supermatrix[key]
end

mt.__eq = function(m, n)
	assert(#m == #n, "tables should be same size")
	local sum = 0
	for i = 1, #m do
		sum = sum + cmath.abs(m[i] - n[i])
	end
	return sum / #m < 1e-12
end

mt.__unm = function(m)
	local c = matrix()
	for i = 1, #m do
		c[i] = -m[i]
	end
	return c
end

mt.__add = function(m, n)
	assert(#m == #n, "tables should be same size")
	local c = matrix()
	for i = 1, #m do
		c[i] = m[i] + n[i]
	end
	return c
end

mt.__sub = function(m, n)
	return mt.__add(m, mt.__unm(n))
end

mt.__mul = function(m, n)
	n = assert(tocomplex(n), "second argument should be 'complex' or 'number'")
	local c = matrix()
	for i = 1, #m do
		c[i] = m[i] * n
	end
	return c
end

mt.__div = function(m, n)
	n = assert(tocomplex(n), "second argument should be 'complex' or 'number'")
	local c = matrix()
	for i = 1, #m do
		c[i] = m[i] / n
	end
	return c
end

mt.__pow = function(m, k)
	if ismatrix(m) and assert_int(k) then
		local n = m:copy()
		for i = 2, k do
			n = n:rightmatrix(m)
		end
		return n
	end
	error("second argument should be 'number'")
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

mt.__tostring = function(m)
	local size = math.sqrt(#m)
	local s = "matrix"
	for i = 1, #m do
		if i % size == 1 then
			s = s .. "\n"
		end
		s = s .. m[i] .. " "
	end
	return s
end

mt.__call = function(m, i, j, value)
	local size = math.sqrt(#m)
	if value then
		m[(i - 1) * size + j] = value
		return m
	end
	return m[(i - 1) * size + j]
end

return mmath
