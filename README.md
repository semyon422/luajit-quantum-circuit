# luajit-quantum-circuit
```lua
-- CNOT |+0> = 1 / sqrt(2) * (|00> + |11>)
assert(qu.CNOT:rightvector(pluszero):ratio((zerozero + oneone) / sq2) == 1) -- entanglement state

-- beam splitter sHV = |H><H| - |V><V|
local sHV = H:dot(H:bra()) - V:dot(V:bra())

-- average value of sHV for |R>: <R|sHV|R> = 0
assert(cmath.abs(qu.sHV:leftvector(R:bra()):dot(R)) < 1e-12)
```
