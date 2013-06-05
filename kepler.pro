function kepler, M, e, max_err, iterations=iterations

;This IDL function returns the solution to Kepler's elliptic equation for the
;eccentric anomaly EA that satisfies f(EA) = EA - e*sin(EA) - M = 0.
;This algorith used Halley's method to iteratively solve for EA until
;|f(EA)|<max_err. Inputs are the mean anomaly M, eccentricity e, and the allowed
;error max_err, and M and e can be scalars or arrays. This algorithm is from Danby's
;textbook, Fundamentals of Celestial Mechanics, 2nd edition. This function does work in
;GDL as well. By Joe Hahn, jhahn@spacescience.org, 4 June 2013.

;Set the maximum number of iterations allowed.
max_iterations = 15

;Get size of inputs and convert any scalar inputs into one-element arrays.
sz = size(e)
sz0 = sz[0]
if (sz0 eq 0) then begin
    ;convert scalar inputs into one element arrays for now.
    M = [M]
    e = [e]
endif
sz = size(e)

;Initial guess
twopi = 2*!dpi
M = M - twopi*fix(M/twopi)
s = make_array(size=sz, value=1)
sinM = sin(M)
j = where(sinM lt 0d, Nj)
if (Nj gt 0) then s[j] = -1d
EA = M + s*(0.85d)*e
f = make_array(size=sz, value=1)

;This array tracks the number of iterations used to solve Kepler's equation.
iterations = make_array(size=sz, type=2)

;Solve Kepler's equation iteratively.
i = 0
Nj = n_elements(e)
j = indgen(Nj)
while((Nj gt 0) and (i lt max_iterations)) do begin $
    es = e[j]*sin(EA[j])
    ec = e[j]*cos(EA[j])
    f[j] = EA[j] - es - M[j]
    df = 1d - ec
    ddf = es
    dddf = ec
    d1 = -f[j]/df
    d2 = -f[j]/(df + d1*ddf/2d)
    d3 = -f[j]/(df + d2*ddf/2d + d2*d2*dddf/6d)
    EA[j] = EA[j] + d3
    i = i + 1
    iterations[j] = i
    j = where(abs(f) gt max_err, Nj)
endwhile
if (i ge max_iterations) then print,'** kep2.pro failed to converge! **'

;Convert back to scalars, as needed.
if (sz0 eq 0) then begin
    M = M[0]
    e = e[0]
    EA = EA[0]
    iterations = iterations[0]
endif

return, EA
end

