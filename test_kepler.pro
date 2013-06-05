;test_kepler.pro
;by Joe Hahn, jhahn@spacescience.org, 4 June 2013.

;This IDL script calls kepler() and makes some plots to show that the
;numerical errors in this call to kepler(), which solves Kepler's equation numerically,
;are within expected limits. All results are plotted
;against x=1-e where log10(x) is uniformly distributed over -6<log10(x)<0
;so that eccentricity e ranges over 10^(-6)<e<0.999999. The figure shows that
;kepler() has to iterate more when e is close to one. In summary,
;when N=10^7 (which is the number of solutions that kepler() must generate)
;and max_error=1e-15, all numerical errors are<max_error as expected.
;The execution time on my 6 year old Dell Precision T3400 is 8.33 seconds.
;To execute, start IDL and enter @keplers_eqn

;Solve Kepler's equation N times
N = 10000000l

;Plot results for Nplot solutions
Nplot = 100000

;Allowed error in the solution to Kepler's equation.
err_max = 1d-15

;Create random orbits with eccentricity e logarithmically distributed over
;-6<log10(1-e)<0 so that e ranges over 10^(-6)<e<1 and mean anomaly M is
;uniformly distributed over -Pi<M<Pi.
seed = 1
x = 10d^(-6*randomu(seed, N))
e = 1 - x
M = !dpi*(2*randomu(seed, N) - 1)

;Use kepler() to solve Kepler's equation. Also get execution time and errors in that solution.
t0_sec = systime(seconds=1)
EA = kepler(M, e, err_max, iterations=iterations)
t1_sec = systime(seconds=1)
dt_sec = t1_sec - t0_sec
err = abs(EA - e*sin(EA) - M)

;Plot the error in the solution for kepler's equation versus x=1-e. Also show
;that the mean anomalies M are uniformly distributed over -Pi<M<Pi,
;and plot the number of iterations that kepler() used to solve Kepler's equation.
!p.multi = [0, 1, 3]
window, xs=550, ys=850, retain=2
j = indgen(Nplot)
plot_oo, x[j], err[j], psym=3, xrange=[1d-6, 1], yrange=[1d-18, 1d-12], xstyle=1, ystyle=1, $
    xtitle='x = 1 - e', ytitle='error', charsize=3.0
plot_oi, x[j], M[j]/!dpi, psym=3, xrange=[1d-6, 1], yrange=[-1, 1], xstyle=1, ystyle=1, $
    xtitle='x = 1 - e', ytitle='mean anomaly M/!7p!3', charsize=3.0
plot_oi, x[j], iterations[j], psym=3, xrange=[1d-6, 1], yrange=[0, 7], xstyle=1, ystyle=1, $
    xtitle='x = 1 - e', ytitle='number of iterations', charsize=3.0
oplot, x[j], iterations[j], psym=6, symsize=0.4, color=128
!p.multi = 0

print, 'maximum error in solution    = ', max(err[j])
print, 'maximum number of iterations = ', max(iterations[j])
print, 'execution time (sec)         = ', dt_sec

stop
stop

;The following confirms that kepler() works for scalar inputs also
e = 0.2d
M = 1.5d
EA = kepler(M, e, err_max, iterations=iterations)
help, e, M, EA, iterations
print, 'max error in scalar solution = ', abs(EA - e*sin(EA) - M)
print, 'maximum number of iterations = ', max(iterations)


