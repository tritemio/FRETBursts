"""
Plot FRET variation when D or A is changed + or - 10%
"""

k = r_[0:0.1:0.01,0.1:10:0.1,10:100, 1000]

E1 = 1/(1+k*1.1)
E2 = 1/(1+k/(1.1))
E3 = 1/(1+k*0.9)
E4 = 1/(1+k/0.9)

plot(E, (E-E1)/E)
plot(E, (E-E2)/E)
plot(E, (E-E3)/E)
plot(E, (E-E4)/E)
legend(["D +10%", "A +10%", "D -10%", "A -10%"])
title("FRET rel. variation on D or A variations")
xlabel("E"); ylabel(r"$\frac{1}{1+D/A} - \frac{1}{1+D^*/A^*}$", fontsize=16)
savefig("FRET Error.png", bbox_inches='tight')

