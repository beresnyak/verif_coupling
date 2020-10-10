unset logscale
a = "circuit.txt"
b = "ver.txt"
plot a u 1:2 t "V_C", a u 1:3 t "I", a u 1:4 t "V_{device}", b u 1:3 t "theory", b u 1:4 t "theory", b u 1:5 t "theory"
pause -1
