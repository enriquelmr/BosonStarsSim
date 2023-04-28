

set xlabel "{/Symbol s}_0"
set ylabel "N_{max}"
set label 1 at  0.30, 0.6 "{/Symbol L}=0" center rotate by 0  front
set label 2 at  0.26, 0.72 "{/Symbol L}=1" center rotate by 0  front
set label 3 at  0.29, 0.85 "{/Symbol L}=3" center rotate by 0  front
set label 4 at  0.26, 0.99 "{/Symbol L}=10" center rotate by 0  front
set label 5 at  0.26, 1.22 "{/Symbol L}=30" center rotate by -20  front
set label 6 at  0.35, 1.25 "{/Symbol L}=100" center rotate by 0  front
set label 7 at  0.51, 2.1 "{/Symbol L}=200" center rotate by 0  front
set label 8 at  0.51, 2.6 "{/Symbol L}=300" center rotate by 0  front
plot "datosauto.dat" using 1:3 with lines notitle lt rgb "black"
