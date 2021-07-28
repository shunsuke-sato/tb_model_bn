#set term postscript eps enhanced color font "Arial, 24"
set term pdf enhanced color font "Arial, 24"

set datafile missing "NaN"

#set palette defined (0 "white", 1 "red")
set palette defined (0 "blue", 1 "red")
#set palette defined (0 "blue" "white", 1 "red")
#set palette defined (0 '#FFFF0000' , 1 '#FFFF00FF')
#set palette defined (0 '#00FF00' , 1 '#FFFF00') 
#set palette defined (0 "white", 0.5 "blue", 1 "red")
set cbrange [0:*]

set style fill transparent solid 0.50 border

set yrange [-10:10]
set xrange [0:*]
set xtics 0.2

set xlabel "Field strength (a.u.)"
set ylabel "Floquet quasienergy (a.u.)"
unset colorbox

unset key
set output "floquet_qeps.pdf"
p for [i=1:2*(2*64+1)] "floquet_band.out" u 1:(column(i*2)*27.2114):(column(i*2+1)) w l lt 1 lc palette lw 1 
#p for [i=1:2*(2*8+1)] "floquet_band.out" u 1:(column(i*2)) w l lt 1 lc rgb "red" lw 1 



unset output