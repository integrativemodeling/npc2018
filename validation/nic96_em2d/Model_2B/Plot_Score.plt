reset
set terminal pdfcairo enhanced color font "Arial-Bold, 40" size 10,10
set border lw 3 lc rgb "#484848"

set boxwidth 1.0 absolute
set style fill transparent solid 0.50 border rgb "#484848"
unset key 


set encoding iso_8859_1

set ylabel "CCC" tc rgb "#484848" offset 0,0 font "Arial-Bold, 51"
set xlabel "2D class average" tc rgb "#484848" offset 0,0 font "Arial-Bold, 57"

set yr [0:1] noreverse nowriteback
set xr [-1:34] noreverse nowriteback

set pointsize 2
set xtics 1,1,33 border in scale 0,0 mirror norotate  offset character 0, 0, 0 autojustify tc rgb "#484848"

set xtics   ()
set ytics 0, .25, 1  nomirror tc rgb "#484848"
set tic scale 0


set format y "%.2f"

set key tc rgb "#484848"
set output "Score.pdf"
plot "C1_logs_35.txt" usi 3:5 w boxes lc rgb "#111111" notitle

set output