# Function to read input parameters
getValue(row,col,filename)=system('awk ''{if (NR == '.row.') print $'.col.'}'' '.filename.'')
# Setting terminal
set terminal pngcairo enhanced size 1200,800 
# Set output of file
set output "./monitor/controller.png"
# Load paths to velocity and stress data
set loadpath "./monitor"
# Height of domain in y
Ly=getValue(2,1,"./gp_input")+0
# Controller constants for title
gain=getValue(2,2,"./gp_input")+0
integral_time=getValue(2,3,"./gp_input")+0
derivative_time=getValue(2,4,"./gp_input")+0
# Set point
SP=getValue(2,5,"./gp_input")+0

# Setting plot area
set size 1,1
set origin 0,0
set multiplot layout 1,2 columnsfirst scale 1.0,0.8 title sprintf("K_c=%1.2f, {/Symbol t}_I=%1.2f, {/Symbol t}_D=%1.2e",gain,integral_time,derivative_time) font "Times Bold,30" offset 0,-3


# Plotting process variable
set xlabel "Time"
set xrange [0:integral_time]
set ylabel "Process Variable"
set y2label "error (%)"
set yrange [0:Ly]
set ytics nomirror
set y2tics nomirror
plot "controller" using 2:(SP) with points pointtype 7 ps 2 lc rgb "red"  title "SP" axis x1y1,\
     "controller" using 2:4    with points pointtype 7 ps 2 lc rgb "blue" title "PV" axis x1y1,\
     "controller" using 2:5    with points pointtype 7 ps 2 lc rgb "black" title "error" axis x1y2
unset xlabel
unset xrange
unset ylabel
unset y2label
unset yrange
unset key

# Plotting control variable
set xlabel "Time"
set xrange [0:integral_time]
set ylabel "Control Variable"
plot "controller" using 2:6 with points pointtype 7 ps 2 lc rgb "green"
unset xlabel
unset xrange
unset ylabel
unset key

unset multiplot
unset output
  