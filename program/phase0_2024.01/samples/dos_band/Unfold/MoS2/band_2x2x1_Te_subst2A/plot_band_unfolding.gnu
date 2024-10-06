set output "unfolded_band.png"
set terminal pngcairo size 320, 480 
set nokey
set xr [0.00:1.13835448983944] 
set yr [-5:5] 
set palette defined ( 0 'white', 1 'light-salmon', 2 'red' ) 
set arrow from 0.0, 0,0 to 1.13835448983944, 0.0 nohead lt 0.2 
set arrow from 0.416666665851287,-5 to 0.416666665851287,5 nohead lt -1
set arrow from 0.657229278177461,-5 to 0.657229278177461,5 nohead lt -1
set ylabel "Energy (eV)"
set xtics ( " G" 0, " M" 0.416666665851287, " K" 0.657229278177461, " G" 1.13835448983944)
plot "plot_band_energy.dat" u 1:2:(0.02*$3) w circles fill transparent solid 0.8 noborder lc rgb "red", "plot_band_energy.dat" u 1:2 w l lw 0.1 lc "black" 
