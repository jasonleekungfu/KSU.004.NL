set xlabel 'Wavelength [nm]'
set ylabel 'Intensity enhancement'
set yrange[0:70]

set title 'Au nanosphere'
plot 'enh_wl-solid0.dat' w l t "8x10^{10} W/cm^2", 'enh_wl-solid1.dat' w l t "6x10^{11} W/cm^2", 'enh_wl-solid2.dat' w l t "1.2x10^{12} W/cm^2", 'enh_wl-solid3.dat' w l t "2.7x10^{12} W/cm^2", 'enh_wl-solid4.dat' w l t "8x10^{12} W/cm^2"

set title 'SiO_2 core - Au Shell'
plot 'enh_wl-shell0.dat' w l t "8x10^{10} W/cm^2", 'enh_wl-shell1.dat' w l t "6x10^{11} W/cm^2", 'enh_wl-shell2.dat' w l t "1.2x10^{12} W/cm^2", 'enh_wl-shell3.dat' w l t "2.7x10^{12} W/cm^2", 'enh_wl-shell4.dat' w l t "8x10^{12} W/cm^2"

