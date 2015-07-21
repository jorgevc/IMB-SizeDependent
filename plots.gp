set terminal pdf
dir="DATOS_TAM/"

set output "MeanSize.pdf"
plot dir.'thinning' u 4:2 title "Mean Size vs T"
set output

set output "MeanDensity.pdf"
plot dir.'thinning' u 4:3 title "Mean Density vs T"
set output

set output "ResourceConsumption.pdf"
plot dir.'thinning' u 4:5 title "Resource Consumtion rate vs T"
set output

set output "dist_200.pdf"
plot dir."DT_2801" u (2.0*4.7**($1**(3.0/8.0))):(($2)*(($1)**(5.0/8.0)*(8.0/(10.0*2.0*4.7*3.0)))) title "dist at 200"
set output


############

reset
#set terminal wxt


set term gif animate delay 100
set output "dist.gif"

dir="DATOS_TAM/"
prefix="DT"

plots=system(sprintf("ls -v ./%s%s* | xargs -n1 basename",dir,prefix))


do for [name in plots]{
		plot [1:500][0:0.4] dir.name u 1:2 title "Size distribution ".name
  pause 1
}

set output

dir="DATOS_TAM/"
set term gif animate delay 100
set output "field.gif"


prefix="T"

plots=system(sprintf("ls -v ./%s%s* | xargs -n1 basename",dir,prefix))

set xrange [0:200]
set yrange [0:200]
set size ratio -1 
do for [name in plots]{
	plot dir.name u ($1):($2):($3==1?(4.7*$4**(3.0/8.0)):1/0) w circles title "Size ".name
 # pause 1
}
set output

