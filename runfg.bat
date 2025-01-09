G:

cd G:\FlightGear

SET FG_ROOT=G:\FlightGear\data

START .\\bin\fgfs.exe --fdm=null --native-fdm=socket,in,30,localhost,5502,udp --native-ctrls=socket,out,30,localhost,5505,udp   --aircraft=747-400 --fog-fastest --disable-clouds --start-date-lat=2004:06:01:09:00:00 --disable-sound --in-air --airport=LOWI --runway=10L --altitude=0 --heading=0 --offset-distance=4.72 --offset-azimuth=0  
