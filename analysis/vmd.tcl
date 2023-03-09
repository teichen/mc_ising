package require topotools
topo readvarxyz traj.xyz

pbc set {15.0 15.0 15.0} -all

display depthcue off
axes location off
mol modstyle 0 0 VDW 0.5000 12.0000
color Display Background white
display resize 1200 1200
scale by 2.0

mol color Name
mol representation VDW 0.5 12.0
mol selection (name O) and user > 0
mol material Opaque
mol addrep 0

mol color Name
mol representation VDW 0.5 12.0
mol selection (name N) and user > 0
mol material Transparent
mol addrep 0

