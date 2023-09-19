MEK=~/Documents/source_nektar/nektar-v5.2.0/build3/solvers/IncNavierStokesSolver;
GMSH=~/Documents/gmsh/bin;
case="cavity";

#$GMSH/gmsh ${case}_mesh.geo -2
#$MEK/../../utilities/NekMesh/NekMesh ${case}_mesh.msh ${case}_mesh.xml
#$MEK/../../utilities/NekMesh/NekMesh ${case}_mesh.xml ${case}_mesh.xml:xml:uncompress
#$MEK/../../utilities/FieldConvert/FieldConvert ${case}_mesh.xml ${case}_mesh.vtu


testcode="re4000";
cd ${testcode}
cd base

mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_makeic.xml
cp -r ${case}_mesh.fld ${case}_mesh_baserst.chk

mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_nonlinear.xml


mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_nonlinear_tan_w.xml

j=0;

while [  $j -le "4000" ]; do

echo $j
aux=$(bc <<< "4000 - $j")
newname=$(printf "%d" $aux)

mv ${case}_mesh_${j}.chk ../adjoint/${case}_nonlinear_${newname}.chk

j=$(($j+1));

done


cd ../adjoint


mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_adjoint0.xml


cd ../..


testcode="re4000";
cd ${testcode}
cd base


mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_nonlinear_wp.xml
cp -r ${case}_mesh.fld ${case}_mesh_rst_wp.chk
mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_nonlinear_wm.xml
cp -r ${case}_mesh.fld ${case}_mesh_rst_wm.chk
mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_nonlinear_ap.xml
cp -r ${case}_mesh.fld ${case}_mesh_rst_ap.chk
mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_nonlinear_am.xml
cp -r ${case}_mesh.fld ${case}_mesh_rst_am.chk


cd ../..


