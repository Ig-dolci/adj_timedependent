#MEK=/home/tomas/nektar/nektar-v5.2.0/build/solvers/IncNavierStokesSolver;
MEK=~/Documents/source_nektar/nektar-v5.2.0/build3/solvers/IncNavierStokesSolver;
GMSH=~/Documents/gmsh/bin;
case="cylinder";

#$GMSH/gmsh ${case}_mesh.geo -2
#$MEK/../../utilities/NekMesh/NekMesh ${case}_mesh.msh ${case}_mesh.xml
#$MEK/../../utilities/NekMesh/NekMesh ${case}_mesh.xml ${case}_mesh.xml:xml:uncompress
#$MEK/../../utilities/NekMesh/NekMesh -m cyl:surf=100:r=0.5:N=3 ${case}_mesh.xml ${case}_mesh.xml
#$MEK/../../utilities/NekMesh/NekMesh -m bl:surf=100:layers=3:r=2:nq=7 ${case}_mesh.xml ${case}_mesh.xml
#$MEK/../../utilities/NekMesh/NekMesh ${case}_mesh.xml ${case}_mesh.xml:xml:uncompress
#$MEK/../../utilities/FieldConvert/FieldConvert ${case}_mesh.xml ${case}_mesh.vtu


testcode="uni_re100";
cd ${testcode}
cd base

mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_makeic.xml
cp -r ${case}_mesh.fld ${case}_mesh_baserst.chk
cp -r ${case}_mesh.fld ../../sin_re100/base/${case}_mesh_baserst.chk

mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_nonlinear.xml

j=0;

while [  $j -le "200000" ]; do

#echo $j
aux=$(bc <<< "200000 - $j")
newname=$(printf "%d" $aux)

mv ${case}_mesh_${j}.chk ../adjoint/${case}_nonlinear_${newname}.chk

j=$(($j+1));

done


cd ../adjoint


mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_adjoint0.xml

#$MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_adjoint0.xml

#cp ${case}_mesh.fld ../grad/${case}_mesh_adj.fld
#cp ${case}_mesh.fld ${case}_rst_adj.chk

#$MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_adjoint.xml


#cp ${case}_nonlinear_1.chk ../grad/${case}_mesh_ut0m.fld
#cp ${case}_nonlinear_39999.chk ../grad/${case}_mesh_ut0p.fld


cd ../base

#mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_nonlinear_wp.xml
#cp -r ${case}_mesh.fld ${case}_mesh_rst_wp.chk
#mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_nonlinear_wm.xml
#cp -r ${case}_mesh.fld ${case}_mesh_rst_wm.chk
#mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_nonlinear_ap.xml
#cp -r ${case}_mesh.fld ${case}_mesh_rst_ap.chk
#mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_nonlinear_am.xml
#cp -r ${case}_mesh.fld ${case}_mesh_rst_am.chk


cd ../..
testcode="sin_re100";
cd ${testcode}
cd base

#mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_makeic.xml
#cp -r ${case}_mesh.fld ${case}_mesh_baserst.chk

mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_nonlinear.xml

j=0;

while [  $j -le "96000" ]; do

#echo $j
aux=$(bc <<< "96000 - $j")
newname=$(printf "%d" $aux)

mv ${case}_mesh_${j}.chk ../adjoint/${case}_nonlinear_${newname}.chk

j=$(($j+1));

done


cd ../adjoint


mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_adjoint0.xml

#$MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_adjoint0.xml

#cp ${case}_mesh.fld ../grad/${case}_mesh_adj.fld
#cp ${case}_mesh.fld ${case}_rst_adj.chk

#$MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_adjoint.xml


#cp ${case}_nonlinear_1.chk ../grad/${case}_mesh_ut0m.fld
#cp ${case}_nonlinear_39999.chk ../grad/${case}_mesh_ut0p.fld


cd ../..

testcode="uni_re100";
cd ${testcode}
cd base


#mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_nonlinear_wp.xml
#cp -r ${case}_mesh.fld ${case}_mesh_rst_wp.chk
#mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_nonlinear_wm.xml
#cp -r ${case}_mesh.fld ${case}_mesh_rst_wm.chk
mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_nonlinear_ap.xml
cp -r ${case}_mesh.fld ${case}_mesh_rst_ap.chk
mpirun -np 6 $MEK/IncNavierStokesSolver ${case}_mesh.xml ${case}_nonlinear_am.xml
cp -r ${case}_mesh.fld ${case}_mesh_rst_am.chk


cd ../..
testcode="sin_re100";
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
