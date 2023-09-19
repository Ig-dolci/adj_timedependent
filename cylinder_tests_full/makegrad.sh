MEK=/home/tomas/nektar/nektar-v5.2.0/build/solvers/IncNavierStokesSolver;
GMSH=~/Documents/gmsh/bin;
case="cavity";

#$GMSH/gmsh ${case}_mesh.geo -2
#$MEK/../../utilities/NekMesh/NekMesh ${case}_mesh.msh ${case}_mesh.xml
#$MEK/../../utilities/NekMesh/NekMesh ${case}_mesh.xml ${case}_mesh.xml:xml:uncompress
#$MEK/../../utilities/NekMesh/NekMesh -m cyl:surf=100:r=0.5:N=3 ${case}_mesh.xml ${case}_mesh.xml
#$MEK/../../utilities/NekMesh/NekMesh -m bl:surf=100:layers=3:r=2:nq=7 ${case}_mesh.xml ${case}_mesh.xml
#$MEK/../../utilities/NekMesh/NekMesh ${case}_mesh.xml ${case}_mesh.xml:xml:uncompress
#$MEK/../../utilities/FieldConvert/FieldConvert ${case}_mesh.xml ${case}_mesh.vtu


testcode="re2000_hp";
cd ${testcode}
cd grad

#cp ../base/${case}_mesh.xml ./

#$MEK/../../utilities/FieldConvert/FieldConvert -m addfld:fromfld=${case}_mesh_ut0m.fld:scale=-1 ${case}_mesh.xml ${case}_mesh_ut0p.fld ${case}_deltau.fld
#$MEK/../../utilities/FieldConvert/FieldConvert -m scaleinputfld:scale=5000 ${case}_mesh.xml ${case}_deltau.fld ${case}_dudt-scal.fld


#$MEK/../../utilities/FieldConvert/FieldConvert -m innerproduct:fromfld=${case}_mesh_adj.fld:fields="0" ${case}_mesh.xml ${case}_dudt-scal.fld out.stdout
#$MEK/../../utilities/FieldConvert/FieldConvert -m innerproduct:fromfld=${case}_mesh_adj.fld:fields="1" ${case}_mesh.xml ${case}_dudt-scal.fld out.stdout

python3 grad_avgdrag.py
python3 grad_fd_a.py
python3 grad_fd_w.py


