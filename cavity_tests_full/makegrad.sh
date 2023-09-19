MEK=~/Documents/source_nektar/nektar-v5.2.0/build3/solvers/IncNavierStokesSolver;
GMSH=~/Documents/gmsh/bin;
case="cavity";



testcode="re4000_h";
cd ${testcode}
cd grad


python3 grad_avgdrag.py
python3 grad_fd_a.py
python3 grad_fd_w.py


