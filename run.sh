WORK_DIR=$(dirname $0)
QE_HOME=/f2dace-qe-source/src_f90_to_parse

python3 -u -m dace.frontend.fortran.tools.create_preprocessed_ast \
    -i $QE_HOME/dft-d3 \
    -i $QE_HOME/external/devxlib/src \
    -i $QE_HOME/external/mbd/src \
    -i $QE_HOME/FFTXlib/src \
    -i $QE_HOME/KS_Solvers \
    -i $QE_HOME/LAXlib \
    -i $QE_HOME/Modules \
    -i $QE_HOME/PW/src \
    -i $QE_HOME/PW/tools \
    -i $QE_HOME/upflib \
    -i $QE_HOME/UtilXlib \
    -i $QE_HOME/XClib \
    -o $WORK_DIR/out/out.f90 \
    -k exx_bp.vexx_bp_k_gpu \
    -d $WORK_DIR/ckpt \
    &> $WORK_DIR/out.log
