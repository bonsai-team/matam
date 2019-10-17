#!/bin/bash
set -e; #Exit immediately if a command exits with a non-zero status
set -u; #Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error

if [ $# -ne 2 ]
then
        echo "Usage : $0 WORKDIR RELEASE"
        exit 1
fi

WORKDIR=$(realpath $1); shift
RELEASE=$1; shift
CPU=$(nproc)

# The matam db has to be present
MATAM_DB="$WORKDIR/db/SILVA_128_SSURef_NR95"
if [ ! -d $(dirname $MATAM_DB) ];then
  echo "Matam db dir must be present (i.e the directory obtains with index_default_ssu_rrna_db.py script)"
  exit 1
fi

MATAM_DIR="$WORKDIR/matam/"
MATAM_BIN="$MATAM_DIR/scripts/matam_assembly.py"
MATAM_INPUT_READS="$MATAM_DIR/examples/16sp_simulated_dataset/16sp.art_HS25_pe_100bp_50x.fq"
MATAM_TRUE_REF="$MATAM_DIR/examples/16sp_simulated_dataset/16sp.fasta"
MATAM_TRUE_REF_TAXO="$MATAM_DIR/examples/16sp_simulated_dataset/16sp.taxo.tab"


mkdir -p $WORKDIR && cd $WORKDIR

echo "=== Erase previous build ==="
rm -rf $MATAM_DIR

echo "=== Clone matam ==="
cd $WORKDIR && git clone https://github.com/bonsai-team/matam.git && cd $MATAM_DIR

echo "=== Checkout the corresponding revision ==="
git checkout $RELEASE

current_revision=$RELEASE
echo "=== Working on: $current_revision ==="

echo "=== Create conda env ==="
# --force to remove previous env
#conda env create -f environment.yml --force
#source activate matam
# activate the env in subshells
set +eu # https://github.com/conda/conda/issues/8186
source ~/anaconda/etc/profile.d/conda.sh
conda activate matam
set -eu

if [ ! -n "$(command -v exonerate)" ]
then
  echo 'EXONERATE 2.2 is explicitly required'
  exit 1
fi

echo "=== Create the output dir ==="
matam_out_dir="$WORKDIR/$current_revision"
if [ -d $matam_out_dir ];then
  echo "The output dir already exists. Ignore this revision"
  exit 1
else
  echo 'outdir' $matam_out_dir
  mkdir $matam_out_dir
fi

echo "=== Build the current release ==="
# quick fix, remove when committed
# sed -i 's|find_library (LIBRT_LIB NAMES rt librt)|find_library (LIBRT_LIB NAMES rt librt HINTS "/usr/lib/x86_64-linux-gnu/")|' ovgraphbuild/CMakeLists.txt
# sed -i 's|find_library (PTHREAD_LIB NAMES pthread)|find_library (PTHREAD_LIB NAMES pthread HINTS "/usr/lib/x86_64-linux-gnu/")|' ovgraphbuild/CMakeLists.txt

./build.py #&> build.log #>/dev/null 2>&1

echo "=== Start iteration ==="
for i in `seq 1 100`;
  do
    echo "Iteration: $i"
    matam_cmd="$MATAM_BIN -i $MATAM_INPUT_READS -d $MATAM_DB -o ${matam_out_dir}/$i --true_references $MATAM_TRUE_REF --true_ref_taxo $MATAM_TRUE_REF_TAXO --cpu $CPU --max_memory 10000 --debug > ${matam_out_dir}/${current_revision}_$i.log 2>&1"
    echo $matam_cmd
    eval $matam_cmd
    rm -rf ${matam_out_dir}/$i
  done

echo "=== Erase final build ==="
rm -rf $MATAM_DIR
