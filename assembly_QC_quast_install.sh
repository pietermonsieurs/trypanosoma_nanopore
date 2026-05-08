## install via pip but with local settings

module load Python/3.14.2-GCCcore-15.2.0

mkdir -p /user/antwerpen/205/vsc20587/data/software/python_lib/lib/python3.14/site-packages
export PYTHONPATH=/user/antwerpen/205/vsc20587/data/software/python_lib/lib/python3.14/site-packages:${PYTHONPATH}

pip install --prefix=/user/antwerpen/205/vsc20587/data/software/python_lib quast

