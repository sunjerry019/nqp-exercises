#!/usr/bin/env fish

# run source init_modules.fish

# install omf using 
# $ curl https://raw.githubusercontent.com/oh-my-fish/oh-my-fish/master/bin/install > install
# $ fish install --path=~/.local/share/omf --config=~/.config/omf

source $MODULESHOME/init/fish

set --export MKLROOT "/software/opt/focal/x86_64/spack/2021.12/spack/opt/spack/linux-ubuntu20.04-x86_64_v2/gcc-9.3.0/intel-oneapi-mkl-2021.4.0-finhpnl2dqry3iuhcvtrpwagvgvngpl5/mkl/2021.4.0"
set --export LIBMKL "$MKLROOT/lib"
set --export INTELMKL "$LIBMKL/intel64"

# https://superuser.com/a/1235985
bass source $MKLROOT/env/vars.sh -arch intel64

module purge
if not module is-loaded spack/2021.12; module load spack/2021.12; end
if not module is-loaded intel-oneapi-mkl/2021.4.0; module load intel-oneapi-mkl/2021.4.0; end
if not module is-loaded python/3.9-2021.11; module load python/3.9-2021.11; end
