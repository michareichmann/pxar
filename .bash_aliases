PXARPATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export PXARPATH="$PXARPATH"

# file operations in the config folder
alias clean_data='$PXARPATH/python/clean_data.py'
alias change_i2c='$PXARPATH/python/change_i2c.py'
alias find_i2c='$PXARPATH/python/find_i2c.py'

# CLIX
alias iclix='ipython -i $PXARPATH/python/iCLIX.py -- @$'