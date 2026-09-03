# .bashrc

# Cookbook for future reference

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific environment
if ! [[ "$PATH" =~ "$HOME/.local/bin:$HOME/bin:" ]]
then
    PATH="$HOME/.local/bin:$HOME/bin:$PATH"
fi
export PATH

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=

# User specific aliases and functions
# export PATH="/project2/andrewferguson/srusti/venvs/miniconda/bin:$PATH"  # commented out by conda initialize

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/project2/andrewferguson/srusti/venvs/miniconda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/project2/andrewferguson/srusti/venvs/miniconda/etc/profile.d/conda.sh" ]; then
        . "/project2/andrewferguson/srusti/venvs/miniconda/etc/profile.d/conda.sh"
    else
        export PATH="/project2/andrewferguson/srusti/venvs/miniconda/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

alias prepgmx='module load openmpi/4.1.2+gcc-10.2.0 gsl/2.5 cuda/11.2 python'
alias ter='exec bash -l'
alias rcc-squeue='squeue -u srustidonapati'
alias ter='exec bash -l'
alias py='module load python/anaconda-2023.09'
alias prepgmx='module load openmpi/4.1.2+gcc-10.2.0 gsl/2.5 cuda/11.2 python'
alias gmx='/project/sachar/gromacs/gromacs/gromacs-2021.6/installed-files/bin/gmx'
alias project2='cd /project2/andrewferguson/srusti/'
alias lrt='ls -lrt'
alias load_gmx='source ~/bin/load_gmx'


# alias prepgmx='module load openmpi/4.1.2+gcc-10.2.0 gsl/2.5 cuda/11.2 python/anaconda-2023.09'
