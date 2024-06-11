wsl --cd $HOME -d VIPER -e bash -c @'
 exec bash --rcfile <(echo 'source ~/.bashrc && conda activate /home/viper/micromamba/envs/VIPERPhylogeny && python ~/viperGUI/viperGUI_Phylo.py && exit')
'@