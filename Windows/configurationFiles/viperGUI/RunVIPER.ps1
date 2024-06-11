wsl --cd $HOME -d VIPER -e bash -c @'
 exec bash --rcfile <(echo 'source ~/.bashrc && conda activate /home/viper/micromamba/envs/VIPERGenomeAssembler && python ~/viperGUI/viperGUI.py && exit')
'@