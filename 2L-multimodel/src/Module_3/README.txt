
One sample command to run the code:
    python infer.py config.yaml


required input file:
    instance_file
    fringe_tree_file
    original_dataset_file
    
    # the following files are model dependent
    ANN:
        {prop}_biases.txt
        {prop}_desc_norm_selected.txt
        {prop}_desc.csv
        {prop}_fringe.txt
        {prop}_values.txt
        {prop}_weights.txt
        {prop}.sdf
    LR:
        {prop}_desc_norm.csv
        {prop}_desc.csv
        {prop}_fringe.txt
        {prop}_values.txt
        {prop}_linreg.txt
        {prop}.sdf

