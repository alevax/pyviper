# Instructions on how to conduct the unit tests for Pyther

### Updating the permissions of the shell script

If you are using a Mac, to update your permissions so that the shell script, `pyther_funcs_test.sh`, that conducts the tests can run, you may need to execute the following lines of code.
```
chmod u+x pyther/test/unit_tests/pyther_funcs_test.sh
xattr -d com.apple.quarantine pyther/test/unit_tests/pyther_funcs_test.sh
```

### Required Python & R Libraries

Python
* sys
* os
* optparse
* datetime
* math

R
* optparse
* tidyverse
* stringr

### Running the test

To execute the test after updating the permissions of the shell script, simply execute the shell script as follow.
```
pyther/test/unit_tests/pyther_funcs_test.sh
```

Three tests will be conducted using 3 separate sets of input files located in
```
pyther/test/unit_tests/test_1
pyther/test/unit_tests/test_2
pyther/test/unit_tests/test_3
```

Outputs containing results of the testing are stored in the subdirectories of these folders. The key outputs are the three `Pyther_VIPER_corr.txt` files with correlations between the results of Python's Pyther and the results of R's VIPER. If the functions in `Pyther/libs` are working correctly, these values should all be greater than 0.99 for all 3 tests.


### Details on the tests

#### Inputs

We will be using *test_1* code as an example. The other 2 tests are contain similar files and executions. Within the executed commands below, all the paths are relative to the `unit_tests` folder, as that is where `pyther_funcs_test.sh` is located.

Within `test_1_inputs`, there are 3 starting files.
1. `LNCaPWT_gExpr_GES.rds` - A GES (10000 genes, 250 samples)
2. `LNCaPWT_gExpr_GES.tsv` - the TSV version of 1
3. `LNCaPWT_pruned.rds` - an ARACNe-AP network created from 250 metacells made from the complete LNCaPWT_gExpr (21515 genes, 2450 samples)

##### Script 1

Script 1 tests `pyther/libs/InteractomeToTable.R`. The run command executed by `pyther_funcs_test.sh` is shown below.
```
Rscript test_scripts/test_InteractomeToTable_func.R
--net_obj=test_1/test_1_inputs/LNCaPWT_pruned.rds
--out_name_project=LNCaPWT
--out_dir=test_1/test_1_outputs/
```

The output of Script 1 is `LNCaPWT_network.tsv`, which is the ARACNe-AP network `LNCaPWT_pruned.rds` saved as a TSV file.

##### Script 2

Script 2 tests `pyther/libs/area_fn.R`. The run command executed by `pyther_funcs_test.sh` is shown below.

```
Rscript test_scripts/test_aREA_func.R
--ges=test_1/test_1_inputs/LNCaPWT_gExpr_GES.rds
--network=test_1/test_1_inputs/LNCaPWT_pruned.rds
--out_name_project=LNCaPWT
--out_dir=test_1/test_1_outputs/
```

The output of Script 2 is `LNCaPWT_aREA_PAct.rds`, the VIPER matrix.

##### Script 3

Script 3 converts the RDS output of `test_aREA_func.R` into a TSV file that can be used by `test_Pyther.py`. The run command executed by `pyther_funcs_test.sh` is shown below.

```
Rscript test_scripts/GENERAL_convert_to_csv_forTesting.R
--file=test_1/test_1_outputs/LNCaPWT_aREA_PAct.rds
--ext=tsv
--out_dir=test_1/test_1_outputs/
```

The output of Script 3 is `LNCaPWT_aREA_PAct.tsv`.

##### Script 4

Script 4 test `pyther/libs/pyther_fn.py`, in other words, Pyther. The run command executed by `pyther_funcs_test.sh` is shown below.

```
/Users/AlexanderWang/opt/miniconda3/bin/python3.8 test_scripts/test_Pyther.py --ges=test_1/test_1_inputs/LNCaPWT_gExpr_GES.tsv
--interactome=test_1/test_1_outputs/LNCaPWT_network.tsv
--viper_tsv=test_1/test_1_outputs/LNCaPWT_aREA_PAct.tsv
--out_name_project=LNCaPWT
--out_dir=test_1/test_1_outputs/
```

The output of this file is `LNCaPWT_Pyther_VIPER_corr.txt`, one of the `Pyther_VIPER_corr.txt` files described above. A correlation of at least 0.99 indicates Pyther is working properly.
