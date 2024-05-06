#!/bin/bash

Rscript unit_test_data.R -i ges.rds -n1 test_net1.rds -n2 test_net1.rds
python unit_test.py -i ges.csv -n1 test_net1.csv -n2 test_net1.csv
Rscript unit_test.R -i ges.rds -n1 test_net1.rds -n2 test_net1.rds
