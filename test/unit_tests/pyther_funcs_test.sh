#!/bin/sh
cd "$(dirname "$0")"
echo "------------------ BEGINNING TEST 1 ------------------"
echo "BASH running test_InteractomeToTable_func.R."
Rscript test_scripts/test_InteractomeToTable_func.R --net_obj=test_1/test_1_inputs/LNCaPWT_pruned.rds --out_name_project=LNCaPWT --out_dir=test_1/test_1_outputs/

echo "BASH running test_aREA_func.R."
Rscript test_scripts/test_aREA_func.R --ges=test_1/test_1_inputs/LNCaPWT_gExpr_GES.rds --network=test_1/test_1_inputs/LNCaPWT_pruned.rds --out_name_project=LNCaPWT --out_dir=test_1/test_1_outputs/

echo "BASH running GENERAL_convert_to_csv_forTesting.R."
Rscript test_scripts/GENERAL_convert_to_csv_forTesting.R --file=test_1/test_1_outputs/LNCaPWT_aREA_PAct.rds --ext=tsv --out_dir=test_1/test_1_outputs/

echo "BASH running test_Pyther.py."
python3 test_scripts/test_Pyther.py --ges=test_1/test_1_inputs/LNCaPWT_gExpr_GES.tsv --interactome=test_1/test_1_outputs/LNCaPWT_network.tsv --viper_tsv=test_1/test_1_outputs/LNCaPWT_aREA_PAct.tsv --out_name_project=LNCaPWT --out_dir=test_1/test_1_outputs/

echo "------------------ BEGINNING TEST 2 ------------------"
echo "BASH running test_InteractomeToTable_func.R."
Rscript test_scripts/test_InteractomeToTable_func.R --net_obj=test_2/test_2_inputs/LCRN1_pruned.rds --out_name_project=LCRN1 --out_dir=test_2/test_2_outputs/

echo "BASH running test_aREA_func.R."
Rscript test_scripts/test_aREA_func.R --ges=test_2/test_2_inputs/LCRN1_gExpr_GES.rds --network=test_2/test_2_inputs/LCRN1_pruned.rds --out_name_project=LCRN1 --out_dir=test_2/test_2_outputs/

echo "BASH running GENERAL_convert_to_csv_forTesting.R."
Rscript test_scripts/GENERAL_convert_to_csv_forTesting.R --file=test_2/test_2_outputs/LCRN1_aREA_PAct.rds --ext=tsv --out_dir=test_2/test_2_outputs/

echo "BASH running test_Pyther.py."
python3 test_scripts/test_Pyther.py --ges=test_2/test_2_inputs/LCRN1_gExpr_GES.tsv --interactome=test_2/test_2_outputs/LCRN1_network.tsv --viper_tsv=test_2/test_2_outputs/LCRN1_aREA_PAct.tsv --out_name_project=LCRN1 --out_dir=test_2/test_2_outputs/

echo "------------------ BEGINNING TEST 3 ------------------"
echo "BASH running test_InteractomeToTable_func.R."
Rscript test_scripts/test_InteractomeToTable_func.R --net_obj=test_3/test_3_inputs/AJ081_pruned.rds --out_name_project=AJ081 --out_dir=test_3/test_3_outputs/

echo "BASH running test_aREA_func.R."
Rscript test_scripts/test_aREA_func.R --ges=test_3/test_3_inputs/AJ081_gExpr_GES.rds --network=test_3/test_3_inputs/AJ081_pruned.rds --out_name_project=AJ081 --out_dir=test_3/test_3_outputs/

echo "BASH running GENERAL_convert_to_csv_forTesting.R."
Rscript test_scripts/GENERAL_convert_to_csv_forTesting.R --file=test_3/test_3_outputs/AJ081_aREA_PAct.rds --ext=tsv --out_dir=test_3/test_3_outputs/

echo "BASH running test_Pyther.py."
python3 test_scripts/test_Pyther.py --ges=test_3/test_3_inputs/AJ081_gExpr_GES.tsv --interactome=test_3/test_3_outputs/AJ081_network.tsv --viper_tsv=test_3/test_3_outputs/AJ081_aREA_PAct.tsv --out_name_project=AJ081 --out_dir=test_3/test_3_outputs/

echo "------------------ TESTING COMPLETE ------------------"