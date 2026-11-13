cd /tmp/claude-0/-home-claude/b3d49189-bf2d-54da-bcef-49acaa46ea63/scratchpad
M=/mnt/user-data/outputs/Ephys_master_table_DRAFT_plus_slice4.xlsx
C="minimlzip_x/miniML-main (includes sEPSC scripts)/core"
export TF_USE_LEGACY_KERAS=1
python3 miniml_batch_master.py --master "$M" --data-root /tmp/ephys_data --results res_s4 --miniml-core "$C" --model smoke/GC_lstm_model.h5 --only 20261112 > miniml_s4.log 2>&1
python3 miniml_batch_master.py --master "$M" --data-root /tmp/ephys_data --results res_cs_old --miniml-core "$C" --model smoke/GC_lstm_model.h5 --names "$(cat cs_names_old.txt)" > miniml_cs_old.log 2>&1
