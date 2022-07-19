Guppy basecall example
```
nohup /opt/ont/guppy/bin/guppy_basecaller -i fast5_pass -s guppy_v_345/fast5_pass --num_callers 4 --config dna_r9.4.1_450bps_hac.cfg --device cuda:all:100% --qscore_filt
ering --min_qscore 7 >& debug.pass.txt&
```
