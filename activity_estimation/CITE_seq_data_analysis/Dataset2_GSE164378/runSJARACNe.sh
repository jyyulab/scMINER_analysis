dir=/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/18_Hao/02_sjaracneDownSample1000_newTFSIGs

## B
sjaracne lsf -e $dir/B_17336_17336_1000/B_17336_17336_1000.exp -g $dir/B_17336_17336_1000/sig/B_6647_6647_1000_sig.txt -o $dir/B_17336_17336_1000/sig_final -n 100 -pc 1e-2 -j ./config_cwlexec.json
sjaracne lsf -e $dir/B_17336_17336_1000/B_17336_17336_1000.exp -g $dir/B_17336_17336_1000/tf/B_1512_1512_1000_tf.txt -o $dir/B_17336_17336_1000/tf_final -n 100 -pc 1e-2 -j ./config_cwlexec.json

## CD4T
sjaracne lsf -e $dir/CD4_T_16458_16458_1000/CD4_T_16458_16458_1000.exp -g $dir/CD4_T_16458_16458_1000/sig/CD4_T_6491_6491_1000_sig.txt -o $dir/CD4_T_16458_16458_1000/sig_final -n 100 -pc 1e-2 -j ./config_cwlexec.json
sjaracne lsf -e $dir/CD4_T_16458_16458_1000/CD4_T_16458_16458_1000.exp -g $dir/CD4_T_16458_16458_1000/tf/CD4_T_1473_1473_1000_tf.txt -o $dir/CD4_T_16458_16458_1000/tf_final -n 100 -pc 1e-2 -j ./config_cwlexec.json

## CD8T
sjaracne lsf -e $dir/CD8_T_16575_16575_1000/CD8_T_16575_16575_1000.exp -g $dir/CD8_T_16575_16575_1000/sig/CD8_T_6487_6487_1000_sig.txt -o $dir/CD8_T_16575_16575_1000/sig_final -n 100 -pc 1e-2 -j ./config_cwlexec.json
sjaracne lsf -e $dir/CD8_T_16575_16575_1000/CD8_T_16575_16575_1000.exp -g $dir/CD8_T_16575_16575_1000/tf/CD8_T_1496_1496_1000_tf.txt -o $dir/CD8_T_16575_16575_1000/tf_final -n 100 -pc 1e-2 -j ./config_cwlexec.json

## DC
sjaracne lsf -e $dir/DC_18578_18578_1000/DC_18578_18578_1000.exp -g $dir/DC_18578_18578_1000/sig/DC_6949_6949_1000_sig.txt -o $dir/DC_18578_18578_1000/sig_final -n 100 -pc 1e-2 -j ./config_cwlexec.json
sjaracne lsf -e $dir/DC_18578_18578_1000/DC_18578_18578_1000.exp -g $dir/DC_18578_18578_1000/tf/DC_1565_1565_1000_tf.txt -o $dir/DC_18578_18578_1000/tf_final -n 100 -pc 1e-2 -j ./config_cwlexec.json

## Mono
sjaracne lsf -e $dir/Mono_17554_17554_1000/Mono_17554_17554_1000.exp -g $dir/Mono_17554_17554_1000/sig/Mono_6717_6717_1000_sig.txt -o $dir/Mono_17554_17554_1000/sig_final -n 100 -pc 1e-2 -j ./config_cwlexec.json
sjaracne lsf -e $dir/Mono_17554_17554_1000/Mono_17554_17554_1000.exp -g $dir/Mono_17554_17554_1000/tf/Mono_1510_1510_1000_tf.txt -o $dir/Mono_17554_17554_1000/tf_final -n 100 -pc 1e-2 -j ./config_cwlexec.json

## NK
sjaracne lsf -e $dir/NK_16509_16509_1000/NK_16509_16509_1000.exp -g $dir/NK_16509_16509_1000/sig/NK_6473_6473_1000_sig.txt -o $dir/NK_16509_16509_1000/sig_final -n 100 -pc 1e-2 -j ./config_cwlexec.json
sjaracne lsf -e $dir/NK_16509_16509_1000/NK_16509_16509_1000.exp -g $dir/NK_16509_16509_1000/tf/NK_1498_1498_1000_tf.txt -o $dir/NK_16509_16509_1000/tf_final -n 100 -pc 1e-2 -j ./config_cwlexec.json

