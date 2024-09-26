dir=/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/02_sjaracne

# B
sjaracne lsf -e $dir/B_13807_13807_304/B_13807_13807_304.exp -g $dir/B_13807_13807_304/sig/B_5856_5856_304_sig.txt -o $dir/B_13807_13807_304/sig_final -n 100 -pc 1e-2 -j ./config_cwlexec.json
sjaracne lsf -e $dir/B_13807_13807_304/B_13807_13807_304.exp -g $dir/B_13807_13807_304/tf/B_1355_1355_304_tf.txt -o $dir/B_13807_13807_304/tf_final -n 100 -pc 1e-2 -j ./config_cwlexec.json

# CD4T
sjaracne lsf -e $dir/CD4T_17399_17399_2878/CD4T_17399_17399_2878.exp -g $dir/CD4T_17399_17399_2878/sig/CD4T_6609_6609_2878_sig.txt -o $dir/CD4T_17399_17399_2878/sig_final -n 100 -pc 1e-2 -j ./config_cwlexec.json
sjaracne lsf -e $dir/CD4T_17399_17399_2878/CD4T_17399_17399_2878.exp -g $dir/CD4T_17399_17399_2878/tf/CD4T_1482_1482_2878_tf.txt -o $dir/CD4T_17399_17399_2878/tf_final -n 100 -pc 1e-2 -j ./config_cwlexec.json

# CD8T
sjaracne lsf -e $dir/CD8T_14501_14501_539//CD8T_14501_14501_539.exp -g $dir/CD8T_14501_14501_539//sig/CD8T_5944_5944_539_sig.txt -o $dir/CD8T_14501_14501_539/sig_final -n 100 -pc 1e-2 -j ./config_cwlexec.json
sjaracne lsf -e $dir/CD8T_14501_14501_539//CD8T_14501_14501_539.exp -g $dir/CD8T_14501_14501_539//tf/CD8T_1381_1381_539_tf.txt -o $dir/CD8T_14501_14501_539/tf_final -n 100 -pc 1e-2 -j ./config_cwlexec.json

# CM
sjaracne lsf -e $dir/CM_12484_12484_127/CM_12484_12484_127.exp -g $dir/CM_12484_12484_127/sig/CM_5577_5577_127_sig.txt -o $dir/CM_12484_12484_127/sig_final -n 100 -pc 1e-2 -j ./config_cwlexec.json
sjaracne lsf -e $dir/CM_12484_12484_127/CM_12484_12484_127.exp -g $dir/CM_12484_12484_127/tf/CM_1254_1254_127_tf.txt -o $dir/CM_12484_12484_127/tf_final -n 100 -pc 1e-2 -j ./config_cwlexec.json

# NK
sjaracne lsf -e $dir/NK_15039_15039_754/NK_15039_15039_754.exp -g $dir/NK_15039_15039_754/sig/NK_6160_6160_754_sig.txt -o $dir/NK_15039_15039_754/sig_final -n 100 -pc 1e-2 -j ./config_cwlexec.json
sjaracne lsf -e $dir/NK_15039_15039_754/NK_15039_15039_754.exp -g $dir/NK_15039_15039_754/tf/NK_1429_1429_754_tf.txt -o $dir/NK_15039_15039_754/tf_final -n 100 -pc 1e-2 -j ./config_cwlexec.json

