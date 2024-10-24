tmpdir=kaldidebug/
mkdir -p "$tmpdir"
echo "speaker1 sampleid" >"$tmpdir"/spk2utt
echo "sampleid output/MP_15_125.41-168.8.wav" >"$tmpdir"/wavscp
echo "sampleid ma vi ste geograf" >"$tmpdir"/text

wavscp="$tmpdir"/wavscp
spk2utt="$tmpdir"/spk2utt
text="$tmpdir"/text

models=models
kaldi=/opt/kaldi

export LC_ALL=C

cut -f2- -d' ' $text | tr ' ' '\n' | sort -u >$tmpdir/wlist

/home/rupnik/anaconda3/envs/p37/bin/python lexicon.py $tmpdir/wlist $models/phonetisaurus-hr/model.fst $tmpdir
$kaldi/src/featbin/compute-mfcc-feats --config=$models/nnet3/conf/mfcc.conf scp:$wavscp ark:$tmpdir/mfcc.ark
$kaldi/src/online2bin/ivector-extract-online2 --config=$models/nnet3/conf/ivector.conf ark:$spk2utt ark:$tmpdir/mfcc.ark ark:$tmpdir/ivec.ark
$kaldi/egs/wsj/s5/utils/sym2int.pl -f 2- $tmpdir/words.txt $text >$tmpdir/text.int

$kaldi/src/bin/compile-train-graphs --read-disambig-syms=$tmpdir/disambig.int $models/nnet3/tdnn1a_sp/tree $models/nnet3/tdnn1a_sp/final.mdl $tmpdir/L.fst ark:$tmpdir/text.int ark:$tmpdir/graphs.fsts
$kaldi/src/nnet3bin/nnet3-latgen-faster --allow_partial=True --online-ivectors=ark:$tmpdir/ivec.ark --online-ivector-period=10 $models/nnet3/tdnn1a_sp/final.mdl ark:$tmpdir/graphs.fsts ark:$tmpdir/mfcc.ark ark:$tmpdir/ali.lat
$kaldi/src/latbin/lattice-align-words $tmpdir/word_boundary.int $models/nnet3/tdnn1a_sp/final.mdl ark:$tmpdir/ali.lat ark:- | $kaldi/src/latbin/lattice-to-ctm-conf ark:- - | $kaldi/egs/wsj/s5/utils/int2sym.pl -f 5 $tmpdir/words.txt - >$tmpdir/ali.ctm

# # # Phoneme alignment:
$kaldi/src/latbin/lattice-align-phones --replace-output-symbols=true $models/nnet3/tdnn1a_sp/final.mdl ark:$tmpdir/ali.lat ark:$tmpdir/phone_aligned.lats
$kaldi/src/latbin/lattice-to-ctm-conf --inv-acoustic-scale=10 --decode-mbr ark:$tmpdir/phone_aligned.lats $tmpdir/trans_phones.ctm
$kaldi/egs/wsj/s5/utils/int2sym.pl -f 5 $tmpdir/phones.txt $tmpdir/trans_phones.ctm >$tmpdir/trans_phones_with_symbols.ctm

cp $tmpdir/ali.ctm $tmpdir/ali_words.ctm
# # # # Grapheme alignment:
# echo "sampleid t a   p u t   s e   p o j a v i l a   l i s i c a   b o g   d a j   b o g   d a j   o v d i   s a n   p o d   j a b u k u n   k a   s i   t i   j a k o   s i   l i p a   j a   s a n   l i s i c a   h o d i   s i m o   d a   s e   i g r a m o   t a k o   s a n   ž a l o s t a n   n e   m o r e n   s e   s   t o b u n   i g r a t   a š   n i s a n   p i t o m a" >"$tmpdir"/text_for_graphemes

# text="$tmpdir"/text_for_graphemes
# cut -f2- -d' ' $text | tr ' ' '\n' | sort -u >$tmpdir/wlist
# /home/rupnik/anaconda3/envs/p37/bin/python lexicon.py $tmpdir/wlist $models/phonetisaurus-hr/model.fst $tmpdir
# $kaldi/src/featbin/compute-mfcc-feats --config=$models/nnet3/conf/mfcc.conf scp:$wavscp ark:$tmpdir/mfcc.ark
# $kaldi/src/online2bin/ivector-extract-online2 --config=$models/nnet3/conf/ivector.conf ark:$spk2utt ark:$tmpdir/mfcc.ark ark:$tmpdir/ivec.ark
# $kaldi/egs/wsj/s5/utils/sym2int.pl -f 2- $tmpdir/words.txt $text >$tmpdir/text.int

# $kaldi/src/bin/compile-train-graphs --read-disambig-syms=$tmpdir/disambig.int $models/nnet3/tdnn1a_sp/tree $models/nnet3/tdnn1a_sp/final.mdl $tmpdir/L.fst ark:$tmpdir/text.int ark:$tmpdir/graphs.fsts
# $kaldi/src/nnet3bin/nnet3-latgen-faster --allow_partial=True --online-ivectors=ark:$tmpdir/ivec.ark --online-ivector-period=10 $models/nnet3/tdnn1a_sp/final.mdl ark:$tmpdir/graphs.fsts ark:$tmpdir/mfcc.ark ark:$tmpdir/ali.lat
# $kaldi/src/latbin/lattice-align-words $tmpdir/word_boundary.int $models/nnet3/tdnn1a_sp/final.mdl ark:$tmpdir/ali.lat ark:- | $kaldi/src/latbin/lattice-to-ctm-conf ark:- - | $kaldi/egs/wsj/s5/utils/int2sym.pl -f 5 $tmpdir/words.txt - >$tmpdir/ali.ctm
# cp $tmpdir/ali.ctm $tmpdir/ali_graphemes.ctm
