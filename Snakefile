from pathlib import Path
import json
import polars as pl


def preprocess_text(s: str) -> str:
    from string import punctuation

    for p in punctuation + "…":
        s = s.replace(p, "")
    s = s.casefold()
    s = s.replace("3208", "tri dva nula osam")
    s = " ".join(s.split())
    return s


data = json.loads(Path("parlaspeech.json").read_text())

expected_audios = [
    f"""output/wav/{Path(i["audio"]).with_suffix(".wav").name }"""
    for i in data]
expected_textgrids = [
    f"""output/stress/{Path(i["audio"]).with_suffix(".TextGrid").name }"""
    for i in data]


rule gather:
    input: expected_textgrids

rule fix_audio:
    input: "data/flac/{file}.flac"
    output: "output/wav/{file}.wav"
    shell:
        "ffmpeg -i {input[0]} -ac 1 -ar 16000 {output[0]}"

rule prep_files:
    input: "output/wav/{file}.wav"
    output:
        s = temp("tmp/{file}.spk2utt"),
        t = temp("tmp/{file}.text"),
        w = temp("tmp/{file}.wav.scp"),
        g = temp("tmp/{file}.text_for_graphemes")
    run:
        for entry in data:
            if Path(entry["audio"]).with_suffix("").name == Path(input[0]).with_suffix("").name:
                Path(output.s).write_text("speaker1 sampleid")
                Path(output.t).write_text("sampleid "+preprocess_text(entry["text"]))
                Path(output.g).write_text("sampleid "+ " ".join(preprocess_text(entry["text"])).replace("   ", " "))
                Path(output.w).write_text("sampleid "+input[0])
                break

rule do_kaldi:
    input:
        spk2utt="tmp/{file}.spk2utt",
        text="tmp/{file}.text",
        wavscp="tmp/{file}.wav.scp",
        text_for_graphemes="tmp/{file}.text_for_graphemes",
    output:
        temp(directory("tmp/{file}")),
        temp("tmp/{file}/trans_phones_with_symbols.ctm"),
        temp("tmp/{file}/ali_words.ctm"),
        temp("tmp/{file}/L.offsets")
    shell:
        """
        wavscp={input.wavscp}
        spk2utt={input.spk2utt}
        text={input.text}
        tmpdir={output[0]}

        models=models
        kaldi=/opt/kaldi

        export LC_ALL=C

        mkdir -p $tmpdir

        cut -f2- -d' ' $text | tr ' ' '\n' | sort -u > $tmpdir/wlist

        /home/rupnik/anaconda3/envs/p37/bin/python lexicon.py $tmpdir/wlist $models/phonetisaurus-hr/model.fst $tmpdir
        $kaldi/src/featbin/compute-mfcc-feats  --config=$models/nnet3/conf/mfcc.conf scp:$wavscp ark:$tmpdir/mfcc.ark
        $kaldi/src/online2bin/ivector-extract-online2 --config=$models/nnet3/conf/ivector.conf ark:$spk2utt ark:$tmpdir/mfcc.ark ark:$tmpdir/ivec.ark
        $kaldi/egs/wsj/s5/utils/sym2int.pl -f 2- $tmpdir/words.txt $text  > $tmpdir/text.int

        $kaldi/src/bin/compile-train-graphs --read-disambig-syms=$tmpdir/disambig.int $models/nnet3/tdnn1a_sp/tree $models/nnet3/tdnn1a_sp/final.mdl $tmpdir/L.fst ark:$tmpdir/text.int ark:$tmpdir/graphs.fsts
        $kaldi/src/nnet3bin/nnet3-latgen-faster --allow_partial=True --online-ivectors=ark:$tmpdir/ivec.ark --online-ivector-period=10 $models/nnet3/tdnn1a_sp/final.mdl ark:$tmpdir/graphs.fsts ark:$tmpdir/mfcc.ark ark:$tmpdir/ali.lat
        $kaldi/src/latbin/lattice-align-words $tmpdir/word_boundary.int $models/nnet3/tdnn1a_sp/final.mdl ark:$tmpdir/ali.lat ark:- | $kaldi/src/latbin/lattice-to-ctm-conf ark:- - | $kaldi/egs/wsj/s5/utils/int2sym.pl -f 5 $tmpdir/words.txt - > $tmpdir/ali.ctm
        cp $tmpdir/ali.ctm $tmpdir/ali_words.ctm

        # # # My additions:
        $kaldi/src/latbin/lattice-align-phones --replace-output-symbols=true $models/nnet3/tdnn1a_sp/final.mdl ark:$tmpdir/ali.lat ark:$tmpdir/phone_aligned.lats
        $kaldi/src/latbin/lattice-to-ctm-conf --inv-acoustic-scale=10 --decode-mbr ark:$tmpdir/phone_aligned.lats $tmpdir/trans_phones.ctm
        $kaldi/egs/wsj/s5/utils/int2sym.pl -f 5 $tmpdir/phones.txt $tmpdir/trans_phones.ctm > $tmpdir/trans_phones_with_symbols.ctm

        touch $tmpdir/ali_graphemes.ctm
        """

rule do_tg_compilation:
    input:
        alignment = "tmp/{file}/trans_phones_with_symbols.ctm",
        wordalignment = "tmp/{file}/ali_words.ctm",
        graphemealignment = "tmp/{file}/L.offsets",
        thedir = "tmp/{file}"
    output: "output/TG/{file}.TextGrid"
    conda:
        "miciprincalign.yml"
    script:
        "scripts/to_textgrid.py"
rule add_stress:
    input:
        tg="output/TG/{file}.TextGrid",
    output:
        tg="output/stress/{file}.TextGrid"
    conda:
        "miciprincalign.yml"
    script:
        "scripts/add_stress.py"