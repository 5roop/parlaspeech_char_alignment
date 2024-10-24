try:
    alignmentfile = snakemake.input.alignment
    wordalignmentfile = snakemake.input.wordalignment
    graphemealignment = snakemake.input.graphemealignment
    outfile = snakemake.output[0]
# except NameError:
#     alignmentfile = "kaldidebug/trans_phones_with_symbols.ctm"
#     wordalignmentfile = "kaldidebug/ali_words.ctm"
#     graphemealignment = "kaldidebug/L.offsets"
#     outfile = "brisi.TextGrid"
except NameError:
    alignmentfile = "tmp/MP_15_125.41-168.8/trans_phones_with_symbols.ctm"
    wordalignmentfile = "tmp/MP_15_125.41-168.8/ali_words.ctm"
    graphemealignment = "tmp/MP_15_125.41-168.8/L.offsets"
    outfile = "brisi.TextGrid"
import polars as pl

# import textgrids
from praatio import textgrid

phon_df = (
    pl.read_csv(
        alignmentfile,
        separator=" ",
        has_header=False,
        new_columns="sampleid channels start duration phoneme idk idk2".split(),
    )
    .with_columns((pl.col("start") + pl.col("duration")).alias("end"))
    .select(pl.exclude("^idk.*$"))
)


_, mins, maxs = (
    phon_df.group_by("sampleid")
    .agg(pl.col("start").min().alias("min"), pl.col("end").max().alias("max"))
    .to_numpy()
    .reshape(-1)
)

tg = textgrid.Textgrid()
intervals = []
for row in phon_df.iter_rows(named=True):
    intervals.append((round(row["start"], 2), round(row["end"], 2), row["phoneme"]))
tier = textgrid.IntervalTier("PhoneAlign", intervals, mins, maxs)
tg.addTier(tier)

# Words:
word_df = pl.read_csv(
    wordalignmentfile,
    separator=" ",
    has_header=False,
    new_columns="sampleid channels start duration word idk idk2".split(),
).with_columns(end=pl.col("start") + pl.col("duration"))

intervals = []
for row in word_df.iter_rows(named=True):
    try:
        if (previous_end := intervals[-1][1]) < row["start"]:
            intervals.append(
                (
                    round(previous_end, 2),
                    round(row["start"], 2),
                    "( )",
                )
            )
    except IndexError:
        pass
    intervals.append(
        (
            round(row["start"], 2),
            round(row["end"], 2),
            row["word"],
        )
    )
tier = textgrid.IntervalTier("WordAlign", intervals, mins, maxs)
tg.addTier(tier)

# Grapheme alignment:
import json
from pathlib import Path

offsets = json.loads(Path(graphemealignment).read_text())
intervals = []


def entries_overlap(
    this: tuple[float, float, str], other: tuple[float, float, str]
) -> bool:
    return (this[0] < other[1]) and (this[1] > other[0])


def percent_overlap(
    this: tuple[float, float, str], other: tuple[float, float, str]
) -> float:
    max_end = max(this[1], other[1])
    min_end = min(this[1], other[1])
    max_start = max(this[0], other[0])
    min_start = min(this[0], other[0])
    try:
        return (min_end - max_start) / (max_end - min_start)
    except ZeroDivisionError:
        return -1.0


for word in tg.getTier("WordAlign").entries:
    if word.label == "( )":
        continue
    # Find overlapping phonemes:
    phonemes = [
        p
        for p in tg.getTier("PhoneAlign").entries
        if entries_overlap(p, word) and (p.label != "sp")
    ]
    offset = [i for i in offsets if i["word"] == word.label][0]
    assert len(offset["pron"]) == len(phonemes), "Mismatch between phonemes and offset"
    for phoneme, offset in zip(phonemes, offset["offset"]):
        intervals.append(
            (phoneme.start, phoneme.end, word.label[offset[0] : offset[1]])
        )
# to_add = []
# for i, current in enumerate(intervals):
#     if i == 0: continue
#     if intervals[i-1].end != inter
tier = textgrid.IntervalTier("GraphAlign", intervals, round(mins, 2), round(maxs, 2))
tg.addTier(tier)
tg.save(outfile, format="long_textgrid", includeBlankSpaces=True)
