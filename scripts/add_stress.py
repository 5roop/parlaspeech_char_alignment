try:
    tg = snakemake.input.tg
    output = snakemake.output.tg
except NameError:
    tg = "data/TG/fgDDtwU4we4_2125.94-2134.66.TextGrid"
    output = "brisi.stressed.TextGrid"

import json
import polars as pl
from pathlib import Path
from praatio import textgrid
from parse import compile


def preprocess_text(s: str) -> str:
    """
    Lowercase the text, spruce up special characters and peculiarities, removes excess whitespace.
    :param str s: Input string
    :return str: Processed string
    """
    from string import punctuation

    for p in punctuation + "…“•":
        s = s.replace(p, "")
    s = s.casefold()
    s = (
        s.replace("x", "ks")
        .replace("y", "i")
        .replace("werthu", "vertu")
        .replace("ȃ", "a")
        .replace("î", "i")
        .replace("ȋ", "i")
        .replace("ȅ", "e")
        .replace("3208", "tri dva nula osam")
    )
    s = " ".join(s.split())
    return s


stress = json.loads(Path("parlaspeech.json").read_text())
stress = [i for i in stress if Path(i["audio"]).with_suffix("").name  == Path(tg).with_suffix("").name ]
tg = textgrid.openTextgrid(
    tg,
    includeEmptyIntervals=True,
)
mins, maxs = (
    tg.getTier("WordAlign").minTimestamp,
    tg.getTier("WordAlign").maxTimestamp,
)
intervals = []
for s in stress:
    text = s["text"]
    pp_text = preprocess_text(text).replace(" ", "")
    # print("Searching for:", pp_text)

    text_in_graphemes = "".join([g.label for g in tg.getTier("GraphAlign").entries])
    assert (pp_text in text_in_graphemes) or (pp_text == text_in_graphemes), "Mismatch!"
    L = len(tg.getTier("GraphAlign").entries)
    for ii in range(L + 1):
        for jj in range(L + 1):
            if ii == jj:
                continue
            if (
                "".join([g.label for g in tg.getTier("GraphAlign").entries[ii:jj]])
                == pp_text
            ):
                i_opt, j_opt = ii, jj
                break

    graphemes = tg.getTier("GraphAlign").entries[i_opt:j_opt]
    text_in_graphemes = "".join(
        [i.label for i in tg.getTier("GraphAlign").entries[i_opt:j_opt]]
    )
    assert text_in_graphemes == pp_text, "Bad code!"
    # print("Found", "".join([g.label for g in graphemes]))
    text_in_graphemes = "".join([i.label for i in graphemes])
    for i in s["stress_positive"]:
        # Which char is this?
        c = text[i].casefold()
        # How many times does c appear in the text?
        n = text[: i + 1].casefold().count(c)
        assert n > 0, "No hits found!"
        # Which grapheme corresponds to this in our TG?
        g = [g for g in graphemes if g.label == c][n - 1]
        intervals.append((round(g.start, 2), round(g.end, 2), "+"))
    for i in s["stress_negative"]:
        # Which char is this?
        c = text[i].casefold()
        # How many times does c appear in the text?
        n = text[: i + 1].casefold().count(c)
        assert n > 0, "No hits found!"
        # Which grapheme corresponds to this in our TG?
        g = [g for g in graphemes if g.label == c][n - 1]
        intervals.append((round(g.start, 2), round(g.end, 2), "-"))
intervals = sorted(set(intervals), key=lambda tpl: tpl[0])

tier = textgrid.IntervalTier("Stress", intervals, mins, maxs)
tg.addTier(tier)

tg.save(output, format="long_textgrid", includeBlankSpaces=True)
