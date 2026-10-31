# Repository About and topics

Paste these into the repository settings (the gear icon beside **About** on the
repo's front page).

## Description

The About field is capped at 350 characters, and it is what appears in GitHub
search results and in Google's snippet. Lead with the capability, not the
method name.

**Suggested:**

```
Genetic sex estimation from BAM files via the X-to-autosome read ratio (Rx). Species-agnostic, needs no Y chromosome, works at low coverage. Includes threshold calibration by simulation.
```

A shorter alternative, if the long one reads as cluttered:

```
Genetic sex estimation from BAM files via the X-to-autosome read ratio. Species-agnostic, needs no Y chromosome.
```

## Topics

GitHub's internal search weights these heavily, and they cost nothing. Topics
must be lowercase and hyphenated.

```
bioinformatics
genomics
bam
sex-determination
karyotype
non-model-organisms
low-coverage-sequencing
conservation-genetics
population-genetics
pysam
```

## Website field

Leave empty until the methods paper is out, then point it at the DOI.

## Why MIT

`LICENSE` is MIT. For a small scientific tool meant to be picked up and cited,
permissiveness removes friction: anyone can vendor it into a pipeline, a
commercial lab can use it, and nothing propagates to their code. The citation
request in `CITATION.cff` is the mechanism for getting credit, not the licence.

GPL-3.0 would be the alternative, and it is what SexInference uses. It would
oblige anyone distributing a modified version to release their changes. That
matters for a large codebase one wants to keep open; for a single script, it
mostly deters casual reuse without protecting much.

The copyright line names you personally. If your institution claims copyright
on software produced there, change it to the institution — worth checking
before the paper goes out, since it is awkward to alter afterwards.

## CITATION.cff

GitHub renders this as a "Cite this repository" button that produces APA and
BibTeX automatically. Two things to fill in:

* **ORCID** — the file carries a placeholder. Replace it with yours, or delete
  the `orcid:` line.
* **Co-authors** — only you are listed. Add anyone who should share credit for
  the tool, which may differ from the manuscript author list.
* **Version and date** — once you tag a release, add `version:` and
  `date-released:` so the citation is specific.

When the methods paper appears, add a `preferred-citation:` block pointing at
it, so the button yields the paper rather than the repository.

## Before the paper goes out

* Tag a release (`v2.0.0`) so there is something citable and immutable.
* Consider archiving it on Zenodo, which mints a DOI per release and connects
  to GitHub automatically. Reviewers increasingly expect a DOI rather than a
  bare repository link, and it guarantees the version used in the paper stays
  available.
