# Published order-23 class-recovery audit

Date: 2026-07-28

## Conclusion

No additional published order-23 frontier H-class was recovered.

The checked primary sources expose two different order-23 factor byte strings.
Both were already in the local corpus, both exact-verify at
`|det| = 2,779,447,296,000,000`, and both belong to the same H-class:

`db2cddf4b8f12da99a32b1563689ad693b89d0a9ae343a4aedde56361fdc5a81`

This does not contradict Orrick's 2006 statement that 14 distinct
record-attaining matrices had been found. The source for that paper preserves
the count but not those factors, starting states, random seeds, or search
code. An abstract class count was not treated as recovered matrix data.

## Sources checked

| source | retained order-23 material | result |
| --- | --- | --- |
| [Orrick et al. 2003 arXiv source](https://export.arxiv.org/e-print/math/0304410v1) | `matData.tex`, `records.tex` | Exactly one order-23 factor; it is the existing post-April reference. The search description has no seed or code. |
| [Orrick 2006 arXiv source](https://export.arxiv.org/e-print/math/0511141v2) | `equivalenceRevised.tex`, bibliography | Reports 14 matrices; contains no order-23 factor or reproducible generation material. |
| [Indiana d23 distinct-capture index](https://web.archive.org/cdx/search/cdx?url=www.indiana.edu%2F~maxdet%2Fd23.html&output=json&fl=timestamp%2Coriginal%2Cmimetype%2Cstatuscode%2Cdigest%2Clength&filter=statuscode%3A200&collapse=digest) | Two distinct HTML payloads | One factor in each payload; both already local and H-equivalent. |
| [Indiana maximal-determinant site-wide index](https://web.archive.org/cdx/search/cdx?url=http%3A%2F%2Fwww.indiana.edu%3A80%2F~maxdet%2F*&output=json&fl=timestamp%2Coriginal%2Cmimetype%2Cstatuscode%2Cdigest%2Clength&filter=statuscode%3A200&collapse=urlkey) | 191 unique URL keys | `d23.html` is the only order-23-specific entry; no ancillary factor or seed files. |
| [Brent maximal-determinant data](https://maths-people.anu.edu.au/~brent/maxdet/) and [archived URL index](https://web.archive.org/cdx/search/cdx?url=http%3A%2F%2Fmaths-people.anu.edu.au%2F~brent%2Fmaxdet%2F*&output=json&fl=timestamp%2Coriginal%2Cmimetype%2Cstatuscode%2Cdigest%2Clength&filter=statuscode%3A200&collapse=urlkey) | Current data for orders 19, 26, 27, 29, 33, and 37; 53 archived URL keys | No order-23 dataset. A direct `order23/` probe returned HTTP 404. |
| [Mendeley machine-readable collection](https://doi.org/10.17632/hzf94h43c5.1) | `NumM[23]=1`, one `Mat[23]` | Sole factor is token-identical to the existing post-April reference. |

The two Indiana payloads are:

- [pre-April capture](https://web.archive.org/web/20030308034054id_/http://www.indiana.edu:80/~maxdet/d23.html),
  CDX digest `S74B4VBUBY6HX6H7MNTK766GB7YFKRIV`
- [post-April capture](https://web.archive.org/web/20200219170713id_/http://www.indiana.edu/~maxdet/d23.html),
  CDX digest `KYCVLU6RJHDE574QVTP6Y4B2OSP3KZB4`

## Exact verification

| source factor | repository path | raw SHA-256 | normalized SHA-256 | receipt SHA-256 |
| --- | --- | --- | --- | --- |
| pre-April archive | `runs/direct-search/reference-data/orrick-pre-april2003.matrix.txt` | `b6a0801d0d1459ec262ae4e640c1b94b4baa3c44ea4b20aafbd1f8ae62f82be1` | `09a498bbef2d992ebf1fda856006fe20743993a4bc3f1998588d38790b87c202` | `4017a943b67e9120a4a539c1c711c955710b72f44fd715e0ab7d33c48b6845c1` |
| post-April / 2003 paper | `references/orrick-et-al-2003/matrix.txt` | `d134c240811076c9f807b98974ca68fda1e7756d1cf7e7ae72cdc044d8743850` | `c0ea58d361945b20dad78bddb3fd93c0810b762e527e8781e87d6db5e86d993a` | `45578d90f1d0e660ee50aa86d3aecd60c859627fdb1bec264024770ae32a32e5` |

Both passed `./arena verify`. Classification used
`research/h_equivalence_audit.py` with pinned `pynauty==2.8.8.1` and found:

- two unique matrix byte strings
- one H-class
- one H-plus-transpose class
- one normalized Gram class

## Local preservation

Downloaded source archives, Wayback indexes and page payloads, payload
checksums, and the complete classification JSON are under:

`runs/direct-search/reference-data/order23-published-class-recovery-20260728/`

That directory is intentionally under the repository's ignored `runs/` tree.
Its `README.md` records the URLs, hashes, extraction findings, and exact
blocker.

## Blocker

None of the checked primary archives retains bytes or reproducible generation
material for the other matrices behind the 14-matrix statement. Recovering
another class requires an actual 23 by 23 sign matrix, or deterministic code
and a seed that reproduce one, followed by exact verification and H-class
classification.
