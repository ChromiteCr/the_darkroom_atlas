# D1 datasheet inventory: Kodak Tri-X 400, T-Max 100, D-76, XTOL, HC-110 and T-MAX Developer

Research date: 2026-11-20. Every Kodak PDF named "opened" below was fetched, its text layer was extracted, and the curve pages were rendered and inspected. The curve figures in the 2016–2021 editions are **vector paths** inside the PDF, not raster images (see section 6). Numbers quoted from Kodak tables are limited to the rows D1 needs.

Consistency pass 2026-11-21 (re-checked against the extracted text of F-4016 2018, F-4017 2021, J-109 2018, J-24 2017, J-78 2017): corrected the reciprocity statement in 5.1 (Tri-X needs +10 % development at 1/1000 s), the T-Max 100 CI time axis in section 6 (4–18 min, not 0–18), the F-4017 figure page range (pp. 8–10), the J-109 Table 1 CI-aim coverage per film in section 0, and marked the T-Max 100 D-76 30-s agitation as coming from the F-4016 procedure rather than the curve caption; added the 24 °C reference of the T-MAX adjustment factors.

Independent spot-check 2026-11-21 (second reviewer). Re-fetched from kodakprofessional.com today: F-4017 (10-21, MD5 94fd0bb9…, 788 169 B), F-4016 (6-18, MD5 4a193fb1…, 547 354 B), J-109 (02/2018), J-24 (12-17) and J-86 (12-17). Re-extracted their text and re-digitised the CI vector paths myself. **Confirmed against those files:** edition dates and page counts; CI 0.56 aim (F-4017 p. 3, J-86 p. 2); all characteristic-curve captions in 2.1 and 3.1; the 20 °C small-tank rows for both films; the reciprocity tables; RMS granularity 17 and 8 with their read conditions; spectral captions (0.5 s and 1.4 s); J-109 Table 1 rows; HC-110 B = stock 1:7 = concentrate 1:31 and the J-24 dilution-B roll rows; the absence of curves in J-109, J-24 and J-86 (the only long paths are small page icons). Also confirmed: every CI span in tables 2.1 and 3.1 (to within 0.01 min / 0.001 CI), including the large-tank legend mapping. **ISO 6:1993** was re-read from a web.archive.org copy of the preview, and its status was read from iso.org. The Contrast Index citation was confirmed through a patent that cites it (see 5.2). **Checked only against the cached text extractions** (publisher copies not re-fetched): F-9 "Daylight, 1/50 second", the F-4017 2007 old-emulsion notes, the T-MAX 120 4.7-mil base (F-32 2002, F-4016 2007) and the J-78 legacy rows. **Not re-derived:** the speed check in 5.3 and the base + fog proxies, which remain the original author's estimates. **Unverified:** the F-4016 byte-identity statement in 1.1 and the arc geometry of CI (5.2). One caveat was added in 2.1, on the Tri-X spectral-exposure provenance.

Confidence labels: **opened** = I read the source myself. **secondary** = taken from a search summary or index, not from the document. **estimate** = my own derivation from opened data. **recalled** = from memory of the literature, not opened.

---

## 0. Answer in brief

* **Planned split as written is not feasible at the characteristic-curve level.** Kodak publishes no D–log H characteristic curve for Tri-X 400 or T-Max 100 in XTOL or in HC-110, in any edition I found (2002–2021). The only developers with characteristic curves for both films are **D-76** and **T-MAX Developer**. T-Max 100 also has T-MAX RS curves.
* **What exists for the XTOL and HC-110 pairs:**
  * Contrast-index-versus-time curves (CI(t)).
  * Time–temperature tables (all aimed at CI 0.56).
  * Push tables.
  * For XTOL only, J-109 Table 1, which gives small-tank development times for several **CI aims** at full strength (18–24 °C; no 27 °C entry for either film) and 1:1 (20–27 °C): 0.52/0.56/0.62/0.72/0.82 for T-Max 100, and 0.56/0.72/0.82 only for Tri-X 400 (no 27 °C 1:1 entry for Tri-X).
* **Consequence:** the split is feasible only if the comparison metric for the XTOL and HC-110 pairs is CI(t) and time-to-CI, not curve shape, toe or speed.
* **Alternative:** the nearest split that keeps a full characteristic-curve held-out test holds out **T-Max 100 × T-MAX Developer**, with Tri-X × T-MAX Developer in training. The caveat is that the T-MAX Developer formula is proprietary (msds-estimated).
* **Agitation mismatch:** Tri-X D-76 curves are large tank (caption: agitation at 1-min intervals), while T-Max 100 D-76 curves are small tank (caption says only "Small Tank"; the 30-s interval comes from F-4016's small-tank procedure, p. 2). "Film parameters from D-76" therefore spans two agitation regimes.

---

## 1. Publication inventory

### 1.1 Films

| Film | Pub. no. | Exact title (as on cover) | Date / revision | URL | Status |
|---|---|---|---|---|---|
| Tri-X 400 (400TX) + 320TXP | F-4017 | KODAK PROFESSIONAL TRI-X 320 and 400 Films (TECHNICAL DATA / BLACK-AND-WHITE FILM) | **October 2021**, back cover "Revised 10-21", © 2021 Kodak Alaris; 11 pp | https://kodakprofessional.com/sites/default/files/wysiwyg/film/f4017_trix_320400.pdf | opened; **newest found** |
| same | F-4017 | same | December 2016, "Revised 12-16"; 9 pp | https://125px.com/docs/film/kodak/f4017-2016.pdf (mirror). Also listed at https://cdn.shopify.com/s/files/1/0339/5113/files/f4017_trix_320400.pdf (not opened) | opened |
| same | F-4017 | same | February 2016, "Revised 2-16"; 9 pp | https://kodakprofessional.com/sites/default/files/wysiwyg/pro/resources/f4017_TriX.pdf (this is the file the official Resources page still links) | opened |
| same | F-4017 | KODAK PROFESSIONAL TRI-X 320 and 400 Films | May 2007, © Eastman Kodak 2007; 16 pp | https://125px.com/docs/film/kodak/f4017-400TX-2007.pdf | opened |
| same | F-4017 | (first 400TX edition) | 11/2005 | https://125px.com/docs/film/kodak/f4017-400TX-2005.pdf | listed in the 125px index, not opened |
| Tri-X Pan (TX, pre-2007 emulsion "5063") | F-9 | KODAK TRI-X Pan and KODAK TRI-X Pan Professional Films | April 2003, "Minor Revision 4-03"; 12 pp | https://125px.com/docs/film/kodak/f9-Tri-X_Pan.pdf | opened; **legacy emulsion** |
| T-Max 100 (100TMX) | F-4016 | KODAK PROFESSIONAL T-MAX 100 Film | **June 2018**, "Revised 6-18", © 2018 Kodak Alaris on p. 1 (the back page still carries © 2016); 9 pp | https://kodakprofessional.com/sites/default/files/wysiwyg/pro/resources/f4016_TMax_100.pdf (original author reported it byte-identical to https://125px.com/docs/film/kodak/f4016_tmax_100-2018.pdf. **Unverified on spot-check:** the file served on 2026-11-21 is 547 354 B, MD5 4a193fb16cbf43cd0e6d40d03c7ccd89, created 2018-06-05. A cached copy from an earlier fetch is 542 705 B with a different MD5. Text, tables and CI vectors are the same in both) | opened; **newest found** |
| same | F-4016 | KODAK PROFESSIONAL T-MAX 100 Film | 02/2016 | https://125px.com/docs/film/kodak/f4016_TMax_100-2016.pdf | listed, not opened |
| T-MAX 100/400/P3200 (combined) | F-4016 | KODAK PROFESSIONAL T-MAX Films | May 2007, © Eastman Kodak 2007; 30 pp | https://125px.com/docs/film/kodak/f4016-TMAX-2007a.pdf | opened |
| same | F-4016 | KODAK PROFESSIONAL T-MAX Films | October 2007 and 02/2004 editions | https://125px.com/docs/film/kodak/f4016-TMAX-2007b.pdf, https://125px.com/docs/film/kodak/f4016-TMAX-2004.pdf | listed, not opened |
| T-MAX Professional (TMX, previous packaging) | F-32 | KODAK T-MAX Professional Films | March 2002; 28 pp | https://125px.com/docs/film/kodak/f32-TMAX.pdf | opened; legacy |
| T-MAX 400 (400TMY) | F-4043 | KODAK PROFESSIONAL T-MAX 400 Film | October 2007, "New 10-07", © Eastman Kodak 2007; 12 pp | https://www.freestylephoto.com/pdf/product_pdfs/kodak/Kodak_Tmax400.pdf (dealer mirror) | opened. F-4043 is the T-MAX 400 sheet, **not** T-Max 100, so it is not relevant to D1 |

Product pages, which have no PDF links of their own:
* https://kodakprofessional.com/photographers/film/black-white/kodak-professional-tri-x-films/515
* https://kodakprofessional.com/photographers/film/black-white/kodak-professional-t-max-100-film/525

The official PDF links live on https://kodakprofessional.com/photographers/resources.

### 1.2 Developers

| Developer | Pub. no. | Exact title | Date / revision | URL | Status |
|---|---|---|---|---|---|
| D-76 | J-78 | KODAK Developer D-76 (TECHNICAL DATA / CHEMICALS) | December 2017, "Revised 12-17"; 7 pp | https://125px.com/docs/techpubs/kodak/j78-2017.pdf (mirror). Not on the kodakprofessional.com Resources page; …/pro/chemistry/j78.pdf returns 404 | opened |
| D-76 | J-78 | same | 10/2002 | https://125px.com/docs/techpubs/kodak/j78-2002_10.pdf | listed, not opened |
| XTOL | J-109 | KODAK PROFESSIONAL XTOL Developer | February 2018, "Revised 02/2018"; 16 pp | https://kodakprofessional.com/sites/default/files/wysiwyg/pro/chemistry/J-109_Feb_2018.pdf (linked from the Resources page; a retry worked and the file is byte-identical, same MD5, to https://125px.com/docs/techpubs/kodak/J-109_Feb_2018.pdf) | opened |
| XTOL | J-109 | KODAK XTOL Developer | April 1998; 24 pp | https://125px.com/docs/techpubs/kodak/j109-1998_04.pdf | opened |
| XTOL | J-109 | — | 09/2004, 03/2008 | 125px index (`j109-2004_09.pdf`, `j109-2008_03.pdf`) | listed, not opened |
| XTOL | J-107 | KODAK XTOL Developer for Small Tank and Tray Processing | September 1996 | https://125px.com/docs/techpubs/kodak/j107-1996_09.pdf | opened |
| XTOL | J-108 | KODAK XTOL Developer for Rotary-Tube Processors | 09/1996 | https://125px.com/docs/techpubs/kodak/j108-1996_09.pdf | listed, not opened |
| HC-110 | J-24 | KODAK PROFESSIONAL HC-110 Developer | December 2017, "Revised 12-17"; 8 pp | https://kodakprofessional.com/sites/default/files/wysiwyg/pro/chemistry/j24.pdf | opened |
| T-MAX / T-MAX RS | J-86 | KODAK PROFESSIONAL T-MAX Developers | December 2017; 13 pp | https://kodakprofessional.com/sites/default/files/wysiwyg/pro/chemistry/j86.pdf | opened |

Archive copies: the original author could not fetch web.archive.org, so no Wayback URL for the Kodak publications was verified. On the 2026-11-21 spot-check web.archive.org was reachable, and the ISO 6 preview was read through it (see 5.2). The Kodak Wayback URLs below were still not checked. Wikipedia's Xtol article cites http://web.archive.org/web/20131013044821/http://wwwuk.kodak.com/global/en/professional/support/techPubs/j109/j109.pdf. That link is unverified.

**None of the developer publications contain characteristic curves or CI curves.** I checked J-78 (2017), J-109 (1998 and 2018), J-107 (1996), J-24 (2017) and J-86 (2017). They contain only time tables, mixing data and process-control notes.

---

## 2. Tri-X 400 (400TX): F-4017, October 2021 (identical content in Feb 2016, Dec 2016 and, for 400TX, May 2007)

### 2.1 Figures (pp. 8–10 of the Oct 2021 edition; 400TX curves on pp. 8–9, 320TXP sheet curves on pp. 9–10)

| Figure | Film/format | Process (verbatim caption content) | Densitometry | Exposure stated? | Curves |
|---|---|---|---|---|---|
| Characteristic curves | 400TX 35 mm | KODAK PROFESSIONAL T-MAX Developer, Small Tank, agitation at 30-second intervals, 20 °C (68 °F) | Diffuse visual | **No** | 6, 7, 9, 11 min |
| Characteristic curves | 400TX 120 | same as above | Diffuse visual | No | 6, 7, 9, 11 min |
| Characteristic curves | 400TX 35 mm | KODAK PROFESSIONAL D-76 Developer, Large Tank, agitation at 1-minute intervals, 20 °C | Diffuse visual | No | 6, 8, 10, 12 min |
| Characteristic curves | 400TX 120 | same as above | Diffuse visual | No | 6, 8, 10, 12 min |
| Contrast Index curves | 400TX 35 mm and 120 | Small Tank at 20 °C | Diffuse visual | — | D-76, D-76 (1:1), T-MAX, T-MAX RS, XTOL, XTOL 1:1, HC-110 (B). Axes: 0–30 min, CI 0–1.2 |
| Characteristic curves (HC-110 B) | **320TXP sheets only** | HC-110 Dilution B, Large Tank 20 °C, 1-min agitation | Diffuse visual | No | 4, 6, 9, 14 min. **Not 400TX** |
| Characteristic curves (D-76) | 320TXP sheets | D-76, Large Tank 20 °C, 1-min | Diffuse visual | No | 6, 8, 10, 12 min |
| CI curves | 320TXP sheets | Large tank, tray, rotary tube (20 °C) | Diffuse visual | — | D-76, T-MAX RS, XTOL, HC-110 (B), and (tray/rotary) XTOL 1:1, D-76 1:1 |
| MTF | Tri-X (not split by format) | KODAK PROFESSIONAL Developer D-76, Large Tank | Diffuse Visual | No illuminant given. The 2007 edition caption literally reads "Exposure: X" | 1–600 cycles/mm log axis, response % (log) |
| Spectral sensitivity | Tri-X | not stated | Diffuse Visual | **Effective exposure 0.5 s** | Two curves: **D = 0.3 and D = 1.0 above gross fog**. Units: log sensitivity, where "sensitivity = reciprocal of exposure (erg/cm²) required to produce specified density"; 250–750 nm |

**CI curves for 400TX, digitised exactly from the vector paths (estimate), small tank, 20 °C.** Each curve is listed as time span → CI span:

| Developer | Time span (min) | CI span |
|---|---|---|
| D-76 | 5.4 → 10.8 | 0.48 → 0.82 |
| D-76 1:1 | 5.4 → 18.4 | 0.37 → 0.90 |
| XTOL | 3.25 → 14.1 | 0.34 → 0.98 |
| XTOL 1:1 | 4.3 → 17.3 | 0.37 → 0.90 |
| HC-110 (B) | 3.25 → 11.9 | 0.53 → 1.06 |
| T-MAX | 5.4 → 10.8 | 0.52 → 0.86 |
| T-MAX RS | 3.25 → 13.0 | 0.51 → 0.99 |

Legends were matched by dash pattern and line width. Linear interpolation reproduces the table times: CI 0.56 at about 6.75 min (D-76) and about 7 min (XTOL).

**Old-film provenance (important for tags).** The May 2007 F-4017, p. 10, states two things:
* The image-structure (granularity) data "was generated using the old versions of these films."
* The MTF and spectral sensitivity "were generated using the older version of this film (KODAK TRI-X Pan Film 5063)."

Later editions drop the note but show the same figures.

*Spot-check caveat (2026-11-21):* the spectral figure in F-4017 (2007 and 2021) gives "Effective Exposure: 0.5 sec." and "D=0.3/1.0 > Gross Fog". The legacy F-9 (2003, TX/5063) spectral figure gives "1.4 seconds", D-76 20 °C, "0.3/1.0 + D-min". If the F-4017 curves really are the 5063 data, as the 2007 note says, then the two captions do not reconcile. Tag the Tri-X spectral data `datasheet` with this ambiguity noted. Do not treat the 0.5 s as a measured exposure time of the current emulsion.

### 2.2 Tables

| Item | Content | Source |
|---|---|---|
| Processing aim | "starting-point recommendations are intended to produce a contrast index of 0.56" | F-4017 2021 p. 3 |
| Time table (normal) | Small tank (30-s agitation) and large tank (1-min) at 18/20/21/22/24 °C. Developers: T-MAX, T-MAX RS, HC-110 (B), D-76, D-76 1:1, XTOL, XTOL 1:1. Also rotary tube (400TX), plus tray, large tank and gaseous-burst tables (320TXP only) | p. 3–5 |
| 20 °C small-tank rows (400TX) | T-MAX 6; T-MAX RS 4½; HC-110 (B) 3¾; D-76 6¾; D-76 1:1 9¾; XTOL 7; XTOL 1:1 9 min | p. 3. The same values appear in the Feb 2016, Dec 2016 and May 2007 editions |
| 20 °C large-tank rows (400TX) | T-MAX NR; T-MAX RS 5; HC-110 (B) 4½; D-76 7¾; D-76 1:1 11; XTOL 8; XTOL 1:1 10½ min | p. 3 |
| Push tables | EI 1600 (2-stop) and EI 3200 (3-stop). Small tank, large tank and rotary tube, 18–24 °C, same developers | p. 6 |
| Rack-and-tank | EI 400 and 800 (one row): T-MAX RS 4–6, D-76 5–7, XTOL 5–7 min at 22 °C | p. 7 |
| Reciprocity | 1/100 000 s: +1 stop, dev +20 %. 1/10 000: +½ stop, +15 %. 1/1000: none, +10 %. 1/100 and 1/10: none. 1 s: +1 stop or 2 s, −10 %. 10 s: +2 stops or 50 s, −20 %. 100 s: +3 stops or 1200 s, −30 %. A graph covers 1–100 s | p. 2 |
| RMS granularity | **17 (fine)** for 400TX (16 for 320TXP). "Read at a net diffuse density of 1.0, using a 48-micrometre aperture, 12x magnification." Development: **HC-110 (Dilution B), 20 °C** | p. 7 |
| Exposure index | Not stated as "ISO 400" anywhere in the 2007–2021 F-4017. Implied by the push tables (EI 1600 = 2-stop), the rack-and-tank "EI 400", and the single-use camera spec "Film Speed 400". J-109 Table 1 lists 400TX at EI 400/800 → CI 0.56. Legacy F-9 states "ISO 400/27°" for TX | — |
| Base | **Not stated for 400TX** in any F-4017. Legacy TX (F-9, 2003): 135 = 5-mil acetate, 120 = 3.6-mil acetate. 320TXP sheets: 7-mil ESTAR Thick Base | F-9 p. 1 |
| Base + fog | **Not tabulated.** Lowest plotted density of the curves (at log H = −3.45): 35 mm D-76 0.27–0.32; 120 D-76 0.16–0.21; 35 mm T-MAX 0.31–0.38; 120 T-MAX 0.20–0.27 (estimate, from vector data). 35 mm runs about 0.1 higher than 120 in both developers | — |
| Filter factors | Wratten 8/11/12/15/25/47/58, polariser, for daylight and tungsten | p. 2 |

**Legacy F-9 (April 2003, TX).** Characteristic curves are "Exposure: Daylight, 1/50 second" for:
* D-76, large tank, 20 °C: 7, 9, 11 min.
* T-MAX Developer, small tank, 24 °C: 5, 7, 9, 11 min.

The CI curves include DK-50 (1:1), D-76, D-76 (1:1), HC-110 (B) and MICRODOL-X; they contain **no XTOL**. The MTF figure is tungsten exposure, D-76 large tank 20 °C. Spectral sensitivity: 1.4 s, D-76 20 °C, 0.3 and 1.0 + D-min. Granularity 17 in HC-110 (B), 7½ min small tank, to CI 0.56.

This is a different emulsion. Use it only as evidence of Kodak's sensitometric exposure convention.

---

## 3. T-Max 100 (100TMX): F-4016, June 2018 (curves identical in the May 2007 combined edition)

### 3.1 Figures (pp. 8–9)

| Figure | Process | Exposure | Densitometry | Curves |
|---|---|---|---|---|
| Characteristic curves | KODAK PROFESSIONAL Developer D-76, **Small Tank**, 20 °C (caption states no agitation; F-4016 p. 2 small-tank procedure: 5–7 inversions every 30 s) | **Daylight** (no time stated) | Diffuse Visual | 6, 7.5, 10 min (density axis 0–4) |
| Characteristic curves | KODAK PROFESSIONAL T-MAX RS Developer and Replenisher, **Large Tank**, 20 °C | Daylight | Diffuse Visual | 8, 10.5, 13, 15 min |
| Characteristic curves | KODAK PROFESSIONAL T-MAX Developer, Small Tank, 20 °C | Daylight | Diffuse Visual | 6, 7, 10, 12 min (density axis 0–3) |
| CI curves | Small Tank, 20 °C | — | Diffuse Visual | D-76, D-76 (1:1), T-MAX, T-MAX RS, XTOL, XTOL 1:1, HC-110 (B). Axes 4–18 min, CI 0.3–0.9 |
| CI curves | Large Tank, 20 °C | — | Diffuse Visual | D-76, T-MAX RS, XTOL, HC-110 (B) |
| CI curves | Rotary Tube, 20 °C | — | Diffuse Visual | XTOL, XTOL 1:1 |
| CI curves | Tray, 20 °C | — | Diffuse Visual | XTOL, XTOL 1:1 |
| MTF ("Modulation Transfer Curves") | KODAK Developer D-76, 20 °C, Small Tank | **Tungsten** | Diffuse visual | 1–600 cycles/mm |
| Spectral sensitivity | KODAK Developer D-76, 20 °C | **Effective exposure 1.4 seconds** | Diffuse visual | **0.3 and 1.0 greater than D-min**. Log sensitivity, sensitivity = 1/exposure in ergs/cm², 250–750 nm |

**CI curves for T-Max 100, digitised from vectors (estimate).** All curves span only **CI 0.56 → 0.82**:

| Tank | Developer | Time span (min) |
|---|---|---|
| Small tank | D-76 | 6.3 → 9.1 |
| Small tank | D-76 1:1 | 9.4 → 13.3 |
| Small tank | XTOL | 7.3 → 10.6 |
| Small tank | XTOL 1:1 | 9.7 → 13.4 |
| Small tank | HC-110 (B) | 5.9 → 14.3 |
| Small tank | T-MAX | 7.7 → 14.5 |
| Small tank | T-MAX RS | 7.8 → 14.6 |
| Large tank | D-76 | 7.1 → 10.3 |
| Large tank | XTOL | 8.2 → 11.9 |
| Large tank | HC-110 (B) | 6.6 → 16.1 |
| Large tank | T-MAX RS | 8.7 → 16.5 |

These curves have only 3–8 vertices, so they are probably drawn through a few test points.

The CI 0.56 end matches the 20 °C time table: D-76 6½, XTOL 7½, HC-110 6, D-76 1:1 9½, XTOL 1:1 9½ min. So F-4016's times also aim at CI 0.56, even though F-4016 never says so (F-4017 and J-86 do).

### 3.2 Tables

| Item | Content | Source |
|---|---|---|
| Processing aim | "contrast appropriate for printing with a diffusion enlarger" (no CI number in F-4016). J-86 states the aim is CI 0.56 for all films except T-MAX 400 (0.60) | F-4016 p. 2; J-86 p. 2 |
| 20 °C small-tank rows | T-MAX (1:4) 7½; T-MAX RS 8; XTOL 7½; XTOL 1:1 9½; D-76 6½; D-76 1:1 9½; HC-110 (B) 6 min. Full table covers 18/20/21/22/24 °C, and T-MAX 1:7 and 1:9 at 24 °C | p. 3 |
| 20 °C large-tank rolls | T-MAX 8½; T-MAX RS 8¾; XTOL 8¼; D-76 7¼; HC-110 (B) 6½ min. Sheets tables, tray, rotary-tube rolls and sheets tables are also given | p. 3–4 |
| Push | EI 200, 400 and 800 tables for small tank, large tank and rotary | p. 5–6 |
| Development-time adjustment factors | Factors for −20 %, normal, +20 % and +40 % contrast at 18–24 °C (T-MAX/T-MAX RS rows: 20/22/24 °C only), for T-MAX/T-MAX RS, D-76 and HC-110 (B). Example for D-76: 20 °C 1.0; 18 °C 1.2; 24 °C 0.7. Example for HC-110 (B): 20 °C 1.0; 18 °C 1.2; 24 °C 0.7. The T-MAX/T-MAX RS factors are referenced to **24 °C** (1.0 there; 20 °C 1.2), the others to 20 °C. This gives a temperature dependence and a CI-ratio dependence. The factors are only approximately consistent with the time tables (e.g. D-76 small tank 18 °C: 1.2 × 6½ = 7.8 min vs 7½ in the table; HC-110 (B): 1.2 × 6 = 7.2 vs 6½), so treat them as rounded guidance, not data | p. 7 |
| Reciprocity | 1/10 000 s: +⅓ stop. 1/1000, 1/100 and 1/10: none. 1 s: +⅓ stop. 10 s: +½ stop or 15 s. 100 s: +1 stop or 200 s. **No development adjustment column** | p. 2 |
| Resolving power | 63 lines/mm (TOC 1.6:1) and 200 lines/mm (TOC 1000:1), "method similar to ISO 6328" | p. 8 |
| RMS granularity | **8**. "Read at a net diffuse density of 1.00, using a 48-micrometre aperture, 12X magnification." Development: **D-76, 20 °C** | p. 8 |
| Exposure index | "nominal speed … is EI 100. It was determined in a manner published in ISO standards"; "ISO 100/21° in most developers" | p. 1 |
| Base | Not stated in the 2018 edition. The 2007 combined F-4016 and F-32 (2002) say the T-MAX **120** is "coated on a thicker (4.7-mil) base than other black-and-white roll films". **135 base: not found** | F-4016 2007 p. 2 |
| Base + fog | Not tabulated. Lowest plotted densities (at log H about −3.15): D-76 ST 0.21–0.22; T-MAX RS LT 0.22–0.26; T-MAX ST 0.24–0.30 (estimate, from vectors) | — |

Legacy F-32 (2002, TMX packaging) has characteristic curves, all with "Exposure: Daylight" and no exposure time:
* T-MAX Developer, small tank, 24 °C: 5, 7, 9, 11 min.
* D-76, small tank, 20 °C: 6, 8, 10, 12 min.
* T-MAX RS, large tank, 24 °C.
* VERSAMAT.

It has no XTOL or HC-110 characteristic curves.

---

## 4. Developer publications: what they add for D1

| Pub | Useful content | Caveat |
|---|---|---|
| J-78 D-76 (Dec 2017) | Agitation procedures; time tables for small tank and large tank at 18–24 °C, full strength and 1:1; push times; D-76R replenishment. **No formula, no curves** | The film rows use **legacy names only** ("T-MAX 100 Professional", "TRI-X Pan"). For example, T-MAX 100 Professional small tank 20 °C is 9 min, against 6½ min for 100TMX in F-4016. **Do not use J-78 times for 100TMX/400TX**; the film publications govern |
| J-109 XTOL (Feb 2018) | Working-tank specific gravity 1.085 ± 0.003 at 25 ± 0.3 °C, "at pH 8.2 ± 0.05" (p. 1, mixing step 5; checked against the rendered page 2026-11-21). Replenishment 70 mL per 135-36 roll. Table 1 (small tank) has a **CI column per EI** (rows below). Tables 2–8 cover rotary tube, large tank seasoned and fresh, and sheet trays, tubes and tanks | Two-part powder; no formula published |
| J-24 HC-110 (Dec 2017) | Dilution table: stock = 1 + 3 concentrate. **Dilution B = stock 1:7 = concentrate 1:31**. Roll and sheet times give both legacy and PROFESSIONAL film rows. 400TX dilution B, small tank, 18/20/21/22/24 °C: 4½ / 3¾ / 3½ / 3 / 2½ min. PROFESSIONAL T-MAX 100 dilution B small tank: 6½ / 6 / 5½ / 5 / 4 min | 400TX at 20 °C is **3¾ min, below Kodak's own 5-min uniformity warning** |
| J-86 T-MAX Developers (Dec 2017) | "intended to produce a contrast index of 0.60 for … T-MAX 400 Film and 0.56 for the other films"; times for both developers, push, rack-and-tank | Proprietary formula |

**J-109 Table 1 (small tank) for the D1 films.** Time in minutes. Full strength is given at 65/68/70/75/80 °F (18/20/21/24/27 °C); 1:1 at 68/70/75/80 °F. A dash means no time is given.

| Film | EI | CI aim | Full strength | 1:1 |
|---|---|---|---|---|
| 400TX 135 | 400/800 | 0.56 | 8 / 7 / 6¼ / 4¾ / — | 9 / 8½ / 7¼ / — |
| 400TX 135 | 1600 | 0.72 | 11¼ / 9¾ / 8¾ / 6¾ / — | 13¼ / 12¼ / 10½ / — |
| 400TX 135 | 3200 | 0.82 | — / 11½ / 10½ / 8 / — | 15½ / 14½ / 12¼ / — |
| 100TMX 135/120 | 25/50 | 0.52 | 8¼ / 6¾ / 6 / 4¾ / — | 9 / 8 / 6¼ / 4¾ |
| 100TMX | 100/200 | 0.56 | 8¾ / 7½ / 6½ / 5 / — | 9½ / 8½ / 6½ / 5 |
| 100TMX | 200 | 0.62 | 9¾ / 8 / 7¼ / 5½ / — | 10½ / 9½ / 7¼ / 5½ |
| 100TMX | 400 | 0.72 | 11½ / 9½ / 8½ / 6½ / — | 12¼ / 10¾ / 8¼ / 6¼ |
| 100TMX | 800 | 0.82 | 12¾ / 10½ / 9½ / 7¼ / — | 13¼ / 11¾ / 9 / 7 |

These give time-to-CI points at several temperatures and two dilutions for both XTOL pairs. That is richer than the F-4016 CI plot, and it is the best held-out material for **T-Max 100 × XTOL**.

There is no equivalent multi-CI table for HC-110. J-24 gives one CI-0.56 time per temperature.

---

## 5. Sensitometric exposure and ISO 6

### 5.1 What Kodak states about the exposure behind the curves

* **T-Max 100 (F-4016 2018, F-4016 2007, F-32 2002):** characteristic curves "Exposure: Daylight". **No exposure time** is given. MTF: tungsten. Spectral: effective exposure 1.4 s.
* **Tri-X 400 (F-4017 2007–2021):** characteristic curves state **no illuminant and no time**. The 2007 MTF caption is a placeholder ("Exposure: X"). Spectral: effective exposure 0.5 s.
* **Legacy TX (F-9, 2003):** "Exposure: Daylight, 1/50 second". This is the only explicit time found for Tri-X curves.
* **Consequence for the simulator's virtual sensitometer:** use ISO sensitometric daylight. Set exposure time to about 1/50–1/100 s, tagged `assumed`, with the F-9 value as the supporting reference. Both reciprocity tables show no exposure correction from 1/1000 s to 1/10 s, and no correction of any kind (exposure or development) from 1/100 s to 1/10 s (Tri-X does call for +10 % development at 1/1000 s). The curves are therefore insensitive to this choice inside 1/100–1/10 s, which contains the F-9 value of 1/50 s.

### 5.2 ISO 6:1993 (opened: text layer of the official preview at https://cdn.standards.iteh.ai/samples/3580/192a0f97667044f3ba51276b0ded056a/ISO-6-1993.pdf)

Title: "Photography — Black-and-white pictorial still camera negative film/process systems — Determination of ISO speed". Edition 2, 1993. ISO lists it as confirmed in 2023. **Opened** on https://www.iso.org/standard/3580.html on 2026-11-21: "This publication was last reviewed and confirmed in 2023"; the life cycle shows stage 90.93 "International Standard confirmed" dated 2023-12-20. The spot-check re-read the preview from the web.archive.org snapshot http://web.archive.org/web/20250617073625/https://cdn.standards.iteh.ai/samples/3580/192a0f97667044f3ba51276b0ded056a/ISO-6-1993.pdf. It confirmed 5.3, 5.6.3, the formulas in 6.1 and the Table 1 bins below.

* **Exposure (5.3):**
  * Exposure time "between 5 s and 1/1 000 s corresponding to the usage practice for the particular film tested". Because of reciprocity failure, the time "should be specified in use instructions".
  * Specimens are held at 23 ± 2 °C and (50 ± 5) % RH.
  * The sensitometer must be a non-intermittent, illuminance-scale type.
  * Illuminant per ISO 7589: "ISO sensitometric daylight, studio tungsten or photoflood".
  * Step increments must be ≤ 0.20 in log H; a continuous wedge must change by ≤ 0.04 log H per mm.
* **Processing (5.4):** the only processing specification is "that the density difference between points m and n shall be 0,80". Processing is done 5–10 days after exposure for general-purpose films, or 4 h to 7 days for professional films.
* **Densitometry (5.5):** ISO standard Type 2 diffuse transmission density, ISO 5-2 geometry and ISO 5-3 spectral conditions. This is Kodak's "diffuse visual".
* **Evaluation (5.6.3):** "Point m is located on the curve at a density of 0,10 above base plus fog density. Point n is located on the curve where the base 10 logarithm of the exposure is 1,30 units greater than that at point m." The development time is chosen so that ΔD(m→n) = 0.80. Base + fog comes from an unexposed sample processed with the strip. Annex A (informative) explains how to interpolate log H_m across a development series to ΔD = 0.80.
* **Formulas (6.1):** S = 0.80 / H_m (H_m in lux·s); S° = 1 + 10·log10(0.80 / H_m). Rounding uses Table 1:
  * ISO 400/27° ↔ log10 H_m from −2.75 to −2.66.
  * ISO 100/21° ↔ log10 H_m from −2.15 to −2.06.
* **Accuracy (6.3):** error in log10 H_m < 0.05. The product speed is taken from the mean log H_m over batches.
* **The contrast condition is not Kodak's CI.** ISO uses ΔD = 0.80 over Δlog H = 1.30, an average gradient of 0.615. Kodak's processing aim is **CI 0.56**, a different fixed-geometry slope measure. The measure was published by C. J. Niederpruem, C. N. Nelson and J. A. C. Yule, *Photographic Science and Engineering* **10**(1), January/February 1966, pp. 35–41. **Citation confirmed (secondary)** via US Patent 3,772,019, which cites it verbatim as "Phot. Sci. & Eng. vol 10, No. 1, January/February 1966, pages 35-41" and says the authors "term their measurement Contrast Index to distinguish it from the older Gamma" (read from https://patents.google.com/patent/US3772019A/en on 2026-11-21). **The paper itself was not opened.** The arc-construction geometry (gauge arc radii, the reference to gross fog) is **unverified**: only forum and training-manual descriptions were found, and they disagree on the construction. Confirm it from the 1966 paper or a Kodak sensitometry publication before implementing. The simulator must compute both quantities separately.

### 5.3 Speed check on the published curves (estimate, from vector paths)

Method:
* Base + fog is proxied by the lowest plotted density (the curves begin at log H ≈ −3.45 or −3.15, so b+f may be slightly overestimated).
* Find m at b+f + 0.10, then compute ΔD to log H_m + 1.30.
* Interpolate across the development series to ΔD = 0.80.
* Coarse solid-line paths (fewer than 20 vertices) are ignored where they gave non-monotonic results.

| Curve family | Bracketing curves (ΔD, log H_m) | log H_m at ΔD = 0.80 | ISO bin |
|---|---|---|---|
| Tri-X 400 35 mm, D-76 LT 20 °C | 8 min (0.71, −2.66), 10 min (0.85, −2.77) | **≈ −2.73** | 400 |
| Tri-X 400 120, D-76 LT 20 °C | 8 min (0.69, −2.60), 10 min (0.82, −2.70) | **≈ −2.69** | 400 |
| T-Max 100, D-76 ST 20 °C | 6 min (0.71, −2.01), 7.5 min (0.91, −2.13) | **≈ −2.06** | 100 (at the 100/80 boundary) |

The absolute log H axes of the published curves are therefore consistent with the box speeds. A virtual-sensitometer speed test against ISO 6 is a meaningful validation of the D-76 calibration.

---

## 6. Practical notes for digitising

* In the 2016–2021 F-4017 and the 2018 F-4016, every curve (characteristic, CI, MTF, spectral) is a **vector path** (PyMuPDF `get_drawings`): characteristic curves have 136 segments, CI curves 3–24. Chart frames are `re` rectangles. Legend samples carry the same dash pattern and line width as their curves.
* Exact digitisation is possible by mapping each frame to its axis ranges. Characteristic curves: log H −4 to +1 (bar notation on the axis), density 0–4, except the T-Max 100 T-MAX Developer chart, which is 0–3. CI: Tri-X 0–30 min and CI 0–1.2; T-Max 4–18 min and CI 0.3–0.9.
* This removes hand-digitising error. Record the PDF edition and page in each reference CSV's provenance.
* Check the axis range of each frame individually. Mixed scales occur on the same page.

---

## 7. Feasibility of the planned D1 split

| Pair | Role (plan) | Characteristic curves | CI(t) curve | Time tables (CI 0.56) | Extra |
|---|---|---|---|---|---|
| Tri-X × D-76 | train | **Yes**: LT 1-min, 20 °C, 6/8/10/12 min, 35 mm and 120 | Yes (ST 20 °C, 5.4–10.8 min) | ST, LT, rotary, 18–24 °C; push | MTF (D-76 LT) |
| Tri-X × XTOL | train | **No** | Yes (ST, 3.25–14 min, CI 0.34–0.98), plus XTOL 1:1 | ST/LT/rotary, 18–24 °C; push | J-109 Table 1: CI 0.56/0.72/0.82, 18–27 °C, 1:1 |
| Tri-X × HC-110 | train | **No** (the HC-110 B curves in F-4017 are 320TXP sheets, a different emulsion) | Yes (ST, 3.25–11.9 min, CI 0.53–1.06) | ST/LT/rotary; push | **Granularity 17** is measured in HC-110 B |
| T-Max 100 × D-76 | train | **Yes**: **ST** (30-s per F-4016 procedure, not caption), 20 °C, 6/7.5/10 min | ST and LT, CI 0.56–0.82 | ST/LT/tray/rotary; push; temperature factors | Granularity 8, MTF, spectral all in D-76 |
| T-Max 100 × XTOL | held out | **No** | ST, LT, rotary, tray, plus XTOL 1:1, CI 0.56–0.82 | Yes; push | J-109 Table 1: CI 0.52/0.56/0.62/0.72/0.82 × 4–5 temperatures × 2 dilutions |
| T-Max 100 × HC-110 | held out | **No** | ST, LT, CI 0.56–0.82 | Yes; push; temperature factors (F-4016 p. 7) | — |
| Tri-X × T-MAX Dev | (not in plan) | **Yes**: ST 30-s, 20 °C, 6/7/9/11 min, 35 mm and 120 | Yes | Yes | — |
| T-Max 100 × T-MAX Dev | (not in plan) | **Yes**: ST, 20 °C, 6/7/10/12 min | Yes | Yes | T-MAX RS: LT curves 8/10.5/13/15 min |
| Tri-X / T-Max × Rodinal | predicted | none (no Kodak data) | none | none | — |

**Verdict.** The planned split is feasible only for a **CI(t) and time-to-CI** validation of the held-out pairs:
* CI at given times.
* Time to reach CI 0.56 at 18–24 °C.
* For XTOL, time to reach CI 0.52–0.82, from J-109.

It is **not** feasible for curve shape, toe, shoulder, D-max or speed, because neither held-out pair, nor Tri-X × XTOL or Tri-X × HC-110 in training, has a characteristic curve.

**Nearest feasible splits:**

1. **Keep the plan and change the observable (minimal change).**
   * Film parameters: D-76 characteristic curves for both films.
   * Chemistry constants: trained on Tri-X × {D-76, XTOL, HC-110} CI(t) curves, time–temperature tables and push tables (J-109 Table 1 for XTOL), plus T-Max 100 × D-76 CI(t).
   * Held out: T-Max 100 × XTOL (CI(t) in four agitation modes plus J-109 multi-CI times) and T-Max 100 × HC-110 (CI(t) small and large tank, time table, temperature factors).
   * State in the validation report that the held-out test is a contrast-kinetics test, not a curve-shape test.
2. **Add one full-curve held-out pair (recommended in addition to option 1).**
   * Train on Tri-X × T-MAX Developer characteristic curves.
   * Hold out **T-Max 100 × T-MAX Developer characteristic curves** (6/7/10/12 min, small tank, 20 °C).
   * This is the only pair that tests toe, shape and speed transfer across films with a developer whose constants came from the other film.
   * Cost: T-MAX Developer must be added to D1 with an `msds-estimated` formula (as HC-110 already is). Its formula uncertainty confounds the test.
   * T-MAX RS (T-Max 100 LT curves) could serve the same purpose but has no Tri-X characteristic curve.
3. **Secondary check.** Once the D-76 film parameters are calibrated on 35 mm, the Tri-X 35 mm and 120 D-76 curves (the same process, two formats) give a format/base consistency check: 120 base + fog is about 0.1 lower.

**Agitation caveat.** "Film parameters from D-76" mixes regimes:
* Tri-X D-76 curves: large tank, 1-min intervals.
* T-Max 100 D-76 curves: small tank, 30-s intervals.

The boundary-layer and agitation model therefore enters the film calibration itself. The Tri-X small-tank D-76 CI(t) curve can be used to link the two regimes.
