# D1 developer formulas and species data: D-76, XTOL, HC-110, Rodinal

Research date: 2026-11-21. Scope:
* Developer formulas in g/L of working stock, with dilutions and working-solution pH.
* Species data for every chemical that appears: molar mass, CAS, pKa, aqueous diffusion coefficient, redox potential.
* Buffer chemistry of each developer.

The aim is to feed `data/developers/*.yaml` and `data/chemistry/species.yaml` for the D1 black-and-white chain.

**Confidence labels used throughout:**
* **opened**: I read the source myself (text extracted or page rendered).
* **secondary**: taken from a search-engine summary or snippet of a document I could not open under the fetch rules.
* **recalled**: from memory of the literature, not opened in this job. Each recalled item says what would confirm it.
* **estimate**: my own computation from opened data. The method is stated.
* **unverified (spot-check)**: re-checked by an independent reviewer (§8) and not found on a permitted host; do not tag as `literature` or `datasheet` provenance on this basis.

Developer formula classes follow the design spec:
* `published`: a public formula.
* `patent`: a formula from a patent example.
* `msds-estimated`: inferred from safety-data-sheet ranges.

**Sources opened and cached** (under `.cache/datasheets/`, git-ignored):

| Cached file | What it is | Host / provenance |
|---|---|---|
| `patents/US5756271A.html/.txt` | US 5,756,271 (Opitz & Zawadzki, Eastman Kodak), "Weakly alkaline ascorbic acid developing composition, processing kit and method using same". Filed 1997-05-27, priority 1996-06-28, granted 1998-05-26. | patents.google.com |
| `patents/US5702875A.html` | US 5,702,875 (same title; parent application 08/674,497) | patents.google.com |
| `patents/DE1118005B_de.txt`, `GB958678A.txt` | DE 1,118,005 (Eastman Kodak; inventors J.C. Barnes, C.W. Zuehlke, G.W. Luckey; priority 1959-09-04, published 1961-11-23) and its family member GB 958,678 (Kodak Ltd, published 1964-05-21): "Flüssiges Konzentrat eines photographischen Entwicklers" / "Photographic developer concentrates" | patents.google.com (German original text used) |
| `books/cinematographica00unse.txt`, `cinematographi02unse.txt` | *Cinematographic Annual* vol. 1 (1930) and vol. 2 (1931), American Society of Cinematographers, scanned full text | archive.org (Media History Digital Library scan) |
| `sds/hc110_psi_v5_2024.pdf/.txt` | SDS "KODAK PROFESSIONAL HC Developer", product no. 1058692, Photo Systems Inc. (Dexter, MI), version 5.0, 7/22/2024 | cdn.shopify.com/s/files/1/0339/5113/… — a Shopify store CDN. The same store hosts Kodak film data sheets (see sibling `datasheets.md`). **I could not tie the store to Kodak Alaris or PSI with certainty**; the document itself names PSI as supplier. |
| `sds/xtolA_psi_v2_2024.pdf/.txt` | SDS "KODAK PROFESSIONAL XTOL Developer Powder A", product 1058338A, PSI, version 2.0, 2/20/2024 | same Shopify CDN |
| `sds/d76_psi_v2_2024.pdf/.txt` | SDS "KODAK PROFESSIONAL D-76 Developer Powder", product 5160304, PSI, version 2.0, 7/24/2024 | same Shopify CDN |
| `sds/adox_rodinal_msds_en.pdf/.txt` | SDS "Rodinal", ADOX (www.adox.de), version 1.0, 6/24/2024 | fotoimpex.com, linked from ADOX's own Rodinal/Adonal product page (ADOX and Fotoimpex are the same Berlin group) |
| `sds/adox_rodinal_datasheet.pdf` (+ PNG renders) | AGFA-branded RODINAL instruction leaflet, a&o imaging solutions GmbH (Koblenz), code "13 144 140 - 1106", revision mark "F1, 2, 4, 5, 13" | fotoimpex.com, linked from adox.de |
| `pubchem/iupac_hc.csv`, `iupac_refs.csv` | IUPAC Digitized pKa Dataset v2.4 (high-confidence set) and its reference-code table. This is the digitization of Perrin (1965, 1972 supplement) and Serjeant & Dempsey (1979). | github.com/IUPAC/Dissociation-Constants |
| `pubchem/cid*.json` | PubChem compound records (molar mass, CAS, HSDB dissociation constants) | pubchem.ncbi.nlm.nih.gov |
| `ds_txt/j78_2017.txt`, `j109_2018.txt`, `j109_1998.txt`, `j24_2017.txt` (pre-existing) | Kodak J-78 (Dec 2017), J-109 (Feb 2018 and Apr 1998), J-24 (Dec 2017) | cited as the Kodak publications |
| `spotcheck_devf/tsmpe37.txt` + page images of pp. 413–414, 419 | *Transactions of the S.M.P.E.* vol. XIII, no. 38 (1929), incl. H.C. Carlton & J.I. Crabtree, "Some properties of fine-grain developers for motion picture film", pp. 406–444 (Communication No. 388, Kodak Research Laboratories) | archive.org scan (Library of Congress copy); added by the spot-check |
| `spotcheck_devf/gp1985.pdf` | R.N. Goldberg & V.B. Parker, "Thermodynamics of solution of SO₂(g) in water and of aqueous sulfur dioxide solutions", *J. Res. NBS* 90(5), 341–358 (1985) | nvlpubs.nist.gov; added by the spot-check |
| `spotcheck_devf/mda1944.pdf` | G.G. Manov, N.J. DeLollis & S.F. Acree, "Ionization constant of boric acid and the pH of certain borax-chloride buffer solutions from 0° to 60° C", *J. Res. NBS* 33, 287 (1944), RP1609 | nvlpubs.nist.gov; added by the spot-check |
| `spotcheck_devf/iupac.csv`, `refs.csv` | fresh copy of the IUPAC pKa dataset v2.4 and its reference-code table | github.com/IUPAC/Dissociation-Constants; added by the spot-check |

**Not fetched because of the host rule:**
* The 125px.com mirror of Kodak J-1 *Processing Chemicals and Formulas* (1977) and of `Kodak_formulary.pdf`.
* SDS copies hosted by universities and dealers (Freestyle, Beau Photo, Agar Scientific, Parallax, digitaltruth).

Figures from those dealer-hosted copies appear only as **secondary**, via search snippets. The Kodak Alaris official SDS portal (3E Online, "KodakAlaris_14831_SDSHosting", reached from corporate.kodakalaris.com → Safety Data) was searched. It lists 180 products (inks, toners, microfilm chemistry) and **no Kodak Professional B&W developer**: the B&W chemistry moved to Photo Systems Inc. web.archive.org was "temporarily offline" for the whole session.

---

## 0. Answer in brief

| Developer | Best formula available | Class | Working pH (measured, with source) | Dilutions (factor) |
|---|---|---|---|---|
| **D-76 (published)** | Metol 2 g, Na₂SO₃ (anhydrous) 100 g, hydroquinone 5 g, borax (decahydrate) 2 g, water to 1 L | published | **8.5** (Kodak Research Laboratories, Carlton & Crabtree, *Trans. SMPE* 13 (1929), Table II, expt. 1; LaMotte colorimetric indicator, ±0.1 pH, measurement temperature not stated; opened by spot-check). Hobby literature: about 8.3 fresh, rising to about 9 on storage (secondary). Ideal-solution estimate 8.94 is an upper bound. | stock (×1), 1:1 (×2) — J-78 2017 |
| **D-76 (packaged, what the Kodak curves use)** | Same actives, plus **diboron trioxide 1–3 %** and **DTPA < 1 %**; borax supplied as pentahydrate. Exact amounts withheld. | **msds-estimated** | not published | as above |
| **XTOL** | US 5,756,271 Ex. 1. Part A: Na₂SO₃ 10 g, DTPA-Na₅ 1 g, sodium metaborate ("8 mol") 4.0 g, Dimezone S 0.2 g. Part B: Na₂SO₃ 75 g, Na₂S₂O₅ 3.5 g, sodium isoascorbate 12 g. Water to 1 L. | patent (the commercial product differs in detail, see §2) | **8.20 ± 0.05** (patent Ex. 1, opened). J-109: "pH 8.2 ± 0.05" with SG 1.085 ± 0.003 at 25 °C (opened). Part A alone in 850 mL: pH 10.19 ± 0.05 (patent). | stock, 1:1 (×2) — J-109 2018. 1:1, 1:2, 1:3 — J-109 1998. |
| **HC-110 (old, about 1960s–2018)** | Liquid amine–SO₂ concentrate (DE 1,118,005 / GB 958,678 type). SDS: DEA 30–35 %, SO₂ 15–20 %, HQ 5–10 %, DEG 5–10 %, ethanolamine 5–10 %, KBr, 4-methyl-1-phenyl-3-pyrazolidinone. | msds-estimated (secondary only; **unverified** by spot-check, no permitted-host copy found) | not found. Patent analogues: pH 9.6–10.0 at 1+20. | J-24: A 1+15 (×16), **B 1+31 (×32)**, C 1+19, D 1+39, E 1+47, F 1+79 from concentrate |
| **HC-110 (new, ≥2019; "KODAK PROFESSIONAL HC Developer")** | PSI SDS 2024: K₂SO₃ solution 45 % at 40–60 %, HQ 10–15 %, DEG 5–10 %, borax pentahydrate 3–5 %, KOH 45 % at 1–3 %, NaBr 1–3 %, HEDTA-Na₃ (40 %) 1–3 %, DEA 1–3 %, Dimezone S < 1 % | msds-estimated | concentrate pH 9.0; "pH in solution 9.3" at an unstated dilution (SDS §9, opened) | same as J-24 |
| **Rodinal / Adonal (predicted case)** | p-aminophenol 3–5 %, KOH 1–3 %, KBr 1–3 % (ADOX SDS 2024). Sulfite present but not listed. | msds-estimated | concentrate pH > 13 (ADOX SDS). Working pH not found. | 1+25 (×26), 1+50 (×51) (Agfa / a&o leaflet). ADOX: 1+25 to 1+500. |

**Three consequences for D1:**
1. **Bracket the formula-version risk.** The Kodak reference data for the HC-110 and D-76 pairs come from **old HC-110** and **packaged D-76**. The film data in F-4016/F-4017 and J-24 (Dec 2017) predate the HC-110 change. Neither of these is the formula with the best documentation.
2. **Calibrate pH to the measured anchors.** Ideal-solution speciation overestimates the pH of the sulfite/borate developers at their ionic strength (about 2 mol/L). The offset is formula-dependent: about −0.7 to −0.8 on the XTOL anchors (patent, glass electrode) but about −0.4 on the D-76 anchor (Carlton & Crabtree 1929, colorimetric). The model needs an activity correction or a fitted effective pKa, calibrated on both anchors, not a single constant offset.
3. **Expect sulfite to drive the HC-110 contrast.** At dilution B, HC-110 carries only **about 40–80 mM sulfite**, against 0.71–0.79 M in D-76 and XTOL. It also contains 1–3 % bromide in the concentrate, about 3–9 mM at working strength. These two differences are what the HC-110 hold-out actually tests.

---

## 1. Kodak D-76

### 1.1 Published formula

| Ingredient | g/L (working stock) | mol/L | Source |
|---|---|---|---|
| Metol (Elon; 4-(methylamino)phenol sulfate) | 2.0 | 5.81 mM salt = **11.6 mM metol** | see below |
| Sodium sulfite, anhydrous | 100 | **0.793** | |
| Hydroquinone | 5.0 | **45.4 mM** | |
| Borax (Na₂B₄O₇·10H₂O) | 2.0 | 5.24 mM = **21.0 mM total B** | |
| Water to | 1 L | | |

* **Primary source opened:** *Cinematographic Annual* vol. 1 (ASC, Hollywood, 1930), article on sensitometry in motion-picture practice, pp. 119–120, "Formula D-76". It gives Elon 120 grains, sodium sulphite (E.K. Co.) 14 oz, hydroquinone 300 grains, borax 120 grains, water 128 oz (1 gallon), at 65 °F. The text calls it the "borax formula", now used by "practically all of the laboratories".
  * Conversion (estimate): 120 gr/gal = 2.05 g/L; 14 oz/gal = 104.9 g/L; 300 gr/gal = 5.14 g/L. This is within rounding of the metric 2/100/5/2 g/L.
  * The same volume, in a later fine-grain section (scan line ≈ 28 996), gives "Formula D-76 A" (as above) and "B": Elon 160 gr, sulphite 14 oz, hydroquinone 160 gr, borax 120 gr.
* **Disagreement:** vol. 2 (1931) prints the "A" sulphite as **11 oz** where vol. 1 prints 14 oz. This may be an OCR or typesetting error; the "B" column reads 14 oz in both. Carlton & Crabtree (1929, below) use 100 g/L sulfite, which supports 14 oz/gal (104.9 g/L); 11 oz (82 g/L) is the outlier.
* **Metric 2/100/5/2 g/L: opened (spot-check).** Carlton, H.C. & Crabtree, J.I., "Some properties of fine-grain developers for motion picture film", *Trans. SMPE* XIII no. 38, pp. 406–444 (1929), Communication No. 388 from the Kodak Research Laboratories (archive.org scan of the Library of Congress copy). Table II (p. 413), expt. 1, "regular borax developer": Elon 2, hydroquinone 5, sodium sulfite 100, borax 2 g per litre. The text (p. 420) repeats "2 grams of Elon, 5 grams of hydroquinone, and 100 grams of sodium sulfite per liter". The earlier attribution of the metric form to Kodak J-1 was recalled and is not needed.
* **Origin:** the paper's first sentence refers to an *earlier* publication of "the Kodak fine-grain developer formula for motion picture film", its ref. 1: *Eastman Duplicating Film*, Eastman Kodak Co., Rochester, 2nd ed. (1927). That booklet was not opened. Carlton & Crabtree 1929 is a follow-up study, not the first publication.
* **Class: published.** Kodak's current J-78 (Dec 2017) gives **no formula and no pH** (opened; the text was searched).

### 1.2 Packaged D-76 (what the Kodak reference curves were made with)

PSI SDS v2.0, 7/24/2024, §3 (opened). Ranges are w/w of the powder; exact values are "withheld as a trade secret".

| Component | CAS | % w/w | g/L if the kit is ≈ 105–125 g/L (estimate) |
|---|---|---|---|
| Sodium sulfite | 7757-83-7 | 80–95 | 84–119 |
| Hydroquinone | 123-31-9 | 3–5 | 3.2–6.3 |
| Borax pentahydrate | 12179-04-3 | 1–3 | 1.1–3.8 (2 g decahydrate ≙ 1.53 g pentahydrate) |
| bis(4-hydroxy-N-methylanilinium) sulphate (metol) | 55-55-0 | 1–3 | 1.1–3.8 |
| **Diboron trioxide (B₂O₃ → boric acid in water)** | 1303-86-2 | **1–3** | **1.1–3.8 (≙ 1.9–6.7 g H₃BO₃)** |
| **Diethylenetriaminepentaacetic acid** | 67-43-6 | **< 1** | < 1.2 |

How packaged D-76 differs from the published formula:
1. It contains **boric acid (as B₂O₃)** in addition to borax. This is a borax/boric-acid buffer of the D-76d kind. It lowers and stabilises the pH, increases buffer capacity, and counters the rise in pH on storage.
2. It contains a **sequestrant (DTPA)**.
3. Borax is supplied as the **pentahydrate**.

The SDS total kit mass per litre is not stated. The g/L column assumes sulfite fixed at 100 g/L. **Class: msds-estimated.**

* Older sheets (secondary, search snippets of dealer-hosted copies):
  * 2001 sheet: sulfite 85–90 %, borax pentahydrate 1–5 %, HQ 4 %, metol 1–5 %.
  * 2015 Kodak Alaris SDS v5.0: sulfite 85–90 %, HQ 1–5 %, metol 1–5 %.
* **Do not merge** the 2002 MSDS "200000263", which shows borax 10–15 % and HQ 6 %. Its filename ("D-76_Dev_Rep") marks it as the **D-76 replenisher** (D-76R) sheet, not the developer.

### 1.3 Dilutions and working pH

* Dilutions: full strength or **1:1** (developer + water, factor 2), mixed just before use and discarded after one batch (J-78 2017, pp. 2 and 5, opened).
* Packaged D-76 with D-76R replenisher is also used; for T-MAX films the modified replenisher is 5 parts D-76 + 1 part D-76R (J-78 p. 6).
* **Working pH: 8.5** for the published formula, Carlton & Crabtree (1929), Table II, expt. 1 (opened by spot-check).
  * Method: LaMotte colorimetric indicator, stated error ±0.1 pH; the authors say "the absolute accuracy of the pH values is open to question", only relative values are reliable. Measurement temperature not stated (the development tests ran at 65 °F).
  * The same table gives a borax series on the 2/5/100 base: **no borax 8.1** (expt. 12), **2 g/L 8.5** (expt. 1), **6 g/L 8.9** (expt. 11), **20 g/L 9.2** (expt. 10); also Elon 4 + HQ 5 + borax 2: 8.1 (expt. 8); Elon 2 + borax 2, no HQ: 8.6 (expt. 14). Use this series as a test of whatever activity model D1 adopts, not to derive a new offset constant.
  * The paper also states (p. 419) that borax and Elon sulfate neutralise each other weight for weight (381 : 344), so the 2 g/2 g solution is "essentially a solution of Elon base ... and hydroquinone in a 10% solution of sodium sulfite".
  * Hobby literature (Anchell & Troop, *The Film Developing Cookbook*, as quoted on photrio.com) gives about 8.3 fresh, rising to about 9.0 on storage. **secondary / forum**; the storage rise has no primary here.
  * The mechanism of the rise is hydroquinone sulfonation, Q + SO₃²⁻ + H₂O → HQ-SO₃⁻ + OH⁻, which releases hydroxide. This is standard in James (ed.), *The Theory of the Photographic Process*, 4th ed. (1977). **recalled**.
* Estimate (§5): ideal-solution speciation of the published formula gives pH **8.94** (reproduced by the spot-check, 8.94). Against the measured 8.5 the offset is about **−0.4**, smaller than the −0.7 to −0.8 seen on XTOL. The earlier statement "≈ 8.2 after the XTOL offset, consistent with 8.3" does not hold against the Kodak measurement.

---

## 2. Kodak XTOL

### 2.1 Patent formula (the standard reconstruction)

US 5,756,271, Example 1, "Preferred Developing Kit" (opened, patent text). Part A is dissolved in 850 mL water, then Part B is added and the volume made up to 1 L.

| Part | Ingredient | g/L | mol/L (estimate) | Notes |
|---|---|---|---|---|
| A | Sodium sulfite | 10 | 79.3 mM | |
| A | DTPA, pentasodium salt | 1 | 1.99 mM (M = 503.3) | |
| A | Sodium metaborate "(8 mol)" | 4.0 | 29.0 mM B | "8 mol" read as the trade grade Na₂B₂O₄·8H₂O ≡ NaBO₂·4H₂O, M = 137.86 per B (**recalled** trade naming; confirm on a borate producer's product sheet) |
| A | 4-hydroxymethyl-4-methyl-1-phenyl-3-pyrazolidone (Dimezone S) | 0.2 | **0.97 mM** | |
| B | Sodium sulfite | 75 | 0.595 M | total sulfite **0.674 M** |
| B | Sodium metabisulfite | 3.5 | 18.4 mM (= 36.8 mM HSO₃⁻) | |
| B | Sodium isoascorbate (erythorbate) | 12 | **55.5 mM** if monohydrate (216.12); 60.6 mM if anhydrous (198.11) | hydration not stated in the patent |

* Part A solution: **pH 10.19 ± 0.05**. Final: **pH 8.20 ± 0.05** (patent, Example 1).
* Ascorbate : Dimezone S ≈ **57–62 : 1** (estimate). The patent's preferred range is 20:1 to 50:1, so Example 1 lies slightly above it.
* Patent ranges for the aqueous composition (claims and description):

| Component | Range (mol/L) | Preferred (mol/L) |
|---|---|---|
| Ascorbic-acid-type developer | 0.01–0.1 | 0.02–0.07 |
| Auxiliary super-additive developer | 2.5×10⁻⁴ to 2.5×10⁻² | 5×10⁻⁴ to 2×10⁻³ |
| Borate, "used as the only buffer" | up to 0.2 | 0.001–0.16 |
| Sulfite preservative | 0.4–1.6 | 0.4–0.8 |
| Sequestrant | 0–0.02 | — |

  * pH range: "at least 7 and up to, but not including 9"; preferred 8–8.5; more preferably 8.0–8.4.
* Two-part logic (patent):
  * One formulation holds the ascorbate plus sulfite, with preservative : developer ≥ 4:1 (preferred 7:1 to 20:1).
  * The other holds the pyrazolidone, the borate and sulfite (≥ 4:1, preferred 50:1 to 200:1).
  * The pyrazolidone/borate part is dissolved first and is "more alkaline (generally from 9.5 to 11)".
  * This matches J-109: Part A is dissolved first and the solution "may appear somewhat tawny or copper-colored"; the tint clears when Part B is added (J-109 Feb 2018, p. 1, opened).
* Patent Examples 2–5 (liquid variants, all "pH = 8.0", g per 1 L):
  * Ex. 2: NaOH 50 % 12.35, **diethanolamine with 16 % SO₂** 75, Dimezone S 0.55, DTPA-Na₅ 40 % 2.5, Na₂S₂O₅ 18.8, ascorbic acid 8.0.
  * Ex. 4: DEA/SO₂ 15.3, Dimezone S 0.67, DTPA-Na₅ 40 % 1.67, ascorbic acid 10.0, Na₂SO₃ 36.0, **NaHCO₃ 13.3**.
  * The patent prints "sodium metabisulfate"; read it as metabisulfite.
* The patent also has **Table I**, comparing Example 1 with D-76 and T-MAX Developer on TMAX 100, TMAX 400, TMAX P3200, TRI-X PAN and PLUS-X PAN. It gives time to CI 0.56, time to CI 0.82, EI, and Dmin. Conditions: 20 °C, 8 L tank, nitrogen-burst agitation 1 s every 15 s, 21-step tablet. For D1 this is a **paired XTOL-type vs D-76 data set** on the same films (§6).

### 2.2 Commercial XTOL: what the SDS say (msds-estimated)

| Part | Source | Composition (% w/w) | Status |
|---|---|---|---|
| A | PSI SDS v2.0, 2/20/2024 | Na₂SO₃ **80–95**; borax pentahydrate **3–5**; DTPA (acid, 67-43-6) **1–3**; Dimezone S (13047-13-7) **< 0.01** | opened |
| A | Kodak Alaris SDS 2015 v3.0 | Na₂SO₃ 90–95; sodium tetraborate pentahydrate 1–5; DTPA pentasodium 1–5; Dimezone S **0.1 to < 1** | secondary |
| B | PSI SDS (Feb 2024) | Na₂SO₃ **60–80**; sodium erythorbate (6381-77-7) **15–25** | secondary (snippet) |
| B | Kodak Alaris SDS 2019 | Na₂SO₃ 75–80; Na₂S₂O₅ 1–5 (erythorbate not shown in snippet) | secondary |
| B | Kodak MSDS 2014 (02/07/2014 v3.0) | lists only Na₂SO₃ and Na₂S₂O₅ as hazardous components | secondary |

Differences from patent Example 1:
1. The commercial Part A buffers with **borax pentahydrate**, not metaborate. Patent Ex. 1 has metaborate at 26 % of Part A.
2. **Disagreement on Dimezone S.** In Ex. 1, Dimezone S is about 1.3 % of Part A (0.2 g of 15.2 g). Kodak 2015 gives 0.1 to < 1 %. PSI 2024 gives **< 0.01 %**. The PSI figure would mean less than about 2 mg/L in working solution, which is implausible for the developer's activity. Treat it as a possible SDS error, or as a sign that the Dimezone S has moved to Part B; the Part B snippet does not show it.
3. The Part B ranges are consistent with Ex. 1: sulfite 82.9 %, metabisulfite 3.9 %, isoascorbate 13.3 % of Part B (estimate).

Kit masses per litre are **not found**, so commercial g/L cannot be computed. **Class: patent (Ex. 1) as the working formula, cross-checked against msds-estimated.**

### 2.3 Dilutions, pH, specific gravity

* J-109 Feb 2018 (opened):
  * Full strength (replenishable).
  * **1:1** one-shot, "slightly greater film speed, enhanced sharpness and shadow detail, and slightly more grain"; minimum 100 mL of stock per 80 in² of film.
  * Working tank solution: SG **1.085 ± 0.003 at 25 °C, pH 8.2 ± 0.05**.
* J-109 Apr 1998 (opened): "may dilute ... 1:1, 1:2, or 1:3" (factors 2, 3, 4), with "pH 8.2 ± 0.05".
* Working pH at 1:1–1:3: **not found** in any Kodak document. The ideal model predicts no change (the buffer ratio is preserved); see §5 for the activity effect.

---

## 3. Kodak HC-110

**Uncertainty warning.** HC-110 is a proprietary liquid concentrate. The formula **changed around 2019**: the trade name on the 2024 SDS is "KODAK PROFESSIONAL HC Developer", same product number 1058692. No public formula exists for either version. All numbers below are **msds-estimated**, and an SDS lists only classified ingredients, so the ranges are **not a closed formula**.

### 3.1 Old formulation (behind the Kodak reference data for D1)

Composition from SDS versions dated 09/11/2009, 05/02/2014 (v4.0, Z17000000795) and 05/26/2016 (v4.1). All **secondary** (search snippets of dealer- and university-hosted copies). **Unverified (spot-check):** no copy on a Kodak, Kodak Alaris or PSI host, and none on web.archive.org, could be found, so none of the numbers below was seen in a source; treat the table as a prior only.

| Component | CAS | % w/w |
|---|---|---|
| Diethanolamine | 111-42-2 | 30–35 (an EU sheet: 35.0) |
| Sulphur dioxide | 7446-09-5 | 15–20 |
| Hydroquinone | 123-31-9 | 5–10 (an EU sheet: 10.0) |
| Diethylene glycol | 111-46-6 | 5–10 |
| Ethanolamine | 141-43-5 | 5–10 |
| Potassium bromide | 7758-02-3 | listed as "contains", % not seen |
| 4-methyl-1-phenyl-3-pyrazolidinone | 2654-57-1 | listed in the 2014/2016 sheets, % not seen |
| Ethylene glycol | 107-21-1 | 0.1 to < 1 (2009 sheet) |
| Chelating agent | — | mentioned, identity not seen |

**Patent basis (opened):** DE 1,118,005 / GB 958,678 (Kodak, 1959 priority). The claims cover an anhydrous concentrate of:
* a hydroxyalkylamine,
* its SO₂ addition product (or sulfite),
* a developing agent soluble in that mixture, such as hydroquinone with 1-phenyl-, 4-methyl-1-phenyl- or 4,4-dimethyl-1-phenyl-3-pyrazolidone,

diluted ≥ 20-fold for use.

Representative examples, g per batch (German text):

| Ex. | Developing agents | Amine | Amine–SO₂ adduct | Other | Dilution | pH |
|---|---|---|---|---|---|---|
| 1 | 4-methyl-1-phenyl-3-pyrazolidone 1.0, HQ 10.0 | 2-aminoethanol 20.0 | 30.0 | benzotriazole 0.8, methylbenzotriazole 0.2 | 1+19 | not stated |
| 6 | 4,4-dimethyl-1-phenyl-3-pyrazolidone 1.0, HQ 10.0 | 2-ethylaminoethanol 16.0 | 34.0 | 2-methylaminoethanol·HBr 6.0 | 1+20 | **9.9** |
| 7 | 4,4-dimethyl analogue 1.0, HQ 10.0 | diethylaminoethanol 22.0 | 28.0 | — | 1+20 | **9.6** |
| 8 | 4-methyl analogue 1.0, HQ 10.0 | 2-aminoethanol 21.5 | 18.5 | — | 1+20 | **10.0** |
| 29 (diethanolamine) | 4,4-dimethyl analogue 0.25, HQ 6.0 | 2,2'-iminodiethanol 26.0 | 137.0 | Na₂S₂O₃ 6.0, water 30 mL | 1+20 | not stated |

The text of Ex. 1/2 says the amine and its SO₂ adduct are present in "equimolar" proportions.

* The DEA/SO₂/HQ/4-methyl-phenidone/KBr profile in the old SDS matches this patent family. That HC-110 is **this** patent is an inference; it is **not stated** anywhere I could open.
* Approximate old-HC-110 working concentrations at **dilution B (1+31, factor 32)**. Estimate from the SDS ranges with an **assumed** concentrate density 1.1–1.2 kg/L; the density was not found:

| Species | Working concentration |
|---|---|
| HQ | 1.7–3.8 g/L (15–34 mM) |
| Total S(IV) | 50–75 mM |
| Amines (DEA + MEA) | 120–170 mM |
| Bromide | unknown % (likely about 1–3 %, i.e. about 3–9 mM; **estimate** by analogy with the new SDS) |

* **pH of the old formulation is not determinable from the SDS.**
  * Ideal speciation of the SDS mid-ranges gives pH 6.4–7.3. That is impossible for an active HQ/pyrazolidone developer.
  * So either an undisclosed alkali is present, or the "SO₂ %" accounting (adduct stoichiometry) differs from my reading.
  * Use the patent analogues (**pH 9.6–10.0 at 1+20**, opened) as the prior. Treat the working pH as `fitted`.

### 3.2 New formulation (2019 onward)

PSI SDS v5.0, 7/22/2024, §3 and §9 (opened):

| Component | CAS | % w/w |
|---|---|---|
| Potassium sulfite solution 45 % | 10117-38-1 | 40–60 (≙ 18–27 % K₂SO₃) |
| Hydroquinone | 123-31-9 | 10–15 |
| 2,2'-oxydiethanol (diethylene glycol) | 111-46-6 | 5–10 |
| Borax pentahydrate | 12179-04-3 | 3–5 |
| Potassium hydroxide 45 % | 1310-58-3 | 1–3 |
| Sodium bromide | 7647-15-6 | 1–3 |
| Dissolvine H-40 (HEDTA trisodium, about 40 % solution) | 139-89-9 | 1–3 |
| Diethanolamine | 111-42-2 | 1–3 |
| Dimezone S | 13047-13-7 | < 1 |
| Pyrocatechol (catechol) | 120-80-9 | < 0.05 (omitted in the first version of this report; added by spot-check from the same SDS §3) |

§9 physical data:
* pH **9.0**; "pH in solution 9.3 (%)", dilution not stated.
* Relative density **1.07**.
* Vapour pressure 18 mbar; boiling point 100 °C; flash point 93.3 °C.

**Flag:** relative density 1.07 looks low for 18–27 % K₂SO₃ plus 10–15 % HQ. My guess is about 1.2–1.3 (**estimate**; not verified).

Other versions (secondary):
* Kodak Alaris SDS 2019 v02: HQ 10 to < 20, K₂SO₃ 10 to < 20, DEG 5 to < 10, borax pentahydrate 1 to < 3, KOH 0.1 to < 1, DEA 0.1–1.
* PSI v3.0, 2/22/2024: same as v5.0, but the snippet omits HQ.

Working concentrations at **dilution B (factor 32)**, estimate, density 1.07–1.25:

| Species | Working concentration |
|---|---|
| HQ | 3.3–5.9 g/L (30–53 mM), comparable to D-76's 5 g/L |
| K₂SO₃ | 6–10.5 g/L (**38–67 mM**) |
| Borate | 2.6–3.9 mM borax ≙ 10–16 mM B |
| NaBr | 0.33–1.17 g/L (**3.3–11 mM Br⁻**) |
| Dimezone S | < 0.33–0.39 g/L (< 1.6–1.9 mM) |

Ideal-solution pH at 1+31 = 9.26 (an upper bound; §5). This is consistent with the SDS's "9.3 in solution".

### 3.3 Dilutions (J-24 Dec 2017, p. 2, opened)

* Stock solution = 1 part concentrate + 3 parts water.
* **The ratios are "part + parts"**: "1:31" means 1+31, a dilution factor of 32.

| Dilution | From stock (stock : water) | From concentrate (concentrate : water) | Factor vs concentrate |
|---|---|---|---|
| A | 1:3 | 1:15 | 16 |
| **B** | **1:7** | **1:31** | **32** |
| C | 1:4 | 1:19 | 20 |
| D | 1:9 | 1:39 | 40 |
| E | 1:11 | 1:47 | 48 |
| F | 1:19 | 1:79 | 80 |

* F-4017 RMS granularity is measured in HC-110 dilution B at 20 °C (sibling `latent-image.md`).
* "Both mixing methods provide the same photographic characteristics" (J-24).
* Reverse-engineered HC-110 clones are discussed on hobby sites (e.g. covingtoninnovations.com/hc110, photrio "New HC-110 Formula" thread). These were **not used** for any number.

---

## 4. Rodinal / ADOX Adonal / R09 (unreferenced "predicted" case)

| Source | Content | Status |
|---|---|---|
| ADOX Rodinal SDS v1.0, 6/24/2024 (US OSHA format) | 4-aminophenol (123-30-8) **3–5 %**; potassium hydroxide (1310-58-3) **1–3 %**; potassium bromide (7758-02-3) **1–3 %**. pH **> 13**. Sulfite (or metabisulfite) not listed: it is non-classified and silently omitted. | opened |
| ADOX product page | "Created by Momme Andresen in 1891 in Berlin"; produced "according to Agfa Leverkusen's latest Rodinal 2004 formula"; dilution 1+25 to 1+500; one-shot; ≥ 5 mL concentrate per 135 or 120 film | opened |
| AGFA RODINAL leaflet (a&o imaging solutions GmbH) | Dilution **1+25 or 1+50** at 20 °C; agitation continuous for the first minute, then tilting every 30 s; contains potassium hydroxide and aminophenol | opened |
| Adonal SDS 2014 (art. 105100), APH 09 SDS 2014 | label-relevant components KOH and 4-aminophenol; "mixture ... with harmless additions" | secondary |
| Andresen's 1891 patent | number and text **not found** | not found |
| Classic "Rodinal formula" | p-aminophenol hydrochloride with potassium metabisulfite, then KOH added until the p-aminophenol base redissolves as phenolate. Amounts circulate in hobby formularies only. | recalled (hobby literature); confirm in Haist, *Modern Photographic Processing* (1979) vol. 1, or Mason, *Photographic Processing Chemistry* (1975) |

* Working concentration (estimate): p-aminophenol 3–5 % × density 1.1–1.3 (assumed) gives:
  * 1+25: **1.3–2.5 g/L (12–23 mM)**;
  * 1+50: 0.65–1.25 g/L.
* KBr 1–3 % gives 0.4–1.5 g/L at 1+25 (**3.5–13 mM Br⁻**).
* Working pH: **not found**.
  * Ideal-solution sweeps over the undisclosed KOH and sulfite amounts give 11.6–12.3 at 1+25 and 11.4–12.1 at 1+50 (§5).
  * Even at 1+100, β is about 5–14 mM/pH. That is not well buffered at high dilution, which is the point of Rodinal's exhaustion and compensation behaviour.
* **Class: msds-estimated.** Validation status in D1: `predicted`.

---

## 5. Buffer chemistry

### 5.1 Method (estimate) and its known error

I solved the charge balance in ideal solution (activities = concentrations) at 25 °C, with the pKa values in §6. β = dC_base/dpH, computed numerically.

**The XTOL check shows the model's bias:**

| Case | Computed (ideal) | Measured (patent) | Offset |
|---|---|---|---|
| Patent Ex. 1, Part A in 850 mL | 11.03 | 10.19 | −0.84 |
| Patent Ex. 1, final | 8.88 | 8.20 | −0.68 |
| D-76 published (2/5/100/2), Carlton & Crabtree 1929 Table II expt. 1 (added by spot-check) | 8.94 | 8.5 (colorimetric, ±0.1) | about −0.4 |

The offset is therefore not a single constant: the two anchors differ by about 0.3 pH, beyond the ±0.1 indicator error, though the 1929 colorimetric values may themselves carry a salt error in 10 % sulfite. The 1929 borax series (no borax 8.1, 2 g/L 8.5, 6 g/L 8.9, 20 g/L 9.2) spans 1.1 pH units, while the ideal model gives only about 8.87 → 9.14 over the same range (spot-check computation), so the ideal model also underpredicts the sensitivity to borax.

The ionic strength is about 2 mol/L (0.67–0.79 M Na₂SO₃), where activity corrections are large and Davies-type equations are invalid. The `HSO₃⁻/SO₃²⁻` and `B(OH)₃/B(OH)₄⁻` apparent pKa values fall by several tenths there.

**Recommendation for D1:**
* Carry an effective (conditional) pKa for the sulfite and borate pairs, fitted to the measured anchors: XTOL 8.20 and 10.19 (patent), XTOL 8.2 ± 0.05 at SG 1.085 (J-109), and D-76 8.5 with the borax series 8.1–9.2 (Carlton & Crabtree 1929, colorimetric, lower weight).
* Alternatively, carry an explicit activity model (Pitzer or SIT) with literature parameters.
* Treat every computed pH below as an **upper bound**. Treat β as order of magnitude only.

### 5.2 Results

| Developer | Buffer pair(s) | Ideal pH (upper bound) | β (mmol L⁻¹ pH⁻¹) | Measured or reported pH |
|---|---|---|---|---|
| D-76 published, stock | HSO₃⁻/SO₃²⁻ (pKa ≈ 7.2) with B(OH)₃/B(OH)₄⁻ (9.24) and HQ/HQ⁻ (9.85) | 8.94 | 54 | **8.5** (Carlton & Crabtree 1929, colorimetric); 8.3 fresh, ~9 aged (secondary) |
| D-76 published, 1:1 | same | 8.94 | 27 | — |
| D-76 packaged-style, B₂O₃ 1 / 2 / 3 g/L added | borax + boric acid | 8.80 / 8.72 / 8.65 | 74 / 91 / 106 | not published |
| D-76d (borax 8 + boric acid 8 g/L, for comparison) | borax/boric | 8.72 | 146 | — (the 8 + 8 g/L buffered variant is described in Carlton & Crabtree 1929, p. 420) |
| XTOL Ex. 1, stock | sulfite/bisulfite + metaborate | 8.88 | 48 | **8.20 ± 0.05** |
| XTOL Ex. 1, 1:1 / 1:2 / 1:3 | same | 8.88 | 24 / 16 / 12 | not found |
| HC-110 new, dilution B (1+31) | sulfite + borate + HQ | 9.26 | 30–35 | "9.3 in solution" (SDS) |
| HC-110 new, 1+63 | same | 9.26 | 15–17 | — |
| HC-110 old | amine-H⁺/amine (DEA 8.88, MEA 9.50) with HSO₃⁻/SO₃²⁻ | **not determinable**: 6.4–7.3 from SDS mid-ranges, which is inconsistent | — | analogues 9.6–10.0 at 1+20 (patent) |
| Rodinal-type 1+25 / 1+50 / 1+100 | p-aminophenolate/phenol (10.3) and OH⁻ | 11.9–12.3 / 11.6–12.1 / 11.4–11.8 | 20–50 / 11–26 / 6–14 | concentrate > 13 (SDS) |

How dilution changes pH:
1. **Sulfite/borate developers (D-76, XTOL, new HC-110).** The ideal model predicts **no pH change** on dilution, because buffer ratios are preserved. β falls in proportion to the dilution factor. The real shift comes from activity: dilution lowers ionic strength and raises the conditional pKa toward the I→0 value. So pH should **rise** slightly (estimate: +0.1 to +0.3 for 1:1). The size has to be measured or fitted. No source gives it.
2. **Rodinal.** The pH is set by excess hydroxide and phenolate, so it **falls** with dilution, by about 0.3 per doubling (estimate).
3. **Development acid load.** Each Ag⁺ reduced releases about 1 H⁺ (HQ → Q + 2H⁺ + 2e⁻). At D_max-type silver development, a few mmol of acid per litre of film-contact volume is released. That is comparable to β of dilute Rodinal or XTOL 1:3. This is why the design spec's transported pH field matters most for diluted developers.

---

## 6. Species data

pKa values are thermodynamic or "mixed" constants at the stated temperature and ionic strength. IUPAC-dataset entries give the original measurement reference in brackets.

**Spot-check correction (citations).** The IUPAC dataset's reference codes (e.g. `B38`, `E20`, `K25`) mean different papers depending on the row's `source` column (`perrin`, `perrin_supp`, `serjeant`); the reference table has one column per source. The first version of this report resolved every code against the Perrin column, which mis-attributed every Serjeant & Dempsey and Perrin-supplement entry. The pKa **values** were right; the brackets below have been re-resolved against the correct column (fresh download of v2.4, 2026-11-21).

### 6.1 Developing agents

| Species | Formula | M (g/mol) | CAS | pKa (T, I) | Source / confidence |
|---|---|---|---|---|---|
| Hydroquinone | C₆H₆O₂ | 110.11 | 123-31-9 | pKa1 9.85 (25 °C, I 0.65 NaClO₄) / 9.91 (20–22 °C, I 0.1) / 9.82 (25 °C, I 0.1 KCl); pKa2 11.39–11.56 (20–26 °C) | IUPAC dataset, Serjeant & Dempsey entries [9.85: Baxendale & Hardy, *Trans. Faraday Soc.* 49, 1140 (1953); 9.91: Staude & Teupel, *Z. Elektrochem.* 61, 181 (1957) and Sunkel & Staude, *Ber. Bunsenges.* 72, 567 (1968); 9.82: Epstein et al., *JACS* 86, 3075 (1964)] — opened; citations corrected by spot-check (previously given as Beale & Liberman 1950 and Elderfield & Rubin 1953, which are the Perrin-column papers for the same codes) |
| | | | | **Disagreement:** HSDB in PubChem quotes pKa 9.96 (AHFS 1990) and **10.85** at 25 °C (Pearce & Simkins, *Can. J. Chem.* 46, 241 (1968)). 10.85 is an outlier against the IUPAC set; use 9.85–9.9. | PubChem/HSDB — opened |
| Metol base, 4-(methylamino)phenol | C₇H₉NO | 123.15 | 150-75-4 | pKaH (ammonium) 5.9 (25 °C, 0.2 M buffers) | IUPAC dataset [Fieser, *JACS* 52, 4915 (1930)] — opened. **Phenol pKa not found**; by analogy with p-aminophenol, 10.1–10.4 (estimate) |
| Metol sulfate (Elon) | (C₇H₉NO)₂·H₂SO₄ | 344.39 | 55-55-0 | — | PubChem — opened |
| p-Aminophenol | C₆H₇NO | 109.13 | 123-30-8 | pKaH 5.29–5.70 (21–25 °C); phenol pKa 10.30–10.46 (22–25 °C, I ≤ 0.1) | IUPAC dataset [5.48: de Courville, *Compt. Rend. C* 262, 1196 (1966); 10.46: Epstein et al., *JACS* 86, 3075 (1964); 10.44: van Hooidonk & Ginjaar, *Rec. Trav. Chim.* 86, 449 (1967); 5.7: Fieser, *JACS* 52, 4915 (1930); 5.5 / 10.3: Eggers et al., *Mitt. Forschungslab. Agfa Leverkusen-München* 1, 95 (1955)] — opened; 5.48 / 10.46 / 10.44 citations corrected by spot-check (previously Datta & Grzybowski 1961, Elderfield & Rubin 1953, Hirt et al. 1959) |
| p-Aminophenol hydrochloride | C₆H₈ClNO | 145.59 | 51-78-5 | — | PubChem — opened |
| Phenidone (1-phenyl-3-pyrazolidinone) | C₉H₁₀N₂O | 162.19 | 92-43-3 | pKaH −0.3 / −0.4 (protonation); **acid pKa: 9.48 (T not stated) vs 7.57 (20 °C) — unresolved conflict** | IUPAC dataset [9.48 and −0.3: Cauzzo & Mazzucato, *Atti Accad. Naz. Lincei* 34, 539 (1963); −0.4: Zaitsev et al., *Zh. Obshch. Khim.* 39, 1835 (1969); **7.57**: Tabak, Grandberg & Kost, *Tetrahedron* 22, 2703 (1966), "extrapolated from values in ethanol-water", Uncertain] — opened. **Spot-check correction:** the first version called 9.48 a single entry and omitted the 7.57 entry. At XTOL's pH 8.2 the two values give an ionised fraction of about 5 % vs about 80 %, which matters for any superadditivity model of Phenidone/Dimezone S. Keep the pyrazolidone acid pKa `fitted` or `assumed` until a direct aqueous measurement is found. |
| 4-Methyl-1-phenyl-3-pyrazolidinone (old HC-110) | C₁₀H₁₂N₂O | 176.21 | 2654-57-1 | pKaH −0.5 (25 °C) | IUPAC dataset [Zaitsev 1969] — opened. Acid pKa not found. |
| Dimezone (4,4-dimethyl-1-phenyl-3-pyrazolidinone) | C₁₁H₁₄N₂O | 190.24 | 2654-58-2 | not found | PubChem — opened |
| **Dimezone S** (4-hydroxymethyl-4-methyl-1-phenyl-3-pyrazolidinone) | C₁₁H₁₄N₂O₂ | 206.24 | 13047-13-7 | **not found** | PubChem — opened |
| L-Ascorbic acid | C₆H₈O₆ | 176.12 | 50-81-7 | pKa1 4.04 / pKa2 11.34 (25 °C, I 0.1 KNO₃) | IUPAC dataset, Serjeant & Dempsey entries [Khan & Martell, *JACS* 89, 7104 (1967) and *JACS* 91, 4668 (1969)] — opened; **citation corrected by spot-check** (previously given as King 1945/1951, the Perrin-column papers for codes K25/K26). The other pairs under the same label are the same compound under other conditions: **4.49 / 12.72 at 0.4 °C**, **3.57 / 11.98 at 40 °C**, and **4.56 / 12.22 at 25 °C in D₂O** (all I 0.1 KNO₃). The earlier remark that they "probably belong to other compounds" was wrong. Commonly quoted values: 4.10 / 11.6–11.8 (CRC; recalled). |
| Erythorbic (D-isoascorbic) acid | C₆H₈O₆ | 176.12 | 89-65-6 | not in the IUPAC set (the ascorbic rows above are all L-ascorbic acid); expect within 0.1–0.2 of ascorbic | recalled; no source opened |
| Sodium erythorbate (isoascorbate) | C₆H₇NaO₆ (·H₂O) | 198.11 anhydrous; 216.12 monohydrate | 6381-77-7 (PubChem resolves it to the anhydrous salt; 7378-23-6 is another registry entry) | — | PubChem — opened. Hydrate M is my computation. |
| Sodium ascorbate | C₆H₇NaO₆ | 198.11 | 134-03-2 | — | recalled CAS |

### 6.2 Preservatives, buffers, alkalis, halides, sequestrants, amines

| Species | Formula | M (g/mol) | CAS | pKa (25 °C unless stated) | Source / confidence |
|---|---|---|---|---|---|
| Sodium sulfite | Na₂SO₃ | 126.04 | 7757-83-7 | SO₂(aq)/HSO₃⁻ **1.86** (K_A = 0.0139 mol/kg); **HSO₃⁻/SO₃²⁻ 7.19** (K_B = (6.5 ± 0.5)×10⁻⁸ mol/kg; 7.17 from the network-consistent 6.72×10⁻⁸); ΔH°(B) = −3.65 kJ/mol; all at 298.15 K, I → 0 | Goldberg & Parker, *J. Res. NBS* 90, 341 (1985), §4 and Table 3 — **opened by spot-check** (was recalled) |
| Potassium sulfite | K₂SO₃ | 158.26 | 10117-38-1 | as above | PubChem (M) — opened |
| Sodium metabisulfite | Na₂S₂O₅ | 190.11 | 7681-57-4 | 2 HSO₃⁻ ⇌ S₂O₅²⁻ + H₂O, **K = 0.032 ± 0.01 kg/mol** (298.15 K, I → 0; from Connick et al.), ΔH° = −4.6 ± 4.0 kJ/mol; 1 mol gives 2 mol HSO₃⁻ | K: Goldberg & Parker 1985, process (C) — opened by spot-check. **Correction:** the first version gave K ≈ 0.07 L/mol (recalled); that value was not found in Goldberg & Parker (it may be a value at non-zero ionic strength). M: PubChem — opened |
| Potassium metabisulfite | K₂S₂O₅ | 222.33 | 16731-55-8 | — | PubChem — opened |
| Sulfur dioxide | SO₂ | 64.07 | 7446-09-5 | — | — |
| Boric acid | H₃BO₃ | 61.83 | 10043-35-3 | **9.235 (25 °C), 9.278 (20 °C)**, from pK = 2237.94/T + 0.016883 T − 3.305 (0–60 °C), i.e. dpKa/dT ≈ −0.008 K⁻¹; polyborates form above about 0.025 M total B | pKa: Manov, DeLollis & Acree, *J. Res. NBS* 33, 287 (1944), abstract equation — **opened by spot-check** (was recalled). Polyborate threshold: still recalled (Ingri 1960s) |
| Borax decahydrate | Na₂B₄O₇·10H₂O | 381.37 | 1303-96-4 | dissolves to 2 B(OH)₃ + 2 B(OH)₄⁻ (pH ≈ pKa) | PubChem (M) — opened |
| Borax pentahydrate | Na₂B₄O₇·5H₂O | 291.29 | 12179-04-3 | same | PubChem — opened |
| Sodium metaborate tetrahydrate ("8 mol") | NaBO₂·4H₂O | 137.86 | 10555-76-7 | gives B(OH)₄⁻ (a base) | PubChem — opened |
| Diboron trioxide | B₂O₃ | 69.62 | 1303-86-2 | gives 2 B(OH)₃ | PubChem — opened |
| Sodium carbonate / bicarbonate | Na₂CO₃ / NaHCO₃ | 105.99 / 84.01 | 497-19-8 / 144-55-8 | 6.35 / 10.33 | recalled (Harned & Davis 1943; Harned & Scholes 1941); **unverified** by spot-check |
| Sodium / potassium hydroxide | NaOH / KOH | 40.00 / 56.11 | 1310-73-2 / 1310-58-3 | pKw 14.00 (25 °C), 14.17 (20 °C) | recalled |
| Potassium bromide | KBr | 119.00 | 7758-02-3 | — | — |
| Sodium bromide | NaBr | 102.89 | 7647-15-6 | — | — |
| EDTA | C₁₀H₁₆N₂O₈ | 292.24 | 60-00-4 | 1.99, 2.67, 6.16, 10.26 (20 °C, I 0.1 KCl) | IUPAC dataset, Serjeant & Dempsey entry [Schwarzenbach & Ackermann, *Helv. Chim. Acta* 30, 1798 (1947)] — opened; citation corrected by spot-check (previously Schwarzenbach et al., *Helv. Chim. Acta* 20, 498 (1937)) |
| DTPA (pentetic acid) | C₁₄H₂₃N₃O₁₀ | 393.35 | 67-43-6 | about 1.8, 2.6, 4.3, 8.6, 10.5 (25 °C, I 0.1) | recalled (Martell & Smith, *Critical Stability Constants*); not in the IUPAC high-confidence set |
| DTPA pentasodium | C₁₄H₁₈N₃Na₅O₁₀ | 503.26 | 140-01-2 | — | M computed; PubChem lookup failed |
| HEDTA | C₁₀H₁₈N₂O₇ | 278.26 | 150-39-0 | about 2.6, 5.4, 9.8 | recalled (Martell & Smith) |
| HEDTA trisodium (Dissolvine H-40 active) | C₁₀H₁₅N₂Na₃O₇ | 344.20 | 139-89-9 | — | M computed; CAS from PSI SDS |
| Diethanolamine | C₄H₁₁NO₂ | 105.14 | 111-42-2 | **8.883 (25 °C)**; 9.133 (15), 9.005 (20), 8.759 (30) → dpKa/dT ≈ −0.025 K⁻¹ | IUPAC dataset, "Reliable" [Bower, Robinson & Bates, *J. Res. NBS* 66A, 71 (1962)] — opened |
| Ethanolamine (MEA) | C₂H₇NO | 61.08 | 141-43-5 | **9.498 (25 °C)**; 9.647 (20) → about −0.030 K⁻¹ | IUPAC dataset, "Reliable" [Bates & Pinching, *J. Res. NBS* 46, 349 (1951)] — opened |
| Diethylene glycol | C₄H₁₀O₃ | 106.12 | 111-46-6 | non-ionizing | PubChem — opened |
| Amine–SO₂ adducts (DEA·SO₂, MEA·SO₂) | R₂NH·SO₂, i.e. R₂NH₂⁺ HSO₃⁻ in water | DEA·SO₂ 169.2; MEA·SO₂ 125.2 (computed) | not registered as such | In water they behave as amine-H⁺ + HSO₃⁻: the sulfite pair (7.2) plus the amine pair (8.9 / 9.5) | DE 1,118,005 (opened) describes the adduct and says it forms "by adding gaseous SO₂ to 2-aminoethanol". Stoichiometry (1:1 vs 2:1) is **not stated** clearly. |

### 6.3 Aqueous diffusion coefficients

**Measured values for the organic developing agents were not found** in any source I could open.

* Handbook value for hydroquinone: 0.83×10⁻⁵ cm²/s at 20 °C. It is itself a Hayduk–Laudie **estimate** (HSDB, secondary); my calculation reproduces it exactly.
* A measured hydroquinone D vs T (25–375 °C) exists (OSTI record 5469741, J. Phys. Chem.), but the 25 °C value was not seen.
* Estimates below (**estimate**):
  * Wilke–Chang (φ = 2.6, η_water = 1.002 cP at 20 °C, 0.890 cP at 25 °C);
  * Hayduk–Laudie (D = 13.26×10⁻⁵ η⁻¹·¹⁴ V⁻⁰·⁵⁸⁹);
  * Le Bas molar volume V: C 14.8, H 3.7, O 7.4, N(sec) 12.0, N(prim) 10.5, N(tert) 10.8, 6-ring −15.0, 5-ring −11.5 cm³/mol (Le Bas increments as tabulated in Poling, Prausnitz & O'Connell, *The Properties of Gases and Liquids*; **recalled**).

| Species (neutral form) | V_LeBas (cm³/mol) | D, 20 °C (WC / HL) ×10⁻⁵ cm²/s | D, 25 °C (WC / HL) ×10⁻⁵ cm²/s |
|---|---|---|---|
| Hydroquinone | 110.8 | 0.88 / 0.83 | 1.01 / 0.95 |
| p-Aminophenol | 117.6 | 0.85 / 0.80 | 0.97 / 0.91 |
| Metol (base) | 141.3 | 0.76 / 0.72 | 0.87 / 0.82 |
| Phenidone | 173.9 | 0.67 / 0.63 | 0.77 / 0.73 |
| 4-Methyl-phenidone | 196.1 | 0.62 / 0.59 | 0.72 / 0.68 |
| Dimezone S | 225.7 | 0.57 / 0.54 | 0.66 / 0.62 |
| Ascorbic / isoascorbic acid | 151.3 | 0.73 / 0.69 | 0.84 / 0.79 |
| Diethanolamine | 126.7 | 0.81 / 0.76 | 0.93 / 0.87 |
| Ethanolamine | 73.4 | 1.13 / 1.05 | 1.29 / 1.21 |
| Diethylene glycol | 118.4 | 0.85 / 0.80 | 0.97 / 0.91 |

Inorganic ions, infinite dilution, 25 °C. D = RTλ/(|z|F²) from limiting equivalent conductivities λ (S cm² equiv⁻¹). The λ values are **recalled** from the CRC Handbook table "Ionic conductivity and diffusion at infinite dilution"; confirm there. **Unverified (spot-check):** no copy of that table on a permitted host could be opened (the CRC/HBCP site is paywalled; the copies found are on course pages). A search-engine summary of one such copy agreed for Br⁻ (78.1; 2.080) and HSO₃⁻ (58; 1.545), which is secondary at best. The D column is internally consistent with D = RTλ/(|z|F²).

| Ion | λ | D (×10⁻⁵ cm²/s) |
|---|---|---|
| Na⁺ | 50.08 | 1.334 |
| K⁺ | 73.48 | 1.957 |
| Br⁻ | 78.1 | 2.080 |
| OH⁻ | 198.0 | 5.27 |
| ½SO₃²⁻ | 72 | 0.959 |
| HSO₃⁻ | 58 | 1.54 |
| ½CO₃²⁻ | 69.3 | 0.923 |
| HCO₃⁻ | 44.5 | 1.185 |
| ½SO₄²⁻ | 80.0 | 1.065 |
| B(OH)₄⁻ | not recalled reliably | not found |

Corrections for the model:
1. **Viscosity.** In 0.7–0.8 M sulfite the viscosity is higher than water (about 1.2–1.4 cP, **estimate**; not found). Scale all D by η_water/η_solution (Stokes–Einstein).
2. **Charge and ion pairing.** The anionic forms (HQ⁻, ascorbate⁻) diffuse with electro-neutrality coupling to Na⁺/K⁺. D1 can use the neutral-form estimates with a fitted gelatin tortuosity factor.
3. **Temperature.** D ∝ T/η: about +2.4 %/K near 20 °C (estimate from the table above).

### 6.4 Redox potentials

**No tabulated developer potentials could be opened.** The values below are **recalled** and must be confirmed. **Unverified (spot-check):** none of these could be opened on a permitted host either. For AgBr/Ag, older literature spans about 0.0710–0.0713 V (e.g. Harned, Keaton & Donelson 0.07103 V, as quoted in a thesis; secondary), so carry 0.071 ± 0.0003 V rather than four significant figures. The standard compilation is James (ed.), *The Theory of the Photographic Process*, 4th ed. (1977), the chapter on developing agents. Fieser, *JACS* 52, 4915 (1930) is the measurement source for the p-aminophenol and metol potentials (the IUPAC set uses that paper's pH-dependence data for their pKaH).

| Couple | E° or E°′ (V vs SHE, 25 °C) | Confidence |
|---|---|---|
| Ag⁺/Ag | +0.7996 | recalled (high) |
| AgBr/Ag, Br⁻ | +0.0713 | recalled (high) |
| p-Benzoquinone / hydroquinone (2e⁻, 2H⁺) | +0.699; E = 0.699 − 0.0592 pH below pKa1 | recalled (high) |
| Dehydroascorbic acid / ascorbate (pH 7) | about +0.06 (E°′) | recalled; confirm in Borsook & Keighley, *PNAS* 19, 875 (1933), or Creutz, *Inorg. Chem.* 20, 4449 (1981) |
| Ascorbyl radical / ascorbate monoanion (pH 7) | about +0.28 | recalled (Williams & Yandell, *Aust. J. Chem.* 35, 1133 (1982)) |
| p-Quinoneimine / p-aminophenol | about +0.73 (E°) | recalled (Fieser 1930) |
| Metol, phenidone, Dimezone S | **not found** | — |

For D1 these potentials only fix the sign and pH trend of the driving force. The design spec already makes the rate constants `fitted`.

---

## 7. Disagreements and open problems (flags)

1. **HC-110 version:**
   * The Kodak reference data were almost certainly made with the pre-2019 formulation:
     * J-24 (Dec 2017) predates it;
     * the F-4016/F-4017 HC-110 tables and CI curves carry over from the 2007–2018 editions;
     * the F-4017 granularity in "HC-110 dilution B" appears in the Oct 2021 edition, but I found no sign that it was re-measured. This is an inference; confirm by comparing the 2016 and 2021 editions.
   * The old formula is documented only by search snippets of SDS copies (DEA/SO₂/MEA; 4-methyl-phenidone; KBr).
   * The new formula (K₂SO₃/HQ/borax/Dimezone S/NaBr) is the opened source.
   * **Keep two separate YAML files.** Calibrate the HC-110 pairs against `hc-110-pre2019`.
2. **D-76 version:** the film calibration curves come from **packaged** D-76 (B₂O₃ and DTPA present, amounts unknown), not the published formula. Record it as msds-estimated with B₂O₃ in 0–3 g/L. The ideal pH sensitivity is −0.3 over that range; β roughly doubles.
3. **XTOL Dimezone S:** < 0.01 % (PSI 2024) vs 0.1 to < 1 % (Kodak 2015, secondary) vs 1.3 % of Part A (patent Ex. 1).
4. **XTOL Part A borate:** metaborate (patent) vs borax pentahydrate (both SDS).
5. **Hydroquinone pKa1:** 9.82–9.91 (IUPAC set) vs 9.96 and 10.85 (HSDB).
6. **Ascorbic acid pKa:** resolved by spot-check. The extra pairs in the IUPAC set (Khan & Martell 1967/1969, not King) are 0.4 °C, 40 °C and D₂O values; use 4.04 / 11.34 at 25 °C, I 0.1.
7. **D-76 sulphite in the Cinematographic Annual:** 14 oz (1930) vs 11 oz (1931) for formula "A". Carlton & Crabtree 1929 (100 g/L) supports 14 oz; 11 oz is the outlier.
8. **New HC-110 relative density 1.07** looks too low for the stated composition.
9. **SDS ranges are not closed formulas.** Unclassified ingredients are silently omitted:
   * sulfite in Rodinal;
   * erythorbate in some XTOL Part B sheets;
   * possibly the undisclosed alkali in old HC-110.
10. **pH model bias:** the ideal-solution pH runs 0.7–0.8 high on the XTOL anchors and about 0.4 high on the D-76 anchor (8.94 vs 8.5, colorimetric, 1929), all at I ≈ 2 M. The offset is formula-dependent; fit it, do not apply one constant.
11. **Patent Table I films:** "TRI-X PAN" and "TMAX 100" in US 5,756,271 Table I (1996) are the **pre-2007 emulsions**. The current 400TX and the 2002+ T-MAX 100 differ. Use the table only for trends (XTOL-type vs D-76 ratio), not as reference curves.
12. **Pyrazolidone acid pKa:** Phenidone 9.48 (Cauzzo & Mazzucato 1963) vs 7.57 (Tabak et al. 1966, extrapolated from ethanol-water); Dimezone S not found. Unresolved; keep `fitted`.
13. **IUPAC reference codes:** must be resolved against the column named by each row's `source`. The first version of this report did not, and mis-cited the HQ, p-aminophenol, ascorbic acid and EDTA entries (values unaffected; fixed in §6).

**US 5,756,271 Table I, selected rows** (opened). 20 °C; time in minutes; N₂-burst agitation.

| Film | Developer | Log E at ANSI 0.62 | EI | t(CI 0.56) min | t(CI 0.82) min | Dmin at CI 0.56 |
|---|---|---|---|---|---|---|
| TMAX 100 | Example 1 (XTOL-type) | −2.18 | 125 | 6.6 | 10.5 | 0.20 |
| TMAX 100 | D-76 | −2.10 | 100 | 6.7 | 10.2 | 0.22 |
| TMAX 100 | T-MAX Dev. | −2.26 | 160 | 6.7 | 11.6 | 0.23 |
| TRI-X PAN | Example 1 | −2.80 | 500 | 7.0 | 11.0 | 0.27 |
| TRI-X PAN | D-76 | −2.65 | 400 | 6.5 | 12.0 | 0.26 |
| TRI-X PAN | T-MAX Dev. | −2.77 | 500 | 6.7 | 12.0 | 0.32 |

---

## 8. Spot-check log (2026-11-21, independent reviewer)

Each item was re-checked against the source named, re-opened where possible (fresh downloads under `.cache/datasheets/spotcheck_devf/`, SDS PDFs re-extracted with PyMuPDF rather than trusting the cached text). Corrections have been applied in place above.

| # | Claim | Verdict | Evidence / correction |
|---|---|---|---|
| 1 | D-76 published formula 2/100/5/2 g/L; Cinematographic Annual vol. 1 pp. 119–120 (120 gr / 14 oz / 300 gr / 120 gr per gal) | confirmed | Annual text re-read. Metric form now opened in Carlton & Crabtree, *Trans. SMPE* 13 (1929), Table II expt. 1 and p. 420. The J-1 attribution is dropped; the formula's first publication is the paper's ref. 1, *Eastman Duplicating Film* (1927), not opened. |
| 2 | D-76 working pH ≈ 8.3 fresh → 9.0 aged, no Kodak primary | wrong | Kodak primary exists: pH **8.5** (Carlton & Crabtree 1929, Table II, colorimetric ±0.1, temperature not stated). Offset vs the ideal 8.94 is about −0.4, not −0.7 to −0.8. |
| 3 | Packaged D-76 SDS (PSI 5160304 v2.0, 7/24/2024) composition | confirmed | PDF re-extracted: Na₂SO₃ 80–95, HQ 3–5, borax pentahydrate 1–3, metol sulphate 1–3, B₂O₃ 1–3, DTPA < 1 % w/w. |
| 4 | XTOL patent Ex. 1 kit, pH 10.19 / 8.20, composition ranges, J-109 SG 1.085 ± 0.003 at pH 8.2 ± 0.05 | confirmed | US 5,756,271 re-fetched from patents.google.com; J-109 (Feb 2018) p. 1 text. Ratio range in the patent is "at least about 10:1, preferably 20:1 to 50:1". |
| 5 | Patent Table I rows (TMAX 100, TRI-X PAN) | confirmed | Re-fetched; all six rows match. |
| 6 | New HC-110 SDS (PSI 1058692 v5.0) composition, pH 9.0 / 9.3, RD 1.07 | confirmed, with omission | PDF re-extracted; all listed values match. **Pyrocatechol (120-80-9) < 0.05 %** was omitted; added. |
| 7 | Old HC-110 SDS composition (DEA 30–35, SO₂ 15–20, HQ 5–10, …) | unverifiable | Only dealer/university copies exist; no Kodak/Kodak Alaris/PSI host or web.archive.org copy found. Marked unverified in §0 and §3.1. |
| 8 | HC-110 dilutions (J-24 Dec 2017 p. 2) | confirmed | Stock A 1:3, B 1:7, C 1:4, D 1:9, E 1:11, F 1:19; concentrate A 1:15, B 1:31, C 1:19, D 1:39, E 1:47, F 1:79. |
| 9 | Hydroquinone pKa1 9.85 / 9.82–9.91, pKa2 11.39–11.56, and their citations | wrong (citations only) | Values confirmed in the IUPAC set; sources are Baxendale & Hardy 1953 (9.85), Staude & Teupel 1957 / Sunkel & Staude 1968 (9.91), Epstein et al. 1964 (9.82), not Beale & Liberman / Elderfield & Rubin. |
| 10 | p-Aminophenol, metol, Phenidone pKa | wrong (citations and an omission) | p-aminophenol values right but 5.48 / 10.46 / 10.44 are de Courville 1966 / Epstein 1964 / van Hooidonk & Ginjaar 1967; pKaH range widens to 5.29–5.70 (Fieser 5.7). Metol 5.9 (Fieser 1930) confirmed. Phenidone has a second acid pKa entry **7.57** (Tabak et al. 1966) conflicting with 9.48. |
| 11 | Ascorbic acid pKa 4.04 / 11.34 (King 1945/1951); other pairs unresolvable | wrong (citation and interpretation) | Values confirmed; source is Khan & Martell 1967/1969. The other pairs are 0.4 °C, 40 °C and D₂O data. |
| 12 | Amine pKa: DEA 8.883 (25 °C), 9.005 (20 °C); MEA 9.498, 9.647 | confirmed | IUPAC set, "Reliable", Bower, Robinson & Bates 1962 and Bates & Pinching 1951 (Perrin column, correctly resolved). |
| 13 | Inorganic pKa: HSO₃⁻ 7.2, boric 9.24, carbonate 6.35 / 10.33 (recalled) | confirmed in part | Sulfite pKa1 1.86 and pKa2 7.19 (Goldberg & Parker 1985) and boric acid 9.235 at 25 °C / 9.278 at 20 °C (Manov et al. 1944) now opened. Carbonate still unverified. Side correction: pyrosulfite K is 0.032 ± 0.01 kg/mol (G&P), not 0.07. |
| 14 | Ion diffusion coefficients at infinite dilution (CRC, recalled) | unverifiable | No permitted-host copy; search-summary agreement for Br⁻ and HSO₃⁻ only. Marked unverified in §6.3. |
| 15 | Redox potentials (Ag⁺/Ag, AgBr/Ag, quinone/HQ, ascorbate, quinoneimine; recalled) | unverifiable | Not opened on a permitted host. AgBr/Ag literature spread 0.0710–0.0713 V noted. Marked unverified in §6.4. |

Also corrected while re-resolving the IUPAC citations: EDTA → Schwarzenbach & Ackermann, *Helv. Chim. Acta* 30, 1798 (1947).
