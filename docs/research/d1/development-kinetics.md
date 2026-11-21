# D1 research: development kinetics and its numerical treatment

Scope: inputs for stage 3 of the physics chain (chemical development as a reaction-diffusion system) for the D1 black-and-white chain. Films: Kodak Tri-X 400 and T-Max 100. Developers: D-76, XTOL, HC-110, and Rodinal as a predicted case. Research date: 2026-11-21. Everything here must be checked again before a value goes into `data/` with a provenance tag.

Sibling reports: `datasheets.md` (Kodak publication inventory, split feasibility) and `latent-image.md` (latent-image constants). This report does not repeat them.

## Confidence labels

| Label | Meaning |
|---|---|
| **opened-source** | I opened the primary document (full text, a publisher abstract page, or a rendered page image) and read the value myself. |
| **secondary-source** | Seen only in a search-result abstract or snippet, a catalogue record, or a later document that cites the original. The original was not opened. |
| **recalled** | From memory of the literature; not opened in this session. Each one says what would confirm it. |
| **estimate** | Computed here from opened values plus standard constants. The arithmetic is shown. |

## Access log

| Source family | Result |
|---|---|
| AIP (J. Chem. Phys.) abstract pages via the publisher's Silverchair host | **opened** (abstracts only) for James 1944 and James 1946 |
| Hokkaido University HUSCAP (publisher of J. Fac. Sci. Hokkaido Univ. Ser. II) | **full text opened**: Sugai & Furuichi 1953; Sugai 1953 |
| J-STAGE (J. Soc. Photogr. Sci. Technol. Japan, free access) | **full text opened**: Iwano 1970 (Part 2) and Iwano 1972 (Part 4). Tables read from rendered page images |
| NBS Journal of Research (nvlpubs.nist.gov) | **full text opened**: Carroll & Hubbard 1930 (RP183) and 1934 (RP659) |
| Optica Applicata (Wrocław digital library) | **full text opened**: Rajkowski & Nowak 1998 |
| Chimia (journal site) | **full text opened**: Eggers 1961 (colour development review) |
| European Patent Office, Google Patents | **opened**: EP 0 573 700 A1 (Agfa), US 5,315,338 A, US 5,837,437 A (Kodak) |
| SMPTE journal abstract pages | **opened** (abstracts only): Ives & Jensen 1943; Snyder 1966 |
| Kodak J-78, J-109 (1998, 2018), J-24, J-86, F-4016, F-4017 | **opened** via the text already extracted to `.cache/datasheets/ds_txt/` |
| Taylor & Francis (J. Photogr. Sci., Imaging Sci. J.) | HTTP 403. Abstracts seen only through search snippets → secondary |
| IS&T library (Photogr. Sci. Eng., J. Imaging Sci. Technol.) | HTTP 403 → secondary |
| ScienceDirect (J. Colloid Sci.) | HTTP 403 → secondary |
| DTIC report ADA151831 | HTTP 403 |
| RIT thesis repository | not fetched (PDFs returned 403 in the sibling job). Abstracts seen through search snippets → secondary |
| archive.org (Nietz 1922, Mees), lantern.mediahist.org, mci.si.edu, scilit | not on the fetch whitelist. Used only as pointers from search results |
| James (ed.) *TPP* 4th ed.; Mason *Photographic Processing Chemistry*; Haist *Modern Photographic Processing* | **not opened**. No accessible full text |

Fetched files are kept under `.cache/datasheets/lit/` (git-ignored), one folder per document.

---

## 1. Rate laws for chemical development

### 1.1 Mechanism: what the primary literature I could open says

| Finding | Source | Confidence |
|---|---|---|
| Development rate of AgBr grains varies approximately as **[HQ]^0.6 at pH 8.0–8.9**, and as **[catechol]^0.56 at pH 7.7**. James reads the fractional order as adsorption of the developing agent by AgBr or by the Ag/AgBr interface before reaction. He states that this contradicts the "electrode" mechanism and supports direct reaction at the triple interface. For silver ions from soluble salts (physical development), the silver-catalysed rate is directly proportional to [HQ]. | T.H. James, "The Site of Reaction in Direct Photographic Development", J. Chem. Phys. 14, 536–539 (1946), DOI 10.1063/1.1724190, abstract | opened-source (abstract) |
| Negatively charged developer ions approaching an AgBr grain meet **two charge barriers**: one from the gelatin and one from adsorbed bromide. Neutral salt speeds development by negative ions mainly by lowering the gelatin barrier; its effect on the bromide barrier is minor. Changes to the bromide barrier account for the **induction period** and for the accelerating action of phenosafranin. | T.H. James, "The Charge Effect in Photographic Development: Role of the Gelatin", J. Chem. Phys. 12, 453–457 (1944), abstract | opened-source (abstract) |
| For hydroquinone at pH 9.25–10.35 (and chlorohydroquinone at pH 8.33–9.16), the gross development rate correlated log-log with both the monovalent and the divalent ion concentration, but **depended only on the divalent ion**, with or without a quaternary-salt accelerator. Diffusion effects mattered in the accelerated developers. Caveat from the same abstract: the thesis's own hypothesis could not be tested, because "separation of the contributions proved not to be possible"; the divalent-only conclusion rests on the log-log correlation. | R.W. Eskes, "Correlation of development acceleration by decyltrimethylammonium bromide with the ionization states of hydroquinone and of chlorohydroquinone", thesis, RIT School of Photographic Arts and Sciences (1982), advisor B.H. Carroll, https://repository.rit.edu/theses/7317 | opened-source (abstract; spot-check 2026-11-21 opened the RIT repository record) |
| The active ascorbate species are **"likely"** (the patent's hedge) the **mono- and dianion**, and their concentration is pH dependent. Both aerial oxidation and development lower the pH of an ascorbate developer. In a polyhydroxybenzene lith developer, aerial oxidation raises the pH while development lowers it. | Agfa-Gevaert, EP 0 573 700 A1 (publ. 15 Dec 1993), p. 2, l. 41–48, citing R.G. Willis & R.B. Pontius, Photogr. Sci. Eng. 14(6), 384 (1970) | opened-source (patent); Willis & Pontius itself secondary |
| Jaenicke treated AgCl/AgBr/AgI as ionic, electronic and mixed conductors and analysed open-circuit dissolution with a mixed-potential condition (equal anion and cation currents). A 1972 paper "Some Basic Concepts of Photographic Development" (author not confirmed in the snippet) argues the electrode theory is valid also for the primary stages. Jaenicke's fuller statement is cited as his review in *Advances in Electrochemistry and Electrochemical Engineering* vol. 10 (1977). | W. Jaenicke, Croat. Chem. Acta (1972) conference paper (hrcak.srce.hr); citation trail | secondary-source |
| Development rate and latent-image size: in a **sulfur-plus-gold-sensitised** emulsion (primitive and S-only emulsions were also measured) a **fast** development component corresponds to grains whose latent-image centres hold **≥ 5 Ag atoms**, explained by the electrode mechanism. A **slow** component corresponds to centres of **< 5 atoms**, explained by a one-by-one reaction mechanism with a long induction period per grain. Centre sizes were estimated by applying a nucleation-and-growth simulation to measured quantum sensitivities. | **K.** Ohzeki (Katsuhisa Ohzeki per Crossref; the initial "T." was wrong), "Rate of Development and Size of Latent Image Centers in Chemically Sensitized Emulsions", J. Imaging Sci. Technol. 40(6), 591–598 (1996), DOI 10.2352/J.ImagingSci.Technol.1996.40.6.art00016. Development was followed by optical density at 1090 nm | opened-source (abstract; spot-check 2026-11-21 read it from the web.archive.org copy of the IS&T page, snapshot 2026-06-13; full text not opened) |

Reading for D1. The literature splits on mechanism. James (1946) found adsorption-controlled direct reaction. Jaenicke and Ohzeki (1996) support an electrode (mixed-potential) mechanism, at least once the silver centre is large enough. A rate law of **Langmuir–Hinshelwood form on the grain or speck surface** is consistent with both views, and its constants are fitted in either case. It has adsorption of the active ionised developer species competing with Br⁻, and a charge-barrier factor. The electrochemical form is an alternative parameterisation, not an extra degree of freedom (§1.10).

### 1.2 Active ionised forms and pH dependence, per agent

| Agent | Active form(s) | Evidence | Confidence |
|---|---|---|---|
| Hydroquinone | **Dianion HQ²⁻** dominates; monoanion contribution small or nil | Eskes 1982 (rate depends only on the divalent ion, pH 9.25–10.35; but the abstract also says the mono- and divalent contributions could not be separated); James's charge-barrier work (the doubly charged ion has the long induction period) | opened-source (Eskes abstract); treat the "monoanion nil" part as a fitted choice |
| Hydroquinone pKa | pKa₁ ≈ 9.85–10.0; pKa₂ ≈ 11.4–12.0, depending on the compilation | Search-result summaries of compilations: EPA/WHO-IPCS 9.9 and 11.6 at 25 °C; Dean (1973) 10.0 and 12.0; a Russian encyclopaedia 9.85 and 11.4 | secondary-source. Spread of 0.6 in pKa₂. Pick one compilation and tag `literature`. It matters because at pH 8–10 the HQ²⁻ fraction is roughly ∝ 10^(2·pH − pKa₁ − pKa₂) |
| Ascorbate | Monoanion HA⁻ and dianion A²⁻ | Agfa EP 0 573 700 citing Willis & Pontius 1970 | secondary-source |
| Ascorbic acid pKa | pKa₁ ≈ 4.1–4.2, pKa₂ ≈ 11.6–11.8 | recalled (standard tables) | recalled. Confirm in CRC Handbook or IUPAC dissociation-constant tables. At XTOL's pH 8.2 (J-109) ascorbate is almost all HA⁻, and the A²⁻ fraction is ~10⁻³·⁵ |
| Metol (N-methyl-p-aminophenol) | Recalled: the uncharged free base and the phenolate; it is a weakly charged or neutral developer with a short induction period, the opposite of HQ²⁻ in the charge-barrier picture | recalled (James *TPP* 4th ed., development chapters; Mason) | recalled. Confirm the pKa of the ammonium (≈ 5–6) and the phenol (≈ 10) in a primary compilation |
| Phenidone (1-phenyl-3-pyrazolidinone) and Dimezone S | Recalled: the anion, formed by deprotonation of the pyrazolidinone N–H, is the active electron donor. A one-electron oxidation gives a coloured, reversible free radical (Axford & Kendall 1954, secondary below) | Axford & Kendall, J. Photogr. Sci. 2, 1–7 (1954) (abstract via search: "a coloured reversible mesomeric free-radical formed by the loss of an electron from the pyrazolidone ion") | secondary-source for the radical; pKa **not found** |
| p-Aminophenol (Rodinal) | Recalled: the phenolate anion; it needs high pH (Rodinal is strongly alkaline, with little sulfite) | recalled | recalled. Confirm pKa in a compilation |

Thesis-level supporting evidence (secondary-source, RIT abstracts via search):
- Brainerd (1974): for all developers tested, the bromide effect **decreased with increasing pH**.
- Durbin (1976): as bromide increased, the induction period increased and the maximum rate decreased.

### 1.3 Order of reaction

| Quantity | Value | Source | Confidence |
|---|---|---|---|
| Order in total HQ, AgBr grains, pH 8.0–8.9 | ≈ 0.6 | James 1946 | opened-source (abstract) |
| Order in catechol, pH 7.7 | ≈ 0.56 | James 1946 | opened-source (abstract) |
| Order in HQ for silver-ion (physical) reduction on Ag | 1 | James 1946 | opened-source (abstract) |
| Dependence on HQ species | divalent ion only (log-log linear) | Eskes 1982 | opened-source (abstract) |
| Order in silver mass during parallel grain development | **zero-order (linear in time)** in silver mass. The abstract argues that earlier "first-order" conclusions were distortions from plotting density, from polydispersity and from layer thickness. The rate constant depends on how fast pre-adsorbed bromide is replaced by developing species. It also reports a **constant growth rate of silver filaments** | Imaging Sci. J. (1998), DOI 10.1080/13682199.1998.11736441 (authors not shown in the snippet) | secondary-source |

### 1.4 Induction period

| Finding | Source | Confidence |
|---|---|---|
| Sugai models the induction period as an Arrhenius-type barrier: τ_ind = C·exp(U/kT), where U is set by Coulomb repulsion from the negative grain surface minus adsorption to the latent image. Adding KBr to a metol developer lengthens τ_ind (Fig. 2 has three KBr levels labelled about 4.5, 9.0 and 13.5; the units are illegible in the OCR), but U from the temperature dependence did **not** rise with KBr. He concludes that the charge effect matters but other factors are involved. | S. Sugai, "On the Induction Period of Photographic Chemical Development", J. Fac. Sci. Hokkaido Univ. Ser. II, 4(4), 268–274 (1953), hdl 2115/34210, Fig. 2 and eq. (1) | opened-source |
| With an experimental plate of 40 mol% AgI, the measured induction periods were: metol developer 70 ± 10 s at 4.5 °C, 30 ± 5 s at 12.2 °C and about 8 ± 3 s at 21.1 °C. Hydroquinone developer about 50 ± 10 s at 14.0 °C and 25–40 s at 21–23 °C; the OCR'd rows are out of order, so read Table 1 again from the page image before use. Against the charge-barrier expectation, HQ had the **smaller** temperature coefficient. On pure AgBr or ~2 mol% AgI plates with metol, τ_ind was too short to measure. | Sugai 1953, Table 1 | opened-source (the values are OCR'd; re-read from the page image) |
| The induction period lengthens with iodide content | Sugai & Furuichi 1953 (below), Fig. 6 | opened-source |
| James linked the induction period to the charge of the dissociated developer ion: it grows with the ionic charge (Sugai's summary of James, J. Phys. Chem. 43, 701 (1939) and 44, 42 (1940); J. Chem. Phys. 11, 338 (1943)) | Sugai 1953 §2 and references | opened-source (as Sugai's account); James originals not opened |
| Latent-image size: centres of ≥ 5 atoms develop fast; centres of < 5 atoms develop slowly with a long per-grain induction period (S+Au-sensitised emulsion) | Ohzeki 1996 | opened-source (abstract) |
| Phenidone lowers the charge barrier: one RIT thesis abstract attributes phenidone's activation of hydroquinone partly to "a charge-barrier depression effect". A quaternary salt (decyltrimethylammonium bromide) shifts density-time curves to shorter times and steepens them after the induction period | RIT theses (Eskes 1982; a thesis on phenidone derivatives), abstracts via search | secondary-source |

### 1.5 Bromide inhibition

| Finding | Source | Confidence |
|---|---|---|
| Br⁻ adsorption on AgBr in gelatin emulsions (potentiometric) follows a **Langmuir isotherm**, with about **10 % of the AgBr surface covered at saturation**. Br⁻ can also sit in the double layer as counter-ion to adsorbed cations such as cyanine dyes. | A.H. Herz & J.O. Helling, "Bromide adsorption in gelatin–silver halide emulsions", J. Colloid Sci. 17, 293–304 (1962), DOI 10.1016/0095-8522(62)90009-0 (authors and pages from Crossref) | **[unverified, spot-check 2026-11-21]** "Langmuir" and "~10 %" are still seen only in a search abstract; the original was not opened |
| Opened, related data from the same group: bromide adsorption at crystallographically defined AgBr surfaces in aqueous gelatin, pH 4.9, 40 °C. Cubic / octahedral: ΔG = −5.9 ± 0.1 / −6.2 ± 0.1 kcal mol⁻¹ Br⁻; ΔH = −10.7 ± 0.4 / −11.5 ± 0.4 kcal mol⁻¹; ΔS = −15 / −17 cal K⁻¹ mol⁻¹; **limiting area per adsorbed Br⁻ = 1.21 / 0.80 nm²**. Nearly twice as much Br⁻ sits on octahedral as on cubic surfaces, and "only a small fraction" of the geometrically available interface is involved in Br⁻ adsorption at lattice sites. Estimate (mine): a full Br⁻ layer would need ~0.17 nm² per Br⁻ on AgBr(100) and ~0.14 nm² on (111) (a = 0.5775 nm, recalled), so saturation coverage is ~14 % (cubic) and ~18 % (octahedral). That is the same order as the "~10 %" above | A.H. Herz, "Surface Reactions of Some Aqueous Silver Halide Dispersions", Croat. Chem. Acta 53(2), 313–318 (1980), p. 317 Table (citing Gardner, Wrathall & Herz, Photogr. Sci. Eng. 21, 325 (1977)), https://hrcak.srce.hr/file/286842 | opened-source (table read from the page image); the coverage % is an estimate |
| KBr adsorption by AgBr (no gelatin) was measured potentiometrically and plotted against log[Br⁻] (−4 to −3 shown). The page image shows the caption as "30.0 C", so "80.0 C" was an OCR error (checked in the 2026-11-21 spot-check). Values are in a figure only and were not digitised here. | B.H. Carroll & D. Hubbard, "Kinetics of reaction between silver bromide and photographic sensitizers", Bur. Stand. J. Res. 12, 329–344 (1934), RP659, §IV.1, Fig. 4 | opened-source (figure only) |
| Induction period ↑ and rate ↓ with bromide; bromide effect weaker at higher pH | Durbin 1976; Brainerd 1974 (RIT theses) | secondary-source |
| Rate in granular development is set by how fast pre-adsorbed bromide is replaced by the developing species | ISJ 1998 abstract | secondary-source |
| Classic "bromide shifts the curve along log H by a constant amount ∝ log[Br⁻]" and Nietz's "reduction potential" scale (the bromide concentration ratio that gives equal depression, relative to hydroquinone) | A.H. Nietz, *The Theory of Development*, Kodak monograph (Van Nostrand, 1922); scan on archive.org (not fetched; not whitelisted) | recalled. The values (metol ≈ 20, HQ = 1 on Nietz's scale) are not quoted here as data. Confirm in the 1922 monograph |

**Isotherm form for D1.** Competitive Langmuir adsorption on the reaction surface:

  θ_dev = K_D·a_D / (1 + K_D·a_D + K_Br·[Br⁻] + K_I·[I⁻])

Here a_D is the activity of the active ionised species. The rate is ∝ θ_dev, or ∝ θ_dev·(1 − θ_Br) in a Langmuir–Hinshelwood dual-site variant. This reproduces three observations: James's fractional order (0.6 at moderate coverage), the weaker bromide effect at higher pH (larger a_D outcompetes Br⁻), and Br⁻ lengthening the induction period. K_Br, K_D and K_I are **fitted**. The ~10–18 % saturation coverage (Herz 1980 limiting areas; the ~10 % figure is unverified) is a sanity bound, not a parameter value. **No published numeric K_Br for a developing grain was found.** The Herz ΔG (−5.9 to −6.2 kcal/mol at 40 °C, pH 4.9, undeveloped AgBr in gelatin) gives an order of magnitude for K_Br only once the standard state of the 1962/1977 papers is known. That standard state was not seen, so keep K_Br fitted.

### 1.6 Activation energies

| System | Ea | Source | Confidence |
|---|---|---|---|
| Commercial plate (Fuji A.I.), metol–sulfite–carbonate developer (metol 1 g, Na₂SO₃ 10 g, Na₂CO₃ 10 g in 300 mL), rate at D 0.1 → 0.4 | **0.64 ± 0.02 eV = 61.8 kJ/mol = 14.8 kcal/mol** | S. Sugai & J. Furuichi, "On the Activation Energy of Photographic Chemical Development", J. Fac. Sci. Hokkaido Univ. Ser. II, 4(4), 253–267 (1953), hdl 2115/34215, Table 1 and §4(I)i | opened-source |
| Same plate, hydroquinone developer (HQ 2 g, Na₂SO₃ 10 g, Na₂CO₃ 10 g in 300 mL), D 0.2 → 0.5, after the induction period | **0.70 ± 0.02 eV = 67.5 kJ/mol = 16.1 kcal/mol** | same, Table II | opened-source |
| Same plate, metol developer + KBr (metol **2 g**, Na₂SO₃ 10 g, Na₂CO₃ 10 g, KBr 4 g in 300 mL; note that this is twice the metol of the plain metol developer, p. 261 footnote) | **0.70 ± 0.02 eV** | same, Table 3 | opened-source |
| Experimental plates (Fuji), metol developer, AgI 0 / 2.5 / 7 / 40 mol% | **0.62 / 0.63 / 0.66 / 0.91 eV** (± 0.02) | same, Tables 4–5 | opened-source. **Confounded:** plates 1–3 used metol 1 g / sulfite 10 g / carbonate 10 g in 300 mL. Plate 4 (40 mol% AgI) used metol 4 g / sulfite 3 g / carbonate 30 g in **600 mL** and a 2-min exposure (100 W, 1 A) instead of 15–20 s. So the jump to 0.91 eV is not an AgI effect alone (spot-check 2026-11-21, Table 4 page image) |
| Uryu (cited by Sugai & Furuichi): commercial plate, metol developer, initial step | 0.66 eV | same, §2 | secondary-source (citation inside an opened paper) |
| My refit of Sugai & Furuichi's tabulated times (least squares on ln P vs 1/T) | metol 0.67 eV; HQ 0.73 eV; metol + KBr 0.80 eV | computed from Tables 1–3 | estimate. Their stated values are 0.64 / 0.70 / 0.70 eV. The spot-check (2026-11-21) re-read Tables 1–3 from the page images (metol 28/54/136/258 s at 26/18/10/2 °C; HQ 19/34/52/106 s at 25.5/19.8/14.3/8.8 °C; KBr 33/53/83/177 s at 24.5/19.0/14.3/9.5 °C) and reproduced 0.67 / 0.73 / 0.80 eV. **The KBr discrepancy is in the paper itself, not OCR.** Plates 1–4 refit: 0.62 / 0.64 / 0.70 / 0.96 eV |
| **Kodak apparent Ea, time to CI 0.56, training pairs only** (least squares over 18/20/21/22/24 °C; ln t = Ea/RT + c) | Tri-X × D-76 small tank **63.6 kJ/mol**, large tank 62.4. Tri-X × XTOL small tank **63.3**, large tank 63.5; J-109 Table 1 CI 0.56 / 0.72 / 0.82 → 63.5 / 62.0 / 65.7. Tri-X × HC-110 B small tank **71.4**, large tank 64.1. T-Max 100 × D-76 (F-4016 p. 7 factors 1.2 / 1.0 / 0.9 / 0.8 / 0.7) **66.2** | F-4017 (Oct 2021) p. 3 time table; J-109 (Feb 2018) Table 1; F-4016 (June 2018) p. 7 | estimate (from opened tables) |
| **Kodak apparent Ea, 1:1 dilutions (Tri-X)** | D-76 1:1 small / large tank **40.2 / 40.6 kJ/mol**. XTOL 1:1 small / large tank **38.9 / 40.4** | F-4017 p. 3 | estimate |
| Diffusion in water, from the viscosity of water (η = 1.053 cP at 18 °C and 0.911 cP at 24 °C, recalled CRC values), D ∝ T/η | **≈ 20 kJ/mol** (17.4 kJ/mol from η alone) | computed | estimate |

**Temperature-label correction (spot-check 2026-11-21).** The Ea values above use Kodak's rounded °C labels (18/20/21/22/24 °C). Kodak's nominal temperatures are 65/68/70/72/75 °F = 18.3/20.0/21.1/22.2/23.9 °C. On those, the same times give: Tri-X × D-76 small / large **67.4 / 65.9**; XTOL small / large **67.0 / 67.3**; HC-110 B small / large **75.6 / 68.3**; T-Max 100 × D-76 factors **70.2**; D-76 1:1 small / large **42.7 / 43.0**; XTOL 1:1 small / large **41.2 / 42.8** kJ/mol. So the full-strength class is **≈ 66–68 kJ/mol** (not ≈ 63) and the 1:1 class **≈ 41–43** (not ≈ 40). Use the °F-based values; the t₁₈/t₂₄ ratios (1.68 and 1.38–1.40) are unaffected.

Rounding check. The Kodak times are rounded to ¼ min. For D-76 (4¾–8 min) that is ±1.6–2.6 % per entry, about ±5 kJ/mol on Ea. For HC-110 B (2½–4½ min) it is about ±8 kJ/mol. So the full-strength versus 1:1 difference (≈ 67 against ≈ 42 kJ/mol on °F temperatures) is larger than the rounding error. The HC-110 B small-tank versus large-tank difference is not.

**Caveat: the Kodak rows are probably not independent measurements.** The Tri-X D-76 and XTOL rows are nearly the same table:
- Small tank: 8 / 6¾ / 6¼ / 5½ / 4¾ for D-76 against 8 / 7 / 6¼ / 5¾ / 4¾ for XTOL.
- t₁₈/t₂₄ = 1.68 for both, in both tank types.
- D-76 1:1 and XTOL 1:1 both give 1.38–1.40.

These are two chemically different systems (metol/HQ against Dimezone S/ascorbate) with identical temperature slopes. That points to a **shared Kodak time–temperature conversion** (one slope for full-strength developers and another for 1:1 dilutions), not to two independent measurements. No Kodak publication I have says how the tables were built. The extracted F-4017, J-78, J-109, J-24 and J-86 text contains no "time–temperature chart" statement. So treat the Kodak-derived values as **one observation per dilution class**, possibly a convention.

What this means:
1. Sugai & Furuichi's chemical Ea for metol and HQ (62–68 kJ/mol) agrees with the apparent Ea of Kodak's full-strength D-76 rows (66–67 kJ/mol on the °F temperatures; 62–64 on the rounded °C labels). That is consistent with full-strength D-76 development to CI 0.56 being largely reaction-controlled. It says **nothing independent about XTOL's agents**, whose rows look like copies.
2. The 1:1 class shows **≈ 41–43 kJ/mol** (°F basis), between the reaction value and the diffusion value (≈ 20). A partly transport-limited (boundary-layer or exhaustion-limited) regime would give this, and so might a Kodak convention.

For D1:
- Use Sugai & Furuichi as the per-agent prior for metol and HQ.
- Do **not** fit the phenidone, Dimezone S or ascorbate Ea to Kodak rows that copy the D-76 slope.
- Keep the 1:1 slope as a check of transport, labelled "possibly a convention".

Not found and held out:
- No per-agent Ea for **phenidone, Dimezone S, ascorbate or p-aminophenol** was found.
- **Held-out data not used here:** the T-Max 100 × XTOL rows of J-109 Table 1 and the HC-110 B factors on F-4016 p. 7 belong to held-out pairs. They were deliberately not fitted.

### 1.7 Growth law: autocatalytic area law versus linear filament growth

| Model | Evidence | Confidence |
|---|---|---|
| Rate ∝ developed-silver surface area. For a compact sphere this is ∝ m^(2/3); it gives sigmoidal single-grain kinetics. | The classic textbook assumption (James *TPP*; Mason) | recalled |
| Two stages: a slow, zero-order first stage (electron transfer to the grain, ∝ grain surface area × developer concentration), then an autocatalytic stage ∝ the area of the silver speck | Narrative in a patent's background section (USPTO, "Method and apparatus for producing a high speed high resolution radiation sensitive article…", seen via search) | secondary-source; weak authority |
| Silver mass grows **linearly in time** in parallel grain development; silver filaments grow at constant rate; apparent first-order kinetics is an artefact of plotting density | ISJ 1998 abstract (DOI above) | secondary-source |

For D1: per grain-size class c, write dm_c/dt = k(T)·θ_dev·g(m_c), where g is either (a) a constant (filament growth, linear in time) or (b) S_c,0 + κ·m_c^(2/3) (area law). Treat the choice as a model-structure test against the D-76 time series (characteristic curves at 6/8/10/12 min for Tri-X; 6/7.5/10 min for T-Max 100). Do not fix it a priori. The literature does not settle it.

### 1.8 Fog

**Not found:** no quantitative fog development rate (rate for grains without a latent image, versus developer, bromide and time) in any source I could open. Constraints that do exist:
- In the Tri-X × D-76 large-tank 35 mm curves (6–12 min), the lowest plotted density (base + fog) lies in **0.27–0.32**, and the 120 curves in 0.16–0.21. So base + fog rises by **≤ 0.05** over a 2× change in development time (datasheets.md §2.2; estimate from vector paths).
- For T-Max 100 × D-76 small tank: 0.21–0.22 (datasheets.md §3.2).

Model fog as the same rate law applied to a "no-latent-image" state, with a much longer induction period (a fitted barrier ΔU_fog) plus a small spontaneous term. Calibrate ΔU_fog on these spreads.

### 1.9 Dependence on grain size

**Not found** quantitatively. Two things depend on size:
- the time to complete a grain, which scales with grain volume divided by the rate of the developing surface (recalled, James *TPP*);
- the induction period, which depends on the latent-image size distribution of the class. latent-image.md §1.3(d) already makes this size dependent.

The Carroll & Hubbard (1930) opened text supports one indirect, qualitative point: solvent ("fine-grain") developers trade sensitivity for resolution. See §2.3.

### 1.10 Proposed rate-law structure for D1 (synthesis; structure from the sources above, constants fitted)

Per depth cell, per size class c and latent-image state j (j = 0 for fog). The developing agent D_i is in its active ionised form:

  r_c,j = N_c,j · Σ_i k_i·exp(−Ea_i/R·(1/T − 1/T_ref)) · θ_i · g(m_c) · H(t − τ_ind,c,j)

  θ_i = K_i·a_i / (1 + Σ_k K_k·a_k + K_Br·[Br⁻] + K_I·[I⁻])

  a_i = ionised fraction(pH; pKa_i) × [D_i]_total

  τ_ind,c,j = τ₀·exp((U_charge(z_i, I_s, [Br⁻]) − U_LI(j))/kT)

| Symbol | Meaning |
|---|---|
| N_c,j | number of grains of size class c in latent-image state j |
| g(m_c) | growth law (§1.7) |
| H(t − τ_ind) | step at the end of the induction period; replace by a smooth hazard for population averaging |
| z_i | charge of the active species i |
| I_s | ionic strength, which enters through the gelatin charge barrier (James 1944) |
| U_LI(j) | adsorption term that grows with latent-image size (Ohzeki 1996) |

Population averaging: the induction times within a class are distributed, so develop the class fraction with a first-order hazard λ_c,j(t) = λ₀·exp(−ΔU/kT) instead of a step.

Status of the constants:
- **Literature**: pKa (§1.2).
- **Literature with a large spread**: Ea for metol and HQ (§1.6).
- **Fitted**: k_i, K_i, K_Br, K_I, τ₀, U terms, ΔU_fog, and Ea for phenidone, Dimezone S, ascorbate and p-aminophenol.

---

## 2. Superadditivity and sulfite

### 2.1 Regeneration mechanism

| Finding | Source | Confidence |
|---|---|---|
| PQ (phenidone–hydroquinone) mechanism: a coloured, reversible, mesomeric free radical forms by loss of one electron from the pyrazolidone **ion**. Hydroquinone regenerates phenidone from this radical. | A.J. Axford & J.D. Kendall, "Phenidone—Part 1; The Mechanism of its Developing Action and its Super-additivity with Hydroquinone", J. Photogr. Sci. 2(1), 1–7 (1954), DOI 10.1080/03700240.1954.11736538 | secondary-source (abstract via search; T&F 403) |
| MQ: hydroquinone regenerates Elon (metol) by reducing oxidised Elon | *British Kinematography* (1948) article, seen via lantern.mediahist.org search snippet | secondary-source (pointer only) |
| Levenson worked out the regeneration mechanism. A later Levenson paper is cited as J. Photogr. Sci. 24, 68 (1976) (title not seen). | reference lists seen via search (informahealthcare.com) | secondary-source (citation only) |
| Willems & Van Veelen are cited as Photogr. Sci. Eng. 6, 39 (1962) and 6, 49 (1962). **That these papers are on superadditivity is not confirmed** from what I saw. | reference lists via search | secondary-source (citation only) |
| Ascorbate superadditive with 3-pyrazolidones: US 2,688,549. Activity of iso-ascorbic acid first reported by Maurer et al., Ber. 67, 1239 (1934). Ascorbate derivatives are weak developers alone in the pH range where HQ works. | Agfa EP 0 573 700 A1, p. 2 | opened-source (patent background) |
| In a regenerative pair, the electron-transfer agent (phenidone, Dimezone S, metol, p-aminophenol) adsorbs strongly to AgX and has a lower reduction potential. The regenerating agent (HQ, ascorbate) has a higher reduction potential and adsorbs poorly. The stability of the ETA's one-electron oxidised form controls how strong the superadditivity is. | Wikipedia "Superadditive developer" (pointer); consistent with Axford & Kendall's abstract and with the RIT phenidone-derivative thesis abstract ("regeneration … partly depending on the stability of the phenidone semiquinone") | secondary-source; Wikipedia is not a citation |
| A forum argument that regeneration of metol by HQ does not operate without sulfite or ascorbate present | Photrio thread "Superadditivity in developers — and beyond" | forum, pointer only. Not used |

**Rate constants and relative activity numbers: not found.** I found no published second-order rate constant for:
- ETA⁺˙ (or ETA_ox) + HQ²⁻ / HQ⁻ → ETA + HQ⁻˙ (or Q);
- phenidone-radical disproportionation;
- metol quinonemonoimine reduction by HQ;
- any reaction involving ascorbate.

Pulse-radiolysis compilations of phenoxyl and semiquinone radical rate constants exist (P. Neta & J. Grodkowski, J. Phys. Chem. Ref. Data, 2005; NIST, seen via search). They hold hydroquinone/semiquinone data but **no phenidone entries** were seen. The "phenidone is 10–18× more active than metol" figure seen on Wikipedia is **not** a citable number.

Proposed reaction set (structure: Levenson; Axford & Kendall). P = phenidone-type ETA (or metol), H = HQ or ascorbate:

| # | Reaction | Rate constant | Status |
|---|---|---|---|
| R1 | P⁻ + Ag⁺(grain or speck) → P˙ + Ag⁰ | k_P (surface, Langmuir θ_P as §1.10) | fitted |
| R2 | H²⁻ (or HA⁻) + Ag⁺ → H˙⁻ + Ag⁰ | k_H (surface; long induction; charge barrier) | fitted |
| R3 | P˙ + H²⁻ → P⁻ + H˙⁻ (regeneration; in solution or adsorbed) | k_reg | fitted (no measured value found) |
| R4 | 2 H˙⁻ ⇌ H²⁻ + Q (semiquinone disproportionation) | fast | QSSA (§5) |
| R5 | 2 P˙ → P⁻ + P_ox (irreversible loss of the ETA, e.g. to hydroxypyrazole) | k_loss | fitted; controls ETA exhaustion |
| R6 | Q + SO₃²⁻ (+H₂O) → HQSO₃⁻ + OH⁻ | k_sulf | **literature value exists but was not obtained**: LuValle, J. Am. Chem. Soc. (1952), cited in later kinetics papers (secondary). Fit or take from LuValle once opened |

The superadditivity factor is then an **output**: the rate with P + H divided by (rate with P alone + rate with H alone). Calibrate k_reg and k_loss on the pairs where the regenerator differs: D-76 (metol + HQ) against XTOL (Dimezone S + ascorbate). HC-110 is msds-estimated (recalled to contain HQ and a pyrazolidinone ETA), which makes it the weakest constraint.

### 2.2 Sulfite

| Role | Finding | Source | Confidence |
|---|---|---|---|
| AgBr solvent | Solubility of AgBr at 30.0 °C: **7.5 ± 0.5 × 10⁻⁵ N** in 0.010 N Na₂SO₃ at pH 8.5. Reduced to **3.5 × 10⁻⁵ N** in a gelatin emulsion at pH 6.2 at the same sulfite level, and to **4 × 10⁻⁵ N** by adding 4 × 10⁻⁴ N KBr at pH 8.5. That is about a 2× reduction, not 20×. Also **2.9 × 10⁻⁵ N** in 0.0030 N and 2.0 × 10⁻⁵ N in 0.0012 N. About double the values of Mees & Piper; the authors add that at the reaction temperatures the solubility was "probably several times as great". **Corrected in the 2026-11-21 spot-check:** the earlier text had 3.5 × 10⁻⁶ and 4 × 10⁻⁶, OCR errors in the exponent. Re-read from the 500-dpi page image | Carroll & Hubbard 1934, RP659, §IV.2, p. 340 | opened-source |
| AgBr solvent → metallic silver | AgBr + sulfite gives metallic silver autocatalytically (relevant to fog and physical development) | Carroll & Hubbard 1934, abstract | opened-source |
| Fine-grain / solvent developers | "Special developers recommended for reducing the graininess … are all based on increasing the solvent action of the solution on silver bromide, and decreasing the reduction potential." Increased resolution was "in general, accompanied by decreased sensitivity". | B.H. Carroll & D. Hubbard, "A comparison of resolving power and sensitivity of photographic plates with varying development", Bur. Stand. J. Res. 5, 1–12 (1930), RP183, abstract | opened-source |
| Scavenging oxidised developer | Hydroquinone sulfonate is a main development by-product (it fluoresces in UV; Kodak used this to visualise the depleted layer at the film surface) | W.C. Snyder, "An Investigation of Agitation in a Continuous Immersion Film Process", J. SMPTE 75(10) (1966), DOI 10.5594/J07079, abstract | opened-source (abstract) |
| Sulfonation chemistry | Between pH 4.0 and 7.5 the benzoquinone + S(IV) reaction yields only hydroquinone monosulfonate. Below pH 4 it also gives HQ + sulfate (LuValle 1952, as summarised in a 2023 kinetics paper) | J.E. LuValle, J. Am. Chem. Soc. (1952) | secondary-source |
| Quinone monosulfonate + sulfite gives a strong reductant in alkali, more negative than HQ, which may add to development | J. Soc. Photogr. Sci. Technol. Japan 42(5), 325–328 (1979) (abstract via search) | secondary-source |
| Inhibiting autoxidation of HQ | Sulfite inhibits HQ autoxidation (James & Weissberger 1939). The initial autoxidation rate of HQ was reported inversely proportional to [sulfite] | citation via search | secondary-source |
| Colour-developer comparison | For p-phenylenediamines, sulfite reacts with the quinonediimine to give the developer sulfonate. Autoxidation consumes about one more mole of sulfite to sulfate (via H₂O₂), giving a sulfite : developer consumption ratio of 2 that EDTA + Fe³⁺ traces lower towards 1 | J. Eggers, "Chemie der photographischen Farbentwicklung", Chimia 15, 499–510 (1961), §3245 | opened-source |
| Buffer | In D-76 the sulfite/bisulfite pair is the **main buffer** near pH 8.4 (estimate below). Borax adds about 5 % of the buffer capacity | computed | estimate |

Buffer estimate for D-76:
- Formula (recalled, the public Kodak formula; confirm in Kodak J-1 or Anchell & Troop): metol 2 g/L, Na₂SO₃ 100 g/L (0.793 M), HQ 5 g/L (45.4 mM), borax decahydrate 2 g/L (5.24 mM, i.e. 21 mM B). pH about 8.2–8.5 (recalled).
- Buffer capacity β = 2.303·C·α(1 − α), with pKa(HSO₃⁻) ≈ 7.2 and pKa(B(OH)₃) ≈ 9.24 (both recalled; confirm in CRC):
  - β_sulfite ≈ 0.10 mol L⁻¹ pH⁻¹;
  - β_borate ≈ 0.005 mol L⁻¹ pH⁻¹.
- Consequence: the "buffer species" that must diffuse with the developer in D-76 is mainly **SO₃²⁻/HSO₃⁻**, not borate.
- **[unverified inputs, spot-check 2026-11-21]** The arithmetic reproduces (β_sulfite 0.101, β_borate 0.0053 mol L⁻¹ pH⁻¹ at pH 8.4). But the formula and both pKa values are recalled. At D-76's ionic strength (about 2.4 mol/L from 0.79 M Na₂SO₃) the effective pKa′ of HSO₃⁻ is likely lower than 7.2, which would cut β_sulfite by several times. For example, with pKa′ = 6.5 it would be ≈ 0.02 M/pH, against ≈ 0.009 for borate at pKa′ ≈ 8.9. Sulfite probably still dominates, but by about 2–3×, not 20×. Confirm with activity corrections before using it.

Effect of sulfite on speed and grain (qualitative): solvent action dissolves fine-grain edges and the smallest latent images, lowering speed and granularity and adding a solution-physical-development component (Carroll & Hubbard 1930; recalled from Haist and Mason). Kodak's own statements:
- "Dilution at 1:1 will provide slightly greater film speed, enhanced sharpness and shadow detail, and slightly more grain" (J-109, Feb 2018, p. 2; opened).
- "Dilute developers … give slightly higher film speed and a slight increase in graininess" (F-4016, June 2018, p. 3 footnote; opened).

Both are qualitative validation targets: lower sulfite at 1:1 means less solvent action, and more depletion at 1:1 means more adjacency.

**Solution-physical development model for D1.** The source is dissolved Ag(SO₃)ₙ complexes, and silver deposits on developing centres. The rate is first-order in dissolved Ag and proportional to the catalytic area (James 1946: proportional to [HQ] for Ag⁺ reduction on silver). The dissolution rate is fitted. The complexation constants of Ag⁺ with SO₃²⁻ are **literature values I did not obtain** (Smith & Martell critical stability constants would have them; recalled).

---

## 3. Transport

### 3.1 Diffusion coefficients in water

| Species | D (25 °C, cm² s⁻¹) | Source | Confidence |
|---|---|---|---|
| Br⁻ | 2.08 × 10⁻⁵ | CRC Handbook, "Ionic conductivity and diffusion at infinite dilution" | recalled |
| OH⁻ | 5.27 × 10⁻⁵ | same | recalled |
| Na⁺ / K⁺ | 1.33 × 10⁻⁵ / 1.96 × 10⁻⁵ | same | recalled |
| SO₄²⁻ / CO₃²⁻ / HCO₃⁻ | 1.07 / 0.92 / 1.19 × 10⁻⁵ | same | recalled |
| SO₃²⁻, B(OH)₄⁻ | not recalled reliably | — | not found. Use the salt values in gel from Iwano (§3.2) as the anchor |
| Hydroquinone | 1.0 × 10⁻⁵ | Wilke–Chang, Le Bas V = 110.8 cm³ mol⁻¹ | estimate (±15 %) |
| p-Aminophenol | 9.7 × 10⁻⁶ | Wilke–Chang, V = 117.6 | estimate |
| Metol (free base) | 8.7 × 10⁻⁶ | Wilke–Chang, V = 141.3 | estimate |
| Ascorbic acid / ascorbate | 8.4 × 10⁻⁶ | Wilke–Chang, V = 151.3 | estimate |
| Phenidone | 7.7 × 10⁻⁶ | Wilke–Chang, V = 173.9 | estimate. Dimezone S is larger, so somewhat lower |
| N,N-diethyl-p-phenylenediamine | 6.8 × 10⁻⁶ | Wilke–Chang, V = 213.5 | estimate (for comparison with Iwano) |

Wilke–Chang: D = 7.4 × 10⁻⁸·(φM)^½·T/(η·V^0.6), with φ = 2.6, M = 18.0, η = 0.890 cP. It ignores ionic charge. For ions, the Nernst–Planck coupling with counter-ions (electroneutrality) matters more than the ±15 % error.

### 3.2 Diffusion in swollen gelatin (primary data found)

H. Iwano (Fuji Photo Film, Ashigara), J. Soc. Photogr. Sci. Technol. Japan. Opened-source; tables read from the rendered pages:
- "The Rate of Penetration of Photographic Chemicals into Gelatin Layers (2): The Diffusion-retarding Factors in Gelatin", 33(3), 143–149 (1970). Cited below as I-1970.
- "Penetration of Photographic Chemicals into Gelatin Layers (IV): Comparison of Penetration Rates of Reacting Species", 35(2), 98–104 (1972), DOI 10.11454/photogrst1964.35.98. Cited below as I-1972.

**Developing agents** (D in cm² s⁻¹):

| Species | Medium | Conditions | D | Ratio to water (estimate, §3.1) | Source |
|---|---|---|---|---|---|
| Hydroquinone | 8 % gelatin gel | 23 °C, pH 9.0, 2 × 10⁻² M | 0.99 × 10⁻⁶ | ≈ 0.10 | I-1970 Table 5; I-1972 Table 3 (25 °C: 0.99 × 10⁻⁶). Spot-check 2026-11-21: I-1972 Table 3 repeats the I-1970 Table 5 values (29.0 / 11.7 / 9.9 × 10⁻⁷) exactly but labels them 25 °C, against 23 °C in 1970. Count it as **one measurement reprinted**, not two that agree, with the temperature uncertain by 2 °C |
| Metol | 8 % gel | 23 °C, pH 9.0 | 1.17 × 10⁻⁶ | ≈ 0.13 | same |
| N,N-diethyl-p-phenylenediamine | 8 % gel | 23 °C, pH 11.0 | 2.90 × 10⁻⁶ | ≈ 0.43 | same |
| Hydroquinone | coated layer, swelling α = 3.03 | 25 °C, pH 10.5, indicator-layer method | 8.0 × 10⁻⁷ (swollen basis); 0.88 × 10⁻⁷ (dry basis) | ≈ 0.08 | I-1970 Table 7 |
| Metol | same layer | pH 9.0 | 11.7 × 10⁻⁷ (swollen); 1.28 × 10⁻⁷ (dry) | ≈ 0.13 | I-1970 Table 7 |
| DEPPD | same layer | pH 11.0 | 16.9 × 10⁻⁷ (swollen); 1.85 × 10⁻⁷ (dry) | ≈ 0.25 | I-1970 Table 7 |
| Metol, true D | coated layer, Na₂SO₄ 100 / 200 g/L (gel 28.6 / 65.6 wt %) | — | 2.1 × 10⁻⁷ / 3.9 × 10⁻⁷ (apparent 0.24 × 10⁻⁷ at 65.6 %) | — | I-1970 Table 4. This is lower than Table 7 by about 5×; the tables are not mutually consistent |
| "Developing agents", summary statement | gelatin layer | — | about **1 × 10⁻⁷ (dry-film basis) ≈ 1 × 10⁻⁶ (normally hardened swollen film)**. Penetration is "greatly reduced" during diffusion-controlled development, depending on AgX and developer concentration | — | I-1972 §3-1 |

**Inorganic species** (D in cm² s⁻¹; 8 % gel, 25 °C, 0.1 M, conductivity method, s.d. 6 %; I-1972 Table 4 and §5):

| Species | D | Note |
|---|---|---|
| Na₂SO₃ | 2.01 × 10⁻⁵ | |
| Na₂CO₃ | 2.15 × 10⁻⁵ | coated layer: 1.02 × 10⁻⁵ at 8 % (I-1970 §5.1) |
| Na₂SO₄ | 0.98 × 10⁻⁵ | |
| KI | 1.49 × 10⁻⁵ | |
| KBr | 1.59 × 10⁻⁵ | Stonham & Kragh: 1.54 × 10⁻⁵ (8 %). Iwano also gives 3.03 × 10⁻⁵ in 4 % gel and 1.10 × 10⁻⁵ in 20 % gel |
| KCl / NaCl | 2.56 / 2.55 × 10⁻⁵ | above the water values, so treat ±30 % as the real scatter |
| Na₂S₂O₃ | 1.61 × 10⁻⁵ | literature values for thiosulfate in the same gel concentration range from 0.10 to 1.8 × 10⁻⁵ (I-1972 Table 5). This spread shows how method-dependent these numbers are |

**OH⁻ and HCHO** (coated layer, true D; I-1970 Table 4 and Fig. 10):

| Gelatin in swollen layer (wt %) | 11.0 | 17.1 | 28.6 | 43.2 | 65.6 |
|---|---|---|---|---|---|
| Swelling α | 12.1 | 7.84 | 4.68 | 3.10 | 2.04 |
| OH⁻, D (× 10⁻⁷ cm² s⁻¹) | 94 | 79 | 48 | 26 | 19 |
| OH⁻, apparent D (× 10⁻⁷ cm² s⁻¹) | 6.0 | 3.15 | 1.21 | 0.38 | 0.18 |
| HCHO, D (× 10⁻⁷ cm² s⁻¹) | 37.8 | 27.4 | 19.0 | 15.9 | 11.7 |

- OH⁻ extrapolates to the water value at infinite swelling.
- The apparent D for OH⁻ is about 16–106× smaller (94/6.0 at 11 % gelatin up to 19/0.18 at 65.6 %; corrected in the spot-check from "5–100×") because OH⁻ reacts with gelatin functional groups (reactive group content 4.8 × 10⁻⁴ mol/g for OH⁻). This matters for **pH transport**: the gelatin itself is a buffer.
- HCHO in gel: 1.5 × 10⁻⁵ (3.9 %), 7.4 × 10⁻⁶ (7.6 %), 2.3 × 10⁻⁶ (15 %) (I-1970 Table 6).

The inorganic values are **salt (ambipolar) diffusion coefficients** from conductivity. Compare them with salt D in water, not with the single-ion values of §3.1. Salt D in water at 25 °C (recalled, CRC): KBr ≈ 2.0 × 10⁻⁵, NaCl ≈ 1.6 × 10⁻⁵, Na₂SO₄ ≈ 1.2 × 10⁻⁵ cm²/s. The resulting gel/water ratios, in 8 % gel:
- KBr ≈ 0.8;
- Na₂SO₄ ≈ 0.8;
- NaCl ≈ 1.6, which is impossible and shows how large the scatter is.

Findings that matter for D1:
1. **Developing agents diffuse about 10× slower in gel than in water.** Small inorganic salts slow down only by a factor of about 1.2–1.3 (ratio about 0.8). A single "gelatin tortuosity factor" applied to all species would be wrong by an order of magnitude for the developer. Use per-species D_gel/D_water: about 0.1 for HQ and metol and about 0.3–0.4 for a phenylenediamine (opened data). For phenidone, Dimezone S and ascorbate the ratio is **not found**; take 0.1–0.3 as an assumed prior.
2. D in gel depends steeply on the swollen gelatin concentration. Swelling (§3.3) is therefore a needed film parameter.
3. "Penetration of HQ into an emulsion under diffusion-controlled development" behaves as a reaction front (I-1972 Fig. 3, theoretical curves for HQ at 1.1–35.2 g/L and pH 9.4–12.8). This is the moving-boundary unit test in §5.7.
4. Gelatin buffers OH⁻. The design's "pH as a transported buffered field" should include gelatin carboxyl and amino groups as an immobile buffer in the emulsion cells. Iwano's 4.8 × 10⁻⁴ mol/g is the opened number for the OH⁻-reactive group content.

**Non-Fickian penetration.** A. Green & G.I.P. Levenson (Kodak Ltd), "The Penetration of Solutes into Swollen Gel Layers", J. Photogr. Sci. 24, 193–196 (1976) (secondary-source; catalogue abstract via search):
- They measured penetration of 1 M thiosulfate and D-163 developer (full strength and diluted) through gelatin of **4.6, 10.5 and 20 µm dry**, with an AgCl indicator underlayer.
- Penetration time ∝ (swollen thickness) × (1/water fraction)³.
- Penetration was **non-Fickian**: the front advanced at roughly constant velocity in layers up to 20 µm.

This conflicts with Iwano's √t law (I-1970 Figs. 6–8: depth ∝ t^½ for HCHO, OH⁻ and metol). The likely reason is that swelling, which proceeds at the same time, controls penetration in Green & Levenson's set-up. For D1: start from Fickian transport with D_gel from Iwano. Add a swelling-front time (T½ of order 10–30 s, see §3.3) as an initial-penetration delay. Do not fit a non-Fickian law with no data behind it.

### 3.3 Swelling and swollen thickness

| Item | Value | Source | Confidence |
|---|---|---|---|
| Swelling ratio of emulsion-side layers in processing | (swollen − dry)/dry = **150–400 %**, i.e. swollen/dry = 2.5–5. Backing layers 150–500 % | colour-film patents (USPTO full-text search results) | secondary-source (colour materials, not B&W camera film) |
| Swollen/dry thickness | "normally adjusted so that maximum swollen thickness is **2–6×** dry" | patent seen via search (Japanese-language summary) | secondary-source |
| Swelling half-time | T½ ≈ 20–30 s or less in colour developer at 30 °C (saturated = 90 % of maximum swollen thickness at 3 min 15 s) | colour-film patents via search | secondary-source |
| "Normally hardened" swollen film | Iwano converts dry-basis to swollen-basis D with α = 3.03 (factor α² ≈ 9.2) in I-1970 Table 7. Separately, I-1972 §3-1 says ~1 × 10⁻⁷ (dry basis) corresponds to ~1 × 10⁻⁶ "in a normally hardened swollen film", i.e. a factor of about 10. Combining the two gives α ≈ 3 for a normally hardened film. That combination is mine; Iwano does not state α for "normally hardened" | I-1970 Table 7; I-1972 §3-1 | opened-source (each statement); estimate (the link) |
| Total hydrophilic colloid in a colour film | ≤ 28 µm dry, preferably ≤ 16 µm | colour-film patent via search | secondary-source |
| Emulsion-layer dry thickness for Tri-X / T-Max 100 | **not found** (see latent-image.md §2.3: coverage-based estimate 4.1 µm for 1.75 g/m² Ag + 4.84 g/m² gelatin) | — | — |
| Overcoat as a quiescent boundary layer | Kodak: a hydrophilic overcoat of about ≥ 3 g/m² "establishes a quiescent boundary layer". Effective D threshold for a "diffusible" compound in swollen gelatin: > 5 × 10⁻⁷ cm² s⁻¹ | US 5,837,437 A (Eastman Kodak; Sharma & Texter) | opened-source |

D1 prior (assumed): dry emulsion 5–10 µm plus a ~1 µm overcoat (assumed); swelling α = 3 (range 2–5) in alkaline developer; swollen emulsion **15–30 µm**. Fit α within this range only if the CI(t) data demand it. Otherwise fix it at 3 and tag it `assumed`, citing Iwano for "normally hardened".

### 3.4 Boundary layer and agitation

| Finding | Source | Confidence |
|---|---|---|
| Mass-transfer coefficient h = D/δ in a forced-flow (slot-jet) processor: max 0.0086 cm/s (u = 100 cm/s), nominal **0.0030 cm/s** (u = 25.4 cm/s), min 0.0018 cm/s. Roller-transport prior art: 0.0022 cm/s, and a "more practical" roller case 0.00048 cm/s. Claimed range 0.001–0.010 cm/s. Calculations use D = 5.0 × 10⁻⁶ cm²/s. | US 5,315,338 A, "Apparatus for enhancing heat and mass transfer in a fluid medium", detailed description and claim 3 | opened-source |
| → equivalent δ = D/h | for D = 5 × 10⁻⁶: δ = **6–28 µm** (slot jets), 23 µm (roller, best), **104 µm** (roller, practical) | computed | estimate |
| Good uniformity needs turbulence strong enough to thin and even out the stagnant layer. With intermittent jets, agitation was raised mainly by applying jets more often. | C.E. Ives & E.W. Jensen, "The Effect of Developer Agitation on Density Uniformity and Rate of Development", J. SMPE 40(2) (1943), DOI 10.5594/J14279, abstract | opened-source (abstract) |
| The boundary layer between moving film and solution is laminar at all practical processing speeds | Snyder 1966, abstract | opened-source (abstract) |
| A 1950s SMPTE analysis: stagnant-layer thickness inversely proportional to a Reynolds-type parameter | SMPTE (search snippet; lantern.mediahist) | secondary-source |
| Agitation in Kodak small tanks: 5–7 inversions in 5 s, then 5 s every 30 s. Large tanks: lift-drain-tilt once per minute | F-4017 (2021) p. 3; F-4016 (2018) p. 2–3; J-78 (2017) p. 1–2 | opened-source |
| Boundary-layer thickness in a **hand-inverted small tank** | **not found** | — | — |
| Rotary-tube (continuous) agitation | No measured δ found. Rotary processing is continuous agitation, so expect the forced-flow regime (δ about 20–100 µm, from the US 5,315,338 values above) without the burst-and-regrowth cycle of hand tanks. Split note: F-4017's 400TX rotary-tube rows are training data. The T-Max 100 rotary-tube XTOL CI curves (F-4016 p. 9) belong to a held-out pair | estimate | estimate |

Estimate for small tanks. After an agitation burst the liquid near the film is quiescent. A depletion layer then grows by diffusion as δ(t) ≈ (πDt)^½:
- D = 1 × 10⁻⁵ cm²/s: δ(25 s) ≈ 280 µm, δ(55 s) ≈ 420 µm.
- For Br⁻ (D = 2 × 10⁻⁵): 400–590 µm.

Natural convection and reel geometry cap the growth (recalled: typical natural-convection diffusion layers in unstirred aqueous cells are a few hundred µm). D1 prior: in a well-mixed burst δ resets to 20–100 µm (from the processor values above), then grows as (πDt)^½ to a cap δ_max = 200–500 µm (fitted, single parameter), shared by all films and developers. Kodak's two agitation regimes (30 s and 60 s) then differ only through the interval, which is the observable. This also handles the agitation mismatch between the Tri-X D-76 curves (large tank) and the T-Max 100 D-76 curves (small tank) that datasheets.md §7 raised.

### 3.5 Adjacency effects

| Finding | Source | Confidence |
|---|---|---|
| Nelson introduced the **chemical spread function** (ChSF) as the analogue of the optical line-spread function. Citation confirmed: **C.N. Nelson, Photogr. Sci. Eng. 15, 82 (1971)** | Rajkowski & Nowak 1998, ref. [2] | opened-source (as citation) |
| ChSF modelled with Frieser's exponential LSF form, L(x) ∝ 10^(−\|x\|/K), where K is the distance at which L falls to 1/10 of its maximum. In this model a restraining (bromide) ChSF is negative and an accelerating ChSF is positive. Rajkowski & Nowak *assumed* chemical ranges K_c = 25, 50 and 100 µm with 30–60 % ChSF weight. These are model choices, not measurements. | B. Rajkowski & P. Nowak, "Model investigations of the influence of adjacency edge effects on the shape and properties of the limiting curve in the silver halide light-sensitive materials", Optica Applicata 28(2), 105–113 (1998); Frieser, *Photographische Informationsaufzeichnung* (Focal Press, 1975), ch. 3 | opened-source |
| Nelson-type small-area density prediction for a fine-grain aerial duplicating film in D-76 (1:4) agreed to within 0.02 D on bar peaks, with slight overshoot at bar edges | Cline, RIT thesis (1982), "Prediction of small-scale densities using the chemical and optical spread functions" | secondary-source (abstract) |
| Chemical MTF from MTF_total = MTF_chem × MTF_optical (Kodak High Contrast Copy 5069 in D-76) | Laitre, RIT thesis (1979) | secondary-source |
| Maximum diffuse density correction for edge effects: 0.14 | Univ. of Arizona dissertation, "The influence of adjacency effects on the restoration of noisy photographic images" | secondary-source |
| Released bromide diffuses "a short distance across borders from highlight areas" | patent background via search | secondary-source |
| Kriss, Nelson & Eisen diffusion model of adjacency effects; Simonds; Eberhard and Kostinsky original length scales | **not found**. No accessible abstract or citation was confirmed. The Wikipedia figure "lines within about 1 mm" is not citable | — |

Length-scale estimate (for checking, not for fitting):
- Inside the swollen layer, developer diffusion over a development time of t = 400 s gives √(D_gel·t) = √(1 × 10⁻⁶ × 400) ≈ 200 µm laterally. For Br⁻ (D_gel ≈ 1.5 × 10⁻⁵) it is about 780 µm.
- But most lateral exchange happens in the liquid boundary layer, so the effective ChSF width is set by δ (100–500 µm, §3.4) and by the agitation interval.
- Rajkowski & Nowak's assumed K_c = 25–100 µm sits at the low end of this range.

The opened MTF figure constrains the effective chemical length. KODAK PROFESSIONAL TRI-X 320 and 400 Films, F-4017 (Oct 2021), p. 8, "Modulation Transfer Function" (D-76, large tank, diffuse visual), read from the rendered page (eyeballed, estimate):
- the response sits slightly **above 100 %, about 105–110 %, from about 2 to 20 cycles/mm**, peaking near 10 cycles/mm;
- it crosses 100 % near 25–30 cycles/mm;
- it falls to about 25 % near 70 cycles/mm.

An overshoot over 5–20 cycles/mm corresponds to half-periods of 25–100 µm, the scale of the chemical spread. Per F-4017 (May 2007) p. 10 (datasheets.md §2.1), this MTF was measured on the older TRI-X Pan 5063 emulsion. Fit δ_max only after the optical PSF is fixed, as the design already says.

### 3.6 Transport budget: why the boundary layer controls compensation

All estimates. Assumptions: coated silver 3 g/m² (assumed, since latent-image.md §2.3 found no B&W value) and a swollen emulsion of 15 µm.
- Silver concentration in the swollen layer: 3/107.87 mol m⁻² ÷ 15 µm ≈ **1.9 M**.
- Developing half of it over 300 s releases J ≈ 5 × 10⁻⁵ mol m⁻² s⁻¹ of Br⁻. It also demands about 2.5 × 10⁻⁵ mol m⁻² s⁻¹ of HQ (0.5 HQ per Ag) and releases 2.5 × 10⁻⁵ mol m⁻² s⁻¹ of net H⁺ (§4).
- The steady concentration drop across a boundary layer δ is ΔC = J·δ/D:

| δ | Br⁻ excess at the film (D = 2 × 10⁻⁹ m²/s) | HQ depletion (D = 1 × 10⁻⁹ m²/s) |
|---|---|---|
| 100 µm | 2.5 mM | 2.5 mM |
| 300 µm | 7.5 mM | 7.5 mM |

For comparison, D-76 contains 45 mM HQ, and 1 g/L KBr = 8.4 mM.
- So in dense areas the film sits in a Br⁻ concentration at the "restrainer" level. It also sees 6–17 % HQ depletion in D-76, about 33 % in D-76 1:1, and more in HC-110 dil. B.
- Across the swollen layer itself, τ = L²/D is about 0.15 s for Br⁻ and about 2 s for HQ (D_gel = 1 × 10⁻⁶ cm²/s). That is much shorter than τ ≈ δ²/D ≈ 10–90 s across the boundary layer.
- So the boundary layer, not the gelatin, sets compensation and adjacency.
- Local acid: about 0.5 M H⁺ released in a fully developed highlight cell, against a sulfite buffer capacity of about 0.1 M/pH in D-76. The local pH would fall by several units if exchange with the bath stopped. pH transport is therefore first-order physics, as the design states.

---

## 4. Stoichiometry, bromide release and capacity

### 4.1 Per-agent stoichiometry (balanced here; standard redox chemistry, recalled)

The proton count is relative to the **fully protonated neutral reference species**. That is what a total-component (TOTH) formulation needs (§5.2). The actual pH change then follows from speciation and buffering.

| Agent | Overall reaction per 2 Ag | e⁻ per agent | H⁺ released per Ag (neutral reference) | With sulfite scavenging | Confidence |
|---|---|---|---|---|---|
| Hydroquinone | H₂Q + 2 AgBr → Q + 2 Ag + 2 H⁺ + 2 Br⁻ | 2 | 1 | Q + SO₃²⁻ + H₂O → H₂Q-SO₃⁻ (hydroquinone monosulfonate) + OH⁻. Net: H₂Q + SO₃²⁻ + 2 AgBr → H₂QSO₃⁻ + 2 Ag + 2 Br⁻ + H⁺, i.e. **0.5 H⁺ per Ag** (atom and charge balance checked, and re-checked in the 2026-11-21 spot-check: C6 H6 O5 S1 Ag2 Br2, charge −2, on both sides) | **[unverified, spot-check 2026-11-21]** recalled reactions, with no primary source opened. The balance is correct, but whether real development follows this path is not confirmed. Contrast: US 4,025,344 A (Du Pont, lith developer replenishment) states "one mole of developing agent, hydroquinone, is needed to develop one mole of silver halide", i.e. 1 HQ per Ag in practice against the 0.5 redox minimum used in §3.6 and §4.3. The "HQ is not capacity-limiting" conclusion of §4.3 holds either way |
| Metol (free base) | HO-C₆H₄-NHCH₃ + 2 AgBr → O=C₆H₄=NCH₃ (quinone monoimine) + 2 Ag + 2 H⁺ + 2 Br⁻ | 2 | 1 | The monoimine is sulfonated (−1 H⁺ via OH⁻ release), or hydrolysed to benzoquinone + CH₃NH₂ (releases a base), or reduced back by HQ (regeneration, no net H⁺). Metol is supplied as the sulfate (C₇H₉NO)₂·H₂SO₄, which consumes 1 OH⁻ per metol on dissolution | recalled |
| Phenidone / Dimezone S | one-electron oxidation of the anion to the radical; regenerated by HQ or ascorbate. Irreversible two-electron end product: 1-phenyl-3-hydroxypyrazole (C₉H₁₀N₂O → C₉H₈N₂O + 2 H⁺ + 2 e⁻) | 1 (cycling); 2 (terminal) | 0 while regenerated; 1 per Ag when consumed | — | recalled. Confirm the terminal product in Mason or Haist |
| Ascorbic acid | H₂A + 2 AgBr → DHA + 2 Ag + 2 H⁺ + 2 Br⁻; DHA + H₂O → 2,3-diketogulonic acid (an acid) | 2 | 1, rising to **1.5** once DHA hydrolyses | sulfite plays no scavenging role for DHA (recalled). The Agfa EP 0 573 700 statement that development lowers the pH of ascorbate developers is consistent | recalled reactions; opened qualitative pH statement |
| p-Aminophenol (Rodinal) | HO-C₆H₄-NH₂ + 2 AgBr → O=C₆H₄=NH + 2 Ag + 2 H⁺ + 2 Br⁻ | 2 | 1 | Quinone imine hydrolysis releases NH₃ (a base). Rodinal has little sulfite, so expect more oxidation-product (stain) and less scavenging | recalled |
| Halide | 1 Br⁻ per Ag from AgBr. For AgBr₁₋ₓIₓ, (1 − x) Br⁻ + x I⁻ per Ag; iodide stays strongly adsorbed on the remaining grains (a stronger inhibitor than Br⁻) | — | — | — | recalled |
| Colour developer (comparison) | quinonediimine (QDI) + coupler anion → dye; QDI deamination is first order with k ∝ [OH⁻], second-order k ≈ 5–10 × 10³ L mol⁻¹ s⁻¹ for p-amino-N,N-diethylaniline; coupling is 10–10⁵× faster than deamination depending on coupler; K₉ (deamination equilibrium) = 0.7 × 10⁻⁴ at pH 11.1 | 2 | — | QDI + SO₃²⁻ → developer sulfonate | opened-source (Eggers 1961, §§3242–3243; citing Tong & Glesmann, J. Am. Chem. Soc. 78, 5827 (1956); 79, 583, 592, 4305 (1957); 82, 1988 (1960); Tong, J. Phys. Chem. 58, 1090 (1954)) |

### 4.2 Manufacturer capacity figures (opened-source; Kodak publications)

Kodak equates one 135-36 or 120 roll to 80 in² = **516 cm²** of film (J-109, Feb 2018, p. 2).

| Developer | Capacity, unreplenished | Time compensation | Minimum stock per roll | Source |
|---|---|---|---|---|
| D-76 full strength | **16 rolls per gallon = about 4 per litre** (tank) and 16 8×10 sheets per gallon (tray) | +15 % development time after every 4 rolls per gallon | — | J-78 (Dec 2017) p. 5, "Storage life and capacity", footnote 2 |
| D-76 1:1 | **8 per gallon (2 per litre)**, one-shot | when processing two 36-exposure rolls in a 16-oz (473 mL) tank, +10 % time. One roll per 473 mL or two per 946 mL is the standard; one roll in 237 mL or two in 473 mL needs +10 % | about 237 mL working solution per roll (≈ 118 mL stock) at the minimum | J-78 p. 2 and p. 5 |
| D-76 replenished | discard after 9600 in² per gallon (120 rolls per gallon); D-76R at 22.2–29.6 mL per roll (70 mL per roll of a modified 5:1 D-76 : D-76R replenisher for T-MAX films) | — | — | J-78 p. 6 |
| XTOL full strength | **≈ 15 rolls per litre** | rolls 1–5 normal; 6–10 +15 %; 11–15 a further +15 % on the adjusted time; discard after 15 | — | J-109 (Feb 2018) p. 2; same in J-109 (Apr 1998) and J-107 (1996) |
| XTOL 1:1 (and higher dilutions in 1998) | one-shot | — | **≥ 100 mL full-strength stock per roll**, "even if the processing equipment will allow the use of less solution"; with less, "may not contain enough active ingredients to develop the film fully" | J-109 (2018) p. 2; J-109 (Apr 1998) p. 2 (1:3 example: 400 mL; page corrected from p. 3 in the spot-check) |
| XTOL replenished | 70 mL per roll | — | — | J-109 (2018) p. 4 |
| HC-110 working solutions (tank, without replenishment) | dil. A 10, **dil. B 5**, C 8, D 4, E 3 sheets (8×10) per litre; tray B 2.5 per litre | — | — | J-24 (Dec 2017) p. 7, "Capacity" |
| HC-110 replenished | 22 mL replenisher per roll; discard after 2 weeks, or after replenisher equal to the original volume, or after 50 8×10 sheets per litre | — | — | J-24 p. 7 |
| T-MAX Developer (context) | 48 rolls per gallon (≈ 12.7 per litre) with time compensation | first 16 rolls normal, then steps | — | J-86 (Dec 2017) p. 9 |

### 4.3 What the capacity numbers constrain (estimates)

- **Silver per roll** = m_Ag × 0.0516 m² / 107.87 g mol⁻¹: 0.48 / 1.44 / 2.39 mmol for 1 / 3 / 5 g/m² coated (assumed range). The developed fraction of an average negative is perhaps 20–50 % (assumed), giving **~0.1–1.2 mmol Ag developed per roll**.
- **D-76, 4 rolls/L:**
  - HQ available is 45 mM, about 11 mmol per roll-share, against ≤ 0.6 mmol HQ consumed per roll. HQ consumption is therefore *not* what sets the capacity.
  - The accumulated Br⁻ is ≈ 0.4–5 mM after 4 rolls per litre, the same order as a deliberate restrainer addition.
  - Sulfite consumed by scavenging is ≤ 0.6 mmol per roll against 790 mmol/L.
  - Kodak's +15 % per 4 rolls per gallon (about 1 roll per litre) is therefore a seasoning correction. It lumps together bromide accumulation, the pH drop, metol loss, aerial oxidation and quarter-minute rounding. The estimate also depends on the assumed silver coverage. The **combined** finite-tank model, not K_Br alone, should reproduce a time increase of roughly 10–20 % per roll per litre. This is a consistency check, not a test that isolates K_Br.
- **XTOL ≥ 100 mL stock per roll** is the sharpest constraint on the finite-tank depletion model. With the (patent or msds) ascorbate concentration c_A of the XTOL data file, 100 mL holds 0.1·c_A mol ascorbate against ~0.05–0.6 mmol needed (0.5 ascorbate per Ag). The model must show visible under-development when the stock per roll falls below about 100 mL (a validation check, not a fit).
- **J-109 and Agfa's statement** that ascorbate developers lose activity as the pH drops ties capacity to the acid budget (1–1.5 H⁺ per Ag, §4.1) and the buffer capacity of the formula.

---

## 5. Numerics

### 5.1 Time scales (estimates; swollen emulsion 15–30 µm, cell Δz ≈ 1 µm, pixel 5–10 µm)

| Process | Time scale | Note |
|---|---|---|
| Acid–base proton transfer | ≪ 1 ms (diffusion-limited, ~10¹⁰ M⁻¹ s⁻¹, recalled) | treat as algebraic equilibrium |
| Semiquinone disproportionation; ETA-radical regeneration | ms to s (constants not found) | QSSA candidates |
| Quinone + sulfite → monosulfonate | fast relative to transport (LuValle; value not obtained) | QSSA or stiff |
| Br⁻ across the swollen layer, L²/D | 0.15–0.6 s | |
| HQ across the swollen layer, L²/D_gel | 2–9 s | |
| Boundary-layer diffusion, δ²/D | 10–250 s | comparable to the agitation interval (30/60 s) |
| Development (CI 0.56) | 200–600 s | |
| Lateral adjacency range in the liquid | 100–600 µm, i.e. 10–100 pixels | |
| Explicit stability limit, vertical (Δz = 1 µm, D = 2 × 10⁻⁵) | Δt < Δz²/2D = 0.25 ms | 10⁶ steps for 5 min: implicit is mandatory |
| Explicit stability limit, lateral (Δx = 6 µm, 2-D, D = 2 × 10⁻⁵) | Δt < Δx²/4D ≈ 4.5 ms | implicit, ADI or spectral |

### 5.2 Formulation: transport components, speciate algebraically

Do not integrate fast acid–base equilibria as ODEs. Use the standard reactive-transport formulation (Yeh & Tripathi, Water Resour. Res. 25, 93 (1989); Steefel & MacQuarrie, Rev. Mineral. 34 (1996); both recalled):
- **Transport totals**: total H (TOTH, proton excess relative to the neutral reference species of §4.1), total sulfite, total borate or carbonate, total of each developing agent, Br⁻, I⁻, oxidised products, and Na⁺/K⁺ as needed for electroneutrality.
- **Speciation**: in each cell after each step, solve the proton balance for pH (one scalar Newton per cell; monotone, bracketable). This gives the ionised fractions a_i needed by §1.10.
- **Gelatin**: an immobile buffer (Iwano's 4.8 × 10⁻⁴ mol/g OH⁻-reactive groups) enters TOTH but not transport.
- **Diffusion**: if each species diffuses with its own D (HQ²⁻ ≠ H₂Q, for example), the totals do not obey a single diffusion equation. Practical fix: use one D per component (species-averaged, or that of the dominant species), tagged `assumed`. Use Nernst–Planck with electroneutrality only if a test shows that the liquid-junction effect changes the CI by more than the target error.

### 5.3 Splitting scheme per time step Δt (Strang)

1. **Lateral half-step (Δt/2):** 2-D diffusion within each depth layer (emulsion cells and boundary-layer cells), each component separately.
   - Linear and constant-coefficient per layer: use a DCT-based spectral solve (Neumann boundaries at tile edges inside the halo), which is exact for the linear operator.
   - Alternatively use Peaceman–Rachford ADI (J. SIAM 3, 28 (1955); recalled), which is unconditionally stable and second order.
2. **Column step (Δt):** per pixel column, vertical diffusion (emulsion ↔ overcoat ↔ boundary-layer cells ↔ tank Dirichlet value), together with all reactions and per-class grain ODEs.
   - This is a block-tridiagonal stiff system: species × cells, about 10 × 40.
   - Integrate it **with the reactions coupled, not split**. The development rate depends on concentrations that diffusion sets on the same time scale near the surface, and splitting reaction from vertical diffusion produces errors at the reaction front.
3. **Lateral half-step (Δt/2).**
4. **Tank coupling:** add the area-integrated flux into the well-mixed tank ODE explicitly once per step (the tank time scale is ≫ Δt).

Splitting error is O(Δt²) (Strang, SIAM J. Numer. Anal. 5, 506 (1968); Hundsdorfer & Verwer, *Numerical Solution of Time-Dependent Advection-Diffusion-Reaction Equations*, Springer 2003, ch. IV; recalled). Lateral coupling is weak, since gradients are over 10–100 pixels, so the splitting error is small. Check it by halving Δt.

### 5.4 Integrator for the column step

| Option | Pros | Cons | Reference (recalled) |
|---|---|---|---|
| **Rosenbrock W-methods, ROS2 or ROS3P** (linearly implicit, one block-tridiagonal LU per stage, embedded error estimate) | No Newton iteration, which makes it easy to **vectorise across 10⁶–10⁷ columns** in NumPy (batched block-Thomas). ROS2 was designed for atmospheric ADR. L-stable (ROS2) or A-stable (ROS3P) | Needs an exact or good Jacobian. ROS3P is not L-stable | Verwer, Spee, Blom & Hundsdorfer, SIAM J. Sci. Comput. 20, 1456 (1999) (ROS2); Lang & Verwer, BIT 41, 731 (2001) (ROS3P) |
| BDF2 / backward Euler with Newton | robust, L-stable (BE) | the Newton iteration count varies per column, which vectorises badly | Hundsdorfer & Verwer 2003 |
| SciPy `solve_ivp` (BDF / Radau / LSODA) per column | reference-quality | one Python call per column is far too slow for full frames. Use it only as the **reference solver** in unit tests on a few columns | SciPy documentation |
| Exponential integrators (ETD, exponential Rosenbrock) | exact for linear diffusion | the nonlinear coupling and block structure make them less attractive than Rosenbrock here | Hochbruck & Ostermann, Acta Numerica 19, 209 (2010) |
| QSSA for radicals (P˙, H˙⁻) and sulfonation | removes the fastest modes | introduces algebraic closures. Validate against the full system on test columns | standard (atmospheric chemistry practice) |

Recommendation: ROS2 with an analytic block-tridiagonal Jacobian, a batched solve across columns, an adaptive Δt from the embedded estimate with a global step for the tile, and float64. Keep the speciation (§5.2) outside the Rosenbrock stage, as a function of the totals. Its derivatives then enter the Jacobian via the chain rule or by finite differences on the scalar pH solve.

### 5.5 Events and boundaries

- **Agitation:** at each agitation time, reset the boundary-layer cells to the tank composition, or reset δ to δ_min (§3.4). Align the time steps to the event times and restart the integrator after each event; this is a discontinuity in the data.
- **Immersion start:** the film enters dry or prewetted. Prewet means the layer starts swollen with water. Dry film means a swelling delay of order T½ (§3.3).
- **Tank:** well mixed, finite volume. The composition evolves from the integrated film flux, plus aerial oxidation (a fitted first-order loss, tagged `assumed`).
- **Stop bath / fixer handoff:** the end state carries the remaining totals into the next step, as the design requires.

### 5.6 Tiling

The halo must exceed the lateral influence length. Two cases:
- **Agitated processes.** Agitation resets the liquid boundary layer, which bounds lateral spread in the liquid. The gel layer exchanges vertically with the boundary layer within about 0.15–2 s, so it does not carry lateral spread past a reset. With 3–5 × √(D·t_int), D = 2 × 10⁻⁵ cm²/s, the length is about 0.7–1.2 mm for t_int = 30 s and 1.0–1.7 mm for t_int = 60 s. This is an estimate that assumes complete reset at each agitation.
- **Stand development** (no resets). 3–5 × √(D·t_dev) with t_dev = 400 s gives **2.7–4.5 mm**.

Recommendation: start with a halo of 2–3 mm (about 300–500 pixels at 6 µm) for agitated processes and about 5 mm for stand development. Shrink it only when the overlapping-tile seam test passes (the design already lists "分块无接缝"). This is larger than the "> diffusion range" the design suggests by eye.

### 5.7 Verification tests with analytic or reference solutions

1. **1-D diffusion into a slab** from a constant-concentration boundary: erfc solution, and series solution for finite thickness (Crank, *The Mathematics of Diffusion*, 2nd ed., 1975; recalled).
2. **Steady reaction-diffusion with first-order consumption** in a slab: cosh profile; the Thiele modulus gives the effectiveness factor.
3. **Moving reaction front** (instantaneous reaction, immobile reactant concentration C_s ≫ supply C₀, stoichiometry ν): X(t) ≈ √(2·D·C₀·t/(ν·C_s)). This is the regime of Iwano's "HQ penetration under diffusion-controlled development" (I-1972 Fig. 3).
4. **Pure buffer titration:** adding acid to a cell gives the pH from the analytic Henderson–Hasselbalch / β integral.
5. **Conservation:** total Ag (developed + halide + dissolved complex), total developer (all oxidation states), total Br, and TOTH to machine precision per step.
6. **Splitting convergence:** O(Δt²) under halving; lateral spectral against ADI.
7. **Agitation limit:** with δ → 0, results must converge to the "infinitely agitated" case. With no agitation, they must show the expected (πDt)^½ depletion.

### 5.8 Comparable published simulations

| Work | What it is | Confidence |
|---|---|---|
| Nelson (1971), Photogr. Sci. Eng. 15, 82 | empirical chemical spread function convolved with exposure/density | opened as citation (Rajkowski & Nowak) |
| Rajkowski & Nowak (1998) | Frieser-form LSF + ChSF superposition to compute edge curves; no reaction-diffusion | opened-source |
| Iwano (1970–72) | analytic penetration curves for HQ under diffusion-controlled development, with an indicator-layer method | opened-source |
| Bayer & Hamilton (1965), JOSA 55, 439 | Monte-Carlo of latent-image formation giving D–log E curves (exposure side, not development) | secondary-source |
| Kriss, Nelson & Eisen (diffusion model of adjacency effects) | **not found** | — |
| Filmulator (open-source raw processor) | simulates developer diffusion with a simple reservoir model. Prior art only; **no data or code is taken from it** | recalled |

No published per-pixel stiff reaction-diffusion simulation of B&W development with pH transport was found. D1 would be new in this respect.

---

## 6. Measured in the literature versus must be fitted

| Quantity | Status | Best source found |
|---|---|---|
| pKa of HQ | literature (spread 0.6 in pKa₂) | compilations (secondary); confirm one |
| pKa of ascorbate, metol, p-aminophenol, sulfite, borate | literature, but not opened here | recalled → CRC / IUPAC |
| pKa of phenidone / Dimezone S | **not found** | — |
| Active ionised form: HQ (dianion), ascorbate (mono- and dianion) | literature (qualitative) | Eskes 1982 (abstract opened); Willis & Pontius 1970 via EP 0 573 700 ("likely") |
| Reaction order in [HQ] ≈ 0.6 (pH 8–8.9) | literature | James 1946 (opened abstract) |
| Ea, metol / HQ development | **literature: 0.62–0.70 eV (60–68 kJ/mol)** | Sugai & Furuichi 1953 (opened) |
| Ea, phenidone / Dimezone S / ascorbate / p-aminophenol | **fitted** (not found) | — |
| Apparent Ea from Kodak tables | derived check (full strength ≈ 66–68 kJ/mol; 1:1 ≈ 41–43, on Kodak's °F temperatures). One observation per dilution class: the D-76 and XTOL rows share a slope, so this is possibly a Kodak convention | F-4017, J-109, F-4016 (estimate) |
| Rate constants k_i, adsorption constants K_i, K_Br, K_I | **fitted** | — (structure only: James 1944/1946; Langmuir isotherm, sec.) |
| Induction-period barrier terms, latent-image-size dependence | **fitted** (structure: Sugai 1953; Ohzeki 1996) | — |
| Regeneration k_reg, ETA loss k_loss | **fitted** (not found) | — |
| Quinone + sulfite k_sulf | literature exists (LuValle 1952), not obtained | — |
| AgBr solubility in sulfite | literature (low-sulfite range only; 7.5 × 10⁻⁵ N at 0.010 N sulfite, 3.5 × 10⁻⁵ N in gelatin emulsion) | Carroll & Hubbard 1934 (opened) |
| Ag–sulfite complex constants | literature, not obtained | Smith & Martell (recalled) |
| D in water (ions) | literature | CRC (recalled) |
| D in water (organic agents) | estimate (Wilke–Chang) | — |
| D in gelatin: HQ, metol, DEPPD, salts, OH⁻, HCHO | **literature (measured)** | Iwano 1970, 1972 (opened) |
| D in gelatin: phenidone, Dimezone S, ascorbate | **assumed** ratio 0.1–0.3 × water | — |
| Gelatin OH⁻-reactive group content | literature: 4.8 × 10⁻⁴ mol/g | Iwano 1970 (opened) |
| Swelling ratio of Tri-X / T-Max 100 | **assumed** α = 3 (2–5) | Iwano's "normally hardened"; colour-film patents (sec.) |
| Emulsion thickness, coated silver | **not found** (fitted or assumed; see latent-image.md) | — |
| Boundary-layer thickness, small hand tank | **fitted** (one δ_max), prior 200–500 µm | estimate; processor values 6–104 µm (US 5,315,338, opened) |
| Adjacency length | emerges from the model; validated on MTF > 100 % | Nelson 1971 (citation); MTF curves (datasheets) |
| Stoichiometry (e⁻, H⁺, Br⁻ per Ag) | textbook redox (recalled), balanced here | Mason / James *TPP* to confirm |
| Capacity, minimum stock per roll | **datasheet** | J-78, J-109, J-24, J-86 (opened) |
| Fog rate | **fitted**; constraint: base + fog spread ≤ 0.05 across 6–12 min (Tri-X D-76) | datasheets.md vectors |

---

## 7. Consequences for D1

1. **Per-species gel diffusion.** Developing agents diffuse about 10× slower in swollen gelatin than in water, while salts slow by about 1.3× (Iwano). The species library needs D_water and D_gel/D_water per species (or per class), not one tortuosity factor.
2. **Arrhenius per agent from literature for metol and HQ only.** Use Sugai & Furuichi's 0.64–0.70 eV as `literature` priors. Their AgI trend (0.62 → 0.66 eV from 0 to 7 mol% I) is usable. The 0.91 eV at 40 mol% comes from a different developer and exposure, so it is not an AgI effect alone. Their own tabulated times refit to 0.67–0.80 eV, so give the prior a range of ±0.1 eV.
   - The phenidone, Dimezone S and ascorbate Ea must be `fitted` or `assumed`.
   - The Kodak XTOL time–temperature rows copy the D-76 slope (§1.6), so they do not identify these Ea independently. Expect weak identifiability, and give the priors wide ranges.
3. **Possible internal check on transport.** Kodak's full-strength class gives an apparent Ea of about 66–68 kJ/mol; the 1:1 class gives about 41–43 (on Kodak's °F temperatures; the rounded °C labels give about 63 and 40).
   - This is one observation per dilution class, and it may be a Kodak convention.
   - If the model's reaction Ea is set by Sugai (metol/HQ), D-76 1:1 should show a lower apparent slope through boundary-layer transport and depletion.
   - Hold the 1:1 temperature rows out of calibration as a within-pair check, as the design's "稀释、温度、时间的变化也留出一部分只做检验" asks. Report a mismatch as "model or convention", not as a failure of the physics.
4. **Split hygiene.** The F-4016 HC-110 temperature/contrast factors and the T-Max 100 rows of J-109 Table 1 belong to held-out pairs. They were not used for any derived number here. Keep them out of every fit.
5. **pH transport is not optional.** The local acid release in dense areas (~0.5 M H⁺ in a 15 µm cell at 3 g/m² Ag, estimate) exceeds the D-76 buffer capacity (~0.1 M/pH, mostly sulfite) by several pH units. The species library must carry sulfite/bisulfite as a buffer, and the gelatin as an immobile buffer. Borate is a minor buffer in D-76.
6. **The boundary layer is the main compensation and adjacency control.** One fitted δ_max (with a growth law (πDt)^½ between agitations) shared by all pairs carries the small-tank versus large-tank difference. This also handles the D-76 agitation mismatch between the Tri-X (large tank) and T-Max 100 (small tank) calibration curves. Fit it on Tri-X D-76 small-tank CI(t) against large-tank time tables (training data only).
7. **Superadditivity constants are wholly `fitted`.** No regeneration rate constant was found. To keep D-76 and XTOL identifiable, calibrate k_reg separately for metol/HQ (D-76) and for Dimezone S/ascorbate (XTOL). Accept that HC-110 (msds-estimated formula) will constrain little.
8. **Capacity data give cheap validation checks.**
   - D-76: +15 % time per ~1 roll per litre (bromide seasoning).
   - XTOL: ≥ 100 mL stock per roll for full development.
   - XTOL: 15 rolls per litre with two +15 % steps.

   All are datasheet values for training developers. They test the combined finite-tank model (bromide, pH, agent loss, aerial oxidation together) without touching the held-out pairs. They do not isolate K_Br or any other single constant. They need the silver coverage, which is still `assumed`. Report them as consistency checks, not validations.
9. **Growth law and induction are model-structure choices.** The literature does not settle area-law autocatalysis against linear filament growth (ISJ 1998). Make g(m) switchable and decide on the Tri-X D-76 time series.
10. **Numerics.**
    - Transport components with algebraic speciation.
    - Strang splitting: lateral spectral or ADI, and a coupled column step (ROS2, batched across pixels).
    - Event-aligned agitation.
    - Tile halo 2–3 mm for agitated processes and about 5 mm for stand development, shrunk only after the seam test passes.
    - Add the verification tests of §5.7 to the design's unit-level checks.
11. **Rodinal (predicted).** Its chemistry (p-aminophenol, high pH, little sulfite) needs a p-aminophenol Ea (not found → `assumed`, borrowed from metol with a flag). It also needs the quinone-imine hydrolysis branch (base-releasing), which no training developer exercises. Its status must stay `predicted`, with that caveat stated.

---

## 7a. Named sources searched but not found or not opened

| Source named in the brief | Result |
|---|---|
| Fortmiller & James (hydroquinone development kinetics) | **not found**. Exact-author searches returned nothing, and no citation was seen in any accessible reference list |
| Lee & James | **not found** (same) |
| Willis & Pontius | citations only (Photogr. Sci. Eng. 14, 141 and 14(6), 384 (1970); 16, 406 (1972); 17, 21 and 157 (1973)). The ascorbate-species result comes via Agfa EP 0 573 700. No numbers |
| Tong & Glesmann (colour developers) | originals not opened. Their kinetic picture and constants are summarised in Eggers 1961 (opened), §4.1 |
| Levenson (superadditivity) | regeneration mechanism attributed to him in secondary sources; J. Photogr. Sci. 24, 68 (1976) cited. **No rate constants found** |
| Willems & Van Veelen | citations only (Photogr. Sci. Eng. 6, 39 and 49 (1962)); subject not confirmed |
| Jaenicke (electrochemical model) | 1972 Croat. Chem. Acta conference paper and a 1977 review cited (secondary). The quantitative mixed-potential model was not obtained |
| James (ed.) *TPP* 4th ed. ch. 13–14; Mason; Haist | not opened (no accessible full text) |
| Kriss, Nelson & Eisen diffusion model of adjacency effects | **not found** |
| Simonds (adjacency) | **not found** |
| Eberhard / Kostinsky original length scales | **not found** in a citable source |
| Nelson 1971 chemical spread function widths | citation confirmed; widths **not found** |
| Boundary-layer thickness for hand-inverted tanks or rotary tubes | **not found** (estimate given in §3.4) |
| Diffusion coefficients of phenidone, Dimezone S and ascorbate in gelatin | **not found** |
| Activation energies for phenidone, Dimezone S, ascorbate and p-aminophenol development | **not found** |
| Regeneration (superadditivity) rate constants | **not found** |
| Fog development rate law; quantitative grain-size dependence of development rate | **not found** |

## 8. What would confirm the recalled or secondary items

| Item | Where to confirm |
|---|---|
| James 1939/1940 (J. Phys. Chem. 43, 701; 44, 42) and J. Chem. Phys. 11, 338 (1943): charge, induction period, HQ kinetics | the ACS/AIP full text; James (ed.) *TPP* 4th ed., development chapters |
| Willis & Pontius 1970 (Photogr. Sci. Eng. 14(6), 384) ascorbate species; Pontius & Willis 14, 141 (1970); 16, 406 (1972); 17, 21 and 157 (1973) | IS&T archive |
| Axford & Kendall 1954; Levenson J. Photogr. Sci. 24, 68 (1976); Green & Levenson 24, 193 (1976) | Taylor & Francis full text |
| Willems & Van Veelen, Photogr. Sci. Eng. 6, 39 and 49 (1962): subject and any rate numbers | IS&T archive |
| K. Ohzeki 1996 rate data versus latent-image size (abstract now opened; the full text is needed for numbers) | J. Imaging Sci. Technol. 40, 591 |
| Langmuir bromide adsorption (~10 % coverage) | Herz & Helling, J. Colloid Sci. 17, 293 (1962), DOI 10.1016/0095-8522(62)90009-0; Gardner, Wrathall & Herz, Photogr. Sci. Eng. 21, 325 (1977) for the standard state of the opened Herz 1980 ΔG values |
| H⁺ and HQ per Ag in real development (0.5 redox minimum against 1 per Du Pont US 4,025,344) | Mason; James *TPP* 4th ed. |
| ISJ 1998 linear-kinetics paper (authors) | DOI 10.1080/13682199.1998.11736441 |
| LuValle quinone–sulfite rate constants | J. Am. Chem. Soc. (1952) |
| Nelson 1971 chemical spread function widths; Kriss, Nelson & Eisen; Simonds | Photogr. Sci. Eng. 15, 82 (1971) and the IS&T archive; James *TPP* 4th ed., chapter on adjacency effects (Rajkowski & Nowak cite "Chapt. 21") |
| pKa values (HQ, ascorbate, metol, p-aminophenol, phenidone, sulfite, borate) | CRC Handbook; IUPAC stability-constant database |
| Ionic D in water (Br⁻, OH⁻, SO₃²⁻, B(OH)₄⁻) | CRC Handbook, "Ionic conductivity and diffusion at infinite dilution" |
| D-76 formula and pH | Kodak J-1 (Processing Chemicals and Formulas); Anchell & Troop |
| Stoichiometry and terminal products (phenidone → hydroxypyrazole; DHA hydrolysis) | Mason, *Photographic Processing Chemistry*; Haist, *Modern Photographic Processing* |
| Nietz reduction-potential scale and bromide-shift law | Nietz, *The Theory of Development* (1922) |
| Numerical references (Strang 1968; Peaceman & Rachford 1955; Verwer et al. 1999; Lang & Verwer 2001; Hochbruck & Ostermann 2010; Hundsdorfer & Verwer 2003; Yeh & Tripathi 1989; Steefel & MacQuarrie 1996; Crank 1975) | the cited journals and books (all standard; recalled bibliographic details) |

## Sources opened in this job (URLs)

- T.H. James, J. Chem. Phys. 14, 536 (1946): https://doi.org/10.1063/1.1724190 (abstract via pubs.aip.org)
- T.H. James, J. Chem. Phys. 12, 453 (1944): https://pubs.aip.org/aip/jcp/article-abstract/12/11/453/183897 (abstract)
- S. Sugai & J. Furuichi (1953): https://hdl.handle.net/2115/34215
- S. Sugai (1953): https://hdl.handle.net/2115/34210
- H. Iwano (1970): https://www.jstage.jst.go.jp/article/photogrst1964/33/3/33_3_143/_article/-char/en
- H. Iwano (1972): https://doi.org/10.11454/photogrst1964.35.98
- B.H. Carroll & D. Hubbard (1934), RP659: https://nvlpubs.nist.gov/nistpubs/jres/12/jresv12n3p329_A2b.pdf
- B.H. Carroll & D. Hubbard (1930), RP183: https://nvlpubs.nist.gov/nistpubs/jres/5/jresv5n1p1_A2b.pdf
- B. Rajkowski & P. Nowak, Optica Applicata 28(2), 105 (1998): https://dbc.wroc.pl/Content/41387/PDF/optappl_2802p105.pdf
- J. Eggers, Chimia 15, 499 (1961): https://www.chimia.ch/chimia/article/download/1961_499/7358/25716
- EP 0 573 700 A1 (Agfa-Gevaert): https://data.epo.org/publication-server/rest/v1.2/patents/EP0573700NWA1/document.pdf
- US 5,315,338 A: https://patents.google.com/patent/US5315338A/en
- US 5,837,437 A (Eastman Kodak): https://patents.google.com/patent/US5837437A/en
- Ives & Jensen, J. SMPE 40(2) (1943): https://journal.smpte.org/periodicals/Journal%20of%20the%20Society%20of%20Motion%20Picture%20Engineers/40/2/7/
- Snyder, J. SMPTE 75(10) (1966): https://journal.smpte.org/periodicals/Journal%20of%20the%20SMPTE/75/10/9/
- Kodak J-78 (2017), J-109 (1998, 2018), J-107 (1996), J-24 (2017), J-86 (2017), F-4016 (2018), F-4017 (2021): see datasheets.md §1 for URLs

## Spot-check log (2026-11-21, independent re-check)

Sources were re-fetched from the publishers (J-STAGE, HUSCAP, NIST, EPO, Crossref, RIT repository, web.archive.org copy of the IS&T page, Hrčak). Tables were read from freshly rendered page images. Copies are in `.cache/datasheets/spotcheck_kin/` (git-ignored).

| Claim | Verdict | Note |
|---|---|---|
| Iwano gel D: HQ 0.99, metol 1.17, DEPPD 2.90 × 10⁻⁶ cm²/s (8 % gel) | confirmed | I-1970 Table 5 (23 °C, 2 × 10⁻² M); I-1972 Table 3 repeats the values but labels them 25 °C |
| Iwano coated-layer D, α = 3.03 (I-1970 Table 7) | confirmed | dry 1.85 / 1.28 / 0.88 × 10⁻⁷; swollen 16.9 / 11.7 / 8.0 × 10⁻⁷; ratio ≈ α² |
| Iwano salt D in 8 % gel (I-1972 Table 4, §5) | confirmed | including KBr 3.03 × 10⁻⁵ (4 %) and 1.10 × 10⁻⁵ (20 %) |
| Iwano OH⁻ true and apparent D, 4.8 × 10⁻⁴ mol/g | confirmed, ratio corrected | true/apparent is 16–106×, not 5–100× |
| Sugai & Furuichi Ea 0.64 / 0.70 / 0.70 eV | confirmed | KBr developer has 2 g metol; the tabulated times refit to 0.80 eV for KBr (a discrepancy in the paper itself) |
| Sugai & Furuichi AgI series 0.62 / 0.63 / 0.66 / 0.91 eV | confirmed, caveat added | the 40 mol% plate used a different developer and exposure |
| Kodak apparent Ea (≈ 63 / ≈ 40 kJ/mol) | wrong (temperature basis) | ≈ 66–68 / ≈ 41–43 kJ/mol on °F temperatures; t₁₈/t₂₄ = 1.68 confirmed |
| James 1946 order 0.6 (HQ) and 0.56 (catechol) | confirmed | Crossref abstract |
| Eskes 1982, divalent ion only | confirmed | RIT abstract opened; contributions "not possible" to separate |
| Carroll & Hubbard 1934 AgBr solubility | wrong (two exponents) | 3.5 × 10⁻⁵ and 4 × 10⁻⁵ N, not 10⁻⁶; 2.9 × 10⁻⁵ confirmed |
| US 5,315,338 h = 0.0018–0.0086 cm/s, D = 5 × 10⁻⁶ | confirmed | Eastman Kodak; 0.00048 is "84 percent less" than 0.0030 |
| Kodak capacity and minimum-stock figures (J-78, J-109, J-24) | confirmed | J-109 (1998) minimum-stock statement is on p. 2 |
| HQ + sulfite stoichiometry, 0.5 H⁺ per Ag | unverified | balance correct; Du Pont US 4,025,344 says 1 HQ per AgX in practice |
| Ohzeki 1996 ≥ 5-atom threshold | confirmed (abstract) | author is Katsuhisa Ohzeki; S+Au emulsion |
| Br⁻ Langmuir, ~10 % coverage (J. Colloid Sci. 1962) | unverified | Herz & Helling 1962 identified; the opened Herz 1980 limiting areas imply ~14–18 % |
| Agfa EP 0 573 700 ascorbate species (not in the list above) | confirmed | wording is "likely" |
| F-4017 (2021) p. 8 MTF | confirmed | re-read at 300 dpi: ~106 % at 2.4 c/mm, ~111 % peak near 10, 100 % at ~27, ~25 % at ~71 c/mm |
