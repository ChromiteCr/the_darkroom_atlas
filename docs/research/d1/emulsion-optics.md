# D1 research: light transport in the emulsion, halation, image structure and density

Scope: stage 1 (optics), stage 4 (silver density) and stage 6 (grain) of the physics chain, for Kodak Tri-X 400 (400TX) and Kodak T-Max 100 (TMX). Date of research: 2026-11-21.

Independent spot-check 2026-11-21 (second reviewer). **Re-done from the cached publisher PDFs:** both MTF curves were re-digitised from the vector paths, with the axes mapped through Kodak's own printed gridlines and minor ticks. The T-Max curve agrees within about 1 point; the identical path, figure code F002_0542AC, is in the 2016 and the June 2018 editions. The Tri-X curve agrees up to its 50 % point, but the 60–70 c/mm tail depends on the axis mapping (corrected below). **Confirmed in the source text:** granularity 8 and 17 with their read conditions; every base-thickness statement in §1.2; the H-24 grey-base quote; the EP 0 513 726 covering-power definition and its bibliographic data; the ISO 5-3:2009 visual-density wording on p. 7, whose S_H is specified in §4.2.3, not §4.2.1; the Newson IPOL model details (Eq. 1–2, σ = 0.8, δ = 1/⌈1/μ_r⌉, 0.999 quantile, 40 s / 2.4 s); the US 12,705,810 B2 front page and Eq. 4, 5, 13 and 16–18. Every DOI in the report was checked against Crossref (Nutting, Kelly I/II, Bayer, Schröter, Gilmore, Wolfe et al., Vendrovsky, Farnell–Solman 1963/1970/1973, Saunders, Zhang, Newson ×2); titles, volumes and pages match. **Recomputed and reproduced:** the AgBr dispersion values (from the database formula, coefficients re-read), the Fresnel/TIR/halation geometry, Siedentopf areas, Selwyn S, κ(D) (1.055/1.18/1.34/1.89/2.86/4.66; lognormal CV 0.5 gives 1.57), the Boolean areas and the covering-power table. **Upgraded:** TAC index (patent source found) and the Frieser form and Nelson citation (Nowak 1998). **Still unverified:** the gelatin index, the interpretation of the 35 mm − 120 base + fog offset as grey-dye density, and the internal equations of Kelly's three-stage model. Changes are marked "[spot-check]".

Confidence tags used below:

- **opened-source**: I opened the source and read the number.
- **secondary-source**: the citation was checked (Crossref metadata, a library record, or a quote in another primary document), but the paper itself was not opened.
- **recalled**: textbook knowledge, not re-checked here. The table says what would confirm it.
- **estimate**: my own calculation from opened numbers. The script logic is described so it can be redone.

Files kept for reading (git-ignored): `.cache/datasheets/pdf_optics/f4016_2018_wb.pdf` (T-Max 100 F-4016, February 2016 edition, from web.archive.org snapshot of imaging.kodakalaris.com), `.cache/datasheets/pdf_optics/f4017_kp.pdf` (Tri-X F-4017, October 2021 edition, © 2021 Kodak Alaris, from kodakprofessional.com), `.cache/datasheets/kodak_mp/h24_module15.pdf` (Kodak H-24 Module 15), `.cache/datasheets/kodak_mp/doublex.pdf`, `.cache/datasheets/papers/newson_ipol2017.pdf`, `.cache/datasheets/papers/us12705810.pdf` (+ page renders), `.cache/datasheets/papers/iso5-3_sample.pdf`, `.cache/datasheets/papers/ep0513726.pdf`.

---

## 0. Summary for the simulator

| Quantity | Value to use | Tag | Where |
|---|---|---|---|
| AgBr refractive index | 2.310 (500 nm), 2.276 (550), 2.257 (589), 2.237 (650), 2.232 (670); formula valid 495–670 nm only | opened-source (data file) | §1.1 |
| Gelatin (dry) refractive index | ≈ 1.53–1.54 | recalled, **unverified** [spot-check: no primary source opened] | §1.1 |
| Cellulose triacetate base index | ≈ 1.47–1.48 | 1.48 confirmed in US 6,398,371 B1 (Fuji Photo Film; wavelength not stated) [spot-check] | §1.1 |
| PET (Estar) index (sheet film only; not used for 135/120) | in-plane ≈ 1.64–1.70, normal ≈ 1.54–1.55 (biaxially oriented); isotropic 1.573 at 550 nm | secondary / opened | §1.1 |
| 135 base, T-Max 100 and legacy Tri-X Pan (TX) | 5 mil acetate = 127 µm (Kodak prints 0.13 mm). Not stated for current 400TX [spot-check: "both films" narrowed] | opened-source | §1.2 |
| T-Max 100 120 base | 4.7 mil acetate = 119 µm (Kodak prints "0.11 mm") | opened-source | §1.2 |
| Tri-X 120 base | 3.6 mil acetate = 91 µm (old TX, F-9 2003); 400TX current base **not stated** | opened-source | §1.2 |
| Grey-dye density of 35 mm base | ≈ 0.10 above 120 (diffuse visual), from the 35 mm vs 120 base + fog difference | estimate. [spot-check] The ≈ 0.1 offset reproduces (D-76 charts: ≈ 0.28–0.34 vs ≈ 0.18–0.22). That the offset is grey dye is **unverified**; no Kodak source states a grey base for still film | §1.3 |
| Halation ring inner radius, 135 TAC | 2·t·tan θc = 234 µm (t = 127 µm, θc = 42.7°) | estimate | §3 |
| Halation ring inner radius, 120 | 168 µm (Tri-X, 91 µm), 220 µm (T-Max, 119 µm) | estimate | §3 |
| Fraction of Lambertian light in gelatin that is TIR at base–air | 1 − 1/n_gel² ≈ 0.58 | estimate | §3 |
| T-Max 100 MTF (D-76, 20 °C) | 104 % at 2.5 c/mm, peak ≈ 113–114 % at ≈ 19 c/mm, 100 % at 50, 64.5 % at 100, 39 % at 150 c/mm | opened-source (vector digitised; spot-check re-digitisation agrees within ≈ 1 point) | §2.4 |
| Tri-X MTF (D-76 large tank) | 106 % at 2.5, peak 111 % at 9–10 c/mm, 100 % at ≈ 23–24, 55 % at 50 c/mm, 50 % at 52.5–53 c/mm; **tail 36–39 % at 60, 24–26 % at 70 c/mm; curve ends at 70–72 c/mm** [spot-check: the tail depends on the axis mapping] | opened-source (vector digitised) | §2.4 |
| RMS diffuse granularity | Tri-X 400 = 17 (HC-110 B, 20 °C, old emulsion); T-Max 100 = 8 (D-76, 20 °C); both at net diffuse D = 1.0, 48 µm aperture, 12× | opened-source | §5.2 |
| Siedentopf effective grain area at D = 1 | Tri-X 1.20 µm², T-Max 0.27 µm² | estimate | §5.3 |
| Boolean (exact overlap) effective grain area at D = 1 | Tri-X 0.64 µm², T-Max 0.14 µm² (constant-radius disks) | estimate | §5.4 |
| Overlap variance factor κ(D) = Var_Boolean / Var_Siedentopf | 1.06 (D 0.1), 1.18 (0.3), 1.34 (0.5), 1.89 (1.0), 2.86 (1.5), 4.66 (2.0) | estimate | §5.4 |
| ISO visual density spectral product | V(λ) × CIE illuminant A (2856 K), ISO 5-3 §6.1 | opened-source | §4.4 |

### Consequences for D1

1. **Statistical mode needs κ(D).** The binomial statistical mode will not match the Boolean micro mode at D ≈ 1 unless it carries the overlap factor κ(D) (§5.4). Without it, σ_D is too low by √κ ≈ 1.37 at D = 1, and the D1 criterion "statistical and micro modes agree on granularity" fails.
2. **The holdouts have no image-structure reference.** All image-structure reference data sits in training pairs: T-Max granularity and MTF in D-76, Tri-X MTF in D-76, Tri-X granularity in HC-110 (B). The holdout pairs T-Max × XTOL and T-Max × HC-110 have no MTF or granularity reference, so those rows can only be `predicted`.
3. **Design decision for the user: Tri-X granularity conflicts with the D-76-only rule.** Tri-X granularity 17 was measured in HC-110 (B). The design says film parameters are calibrated on D-76 only, so the Tri-X grain classes cannot be fitted to 17 without using HC-110 data. If the rule stands, 17 becomes a cross-developer check rather than a calibration target. This report does not resolve that choice.
4. **Some Tri-X data is from older emulsions.** The Tri-X MTF comes from the older TRI-X Pan 5063 film, and the granularity from the old TX/TXP emulsions (both stated in F-4017).
5. **Base data per format.** 135 and 120 differ in base thickness (127 µm vs 91 µm Tri-X / 119 µm T-Max) and in grey dye (≈ 0.10). The material model therefore needs base data per format. The ≈ 0.1 base + fog offset between 35 mm and 120 in the datasheets becomes a check on the base model.
6. **The optical PSF is not fitted to the MTF directly.** The optical PSF must not be fitted to the datasheet MTF directly (§2.4, as the design already says).
7. **The MTF measurement needs an assumed mean density.** Kodak does not state the mean density or exposure of its MTF measurement. The virtual sine-target measurement has to copy Kodak's effective-exposure method, with the mean density tagged `assumed`.

---

## 1. Optical properties of the layers

### 1.1 Refractive indices

| Material | n | λ | Source | Tag |
|---|---|---|---|---|
| AgBr | 2.3095 | 500 nm | H. Schröter, "Über die Brechungsindizes einiger Schwermetallhalogenide im Sichtbaren …", *Z. Phys.* **67**, 24–36 (1931), doi:10.1007/BF01391040. Values computed from the dispersion fit in the refractiveindex.info database file `main/AgBr/nk/Schroter.yml`: (n²−1)/(n²+2) = 0.452505 + 0.09939 λ²/(λ²−0.070537) − 0.000150 λ² (λ in µm), valid 0.495–0.67 µm | opened-source (data file); paper not opened |
| AgBr | 2.2756 | 550 nm | same | same |
| AgBr | 2.2572 | 589.3 nm | same | same |
| AgBr | 2.2531 | 600 nm | same | same |
| AgBr | 2.2372 | 650 nm | same | same |
| AgBr | 2.2322 | 670 nm | same | same |
| AgBr, 400–490 nm | 2.37 (450 nm), 2.51 (400 nm) | — | Polyanskiy et al., *Opt. Mater. Express* **14**, 696 (2024), a fit to Schröter + McCarthy 1973 data (database file `Polyanskiy.yml`). Its stated range starts at 0.495 µm. Below 495 nm this is an **extrapolation** towards the absorption edge (~2.7 eV), not a measurement | estimate |
| Gelatin, dry film | ≈ 1.53–1.54 (visible) | — | Commonly quoted value in emulsion optics. Confirm with a measured gelatin film dispersion (e.g. prism coupler data). [spot-check] **Unverified.** No primary source was opened. A forum post attributes "about 1.54" for dry pure gelatin to H. Bjelkhagen, *Silver-Halide Recording Materials* (Springer); opening that book would confirm it. Dichromated gelatin (1.58 at 633 nm, NASA report) is not pure gelatin | recalled |
| Gelatin, swollen in processing | falls towards water (1.33) with swelling | — | physical reasoning; matters only if exposure through wet layers is ever modelled (it is not) | estimate |
| Cellulose triacetate (TAC) | ≈ 1.47–1.48 | 589 nm | Handbook value. Not opened. The refractiveindex.info entry "cellulose" (Sultanova, Kasarova & Nikolov, *Acta Phys. Pol. A* **116**, 585 (2009)) gives 1.4701 at 589 nm, 1.4752 at 500 nm and 1.4663 at 700 nm. That entry is for a "cellulose" film whose exact material (cellulose vs an acetate) the database does not state. [spot-check] US 6,398,371 B1 (Fuji Photo Film Co.) states: "the refractive index of the triacetyl cellulose used as the transparent support is 1.48". No wavelength is given; read on Google Patents | 1.48 confirmed (patent); Sultanova numbers opened (data file) |
| PET, isotropic sample | 1.5811 (500), 1.5733 (550), 1.5690 (590), 1.5637 (650), 1.5605 (700 nm); k ≈ 1×10⁻⁶ | — | X. Zhang, J. Qiu, X. Li, J. Zhao, L. Liu, *Appl. Opt.* **59**, 2337–2344 (2020), doi:10.1364/AO.383831 (tabulated in the database file `Zhang.yml`) | opened-source (data file) |
| PET, biaxially oriented film (Estar is this kind of film) | n_x ≈ 1.66–1.70, n_y ≈ 1.54–1.65, n_z ≈ 1.54–1.55 at 633 nm, depending on stretch ratio (lower n_z near 1.49–1.50 for strongly biaxial film: recalled) | — | Values quoted in optical-film patents found by search (USPTO hits such as US 10,350,818). The patents were not opened, so which patent gave which value is unverified. KLA/Filmetrics lists 1.6357 at 632.8 nm for "PET, Estar, Melinex, Mylar". Not checked against an Estar-specific measurement. Not needed for D1: the 135 and 120 bases are acetate; Estar only matters for sheet film | secondary-source |

Notes:

- The AgBr/gelatin index ratio is m ≈ 1.45–1.50 across 500–650 nm. That puts sub-micron grains in the strong-Mie regime (§2.3).
- Gelatin (1.54) over TAC (1.475) gives normal-incidence reflectance 4.6 × 10⁻⁴ (Fresnel, computed). The emulsion–base interface is optically negligible except at grazing angles. Total internal reflection gelatin→TAC needs θ > 73.3°.
- Only the **base–air** interface matters for halation: R = 3.7 % at normal incidence, 5.0 % at 30°, 18.6 % at 40°, 41 % at 42°, and 1 above θc = 42.7° (TAC, n = 1.475; Fresnel, unpolarised, computed).

### 1.2 Base thickness (Kodak publications)

| Film / format | Base | Source | Tag |
|---|---|---|---|
| T-Max 100 (TMX), 135-24/36, 35 mm × 100 ft | "5-mil (0.13 mm) acetate" | Kodak F-32, *KODAK T-MAX Professional Films* (2002), p. 26, "Sizes available" | opened-source |
| T-Max 100, 120 | "4.7-mil (0.11 mm) acetate". 4.7 mil = 119 µm; Kodak's 0.11 mm is a rounding down. Use 119 µm | F-32 (2002) p. 26; F-4016 (May 2007) p. 2: "120-size film coated on a thicker (4.7-mil) base than other black-and-white roll films" | opened-source |
| T-Max 100, sheets | "7-mil (0.18 mm) ESTAR Thick Base" | F-32 (2002) p. 26 | opened-source |
| Tri-X Pan (TX, the predecessor of 400TX), 135, 35 mm long rolls, 70 mm | "5-mil acetate" | Kodak F-9 (April 2003), p. 1, sizes table | opened-source |
| Tri-X Pan (TX), 120 | "3.6 mil acetate" (91 µm) | F-9 (2003) p. 1 | opened-source |
| Tri-X 320 (320TXP), 120/220 | "3.9-mil acetate" (99 µm) | F-4017 (May 2007) p. 1 | opened-source |
| Tri-X 320, sheets | "7-mil KODAK ESTAR Thick Base" | F-4017 (May 2007, Dec 2016) p. 1 | opened-source |
| Tri-X 400 (400TX), current 135 and 120 | **not stated** in any F-4017 edition read | — | not found |

1 mil = 25.4 µm: 3.6 mil = 91.4 µm, 3.9 mil = 99.1 µm, 4.7 mil = 119.4 µm, 5 mil = 127 µm, 7 mil = 177.8 µm.

### 1.3 Grey base and antihalation

- Kodak H-24 *Processing KODAK Motion Picture Films*, Module 15 "Processing Black-and-White Films", p. 15-5: "All black-and-white Eastman cine camera films, and some laboratory films, have a gray antihalation safety base. The gray base reduces light piping during exposure, minimizing the halo effect and producing a sharper image." Fig. 15-1 shows the stack OVERCOAT / SENSITIZED LAYER / SUPPORT and is marked "not drawn to scale". It gives **no thicknesses**. (opened-source)
- Kodak EASTMAN DOUBLE-X 5222 technical data (2022 PDF): "This film has a grey acetate safety base." (opened-source; motion-picture stock, cited only as evidence of Kodak practice.)
- **Still-camera Tri-X / T-Max datasheets say nothing** about a grey base, antihalation underlayer or backing. Grep of all F-4016, F-4017, F-32 and F-9 texts for "antihalation / backing / undercoat" returned nothing. (not found)
- Indirect evidence, from the sibling report `datasheets.md`. Lowest plotted density (≈ base + fog) is 0.27–0.32 for Tri-X 35 mm in D-76 against 0.16–0.21 for 120. In T-MAX developer it is 0.31–0.38 against 0.20–0.27. So 35 mm runs about **0.10–0.11 higher than 120** in both developers. The emulsion is nominally the same, and fog would not be expected to differ by format, so the **grey-dye density of the 35 mm acetate is about 0.10 (diffuse visual)**. This is an estimate. Part of the difference could be a format fog or base difference. It is lower than the 0.2–0.3 figure in the task brief; the brief's figure looks like total base + fog. **Use D_grey ≈ 0.10 as the prior and tag it `estimate`.** [spot-check] I re-read the lowest points of the four D-76 curves per format from the F-4017 (Oct 2021) p. 9 vectors, using the density-label positions; the gridlines were not fitted. They give ≈ 0.28–0.34 (35 mm) and ≈ 0.18–0.22 (120), so the ≈ 0.1 offset reproduces. Its attribution to grey dye remains an unverified interpretation.
- Why a grey base works (estimate, geometry). Image-forming light passes the base only after exposure, so it is never attenuated. A halation ray crosses the base twice at angle θ, giving optical path density D_grey·2/cos θ. At θc = 42.7° (cos θc = 0.735) that is 2.72·D_grey. With D_grey = 0.10 the ring is attenuated to 10^(−0.272) = 0.53, and further out (θ = 60°) to 10^(−0.40) = 0.40. Light piping along the film plane (θ → 90°) is suppressed almost completely, which is the stated purpose.
- 120 roll film. Classic practice is a removable, dyed gelatin anti-curl/antihalation backing on the back of the base, decolourised in processing; T-Max is also said to carry dyes in the emulsion (the "pink" T-Max stain). **Kodak's still-film datasheets do not document either** (recalled; would be confirmed by a Kodak patent for the TMX/400TX layer structure or a cross-section image). One related datasheet sentence: F-4016 (2018) and F-32 (2002) say that a "magenta (pink) stain after fixing" means the fixer may be near exhaustion or the fixing time too short. This shows that dyes which should be removed in processing are present in T-Max. (opened-source)
- Recommendation: model 135 as clear-index TAC + grey dye (D ≈ 0.10, spectrally neutral, `estimate`); model 120 as TAC + an index-matched backing with an absorbing dye whose density is `assumed` (fitted only if a halation measurement ever appears). Without data, the halation strength is not a validated quantity in D1.

### 1.4 Emulsion and overcoat thickness

| Quantity | Value | Tag |
|---|---|---|
| Tri-X / T-Max emulsion layer thickness | **not found** (no Kodak figure) | — |
| Dry layer from coverage | t ≈ m_Ag·(187.77/107.87)/ρ_AgBr + m_gel/ρ_gel; 1.75 g/m² Ag + 4.84 g/m² gelatin gives 4.1 µm (sibling report `latent-image.md` §2.3, using a Kodak tabular-grain patent coating) | estimate |
| Overcoat (protective gelatin) | ≈ 1 µm is typical | recalled |
| Number of emulsion layers | T-Max 400 is said to have a two-layer structure (forum description); T-Max 100 and 400TX unknown | not found |

---

## 2. Scattering in the emulsion and the MTF

### 2.1 Models in the literature

| Work | What it gives | Citation | Tag |
|---|---|---|---|
| Frieser's exponential spread function | Line spread function L(x) ∝ exp(−|x|/a), so MTF(ν) = 1/(1 + (2πaν)²). Frieser expressed the constant as a "spread" length; for typical camera negatives a is a few µm. Frieser's own constants **were not found** | H. Frieser, "Spread function and contrast transfer function of photographic layers", *Photogr. Sci. Eng.* **4**, (1960) (volume/year confirmed from a reference list; page not confirmed). [spot-check] Confirmed in a peer-reviewed source. P. Nowak, "Modelling of the modulation transfer function of silver halide light-sensitive layer including DIR and DAR dye couplers", *Optica Applicata* **28**(2), 95 (1998), Eq. 3–5, gives the Frieser LSF as exponential in |x| with "Frieser constant" K and MTF(ν) = 1/[1 + (Kπν/2.303)²]. That equals the form here with a = K/(2 ln 10). Nowak cites Frieser's book *Photographische Informationsaufzeichnung* (Focal Press, 1975) and Frieser, *Kinotechnik* 17, 167 (1935) | functional form: confirmed (Nowak 1998); 1960 citation: secondary-source |
| Gilmore 1967 | Point spread function from a simplified diffusion (radiative transfer) calculation in a thin emulsion; PSF singular at the origin; limits on PSF and MTF | H. F. Gilmore, "Models of the point spread function of photographic emulsions based on a simplified diffusion calculation", *JOSA* **57**, 75 (1967), doi:10.1364/JOSA.57.000075 | secondary-source |
| Wolfe, Marchand & DePalma 1968 | Volume reflectance and transmittance of fine-grain AgBr layers, analysed with Kubelka–Munk two-flux theory, give an MTF. The calculated MTF is somewhat higher than the photographically measured one, partly because of halation from base–air reflections | R. N. Wolfe, E. W. Marchand, J. J. DePalma, *JOSA* **58**(9), 1245–1256 (1968), doi:10.1364/JOSA.58.001245 | secondary-source (abstract via library record) |
| DePalma & Gasper 1972 | Monte Carlo photon transport in emulsions, single-particle (Mie) phase functions, giving PSF/MTF versus grain size and silver loading | J. J. DePalma, J. Gasper, "Determining the optical properties of photographic emulsions by the Monte Carlo method", *Photogr. Sci. Eng.* **16**(3), 181 (1972) | secondary-source (from reference lists only) |
| Radiative-transfer MTF expression | MTF as a function of emulsion thickness, crystal size, concentration and refractive index | *Appl. Opt.* **29**, 151 (1990), doi:10.1364/AO.29.000151 (seen only as a search hit) | secondary-source |
| Kelly three-stage model | Optical spread (linear in exposure) → nonlinear D–log H → linear chemical spread (adjacency). This is the standard framework that explains MTF > 1 | D. H. Kelly, "Systems analysis of the photographic process. I. A three-stage model", *JOSA* **50**(3), 269 (1960), doi:10.1364/JOSA.50.000269; "II. Transfer function measurements", *JOSA* **51**(3), 319 (1961), doi:10.1364/JOSA.51.000319. [spot-check] The Optica abstracts were read. Part I describes a nonlinear, isotropic model that predicts "adjacency" effects. Part II reports measured optical and chemical transfer functions for several emulsions and processes, supporting the Part I model. The exact stage equations were not seen | secondary-source (Crossref + abstracts) |
| Nelson's adjacency model | Prediction of densities in fine detail using a chemical spread function | C. N. Nelson, *Photogr. Sci. Eng.* **15**, 82 (1971). [spot-check] Cited as the source of the "chemical scattering function" in Nowak, *Optica Applicata* 28, 95 (1998), ref. [13] | secondary-source |

### 2.2 How optical spread and chemical adjacency combine

- Exposure in the layer is H_eff(x) = (H ∗ PSF_opt)(x). This step is linear.
- Development is local and nonlinear: D = f(log H_eff). Developer exhaustion and bromide/iodide release by development diffuse laterally over tens of µm. A dense area depletes developer and releases inhibitor, which slows development nearby, and the reverse happens next to a thin area. In the simulator this comes out of the reaction–diffusion solver; nothing extra needs adding.
- Kelly's linearised view: D ≈ γ·log(H ∗ PSF_opt) − k·(D ∗ CSF), with CSF the chemical spread function (normalised). In the frequency domain the chemical term multiplies the response by 1/(1 − k·(1 − CSF̂(ν))) or similar, up to normalisation. At ν → 0 the chemical term is a uniform suppression that is already in the large-area D–log H curve, so the response is 1. At intermediate ν (spatial periods of 20–200 µm, comparable to the developer diffusion length in the layer and boundary layer), small areas escape the suppression and the response rises **above 1**. At high ν the optical PSF dominates and the response falls.
- Measured MTFs above 100 % at low frequency are therefore a development (adjacency) effect, not an optical one. Their size depends on agitation, dilution and developer. That is why the datasheets state "Small Tank" or "Large Tank" and the developer.
- **Kodak's MTF definition.** The curves are labelled "Response (%)". The standard Kodak method (recalled; confirm against a Kodak MTF method note or ISO 6328 / ANSI PH2.39 text) is to expose a sinusoidal pattern, measure the density modulation with a microdensitometer, convert the densities back to effective exposure through the large-area D–log H curve, and divide the effective-exposure modulation by the input exposure modulation. The **mean density / exposure level of the measurement is not stated** on either datasheet (not found).
- Consequence (already in the design doc): fitting PSF_opt directly to the datasheet curve double-counts adjacency. The calibration must run the virtual sine-target exposure → full development → microdensitometer → same D–log H inversion, and compare the resulting response with the Kodak curve.

### 2.3 Mie regime for AgBr microcrystals (computed)

A Bohren–Huffman Mie calculation for AgBr spheres (Schröter index) in gelatin (n = 1.54), with the sphere assumption and **independent scattering**, at AgBr volume fraction f = 0.115 (from the 4.1 µm / 0.47 µm AgBr example above). μ_s = 1.5·f·Q_sca/d; ℓ_s = 1/μ_s; transport mean free path ℓ* = ℓ_s/(1 − g). Tag: **estimate**.

| λ (nm) | d (µm) | x = πd·n/λ | Q_sca | g | ℓ_s (µm) | ℓ* (µm) |
|---|---|---|---|---|---|---|
| 550 | 0.2 | 1.76 | 1.12 | 0.62 | 1.0 | 2.7 |
| 550 | 0.4 | 3.52 | 3.89 | 0.76 | 0.60 | 2.5 |
| 550 | 0.6 | 5.28 | 3.87 | 0.74 | 0.90 | 3.4 |
| 550 | 1.0 | 8.80 | 1.81 | 0.56 | 3.2 | 7.3 |
| 550 | 1.5 | 13.2 | 2.35 | 0.75 | 3.7 | 14.7 |
| 500 | 0.4 | 3.87 | 4.08 | 0.76 | 0.57 | 2.4 |
| 650 | 0.4 | 2.98 | 3.01 | 0.75 | 0.77 | 3.1 |
| 650 | 1.0 | 7.44 | 2.10 | 0.52 | 2.8 | 5.8 |

Reading: for grains of 0.3–1 µm the scattering mean free path is 0.6–3 µm, below a 4–10 µm layer thickness. The emulsion is several scattering lengths thick, strongly forward scattering (g ≈ 0.6–0.8), and in the multiple-scattering regime: a diffusion-like PSF with a sharp core and exponential-like wings, consistent with Frieser's empirical form. Caveats:

1. At f ≈ 0.1 the grains are close-packed enough for dependent scattering, which lowers μ_s.
2. Monodisperse results oscillate; polydisperse populations smooth them.
3. Tabular T-grains (ECD 1–2 µm, thickness ~0.1 µm, lying flat) are strongly anisotropic scatterers. Mie spheres do not describe them, and their lower forward-scattering into large angles is the stated sharpness advantage. A D1 micro-mode Monte Carlo could use a T-matrix or simply calibrate a two-parameter PSF.
4. Absorption by AgBr (blue) and by sensitising dye changes the depth profile and also the spread.

Practical model for D1 (recommendation): per layer and per wavelength band, a radially symmetric PSF = (1 − w)·δ + w·(exponential kernel with length a_λ), two parameters per film, `fitted` through the full forward model. Priors from §2.4 (a ~1–2 µm T-Max, ~3–5 µm Tri-X), with a_λ growing towards blue for conventional grains.

### 2.4 Kodak MTF curves, digitised

Method: both PDFs carry the MTF curve as a vector path. I took the path points (Bezier segments sampled at 6 points each) and mapped them with a **piecewise log interpolation** through the printed gridlines (the printed decades are not exactly equal in length; the mismatch is up to ~1.5 pt, about 0.02 in log ν). Script logic: locate the chart frame and gridlines in `page.get_drawings()`, pick the path with line width 0.72 pt inside the frame, interpolate R at chosen ν. Estimated reading error ±1.5 % response, ±3 % frequency.

**T-Max 100**: Kodak F-4016 (February 2016 edition; the same figure, code F002_0542AC, is in the 2018 edition), p. 8, "Modulation Transfer Curves". Exposure: tungsten. Process: KODAK Developer D-76, 68 °F (20 °C), small tank. Densitometry: diffuse visual. (opened-source)

| ν (c/mm) | 2.5 | 3 | 5 | 7 | 10 | 15 | 20 | 30 | 40 | 50 | 60 | 70 | 80 | 100 | 120 | 150 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Response (%) | 104.0 | 105.4 | 108.9 | 110.7 | 112.2 | 113.5 | 113.8 | 111.6 | 106.2 | 99.6 | 92.2 | 84.7 | 78.0 | 64.5 | 52.5 | 39.3 |

The curve spans 2.44–152 c/mm. 50 % response at 126 c/mm.

[spot-check] I re-digitised the curve with every printed gridline and minor tick, using the 10–100 % decade length (66.9 pt) above 100 %. Results: 103.7 (2.5), 108.4 (5), 111.6 (10), 113.1 (20), 99.6 (50), 64.4 (100), 52.6 (120), 38.9 (150); peak 113.1 % at 19 c/mm; 50 % at 125 c/mm. That agrees with the table within ≈ 1 point. The peak is ≈ 113 %, not 114 %.

**Tri-X 400**: Kodak F-4017 (October 2021 edition), p. 8, "Modulation Transfer Function" (p. 7 in the December 2016 edition, p. 10 in the May 2007 edition). Process: KODAK PROFESSIONAL Developer D-76, large tank. Densitometry: diffuse visual. No exposure source stated (the 2007 edition prints "Exposure: X"). The 2007 edition notes the MTF "was generated using the older version of this film (KODAK TRI-X Pan Film 5063)". (opened-source)

| ν (c/mm) | 2.5 | 3 | 5 | 7 | 10 | 15 | 20 | 30 | 40 | 50 | 60 | 70 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Response (%) | 106.4 | 107.1 | 109.4 | 110.8 | 111.1 | 108.0 | 103.8 | 90.6 | 73.1 | 55.3 | 38.8 | 26.3 |

The curve spans 2.45–72 c/mm. 50 % response at 53 c/mm.

[spot-check] **Tail corrected to a range.** I re-digitised the curve with Kodak's printed minor ticks: x for 2–9, 20–90 and 200 c/mm; y for 2–9, 20–90 % and the 100 % gridline. Results: 106.1 (2.5), 109.0 (5), 110.6 (10), 103.6 (20), 90.7 (30), 73.9 (40), 55.0 (50), **36.4 (60)**; 100 % at 23.7 c/mm; 50 % at 52.5 c/mm. The path's last point (x = 207.41 pt) sits on the printed 70 c/mm gridline (x = 207.49 pt), at ≈ **24 %**. So the curve ends at ≈ 70 c/mm, not 72. Interpolating between major gridlines only gives the table's values (38.8 % at 60, 26.3 % at 70, end 72 c/mm). Kodak's minor ticks sit ≈ 1 pt off exact log spacing, so neither mapping is authoritative. Use 60 c/mm = 36–39 %, 70 c/mm = 24–26 %, with reading error ±3 points beyond 50 c/mm. Values up to 50 c/mm agree within ≈ 1 point.

**Descriptive tail fits** (estimate; for orientation only, **must not seed the optical PSF directly**). R(ν) = g/(1 + (2πaν)²) fitted to ν ≥ 60 c/mm (T-Max): g = 1.24, a = 1.54 µm, rms residual 0.004. For Tri-X (ν ≥ 30) the Lorentzian fits poorly (g = 1.65, a = 4.7 µm, rms 0.034); a Gaussian g·exp(−2π²σ²ν²) fits better (g = 1.20, σ = 3.97 µm, rms 0.002). The gains g > 1 are the adjacency boost folded into the tail, and a/σ mix optical and chemical spread, so they are upper bounds on the purely optical spread length. Usable as prior ranges: optical exponential length ≲ 1.5 µm (T-Max), ≲ 4–5 µm (Tri-X). [spot-check] I refitted to the re-digitised points. T-Max Lorentzian (ν ≥ 60): g = 1.22, a = 1.52 µm, rms 0.005, which reproduces. Tri-X (30–70 c/mm, tick mapping): Lorentzian g = 1.76, a = 4.97 µm, rms 0.044; Gaussian g = 1.23, σ = 4.09 µm, rms 0.010. The Tri-X fit constants carry the tail uncertainty above: the folded-in gain spans ≈ 1.2–1.8, not 1.2–1.65. The prior ranges stand.

Resolving power (T-Max 100, F-4016 p. 8): 63 lines/mm at TOC 1.6:1, 200 lines/mm at TOC 1000:1, "method similar to … ISO 6328". (opened-source) Tri-X: not given in the current F-4017.

---

## 3. Halation

Geometry (standard ray optics; estimate where numbers are computed):

- Light scattered downward in the emulsion enters the base. By Snell invariance n·sin θ is conserved through parallel layers, so total internal reflection at the base–air face happens when n_emulsion·sin θ_em ≥ 1, i.e. n_base·sin θ_b ≥ 1, independent of the intermediate indices. Inside the base the critical angle is θc = arcsin(1/n_base): 42.7° (TAC, n = 1.475), 39.6° (isotropic PET 1.57), 37.3° (PET in-plane 1.65).
- A ray leaving the emulsion at base angle θ returns to the emulsion at lateral distance r = 2t·tan θ. All rays with θ ≥ θc return with R = 1, so the halation ring has a sharp inner edge at r_c = 2t·tan θc. Inside it, there is a weak disc from partial Fresnel reflection (3.7–19 % for θ = 0–40°).
- Ring inner radius r_c (computed): TAC 135 (t = 127 µm): **234 µm**; Tri-X 120 (91 µm): 168 µm; T-Max 120 (119 µm): 220 µm; Estar sheets (178 µm, n 1.57–1.65): 271–294 µm. Add the emulsion contribution (≈ 2·t_em·tan θ_em, a few µm), which is negligible.
- Light fraction (estimate). For a Lambertian radiance distribution in gelatin (n = 1.54), the flux fraction with sin θ_em > 1/n_em is 1 − 1/n_em² = **0.58**. So more than half of the diffusely scattered light that reaches the back surface is trapped by TIR, if nothing absorbs it.
- Usable kernel. Write I(θ) for the radiant intensity per steradian in the base of the downward light leaving a point at the emulsion bottom (from the layer transport). The returned irradiance at radius r is

  E(r) ∝ I(θ)·R(θ)·10^(−2·D_grey/cos θ)·cos³θ / (2t)², with θ = arctan(r/2t),

  where R(θ) is the unpolarised Fresnel reflectance at base–air (1 for θ ≥ θc), D_grey is the base dye density at normal incidence, and the cos³θ/(2t)² factor is the Jacobian of r = 2t·tan θ for a ring of radius r. For a Lambertian source, I(θ) ∝ cos θ, so the kernel goes as cos⁴θ·R(θ)·10^(−2D_grey/cos θ). Higher orders land at r_k = 2k·t·tan θ with weight [R·10^(−2D_grey/cos θ)·ρ_em]^k, where ρ_em is the emulsion's diffuse reflectance back into the base. They are usually negligible with a grey base or backing.
- With D_grey = 0.10 the ring is attenuated to about 53 % at θc; with an absorbing backing it effectively disappears.
- Colour films: the red-sensitive layer is nearest the base, so halation is red-tinted. That is a D2 issue.

Published quantitative treatments:

- K. Vendrovsky, I. Pakoushko, "Halation by reflexion in photographic layers", *J. Photogr. Sci.* **12**, 71–75 (1964), doi:10.1080/00223638.1964.11737228 (secondary-source, Crossref; not opened). This is the obvious paper to open for measured ring profiles.
- Wolfe, Marchand & DePalma 1968 attribute part of the difference between computed and measured MTF to halation from base–air reflections (secondary-source).
- Kodak H-24 Module 15 (opened) gives the qualitative grey-base statement quoted in §1.3.

---

## 4. Density from developed silver

### 4.1 Nutting's formula

P. G. Nutting, "On the absorption of light in heterogeneous media", *Phil. Mag.* (6) **26**(153), 423–426 (1913), doi:10.1080/14786441308634988 (secondary-source, Crossref). For N opaque particles per unit area, each of projected area a, randomly placed (a Poisson / Boolean layer), the transmittance is T = exp(−N·a), so

D = log10(e)·N·a = 0.4343·N·a.

[spot-check] Nutting's paper itself was not opened; the Crossref metadata match. The Boolean coverage formula P = 1 − exp(−λE[πr²]), which gives the same T, is confirmed in Newson et al. IPOL 2017 Eq. 1. The overlap of particles is included exactly by the exponential (this is the Boolean-model coverage formula, §6). Summed over size classes: D = 0.4343·Σ_c N_c·a_c. This is the same equation as the Boolean model's mean, so the statistical mode and the micro mode can share size classes and covering power directly.

### 4.2 Covering power

- Definition (opened-source): "Covering power is commonly defined as the maximum image density divided by developed silver per unit area, typically reported in units of g/m² or mg/dm²" (EP 0 513 726 A1, applicant Eastman Kodak Company, inventors Loblaw, Tsaur, Kam-Ng; filed 12.05.92, published 19.11.92; quote on p. 2 of the EPO publication PDF). More generally, CP(D) = D/M_Ag at any density.
- From Nutting: CP = 0.4343·a/m_grain = 0.4343/(ρ_Ag·t_eff), with ρ_Ag = 10.49 g/cm³ and t_eff = V_Ag/a (silver volume per projected area). Computed (estimate):

| t_eff (µm) | 0.05 | 0.07 | 0.10 | 0.20 | 0.33 | 0.67 |
|---|---|---|---|---|---|---|
| CP (density per g/m² Ag) | 0.83 | 0.59 | 0.41 | 0.21 | 0.13 | 0.06 |

A compact sphere of d = 1 µm has t_eff = 2d/3 = 0.67 µm, giving CP ≈ 0.06. Filamentary silver from chemical development spreads the same mass over a much larger projected area (t_eff ~0.05–0.1 µm), giving CP ≈ 0.4–0.8 per g/m². This tenfold difference is why the morphology of developed silver, and not only its mass, sets density.
- Literature (titles and metadata checked through Crossref; texts not opened, so no numbers taken):
  - D. M. Farnell & L. R. Solman, "The relationship between covering power at saturation density and undeveloped grain size", *J. Photogr. Sci.* **11**, 347–350 (1963), doi:10.1080/00223638.1963.11736943. CP decreases with increasing grain size (consistent with an RIT thesis finding grain size the only significant factor among twelve).
  - G. C. Farnell & L. R. Solman, "The covering power of photographic silver deposits I. Chemical development", *J. Photogr. Sci.* **18**, 94–101 (1970), doi:10.1080/00223638.1970.11737534; L. R. Solman, "… II. Physical development", **18**, 136–147 (1970); "… III. Comparison of chemical and physical development", **18**, 179–184 (1970).
  - L. R. Solman & G. C. Farnell, "The relation between covering power and density for model emulsion systems", *J. Photogr. Sci.* **21**, 78–84 (1973), doi:10.1080/00223638.1973.11737716. CP is not constant along the curve.
  - A Kodak radiographic patent found by search claims D_max ≥ 3.0 at < 12 mg/dm² (1.2 g/m²) silver for tabular grains. This is double-coated radiographic film, so it is not directly comparable (secondary-source).
  - T. H. James (ed.), *The Theory of the Photographic Process*, 4th ed. (1977), the chapter on optical properties and covering power (ch. 18 or 21 depending on edition; I could not confirm which) is the general reference (recalled).
- Dependences to carry in the model (recalled qualitative facts; the Farnell–Solman titles support 1 and 4):
  1. CP falls with grain size.
  2. Solvent developers (high sulphite, e.g. D-76) produce more compact silver (lower CP); non-solvent, high-activity developers produce more filamentary silver.
  3. Physical development (silver deposited from solution) gives compact silver with lower CP.
  4. CP varies with the degree of development and with density.
  5. Tabular grains have high CP because the developed silver keeps the large projected area.
- For D1: CP per size class and per film is a `fitted` parameter, with priors 0.3–1.0 density per g/m² for conventional chemical development. Silver coverage of Tri-X / T-Max is **not published** (not found), so CP and coated silver are degenerate unless D_max at long development pins them.

### 4.3 Diffuse versus specular density; Callier Q

- Q = D_specular/D_diffuse ≥ 1. Silver grains scatter as well as absorb, so a narrow-acceptance (specular) system loses scattered light and reads higher density.
- Typical values: Q ≈ 1.3–1.5 for coarse-grained developed negatives and ≈ 1.0–1.1 for fine-grained images (recalled; confirm in James ch. on Callier effect, or Streiffert, J. SMPTE 1947, on Q of motion-picture emulsions).
- A. E. Saunders, "A discussion of the factors influencing Callier coefficient", *J. Photogr. Sci.* **34**, 26–32 (1986), doi:10.1080/00223638.1986.11738421: Q varies with diffuse density; base density and surface scatter (matting agent, reticulation) change the shape of the Q–D curve (secondary-source, abstract seen in search).
- For the simulator: Nutting's opaque-particle formula gives the absorbed or blocked fraction, which is closest to **diffuse** density. The datasheet curves are all "Diffuse Visual", so calibrate against diffuse density. The "standard scan" exit is a near-specular condenser or LED scanner and should apply Q(D) > 1, as an `assumed` function, documented, since no per-film Q is published.

### 4.4 What "diffuse visual density" means

- Spectral: ISO 5-3:2009, §6.1, "ISO 5 standard visual density", notation D_T(S_H:s_V). Read from the sample preview of ISO 5-3:2009 hosted by iTeh Standards (a reseller, cdn.standards.iteh.ai), not from iso.org; cite the standard itself. The spectral products "are chosen to match the product of the spectral luminous efficiency function for photopic vision, V_λ … and the relative spectral power distribution of the influx spectrum specified for reflection measurements, S_A. This is essentially the CIE tristimulus Y function of illuminant A" (opened-source, ISO 5-3:2009 p. 7, via the iTeh preview). The influx is CIE illuminant A (Planckian, ≈ 2856 K); for transmittance, S_H is illuminant A modified by a heat-absorbing filter, in the region above 560 nm (§4.2.3; [spot-check] the earlier "§4.2.1" was wrong). Weighting factors for computation from spectral data at 10/20 nm intervals are in ISO 5-3 Table 8 and the files "10nmWeights.csv" etc. **The numerical table was not opened** (behind the paywall).
- Geometry: ISO 5-2:2009 specifies "ISO 5 standard diffuse" and f/4.5 and f/1.6 projection transmittance densities. The 1985 edition replaced the integrating sphere with an opal-type diffuser (secondary-source, ANSI preview description). The diffuse condition is, as far as I recall, a fully diffuse influx (opal) with collection in a narrow cone around the normal (half-angle about 10°), or the reciprocal arrangement. Recalled; confirm in ISO 5-2 §6.
- For the simulator's densitometer: compute T(λ) from the layer model, weight with V(λ)·S_A(λ) (or S_H; equivalent for non-fluorescent film), integrate, and use hemispherical (diffuse) transmittance, i.e. include all scattered transmitted light.

---

## 5. Granularity

### 5.1 Selwyn's law

E. W. H. Selwyn, "A theory of graininess", *Photogr. J.* **75**, 571 (1935) (recalled citation; not found in Crossref). For apertures much larger than the grain correlation length, σ_D·√A = constant. Selwyn granularity is S = σ_D·√(2A). With the Kodak aperture (A = π·24² = 1810 µm²): Tri-X S = 0.017·√3619 = **1.02 µm**; T-Max S = **0.48 µm** (estimate). The Zhang et al. patent (§6.2) uses the same definition: G_Selwyn = G_RMS·√(2A(S)) (opened-source).

### 5.2 Kodak RMS diffuse granularity

- Definition on the datasheets (opened-source): "Read at a net diffuse density of 1.00, using a 48-micrometre aperture, 12X magnification" (F-4016 June 2018 p. 8; F-4017 Oct 2021 p. 7, Dec 2016 p. 6, May 2007 p. 10). The number is σ_D × 1000. "Net" means density above base + fog. The 12× refers to the viewing magnification at which the 48 µm aperture correlates with graininess; it is not an instrument setting.
- Values:

| Film | RMS granularity | Developer | Note | Source |
|---|---|---|---|---|
| T-Max 100 | 8 | D-76, 20 °C | — | F-4016 (2018) p. 8 |
| Tri-X 400 | 17 ("fine") | HC-110 dilution B, 20 °C | "generated using the old versions of these films"; sibling report: 7½ min small tank, CI 0.56 | F-4017 p. 7 (Oct 2021), p. 6 (Dec 2016), p. 10 (May 2007) |
| Tri-X 320 | 16 | HC-110 (B) | same note | F-4017 |

- So the Tri-X granularity validation belongs to the **Tri-X × HC-110** training pair, not to D-76.
- Kodak does not state the microdensitometer's optical numerical aperture. "Diffuse" granularity presumably means the readings were made or converted to diffuse density (not found; recalled practice).

### 5.3 From a grain model to RMS granularity (Siedentopf)

For a Poisson number of opaque grains of area a in an aperture A, with density read through Nutting linearly in the count: σ_D² = 0.4343·D·a/A. H. Siedentopf, *Phys. Z.* **38**, 454 (1937) (recalled citation). For mixed sizes: σ_D² = 0.4343²·Σ_c N_c a_c²/A = 0.4343·D·⟨a²⟩/⟨a⟩/A. This is a weighted mean area, dominated by the largest grains.

Implied effective (developed-silver) areas at D = 1, A = 1810 µm² (estimate):

| Film | σ_D | a_eff (Siedentopf) | Equivalent disk radius | Implied N at D = 1 |
|---|---|---|---|---|
| Tri-X 400 | 0.017 | 1.20 µm² | 0.62 µm | 1.9 µm⁻² |
| T-Max 100 | 0.008 | 0.27 µm² | 0.29 µm | 8.6 µm⁻² |

A binomial count (n of N_total grains developed, per class) gives Var(n) = N·p(1 − p), smaller than Poisson by (1 − p). At D = 1 for the largest classes (p near 1) this matters, and it should be kept, since D1 samples binomially.

### 5.4 Correction for grain overlap (Boolean model) — important for D1

The transmittance of a Boolean layer of opaque disks is not a linear function of the count, because overlapping grains shadow each other. The exact covariance of the indicator (the same for coverage and transmittance) for constant mean coverage is

C(τ) = q²·(exp(λ·γ(τ)) − 1), q = exp(−λ·E[πr²]) = mean transmittance,

γ(τ) = E[ intersection area of two disks of radius r at distance τ ] = E[ 2r² arccos(τ/2r) − (τ/2)√(4r² − τ²) ] for τ < 2r.

This is the textbook Boolean-model result; it is printed as Eq. 4–5 in the Zhang et al. patent, opened-source. For a large aperture Var_T ≈ (1/A)∫C(τ)d²τ and σ_D ≈ 0.4343·σ_T/T̄. Expanding exp(λγ) − 1 ≈ λγ gives back Siedentopf. The exact value is larger, by κ(D) (computed by numerical integration, constant-radius disks; estimate):

| Net D | 0.1 | 0.3 | 0.5 | 1.0 | 1.5 | 2.0 |
|---|---|---|---|---|---|---|
| κ = Var_Boolean/Var_Siedentopf | 1.055 | 1.18 | 1.34 | 1.89 | 2.86 | 4.66 |

For lognormal radii with coefficient of variation 0.5, κ(1.0) = 1.58. B. E. Bayer, "Relation between granularity and density for a random-dot model", *JOSA* **54**(12), 1485 (1964), doi:10.1364/JOSA.54.001485, gives the analytic treatment of this effect (secondary-source).

Consequences:

- The effective disk area that reproduces Kodak's numbers under the exact Boolean model is about half the Siedentopf value: Tri-X **0.64 µm²** (r = 0.45 µm), T-Max **0.14 µm²** (r = 0.21 µm) for constant radii.
- A statistical mode that draws binomial counts and maps them **linearly** to density (Nutting per pixel) reproduces Siedentopf, so it will under-predict the micro mode's σ_D by √κ ≈ 1.37 at D = 1. Remedy: compute the per-pixel transmittance noise with the Boolean covariance. Either add independent extra variance (κ(D) − 1)·σ_Siedentopf² with the Boolean correlation shape, or synthesise the noise field with the covariance C(τ) (§6.2). The D1 check "statistical and micro modes agree on granularity" will fail otherwise.
- Real granularity also includes the clumping of developed silver, the microdensitometer optics, and the 3D layer. For opaque grains in transmission, overlap in projection is exactly what Nutting and the Boolean model count, whatever the depth of each grain. Depth enters instead through the diffuse-illumination geometry (oblique rays see different projections) and through light scattered by the silver, i.e. the Callier effect. The "a" fitted to granularity is an effective parameter (`fitted`), not the crystal's projected area.

### 5.5 Noise power spectrum shape

- For a stationary grain field, NPS(ν) = Fourier transform of the density autocovariance. NPS(0) = σ_D²·A for large A, i.e. NPS(0) = S²/2 (Selwyn).
- For a Boolean field of disks, the covariance has support 2r_max (plus the exp(λγ) nonlinearity), so the NPS is flat ("white") up to roughly 1/(2πr) and rolls off like the squared disk transform |2J₁(2πrν)/(2πrν)|². For r ≈ 0.3–0.6 µm that roll-off sits at hundreds of c/mm, far beyond the 48 µm aperture and any scanner.
- Measured film NPS also falls at tens of c/mm. This comes from (i) grain clumping and the depth distribution of grains, and (ii) the measuring microdensitometer's MTF², not from the exposure PSF (grain noise is not filtered by the exposing optics). Published film NPS curves (Dainty & Shaw, *Image Science*, 1974, ch. 8–9; recalled) look roughly flat to ~10–20 c/mm and then fall. No Kodak NPS curves exist for these films (not found).

---

## 6. Boolean grain model for rendering

### 6.1 Newson, Delon, Galerne (2017)

Citations (Crossref-checked):

- A. Newson, J. Delon, B. Galerne, "A Stochastic Film Grain Model for Resolution-Independent Rendering", *Computer Graphics Forum* **36**(8), 684–699 (2017), doi:10.1111/cgf.13159 (secondary-source; the paper itself was blocked).
- A. Newson, N. Faraj, J. Delon, B. Galerne, "Realistic Film Grain Rendering", *Image Processing On Line* **7**, 165–183 (2017), doi:10.5201/ipol.2017.192 (opened-source; this is the implementation paper, with C++ code on the IPOL page).
- Preprint: "Stochastic modeling and resolution-free rendering of film grain", MAP5 tech. report 2016-23, HAL hal-01358392.

Model definition (IPOL §3, opened-source):

1. Z = ∪ᵢ B(xᵢ, rᵢ): centres xᵢ form a Poisson process of intensity λ in the plane; radii rᵢ are i.i.d., independent of positions, either constant or log-normal with mean μ_r and standard deviation σ_r.
2. Coverage: P(1_Z(y) = 1) = 1 − exp(−λ·E[πr²]) (IPOL Eq. 1).
3. Inhomogeneous model: the input image u is normalised, ũ = u/(u_max + ε) ∈ [0, 1), with ε = 0.1 grey levels, and λ(y) = (1/E[πr²])·log(1/(1 − ũ(⌊y⌋))) (Eq. 2), piecewise constant per input pixel. E[πr²] = π(μ_r² + σ_r²). Grain units: one input pixel = one unit length.
4. Filtering: the observed image is v = φ ∗ 1_Z, with φ a Gaussian of standard deviation σ in output pixels (default **σ = 0.8 output pixels**). It is evaluated by Monte Carlo, v(y) = (1/N)·Σ_k 1_Z((y − ξ_k)/s) with ξ_k ~ N(0, σ²I) and s the zoom factor. The same offsets ξ_k are used for every pixel (Eq. 4).
5. Two algorithms:
   - Grain-wise: for each input pixel draw Q ~ Poisson(λ), place Q centres uniformly in the pixel, draw radii, rasterise into N binary images, average.
   - Pixel-wise: partition the plane into cells of size δ = 1/⌈1/μ_r⌉, regenerate the grains of any cell reproducibly from a cell-seeded PRNG, and for each output sample test cells within r_max. r_max = the 0.999 quantile of the radius distribution.
   - Complexities are given in Eq. 8–10; pixel-wise is faster for small σ_r/μ_r.
6. Reported parameter ranges: μ_r ∈ {0.025 … 0.3} input pixels, σ_r/μ_r from 0 to 0.9. Timing for a 512 × 512 constant grey (128): 40 s grain-wise and 2.4 s pixel-wise (parallel).
7. Stated limitation: it cannot reproduce emulsions that look like "white grains" (their Ilford Delta 3200 example), and it is sparse in dark areas.

Mapping onto the Darkroom Atlas micro mode (recommendation):

- In the photographic version the disks are the developed silver of each grain, and 1 − 1_Z is the transmittance. **The image is a negative**: λ comes from the simulated developed-grain population, not from an input grey level. Use physical units: µm, disk area = the class's effective developed-silver area a_c (§5.4), λ_c = N_c(x)·p_dev,c(x) per µm².
- Mean transmittance = exp(−Σ_c λ_c a_c), which is exactly Nutting. So per-class counts and covering powers carry over unchanged between the statistical and micro modes.
- The Gaussian filter φ plays the role of the microdensitometer aperture or the scanner and enlarger optics. For a granularity check use a 48 µm circular aperture instead.
- For a 3D micro mode, project each developed grain's silver onto the film plane. The union of the projected disks is again a Boolean set, provided positions are Poisson in 3D.

### 6.2 Zhang, Wang, Tian, Pappas (2023): analytic acceleration

- Citation (Crossref-checked): Kaixuan Zhang, Jingxian Wang, Daizong Tian, Thrasyvoulos N. Pappas, "Film Grain Rendering and Parameter Estimation", *ACM Transactions on Graphics* **42**(4), pp. 1–14 (2023), doi:10.1145/3592127, SIGGRAPH 2023. The article number was not visible in the metadata (not found). The paper itself was not opened (ACM blocked). The abstract claims a speed-up of six orders of magnitude over Newson's Monte Carlo and real-time rendering when individual grains are not visible (secondary-source).
- The same method is in **US 12,705,810 B2**, "Real time film grain rendering and parameter estimation", Northwestern University; inventors D. Tian, J. Wang, K. Zhang, T. N. Pappas; filed 18 Mar 2024, provisional 63/490,545 of 16 Mar 2023; granted 11 Aug 2026 (opened-source; USPTO image PDF, columns 5–12 read). The equations below are from the patent.
  - Eq. 3: λ(p) = log(1/(1 − u[p̄]))/(π·E[r²]), with u the normalised grey value (coverage).
  - Eq. 4–5: C_B(Δ) = Cov{1_B(p), 1_B(p + Δ)} = (1 − u)²·(exp(λ·γ(‖Δ‖)) − 1), with γ(‖Δ‖) = ∫_{‖Δ‖/2}^∞ (2x²·arccos(‖Δ‖/2x) − (‖Δ‖/2)·√(4x² − ‖Δ‖²))·f_R(x) dx.
  - Lemma 2.2: scale invariance, C_B(s‖Δ‖, s·r) = C_B(‖Δ‖, r).
  - Eq. 13: Var{B̃(p)} = (1/2σ²)·∫₀^∞ exp(−x²/4σ²)·C_B(x)·x dx, for B̃ = φ ∗ 1_B with a Gaussian φ of standard deviation σ.
  - Eq. 16: for σ ≳ 3μ_r (grains not individually visible) the correlation of B̃ is ≈ exp(−‖Δ‖²/4σ²), independent of grey level. [spot-check] The patent (col. 9, ll. 16–19) literally prints "T(p, ‖Δ‖)≈1 if σ≤3μ_r", allowing 3 % approximation error. That contradicts its own premise a few lines earlier ("σ_r=0 and 2μ_r<<σ") and its "individual grains are not visible" remark, so it looks like a typo. Read it as σ ≥ 3μ_r, and confirm against the ACM TOG paper.
  - Eq. 17–18: RMS granularity ≈ √(Var{B̃}·A(S))/((1 − E{B̃})·A(S))·log e; G_Selwyn = G_RMS·√(2A(S)).
  - Rendering, Eq. 19–24: y[p] = u[p] + (n ∗ h)[p], with n[p] ~ N(0, Var{B̃[p]}) signal-dependent white noise and h an isotropic filter normalised so that Σh² = 1. h is chosen by least squares so that (h ∗ h)[Δ]/Σh² ≈ exp(−‖Δ‖²/4σ²); a sampled Gaussian suffices when σ > 1 pixel. Clipping to [0, 1] is harmless in practice; a bounded distribution (truncated Gaussian, binomial, beta) is an alternative.
  - Parameter estimation, Eq. 25–37: assume σ_r = 0. Then Var{B̃} ≈ (μ_r²/2σ²)·∫₀² C_B(x; r = 1)·x dx (Eq. 32). σ is estimated from the residual correlation (Eq. 34), μ_r from the residual variance per denoised grey level (Eq. 36), and q-bit quantisation is modelled as uniform noise of variance 2^(−2q)/12 (Eq. 25–29, 37).
- Use in D1:
  1. Eq. 4, 5 and 13 give the **exact** aperture-filtered variance of a Boolean layer. They are what §5.4 needs to make the binomial statistical mode agree with the micro mode, with a Gaussian replaced by the 48 µm disk aperture for granularity: Var = (1/A²)∫C_B(τ)·K_A(τ)d²τ, with K_A the aperture autocorrelation (overlap area of two disks of radius 24 µm at distance τ).
  2. Zhang's correlated-Gaussian synthesis is a cheap full-frame renderer whose variance and correlation match the Boolean model when grains are below the display resolution. It is a candidate for the "statistical mode" output stage. It is not a physical model of development, so it should take its per-pixel mean from the simulated density and its variance from Eq. 13 / §5.4.
  3. The micro-mode reference renderer for a small patch is Newson's pixel-wise algorithm, in µm units, with a disk aperture or Gaussian φ.

---

## 7. Not found

- Tri-X 400 (current 400TX) base type and thickness for 135 and 120; any Kodak statement of grey base or antihalation layer on still B&W films.
- Emulsion and overcoat thicknesses, silver coverage and number of layers for Tri-X 400 and T-Max 100.
- Frieser's spread constants for named films, and the page number of Frieser 1960.
- Numerical values from DePalma & Gasper 1972, Gilmore 1967, Vendrovsky & Pakoushko 1964, Farnell & Solman 1963/1970/1973, Saunders 1986 (texts behind paywalls; only metadata and abstracts seen).
- Callier Q values specific to Tri-X or T-Max.
- ISO 5-2 exact angular tolerances and the ISO 5-3 Table 2/8 weights (paywalled).
- Kodak's MTF measurement method, mean density and exposure level for the published curves; the Tri-X MTF exposure source.
- Kodak's microdensitometer optics (numerical aperture) for the "diffuse" RMS granularity.
- Noise power spectra for Tri-X or T-Max.
- Zhang et al. 2023 article number and full text (the patent text was read instead).
- A gelatin dispersion curve from a measured source (the gelatin index is still unverified), and a TAC dispersion curve from an opened measured source. [spot-check] A single TAC value of 1.48, with no wavelength, is now sourced to US 6,398,371 B1.
