# D1 research: calibration, virtual sensitometry and validation methodology

Research date: 2026-11-21. Scope: what the virtual sensitometer and densitometer must implement, how to digitise the Kodak curves, how to fit a mechanistic model with roughly 10–20 film parameters and shared chemistry constants, which accuracy is realistic, and which open implementations exist.

Confidence labels (same convention as the sibling reports `datasheets.md` and `latent-image.md`):
* **opened**: I read the source myself (for ISO standards: the publicly hosted preview pages, see 1.0).
* **secondary**: from a search-engine summary, abstract page or index, not from the document itself.
* **estimate**: my own computation from opened data (scripts were run on files in `.cache/datasheets/`).
* **recalled**: from memory of the literature; not opened. The entry says what would confirm it.

Files fetched for this report and kept, git-ignored, under `.cache/datasheets/`: `iso6/ISO-6-1993-preview.pdf`, `iso5_2/ISO-5-2-2009-preview.pdf`, `iso5_2/ISO-5-2-1991-preview.pdf`, `iso5_3/ISO-5-3-1995-preview.pdf`, `iso7589/ISO-7589-1984-preview.pdf`, `smpte/rp0014-1997.pdf`, `nist/SP250-99.pdf`, `aapm/RPT_94.pdf`, `z133/z133.pdf` (Kodak AJ-3 + Z-133E bundle), `f4016/f4016_TMax_100_2018.pdf` (MD5 4a193fb16cbf43cd0e6d40d03c7ccd89). Added by the spot-check: `spotcheck_cv/h740.pdf` (Kodak H-740).

> **Independent spot-check, 2026-11-21.** 15 key facts were re-checked against their sources. See section 9 for the verdict table.
>
> Three corrections were made in place:
> * Z-133 Formula 2 is printed two ways by Kodak.
> * The AAPM §3.2.3 page number.
> * The Radiology 1983 citation and its context.
>
> Two upgrades were made:
> * The CI 0.2/2.2 construction is now confirmed in Kodak H-740.
> * The ISO Type 2 and lens-τ tables are legible in the cached previews.
>
> Edits are marked "*Spot-check*".

---

## 0. Answer in brief

1. **ISO 6 does not use visual density.** ISO 6:1993 §5.5 specifies "ISO standard Type 2 diffuse transmission density". In ISO 5-3, Type 2 is the printing density for silver-halide print material, a blue/UV-weighted response. It is not ISO visual density. The sibling report `datasheets.md` §5.2 says ISO 6 densitometry "is Kodak's diffuse visual"; **that statement is wrong** and should be corrected. Kodak's curves are diffuse **visual**. The virtual densitometer therefore needs both:
   * Type 2 for the ISO speed check.
   * Visual for comparison with the Kodak curves.
2. **The visual spectral product needs no ISO table.** ISO 5-3 §8.1 requires S_H·V_T = S_A·V, so the visual product equals CIE illuminant A × CIE V(λ). I checked this numerically against the ISO 5-3:1995 Table 2 values; it agrees within 0.002 in log₁₀ (section 1.4).
   * Two items have no public substitute: the **Type 2 spectral product** and the **ISO standard camera-lens transmittance** inside the ISO 7589 sensitometric daylight.
     * *Spot-check correction (2026-11-21):* both **can be read in full** in the freely hosted iTeh previews already cached. They are the ISO 5-3:1995 Table 2 Type 2 column, 340–540 nm (preview p. 4), and the ISO 7589:1984 Table 1 τ(λ), D₅₅ and S_λ columns, 350–690 nm (preview p. 4).
     * Availability is therefore not the problem. Whether the project may embed ISO's tabulated values in its data files is a licensing question for the user.
     * Caveat: ISO 6:1993 cites ISO 5-3:**1984**, and the anchors here come from the 1995 edition. Whether the Type 2 column changed between those editions is not known.
3. **The general contrast-index construction was not opened.**
   * The defining paper is Niederpruem, Nelson & Yule, Phot. Sci. Eng. 10(1) (1966) 35–41. The citation is confirmed through a patent; the paper itself was not opened.
   * Kodak's Z-133E gives two control-strip formulas, CI = (HD − LD)/2.26 + 0.10 and a five-step regression. *Spot-check:* Kodak prints the first one as (HD − TD) on p. 5 and on Worksheet 1 (p. 45), and as (HD − LD) on pp. 6, 8 and 9 (see 1.3.3). These are **regressions specific to the T-MAX 400 control strip**, not the definition.
   * *Spot-check (2026-11-21):* the 0.2/2.2 construction is **confirmed by a Kodak source**, *Basic Photographic Sensitometry Workbook*, Kodak publication H-740 (© 2006), p. 9. It uses a straight edge marked at 0.0, 0.2 and 2.2 log E units: the 0.0 mark sits on the D-min line, the 0.2 and 2.2 marks touch the curve, and CI is the slope of that line. H-740 is a Kodak teaching text, not the 1966 definition, and the 1966 paper is still unopened.
   * I tested the 0.2/2.2-arc construction and the 0.10 + 2.0 approximation on the vector T-Max 100 × D-76 curves (section 1.3.4).
     * The two constructions agree with each other to within 0.02.
     * Against Kodak's own CI(t) curve, only one time can be interpolated: at 7.5 min the gap is 0.03. At 6 min the agreement is better, and at 10 min the gap is about 0.04–0.05, but both of those are extrapolations.
     * A second check is consistent: the characteristic curves reach CI 0.56 at about 6.2 min, Kodak's CI curve at 6.30 min, and the time table gives 6½ min.
   * **About 0.03 CI is an upper bound from a single point**, not a measured property of the datasheet. It mixes the datasheet's own inconsistency with the unknown error of the (recalled) CI construction.
4. **Digitisation.**
   * Vector extraction from the 2016–2021 Kodak PDFs is exact with respect to the plotted line.
   * The remaining error lies in two other places:
     * how closely the plotted line follows Kodak's measurements, which is unknown (the CI curves have 3–8 vertices);
     * how representative production coatings are. Kodak's own NOTICE says the curves are "representative of production coatings" and are not specifications.
5. **Realistic accuracy scale, all from opened sources:**

   | Source | Scale |
   |---|---|
   | NIST density calibration | 0.007 D (k = 2) |
   | ISO 5-2 instrument realisation | < 0.01 D |
   | Densitometer repeatability | ±0.01 D |
   | AAPM repeat-exposure flags | 0.04 D (base + fog), 0.10 D (straight line) |
   | ISO 6 speed accuracy | 0.05 in log H_m, which is half of one 1/3-stop speed bin |
   | Kodak CI tuning tolerance | ±0.02 |
   | Kodak CI action limits | ±0.07 |
   | Kodak CI control limits | +0.20 / −0.12 |

   Section 4 turns these into proposed thresholds for D1.
6. **Fitting.**
   * Approach: weighted, prior-regularised least squares, using `scipy.optimize.least_squares(method="trf")` in log-parameter space with multi-start.
   * Add a global stage (DE or CMA-ES) only on a reduced, identifiable parameter set.
   * Check identifiability with the profile likelihood and the SVD of the weighted Jacobian.
   * Expect plateaus for the XTOL and HC-110 constants: they are constrained only by CI(t) curves and time tables, about 10–20 scalars per developer.
   * Freeze code and data hashes before the single held-out evaluation.
7. **Not found:**
   * a standard that specifies how reciprocity tests are run;
   * interlaboratory or batch data for pictorial black-and-white film;
   * the Status A/M tables (beyond the preview);
   * an open-source package that fits mechanistic parameters to H&D curves.

---

## 1. Standards and definitions for the virtual sensitometer and densitometer

### 1.0 Where the ISO text came from

* **ISO 6, 5-2, 5-3 and 7589 previews.** ISO sells these standards. I read only the **preview PDFs hosted by iTeh Standards** (cdn.standards.iteh.ai), a reseller, not iso.org. The ANSI webstore preview of ISO 7589:2002 returned HTTP 403.
* **Scope of the previews.** They cover the scope, definitions and, for most of these standards, the normative core, but not every table. Clause numbers below refer to the edition named.
* **Tables not reproduced.** ISO tables are copyright ISO, so this report quotes only a few anchor values. The implementation must take full tables from a licensed copy, or use the public CIE substitutes named below.

**Editions.** ISO 6:1993 (2nd edition) is current; ISO confirmed it in 2023 (secondary, iso.org listing). Its normative references are dated, and are older than the current editions:

| ISO 6:1993 cites | Current edition |
|---|---|
| ISO 5-2:1991 | ISO 5-2:2009 (5th edition, confirmed 2025) |
| ISO 5-3:1984 | ISO 5-3:2009 |
| ISO 7589:1984 | ISO 7589:2002 (2nd edition) |

Source: ISO 6:1993 p. 1 §2 (opened).

**Recommendation:** the simulator follows the editions ISO 6:1993 cites for the speed computation. Kodak's quoted speeds were determined under those editions, or earlier ones. Record the edition in the run manifest.

### 1.1 ISO 6:1993, black-and-white pictorial negative film speed

Opened: preview pp. 1–4 (https://cdn.standards.iteh.ai/samples/3580/192a0f97667044f3ba51276b0ded056a/ISO-6-1993.pdf). The sibling report §5.2 covers this too; the items below are re-verified, and several are added.

| Clause | Requirement | Simulator consequence |
|---|---|---|
| §3.1 | Exposure H is the time integral of illuminance on the film, in lux·s | Convert spectral irradiance × time to photometric H with 683 lm/W × ∫E_e(λ)V(λ)dλ × t, using the SPD **at the film plane**, after lens and modulator |
| §4 | Samples must represent several batches, aged as users would use them; "several independent evaluations" | The datasheet curve is a batch average; there is no single-batch truth to fit |
| §5.3.1 | Specimens held at 23 ± 2 °C and (50 ± 5) % RH during exposure | Exposure temperature: `assumed` 23 °C |
| §5.3.2 | Non-intermittent, illuminance-scale sensitometer | Virtual sensitometer: one exposure time, modulated illuminance |
| §5.3.3 | Illuminant per ISO 7589: sensitometric daylight, studio tungsten or photoflood. The illuminant used should be stated in use instructions | Kodak states "Daylight" (F-4016) and nothing for Tri-X (F-4017); use ISO sensitometric daylight (1.5) |
| §5.3.5 | Modulator neutrality: over 400–700 nm, total range of spectral density ≤ 5 % of the mean or 0.03, whichever is greater. Over 360–400 nm: ≤ 10 % or 0.06. Steps ≤ 0.20 log H. Continuous wedge ≤ 0.04 log H per mm | Virtual modulator is spectrally neutral by construction; use 0.15 steps (1.6) |
| §5.3.6 | Exposure time between 5 s and 1/1000 s, "corresponding to the usage practice" | Kodak gives no time for the D1 curves. F-9 (legacy TX) gives 1/50 s. Use 1/50 s, tagged `assumed` (sibling §5.1) |
| §5.4.1 | Exposure-to-processing interval: 5–10 days (general-purpose films) or 4 h–7 days (professional films), at 23 °C / 50 % RH | **Latent-image keeping is part of the standard condition.** If the latent-image model has regression, the virtual test must include the delay or state that it omits it |
| §5.4.2 | The only processing requirement: ΔD(m→n) = 0.80. Process information must be published | Speed is defined at one contrast, so it is a development-series interpolation (Annex A) |
| §5.5 | "ISO standard **Type 2** diffuse transmission density", geometry per ISO 5-2, spectral per ISO 5-3 | See 1.4: **Type 2 is printing density, not visual** |
| §5.6.2 | Base + fog from an unexposed sample processed with the strip | Virtual: a zero-exposure patch in every run (the design already lists "zero exposure → fog only") |
| §5.6.3 | m: 0.10 above base + fog. n: log H_m + 1.30. Time chosen so ΔD = 0.80 | Implement exactly. Note the speed point is on the **Type 2** curve |
| §6.1 | S = 0.80/H_m; S° = 1 + 10·log₁₀(0.80/H_m) | — |
| §6.1 Table 1 | ISO speed bins are **0.10 wide in log₁₀ H_m** (1/3 stop). ISO 400 ↔ −2.75 … −2.66; ISO 100 ↔ −2.15 … −2.06 | Compare log H_m numerically, not by bin. T-Max 100 in D-76 sits at about −2.06, right on the 100/80 bin edge (sibling §5.3) |
| §6.2 | Product speed comes from the **arithmetic mean of log₁₀ H_m over batches** | — |
| §6.3 | Calibration must keep the **error in log₁₀ H_m below 0.05** | The standard's own accuracy target is 1/6 stop (section 4) |

**Contrast condition versus Kodak's aim.** ISO's ΔD = 0.80 over Δlog H = 1.30 is an average gradient of 0.615 between fixed points. Kodak develops to CI 0.56 for these films (F-4017 2021 p. 3; J-86). The simulator must compute G-bar_ISO and CI separately (1.3).

### 1.2 Correction to the sibling report: ISO 6 densitometry is Type 2, not visual

* **What ISO 6 says.** ISO 6:1993 §5.5 and the Figure 1 y-axis label ("ISO standard Type 2 transmission density") both specify **Type 2** (opened, p. 3).
* **What Type 2 is.** ISO 5-3:1995 §8.2.2 defines ISO Type 2 printing density, D_T(S_H : s₂), for printing onto non-colour-sensitised silver-halide material such as black-and-white paper (opened, p. 5). Its spectral product combines:
  * the average sensitivity of such a print material;
  * an ultraviolet filter with a sharp cut-off at 360 nm.
* **Visual is separate.** ISO visual density, D_T(S_H : V_T), is §8.1.
* **Not a geometric variant.** ISO 5-2:1991, the geometry part, has no "Type" classes; it defines diffuse density and f/4.5 and f/1.6 projection density (opened, p. 1–2). Type 2 is therefore a **spectral** type.

Anchor values for the Type 2 product from ISO 5-3:1995 Table 2 (log₁₀ Π₂, normalised to a 5.000 peak; opened, p. 4):

| Wavelength | log₁₀ Π₂ |
|---|---|
| 430 nm | 5.000 (peak) |
| 360 nm | 4.280 |
| 400 nm | 4.916 |
| 500 nm | 3.993 |
| 540 nm | < 1.000 |

**Consequence:**
* **Silver images.** A developed silver image is close to neutral, so its Type 2 and visual densities differ only slightly. Silver is not perfectly neutral, though, and fine-grained silver is warmer; the size of the difference for T-grain versus conventional silver is **not found**.
* **Base and fog tints.** Base dyes, antihalation residue and sensitising-dye stain are not neutral. ISO speed computed from a Type 2 curve can therefore differ slightly from a speed computed on Kodak's visual curves.
* **Validation convention.** ISO speed validation should use Type 2 density from the simulator's own spectral transmittance. Curve-shape validation against Kodak uses visual density. Any spectrally flat base model is `assumed`.

### 1.3 Contrast index versus average gradient

#### 1.3.1 Definitions

| Quantity | Definition | Source / confidence |
|---|---|---|
| Gamma γ | Slope of the straight-line part of the D–log H curve | Standard; recalled (James 1977, ch. 17, or Todd & Zakia). Not used for validation |
| ISO G-bar | [D(log H_m + 1.30) − D(log H_m)] / 1.30, with m at 0.10 above base + fog | ISO 6 §5.6.3 (opened). It equals 0.615 at the ISO speed condition |
| Kodak contrast index (CI) | Average slope over the part of the curve "used" by a normally exposed negative of an average scene. It was introduced to replace gamma because it includes the toe | Niederpruem, Nelson & Yule, "Contrast Index", Phot. Sci. Eng. **10**(1) (Jan/Feb 1966) 35–41. Citation **secondary**, confirmed via the text of US 3,772,019, which also describes CI as using "a longer portion of the curve" than gamma; paper **not opened** |

#### 1.3.2 The CI construction (confirmed in Kodak H-740 p. 9 at the spot-check; not verified against the 1966 paper)

> **Spot-check, 2026-11-21.** Kodak's *Basic Photographic Sensitometry Workbook*, publication H-740 (© Eastman Kodak Company, November 2006), p. 9 ("Contrast Index"), is at https://www.kodak.com/content/products-brochures/Film/Basic-Photographic-Sensitometry-Workbook.pdf. It was fetched to `.cache/datasheets/spotcheck_cv/h740.pdf`.
>
> It gives the procedure:
> 1. Mark a straight edge at 0.0, at 0.2 log E units, and 2.0 log E units further on at 2.2.
> 2. Draw a horizontal line at D-min.
> 3. Move the edge until the 0.0 mark touches the D-min line and the 0.2 and 2.2 marks both touch the curve.
> 4. CI is the slope of the line along the edge.
>
> This is the collinear 0.2/2.2 construction below. The radii are therefore **confirmed by a Kodak source**, so they no longer rest on a forum post. The "URL dead" remark below is outdated.
>
> Two points remain open:
> * The 1966 defining paper itself was still not opened.
> * The 0.10 + 2.0 chord approximation remains sourced only from the dofmaster course text.

The geometric construction found in secondary sources:

1. Plot D and log H on **equal scales** (1 density unit = 1 log H unit).
2. Take a centre point C on the base + fog line.
3. Draw arcs of radius **0.2** and **2.2** about C. They cut the curve at P₁ (toe) and P₂.
4. Slide C along the base + fog line until **C, P₁ and P₂ are collinear**.
5. CI is the slope of P₁P₂. The segment P₁P₂ is 2.0 units long.

Provenance of this construction:
* **Overlay meter.** The procedure with a transparent meter ("slide until the low-density arc meets the toe and the high-density arc reads the same value") is in the *Advanced Photography* course text at dofmaster.com, p. 2-17 (© 2006 Don Fleming). That page is a hobby site with no stated primary source. It is cited here only because nothing better was found, and it does not give the radii.
* **Radii.** The values 0.2/2.2 come only from a Photrio forum post seen in a search summary. **Treat the radii as unverified.**
* **Approximation from the same course text:** take the toe point at 0.10 above base + fog, swing a compass of span 2.0 to the curve, and take the slope of the chord.
* **How to confirm:** open Phot. Sci. Eng. 10, 35 (1966), or Kodak's *Basic Photographic Sensitometry Workbook*. Its old kodak.com URL is dead, and web.archive.org could not be fetched by my tools.

#### 1.3.3 What Kodak Z-133E gives, and what it is not

Opened: *Monitoring and Troubleshooting KODAK Black-and-White Film Processes*, Kodak publication Z-133E, © 2003, "Minor Revision 3/03". It was fetched from imaging.kodakalaris.com, which redirects to business.kodakmoments.com. The PDF bundles it with AJ-3 and R-20.

* **Strip and readings.**
  * The control strip is pre-exposed **T-MAX 400** film with five steps: D-min, TD (toe), LD, HD and D-max.
  * The steps are read in densitometer "visual" mode (Z-133 p. 4).
* **Formula 1:** CI = 0.128 + 0.267·D-min − 0.969·TD + 0.454·LD + 0.183·HD + 0.039·D-max.
* **Formula 2:** CI = (HD − LD)/2.26 + 0.10 on Z-133 pp. 6, 8 and 9.
  * **Spot-check correction:** Kodak's publication is internally inconsistent, and the p. 5 text layer is **not** garbled. The rendered p. 5 prints CI = (HD − **TD**)/2.26 + 0.10, and so does Worksheet 1 on p. 45, whose boxes are labelled "HD" and "TD".
  * Which version Kodak intended is not settled here. It does not matter for D1, because the strip formulas must not be applied to simulated curves (below).
* **Speed value:** 140.9 × (TD − D-min) + 335. This is "not an ISO/ASA speed" (p. 5).
* **Aims:** CI 0.58 for a diffusion enlarger, 0.43 for a condenser enlarger; speed 355; D-min 0.06.
* **What these formulas are.** They are **regressions specific to that strip's step exposures**. They show that Kodak computes CI from a full curve, and they give the CI tolerance scale used in section 4. They are **not** a general CI algorithm and must not be applied to simulated curves.
* **The Z-133 film table is not a film CI target.** The table on p. 43 lists, for example, D-76 + Tri-X 400 normal 0.59 ± 0.03 and D-76 + T-Max 100 normal 0.62 ± 0.03. These are the **control-strip** CI values at which the process suits that film. They are not the CI of the film itself, and they must **never** be used as a Tri-X or T-Max target. The film targets remain CI 0.56 (F-4017, J-86) and the datasheet CI(t) curves.

#### 1.3.4 Numerical test on Kodak's vector curves (estimate)

* **Data.** F-4016 (T-Max 100), June 2018, p. 8, top-right chart: "KODAK PROFESSIONAL Developer D-76; Small Tank; 20 °C (68 °F)", Daylight, Diffuse Visual.
* **Paths.** The curves were extracted with PyMuPDF `get_drawings()`:
  * path 114, dash [6.48 3.24], 6 min;
  * path 112, dash [3.24 1.62], 7.5 min;
  * path 113, solid, 15 Bézier and line items, 10 min.
* **Curve-to-time assignment.** Made by final density order, because the legend labels sit beside unlabelled samples.
* **Frame calibration.** The frame rectangle runs x 355.6 → 540.0 pt (log H −4 → +1) and y 241.5 → 57.0 pt (D 0 → 4). This was checked against the tick-label positions.
* **Kodak CI(t) for comparison.** Taken from the "D-76" small-tank CI curve on p. 9 (dash [6.48 3.24]; frame 4–18 min, CI 0.3–0.9). Its vertices run from (6.30 min, 0.560) to (9.10 min, 0.820), and interpolation gives 0.658 at 7.5 min.

| Time | log H_m | ΔD over 1.30 | G-bar_ISO | CI, 0.10 + 2.0 chord | CI, 0.2/2.2 collinear | Kodak CI(t) curve |
|---|---|---|---|---|---|---|
| 6 min | −2.011 | 0.712 | 0.547 | 0.543 | 0.542 | ≈ 0.535 (extrapolated; the curve starts at 6.30 min) |
| 7.5 min | −2.129 | 0.911 | 0.701 | 0.689 | 0.693 | 0.658 |
| 10 min | −2.212 | 1.142 | 0.879 | 0.876 | 0.895 | ≈ 0.93 (extrapolated beyond 9.10 min) |

**Reading the table:**
* The two CI constructions differ by at most 0.02, so this test cannot tell them apart on these curves.
* **7.5 min is the only interpolated comparison.** There the characteristic-curve CI is 0.03 above Kodak's CI(t) curve. A plausible cause is that the two figures come from different test runs or were smoothed separately.
* **6 and 10 min are extrapolations.** At 6 min the agreement is within 0.01. At 10 min the gap is about 0.04–0.05 (0.876 or 0.895 against ≈ 0.93), but Kodak's curve ends at 9.10 min, so this point is weak evidence.
* **Second check: time to CI 0.56.** Interpolating the characteristic-curve CI between 6 and 7.5 min gives CI 0.56 at about 6.2 min. Kodak's CI curve gives 6.30 min, and the F-4016 time table gives 6½ min.
* **Consequence.** About 0.03 CI is an **upper bound from a single point**. It combines the datasheet's internal inconsistency with the unknown error of the (recalled) CI construction. Training thresholds on CI should be stated as RMS ≤ 0.03 with a maximum ≤ 0.05 (section 4.2). A bare "within 0.03" would sit on this noise floor.
* **Speed check.** The log H_m values reproduce the sibling's speed check, which lands at the ISO 100/80 boundary.
* **Status.** Label the 0.2/2.2 algorithm "estimate; consistent with Kodak data to about 0.03 at the one interpolable point". It stays `recalled` until the 1966 paper is opened.
  * *Spot-check update:* the construction is now sourced to Kodak H-740 (2006) p. 9 (see 1.3.2), so tag it `literature` (Kodak H-740), not `recalled`.
  * *Spot-check reproduction:* the table above was recomputed independently from the same PDF (MD5 4a193fb16cbf43cd0e6d40d03c7ccd89, p. 8 paths 112–114, same frame) and matches to all three decimals.
  * The p. 8 legend samples do sit beside their labels: solid = 10 min, [3.24 1.62] = 7.5 min, [6.48 3.24] = 6 min. This agrees with the D-max ordering.
  * The p. 9 D-76 small-tank CI vertices (6.30 min, 0.560) → (9.10 min, 0.820) and the 0.658 value at 7.5 min were also reproduced.

The script was a scratch file and was not kept. The method is fully described above, so it can be reimplemented in `validate/`.

### 1.4 ISO 5-2 geometry and ISO 5-3 spectral conditions

**ISO 5-2:2009 diffuse transmittance density** (opened, clauses 5–6, pp. 2–5):

| Item | Value |
|---|---|
| Diffuse-influx mode | θᵢ = 0°, κᵢ = 90° (diffuse) and θₜ = 0°, κₜ = 10°. Realised tolerances: ±2° on the angles and κₜ = 10° ± 2° |
| Diffuse-efflux mode | The reverse arrangement. It is stated to be equivalent |
| Realisation | Parameter tolerances "for many applications have an effect of less than 0,01 on the density values" (§5) |
| Diffuser | Diffusion coefficient 0.92 ± 0.02. Spectral reflectance factor 0.55 ± 0.05. Emulsion faces the diffuser; diffuser side polished (§6.2–6.3) |
| Effect of opal vs. integrating sphere | Opal densities come out "slightly smaller" than integrating-sphere densities because of inter-reflection between opal and specimen (Introduction) |
| Sampling aperture | Measurements on areas below **0.5 mm** diameter count as microdensitometry and are out of scope. Illumination uniformity within 10 % (§6.2) |
| Projection density (not needed for D1) | f/4.5: κ = 6.4°; f/1.6: κ = 18.2° |

**NIST SP 250-99** (Zarobila, Vest & Nadal, August 2022, *Diffuse Visual Transmission Density*; opened):
* **Expanded uncertainty.** NIST's step-tablet calibration has **U(k = 2) = 0.007** in density (Table 6, p. 29). The dominant component is step non-uniformity, 0.003.
* **Diffusion coefficient.** Moving the diffusion coefficient within its realisable range (0.90–0.95) changes density by **+0.005 to −0.019** (Eq. 45, p. 27).

Consequence for the virtual densitometer:
* **Ideal geometry.** Compute the ideal diffuse density: D = −log₁₀ of the hemispherically collected transmittance for collimated input, or the equivalent diffuse-influx/directional-efflux quantity. Do not model the opal inter-reflection.
* **What the model must supply.** The Nutting / covering-power model must output **diffuse** density. Literature covering-power values are diffuse; the specular density (Callier) is larger and is not needed for D1.

**ISO 5-3:1995 spectral conditions** (opened pp. iv–5; the 2009 edition is current, clause 6 in 2009 corresponds to clause 8 in 1995):

| Item | Content |
|---|---|
| Influx for transmission, S_H | CIE illuminant A, modified in the infrared by a heat-absorbing filter. Table 1, normalised to 100 at 560 nm; it equals CIE A up to about 570 nm and falls below it beyond (e.g. 62 vs. 198 at 700 nm) |
| Spectral product | Π = S·s at 10 nm intervals, normalised to a peak of 100 000 and tabulated as log₁₀ (peak 5.000) (§6, §8) |
| ISO visual transmission density D_T(S_H : V_T) | The receiver compensates so that **S_H·V_T = S_A·V** at every 10 nm (§8.1). The visual product is therefore CIE A (Planck, 2856 K) × CIE V(λ) |
| Check (estimate) | I computed log₁₀[A(λ)V(λ)], normalised to 5.000 at 570 nm, with CIE V(λ) and Planck A (c₂ = 1.4388 × 10⁻² m·K). Values: 450 nm 3.091, 500 nm 4.278, 560 nm 4.989, 600 nm 4.902, 650 nm 4.238. ISO 5-3:1995 Table 2 gives 3.090, 4.276, 4.989, 4.902, 4.238. Agreement is within **0.002** |
| Type 1 printing (diazo/vesicular) and Type 2 printing (silver halide) | §8.2.1–8.2.2; Type 2 is described in 1.2 |
| Status A, Status M | §8.3–8.4, Tables 3 and 4. **Not in the preview: not found.** They are not needed for black-and-white D1; they are needed in D2 |

**Virtual densitometer formula.** For a spectral product Π(λ) and specimen spectral diffuse transmittance T(λ):

D = −log₁₀ [ Σ Π(λ)·T(λ) / Σ Π(λ) ], summed at ≤ 10 nm.

* **Visual:** Π = A·V. Both factors are public CIE data, so tag them `literature`.
* **Type 2:** Π from ISO 5-3 Table 2, from a licensed copy. Tag it `literature` with the ISO citation. Until a licensed copy is available, tag it `assumed`, and mark the gap in the provenance plan.
  * *Spot-check note:* the full Type 2 column (340–540 nm) is legible in the cached iTeh preview (p. 4). Licensing, not availability, is the open question; the user decides.

### 1.5 ISO 7589, illuminants for sensitometry

Opened: ISO 7589:1984 preview (cover, §0–6, Table 1, the figure). The current edition is ISO 7589:2002, which adds a printer illuminant for black-and-white paper. Its preview was refused (HTTP 403), so the 2002 values are **not verified**.

| Item | Content (ISO 7589:1984) |
|---|---|
| Photographic daylight (§3.5) | Typical daylight of correlated colour temperature **5500 K** (skylight plus sun at about 40° elevation), designated **D₅₅**. Data from Judd, MacAdam & Wyszecki, JOSA 54 (1964) 1031, via CIE Publication 15, Table 1.1.4 |
| ISO sensitometric daylight illuminant (§4.1) | S_λ = D₅₅(λ) × τ(λ). τ(λ) is the spectral transmittance of the **ISO standard camera lens (ISO 6728)**. Table 1 runs 350–690 nm at 10 nm |
| Lens anchors (Table 1) | τ = 0.00 at 350 nm; 0.42 at 380 nm; 0.74 at 400 nm; 0.94 at 440 nm; ≥ 0.97 from 460 to 640 nm; 0.94 at 690 nm |
| Studio tungsten (§4.2) | For 3200 K lamps; the effective distribution temperature is about 3050 K in studios. Values in Table 2 (present in the cached preview, PDF p. 7; spot-check saw only the table heading) |
| Photoflood (§4.3) | Black body at 3400 K, modulated by the standard lens. Values in Table 3 (present in the cached preview, PDF p. 8; spot-check saw only the table heading) |
| Tolerance (§5–6) | Each illuminant is judged by its ISO/SDI: three numbers, 100 × log₁₀ of blue, green and red responses weighted by W_B, W_G, W_R, with the smallest set to 0. Tolerances: **red minus green within ±3** and **blue minus green within ±4** |
| Scope | Applies to black-and-white as well as colour camera films (§1) |

Simulator consequences:
* **Exposure SPD.** Use S_λ = D₅₅ × τ_lens as the exposure SPD at the film plane.
  * D₅₅ can be generated from the public CIE daylight basis functions. The 1984 table cites CIE 15; small differences from a 5 nm D55 are below the SDI tolerance.
  * τ_lens has no public source except the ISO tables. *Spot-check note:* the complete 350–690 nm τ column (and D₅₅ and S_λ) is legible in the cached ISO 7589:1984 preview, Table 1, p. 4, so the values are available; whether they may be embedded is a licensing decision for the user. If not used, use a flat τ = 1 above 460 nm with a UV cut-off, tagged `assumed`. A panchromatic film's blue/UV response makes this matter only at the level of a few hundredths in log H. That is an estimate, to be checked once the table is in hand.
* **Range.** The 1984 table ends at **690 nm**. The project's grid runs to 780 nm. Extend with D₅₅ × 0.94, `assumed`; panchromatic sensitivity is near zero beyond 690 nm.
* **Unit chain.** H in lux·s must be computed from the **same** S_λ (after the lens) that drives photon absorption. The film's speed in lux·s then depends on that SPD, which is why Kodak and ISO name the illuminant.

### 1.6 Step tablets

| Source | Content | Confidence |
|---|---|---|
| SMPTE RP 14-1997, *Plotting Data from Sensitometric Strips Exposed on Type Ib (Intensity Scale) Sensitometers*, reaffirmed 2004, "Stable" (pub.smpte.org) | §3.3: step increments "should typically" be **0.15** diffuse visual density. Single-step departures of about 0.015 are not significant, but cumulative errors are. §4.3: data **must be plotted against actual step densities** unless (1) every adjacent increment is 0.150 ± 0.015, (2) every step is within 0.02 density or 2 % of nominal, whichever is greater, and (3) every step is within 0.02 of the best linear fit. §6.2: tablet densities change little over periods of up to two years | opened |
| Kodak Photographic Step Tablet No. 2 | 21 steps, density about 0.05–3.05, nominal increment 0.15 | secondary (retailer listing for the Tiffen/Kodak tablet; Kodak's own catalogue page does not give densities) |
| AAPM Report 94 (TG-22, 2006) | Sensitometers use 0.15 log relative exposure per step; "light gamma" is the step-to-step slope over 0.15 | opened (p. 7–9) |
| ISO 6 §5.3.5 | Steps ≤ 0.20 log H | opened |

Virtual sensitometer:
* **Increment.** Use exact 0.15 increments, so there is no tablet error.
* **Range.** Do not stop at 21 steps: Kodak's axes span log H −4 to +1 (5 log units, 34 steps at 0.15). Make the virtual tablet cover the figure's full log H range, plus a zero-exposure patch.
* **Spatial scale.** Each step is millimetres wide, while diffusion lengths are tens of µm. Each step can therefore be simulated as an independent 1-D depth problem, coupled only through tank depletion (section 3.4).
* **Edges.** Edge effects at step boundaries are irrelevant to the curve if densities are read at step centres. ISO 5-2 says to read in "a uniform area"; Z-133 says to read "in the center" of each step.

### 1.7 Reciprocity tests

* **No standard found** (searched ISO and SMPTE titles; **not found**). ISO 6 §5.3.6 only says the exposure time must be stated because speed depends on it.
* **How Kodak specifies reciprocity.** Both D1 datasheets publish reciprocity as tables of "exposure adjustment" (stops, or adjusted time) and, for Tri-X, a "development adjustment" (sibling report §2.2; `latent-image.md` §4). The **criterion density** at which the adjustment holds is **not stated** by Kodak (not found).
* **Classical method (recalled).** Plot log(I·t) needed to reach a fixed net density against log I or log t; this is the reciprocity-failure curve in James 1977 (ch. 5). Confirm against James.
* **Proposed virtual test**, criterion tagged `assumed`:
  1. Run the sensitometer at t = 10⁻⁵ … 10² s.
  2. Find the extra log H needed to hold net density at **0.10** (the speed point) and at about **0.6–1.0** (mid-tone).
  3. Report both, because Kodak's adjustment may correspond to either.
* **The development column.** Tri-X's +10 % … +20 % development column means the criterion is contrast as well as speed. Compare the simulated CI change at each exposure time with the implied correction.

---

## 2. Digitising datasheet curves

### 2.1 Vector PDFs: extract the paths, do not trace pixels

Libraries:
* **PyMuPDF**, `Page.get_drawings()` (opened docs: https://pymupdf.readthedocs.io/en/latest/page.html).
  * It returns one dict per path, with keys `items`, `rect`, `color`, `fill`, `width`, `dashes`, `type` ("s", "f" or "fs"), `seqno` and `layer`.
  * Item tuples are ("l", p1, p2) lines, ("c", p1, p2, p3, p4) cubic Béziers with control points p2 and p3, ("re", rect, orientation) and ("qu", quad).
  * `get_cdrawings()` is equivalent but returns plain tuples. The docs cite a 10× speed-up (2 s → 200 ms on one page).
  * Coordinates use MuPDF's **top-left origin, in points**, on the unrotated page. **Flip y.**
  * `extended=True` (v1.22) also exposes clip paths and groups, which helps discard curve segments clipped by the frame.
  * Versions: `get_drawings` dates from v1.18.0. PyMuPDF 1.27.2 is available on this machine; it is not a project dependency.
* **pdfplumber** (opened README: https://github.com/jsvine/pdfplumber), built on pdfminer.six.
  * It exposes `.lines`, `.rects` and `.curves`.
  * `curve["pts"]` holds the on-curve points; `curve["path"]` holds the full path with Bézier control points.
  * Other attributes: `stroking_color`, `dash`, `linewidth`.
  * Coordinates: `top`/`bottom` are measured from the page top, `y0`/`y1` from the page bottom.
* **pdfminer.six**, the underlying parser (LTCurve / LTLine / LTRect; recalled). Use it directly only if a PDF defeats the other two.

Recommended pipeline for `tools/digitize/`. Each step was checked on F-4016 p. 8–9 (section 1.3.4).

1. **Identify the frame.** The plot box is a stroked `re` path, width 0.432 pt in these Kodak files. The grey filled `re` panels are backgrounds; ignore them.
2. **Calibrate the axes from tick labels, not assumed ranges.**
   * Read text spans and their bounding boxes with `page.get_text("dict")`.
   * Pair each numeric label with its tick.
   * Fit x = a + b·value by least squares using **all** labelled ticks. The residuals give a self-check and should be < 0.1 pt.
   * **Kodak's log H axis uses bar notation.** The labels extract as "4.0 3.0 2.0 1.0 0.0 1.0", so the negative sign is lost. Encode it per figure (`−4 … +1`), and assert that the labels are monotonic after the sign fix.
3. **Log axes.** MTF (response %, cycles/mm) and some spectral plots are logarithmic. Map the labels as log₁₀(value) and fit linearly in that space; a two-tick check suffices. Density, log H, CI and time axes are linear.
4. **Select curve paths.**
   * Keep stroked paths whose rect lies inside the frame.
   * Match each to its legend by (dash array, line width).
   * The 2018 F-4016 uses widths 0.72 and 1.08 pt and dashes such as [6.48 3.24], [3.24 1.62], [6.48 3.24 1.00 3.24] and [16.19 3.24 3.24 3.24 3.24 3.24].
   * Where the legend gives only labels beside identical samples, as in the characteristic-curve charts, assign by the physically required order (longer development → higher D-max) and record that this was done.
5. **Sample the paths.**
   * Expand Béziers by evaluating the cubic (about 30 points per segment is ample).
   * De-duplicate shared endpoints.
   * Sort by x for functions of log H or time.
   * Keep the raw vertex list too: the CI curves have only 3–8 vertices, which says something about how Kodak drew them.
6. **Record provenance in each reference CSV:** publication, edition date, page, file MD5, path `seqno`, dash and width, frame coordinates, tick-fit residuals, and the curve-assignment rule.

### 2.2 Raster scans: older editions and any figure without vectors

Method:
* Use 600 dpi or better. Deskew, then fit a **projective** (8-parameter) map from at least 4 ticks per axis; scans carry skew and keystone.
* Extract the centreline column by column: take the intensity-weighted centroid of dark pixels on each line, or skeletonise.
* Bridge dashes by fitting a smoothing spline through the dash fragments.
* Report the tick-fit residuals as the calibration error.
* WebPlotDigitizer (Rohatgi; GPL-3 per its citation page, secondary) implements this interactively. A small OpenCV/scikit-image script is reproducible and fits the project's no-copying rule better.

Resolution budget for the page images already in `ds_img/` (estimate): they are 850 × 1100 px for a 612 pt page, i.e. 1.39 px/pt, about 100 dpi.
* *Spot-check correction:* this holds for `f4016_2018_p8/p9` and `kpf4017_p8–p10`. The `f4017_2007_*` images are **723 × 935 px**, about 85 dpi (1.18 px/pt). Their budget is coarser, and the F-4016 chart scale below must not be reused for that edition without measuring its own frame.

| Chart | Axis | Scale at 100 dpi | Half-pixel uncertainty |
|---|---|---|---|
| Characteristic curves | density | 46.1 pt per D → 64 px per D (6.4 px per 0.1 D) | ±0.008 D |
| Characteristic curves | log H | 36.9 pt per log unit → 51 px per log unit | ±0.01 log H |
| CI charts | CI | 256 px per 0.6 CI (≈ 0.0023 CI per px) | — |
| CI charts | time | 0.055 min per px | — |

* **Line width.** The curves are 0.72 pt wide, about 1 px at this resolution. Rendering at 600 dpi brings the quantisation error below 0.002 D, and the remaining error is then the printed line itself.

Published digitisation studies (these cover only extraction error, from raster images, in non-photographic fields):
* Jelicic Kadic et al., *J. Clin. Epidemiol.* 74 (2016) 119–123 (secondary, abstract):
  * Software extraction was 47 % faster than manual estimation.
  * Agreement with the original data was 75 % and 73 % for software versus 66 % and 69 % for manual estimation (two raters).
* Drevon, Fursa & Malcolm, *Behavior Modification* 41(2) (2017) 323–339 (secondary, abstract): high intercoder reliability and validity with WebPlotDigitizer over 3596 points.
* Neither gives an error in physical axis units applicable here, so section 2.2's own resolution budget is the better number.

### 2.3 Three layers of curve error, to be kept separate in the reports

| Layer | Size | Basis |
|---|---|---|
| (a) Extraction | Vector: about 0 (floating-point, sub-0.001 D). Raster at 100 dpi: ±0.01 D / ±0.01 log H. At 600 dpi: < 0.002 | 2.1–2.2 (estimate) |
| (b) Fidelity of Kodak's plotted line to Kodak's measurement | **Unknown**. Characteristic curves are 136-segment polylines, smooth like plotting software, so they are fitted curves rather than measured points. CI curves have 3–8 vertices (sibling §3.1). The characteristic-curve and CI figures disagree by up to 0.03 CI (1.3.4) | estimate |
| (c) Representativeness of the curve for any given film | Kodak's NOTICE: the curves are "representative of production coatings" and "do not represent standards or specifications" (F-4016 2007, F-4017 2007, F-9 2003; the 2016–2021 editions omit the notice). ISO 6 §4 and §6.2 likewise define speed as a mean over batches | opened |

Practical σ for fitting (section 3.2): σ_D ≈ 0.02–0.03 D, combining (b) and (c). Raster extraction (a) only adds to this if raster inputs are used.

---

## 3. Parameter estimation for the mechanistic model

### 3.1 Structure of the D1 fit

* **Parameter blocks:**
  * film block θ_F for each film, about 10–20 parameters (grain classes, nucleation and growth, covering power, fog rate, …);
  * developer-specific block θ_B for each bath, covering formula uncertainties such as the HC-110 msds-estimated composition;
  * shared chemistry block θ_C (rate constants, adsorption constants, regeneration rates, activation energies).
* **Observables, with very different information content:**
  * D–log H curves: Tri-X × D-76 large tank in 35 mm and 120; T-Max 100 × D-76 small tank.
  * CI(t) curves.
  * Time-to-CI tables at 18–24 °C.
  * J-109 multi-CI times for XTOL.
  * Push tables, and F-4016 temperature factors.
* **Per sibling §7:** only D-76 provides curve shape for both films. The XTOL and HC-110 constants are constrained by CI(t) plus tables, roughly **10–20 scalar targets per developer**. Some of those are not independent: the time table entries are read off the same CI curve.

### 3.2 Objective function

* **Weighted least squares with priors (MAP):**

  χ²(θ) = Σᵢ [(yᵢ − fᵢ(θ))/σᵢ]² + Σⱼ [(gⱼ(θ) − μⱼ)/τⱼ]²

  The second sum holds the literature or plausibility priors on transformed parameters gⱼ (log rate constants, logit fractions). This is Tikhonov regularisation toward the priors, and it is the Gaussian-prior MAP estimate (recalled standard result, e.g. Tarantola, *Inverse Problem Theory*, SIAM 2005, ch. 3; confirm there).
* **Give each observable type its own σ before concatenating residuals:**

  | Observable | Residual unit | Suggested σ | Basis |
  |---|---|---|---|
  | Density on a characteristic curve | D | 0.02–0.03 | 2.3 (b) + (c). Weight the toe more lightly if b + f is uncertain |
  | CI at a given time | CI | 0.02–0.03 | Z-133 tuning tolerance ±0.02; datasheet consistency 0.03 (1.3.4) |
  | Time to reach a CI aim | log₁₀ t | ≈ 0.01–0.02 | Kodak rounds times to ¼ min. ±0.125 min at 6 min is ±2 %, i.e. a uniform rounding σ of about 0.005 in log₁₀ t. Add the CI-aim tolerance mapped through dCI/dt |
  | Speed point log H_m | log H | 0.03–0.05 | ISO 6 §6.3 |

* **Fit time in log t.** The residual is then a relative error, comparable across 4-minute and 15-minute developments.
* **Do not fit time-table entries and the CI curve they were read from as independent data.** Down-weight them, or keep one. Otherwise the fit double-counts one measurement.

### 3.3 Optimisers

* **Local, the main workhorse:** `scipy.optimize.least_squares` (opened docs).
  * `method="trf"` (Branch, Coleman & Li, SIAM J. Sci. Comput. 1999) supports bounds and sparse Jacobians. `"lm"` does not support bounds.
  * `x_scale="jac"`.
  * Use a robust `loss="soft_l1"` with `f_scale` of about 1 (in σ units) during exploration. Report final results with `loss="linear"` so χ² keeps its meaning.
  * Parametrise positive quantities in log space so the bounds become box constraints on log values.
  * Finite-difference Jacobians (`"2-point"` or `"3-point"`, `diff_step`) need a forward model whose numerical noise is far below the step. Use fixed time steps and grids, or tolerances 10–100× tighter than for production runs, while fitting. Adaptive-step noise produces spurious Jacobian columns.
* **Global, on a reduced set only:**
  * differential evolution (Storn & Price, *J. Global Optim.* 11 (1997) 341–359, confirmed; `scipy.optimize.differential_evolution`);
  * CMA-ES (Hansen, "The CMA Evolution Strategy: A Tutorial", arXiv:1604.00772, confirmed; the `cma` package).
  * Both need about 10²–10⁴ evaluations per dimension-scale, so restrict them to the ~5–8 identifiable combinations found in 3.4.
  * For the full set, multi-start TRF from Sobol points (`scipy.stats.qmc`) inside the prior box is cheaper and usually sufficient.
* **Cost control:**
  * Treat each sensitometer step as a 1-D depth column; sibling geometry gives mm steps against tens of µm of diffusion. Vectorise all steps of all curves of one bath into one array.
  * The only cross-step coupling is the finite-tank bulk concentration. Bromide release and developer consumption from a 35 mm strip in a small tank are small but not zero, so keep it.
  * Cache the latent-image look-up table per film and illuminant; it is independent of chemistry.
* **Staging:**
  1. Fit θ_F(Tri-X) and θ_F(T-Max 100) on the D-76 curves, with θ_C at prior values.
  2. Fit θ_C (+ θ_B) jointly on all training CI(t) and table data, with θ_F frozen.
  3. Do a final joint polish of all training data, starting from stages 1 and 2.
  4. Compare stage 3 with stage 1+2. A large shift in θ_F during the joint polish shows that θ_F is absorbing chemistry misfit (see 3.6).

### 3.4 Identifiability and sloppiness

* **Structural non-identifiability.** Some parameters enter only as products or ratios. Candidates in this model:
  * an Arrhenius pre-factor and an activation energy, if all data are at one temperature (use the 18–24 °C tables);
  * a boundary-layer thickness and an effective diffusion coefficient in the mass-transfer term;
  * developer-ion concentration times rate constant, where pKa and pH are fixed;
  * grain number density times covering power, for D-max;
  * regeneration rate and the semiquinone/oxidised-developer lifetime.
* **Practical non-identifiability: the profile likelihood.** Raue et al., *Bioinformatics* 25(15) (2009) 1923–1929 (confirmed).
  * For each parameter, re-optimise all others on a grid and plot χ²_min(θₖ).
  * A profile that crosses Δχ² = 3.84 on both sides gives a 95 % pointwise confidence interval.
  * A flat side means the parameter is practically non-identifiable: fix it to its prior, or report a one-sided bound.
  * Each profile costs a series of local fits. Do it only for parameters you intend to report as `fitted`.
* **Sloppiness: the eigenspectrum of the weighted Fisher information Jᵀ W J at the optimum.**
  * Gutenkunst et al., *PLoS Comput. Biol.* 3(10) (2007) e189 (confirmed), found eigenvalues spread evenly over many decades in 17 biology models.
  * They argue predictions can be constrained even when parameters are not.
  * Transtrum et al., *J. Chem. Phys.* 143 (2015) 010901 (confirmed) review the topic.
  * Procedure: SVD of W^{1/2}J in log-parameter space; inspect the singular values and right singular vectors.
  * The stiff directions are what the data determine. Report them as "fitted combination" parameters, if they are reported at all.
* **What to expect in D1:**
  * Stiff: the film parameters from D-76 curves.
  * The XTOL and HC-110 chemistry blocks are expected to show plateaus. For HC-110 this is made worse by the `msds-estimated` formula, whose composition uncertainty overlaps the rate constants. This is an estimate based on the data inventory, to be verified by the first profile run.
  * Recommendation: decide the fitted set **after** a profile and SVD pass on the training data. Fix the rest at priors, tagged `assumed` or `literature`. Re-fit, and only then touch the held-out data.

### 3.5 Reporting fitted values with uncertainty

* **Covariance.** Estimate it as C ≈ s²(JᵀWJ)⁻¹ with s² = χ²_min/(N − p) (recalled: Bates & Watts, *Nonlinear Regression Analysis and Its Applications*, Wiley 1988; Seber & Wild 1989). The scipy documentation states only that JᵀJ is the Gauss–Newton Hessian approximation; the s² scaling must come from a statistics text.
  * Report it in log space: "k = 3.1 × 10⁻² (×/÷ 1.4)".
  * Where the profile likelihood is asymmetric or open, report the profile interval instead.
* **Model discrepancy.** With an imperfect mechanistic model, residuals contain model discrepancy, and parametric intervals are too narrow and biased.
  * Kennedy & O'Hagan, *J. R. Stat. Soc. B* 63(3) (2001) 425–464 (confirmed) introduced the discrepancy term.
  * Brynjarsdóttir & O'Hagan, *Inverse Problems* 30 (2014) 114007 (confirmed) show that physical parameters are recovered only with realistic discrepancy priors.
  * Practical consequence: report fitted constants as "calibrated under model M version X". Treat their uncertainty as a lower bound. The meaningful uncertainty is the **held-out prediction error**, which the design already treats as the credibility measure.
* **Prediction uncertainty.** Propagate parameter uncertainty into predicted curves by sampling θ from C (or from profile-based intervals) and re-running. Report bands on held-out predictions alongside the point prediction.
* **Language.** Use GUM terms: standard uncertainty, expanded uncertainty with k. JCGM 100:2008 (*Evaluation of measurement data — Guide to the expression of uncertainty in measurement*) is recalled, not opened, and is freely available from BIPM.
* **Data fields** for each `fitted` value: value, σ or interval (with method: covariance or profile), fit ID, training-set hash, and whether the value lay on a bound (`active_mask` from least_squares).

### 3.6 Train / held-out protocol and what leakage looks like here

Protocol:
1. **Pre-register.** Before any fit, write down:
   * the held-out data: T-Max 100 × XTOL (all CI(t) curves, J-109 rows for 100TMX, F-4016 tables) and T-Max 100 × HC-110 (CI(t) curves, J-24 and F-4016 rows, F-4016 temperature factors);
   * the metrics (section 4);
   * the thresholds.
2. **Develop and select the model using internal cross-validation inside the training set.** For example:
   * fit chemistry on Tri-X × {D-76, XTOL} and predict Tri-X × HC-110;
   * leave out one temperature row (e.g. 24 °C), or the 1:1 dilution, and predict it.
   Use **block** splits, not random points (Roberts et al., *Ecography* 40 (2017) 913–929, confirmed: random CV underestimates error when data are dependent). Points within one CI curve are strongly dependent.
3. **Freeze:** code commit, data-file hashes and fitted-parameter file. The run manifest already records these.
4. **Run the held-out evaluation once.** Report every metric, including failures, and the `validated`/`failed` status as the design specifies.
5. **Any later change to the model is a new model version.** Its held-out numbers must be labelled "after unblinding".

Leakage routes specific to D1:
* Using any **T-Max 100 row of J-109 Table 1 or J-24** (time to CI for 100TMX in XTOL or HC-110) anywhere in training, including as a "sanity check" during development.
* **Re-tuning shared physics after seeing held-out residuals.** This includes the agitation/boundary-layer model, the tank-volume assumption, the msds-estimated HC-110 composition, the XTOL patent-derived formula and the developer pKa choices. These parameters are shared, so a held-out-motivated change leaks.
* **Using the D-76 small-tank CI(t) curve of T-Max 100** to calibrate agitation, then claiming the held-out small-tank XTOL curve as independent. It is not fully independent; it shares the agitation regime. State this, and check large-tank held-out curves separately.
* **Choosing which parameters to fix** (3.4) by looking at held-out fit quality.
* **Holding out random time points inside a curve** instead of whole curves, temperatures or dilutions.
* **Fitting temperature dependence** on T-Max 100 temperature factors (F-4016 p. 7). Those for D-76 are training; those for HC-110 are held out.
* **Selecting the CI construction or the densitometer spectral type** after seeing which gives better held-out agreement.
* **Film-parameter leakage through the joint fit.** In a joint fit, Tri-X film parameters see XTOL and HC-110 data, but T-Max 100 film parameters see only D-76. This asymmetry is intended. It also means held-out success depends on whether the chemistry generalises across grain types (T-grain versus conventional) — that is the hypothesis under test.

---

## 4. What accuracy is realistic, and proposed thresholds

### 4.1 Reference scales

| Quantity | Scale | Source | Confidence |
|---|---|---|---|
| Density, primary calibration | U(k = 2) = 0.007 D | NIST SP 250-99 (2022), Table 6 | opened |
| Density, instrument realisation | < 0.01 D for ISO 5-2 tolerances; +0.005/−0.019 D across the allowed diffusion coefficients | ISO 5-2:2009 §5, §6.3 note 1; NIST SP 250-99 Eq. 45 | opened |
| Densitometer short-term repeatability | ±0.01 D (0–5 D, visual, 2–3 mm aperture) | X-Rite 361T support page | secondary |
| Step-tablet conformance | ±0.015 per increment; ±0.02 D per step | SMPTE RP 14-1997 §4.3 | opened |
| Repeat sensitometric exposures (same film batch, same processor) | Re-test if they differ by > 0.04 D (base + fog) or > 0.10 D (straight line) | AAPM Report 94 (2006) §3.2.3, printed p. 5 (spot-check page correction); the same 0.04 / 0.10 thresholds reappear in §3.3.9, p. 7, there for "speed density" | opened |
| Base + fog drift worth investigating | > 0.04 D, same film and batch | AAPM Report 94 p. 9 | opened |
| ISO speed accuracy | Error in log H_m < 0.05 (1/6 stop); speed bins 0.10 wide (1/3 stop) | ISO 6:1993 §6.3, Table 1 | opened |
| Kodak black-and-white process control (T-MAX 400 control strip) | CI: tuning ±0.02 of aim; action limits ±0.07 (upper action region +0.07 to +0.20, lower −0.07 to −0.12); control limits +0.20 / −0.12. Speed value: action −9/+15, control −17/+22 units of 140.9 × (TD − D-min), i.e. toe-step density about −0.06/+0.11 D (action) and −0.12/+0.16 D (control). D-min: action +0.02 | Kodak Z-133E (2003) pp. 6, 9 | opened (conversion of speed units to density: estimate) |
| Kodak film-to-strip CI mapping scatter | ±0.03 to ±0.08 CI per film ("Normal" column; the "Push 2" column reaches ±0.11, PLUS-X 125 in T-MAX RS) | Z-133E p. 43 | opened |
| Datasheet internal consistency (characteristic-curve figure versus CI(t) figure) | 0.03 CI at 7.5 min, the only interpolable point (T-Max 100, D-76 small tank). This is an upper bound that includes the construction's error. Time to CI 0.56: about 6.2 min (characteristic curves), 6.30 min (CI curve), 6½ min (table) | 1.3.4 | estimate |
| Datasheet speed vs box speed | Tri-X 35 mm D-76: log H_m −2.73 (inside the ISO 400 bin). T-Max 100: −2.06 (edge of the 100 bin) | sibling §5.3; re-computed in 1.3.4 | estimate |
| Batch-to-batch film variation (radiographic film, not pictorial) | For films of the same type from different batches, the X-ray exposure (through intensifying screens) needed for net OD 1.0 varied by ±20 %, i.e. +0.08 / −0.10 log H. The study covered 50 batches from 20 institutions in Ontario, of several types, so 50 is not a within-type count. The abstract also says white-light (Wejex) sensitometry was not always a good indicator of X-ray sensitivity, which weakens this as a pictorial analogue | C. Daniels, A. F. Holloway, "Comparison and variations of the speed of radiographic film", *Radiology* 146(1) (1983) 203–207, DOI 10.1148/radiology.146.1.6849046 (PMID 6849046) | abstract opened (PubMed E-utilities) at spot-check; full text not opened |
| Interinstitutional variation (radiographic dosimetry film) | OD differences of 17 % (batch) and 25 % (processing) at 50 cGy | Bos et al., *Med. Phys.* 2002 (PubMed 12201424) | secondary (abstract) |
| Pictorial black-and-white batch or interlaboratory reproducibility | — | — | **not found** |

### 4.2 Proposed thresholds for the first D1 report

The design leaves thresholds to be agreed with the user after the first report. The table below is a starting proposal derived from 4.1, not a decision.

| Metric | Training pairs | Held-out pairs | Derivation |
|---|---|---|---|
| RMS density error on a D–log H curve, between b + f + 0.1 and the last plotted point | ≤ 0.03 D | (no held-out curves in the plan) | Above metrology (0.007–0.01), at the datasheet layer (b) + (c) scale, well inside the AAPM 0.10 straight-line flag |
| Max density error on a curve | ≤ 0.08 D | — | Below the AAPM 0.10 flag |
| Base + fog | within 0.03 D | — | AAPM 0.04; Z-133 D-min action +0.02 |
| Speed, log H_m at ΔD = 0.80 (Type 2) | within 0.05 (ISO 6 §6.3) | within 0.10 (one 1/3-stop bin) | ISO 6 |
| CI at a given time | RMS ≤ 0.03 and max ≤ 0.05 | each point within 0.07 | 0.03 = single-point upper bound on datasheet consistency (1.3.4); 0.07 = Kodak action limit |
| Time to CI 0.56 | equivalent of RMS ±0.03 CI via the local dCI/dt | equivalent of ±0.07 CI | See below |

Converting the CI tolerance to time uses the local slope dCI/dt of the digitised CI(t) curve (vertex spans from sibling §2.1 and §3.1; estimate):

| Pair (small tank, 20 °C) | dCI/dt (per min) | Time equivalent of ±0.07 CI | Time equivalent of ±0.03 CI |
|---|---|---|---|
| T-Max 100 × D-76 | 0.26 / 2.8 min = 0.093 | ±0.75 min (±12 % of 6.5 min) | ±0.32 min |
| T-Max 100 × XTOL | 0.26 / 3.3 min = 0.079 | ±0.9 min (±12 % of 7.5 min) | ±0.38 min |
| T-Max 100 × HC-110 (B) | 0.26 / 8.4 min = 0.031 | **±2.3 min** (±38 % of 6 min) | ±1.0 min |
| Tri-X × D-76 | 0.34 / 5.4 min = 0.063 | ±1.1 min | ±0.48 min |
| Tri-X × XTOL | 0.64 / 10.85 min = 0.059 | ±1.2 min | ±0.51 min |
| Tri-X × HC-110 (B) | 0.53 / 8.65 min = 0.061 | ±1.1 min | ±0.49 min |

* **HC-110 and T-Max 100.** Kodak's CI(t) curve is very flat, so a time-to-CI test is weak; the time tolerance becomes large. Judge this pair by **CI at the tabulated times**, not by time-to-CI.
  * *Spot-check:* the span was re-measured from F-4016 (2018) p. 9, path with dash [16.19 3.24 3.24 3.24 3.24 3.24]. It runs from 5.90 min at CI 0.560 to 14.30 min at 0.820, confirming 0.031 per min.
  * The **local** slope on the first vertex segment, at CI 0.56, is lower: about 0.028 per min. That gives about **±2.5 min** for ±0.07 CI.
  * The T-Max × D-76 (0.093) and × XTOL (0.079, 7.30 → 10.60 min) spans were also reproduced.
* **Slopes are averages.** These slopes are span averages. Use the local slope at CI 0.56 from the vector vertices in the implementation.

---

## 5. Existing open implementations

| Project | What it does | Approach to curves/calibration | Limits for this project | Confidence |
|---|---|---|---|---|
| **spektrafilm** (formerly agx-emulsion), A. Volpato, https://github.com/andreavolpato/spektrafilm | Spectral simulation of colour negative film and print: layer spectral sensitivities, per-layer D–log H curves, dye spectra, masking couplers as negative absorption, DIR couplers with diffusion, halation and scatter as Gaussian blurs, stochastic grain | Uses curves "extracted" from manufacturer datasheets; the extraction tool is not named. Generic CMY diffuse densities where datasheets omit them. Print paper uses "unmixed" density curves. No development chemistry. Code GPL-3; profiles and LUTs CC BY-SA 4.0. Used by darktable's spektrafilm module and by ART | Curve-driven (approach A in the design), so no developer prediction. The README says datasheet curves alone do not give a convincing look. **Do not copy its data** (project rule; also share-alike) | opened (README) |
| **Filmulator**, C. Vaculik, https://github.com/CarVac/filmulator-gui | Raw editor that simulates developer depletion and diffusion to produce local-contrast and highlight-compression effects | Exposes "film area" and "drama" controls, not film stocks or chemistry; there is no characteristic-curve calibration | No film or developer data; not a sensitometric tool. Licence file present but its type is unverified | opened (README); licence not verified |
| **WebPlotDigitizer**, A. Rohatgi, https://automeris.io | Interactive raster or vector-image digitiser with automatic curve extraction | Manual axis calibration; colour- or mask-based line extraction | Raster-only workflow; for Kodak vector PDFs, direct path extraction (2.1) is exact and scriptable | secondary |
| Python sensitometry or H&D-fitting package with a documented method | — | — | **Not found.** Forum posts describe ad hoc SciPy spline fits; commercial tools (e.g. BTZS "Plotter", ExpoDev) are closed and were not examined | not found |

Related classical method, useful for compact curve families and residual analysis: Simonds, "Application of characteristic vector analysis to photographic and optical response data", *JOSA* 53 (1963) 968. This is principal-component analysis of families of response curves. Secondary; citation confirmed via reference lists, paper not opened. It could serve as a sanity check that simulated curve families span the same low-dimensional space as Kodak's.

---

## 6. Consequences for D1

1. **Virtual densitometer:**
   * Implement spectral diffuse transmittance density with two products: **visual** (CIE A × V(λ), public data, verified against ISO 5-3) and **ISO Type 2** (ISO 5-3 Table 2, needs a licensed copy).
   * Compare with Kodak curves in visual density.
   * Compute ISO speed in Type 2 density.
   * Correct sibling `datasheets.md` §5.2 accordingly.
2. **Virtual sensitometer:**
   * Illuminant: S_λ = D₅₅ × τ_ISO-lens. τ is `assumed` until a licensed ISO 7589 copy is available; extend the 690 nm table end to 780 nm, `assumed`.
   * Exposure time: 1/50 s, `assumed` from F-9.
   * Steps: 0.15 log H spanning the full figure axis (−4 … +1), plus a zero-exposure patch.
   * Model each step as an independent 1-D depth column coupled only through the tank.
   * Optionally simulate the ISO 6 latent-image keeping interval.
3. **Metrics:**
   * Compute G-bar_ISO and log H_m (ISO 6).
   * Compute CI with the 0.2/2.2 collinear construction. Tag it `literature` (Kodak H-740, 2006, p. 9; spot-check), "consistent with Kodak curves to 0.03".
   * Also report the 0.10 + 2.0 chord as a cross-check.
   * Never use the Z-133 strip formulas or strip targets as film targets.
4. **Digitisation tool:**
   * Vector path extraction with tick-label axis fitting, the bar-notation sign fix, legend matching by dash and width, and a full provenance record.
   * Raster only for editions without vectors, at ≥ 600 dpi with a projective calibration.
5. **Fit design:**
   * Objective: per-observable σ; time-to-CI fitted in log t; no double counting of tables and curves.
   * Optimisation: log-parameter TRF with bounds and multi-start; DE or CMA-ES only on the identified stiff subset; staged fit (film on D-76, then chemistry, then joint polish).
   * Before any held-out look: a profile-likelihood and Jacobian-SVD pass to fix non-identifiable parameters to priors.
6. **Held-out test:**
   * Pre-register metrics and thresholds; freeze hashes; evaluate once.
   * Use block internal CV inside training (leave out a developer, a temperature row or a dilution).
   * Respect the leakage list in 3.6, in particular the T-Max 100 rows of J-109 and J-24 and any retuning of shared agitation or formula assumptions.
7. **Thresholds** (to agree with the user):
   * Curves: RMS ≤ 0.03 D, maximum ≤ 0.08 D.
   * Speed: log H_m within 0.05 (training) or 0.10 (held-out).
   * CI: RMS ≤ 0.03 with maximum ≤ 0.05 (training); each point within 0.07 (held-out).
   * Express T-Max 100 × HC-110 as CI at time, because its CI(t) is too flat for a time-to-CI test.
   * The datasheet's self-consistency is bounded by about 0.03 CI at its one interpolable point.
8. **Reported uncertainties:**
   * Parametric σ and profile intervals for `fitted` values, stated as lower bounds under model discrepancy.
   * The held-out prediction error is the headline credibility number.

---

## 7. Not found

* Text of Niederpruem, Nelson & Yule (1966). *Spot-check:* Kodak's *Basic Photographic Sensitometry Workbook* (H-740, 2006) **was found** at kodak.com and confirms the 0.2/2.2 construction (p. 9), so the radii are no longer unverified.
* Status A and Status M spectral products (ISO 5-3 Tables 3–4; not in the preview).
* ISO 7589:2002 values; the preview was refused. (*Spot-check:* the 1984 Tables 2–3, studio tungsten and photoflood, **are** in the cached 1984 preview, PDF pp. 7–8, so they are no longer "not found"; their values were not read here.)
* ISO 6728 standard-lens transmittance beyond the anchors quoted from the ISO 7589:1984 Table 1.
* A standard or Kodak statement of the criterion density behind reciprocity tables.
* Published batch-to-batch or interlaboratory reproducibility for pictorial black-and-white film. Only radiographic-film analogues were found, from abstracts.
* Size of the visual versus Type 2 density difference for developed silver in Tri-X and T-Max 100.
* An open-source package that fits mechanistic or H&D-curve parameters with a documented method.
* Filmulator's licence type.

---

## 8. Sources

Standards and metrology:
* ISO 6:1993, *Photography — Black-and-white pictorial still camera negative film/process systems — Determination of ISO speed*. Preview via iTeh: https://cdn.standards.iteh.ai/samples/3580/192a0f97667044f3ba51276b0ded056a/ISO-6-1993.pdf. Opened: §1–7, Fig. 1, Table 1.
* ISO 5-2:2009 preview: https://cdn.standards.iteh.ai/samples/52914/b670fdd035aa4ae48069396eae658ba7/ISO-5-2-2009.pdf (opened, §1–6.4). ISO 5-2:1991 preview: https://cdn.standards.iteh.ai/samples/20037/70520650713f43c2af7955400729998e/ISO-5-2-1991.pdf (opened, §1–5.2).
* ISO 5-3:1995 preview: https://cdn.standards.iteh.ai/samples/20036/904283a1baf2405bb2324932b7523465/ISO-5-3-1995.pdf (opened, Introduction, §1–8.4, Tables 1–2).
* ISO 7589:1984 preview: https://cdn.standards.iteh.ai/samples/14377/92d4193ebb264e2c8dd773bad38ef097/ISO-7589-1984.pdf (opened, §0–6, Table 1, the figure). ISO 7589:2002 catalogue page: https://www.iso.org/standard/33979.html (secondary).
* SMPTE RP 14-1997: https://pub.smpte.org/doc/rp14/19970122-pub/rp0014-1997_stable2015.pdf (opened).
* C. J. Zarobila, R. E. Vest, M. E. Nadal, *Diffuse Visual Transmission Density*, NIST SP 250-99 (2022), https://doi.org/10.6028/NIST.SP.250-99 (opened).
* AAPM Report No. 94, *Validating Automatic Film Processor Performance*, TG-22 (2006): https://www.aapm.org/pubs/reports/RPT_94.pdf (opened).
* X-Rite 361T support specifications: https://xrite.com/361t-transmission-densitometer/Support (secondary).
* C. Daniels, A. F. Holloway, *Radiology* 146(1) (1983) 203–207, doi:10.1148/radiology.146.1.6849046; abstract via PubMed E-utilities, PMID 6849046 (opened at spot-check).

Kodak:
* *Monitoring and Troubleshooting KODAK Black-and-White Film Processes*, Z-133E, © 2003, Minor Revision 3/03: https://business.kodakmoments.com/sites/default/files/wysiwyg/pro/chemistry/z133.pdf (opened, pp. 3–9, 42–43).
* F-4016 *KODAK PROFESSIONAL T-MAX 100 Film*, Revised 6-18: https://kodakprofessional.com/sites/default/files/wysiwyg/pro/resources/f4016_TMax_100.pdf (opened; vector paths pp. 8–9).
* F-4016 (2007), F-4017 (2007), F-9 (2003) NOTICE text, via the cached text extracts in `.cache/datasheets/ds_txt/` (opened; publication URLs are in `datasheets.md` §1).
* *Basic Photographic Sensitometry Workbook*, Kodak publication H-740, © 2006: https://www.kodak.com/content/products-brochures/Film/Basic-Photographic-Sensitometry-Workbook.pdf (opened at spot-check, p. 9 Contrast Index; cached `.cache/datasheets/spotcheck_cv/h740.pdf`).
* US Patent 3,772,019 (citation of Niederpruem, Nelson & Yule, Phot. Sci. Eng. 10(1) 1966 35–41): https://patents.google.com/patent/US3772019A/en (opened).
* *Advanced Photography* course text, p. 2-17, https://dofmaster.com/courses/advanced/advancedphotography-55.html (opened; a secondary source for the CI approximation).

Digitisation:
* PyMuPDF `Page` documentation: https://pymupdf.readthedocs.io/en/latest/page.html (opened).
* pdfplumber README: https://github.com/jsvine/pdfplumber (opened).
* B. Jelicic Kadic et al., *J. Clin. Epidemiol.* 74 (2016) 119–123 (secondary).
* S. Drevon, S. R. Fursa, A. L. Malcolm, *Behav. Modif.* 41(2) (2017) 323–339, doi:10.1177/0145445516673998 (secondary).

Estimation and validation methodology:
* `scipy.optimize.least_squares` documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html (opened).
* R. Storn, K. Price, *J. Global Optim.* 11 (1997) 341–359, doi:10.1023/A:1008202821328 (confirmed via search).
* N. Hansen, "The CMA Evolution Strategy: A Tutorial", arXiv:1604.00772 (confirmed via search).
* A. Raue et al., *Bioinformatics* 25(15) (2009) 1923–1929, doi:10.1093/bioinformatics/btp358 (confirmed via search).
* R. N. Gutenkunst et al., *PLoS Comput. Biol.* 3(10) (2007) e189, doi:10.1371/journal.pcbi.0030189 (confirmed via search).
* M. K. Transtrum et al., *J. Chem. Phys.* 143 (2015) 010901, doi:10.1063/1.4923066 (confirmed via search).
* M. C. Kennedy, A. O'Hagan, *J. R. Stat. Soc. B* 63(3) (2001) 425–464, doi:10.1111/1467-9868.00294 (confirmed via search).
* J. Brynjarsdóttir, A. O'Hagan, *Inverse Problems* 30 (2014) 114007, doi:10.1088/0266-5611/30/11/114007 (confirmed via search).
* D. R. Roberts et al., *Ecography* 40(8) (2017) 913–929, doi:10.1111/ecog.02881 (confirmed via search).
* Recalled, not opened:
  * Bates & Watts (1988), *Nonlinear Regression Analysis and Its Applications*;
  * Tarantola (2005), *Inverse Problem Theory*;
  * JCGM 100:2008 (GUM);
  * T. H. James (ed.), *The Theory of the Photographic Process*, 4th ed. (1977).

Open implementations:
* spektrafilm: https://github.com/andreavolpato/spektrafilm (README opened).
* Filmulator: https://github.com/CarVac/filmulator-gui (README opened).
* J. L. Simonds, *JOSA* 53 (1963) 968 (secondary).

---

## 9. Spot-check log (independent verification, 2026-11-21)

How the check was done:
* Each fact was re-checked by opening the source myself.
* For the ISO, SMPTE, NIST, AAPM and Z-133 sources, I read the PDFs cached under `.cache/datasheets/`, extracting the text with PyMuPDF and rendering tables or formulas to images where layout mattered.
* Kodak H-740 was newly fetched from kodak.com.
* The Radiology abstract came from PubMed E-utilities.
* The F-4016 numbers were recomputed from the vector paths with a fresh script, which was not kept.

| # | Claim | Verdict | Evidence / correction |
|---|---|---|---|
| 1 | CI = slope of the collinear 0.2/2.2 construction from a centre on the b + f line (was `recalled`) | **confirmed** | Kodak H-740 (2006) p. 9: a straight edge marked 0.0 / 0.2 / 2.2 log E, with 0.0 on the D-min line and 0.2 and 2.2 on the curve. The 1966 paper is still not opened, and the 0.10 + 2.0 chord is still sourced only to the course text |
| 2 | Batch speed variation ±20 % (radiographic; was `secondary`) | **confirmed, context corrected** | PubMed abstract: Daniels & Holloway, *Radiology* 146(1):203–207. The ±20 % is X-ray exposure through screens to net OD 1.0, for same-type film across batches. The 50 batches span several types and 20 institutions. ±20 % = +0.08/−0.10 log H |
| 3 | F-4016 (2018) T-Max 100 × D-76 small tank: log H_m −2.011/−2.129/−2.212; CI(0.2/2.2) 0.542/0.693/0.895; chord 0.543/0.689/0.876; Kodak CI(t) 0.658 at 7.5 min | **confirmed** | Recomputed from paths 112–114 (p. 8) and the p. 9 D-76 path. Identical to 3 decimals |
| 4 | T-Max × HC-110 (B) dCI/dt ≈ 0.031 per min (±2.3 min for ±0.07) | **confirmed** | p. 9 path: 5.90 → 14.30 min for CI 0.56 → 0.82. The local slope at CI 0.56 is about 0.028 per min (≈ ±2.5 min); added in 4.2 |
| 5 | Visual product = CIE A × V(λ), within 0.002 of ISO 5-3 Table 2 | **confirmed** | Recomputed with Planck 2856 K: 450 3.091, 500 4.278, 560 4.989, 600 4.902, 650 4.238, against the table's 3.090/4.276/4.989/4.902/4.238. It also matches at 410–440 nm |
| 6 | Page images 850 × 1100 px ≈ 100 dpi; 6.4 px per 0.1 D; ±0.008 D at half a pixel | **confirmed, scope corrected** | True for the F-4016 2018 and kp F-4017 images. The F-4017 2007 images are 723 × 935 px (≈ 85 dpi) |
| 7 | ISO 6: m at 0.10 over b + f; n at +1.30; ΔD = 0.80; S = 0.80/H_m | **confirmed** | ISO 6:1993 §5.6.3, §6.1, preview p. 3 |
| 8 | ISO 6 densitometry is Type 2. Type 2 anchors: 360 nm 4.280, 400 nm 4.916, 430 nm 5.000, 500 nm 3.993, 540 nm < 1.000 | **confirmed** | ISO 6 §5.5 and §5.6.1. The ISO 5-3:1995 §8.2.2 and Table 2 render checked column by column |
| 9 | Speed bins 0.10 wide (ISO 400 −2.75…−2.66; ISO 100 −2.15…−2.06); accuracy < 0.05; product speed = mean of log H_m | **confirmed** | ISO 6 Table 1, §6.2, §6.3, preview p. 4 |
| 10 | Exposure-to-processing interval 5–10 days (general purpose) or 4 h–7 days (professional), at 23 ± 2 °C and 50 ± 5 % RH | **confirmed** | ISO 6 §5.4.1, p. 2 |
| 11 | ISO 7589:1984 lens τ: 0.00 at 350, 0.42 at 380, 0.74 at 400, 0.94 at 440, ≥ 0.97 from 460 to 640, 0.94 at 690 nm; D₅₅ = 5500 K | **confirmed** | Table 1, p. 4 (rendered). The whole table is legible in the preview, so "needs a licensed copy" was changed to a licensing note |
| 12 | Kodak Z-133 Formula 2: CI = (HD − LD)/2.26 + 0.10 on pp. 5, 6, 9; the p. 5 text layer is garbled | **wrong** | The rendered p. 5 prints (HD − **TD**)/2.26 + 0.10, and Worksheet 1 (p. 45) also uses HD and TD. Only pp. 6, 8 and 9 print (HD − LD). Kodak is inconsistent; the text layer is faithful |
| 13 | Z-133 tolerances: tuning ±0.02; action ±0.07; control +0.20/−0.12; D-min action +0.02 | **confirmed** | p. 9 table (CI action +0.07 to +0.20 / −0.07 to −0.12; D-min action +0.02 to +0.03); tuning ±0.02 on p. 6 and p. 43. The p. 43 normal-column scatter is ±0.03…±0.08, and Push 2 reaches ±0.11 (row qualified) |
| 14 | AAPM 94: repeat exposures differing by > 0.04 D (b + f) or > 0.10 D (linear) should be repeated; b + f drift > 0.04 investigated | **confirmed, page corrected** | §3.2.3 is on printed p. 5, not p. 7. §3.3.9 on p. 7 repeats 0.04 / 0.10 ("speed density"). The b + f drift item is on p. 9 |
| 15 | NIST SP 250-99: U(k = 2) = 0.007 D; d = 0.90–0.95 gives +0.005 to −0.019 D | **confirmed** | Table 6 (printed pp. 28–29); Eq. 45 (printed p. 27) |

Also re-read and consistent with the report; these were not counted in the 15:
* ISO 5-2:2009 geometry: κ 90°/10° ± 2°, d = 0.92 ± 0.02, R = 0.55 ± 0.05, < 0.5 mm aperture, < 0.01 D.
* ISO 7589 SDI tolerances: ±3 red–green, ±4 blue–green (rendered p. 3).
* SMPTE RP 14 §3.3 and §4.3.
* The Kodak NOTICE wording in F-4016 (2007), F-4017 (2007) and F-9 (2003); also present in F-32 (2002).
