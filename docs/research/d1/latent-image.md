# D1 research: latent-image formation and grain population data

Scope: inputs for stage 2 of the physics chain (latent-image formation, per grain size class) for Kodak Tri-X 400 and T-Max 100. Research date: 2026-11-20. Everything here is meant to be checked again before a value goes into `data/` with a provenance tag.

Revision 2026-11-21 (consistency pass, no new research): aligned the τ₀ prior with its own derivation (§1.2, §1.4); corrected the reading of the Hailstone & De Keyzer 2004 abstract (§1.2); corrected the grain-size scaling of intrinsic absorption (§1.3d); stated that the spectral-curve reference exposures lie at the start of LIRF, not inside the reciprocity-normal range (§3.3, §5.2); corrected §5.6, which wrongly said F-4017 has no D-76 characteristic curves; added the Tri-X 5063 provenance note and a tabular-AR caveat on the T-Max 100 prior.

Revision 2026-11-21 (independent spot-check of the key facts). Abstracts were re-read from publisher metadata (Crossref) and OpenAlex, and the Kodak sheets and the Kofron patent were re-opened. Findings:
- **Upgraded:** Fayet et al. 1985 (Ag₄⁺ is the smallest developable deposited cluster), Katz 1949 (Ag₂), both Hailstone 2004 abstracts, and the 12-photon residual-size result now count as opened abstracts.
- **Corrected:**
  - The 12-photon paper is Broadhead & Hymers (1987), not "Hailstone et al.".
  - The hole mobility is a measured Hall value of 2.0 ± 0.5 cm² V⁻¹ s⁻¹, not ~1.
  - The HIRF Schwarzschild exponents are recomputed for aperture (intensity) adjustment: 1.15 / 1.10, not 1.18 / 1.11.
  - The T-Max 100 sheet actually re-checked is the June 2018 F-4016.
- **Added from the opened abstracts:**
  - Katz's requirement of a distribution of trap depths.
  - Experimental LIRF onset at 0.1–1 s, a 3-atom minimum developable size and measured sensitizer trap depths (Hailstone, French & De Keyzer 2004).
  - A φ_r caveat (Tani 1980, near unity).
- Items still unverified are marked **UNVERIFIED** in place.

## Confidence labels

| Label | Meaning |
|---|---|
| **opened-source** | I opened the primary document and read the value myself, either as extracted PDF text or as a rendered page image. |
| **secondary-source** | I saw the value only in a search-result abstract or snippet, or in a later document that cites the original. The original paper was not opened. |
| **recalled** | From my memory of the literature; not opened in this session. Each one says what would confirm it. |
| **estimate** | Computed here from opened values plus standard constants. The arithmetic is shown. |

Access problems: Taylor & Francis (J. Photogr. Sci., Imaging Sci. J.), the IS&T library (J. Imaging Sci. Technol.), the RIT repository PDFs and the Katz 1949 PDF on Deep Blue all returned HTTP 403 or a JavaScript wall. As a result, **none of the kinetic constants of the Hamilton / Hailstone models could be read from a primary source**. They are listed below as recalled or as fit targets. (Spot-check, 2026-11-21: the full abstracts of Katz 1949, Hamilton 1988, Fayet et al. 1985, Broadhead & Hymers 1987 and both Hailstone 2004 papers were read from publisher-deposited metadata (Crossref) and OpenAlex. No paper body and no kinetic constant could be read. The Deep Blue copy of Katz sits behind a Cloudflare challenge, which was not bypassed.)

---

## 1. The nucleation-and-growth model of latent-image formation

### 1.1 Sources and what was actually accessible

| Source | What it is | Access |
|---|---|---|
| Gurney, R.W. & Mott, N.F., Proc. R. Soc. A 164, 151 (1938) | Original two-stage mechanism: an electronic step (photoelectron trapped at a sensitivity speck), then an ionic step (an interstitial Ag⁺ migrates to the trapped electron). | recalled |
| [Anonymous in snippet; probably J.H. Webb], "Latent Image Theory and its Experimental Application to Motion Picture Sound-film Emulsion", J. SMPE 32(1) (1939) | Earliest quantitative reciprocity model. The abstract attributes high-intensity failure to the minimum time needed to attract a silver ion. It attributes low-intensity failure to a speck holding one silver atom being unstable, with "activated" grains decaying exponentially in time. | secondary-source (abstract via search) |
| Katz, E., "On the Photographic Reciprocity Law Failure and Related Effects. I. The Low Intensity Failure", J. Chem. Phys. 17(11), 1132 (1949). [Deep Blue](https://deepblue.lib.umich.edu/bitstream/handle/2027.42/71185/JCPSA6-17-11-1132-1.pdf;sequence=2) | Quantitative low-intensity-failure (LIRF) theory. A monatomic speck is "unstable against thermal agitation"; a diatomic speck is "stable for most practical purposes"; a triatomic speck is "completely stable". Assumption (b) of the abstract: the electron trap depth of monatomic specks "is not the same for all traps"; a quasi-continuous distribution of depths is needed to fit l.i.f. curves. DOI 10.1063/1.1747127, pp. 1132–1141. | **abstract opened** (spot-check, via Crossref publisher metadata; the PDF is behind a Cloudflare challenge) |
| Bayer, B.E. & Hamilton, J.F., "Computer Investigation of a Latent-Image Model", J. Opt. Soc. Am. 55(4), 439–452 (1965), and "Investigation of a Latent-Image Model: Recombination of Trapped Electrons and Free Holes", JOSA 55, 528 (1965) | First Monte-Carlo simulation of the nucleation-and-growth model. | bibliographic record only (secondary-source) |
| Hamilton, J.F., "The silver halide photographic process", Adv. Phys. 37(4), 359–441 (1988), DOI 10.1080/00018738800101399 | The standard quantitative review. | abstract opened (spot-check, OpenAlex); it states the nucleation-and-growth mechanism but gives no constants |
| Hailstone, R.K., J. Phys. Chem. 99, 4414–4428 (1995), "Computer simulation studies of silver cluster formation on AgBr microcrystals" | Modern Monte-Carlo nucleation-and-growth model with recombination. Parameters: electron and hole diffusion coefficients, recombination radius, time for an interstitial Ag⁺ to be captured by a trapped electron, and Ag-atom lifetime. | secondary-source (abstract via search; parameter values not seen) |
| Hailstone, R.K. & De Keyzer, R., "The effect of trap depth on latent-image formation in an AgBr tabular grain: a computer simulation study", J. Imaging Sci. Technol. 45(4), 388–392 (2001) | | record only (403) |
| Hailstone, R. & DeKeyzer, R., "Latent-image formation in tabular AgBr grains: simulation studies", Imaging Sci. J. 52(3), 164–175 (2004), DOI 10.1179/136821904225011645. [RIT record](https://repository.rit.edu/article/841) | | **abstract opened** (the PDF returned 403; the full abstract was re-read on OpenAlex in the spot-check) |
| Hailstone, R.K., French, J. & De Keyzer, R., "Latent-image formation in tabular AgBr grains: experimental studies", Imaging Sci. J. 52(3), 151–163 (2004), DOI 10.1179/136821904225011654 | Five tabular AgBr emulsions, reciprocity data from 10⁻⁴ to 10³ s. Low-irradiance reciprocity failure "commenced at 0.1–1 s". Dye-present sensitized emulsions showed little HIRF, "suggesting the minimum developable size of the latent image was three atoms for the development conditions used". Sensitizer-centre electron trap depths: 0.1 eV maximum (550 nm centre), 0.225–0.425 eV (650 nm centre), 0.45–0.65 eV (750 nm centre). | **abstract opened** (spot-check, full abstract via OpenAlex) |
| Fayet, P., Granzer, F., Hegenbart, G., Moisar, E., Pischel, B. & Wöste, L., "Latent-image generation by deposition of monodisperse silver clusters", Phys. Rev. Lett. 55(27), 3002–3004 (1985), DOI 10.1103/PhysRevLett.55.3002 | Size-selected Ag₁⁺, Ag₃⁺, Ag₄⁺, Ag₅⁺, Ag₇⁺, Ag₉⁺ deposited on binder-free AgBr microcrystals. Development occurred with Ag₄⁺ or larger, and not with smaller clusters. The critical size "may depend on the development conditions (e.g., the redox potential)". Ag₂⁺ was not tested. | **abstract opened** (spot-check, APS abstract page and OpenAlex) |
| Broadhead, P. & Hymers, C.E., "The Quantum Sensitivity of an Emulsion and the Size of the Grains Producing the Dye Image in a Coupler Starved Format", J. Photogr. Sci. 35(1), 2–9 (1987), DOI 10.1080/00223638.1987.11738425 | Residual-size analysis of a polydisperse **colour** emulsion binned into 12 size classes (√2 steps in projected area). | **abstract opened** (spot-check, OpenAlex) |
| Hanson, R.C. (Ph.D. thesis, Univ. of Illinois, 1960), "Hall Mobility of Holes in Silver-Bromide"; journal version Hanson, R.C. & Brown, F.C., J. Appl. Phys. 31(1), 210–211 (1960), DOI 10.1063/1.1735400 | Hole Hall mobility in AgBr. | thesis abstract opened (spot-check, OpenAlex); the J. Appl. Phys. note was not opened |
| Tani, T., "Physics of the photographic latent image", Physics Today 42(9) (Sept 1989). [link](https://physicstoday.aip.org/features/physics-of-the-photographic-latent-image) | | **intro opened** |
| James, T.H. (ed.), *The Theory of the Photographic Process*, 4th ed., Macmillan, New York (1977), 714 pp., ISBN 0023601905 | Ch. 4–5 not accessible; no full text online. | not opened |
| Tani, T., *Photographic Sensitivity: Theory and Mechanisms*, Oxford Univ. Press (1995) | | not opened |

### 1.2 Model structure (for the Markov-chain implementation)

The *structure* below is the consensus Gurney–Mott / Hamilton nucleation-and-growth picture as summarised in Bayer & Hamilton (1965), Hamilton (1988) and Hailstone (1995). The step list matches the Hailstone 1995 abstract (opened via snippet), but the detailed equations are **recalled**. Confirm against Hamilton 1988 §§4–6 and Hailstone 1995.

Per grain, events after a photon is absorbed in the grain or by an adsorbed dye:

| # | Event | Notation | Typical time scale | Literature status |
|---|---|---|---|---|
| E1 | Photon absorbed → free electron e⁻ + free hole h⁺. Dye excitation injects e⁻; the dye hole may or may not return. | rate Φ_g = σ_g(λ)·q(λ) [photons s⁻¹] | — | σ_g: see §3 |
| E2 | e⁻ trapped, shallow or deep. With sulfur (+gold) sensitization the deep traps are the Ag₂S / AuAgS centres, mostly at grain edges and corners in tabular grains. | trap depth E_t; detrapping rate ν₀·exp(−E_t/kT) | trapping ~ns (recalled) | E_t is fitted in simulations: 0.2–0.6 eV edge traps (Hailstone & De Keyzer 2004, opened abstract). It was also deduced experimentally for S(+Se)+Au sensitizer centres: ≤ 0.1, 0.225–0.425 and 0.45–0.65 eV (Hailstone, French & De Keyzer 2004, opened abstract). |
| E3 | Ionic step: an interstitial Ag_i⁺ reaches the trapped e⁻ → Ag⁰ atom | τ_i | µs order (recalled) | from ionic conductivity / interstitial concentration (measured quantities); used as an input in the simulations |
| E4 | Single Ag⁰ is thermally unstable: Ag⁰ → Ag_i⁺ + e⁻ | lifetime τ₀ | order-of-magnitude prior ~0.05–5 s at room temperature (estimate, §1.4) | inferred and **fitted** from LIRF; never directly measured for photolytic Ag⁰ on S+Au emulsions in anything I could open |
| E5 | Nucleation: a second e⁻ is trapped at Ag⁰ (Ag⁰ is a deeper trap), followed by Ag_i⁺ → Ag₂. Ag₂ is the **minimum stable size**. | — | — | consensus (Katz 1949 abstract; Hailstone 1995 abstract: "trapping of an electron at a Ag atom, followed by the migration of an interstitial Ag ion to yield a stable two-atom cluster") |
| E6 | Growth: Ag_n + e⁻ + Ag_i⁺ → Ag_{n+1}. Capture cross-section rises with n. | — | — | structure |
| E7 | Developability: a grain is developable if any cluster has n ≥ n_dev | n_dev | — | n_dev = 3–4 (Tani 1989, opened). Fayet et al. 1985 (abstract opened in the spot-check): deposited Ag₄⁺ and larger developed, Ag₃⁺ did not, and the critical size "may depend on the development conditions (e.g., the redox potential)". Hailstone, French & De Keyzer 2004 (abstract opened) inferred 3 atoms for their development conditions. That weaker developers need larger clusters remains recalled (Tani 1995; James 1977). |
| E8 | Recombination: free h⁺ + trapped e⁻ → lost. Also "hole attack": h⁺ + Ag_n → Ag_{n−1} + Ag_i⁺. | recombination radius r_c, hole diffusion D_h | — | Bayer & Hamilton 1965b treat trapped-e⁻ / free-h⁺ recombination explicitly. Hailstone & De Keyzer 2001 found enhanced recombination at corner traps (secondary). |
| E9 | Hole removal: h⁺ trapped by a gold centre, by gelatin / halogen acceptors, or escapes. This reduces E8. | — | — | qualitative; fitted |
| E10 | Dispersity (multiple competing nuclei): if more than one Ag₂ forms on separate sites, the photons are split between sub-developable clusters. | number of sites S per grain | — | Explains high-intensity reciprocity failure (HIRF) and the drop in quantum sensitivity of large grains (1969 J. Photogr. Sci. abstract, secondary: "latent-image silver is dispersed over a number of sites"). |

**How reciprocity failure arises in this structure**

- **LIRF (long exposures).** When the mean interval between electron arrivals at the speck, 1/(η·Φ_g), becomes comparable to τ₀, a single Ag⁰ decays (E4) before the second electron arrives (E5). Nucleation efficiency falls, so more photons are needed. Growth beyond Ag₂ is unaffected. LIRF therefore lowers speed and **raises contrast**. This is consistent with Kodak's instruction to *reduce* development by 10–30 % for long exposures (Tri-X table, §4).
- **HIRF (short exposures).** When many photons arrive within τ_i, several electrons are trapped at different sites before any ionic step completes. Several Ag⁰/Ag₂ nuclei form at once (dispersity, E10), and the high instantaneous hole density raises recombination (E8). The result is many sub-developable clusters, lower speed and **lower contrast**. This is consistent with Kodak's instruction to *increase* development by 10–20 % for short exposures (Tri-X table, §4). The 1939 SMPE abstract attributes HIRF to "the minimum time required for the attraction of a silver ion" (E3).
- The reciprocity-optimum exposure time lies between the two regimes. Hailstone & De Keyzer (2004, opened abstract) report a best simulated efficiency of 6–8 absorbed photons/grain for a 0.01 s exposure. The grain-diameter and thickness effects on efficiency that they saw at 0.01 s disappeared in simulations with a 10⁻⁶ s exposure, which they take to show that the efficiency loss at 0.01 s came from differences in the onset of low-intensity reciprocity failure. In other words, in their simulated tabular grains LIRF already begins near 0.01 s; the abstract does not state a reciprocity-optimum time. (Kodak's tables show no LIRF for either film until after 1/10 s, so the D1 constants must place LIRF onset later than in those simulations. The companion experimental paper by Hailstone, French & De Keyzer (2004; abstract opened in the spot-check) measured LIRF onset at 0.1–1 s on real tabular emulsions. That agrees with the Kodak tables, so the early onset at 0.01 s is a feature of the simulation parameters, not of tabular grains.)

### 1.3 Equations to implement

All equations are **recalled / estimate**. They are standard kinetic-theory forms consistent with the structure above, and their constants are D1 fit parameters.

**(a) Reciprocity-obeying limit: the quantum-hit (Poisson) model.** No decay, no dispersity. For a size class with mean absorbed photons m = σ_g·H_q (H_q = photon exposure, photons/area) and per-photon efficiency η:

  P_dev(m) = 1 − Σ_{k=0}^{r−1} e^{−ηm} (ηm)^k / k!

where r is the effective number of "useful" electrons required (≥ n_dev). This gives the characteristic-curve toe per size class and is the right unit test for the Markov chain: with τ₀→∞ and τ_i→0 the chain must reduce to this.

**(b) Single-site Markov chain (concentration principle).** One speck per grain, state k = 0, 1, …, n_dev (n_dev absorbing). For constant photon rate Φ (per grain), with η = probability that an absorbed photon's electron reaches the speck without recombining:

| Transition | Rate |
|---|---|
| 0 → 1 | η₁·Φ |
| 1 → 0 | 1/τ₀ + ρ₁·Φ (thermal decay + hole attack) |
| 1 → 2 | η₂·Φ |
| k → k+1 (k ≥ 2) | η_g·Φ |
| k → k−1 (k ≥ 2) | ρ_g·Φ (small with gold / hole traps) |

The master equation is dP/dt = Q(Φ)·P with tridiagonal Q. Constant irradiance gives P(t) = exp(Q t)·P(0). Precompute P_dev(Φ, t) per size class on a log grid; this is the lookup table in the design. The model is linear in Φ apart from the 1/τ₀ term, so LIRF appears naturally. In the chain-of-two-steps limit (Katz-type), the nucleation probability per first atom is η₂Φ/(η₂Φ + 1/τ₀). A reciprocity curve at fixed density then rises as 1 + 1/(η₂Φτ₀) at low Φ. Katz (1949, abstract opened in the spot-check) found that a single trap depth does not reproduce measured LIRF curves. A quasi-continuous distribution of monatomic-speck trap depths (that is, of τ₀) is required. So τ₀ should be implemented as a distribution (for example log-normal in τ₀, i.e. normal in E_t) from the start.

**(c) HIRF extension.** Either:

- (i) a multi-site version, S sites per grain. Each electron goes to a site chosen ∝ trap strength. The ionic step is a separate state with rate 1/τ_i, so electrons trapped while the ionic step is pending can nucleate elsewhere. This is the Bayer–Hamilton approach.
- (ii) the cheaper option: make η₁, η₂ functions of Φτ_i, for example η(Φ) = η₀/(1 + Φτ_i·(S−1)/S). This is an estimate-level closure, not from a source.

Option (i) is what the literature uses. Option (ii) must be marked `assumed`.

**(d) Grain-size dependence.**

- Absorbed photons scale with σ_g. For intrinsic absorption this is ∝ grain volume where α·d ≪ 1 (weakly absorbed light, thin grains), tending to ∝ projected area only where α·d ≳ 1 (strongly absorbed UV / short blue); for dye absorption it is ∝ surface area (see §3.1).
- The *requirement* also grows with size, because of more competing sites and longer electron paths. This is shown by:
  - the 1969 J. Photogr. Sci. abstract (secondary): "number of quanta per grain required for developability steadily increases with size within the normal size range";
  - Farnell, J. Photogr. Sci. 17, 116–125 (1969), as cited in Kofron US 4,439,520 (opened patent text): for unsensitized AgBr/AgBrI, blue speed stops rising above about 0.5 µm² projected area (about 0.8 µm diameter);
  - Broadhead & Hymers, J. Photogr. Sci. 35(1), 2–9 (1987) (abstract opened in the spot-check; an earlier draft wrongly gave the authors as "Hailstone et al."): the most sensitive grains were among the largest, needed on average 12 absorbed photons for 50 % developability, and quantum sensitivity fell off with decreasing grain size. In that polydisperse **colour** emulsion small grains were *less* efficient.
- So the dependence is not monotonic across emulsions. For D1, make η and S functions of size class with 1–2 fitted parameters, not a fixed law.

### 1.4 Parameter values found

| Parameter | Value | Source | Confidence | Measured / fitted in literature |
|---|---|---|---|---|
| Minimum stable cluster | Ag₂ (Katz: diatomic "stable for most practical purposes", triatomic "completely stable") | Katz 1949, DOI 10.1063/1.1747127 (abstract); Hailstone 1995 (abstract via search snippet only) | opened-source (abstract, spot-check) | inferred from reciprocity and two-flash experiments |
| Minimum developable cluster | 3–4 atoms ("as few as three or four silver atoms") | Tani, Physics Today 42(9), 1989, DOI 10.1063/1.881181 | opened-source | inferred. Supported by cluster deposition: Fayet et al., Phys. Rev. Lett. 55(27), 3002 (1985), abstract opened in the spot-check, found deposited Ag₄⁺ and larger developable, Ag₃⁺ and Ag₁⁺ not (Ag₂⁺ not tested), with a size that may depend on redox potential. Hailstone, French & De Keyzer 2004 (abstract) inferred 3 atoms under their development conditions. |
| Development gain | > 10⁹ | Tani 1989 | opened-source | — |
| Ag ions in a 1.5 µm × 0.08 µm tabular grain | 2.9 × 10⁹ | computed (ρ_AgBr = 6.473 g cm⁻³, recalled CRC value) | estimate | — |
| Quantum sensitivity, best simulated (tabular 1.2 × 0.2 µm, edge traps) | 6–8 absorbed photons/grain at 0.01 s | Hailstone & De Keyzer, Imaging Sci. J. 52, 164 (2004) | opened-source (abstract) | simulation (fitted traps) |
| Quantum sensitivity, measured, 5 tabular AgBr emulsions | 10–19 absorbed photons/grain (400 nm); 13–27 absorbed photons/grain (spectral, dye-sensitized); midpoint of D–logE, 0.1 s exposures; the authors computed these "neglecting the polydisperse nature of these emulsions" | Hailstone, French & De Keyzer, Imaging Sci. J. 52(3), 151–163 (2004), DOI 10.1179/136821904225011654 | opened-source (full abstract via OpenAlex, spot-check) | measured |
| Quantum sensitivity, residual-size analysis | 12 absorbed photons/grain for 50 % developable (largest, most sensitive grains; colour emulsion, 12 size classes) | Broadhead, P. & Hymers, C.E., J. Photogr. Sci. 35(1), 2–9 (1987), DOI 10.1080/00223638.1987.11738425 | opened-source (abstract via OpenAlex, spot-check) | measured |
| Most sensitive grains of "modern" emulsions | ≥ 4 absorbed photons | Farnell & Chanter, "The Quantum Sensitivity of Photographic Emulsion Grains", J. Photogr. Sci. 9, 73–83 (1961), as quoted in a patent | secondary-source | measured |
| Edge trap depth (S+Au sensitization analogue) | 0.2–0.6 eV | Hailstone & De Keyzer 2004 | opened-source (abstract) | fitted |
| Ag-atom lifetime from Ag₂ donor centres (reduction-sensitized) | ~1 s at room temperature | Imaging Sci. J. (1997), DOI 10.1080/13682199.1997.11736362 | secondary-source | **Not the photolytic τ₀ of E4.** Do not use as τ₀; at most an order-of-magnitude prior. |
| τ₀ (photolytic Ag⁰, S+Au emulsion) | not found. Order-of-magnitude prior 0.05–5 s (Tri-X ~0.05–1 s, T-Max 100 ~0.2–4 s). | see derivation below | estimate | fitted in all simulations seen |
| τ_i (ionic step) | not found; µs order | recalled (Hamilton 1988) | recalled; **UNVERIFIED** (the Hamilton 1988 abstract, opened in the spot-check, gives no time constant) | input from ionic-conductivity data |
| Hole mobility in AgBr, room temperature | **2.0 ± 0.5 cm² V⁻¹ s⁻¹** (Hall, 27 °C), falling to 0.5 ± 0.15 at 150 °C, about T⁻⁴. The earlier "~1" was **wrong** for the Hall value. A room-temperature hole *drift* mobility was not found; it may differ from the Hall value. | Hanson, Ph.D. thesis, Univ. of Illinois (1960), abstract; journal version Hanson & Brown, J. Appl. Phys. 31, 210 (1960) | opened-source (thesis abstract, spot-check) | measured |
| Electron mobility in AgBr, room temperature | ~60 cm² V⁻¹ s⁻¹ | Hanson (1960) abstract: the hole mobility "is about 30 times smaller than the electron mobility determined by other workers at room temperature". ~60 is 30 × 2.0, computed here; the original electron measurement was not opened. Burnham, Brown & Knox, Phys. Rev. 119, 1560 (1960) only cover 4–125 K (150 cm² V⁻¹ s⁻¹ at 125 K). | estimate (from a ratio in an opened abstract); **UNVERIFIED as a direct measurement** | measured (Hall effect) |
| Recombination radius, trap densities | not found | Hailstone 1995 lists them as model parameters | — | fitted |

**Derivation of the τ₀ prior (estimate).** The Tri-X table needs +1 stop at an indicated 1 s; T-Max 100 needs +1/3 stop at 1 s. With roughly 10–20 absorbed photons per grain near mid-curve (the quantum-sensitivity rows above), an indicated 1 s exposure means a per-grain photon rate Φ ≈ 10–20 s⁻¹, i.e. photons every 0.05–0.1 s. In the two-step form of §1.3(b), the exposure multiplier is 1 + 1/x with x = η₂Φτ₀.
- Tri-X, ×2 at 1 s → x = 1 → τ₀ = 1/(η₂Φ) ≈ (0.05–0.1 s)/η₂, i.e. about 0.05–1 s for η₂ = 0.1–1. Check at 1/10 s: x = 10 → ×1.1 (0.14 stop), consistent with "none" in the table.
- T-Max 100, ×2^(1/3) = 1.26 at 1 s → x ≈ 3.8 → τ₀ ≈ (0.19–0.38 s)/η₂, i.e. about 0.2–4 s (or a deeper first trap).

Caveat: with one fixed τ₀ this simple form is too steep. Fixing x = 1 at 1 s gives ×11 at an indicated 10 s and ×101 at 100 s, against Tri-X table values of ×5 (50 s) and ×12 (1200 s); its asymptotic Schwarzschild exponent is 0.5, against about 0.72 derived in §4.3. The full chain (growth competition, a spread of trap depths or τ₀ across sites/size classes) has to flatten this. Katz (1949, abstract opened in the spot-check) reached the same conclusion: a quasi-continuous distribution of electron trap depths at monatomic specks was required to reproduce measured LIRF curves, and the l.i.f. slope varied systematically with grain size. Treat all of this as a prior range for fitting only.

**Measured vs fitted, summary.** Measured in the literature: quantum sensitivity (photons per grain), grain absorption, electron/hole mobilities, ionic conductivity, reciprocity curves, and sensitizer-centre trap depths (from long-wavelength sensitivity, Hailstone, French & De Keyzer 2004). Inferred from experiments: n_stable = 2 and n_dev ≈ 3–4 (Ag₄⁺ is the smallest developable size in cluster deposition). Fitted in every simulation I could identify: trap depths and densities, τ₀, recombination probability, and number of sites.

---

## 2. Grain population data

### 2.1 Kodak Tri-X 400

| Item | Value | Source | Confidence |
|---|---|---|---|
| Grain size, distribution, shape | **not found** in any Kodak publication or accessible paper | — | — |
| Emulsion type | conventional (three-dimensional), polydisperse silver bromoiodide | general knowledge; the datasheet says nothing on emulsion | recalled |
| Kodak's own "conventional" AgBrI comparison emulsions (not Tri-X) | Control 5: 4.7 mol% I, ECD 1.4 µm, t = 0.42 µm, aspect ratio 3.3:1. Control 6: 10 mol% I, 1.1 µm, ≅0.40 µm, 2.8:1. Control 7: 5 mol% I, 1.0 µm, ≅0.40 µm, 2.5:1. | Kofron et al., US 4,439,520, table comparing Examples 1–4 with Controls 5–7 ([Google Patents](https://patents.google.com/patent/US4439520A/en)) | opened-source |
| Kodak low-aspect-ratio control "Emulsion A" | 0.8 µm average grain size, aspect ratio < 3:1, AgBrI 1 mol% I, double-jet. Precipitated with 0.12 mg/Ag mol ammonium hexachlororhodate(III); a Rh dopant typically alters reciprocity and contrast. Chemically sensitized with 4.4 mg/Ag mol sodium thiosulfate pentahydrate + 1.75 mg/Ag mol potassium tetrachloroaurate, 23 min at 60 °C. Dye 87 mg/Ag mol. Coated at 1.75 g/m² Ag, 4.84 g/m² gelatin **on a TiO₂–gelatin layer on paper support** (not film). | US 4,439,520 | opened-source (verbatim) |
| Diffuse RMS granularity | 17 (TRI-X 400); 16 (TRI-X 320). 48 µm aperture, net diffuse D = 1.0, HC-110 dil. B, 20 °C. | Kodak F-4017, Oct 2021, p. 7 | opened-source |
| Generic textbook range for conventional camera-speed negative emulsions | ECD about 0.2–2 µm, lognormal with σ_ln about 0.3–0.6, iodide about 2–10 mol% | recalled (James 1977, ch. 3) | recalled; **UNVERIFIED**; confirm in James ch. 3 or Tani 1995 ch. 1 |

### 2.2 Kodak T-Max 100

| Item | Value | Source | Confidence |
|---|---|---|---|
| ECD, thickness, aspect ratio, size distribution | **not found**. Kodak has not published them for T-Max 100. | — | — |
| Emulsion type | T-GRAIN (tabular) | Kodak F-4016 | opened-source |
| Diffuse RMS granularity | 8 (D-76, 20 °C; 48 µm aperture, net diffuse D = 1.00, 12× magnification) | Kodak F-4016, Feb 2016, p. 8; re-checked on the rendered p. 8 of the June 2018 F-4016 ("Image Structure: … based on development in KODAK Developer D-76, at 68°F (20°C)") | opened-source |
| Resolving power | 63 lines/mm at TOC 1.6:1; 200 lines/mm at TOC 1000:1 | F-4016 (2016) p. 8 | opened-source |
| 120 base | 4.7-mil (≈119 µm) | F-4016 (Jan 2002) p. 2 | opened-source |
| T-Max 400 cross-section ("tabular grains over a layer of conventional grains") | — | forum description of the Kodak book *Making KODAK Film*; T-Max 400 only | forum, pointer only |

**Tabular-grain definitions and example dimensions.** Kofron, J.T., Booms, R.E., Jones, C.G., Haefner, J.A., Wilgus, H.S., Evans, F.J., "Sensitized high aspect ratio silver halide emulsions and photographic elements", US 4,439,520 (filed 30 Sep 1982, priority 12 Nov 1981, granted 27 Mar 1984; Eastman Kodak).

Definitions:
- thickness < 0.3 µm (0.5 µm permitted for blue-sensitive emulsions);
- diameter ≥ 0.6 µm;
- average aspect ratio > 8:1, preferred ≥ 12:1, optimally ≥ 20:1;
- tabular grains ≥ 50 % of projected area, preferred ≥ 70 %, optimally ≥ 90 %;
- iodide preferably ≥ 0.1 mol%, up to about 15–20 mol%.

| Patent emulsion | I (mol%) | ECD (µm) | t (µm) | AR | % proj. area | Coating | Confidence |
|---|---|---|---|---|---|---|---|
| Example 1 | 6 | ≅3.8 | 0.14 | 27:1 | >50 | — | opened (HTML table) |
| Example 2 | 1.2 | ≅3.8 | 0.14 | 27:1 | 75 | — | opened |
| Example 3 | 12.0 | 2.8 | 0.15 | 19:1 | >90 | — | opened |
| Example 4 | 12.3 | 1.8 | 0.12 | 15:1 | >50 | — | opened |
| Emulsion B (Group VIII noble-metal-doped example) | 1 | 1.5 | 0.08 | 19:1 | 90 | coated as Emulsion A: 1.75 g/m² Ag, 4.84 g/m² gel, on paper. S+Au: 5 mg/Ag mol sodium thiosulfate pentahydrate + 5 mg/Ag mol potassium tetrachloroaurate, 30 min at 65 °C. 150 mg dye per Ag mol. Exposed at 10⁻⁴ s (EG&G sensitometer). | opened (verbatim) |
| AgBr series (no iodide) | 0 | 1.43–2.96 | 0.07–0.12 | 16.7–37:1 | >95 | — | opened via fetch summary; the 2.96 µm / 0.08 µm / 37:1 row was verified verbatim |
| AgBrI 3 mol% example (blue) | 3 | 1.0 | 0.10 | ~10:1 | >85 | 1.08 g/m² Ag, 2.9 g/m² gel | fetch summary only (not verbatim-checked) |
| Other single-layer coatings in the patent | — | — | — | — | — | 1.07, 1.35, 2.15 g/m² Ag | fetch summary only |

Use for D1: these are Kodak's own tabular emulsions from the T-GRAIN programme, not T-Max 100 itself. For T-Max 100 the size class grid has to be **fitted**: ECD and thickness enter through σ_g and the Nutting covering power, and the fit is constrained by RMS granularity 8 and the D-76 curve. A defensible prior is ECD 0.6–1.5 µm and t 0.08–0.15 µm, but this is **estimate** only; no source gives it for T-Max 100. Note that the small-ECD / thick corner of this box (0.6 µm / 0.15 µm, AR 4) falls below the patent's > 8:1 tabular definition, so the prior should be sampled jointly with AR ≥ 8 if it is meant to stay "tabular".

### 2.3 Coated silver and layer thickness

| Item | Value | Source | Confidence |
|---|---|---|---|
| Silver in pictorial B&W films (Tri-X, T-Max) | **not found** in any primary source. Kodak J-210 says per-product values were available only from Kodak Environmental Services on request. | Kodak J-210, "Sources of Silver in Photographic Processing Facilities" (1998), pp. 2–3 | opened-source (the absence) |
| Kodak colour negative films (context) | GOLD 100/200: 16; GOLD 400: 24; ROYAL GOLD 100/200/400: 17/18/22 troy oz / 1000 ft². That is 5.4 / 8.0 / 5.7 / 6.0 / 7.4 g/m². | J-210, p. 4 (×0.3348 g m⁻² per troy oz/1000 ft²) | opened-source; conversion is an estimate |
| Kodak graphic-arts films (context) | Scanner 2000 SAI and Camera 2000 CGP: 9.6 troy oz / 1000 ft² = 3.2 g/m² | J-210, p. 5 | opened-source |
| Kodak colour paper (context) | Ektacolor Royal III: 2 troy oz / 1000 ft² = 0.67 g/m² | J-210, p. 4 | opened-source |
| "Typical B&W film 3–5 g/m²" | forum statement only; never cite as data | Photrio thread "Silver content of BW film in g/cm2" | forum (no better source found) |
| Single-layer tabular coatings in the Kofron patent | 1.07–2.15 g/m² Ag; gelatin 2.15–10.4 g/m² | US 4,439,520 | opened-source (1.75/4.84 verified verbatim) |
| Emulsion-layer thickness for Tri-X / T-Max | **not found** | — | — |
| Dry layer thickness from coverage (formula) | t_dry [µm] = m_Ag·(187.77/107.87)/ρ_AgBr + m_gel/ρ_gel, with ρ_AgBr = 6.473 g cm⁻³ and ρ_gel ≈ 1.35 g cm⁻³ (both recalled), m in g m⁻². For 1.75 g/m² Ag + 4.84 g/m² gel: 0.47 + 3.59 = **4.1 µm**. | computed | estimate |
| Grain number density, Emulsion B | V = 0.141 µm³, 5.26 × 10⁻¹³ g Ag per grain → 3.33 grains/µm² of film at 1.75 g/m² Ag. Projected area 1.77 µm² per grain, so about 5.9 grain layers deep. | computed | estimate |

---

## 3. Light absorption by grains and the datasheet spectral sensitivity curve

### 3.1 Intrinsic AgBr / AgBrI absorption

- The AgBr edge is indirect (phonon-assisted), with an exponential Urbach tail. Absorption is weak beyond about 480–500 nm, and iodide extends it by tens of nm.
- Absolute absorption spectra: Moser, F. & Urbach, F., "Optical Absorption of Pure Silver Halides", Phys. Rev. 102, 1519–1523 (1956), accurate to better than 10 % at room temperature (abstract, secondary). **Absorption coefficients α(λ) in cm⁻¹ were not obtained.** Open Moser & Urbach 1956 or Palik's handbook.
- Surface / size effect: 65 Å-radius AgBr has about half the bulk α (Phys. Rev. B, 1967; abstract, secondary). This is negligible for camera-film grain sizes.
- For grains with α·d ≪ 1, σ_g,intrinsic(λ) ≈ α(λ)·V_g. For thicker grains, use 1 − exp(−α·d_eff) times the projected area (estimate, standard optics). Tabular grains at the same volume have a larger projected area but absorb about the same intrinsic light. They gain on dye absorption, because a dye monolayer covers the larger surface area. That is the T-GRAIN rationale.

### 3.2 Spectral sensitising dyes

- Dyes saturate at about one adsorbed monolayer. Fuji patents take a molecular occupancy area of **80 Å²** per dye molecule to define the saturation coverage, and use multilayer adsorption to exceed it. The patent defines "light absorption strength" as ∫ OD dν over 5,000–35,000 cm⁻¹, per unit grain surface, with ≥ 100 as the multilayer target (Suzumoto, Urabe, Yamashita, US 6,048,681, Fuji Photo Film; opened via Google Patents). It gives no numeric monolayer absorption fraction.
- **Per-molecule cross-section (standard conversion):** σ = (1000·ln10/N_A)·ε = 3.82 × 10⁻²¹ cm² × ε [M⁻¹ cm⁻¹].
- **Monolayer absorption estimate:** for J-aggregate peak ε ~ 2 × 10⁵ M⁻¹ cm⁻¹ (recalled typical value), σ ≈ 7.6 Å². Divided by 80 Å² per molecule, that gives ≈ 10 % absorbed per covered face at the band peak. A tabular grain has two faces, so ≈ 15–20 % of light on its projected area is absorbed at peak. Treat this as a prior only.
- Efficiency: the relative quantum yield of spectral sensitization φ_r (electrons injected per photon absorbed by the dye, relative to intrinsic absorption at 400 nm) is the standard measure (Tani). The measured-vs-spectral numbers in Hailstone et al. 2004 (10–19 photons at 400 nm vs 13–27 photons for spectral exposure) imply an effective φ_r of about 0.7 for those emulsions (estimate; the abstract was opened in the spot-check). Caveats:
  - The two ranges are not paired emulsion by emulsion.
  - The spectral exposures were at the dye band, not at a matched absorbed fraction.
  - The ratio of photons per grain is only a crude proxy for φ_r.
  - Tani (1980), as summarised in Kofron US 4,439,520 (opened), reports that the relative quantum yield of spectral sensitization "has been observed to be near unity".
  - Use a prior φ_r ≈ 0.7–1.0, not a point value of 0.7.

### 3.3 From the datasheet curve to a per-grain absorbed-photon rate

What the curves are (both **opened-source**, rendered pages):

| Film | Datasheet | Definition | Effective exposure | Process | Criteria | Read-off level (eyeballed from the plot; to be digitised) |
|---|---|---|---|---|---|---|
| Tri-X 400 | F-4017, Oct 2021, p. 8 | log S, S = 1/(exposure in erg/cm² for the specified density) | 0.5 s | not stated on the graph | D = 0.3 and D = 1.0 above gross fog; densitometry diffuse visual | solid curve (legend "D=1.0") ≈ 2.2–2.7; dashed (legend "D=0.3") ≈ 1.2–1.6; cut-off about 650–670 nm |
| T-Max 100 | F-4016, Feb 2016, p. 8 | same | 1.4 s | D-76, 20 °C | 0.3 and 1.0 above D-min | 0.3-curve ≈ 1.0–1.5; 1.0-curve ≈ 0.0–0.5; cut-off about 650–700 nm |

**Warning: the Tri-X legend appears inverted** (estimate, from a consistency check).
- In F-4017 the curve labelled D = 1.0 lies *above* the one labelled D = 0.3. That is impossible, because higher density needs more exposure and therefore has lower sensitivity.
- An ISO cross-check:
  - ISO 400 means H_m = 0.8/400 = 0.002 lx·s at D_min + 0.1.
  - 1 lx·s at 555 nm = 1.464 erg/cm².
  - That predicts Tri-X log S(D = 0.3) ≈ 1.7–2.3 and log S(D = 1.0) ≈ 0.5–1.1.
  - The same check reproduces the T-Max 100 levels well.
- The Tri-X solid curve (2.2–2.7) therefore behaves like the D = 0.3 curve, and matches the upper end of the prediction. Even with the legend swapped, the dashed curve (1.2–1.6) sits about 0.5 log above the predicted D = 1.0 level, so the absolute level does not fit ISO 400 either way. (I re-checked the legend on the rendered p. 8 of the Oct 2021 edition: the solid line is labelled "D=1.0>Gross Fog".)
- The May 2007 F-4017 (p. 10) states that the MTF and spectral sensitivity data were generated with the older KODAK TRI-X Pan Film 5063 (per the sibling report `datasheets.md`, §2.1), which may account for part of the mismatch.
- Spot-check (2026-11-21): the rendered p. 8 of the Oct 2021 F-4017 was viewed again. The solid curve is labelled "D=1.0>Gross Fog" and lies at about 2.2–2.7, above the dashed "D=0.3>Gross Fog" curve at about 1.15–1.6; the effective exposure is 0.5 s. The rendered p. 8 of the June 2018 F-4016 shows the T-Max 100 curves at 1.4 s, D-76 20 °C, "0.3 greater than D-min" at about 1.0–1.5 and "1.0 greater than D-min" at about −0.2 to 0.5, in the expected order. The T-Max 100 levels in the table above hold for that edition.
- Use only the **relative spectral shape** from F-4017. Do not use its absolute level or the legend for calibration without flagging the problem.

**Procedure** (estimate-level method; standard radiometry):

1. **Photon exposure per erg.** 1 erg/cm² at wavelength λ equals 5.034 × λ[nm] photons/µm². That is 2,014 at 400 nm, 2,769 at 550 nm and 3,272 at 650 nm. At 555 nm, 1 lx·s = 4.09 × 10³ photons/µm². The ISO 400 speed point (0.002 lx·s) is therefore about 8 photons/µm² incident (555 nm equivalent), and ISO 100 (0.008 lx·s) about 33 photons/µm².
2. **Quantum spectral sensitivity.** The datasheet S(λ) is per unit energy. Per incident photon, S_q(λ) ∝ S(λ)/λ.
3. **Assumption.** The per-absorbed-photon latent-image efficiency is independent of λ (φ_r constant across the dye band, relative to intrinsic absorption). The reference exposure times (0.5 s / 1.4 s) are **not** inside the reciprocity-normal range of the tables (no correction from 1/1000 to 1/10 s): both lie at the start of LIRF, roughly +1/2 stop for Tri-X at 0.5 s (interpolated between none at 1/10 s and +1 stop at 1 s) and about +1/3 stop for T-Max 100 at 1.4 s. In this model LIRF acts after electron generation and so is wavelength-independent; it lowers the absolute level but leaves the relative spectral shape unchanged (assumed). Then the film-level sensitivity is proportional to the population-weighted absorption: S_q(λ) ∝ Σ_classes w_c·σ_c(λ), up to optical depth effects in the layer.
4. **Parameterisation.** Write σ_c(λ) = A_c·[a_intr(λ)·(V_c/A_c) + a_dye(λ)·(S_c/A_c)], where A_c is the projected area, V_c the volume and S_c the surface area of size class c. Fit the two shape functions a_intr(λ) and a_dye(λ), the dye term with its known J-band. Constrain the shape by S_q(λ) normalised at one wavelength.
5. **Absolute scale.** Fix the scale by forward simulation: the full chain (absorption → Markov chain → development → Nutting density) must reproduce the D-76 characteristic curve and ISO speed. **Do not** take the absolute level from the spectral curve. Also check that the D = 0.3 and D = 1.0 curves have the same shape. For T-Max 100 they are offset by about 1.0–1.2 log units with a similar shape; shape differences indicate layer optical-density effects in the blue, which the layer transport of stage 1 should reproduce.
6. **Per-grain photon rate.** Φ_g,c = ∫ σ_c(λ)·E_q(λ, z) dλ, where E_q is the local photon irradiance at depth z, from stage 1. Absorbed count per exposure is Poisson with mean Φ_g,c·t. This is what feeds the Markov chain.

---

## 4. Kodak reciprocity tables

### 4.1 Tri-X 320 and 400: KODAK PROFESSIONAL TRI-X 320 and 400 Films, F-4017, October 2021

Source: "Exposure and Development Adjustments for Long and Short Exposures", p. 2. One table covers both films. [kodakprofessional.com PDF](https://kodakprofessional.com/sites/default/files/wysiwyg/film/f4017_trix_320400.pdf). Opened-source.

| Indicated exposure time (s) | Aperture adjustment | OR adjusted time (s) | AND development adjustment |
|---|---|---|---|
| 1/100,000 | +1 stop | change aperture | +20 % |
| 1/10,000 | +1/2 stop | change aperture | +15 % |
| 1/1,000 | none | none | +10 % |
| 1/100 | none | none | none |
| 1/10 | none | none | none |
| 1 | +1 stop | 2 | −10 % |
| 10 | +2 stops | 50 | −20 % |
| 100 | +3 stops | 1200 | −30 % |

The sheet adds that interpolation between 1 and 100 s is difficult and refers to graphs on the same page (calculated vs adjusted time, 1–100 s).

### 4.2 T-Max 100

F-4016, February 2016 (© 2016 Kodak Alaris), p. 2, "Adjustments for Long and Short Exposures". Identical values appear in the Jan 2002 F-4016 (T-MAX family), p. 5. Opened-source. There is no development-adjustment column. Spot-check: the Feb 2016 edition was not re-opened. The June 2018 F-4016 (© 2018 Kodak Alaris), p. 2, was opened and gives exactly the values below. Cite that edition if the Feb 2016 copy cannot be produced.

| Indicated exposure time (s) | Aperture adjustment | OR adjusted time (s) |
|---|---|---|
| 1/10,000 | +1/3 stop | change aperture |
| 1/1,000 | none | none |
| 1/100 | none | none |
| 1/10 | none | none |
| 1 | +1/3 stop | change aperture |
| 10 | +1/2 stop | 15 |
| 100 | +1 stop | 200 |

For comparison, the 2002 sheet gives T-Max 400: none from 1/10,000 to 1/10 s; 1 s +1/3 stop; 10 s +1/2 stop → 15 s; 100 s +1 1/2 stops → 300 s.

### 4.3 Schwarzschild exponents

**Kodak publishes none** in either sheet. Derived local values are below (estimate). They assume equal density for I·t_adj^p = const between table rows, so p = Δlog t_ind / Δlog t_adj.

| Film | Segment | Local p |
|---|---|---|
| Tri-X | 1 → 10 s | 0.715 |
| Tri-X | 10 → 100 s | 0.725 |
| Tri-X | fit over 1–100 s | t_adj ≈ 2.0 · t_ind^1.39 (reproduces 49 s and 1200 s) |
| Tri-X | 10⁻³ → 10⁻⁴ s and 10⁻⁴ → 10⁻⁵ s (HIRF) | **1.15** (corrected in the spot-check; was 1.18) |
| T-Max 100 | 1 → 10 s | 0.93 (0.94 if the 1 s +1/3-stop row is treated as an intensity change) |
| T-Max 100 | 10 → 100 s | 0.89 |
| T-Max 100 | 10⁻³ → 10⁻⁴ s | **1.10** (corrected in the spot-check; was 1.11) |

**Correction to the HIRF rows (spot-check).** Kodak prescribes an *aperture* change for these rows, so the exposure time is fixed and the intensity is multiplied by the stop factor f. Under I·t^p = const at fixed t, f ∝ t^(1−p), so p = 1 + log₁₀ f per decade. For Tri-X, f = 2^(1/2) per decade gives p = 1 + 0.1505 = 1.15; for T-Max 100, f = 2^(1/3) gives p = 1.10. The earlier values (1.18, 1.11) came from treating f as a lengthened time, t_adj = f·t, at unchanged intensity. That models a different exposure from the one Kodak specifies. The long-exposure rows are unaffected, because Kodak gives adjusted times there. Arithmetic re-checked: 1/log₁₀(50/2) = 0.715, 1/log₁₀(1200/50) = 0.725, 1/log₁₀(200/15) = 0.889, and 2.0·10^1.39 = 49 s, 2.0·100^1.389 = 1200 s.

These are fit or validation targets for the Markov chain, not inputs. The design says reciprocity failure must emerge from the model.

---

## 5. Consequences for D1 (summary)

1. **Circularity risk in reciprocity validation.** The latent-image constants (τ₀, τ_i, η, ρ, number of sites) are not available from opened literature and must be fitted. If they are fitted to the Kodak reciprocity tables, the "倒易律修正 vs 手册表格" check becomes circular. Proposed split, either:
   - (a) fit on the 1/1000–10 s rows of each film and hold out the 100 s rows and the 10⁻⁵ s Tri-X row; or
   - (b) fit the latent-image constants on Tri-X only and predict T-Max 100's table. Note that grain shape and trap distribution differ, so this is a hard test.

   In both cases the Tri-X development adjustments (+20 / +15 / +10 % short, −10 / −20 / −30 % long) are an independent contrast-shape check. The model must reproduce *contrast* changes with exposure time, not only speed loss.
2. **Spectral reference conditions sit at the start of LIRF.** Spectral sensitivity curves use 0.5 s (Tri-X) and 1.4 s (T-Max 100) effective exposure. Both films need no exposure correction between 1/1000 and 1/10 s, but +1 stop (Tri-X) and +1/3 stop (T-Max 100) at 1 s. So the Tri-X 0.5 s reference is roughly +1/2 stop into LIRF (interpolated) and the T-Max 1.4 s reference about +1/3 stop. The forward model should simulate the spectral curves at those exposure times, not at the 1/50–1/100 s sensitometric time.
3. **Tri-X spectral curve legend is inverted, and its absolute level does not match ISO 400.** It was also measured on the older TRI-X Pan 5063 film. Use the shape only and flag it in the data provenance.
4. **No published grain data for either film.** The size-class grids (ECD, thickness, lognormal width, silver coverage) are `fitted` parameters with `assumed` priors. Priors:
   - Tri-X: Kofron controls / textbook ranges.
   - T-Max 100: Kofron tabular examples, scaled down for RMS 8.
   - Silver coverage: no B&W primary value exists. Use a prior and let D_max and covering power constrain it.
5. **n_dev links to development.** The minimum developable size depends on developer potential. The latent-image output should be the full cluster-size distribution per size class, as the design already says, not a binary developable flag. Then D-76 vs XTOL vs HC-110 speed differences can emerge from developer-dependent n_dev and induction time rather than from per-pair fudge factors.
6. **Characteristic-curve developer.** F-4017 (Oct 2021) has 400TX characteristic curves for both **T-MAX Developer** (p. 8, small tank, 6/7/9/11 min) and **D-76** (p. 9, large tank, 1-min agitation, 20 °C, 6/8/10/12 min), for 35 mm and 120; re-checked on the rendered pages. Granularity is HC-110 dil. B; MTF is D-76. HC-110 curves exist only for 320TXP sheets, and there are no XTOL curves (see `datasheets.md`). F-4016 (2016) has D-76, T-MAX RS and T-MAX Developer curves for T-Max 100.
7. **Unit tests suggested by the literature.** With τ₀ → ∞ and τ_i → 0 the chain must reduce to the Poisson r-hit model. Quantum sensitivity at the optimum exposure time should land in the 6–27 photons/grain range reported by Hailstone and colleagues.

## 6. What would confirm the recalled items

| Item | Where to confirm |
|---|---|
| Ionic-step time, room-temperature electron mobility (direct measurement), hole drift mobility, recombination radius, trap densities (all **UNVERIFIED**) | Hamilton, Adv. Phys. 37, 359 (1988), §§3–6; Hailstone, J. Phys. Chem. 99(13), 4414–4428 (1995), DOI 10.1021/j100013a009, parameter table; Hanson & Brown, J. Appl. Phys. 31, 210 (1960) for the electron reference they cite |
| Lognormal ranges for conventional emulsions | James (ed.), *TPP* 4th ed. (1977), ch. 3; Tani (1995) ch. 1 |
| AgBr α(λ) | Moser & Urbach, Phys. Rev. 102, 1519 (1956) |
| J-aggregate ε | any cyanine J-aggregate spectroscopy paper (e.g. in Tani 1995, ch. 6) |
| ρ_AgBr = 6.473 g cm⁻³, gelatin ≈ 1.35 g cm⁻³ | CRC Handbook |
| ISO speed definition (S = 0.8/H_m at D_min + 0.10) | ISO 6:1993 |
