"""The shared kinetics constants file ``data/chemistry/kinetics.yaml`` (spec D1 2.3.6, 2.3.7,
2.4.11, 5.3, 5.8, 7.1)."""

from pathlib import Path

import pytest
from typer.testing import CliRunner

from darkroom_atlas.cli import app
from darkroom_atlas.data.library import Library
from darkroom_atlas.data.models import KineticsFile, iter_quantities
from darkroom_atlas.data.provenance import Source

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
FIXTURES = ROOT / "tests" / "fixtures" / "library"

LOG = "lognormal"
UNI = "uniform"
NORM = "normal"

# field path -> (value, unit, source, bounds, prior). Values the spec gives are pinned; the
# placeholder values of fitted parameters follow the spec D1 5.8 example. Bounds and priors
# are the "prior/bounds" columns of spec D1 2.3.6 and 2.4.11 ("(log)" -> lognormal, linear
# box of a fitted parameter -> uniform, N(mu, sigma^2) -> normal).
_AGENT_COMMON = {
    "k_et": (0.05, "1/s", "assumed", (1e-4, 10.0), LOG),
    "K_ads": (100.0, "L/mol", "assumed", (1.0, 1e4), LOG),
    "lambda0": (0.01, "1/s", "assumed", (1e-4, 10.0), LOG),
}
_K_DES = (0.1, "1/s", "assumed", (1e-3, 10.0), LOG)
_K_REG = (100.0, "L/(mol*s)", "assumed", (1e-2, 1e4), LOG)
_EA_FIXED = (65.0, "kJ/mol", "assumed", None, None)


def _agent(prefix, extra, fixed=False):
    rows = {**_AGENT_COMMON, **extra}
    out = {}
    for name, row in rows.items():
        if fixed and row[4] in (LOG, UNI):
            row = (*row[:4], None)  # p-aminophenol: metol values, not fitted (spec D1 7.1)
        out[f"agents.{prefix}.{name}"] = row
    return out


EXPECTED = {
    # Spec D1 2.3.3, 2.3.6.
    "activity.debye_huckel_A_20C": (0.5050, "1", "assumed", None, None),
    "activity.debye_huckel_A_25C": (0.5091, "1", "assumed", None, None),
    "activity.sit_b": (1.5, "1", "assumed", None, None),
    **{
        f"activity.delta_eps.{c}": (0.0, "kg/mol", "assumed", (-0.45, 0.45), NORM)
        for c in ("sulfite", "borate", "organic", "amine")
    },
    # Spec D1 2.4.11 agent table, 2.4.5, 7.1.
    **_agent(
        "metol",
        {
            "activation_energy": (61.8, "kJ/mol", "literature", None, NORM),
            "k_des": _K_DES,
        },
    ),
    **_agent(
        "hydroquinone",
        {
            "activation_energy": (67.5, "kJ/mol", "literature", None, NORM),
            "k_reg": _K_REG,
        },
    ),
    **_agent(
        "pyrazolidone",
        {
            "activation_energy": _EA_FIXED,
            "k_des": _K_DES,
            "pka": (8.5, "1", "assumed", (7.5, 9.5), UNI),
        },
    ),
    **_agent("ascorbate", {"activation_energy": _EA_FIXED, "k_reg": _K_REG}),
    **_agent("p-aminophenol", {"activation_energy": _EA_FIXED, "k_des": _K_DES}, fixed=True),
    # Spec D1 2.4.11 shared rows.
    "halides.K_Br": (300.0, "L/mol", "assumed", (10.0, 1e4), LOG),
    "halides.adsorption_enthalpy": (-45.0, "kJ/mol", "literature", (-60.0, -30.0), None),
    "halides.iodide_ratio": (30.0, "1", "assumed", (10.0, 1000.0), None),
    "barriers.u_gel": (2.0, "1", "assumed", (0.0, 10.0), UNI),
    "barriers.u_halide": (5.0, "1", "assumed", (0.0, 20.0), UNI),
    "barriers.ionic_strength_ref": (1.0, "mol/L", "assumed", None, None),  # spec D1 2.4.2
    "latent.u_li": (1.0, "1", "assumed", (0.0, 5.0), UNI),
    "latent.n_sat": (5.0, "1", "literature", None, None),
    "solvent.latent_attrition": (1.0, "1", "assumed", (0.0, 5.0), UNI),
    "solvent.k_diss": (1e-4, "um/s", "assumed", (1e-6, 1e-2), LOG),
    "solvent.c_sat_ref": (1e-3, "mol/L", "assumed", (1e-5, 1e-2), LOG),
    "solvent.sulfite_exponent": (1.0, "1", "assumed", (0.6, 2.0), UNI),
    "solvent.bromide_exponent": (0.4, "1", "assumed", (0.3, 1.0), None),
    "solvent.ligands": (2.0, "1", "assumed", (1.0, 3.0), None),
    "solvent.k_pd": (1e-4, "um/s", "assumed", (1e-6, 1e-2), LOG),
    "solvent.sulfonation_switch": (1e-4, "mol/L", "assumed", None, None),  # c_s, 2.4.6
    "transport.delta_max": (300.0, "um", "assumed", (200.0, 500.0), UNI),
    "transport.delta_agit": (50.0, "um", "assumed", (20.0, 100.0), UNI),
    "transport.swelling_half_time": (20.0, "s", "assumed", (10.0, 40.0), None),
    "transport.viscosity_ionic_coefficient": (0.125, "L/mol", "assumed", (0.0, 0.2), None),
    "transport.water_viscosity.temperature": (
        [18.0, 20.0, 24.0, 25.0],
        "C",
        "assumed",
        None,
        None,
    ),
    "transport.water_viscosity.viscosity": (
        [1.053, 1.002, 0.911, 0.890],
        "cP",
        "assumed",
        None,
        None,
    ),
    "transport.proton_diffusion": (8130.0, "um2/s", "assumed", None, None),
    "transport.hydroxide_diffusion": (4600.0, "um2/s", "assumed", None, None),
    "gelatin.groups_per_gram": (4.8e-4, "mol/g", "literature", None, None),
    "gelatin.pka": (10.5, "1", "assumed", (9.5, 11.5), None),
    "gelatin.initial_ph": (5.5, "1", "assumed", (5.0, 6.5), None),
    "gelatin.dry_density": (1.35, "g/cm3", "assumed", None, None),
}

# Parameters calibration fits (spec D1 7.1 blocks A activity, C agents, C shared).
FITTED = {
    *(f"activity.delta_eps.{c}" for c in ("sulfite", "borate", "organic", "amine")),
    *(f"agents.metol.{n}" for n in ("k_et", "K_ads", "lambda0", "k_des", "activation_energy")),
    *(
        f"agents.hydroquinone.{n}"
        for n in ("k_et", "K_ads", "lambda0", "k_reg", "activation_energy")
    ),
    *(f"agents.pyrazolidone.{n}" for n in ("k_et", "K_ads", "lambda0", "k_des", "pka")),
    *(f"agents.ascorbate.{n}" for n in ("k_et", "K_ads", "lambda0", "k_reg")),
    "halides.K_Br",
    "barriers.u_gel",
    "barriers.u_halide",
    "latent.u_li",
    "solvent.latent_attrition",
    "solvent.k_diss",
    "solvent.c_sat_ref",
    "solvent.sulfite_exponent",
    "solvent.k_pd",
    "transport.delta_max",
    "transport.delta_agit",
}


@pytest.fixture(scope="module")
def kinetics() -> KineticsFile:
    return Library().kinetics()


@pytest.fixture(scope="module")
def quantities(kinetics) -> dict:
    return dict(iter_quantities(kinetics))


def test_file_header(kinetics):
    assert kinetics.schema_version == 1
    assert kinetics.model_version == "D1c"
    assert kinetics.calibration is None  # prior values: no calibration has run
    assert sorted(kinetics.agents) == sorted(
        ["metol", "hydroquinone", "pyrazolidone", "ascorbate", "p-aminophenol"]
    )


def test_every_constant_is_listed(quantities):
    assert sorted(quantities) == sorted(EXPECTED)


@pytest.mark.parametrize("path", sorted(EXPECTED))
def test_constant_matches_spec_tables(quantities, path):
    value, unit, source, bounds, prior = EXPECTED[path]
    q = quantities[path]
    assert q.value == value
    assert q.unit == unit
    assert q.source == Source(source)
    assert q.bounds == bounds
    assert q.prior == prior
    assert q.fit is None


def test_agent_optional_fields(kinetics):
    # ETAs carry k_des, regenerators k_reg; only the pyrazolidone class has a pKa.
    for name, entry in kinetics.agents.items():
        eta = name in ("metol", "pyrazolidone", "p-aminophenol")
        assert (entry.k_des is not None) is eta, name
        assert (entry.k_reg is not None) is (not eta), name
        assert (entry.pka is not None) is (name == "pyrazolidone"), name


def test_literature_only_where_the_spec_says(quantities):
    literature = {p for p, q in quantities.items() if q.source == Source.LITERATURE}
    assert literature == {
        "agents.metol.activation_energy",
        "agents.hydroquinone.activation_energy",
        "halides.adsorption_enthalpy",
        "latent.n_sat",
        "gelatin.groups_per_gram",
    }
    for path in literature:
        assert quantities[path].ref, path
    for path, q in quantities.items():
        assert q.source in (Source.LITERATURE, Source.ASSUMED), path
        if q.source == Source.ASSUMED:
            assert q.note, path  # every assumed value says where it comes from


def test_priors_exactly_on_the_fitted_parameters(quantities):
    with_prior = {p for p, q in quantities.items() if q.prior is not None}
    assert with_prior == FITTED
    for path in FITTED - {
        "agents.metol.activation_energy",
        "agents.hydroquinone.activation_energy",
    }:
        assert "PLACEHOLDER" in quantities[path].note, path


def test_normal_priors_have_their_sigma(quantities):
    assert quantities["agents.metol.activation_energy"].uncertainty == 10.0
    assert quantities["agents.hydroquinone.activation_energy"].uncertainty == 10.0
    for c in ("sulfite", "borate", "organic", "amine"):
        assert quantities[f"activity.delta_eps.{c}"].uncertainty == 0.15


def test_values_inside_their_bounds(quantities):
    for path, q in quantities.items():
        if q.bounds is not None:
            lo, hi = q.bounds
            assert lo <= q.value <= hi, path
            if q.prior == LOG:
                assert lo > 0.0, path


def test_sugai_activation_energies_from_ev():
    # 0.64 and 0.70 eV (spec D1 2.4.5); 1 eV = 96.4853 kJ/mol (spec D1 5.1).
    assert 0.64 * 96.4853 == pytest.approx(61.8, abs=0.05)
    assert 0.70 * 96.4853 == pytest.approx(67.5, abs=0.05)


def test_water_ion_diffusion_from_tracer_values(kinetics):
    # 25 C tracer values scaled to 20 C by T / eta (spec D1 2.3.7).
    scale = (293.15 / 298.15) * (0.890 / 1.002)
    t = kinetics.transport
    assert t.proton_diffusion.to("um2/s") == pytest.approx(9.31e-5 * 1e8 * scale, rel=1e-3)
    assert t.hydroxide_diffusion.to("um2/s") == pytest.approx(5.27e-5 * 1e8 * scale, rel=1e-3)


def test_halide_enthalpy_matches_herz():
    # -10.7 / -11.5 kcal/mol (Herz 1980, spec D1 2.4.11) bracket -45 kJ/mol.
    assert -11.5 * 4.184 < -45.0 < -10.7 * 4.184


def test_c_sat_anchor_extrapolation():
    # Carroll & Hubbard 7.5e-5 M at 0.005 M sulfite, extrapolated to 0.1 M (spec D1 2.4.11).
    assert 7.5e-5 * (0.1 / 0.005) ** 0.6 == pytest.approx(4.5e-4, rel=0.03)
    assert 7.5e-5 * (0.1 / 0.005) ** 1.0 == pytest.approx(1.5e-3, rel=1e-12)


def test_check_data_accepts_the_data_library():
    result = CliRunner().invoke(app, ["check-data", "--data", str(DATA)])
    assert result.exit_code == 0, result.output
    assert "OK:" in result.output


def test_fixture_copy_has_the_prior_values(kinetics):
    # The fixture library carries the same prior file (spec D1 10 D1d).
    fixture = Library(FIXTURES).kinetics()
    assert fixture == kinetics
