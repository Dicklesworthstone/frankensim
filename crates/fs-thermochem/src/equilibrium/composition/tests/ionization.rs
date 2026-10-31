//! Independent ideal atomic-plasma oracles, not fitted plasma validation.
//!
//! The source cards use the monatomic Sackur--Tetrode entropy and constant
//! Cp=5R/2. Their Gibbs minimizer must reproduce Saha without invoking the
//! production equilibrium arithmetic in the oracle. Electron spin is counted
//! once, in its degeneracy. Actual translational masses retain the small
//! (m_proton/m_atom)^(3/2) correction to the usual heavy-mass approximation.
//!
//! References: NASA RP-1311, section 3.7 (charge-constrained ideal equilibrium),
//! <https://ntrs.nasa.gov/citations/19950013764>; the Saha/pressure/energy/Cv
//! derivation <https://www.astro.princeton.edu/~gk/A403/ioniz.pdf>, i.3--i.9;
//! NIST ASD hydrogen ionization energy 13.598434599702 eV,
//! <https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=H+I>.
//! Excited levels, nonideal plasma corrections and nonequilibrium rates are
//! intentionally absent from this independently soluble test model.

use super::*;

const K_B: f64 = 1.380_649e-23;
const PLANCK: f64 = 6.626_070_15e-34;
const AVOGADRO: f64 = 6.022_140_76e23;
const ELECTRON_MASS: f64 = 9.109_383_713_9e-31;
const PROTON_MASS: f64 = 1.672_621_925_95e-27;
const ATOM_MASS: f64 = PROTON_MASS + ELECTRON_MASS;
const ELECTRON_MOLAR_MASS: f64 = ELECTRON_MASS * AVOGADRO;
const ATOM_MOLAR_MASS: f64 = ATOM_MASS * AVOGADRO;
const IONIZATION_J: f64 = 13.598_434_599_702 * 1.602_176_634e-19;
const IONIZATION_J_MOL: f64 = IONIZATION_J * AVOGADRO;
const P0: f64 = 100_000.0;

fn relative_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite()
            && expected.is_finite()
            && (actual - expected).abs() <= tolerance * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}, relative tolerance={tolerance}"
    );
}

fn monatomic(
    name: &str,
    mass_kg: f64,
    degeneracy: f64,
    formation_j: f64,
) -> Nasa9StandardStateModelV1 {
    let entropy_constant = 2.5
        + (degeneracy * K_B / P0
            * (2.0 * core::f64::consts::PI * mass_kg * K_B / PLANCK.powi(2)).powf(1.5))
        .ln();
    let coefficients = [
        0.0,
        0.0,
        2.5,
        0.0,
        0.0,
        0.0,
        0.0,
        formation_j / K_B,
        entropy_constant,
    ];
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(name.to_string(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: P0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 50.0, 100_000.0),
        sources: Vec::new(),
        provenance: Provenance {
            source: "ideal monatomic statistical mechanics; NIST H ionization energy".into(),
            license: "analytic test fixture".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        id(name),
        MolarMass::new(mass_kg * AVOGADRO),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(P0),
            ElementalReferenceIdV1::new("atomic-ground-state-test").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}

fn charged_owner(
    models: Vec<Nasa9StandardStateModelV1>,
    nuclei: Vec<i128>,
    charges: Vec<i128>,
) -> Result<IdealGasEquilibriumV1, CompositionEquilibriumErrorV1> {
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![nuclei],
        )
        .unwrap(),
        ChargeVector::new(ids, charges).unwrap(),
        models,
    )
}

fn hydrogen(reversed: bool) -> IdealGasEquilibriumV1 {
    let mut models = vec![
        monatomic("H", ATOM_MASS, 2.0, 0.0),
        monatomic("H+", PROTON_MASS, 1.0, IONIZATION_J),
        monatomic("e-", ELECTRON_MASS, 2.0, 0.0),
    ];
    let mut nuclei = vec![1, 1, 0];
    let mut charges = vec![0, 1, -1];
    if reversed {
        models.reverse();
        nuclei.reverse();
        charges.reverse();
    }
    charged_owner(models, nuclei, charges).unwrap()
}

fn saha(t: f64) -> f64 {
    let translational_mass = ELECTRON_MASS * (PROTON_MASS / ATOM_MASS);
    (2.0 * core::f64::consts::PI * translational_mass * K_B * t / PLANCK.powi(2)).powf(1.5)
        * (-IONIZATION_J / (K_B * t)).exp()
}

fn fraction_tp(t: f64, p: f64) -> f64 {
    let a = saha(t) * K_B * t / p;
    (a / (1.0 + a)).sqrt()
}

fn fraction_tv(t: f64, volume: f64) -> f64 {
    let r = saha(t) * volume / AVOGADRO; // One mole of nuclei.
    let root = r.sqrt();
    2.0 * root / ((r + 4.0).sqrt() + root)
}

fn accepted_charge(solution: &IdealGasEquilibriumSolutionV1) {
    let positive = species(solution, "H+").amount.value();
    let negative = species(solution, "e-").amount.value();
    relative_close(negative, positive, 2e-9);
    assert!(solution.residuals().relative_charge_residual.abs() <= 1e-10);
}

#[test]
fn g1_ionization_tp_matches_saha_from_a_neutral_unseeded_feed() {
    let chemistry = hydrogen(false);
    let feed = [(id("H"), Amount::new(1.0))];
    for (t, p) in [(4000.0, P0), (8000.0, 1000.0), (15000.0, P0), (25000.0, P0)] {
        let state = chemistry
            .equilibrate_tp(
                Temperature::new(t),
                Pressure::new(p),
                &feed,
                EquilibriumSolveV1::default(),
                || false,
            )
            .unwrap();
        let x = fraction_tp(t, p);
        relative_close(species(&state, "H+").amount.value(), x, 2e-8);
        relative_close(species(&state, "e-").amount.value(), x, 2e-8);
        relative_close(species(&state, "H").amount.value(), 1.0 - x, 2e-8);
        relative_close(state.mass().value(), ATOM_MOLAR_MASS, 2e-10);
        relative_close(
            state.enthalpy().value(),
            2.5 * R * t * (1.0 + x) + x * IONIZATION_J_MOL,
            2e-8,
        );
        accepted_charge(&state);
        assert!(state.species().iter().all(|s| s.source.is_some()));
    }
}

#[test]
fn g1_ionization_tv_includes_electron_pressure_and_ionization_energy() {
    let chemistry = hydrogen(false);
    let feed = [(id("H"), Amount::new(1.0))];
    let t = 15000.0;
    for volume in [0.01, 1.0, 100.0] {
        let state = chemistry
            .equilibrate_tv(
                Temperature::new(t),
                Volume::new(volume),
                &feed,
                EquilibriumVolumeSolveV1::default(),
                || false,
            )
            .unwrap();
        let x = fraction_tv(t, volume);
        relative_close(species(state.equilibrium(), "H+").amount.value(), x, 2e-8);
        relative_close(
            state.equilibrium().pressure().value(),
            (1.0 + x) * R * t / volume,
            2e-8,
        );
        relative_close(
            state.internal_energy().value(),
            1.5 * R * t * (1.0 + x) + x * IONIZATION_J_MOL,
            2e-8,
        );
        relative_close(
            state.thermodynamic().response().density_kg_m3,
            ATOM_MOLAR_MASS / volume,
            2e-9,
        );
        accepted_charge(state.equilibrium());
    }
}

#[test]
fn g1_ionization_equilibrium_caloric_and_acoustic_derivatives_match_saha() {
    let chemistry = hydrogen(false);
    let feed = [(id("H"), Amount::new(1.0))];
    for t in [8000.0, 15000.0, 25000.0] {
        let state = chemistry
            .equilibrate_tp_with_response(
                Temperature::new(t),
                Pressure::new(P0),
                &feed,
                EquilibriumSolveV1::default(),
                || false,
            )
            .unwrap();
        let x = fraction_tp(t, P0);
        let d = 1.5 + IONIZATION_J / (K_B * t);
        let q = x * (1.0 - x) / (2.0 - x);
        let cv_r = 1.5 * (1.0 + x) + q * d * d;
        let cp_r = 2.5 * (1.0 + x) + 0.5 * x * (1.0 - x * x) * (d + 1.0).powi(2);
        let sound_squared =
            R * t / ATOM_MOLAR_MASS * (1.0 + x - q + (1.0 + x + q * d).powi(2) / cv_r);
        let response = state.response();
        relative_close(response.heat_capacity_p_j_k, R * cp_r, 2e-8);
        relative_close(response.heat_capacity_v_j_k, R * cv_r, 2e-8);
        relative_close(response.sound_speed_m_s.powi(2), sound_squared, 2e-8);
        relative_close(
            response.thermal_expansion_per_k,
            (1.0 + 0.5 * x * (1.0 - x) * (d + 1.0)) / t,
            2e-8,
        );
        relative_close(
            response.isothermal_compressibility_per_pa,
            (1.0 + 0.5 * x * (1.0 - x)) / P0,
            2e-8,
        );
    }
}

fn bisect_energy(target: f64, initial_t: f64, energy: impl Fn(f64) -> f64) -> f64 {
    let (mut lo, mut hi) = (500.0, initial_t);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if energy(mid) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn g1_ionization_hp_and_uv_spend_thermal_energy_on_ionization() {
    let chemistry = hydrogen(false);
    let feed = [(id("H"), Amount::new(1.0))];
    let initial_t = 30000.0;
    let enthalpy = 2.5 * R * initial_t;
    let hp_expected = bisect_energy(enthalpy, initial_t, |t| {
        let x = fraction_tp(t, P0);
        2.5 * R * t * (1.0 + x) + x * IONIZATION_J_MOL
    });
    let hp = chemistry
        .equilibrate_hp(
            Pressure::new(P0),
            &feed,
            Temperature::new(initial_t),
            EquilibriumEnthalpySolveV1 {
                minimum_temperature_k: 500.0,
                maximum_temperature_k: initial_t,
                enthalpy_tolerance_j: 1e-5,
                ..hp_controls()
            },
            || false,
        )
        .unwrap();
    relative_close(hp.temperature().value(), hp_expected, 2e-8);
    assert!(hp.temperature().value() < 0.75 * initial_t);
    assert!(hp.enthalpy_residual().unwrap().value().abs() <= 1e-5);
    accepted_charge(&hp);

    let volume = 0.5;
    let energy = 1.5 * R * initial_t;
    let uv_expected = bisect_energy(energy, initial_t, |t| {
        let x = fraction_tv(t, volume);
        1.5 * R * t * (1.0 + x) + x * IONIZATION_J_MOL
    });
    let uv = chemistry
        .equilibrate_uv(
            Volume::new(volume),
            &feed,
            Temperature::new(initial_t),
            EquilibriumInternalEnergySolveV1 {
                minimum_temperature_k: 500.0,
                maximum_temperature_k: initial_t,
                energy_tolerance_j: 1e-5,
                ..uv_controls()
            },
            || false,
        )
        .unwrap();
    relative_close(uv.equilibrium().temperature().value(), uv_expected, 2e-8);
    assert!(uv.energy_residual().unwrap().value().abs() <= 1e-5);
    accepted_charge(uv.equilibrium());
}

#[test]
fn g1_ionization_signed_negative_and_multiple_charge_constraints() {
    // At T_ref and p0, author g_i/RT=-ln(y_i). Every chemical potential is
    // exactly zero at the chosen feasible composition, hence this is the
    // independently known strict Gibbs minimum. Both negative ions and free
    // electrons must enter charge balance, with +2 counted twice.
    let expected: [f64; 4] = [0.5, 0.2, 0.3, 0.4];
    let chemistry = charged_owner(
        vec![
            model("X", 0.01, 2.5, 0.0, -(expected[0] / 1.4).ln()),
            model(
                "Xm",
                0.01 + ELECTRON_MOLAR_MASS,
                2.5,
                0.0,
                -(expected[1] / 1.4).ln(),
            ),
            model(
                "Xpp",
                0.01 - 2.0 * ELECTRON_MOLAR_MASS,
                2.5,
                0.0,
                -(expected[2] / 1.4).ln(),
            ),
            model(
                "e-",
                ELECTRON_MOLAR_MASS,
                2.5,
                0.0,
                -(expected[3] / 1.4).ln(),
            ),
        ],
        vec![1, 1, 1, 0],
        vec![0, -1, 2, -1],
    )
    .unwrap();
    let result = chemistry
        .equilibrate_tp(
            Temperature::new(REFERENCE_T),
            Pressure::new(P0),
            &[(id("X"), Amount::new(1.0))],
            EquilibriumSolveV1::default(),
            || false,
        )
        .unwrap();
    for (name, expected) in ["X", "Xm", "Xpp", "e-"].into_iter().zip(expected) {
        relative_close(species(&result, name).amount.value(), expected, 2e-9);
    }
    relative_close(result.mass().value(), 0.01, 2e-10);
    assert!(result.residuals().relative_charge_residual.abs() <= 1e-10);
}

#[test]
fn g1_ionization_trace_balance_is_scaled_by_actual_charge_carriers() {
    let chemistry = hydrogen(false);
    let t = 2000.0;
    let state = chemistry
        .equilibrate_tp(
            Temperature::new(t),
            Pressure::new(P0),
            &[(id("H"), Amount::new(1.0))],
            EquilibriumSolveV1::default(),
            || false,
        )
        .unwrap();
    let expected = fraction_tp(t, P0);
    assert!(expected > 0.0 && expected < 1e-14);
    relative_close(species(&state, "H+").amount.value(), expected, 2e-8);
    relative_close(species(&state, "e-").amount.value(), expected, 2e-8);
    accepted_charge(&state);
}

#[test]
fn g1_ionization_retains_charge_balance_beyond_direct_amount_underflow() {
    // A mathematical ideal-atomic-gas fixture, not a claim that cold physical
    // hydrogen stays atomic. Neither charged amount is representable here.
    let chemistry = hydrogen(false);
    let t = 50.0;
    let state = chemistry
        .equilibrate_tp(
            Temperature::new(t),
            Pressure::new(P0),
            &[(id("H"), Amount::new(1.0))],
            EquilibriumSolveV1::default(),
            || false,
        )
        .unwrap();
    let mass = ELECTRON_MASS * (PROTON_MASS / ATOM_MASS);
    let log_saha = 1.5 * (2.0 * core::f64::consts::PI * mass * K_B * t / PLANCK.powi(2)).ln()
        - IONIZATION_J / (K_B * t);
    let log_x = 0.5 * (log_saha + (K_B * t / P0).ln());
    for name in ["H+", "e-"] {
        let charged = species(&state, name);
        assert_eq!(charged.amount.value(), 0.0);
        assert!((charged.log_mole_fraction.unwrap() - log_x).abs() < 1e-8);
    }
    assert!(state.residuals().relative_charge_residual < 1e-10);
    assert!(state.is_bound_to(&chemistry));
    assert!(state.is_bound_to(&hydrogen(true)));
    let different_charges =
        charged_owner(chemistry.models().to_vec(), vec![1, 1, 0], vec![0, 2, -1]).unwrap();
    assert!(!state.is_bound_to(&different_charges));
}

#[test]
fn g3_ionization_preserves_axis_permutation_and_feed_scale() {
    let chemistry = hydrogen(false);
    assert_eq!(chemistry, hydrogen(true));
    let expected = fraction_tp(15000.0, P0);
    for scale in [1e-12, 1.0, 1e12] {
        let state = chemistry
            .equilibrate_tp(
                Temperature::new(15000.0),
                Pressure::new(P0),
                &[(id("H"), Amount::new(scale))],
                EquilibriumSolveV1::default(),
                || false,
            )
            .unwrap();
        relative_close(species(&state, "H+").amount.value(), expected * scale, 2e-8);
        relative_close(state.mass().value(), ATOM_MOLAR_MASS * scale, 2e-10);
        accepted_charge(&state);
    }
}

#[test]
fn g0_ionization_refuses_nonphysical_nucleus_free_species_and_net_charge() {
    for charge in [0, 1, -2] {
        assert!(matches!(
            charged_owner(
                vec![
                    model("H", 0.001, 2.5, 0.0, 0.0),
                    model("empty", ELECTRON_MOLAR_MASS, 2.5, 0.0, 0.0)
                ],
                vec![1, 0],
                vec![0, charge],
            ),
            Err(CompositionEquilibriumErrorV1::Invalid(_))
        ));
    }
    assert!(matches!(
        charged_owner(
            vec![
                model("H", 0.001, 2.5, 0.0, 0.0),
                model("e1", ELECTRON_MOLAR_MASS, 2.5, 0.0, 0.0),
                model("e2", ELECTRON_MOLAR_MASS, 2.5, 0.0, 0.0)
            ],
            vec![1, 0, 0],
            vec![0, -1, -1],
        ),
        Err(CompositionEquilibriumErrorV1::Invalid(_))
    ));
    assert!(matches!(
        hydrogen(false).equilibrate_tp(
            Temperature::new(15000.0),
            Pressure::new(P0),
            &[(id("H+"), Amount::new(1.0)), (id("e-"), Amount::new(0.1))],
            EquilibriumSolveV1::default(),
            || false,
        ),
        Err(CompositionEquilibriumErrorV1::Invalid(_))
    ));
}
