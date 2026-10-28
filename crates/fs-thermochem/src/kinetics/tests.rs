use super::*;
use crate::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, InitialStatePolicy,
    NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, StandardStateConventionV1,
    StandardStatePhaseV1, required_nasa9_parameter_dims,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, LawId, LawParameter, Provenance};
use fs_qty::{Amount, MolarMass, Volume};
use std::collections::BTreeMap;

#[path = "pressure/tests.rs"]
mod pressure_rates;

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn source(
    name: &str,
    mass: f64,
    cp_r: f64,
    formation: f64,
    entropy_offset: f64,
) -> Nasa9StandardStateModelV1 {
    let coefficients = [
        0.0,
        0.0,
        cp_r,
        0.0,
        0.0,
        0.0,
        0.0,
        formation / R,
        entropy_offset,
    ];
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(name.into(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 6000.0),
        sources: vec![],
        provenance: Provenance {
            source: "manufactured kinetic conformance".into(),
            license: "fixture".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        id(name),
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("kinetics-fixture").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
fn thermo(models: Vec<Nasa9StandardStateModelV1>, atoms: Vec<i128>) -> IdealGasEquilibriumV1 {
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    IdealGasEquilibriumV1::new(
        ElementalMatrix::new(vec![ElementId::new("H").unwrap()], ids.clone(), vec![atoms]).unwrap(),
        ChargeVector::new(ids.clone(), vec![0; ids.len()]).unwrap(),
        models,
    )
    .unwrap()
}
fn rate(scale: f64, b: f64, ea: f64) -> ArrheniusRateV1 {
    ArrheniusRateV1::new(
        Temperature::new(500.0),
        ReactionProgressRateQuantityV1::new(scale),
        b,
        MolarEnergyQuantityV1::new(ea),
    )
    .unwrap()
}
fn reaction(
    name: &str,
    lhs: &[(&str, u32)],
    rhs: &[(&str, u32)],
    law: ArrheniusRateV1,
) -> ElementaryReactionV1 {
    ElementaryReactionV1::new(
        ReactionId::new(name).unwrap(),
        lhs.iter().map(|(s, n)| (id(s), *n)).collect(),
        rhs.iter().map(|(s, n)| (id(s), *n)).collect(),
        law,
    )
    .unwrap()
}
fn isomers(q: f64, law: ArrheniusRateV1) -> IdealGasKineticsV1 {
    IdealGasKineticsV1::new(
        thermo(
            vec![
                source("A", 0.01, 4.0, 0.0, 0.0),
                source("B", 0.01, 4.0, -q, 0.0),
            ],
            vec![1, 1],
        ),
        "isomer-fixture",
        1,
        vec![reaction("A-to-B", &[("A", 1)], &[("B", 1)], law)],
        1e-12,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tol={tolerance:e}"
    );
}

#[test]
fn g1_arrhenius_reverse_rate_temperature_response_and_energy_share_sources() {
    let m = isomers(12000.0, rate(7.0, 0.7, 25000.0));
    let mut w = KineticWorkspaceV1::new(&m).unwrap();
    for t in [300.0, 500.0, 1100.0] {
        let c = [2.0, 0.5];
        let e = m
            .evaluate(Temperature::new(t), &c, &mut w, || false)
            .unwrap();
        let kf = 7.0 * (t / 500.0).powf(0.7) * (25000.0 / R * (1.0 / 500.0 - 1.0 / t)).exp();
        let qr = kf * c[1] * (-12000.0 / (R * t)).exp();
        let net = kf * c[0] - qr;
        close(w.reaction_rates()[0], net, net.abs() * 2e-13);
        close(w.species_rates()[0], -net, net.abs() * 2e-13);
        close(w.species_rates()[1], net, net.abs() * 2e-13);
        close(w.rate_jacobian()[0], kf, kf * 2e-13);
        close(
            w.rate_jacobian()[1],
            -kf * (-12000.0 / (R * t)).exp(),
            kf * 2e-13,
        );
        close(
            w.temperature_derivatives()[0],
            net * (0.7 / t + 25000.0 / (R * t * t)) - qr * 12000.0 / (R * t * t),
            kf * 2e-13,
        );
        close(
            e.adiabatic_temperature_rate_k_s,
            12000.0 * net / (2.5 * 3.0 * R),
            kf * 1e-10,
        );
        close(
            e.entropy_production_w_m3_k,
            R * net * (12000.0 / (R * t) + (c[0] / c[1]).ln()),
            kf * 1e-10,
        );
        assert!(e.entropy_production_w_m3_k > 0.0);
    }
}

#[test]
fn g1_dissociation_uses_concentration_standard_state_and_matches_gibbs_equilibrium() {
    let tr = 1000.0_f64;
    let th = thermo(
        vec![
            source("A", 0.001, 2.5, 0.0, 2.5 * (1.0 - tr.ln())),
            source("A2", 0.002, 3.5, 0.0, 3.5 * (1.0 - tr.ln())),
        ],
        vec![1, 2],
    );
    let m = IdealGasKineticsV1::new(
        th.clone(),
        "dissociation-fixture",
        1,
        vec![reaction(
            "dissociation",
            &[("A2", 1)],
            &[("A", 2)],
            rate(4.0, 0.0, 0.0),
        )],
        1e-12,
    )
    .unwrap();
    let mut w = KineticWorkspaceV1::new(&m).unwrap();
    let kc = 100000.0 / (R * tr);
    let c = [kc.sqrt(), 1.0];
    m.evaluate(Temperature::new(tr), &c, &mut w, || false)
        .unwrap();
    close(w.log_equilibrium_constants()[0], kc.ln(), 3e-14);
    close(w.reaction_rates()[0], 0.0, 2e-13);
    for v in [0.001, 0.02, 1.0] {
        let equilibrium = th
            .equilibrate_tv(
                Temperature::new(tr),
                Volume::new(v),
                &[(id("A2"), Amount::new(1.0))],
                crate::equilibrium::composition::volume::EquilibriumVolumeSolveV1::default(),
                || false,
            )
            .unwrap();
        let c: Vec<_> = equilibrium
            .equilibrium()
            .species()
            .iter()
            .map(|s| s.amount.value() / v)
            .collect();
        m.evaluate(Temperature::new(tr), &c, &mut w, || false)
            .unwrap();
        assert!(w.reaction_rates()[0].abs() < 1e-7 * c[1].max(1.0));
    }
    let c = [2.0, 3.0];
    m.evaluate(Temperature::new(tr), &c, &mut w, || false)
        .unwrap();
    let derivative = w.temperature_derivatives()[0];
    let h = 0.001;
    m.evaluate(Temperature::new(tr + h), &c, &mut w, || false)
        .unwrap();
    let plus = w.reaction_rates()[0];
    m.evaluate(Temperature::new(tr - h), &c, &mut w, || false)
        .unwrap();
    let minus = w.reaction_rates()[0];
    close(derivative, (plus - minus) / (2.0 * h), 1e-10);
}

#[test]
fn g0_zero_products_and_catalysts_have_correct_inward_rates_and_derivatives() {
    let m = isomers(0.0, rate(3.0, 0.0, 0.0));
    let mut w = KineticWorkspaceV1::new(&m).unwrap();
    let e = m
        .evaluate(Temperature::new(500.0), &[2.0, 0.0], &mut w, || false)
        .unwrap();
    close(w.reaction_rates()[0], 6.0, 1e-14);
    close(w.rate_jacobian()[0], 3.0, 1e-14);
    close(w.rate_jacobian()[1], -3.0, 1e-14);
    assert_eq!(e.entropy_production_w_m3_k, f64::INFINITY);
    let th = thermo(
        vec![
            source("A", 0.01, 4.0, 0.0, 0.0),
            source("B", 0.01, 4.0, 0.0, 0.0),
            source("X", 0.01, 3.0, 0.0, 0.0),
        ],
        vec![1, 1, 1],
    );
    let m = IdealGasKineticsV1::new(
        th,
        "catalyst-fixture",
        1,
        vec![reaction(
            "catalyzed",
            &[("A", 1), ("X", 1)],
            &[("B", 1), ("X", 1)],
            rate(3.0, 0.0, 0.0),
        )],
        1e-12,
    )
    .unwrap();
    let mut w = KineticWorkspaceV1::new(&m).unwrap();
    for delta in [1.0, 1e-6, 1e-12] {
        let e = m
            .evaluate(
                Temperature::new(500.0),
                &[1.0 + delta, 1.0, 0.0],
                &mut w,
                || false,
            )
            .unwrap();
        assert_eq!(w.reaction_rates()[0], 0.0);
        assert_eq!(w.species_rates()[2], 0.0);
        close(
            w.rate_jacobian()[2],
            3.0 * ((1.0 + delta) - 1.0),
            1e-14 * delta.max(1e-12),
        );
        assert_eq!(e.entropy_production_w_m3_k, 0.0);
        close(
            w.dimensionless_affinities()[0].unwrap(),
            (1.0 + delta).ln(),
            1e-15,
        );
    }
}

#[test]
fn g3_near_equilibrium_progress_and_entropy_retain_the_chemical_driving_force() {
    for scale in [1e-150, 1.0, 1e150] {
        let m = isomers(0.0, rate(scale, 0.0, 0.0));
        let mut w = KineticWorkspaceV1::new(&m).unwrap();
        for direction in [-1.0, 1.0] {
            let a = 1.0 + direction * 1e-12;
            let e = m
                .evaluate(Temperature::new(500.0), &[a, 1.0], &mut w, || false)
                .unwrap();
            close(w.reaction_rates()[0] / scale, a - 1.0, 1e-25);
            assert!(e.entropy_production_w_m3_k > 0.0);
        }
    }
}

#[test]
fn g1_network_jacobian_matches_independent_concentration_and_uv_perturbations() {
    let th = thermo(
        vec![
            source("A", 0.01, 4.0, 0.0, 0.0),
            source("B", 0.01, 5.0, -5000.0, 1.0 - 500.0_f64.ln()),
            source("C", 0.01, 3.0, -8000.0, 500.0_f64.ln() - 1.0),
        ],
        vec![1, 1, 1],
    );
    let m = IdealGasKineticsV1::new(
        th,
        "chain",
        1,
        vec![
            reaction("A-B", &[("A", 1)], &[("B", 1)], rate(2.0, 0.3, 12000.0)),
            reaction("B-C", &[("B", 1)], &[("C", 1)], rate(0.7, -0.1, 5000.0)),
        ],
        1e-12,
    )
    .unwrap();
    let mut w = KineticWorkspaceV1::new(&m).unwrap();
    let t = 700.0;
    let c = [2.0, 1.0, 0.5];
    let e = m
        .evaluate(Temperature::new(t), &c, &mut w, || false)
        .unwrap();
    let jac = w.rate_jacobian().to_vec();
    let dt = w.temperature_derivatives().to_vec();
    let u = w.species_internal_energies().to_vec();
    for i in 0..3 {
        let h = 1e-5;
        let mut plus = c;
        let mut minus = c;
        plus[i] += h;
        minus[i] -= h;
        m.evaluate(Temperature::new(t), &plus, &mut w, || false)
            .unwrap();
        let rp = w.reaction_rates().to_vec();
        m.evaluate(Temperature::new(t), &minus, &mut w, || false)
            .unwrap();
        let rm = w.reaction_rates().to_vec();
        for j in 0..2 {
            close(
                jac[j * 3 + i],
                (rp[j] - rm[j]) / (2.0 * h),
                2e-9 * jac[j * 3 + i].abs().max(1.0),
            );
        }
        let temperature = |amounts: &[f64; 3]| {
            (e.thermodynamics.internal_energy_density_j_m3
                + 5000.0 * amounts[1]
                + 8000.0 * amounts[2])
                / (R * (3.0 * amounts[0] + 4.0 * amounts[1] + 2.0 * amounts[2]))
        };
        m.evaluate(Temperature::new(temperature(&plus)), &plus, &mut w, || {
            false
        })
        .unwrap();
        let rp = w.reaction_rates().to_vec();
        m.evaluate(
            Temperature::new(temperature(&minus)),
            &minus,
            &mut w,
            || false,
        )
        .unwrap();
        let rm = w.reaction_rates().to_vec();
        for j in 0..2 {
            close(
                jac[j * 3 + i] - dt[j] * u[i] / e.thermodynamics.heat_capacity_v_density_j_m3_k,
                (rp[j] - rm[j]) / (2.0 * h),
                2e-9 * (jac[j * 3 + i]
                    - dt[j] * u[i] / e.thermodynamics.heat_capacity_v_density_j_m3_k)
                    .abs()
                    .max(1.0),
            );
        }
    }
}

#[test]
fn g0_invalid_mechanisms_sources_range_and_cancellation_refuse_explicitly() {
    let good = isomers(0.0, rate(1.0, 0.0, 0.0));
    let bad = reaction("unbalanced", &[("A", 1)], &[("B", 2)], rate(1.0, 0.0, 0.0));
    assert!(
        IdealGasKineticsV1::new(good.thermodynamics().clone(), "bad", 1, vec![bad], 1e-12).is_err()
    );
    let th = thermo(
        vec![
            source("A", 0.01, 4.0, 0.0, 0.0),
            source("B", 0.01, 4.0, 0.0, 0.0),
        ],
        vec![1, 2],
    );
    assert!(matches!(
        IdealGasKineticsV1::new(th, "bad-elements", 1, good.reactions().to_vec(), 1e-12),
        Err(KineticsErrorV1::Conservation(_))
    ));
    let mut w = KineticWorkspaceV1::new(&good).unwrap();
    assert!(
        good.evaluate(Temperature::new(500.0), &[-1.0, 2.0], &mut w, || false)
            .is_err()
    );
    assert!(
        good.evaluate(Temperature::new(500.0), &[0.0, 0.0], &mut w, || false)
            .is_err()
    );
    assert_eq!(
        good.evaluate(Temperature::new(500.0), &[1.0, 0.0], &mut w, || true),
        Err(KineticsErrorV1::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        good.evaluate(Temperature::new(500.0), &[1.0, 0.0], &mut w, || {
            polls += 1;
            polls == 4
        }),
        Err(KineticsErrorV1::Cancelled)
    );
    let tiny = isomers(0.0, rate(1e-200, 0.0, 0.0));
    let mut w = KineticWorkspaceV1::new(&tiny).unwrap();
    assert!(matches!(
        tiny.evaluate(Temperature::new(500.0), &[1e-200, 1e-200], &mut w, || false),
        Err(KineticsErrorV1::Arithmetic(_))
    ));
    let disabled = isomers(0.0, rate(0.0, 0.0, 0.0));
    let mut w = KineticWorkspaceV1::new(&disabled).unwrap();
    assert_eq!(
        disabled
            .evaluate(Temperature::new(500.0), &[1.0, 0.0], &mut w, || false)
            .unwrap()
            .entropy_production_w_m3_k,
        0.0
    );
    assert!(
        disabled
            .evaluate(Temperature::new(100.0), &[1.0, 0.0], &mut w, || false)
            .is_err()
    );
}

#[test]
fn g5_canonical_mechanism_and_workspace_reuse_replay_exactly() {
    let m = isomers(7000.0, rate(3.0, 0.4, 6000.0));
    let mut w = KineticWorkspaceV1::new(&m).unwrap();
    let first = m
        .evaluate(Temperature::new(700.0), &[1.5, 0.2], &mut w, || false)
        .unwrap();
    let rates = w.reaction_rates().to_vec();
    let jac = w.rate_jacobian().to_vec();
    m.evaluate(Temperature::new(300.0), &[0.0, 0.5], &mut w, || false)
        .unwrap();
    let repeated = m
        .evaluate(Temperature::new(700.0), &[1.5, 0.2], &mut w, || false)
        .unwrap();
    assert_eq!(first, repeated);
    assert_eq!(rates, w.reaction_rates());
    assert_eq!(jac, w.rate_jacobian());
    assert!(m.conservation().matches(
        m.thermodynamics().elemental_matrix(),
        m.stoichiometry(),
        m.thermodynamics().charge_vector()
    ));
}
