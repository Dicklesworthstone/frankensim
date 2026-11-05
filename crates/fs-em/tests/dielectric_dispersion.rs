//! G1/G3/G4 spatial Debye memory: analytic dynamics, physical power balance,
//! heterogeneous fields and transactional refusal. Synthetic source fixtures.
mod common;

use common::{config, ends, layered, limits};
use fs_em::electrostatic::{
    ElectrostaticControls, ElectrostaticError, dispersion::PlanarDebyeDielectric,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS,
        dielectric::{DebyePoleProperties, bulk::ResolvedDebyeMaterial},
    },
    state_point::MaterialPropertySelection,
};
use fs_qty::Time;
use fs_rep_mesh::{Metric2, TriComplex2};

fn law(epsilon: f64, poles: &[(f64, f64)]) -> ResolvedDebyeMaterial {
    let mut claims = ClaimSet::new();
    let mut values = vec![(
        "epsilon_inf".to_owned(),
        ELECTRIC_PERMITTIVITY_DIMS,
        epsilon,
    )];
    let names: Vec<_> = poles
        .iter()
        .enumerate()
        .map(|(i, &(delta, tau))| {
            let permittivity_increment = format!("delta_{i}");
            let relaxation_time = format!("tau_{i}");
            values.push((
                permittivity_increment.clone(),
                ELECTRIC_PERMITTIVITY_DIMS,
                delta,
            ));
            values.push((relaxation_time.clone(), Time::DIMS, tau));
            DebyePoleProperties {
                permittivity_increment,
                relaxation_time,
            }
        })
        .collect();
    for (name, dims, value) in values {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic bulk Debye regression".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    ResolvedDebyeMaterial::resolve(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &names,
    )
    .unwrap()
}

fn controls() -> ElectrostaticControls {
    ElectrostaticControls {
        linear: config(),
        charge_tolerance_c: 1e-17,
        energy_tolerance_j: 1e-17,
    }
}
fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-22),
        "{actual:e} != {expected:e}"
    );
}
fn initial(materials: Vec<ResolvedDebyeMaterial>, boundary: &[usize]) -> PlanarDebyeDielectric {
    let polarization = materials
        .iter()
        .map(|law| vec![[0.0; 2]; law.poles().len()])
        .collect();
    PlanarDebyeDielectric::new(
        layered(0.5, 1.0, 0.2),
        materials,
        vec![0.0, 1.0, 2.0, 0.0, 1.0, 2.0],
        polarization,
        boundary,
        0.0,
        limits(),
        controls(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn held_voltage_polarization_converges_at_second_order_with_separate_heat_and_storage() {
    fn run(n: usize) -> (f64, f64, f64, f64) {
        let material = law(2e-11, &[(3e-11, 0.2)]);
        let mut state = initial(vec![material; 4], &[0, 3, 2, 5]);
        let prepared = state.prepare_step(0.2 / n as f64, &mut || false).unwrap();
        let original = state.stored_energy_j();
        let (mut heat, mut work) = (0.0, 0.0);
        for _ in 0..n {
            let (next, step) = state
                .step_candidate(&prepared, &ends(2.0), controls(), &mut || false)
                .unwrap()
                .into_parts();
            assert!(step.heat_j > 0.0);
            assert!(step.energy_residual_j.abs() < 1e-23);
            near(
                step.supplied_work_j,
                2.0 * (step.nodal_boundary_charge_increment_c[2].unwrap()
                    + step.nodal_boundary_charge_increment_c[5].unwrap()),
                2e-12,
            );
            heat += step.heat_j;
            work += step.supplied_work_j;
            state = next;
        }
        for e in state.electric_field_v_m() {
            near(e[0], -2.0, 1e-12);
            assert!(e[1].abs() < 1e-12);
        }
        near(state.stored_energy_j() - original + heat, work, 1e-12);
        (
            state.polarization_c_m2()[0][0][0],
            heat,
            work,
            state.stored_energy_j(),
        )
    }
    let coarse = run(20);
    let fine = run(40);
    let expected = -6e-11 * (1.0 - (-1.0_f64).exp());
    let ratio = (coarse.0 - expected).abs() / (fine.0 - expected).abs();
    assert!((3.95..4.05).contains(&ratio), "second-order ratio {ratio}");
    near(fine.0, expected, 4e-5);
    // Independent continuous pole oracles, volume 0.2 m3 and E²=4.
    near(
        fine.1,
        0.5 * 0.2 * 3e-11 * 4.0 * (1.0 - (-2.0_f64).exp()),
        5e-5,
    );
    near(fine.2, 0.2 * 3e-11 * 4.0 * (1.0 - (-1.0_f64).exp()), 5e-5);
}

#[test]
fn isolated_precharge_relaxes_voltage_without_inventing_source_work_or_losing_charge() {
    let material = law(2e-11, &[(3e-11, 0.2)]);
    let mut state = initial(vec![material; 4], &[0, 3]);
    let charge = state.nodal_charge_c().to_vec();
    let original = state.stored_energy_j();
    let prepared = state.prepare_step(0.001, &mut || false).unwrap();
    let mut heat = 0.0;
    for _ in 0..160 {
        let (next, step) = state
            .step_candidate(&prepared, &[(0, 0.0), (3, 0.0)], controls(), &mut || false)
            .unwrap()
            .into_parts();
        assert_eq!(step.supplied_work_j, 0.0);
        assert!(step.heat_j >= 0.0);
        heat += step.heat_j;
        state = next;
    }
    // D stays fixed. Effective relaxation time is tau*epsilon_inf/epsilon_static.
    let expected = 0.8 + 1.2 * (-2.0_f64).exp();
    near(state.voltage_v()[2], expected, 7e-6);
    for &v in &[1, 2, 4, 5] {
        assert!((state.nodal_charge_c()[v] - charge[v]).abs() < 1e-23);
    }
    near(state.stored_energy_j() + heat, original, 2e-12);
}

#[test]
fn layered_poles_redistribute_fields_and_keep_free_interface_charge() {
    let left = law(2e-11, &[(6e-11, 0.1)]);
    let right = law(2e-11, &[]);
    let mut state = initial(
        vec![left.clone(), left, right.clone(), right],
        &[0, 3, 2, 5],
    );
    let free = [state.nodal_charge_c()[1], state.nodal_charge_c()[4]];
    let prepared = state.prepare_step(0.005, &mut || false).unwrap();
    for _ in 0..400 {
        state = state
            .step_candidate(&prepared, &ends(2.0), controls(), &mut || false)
            .unwrap()
            .into_parts()
            .0;
    }
    // Static series slabs: epsilon_left=8e-11, epsilon_right=2e-11.
    for cell in 0..4 {
        near(
            state.electric_field_v_m()[cell][0],
            if cell < 2 { -0.8 } else { -3.2 },
            1e-9,
        );
        near(state.displacement_c_m2()[cell][0], -6.4e-11, 1e-9);
    }
    assert!((state.nodal_charge_c()[1] - free[0]).abs() < 1e-23);
    assert!((state.nodal_charge_c()[4] - free[1]).abs() < 1e-23);
    assert_eq!(state.materials()[2].poles().len(), 0);
}

#[test]
fn multipole_prehistory_changes_field_without_becoming_ohmic_leakage() {
    let material = law(2e-11, &[(3e-11, 0.05), (5e-11, 0.3)]);
    let mut state = initial(vec![material; 4], &[0, 3, 2, 5]);
    let prepared = state.prepare_step(0.002, &mut || false).unwrap();
    for _ in 0..100 {
        state = state
            .step_candidate(&prepared, &ends(2.0), controls(), &mut || false)
            .unwrap()
            .into_parts()
            .0;
    }
    let p = state.polarization_c_m2()[0].clone();
    near(p[0][0], -6e-11 * (1.0 - (-4.0_f64).exp()), 2e-5);
    near(p[1][0], -1e-10 * (1.0 - (-2.0_f64 / 3.0).exp()), 3e-6);
    // A shorted, previously polarized cell dissipates its remaining memory.
    let first = state
        .step_candidate(&prepared, &ends(0.0), controls(), &mut || false)
        .unwrap();
    let next = first.into_parts().0;
    let relaxation = next
        .step_candidate(&prepared, &ends(0.0), controls(), &mut || false)
        .unwrap();
    assert!(relaxation.step().heat_j > 0.0);
    assert_eq!(relaxation.step().supplied_work_j, 0.0);
    assert!(relaxation.field().stored_energy_j() < next.stored_energy_j());
}

#[test]
fn zero_pole_limit_and_independent_component_gauges_preserve_physics() {
    let mesh = TriComplex2::from_indexed_triangles(
        "disconnected-dielectric",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [3.0, 0.0],
            [3.0, 1.0],
            [2.0, 1.0],
        ],
        vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
        Metric2::planar(0.2).unwrap(),
    )
    .unwrap();
    let initial_v = vec![
        1e9,
        1e9 + 2.0,
        1e9 + 2.0,
        1e9,
        -1e8,
        -1e8 + 4.0,
        -1e8 + 4.0,
        -1e8,
    ];
    let state = PlanarDebyeDielectric::new(
        mesh,
        vec![law(2e-11, &[]); 4],
        initial_v,
        vec![vec![]; 4],
        &[0, 1, 2, 3, 4, 5, 6, 7],
        0.0,
        limits(),
        controls(),
        &mut || false,
    )
    .unwrap();
    let prepared = state.prepare_step(0.1, &mut || false).unwrap();
    let boundary = [
        (0, 2e9),
        (1, 2e9 + 3.0),
        (2, 2e9 + 3.0),
        (3, 2e9),
        (4, -3e8),
        (5, -3e8 + 5.0),
        (6, -3e8 + 5.0),
        (7, -3e8),
    ];
    let candidate = state
        .step_candidate(&prepared, &boundary, controls(), &mut || false)
        .unwrap();
    assert_eq!(candidate.step().heat_j, 0.0);
    near(
        candidate.step().supplied_work_j,
        0.5 * 0.2 * 2e-11 * ((9.0 - 4.0) + (25.0 - 16.0)),
        1e-12,
    );
    near(candidate.field().electric_field_v_m()[0][0], -3.0, 1e-12);
    near(candidate.field().electric_field_v_m()[2][0], -5.0, 1e-12);
}

#[test]
fn prepared_operator_normalizes_cancellation_inside_native_hodge_assembly() {
    let material = law(2e-11, &[(3e-11, 0.2)]);
    let state = initial(vec![material; 4], &[0, 3, 2, 5]);
    let mut total_polls = 0;
    state
        .prepare_step(0.01, &mut || {
            total_polls += 1;
            false
        })
        .unwrap();
    // Exercise every cancellation seam in this small prepared operator,
    // including the native Hodge's coefficient/cell/sparse-assembly polls.
    // No hard-coded poll number can silently stop exercising that owner when
    // preparation gains or loses an earlier cancellation check.
    for stop in 1..=total_polls {
        let mut polls = 0;
        let result = state.prepare_step(0.01, &mut || {
            polls += 1;
            polls == stop
        });
        assert!(
            matches!(result, Err(ElectrostaticError::Cancelled)),
            "preparation cancellation at poll {stop} was not normalized"
        );
    }
    let prepared = state.prepare_step(0.01, &mut || false).unwrap();
    assert_eq!(prepared.time_step_s(), 0.01);
    assert_eq!(state.time_s(), 0.0);
}

#[test]
fn cancellation_failed_budget_and_foreign_operator_leave_accepted_history_unchanged() {
    let material = law(2e-11, &[(3e-11, 0.2)]);
    let state = initial(vec![material.clone(); 4], &[0, 3, 2, 5]);
    let foreign = initial(vec![material; 4], &[0, 3, 2, 5]);
    let prepared = state.prepare_step(0.01, &mut || false).unwrap();
    let wrong = foreign.prepare_step(0.01, &mut || false).unwrap();
    assert!(
        state
            .step_candidate(&wrong, &ends(2.0), controls(), &mut || false)
            .is_err()
    );
    let before = state.polarization_c_m2().to_vec();
    let voltages = state.voltage_v().to_vec();
    let mut polls = 0;
    let failed = state.step_candidate(&prepared, &ends(2.0), controls(), &mut || {
        polls += 1;
        polls > 20
    });
    assert!(matches!(failed, Err(ElectrostaticError::Cancelled)));
    assert!(
        state
            .step_candidate(&prepared, &[(0, 0.0), (3, 0.0)], controls(), &mut || false)
            .is_err()
    );
    let mut invalid = controls();
    invalid.energy_tolerance_j = f64::NAN;
    assert!(
        state
            .step_candidate(&prepared, &ends(2.0), invalid, &mut || false)
            .is_err()
    );
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.polarization_c_m2(), before);
    assert_eq!(state.voltage_v(), voltages);
    let a = state
        .step_candidate(&prepared, &ends(2.0), controls(), &mut || false)
        .unwrap();
    let b = state
        .step_candidate(&prepared, &ends(2.0), controls(), &mut || false)
        .unwrap();
    assert_eq!(a.field().voltage_v(), b.field().voltage_v());
    assert_eq!(a.field().polarization_c_m2(), b.field().polarization_c_m2());
    assert_eq!(a.step().heat_j.to_bits(), b.step().heat_j.to_bits());
}
