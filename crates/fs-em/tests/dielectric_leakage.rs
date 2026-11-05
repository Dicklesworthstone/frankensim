//! G1/G3 source-backed field leakage and independent charge/energy oracles.
//! Synthetic material cards verify the admitted EQS law, not measured accuracy.
#![cfg(feature = "resolved-leaky-dielectric")]

mod common;

use common::{config, ends, layered, limits};
use fs_em::electrostatic::{
    ElectrostaticControls, ElectrostaticError,
    dispersion::{AxisymmetricDebyeDielectric, PlanarDebyeDielectric},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::ELECTRICAL_CONDUCTIVITY_DIMS,
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS,
        dielectric::{DebyePoleProperties, leaky::ResolvedLeakyDielectric},
    },
    state_point::MaterialPropertySelection,
};
use fs_qty::Time;
use fs_rep_mesh::{Metric2, TriComplex2};

const BOTH_ENDS: [usize; 4] = [0, 3, 2, 5];
const GROUNDED_LEFT: [usize; 2] = [0, 3];
const AREA_M2: f64 = 0.2;
const VOLUME_M3: f64 = 0.2;

fn law(epsilon: f64, sigma: f64, poles: &[(f64, f64)]) -> ResolvedLeakyDielectric {
    let mut values = vec![
        (
            "epsilon_inf".to_owned(),
            ELECTRIC_PERMITTIVITY_DIMS,
            epsilon,
        ),
        (
            "dc_conductivity".to_owned(),
            ELECTRICAL_CONDUCTIVITY_DIMS,
            sigma,
        ),
    ];
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
    let mut claims = ClaimSet::new();
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
                    source: "synthetic EQS leakage regression".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic leaky dielectric".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    ResolvedLeakyDielectric::resolve(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &names,
        "dc_conductivity",
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
        (actual - expected).abs() <= relative * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}",
    );
}

fn voltages(interface: f64) -> Vec<f64> {
    vec![0.0, interface, 2.0, 0.0, interface, 2.0]
}

fn initial(
    materials: Vec<ResolvedLeakyDielectric>,
    voltage_v: Vec<f64>,
    polarization: Vec<Vec<[f64; 2]>>,
    boundary: &[usize],
) -> PlanarDebyeDielectric {
    PlanarDebyeDielectric::new_leaky(
        layered(0.5, 1.0, 0.2),
        materials,
        voltage_v,
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
fn held_voltage_supplies_ohmic_charge_without_changing_equilibrated_storage() {
    let (epsilon, delta, sigma, dt) = (2e-11, 3e-11, 1e-10, 0.075);
    let material = law(epsilon, sigma, &[(delta, 0.2)]);
    let state = initial(
        vec![material; 4],
        voltages(1.0),
        vec![vec![[-2.0 * delta, 0.0]]; 4],
        &BOTH_ENDS,
    );
    let original_charge = state.nodal_charge_c().to_vec();
    let original_energy = state.stored_energy_j();
    let prepared = state.prepare_step(dt, &mut || false).unwrap();
    let (next, step) = state
        .step_candidate(&prepared, &ends(2.0), controls(), &mut || false)
        .unwrap()
        .into_parts();

    let current = sigma * AREA_M2 * 2.0; // A sigma V/L, L = 1 m.
    let supplied_right = step.nodal_boundary_supplied_charge_c[2].unwrap()
        + step.nodal_boundary_supplied_charge_c[5].unwrap();
    near(supplied_right, current * dt, 3e-9);
    let stored_right = step.nodal_boundary_charge_increment_c[2].unwrap()
        + step.nodal_boundary_charge_increment_c[5].unwrap();
    assert!(stored_right.abs() < 1e-23);
    for (&before, &after) in original_charge.iter().zip(next.nodal_charge_c()) {
        assert!((after - before).abs() < 1e-23);
    }
    near(
        original_charge[2] + original_charge[5],
        (epsilon + delta) * AREA_M2 * 2.0,
        3e-9,
    );
    near(next.stored_energy_j(), original_energy, 3e-9);
    assert!(step.stored_change_j.abs() < 1e-23);
    assert!(step.polarization_heat_j.abs() < 1e-28);
    near(step.joule_heat_j, dt * sigma * 4.0 * VOLUME_M3, 3e-9);
    near(step.heat_j, step.joule_heat_j, 3e-9);
    near(step.supplied_work_j, 2.0 * current * dt, 3e-9);
    near(step.supplied_conduction_work_j, step.joule_heat_j, 3e-9);
    assert!(step.supplied_displacement_work_j.abs() < 1e-23);
    assert!(step.energy_residual_j.abs() < 1e-22);
    for cell in 0..4 {
        near(
            step.conduction_current_density_a_m2[cell][0],
            -2.0 * sigma,
            3e-9,
        );
        assert!(step.conduction_current_density_a_m2[cell][1].abs() < 1e-20);
        assert!(step.polarization_current_density_a_m2[cell][0].abs() < 1e-22);
        near(
            step.cell_joule_heat_j[cell],
            dt * sigma * 4.0 * next.cell_volumes_m3()[cell],
            3e-9,
        );
    }
}

#[test]
fn isolated_maxwell_charge_relaxation_is_second_order_and_converts_storage_to_heat() {
    let (epsilon, sigma, duration) = (2e-11, 1e-10, 0.3);
    let run = |n: usize| {
        let mut state = initial(
            vec![law(epsilon, sigma, &[]); 4],
            voltages(1.0),
            vec![vec![]; 4],
            &GROUNDED_LEFT,
        );
        let original = state.stored_energy_j();
        let original_right = state.nodal_charge_c()[2] + state.nodal_charge_c()[5];
        let prepared = state
            .prepare_step(duration / n as f64, &mut || false)
            .unwrap();
        let mut heat = 0.0;
        for _ in 0..n {
            let (next, step) = state
                .step_candidate(&prepared, &[(0, 0.0), (3, 0.0)], controls(), &mut || false)
                .unwrap()
                .into_parts();
            assert_eq!(step.supplied_work_j, 0.0);
            assert_eq!(step.polarization_heat_j, 0.0);
            assert!(step.joule_heat_j > 0.0);
            assert!(step.maximum_free_charge_residual_c < 1e-22);
            heat += step.joule_heat_j;
            state = next;
        }
        near(state.stored_energy_j() + heat, original, 3e-9);
        let right_charge = state.nodal_charge_c()[2] + state.nodal_charge_c()[5];
        assert!(right_charge < 0.3 * original_right);
        near(right_charge, epsilon * AREA_M2 * state.voltage_v()[2], 3e-9);
        state.voltage_v()[2]
    };
    // For zero Debye poles, epsilon Edot + sigma E = 0.
    let exact = 2.0 * (-duration * sigma / epsilon).exp();
    let coarse = run(40);
    let fine = run(80);
    let ratio = (coarse - exact).abs() / (fine - exact).abs();
    assert!((3.95..4.05).contains(&ratio), "second-order ratio {ratio}");
    near(fine, exact, 5e-5);
}

#[test]
fn maxwell_wagner_interface_charge_evolves_toward_the_conductivity_divider() {
    let (epsilon_left, epsilon_right) = (2e-11, 6e-11);
    let (sigma_left, sigma_right) = (1e-10, 4e-10);
    let left = law(epsilon_left, sigma_left, &[]);
    let right = law(epsilon_right, sigma_right, &[]);
    let mut state = initial(
        vec![left.clone(), left, right.clone(), right],
        voltages(1.5),
        vec![vec![]; 4],
        &BOTH_ENDS,
    );
    let initial_charge = state.nodal_charge_c()[1] + state.nodal_charge_c()[4];
    assert!(initial_charge.abs() < 1e-23);
    let original_energy = state.stored_energy_j();
    let (duration, n) = (0.32, 80);
    let prepared = state
        .prepare_step(duration / n as f64, &mut || false)
        .unwrap();
    let (mut heat, mut work) = (0.0, 0.0);
    for _ in 0..n {
        let before = state.nodal_charge_c()[1] + state.nodal_charge_c()[4];
        let (next, step) = state
            .step_candidate(&prepared, &ends(2.0), controls(), &mut || false)
            .unwrap()
            .into_parts();
        let after = next.nodal_charge_c()[1] + next.nodal_charge_c()[4];
        let outflow = step.nodal_conduction_outflow_c[1] + step.nodal_conduction_outflow_c[4];
        assert!((after - before + outflow).abs() < 1e-22);
        assert_eq!(step.polarization_heat_j, 0.0);
        assert!(step.joule_heat_j > 0.0);
        heat += step.joule_heat_j;
        work += step.supplied_work_j;
        state = next;
    }
    // Equal-width layers obey:
    // (epsL+epsR) Vdot + (sigmaL+sigmaR) V = 2 sigmaR.
    let tau = (epsilon_left + epsilon_right) / (sigma_left + sigma_right);
    let steady = 2.0 * sigma_right / (sigma_left + sigma_right);
    near(tau, 0.16, 1e-14);
    near(steady, 1.6, 1e-14);
    let expected = steady + (1.5 - steady) * (-duration / tau).exp();
    assert!((state.voltage_v()[1] - expected).abs() < 2e-6);
    near(state.voltage_v()[4], state.voltage_v()[1], 3e-9);
    let interface_charge = state.nodal_charge_c()[1] + state.nodal_charge_c()[4];
    assert!(interface_charge > 2e-12, "interface free charge was frozen");
    let expected_charge = 2.0 * AREA_M2 * (epsilon_left + epsilon_right) * (expected - 1.5);
    near(interface_charge, expected_charge, 2e-5);
    near(
        interface_charge,
        AREA_M2 * (state.displacement_c_m2()[2][0] - state.displacement_c_m2()[0][0]),
        3e-9,
    );
    near(state.stored_energy_j() + heat, original_energy + work, 3e-9);
}

// Independent exponential of the two physical ODEs, using x = P/delta:
// eps Edot + delta xdot + sigma E = 0; tau xdot + x = E.
// Initial state E=-2 V/m, x=0. No field solver or midpoint recurrence is used.
fn open_debye_ohmic_oracle(time: f64, epsilon: f64, delta: f64, tau: f64, sigma: f64) -> [f64; 2] {
    let a = -(sigma + delta / tau) / epsilon;
    let b = delta / (tau * epsilon);
    let c = 1.0 / tau;
    let d = -1.0 / tau;
    let half_trace = 0.5 * (a + d);
    let gap = ((0.5 * (a - d)).powi(2) + b * c).sqrt();
    let factor = (half_trace * time).exp();
    let cosh = (gap * time).cosh();
    let sinh_over_gap = (gap * time).sinh() / gap;
    [
        factor * (-2.0 * cosh - 2.0 * (a - half_trace) * sinh_over_gap),
        delta * factor * (-2.0 * c * sinh_over_gap),
    ]
}

#[test]
fn combined_debye_and_ohmic_relaxation_matches_independent_two_state_exponential() {
    let (epsilon, delta, tau, sigma, duration) = (2e-11, 3e-11, 0.2, 1e-10, 0.4);
    let run = |n: usize| {
        let mut state = initial(
            vec![law(epsilon, sigma, &[(delta, tau)]); 4],
            voltages(1.0),
            vec![vec![[0.0; 2]]; 4],
            &GROUNDED_LEFT,
        );
        let original = state.stored_energy_j();
        let prepared = state
            .prepare_step(duration / n as f64, &mut || false)
            .unwrap();
        let (mut polarization_heat, mut joule_heat, mut total_heat) = (0.0, 0.0, 0.0);
        for _ in 0..n {
            let (next, step) = state
                .step_candidate(&prepared, &[(0, 0.0), (3, 0.0)], controls(), &mut || false)
                .unwrap()
                .into_parts();
            assert_eq!(step.supplied_work_j, 0.0);
            assert!(step.polarization_heat_j > 0.0);
            assert!(step.joule_heat_j > 0.0);
            polarization_heat += step.polarization_heat_j;
            joule_heat += step.joule_heat_j;
            total_heat += step.heat_j;
            state = next;
        }
        near(total_heat, polarization_heat + joule_heat, 3e-12);
        near(
            state.stored_energy_j() + polarization_heat + joule_heat,
            original,
            3e-9,
        );
        let e = state.electric_field_v_m()[0][0];
        let p = state.polarization_c_m2()[0][0][0];
        near(
            state.nodal_charge_c()[2] + state.nodal_charge_c()[5],
            -AREA_M2 * (epsilon * e + p),
            3e-9,
        );
        [e, p]
    };
    let exact = open_debye_ohmic_oracle(duration, epsilon, delta, tau, sigma);
    let coarse = run(80);
    let fine = run(160);
    for i in 0..2 {
        let ratio = (coarse[i] - exact[i]).abs() / (fine[i] - exact[i]).abs();
        assert!(
            (3.95..4.05).contains(&ratio),
            "state {i}: second-order ratio {ratio}"
        );
    }
    near(fine[0], exact[0], 1.1e-5);
    near(fine[1], exact[1], 3e-6);
}

#[test]
fn sourced_zero_leakage_parity_gauge_identity_and_transactional_refusal() {
    let zero = law(2e-11, 0.0, &[(3e-11, 0.2)]);
    let original_bundle = zero.material().clone();
    let bulk = zero.clone().into_parts().0;
    let leaky = initial(
        vec![zero; 4],
        voltages(1.0),
        vec![vec![[0.0; 2]]; 4],
        &BOTH_ENDS,
    );
    let pure = PlanarDebyeDielectric::new(
        layered(0.5, 1.0, 0.2),
        vec![bulk; 4],
        voltages(1.0),
        vec![vec![[0.0; 2]]; 4],
        &BOTH_ENDS,
        0.0,
        limits(),
        controls(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(leaky.materials(), pure.materials());
    for cell in 0..4 {
        assert_eq!(leaky.materials()[cell].material(), &original_bundle);
        assert_eq!(
            leaky.conductivity_properties()[cell].as_deref(),
            Some("dc_conductivity")
        );
        assert_eq!(pure.conductivity_properties()[cell], None);
        assert_eq!(leaky.conductivity_s_m()[cell], 0.0);
        let key = PropertyKey::new("dc_conductivity", ELECTRICAL_CONDUCTIVITY_DIMS);
        assert_eq!(
            leaky.materials()[cell].material().property_by_key(&key),
            original_bundle.property_by_key(&key),
        );
    }
    let zero_prepared = leaky.prepare_step(0.05, &mut || false).unwrap();
    let pure_prepared = pure.prepare_step(0.05, &mut || false).unwrap();
    let zero_step = leaky
        .step_candidate(&zero_prepared, &ends(2.25), controls(), &mut || false)
        .unwrap();
    let pure_step = pure
        .step_candidate(&pure_prepared, &ends(2.25), controls(), &mut || false)
        .unwrap();
    assert_eq!(
        zero_step.field().relative_voltage_v(),
        pure_step.field().relative_voltage_v()
    );
    assert_eq!(
        zero_step.field().polarization_c_m2(),
        pure_step.field().polarization_c_m2()
    );
    assert_eq!(zero_step.step().heat_j, pure_step.step().heat_j);
    assert_eq!(
        zero_step.step().supplied_work_j,
        pure_step.step().supplied_work_j
    );
    assert_eq!(zero_step.step().joule_heat_j, 0.0);

    let material = law(2e-11, 1e-10, &[(3e-11, 0.2)]);
    let source = material.material().clone();
    let make = |shift: f64| {
        initial(
            vec![material.clone(); 4],
            voltages(1.0).iter().map(|v| v + shift).collect(),
            vec![vec![[0.0; 2]]; 4],
            &BOTH_ENDS,
        )
    };
    let state = make(0.0);
    let shift = 2.0_f64.powi(40);
    let shifted = make(shift);
    let prepared = state.prepare_step(0.02, &mut || false).unwrap();
    let shifted_prepared = shifted.prepare_step(0.02, &mut || false).unwrap();
    let endpoint = [(0, 0.25), (3, 0.25), (2, 2.5), (5, 2.5)];
    let shifted_endpoint: Vec<_> = endpoint.iter().map(|&(v, e)| (v, e + shift)).collect();
    let before_charge = state.nodal_charge_c().to_vec();
    assert!(
        state
            .step_candidate(&shifted_prepared, &endpoint, controls(), &mut || false)
            .is_err()
    );
    assert!(state.prepare_step(0.0, &mut || false).is_err());
    let mut polls = 0;
    let cancelled = state.step_candidate(&prepared, &endpoint, controls(), &mut || {
        polls += 1;
        polls > 20
    });
    assert!(matches!(cancelled, Err(ElectrostaticError::Cancelled)));
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.nodal_charge_c(), before_charge);
    assert_eq!(state.polarization_c_m2(), vec![vec![[0.0; 2]]; 4]);
    let candidate = state
        .step_candidate(&prepared, &endpoint, controls(), &mut || false)
        .unwrap();
    let gauge = shifted
        .step_candidate(
            &shifted_prepared,
            &shifted_endpoint,
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        candidate.field().relative_voltage_v(),
        gauge.field().relative_voltage_v()
    );
    assert_eq!(
        candidate.field().electric_field_v_m(),
        gauge.field().electric_field_v_m()
    );
    assert_eq!(candidate.step().heat_j, gauge.step().heat_j);
    assert_eq!(
        candidate.step().supplied_work_j,
        gauge.step().supplied_work_j
    );
    assert_eq!(candidate.field().materials()[0].material(), &source);
}

fn axial_mesh(angle: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "axial leaky dielectric fixture",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [0.0, 0.2],
            [1.0, 0.2],
            [2.0, 0.2],
            [0.0, 0.4],
            [1.0, 0.4],
            [2.0, 0.4],
        ],
        vec![
            [0, 1, 4],
            [0, 4, 3],
            [1, 2, 5],
            [1, 5, 4],
            [3, 4, 7],
            [3, 7, 6],
            [4, 5, 8],
            [4, 8, 7],
        ],
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap()
}

#[test]
fn axisymmetric_axial_leakage_uses_physical_volume_span_and_voltage_traces() {
    let (epsilon, delta, tau, sigma, dt) = (2e-11, 3e-11, 0.2, 1e-10, 0.04);
    let material = law(epsilon, sigma, &[(delta, tau)]);
    let voltage = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0];
    let electrodes = [0, 1, 2, 6, 7, 8];
    let endpoint = [(0, 0.0), (1, 0.0), (2, 0.0), (6, 2.0), (7, 2.0), (8, 2.0)];
    let mut scaled_results = Vec::new();
    for angle in [core::f64::consts::TAU, 0.7] {
        let state = AxisymmetricDebyeDielectric::new_leaky(
            axial_mesh(angle),
            vec![material.clone(); 8],
            voltage.clone(),
            vec![vec![[0.0; 2]]; 8],
            &electrodes,
            0.0,
            limits(),
            controls(),
            &mut || false,
        )
        .unwrap();
        let physical_volume = 0.5 * angle * 2.0_f64.powi(2) * 0.4;
        let electrode_area = 0.5 * angle * 2.0_f64.powi(2);
        near(state.cell_volumes_m3().iter().sum(), physical_volume, 3e-12);
        near(
            state.stored_energy_j(),
            0.5 * epsilon * 25.0 * physical_volume,
            3e-9,
        );
        let prepared = state.prepare_step(dt, &mut || false).unwrap();
        let (next, step) = state
            .step_candidate(&prepared, &endpoint, controls(), &mut || false)
            .unwrap()
            .into_parts();
        let b = dt / (2.0 * tau + dt);
        let polarization_z = 2.0 * b * delta * -5.0;
        let polarization_heat_density = tau * polarization_z.powi(2) / (delta * dt);
        let joule_heat_density = dt * sigma * 25.0;
        near(
            step.polarization_heat_j,
            physical_volume * polarization_heat_density,
            3e-9,
        );
        near(
            step.joule_heat_j,
            physical_volume * joule_heat_density,
            3e-9,
        );
        near(
            step.heat_j,
            physical_volume * (polarization_heat_density + joule_heat_density),
            3e-9,
        );
        near(
            step.stored_change_j,
            physical_volume * polarization_z.powi(2) / (2.0 * delta),
            3e-9,
        );
        let supplied_top: f64 = [6, 7, 8]
            .iter()
            .map(|&v| step.nodal_boundary_supplied_charge_c[v].unwrap())
            .sum();
        near(
            supplied_top,
            electrode_area * (-polarization_z + dt * sigma * 5.0),
            3e-9,
        );
        near(step.supplied_work_j, 2.0 * supplied_top, 3e-9);
        assert!(step.energy_residual_j.abs() < 1e-21);
        for cell in 0..8 {
            assert!(next.electric_field_v_m()[cell][0].abs() < 1e-11);
            near(next.electric_field_v_m()[cell][1], -5.0, 3e-9);
            near(next.polarization_c_m2()[cell][0][1], polarization_z, 3e-9);
            near(
                step.conduction_current_density_a_m2[cell][1],
                -5.0 * sigma,
                3e-9,
            );
            near(
                step.cell_joule_heat_j[cell],
                next.cell_volumes_m3()[cell] * joule_heat_density,
                3e-9,
            );
        }
        scaled_results.push([
            next.stored_energy_j() / angle,
            step.heat_j / angle,
            supplied_top / angle,
        ]);
    }
    for i in 0..3 {
        near(scaled_results[0][i], scaled_results[1][i], 3e-9);
    }

    assert!(
        AxisymmetricDebyeDielectric::new_leaky(
            layered(0.5, 1.0, 0.2),
            vec![material.clone(); 4],
            voltages(1.0),
            vec![vec![[0.0; 2]]; 4],
            &BOTH_ENDS,
            0.0,
            limits(),
            controls(),
            &mut || false,
        )
        .is_err()
    );
    assert!(
        PlanarDebyeDielectric::new_leaky(
            axial_mesh(core::f64::consts::TAU),
            vec![material.clone(); 8],
            voltage.clone(),
            vec![vec![[0.0; 2]]; 8],
            &electrodes,
            0.0,
            limits(),
            controls(),
            &mut || false,
        )
        .is_err()
    );
    // The symmetry axis sweeps zero physical electrode area.
    assert!(
        AxisymmetricDebyeDielectric::new_leaky(
            axial_mesh(core::f64::consts::TAU),
            vec![material; 8],
            voltage,
            vec![vec![[0.0; 2]]; 8],
            &[0, 3, 6],
            0.0,
            limits(),
            controls(),
            &mut || false,
        )
        .is_err()
    );
}
