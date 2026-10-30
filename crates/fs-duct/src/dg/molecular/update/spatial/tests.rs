//! Physical weak EOS, entropy/momentum retention, and continued native dynamics.
use super::super::tests::{declaration, excited};
use super::*;
use crate::profile::gas_field::tests::field;
use fs_matdb::QueryPoint;
use fs_material::gas_field::PrimitiveGasState;
use fs_material::visco::GeneralizedMaxwell;
use fs_phs::WallPin;

const RADIUS: f64 = 0.0005;

fn base(field: &AxialGasField, active: bool) -> DgDuct {
    DgDuct::new_with_gas_field(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: field.length_m(),
            }],
        },
        field,
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap()
    .with_spatial_viscothermal(DgWallResolution {
        poles_per_field: 8,
        quadrature_points: 4,
        band_rad_s: [50.0, 2000.0].map(|f| core::f64::consts::TAU * f),
        maximum_relative_error: 0.02,
    })
    .unwrap()
    .with_molecular_relaxation_field(declaration(active))
    .unwrap()
}

fn close(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() <= 2e-9 * scale.max(1e-25),
        "actual {a:e}, reference {b:e}, scale {scale:e}"
    );
}

// Independent physical ideal-gas equation. Frozen compression and axial
// momentum are integrated with a separate 64-point Gauss rule, not the
// production matrix roots or transfer. The original wall rule supplies only
// physical quadrature volumes, and the Bessel producer supplies radial weights.
// At finite wall quadrature, its thermal correction adds to the exact frozen
// pressure mass; that independent mass axis must not be silently lumped.
fn physical(duct: &DgDuct, state: &DgState) -> Vec<f64> {
    let field = duct.gas_field().unwrap();
    let law = duct.molecular_relaxation_law().unwrap();
    let local = duct.mesh.order + 1;
    let m = duct.field_dofs;
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    let (p, q) = duct.from_energy_coords(&state.tilde);
    let resolution = duct.spatial_viscothermal_admission().unwrap().resolution;
    let (mass_nodes, mass_weights) = gauss_legendre_rule(64);
    let (wall_nodes, wall_weights) = gauss_legendre_rule(resolution.quadrature_points);
    let modes = &duct.molecular.as_ref().unwrap().modes;
    let mut compression = vec![0.0; m];
    let mut momentum = vec![0.0; m];
    let mut entropy = Vec::new();
    let mut viscous = Vec::new();
    let mut at = 0;
    let mut viscous_at = 2 * m;
    for (e, element) in duct.elements.iter().enumerate() {
        let jac = 0.5 * element.length;
        for (&xi, &weight) in mass_nodes.iter().zip(&mass_weights) {
            let x = element.x_left + jac * (xi + 1.0);
            let material = law
                .evaluate(
                    &RelaxingMedium::from_thermochemical(&field.state_at(x).unwrap()).unwrap(),
                )
                .unwrap();
            let gas = material.gas();
            let bulk = gas.density * material.frozen_sound_speed_m_s().powi(2);
            let phi = duct.lgl.eval_basis(xi);
            let pressure: f64 = (0..local).map(|i| phi[i] * p[e * local + i]).sum();
            let flow: f64 = (0..local).map(|i| phi[i] * q[e * local + i]).sum();
            for (i, phi) in phi.iter().enumerate() {
                compression[e * local + i] += jac * weight * area * phi * pressure / bulk;
                momentum[e * local + i] += jac * weight * gas.density / area * phi * flow;
            }
        }
        for (&xi, &weight) in wall_nodes.iter().zip(&wall_weights) {
            let x = element.x_left + jac * (xi + 1.0);
            let material = law
                .evaluate(
                    &RelaxingMedium::from_thermochemical(&field.state_at(x).unwrap()).unwrap(),
                )
                .unwrap();
            let gas = material.gas();
            let cp = material.frozen_cv_j_kg_k() + gas.specific_gas_constant;
            let bulk = gas.density * material.frozen_sound_speed_m_s().powi(2);
            let volume = jac * weight * area;
            let phi = duct.lgl.eval_basis(xi);
            let pressure: f64 = (0..local).map(|i| phi[i] * p[e * local + i]).sum();
            let capacities: Vec<_> = material
                .processes()
                .iter()
                .filter(|process| process.heat_capacity_j_kg_k > 0.0)
                .map(|process| process.heat_capacity_j_kg_k)
                .collect();
            let spectral = fs_phs::ZkDiffusive::new(
                RADIUS,
                gas.density,
                gas.dynamic_viscosity,
                gas.gamma,
                gas.prandtl,
                resolution.poles_per_field,
                resolution.band_rad_s,
                resolution.maximum_relative_error,
            )
            .unwrap();
            let mut temperature_sum = 0.0;
            for term in spectral.shunt_terms() {
                let radial = term.gain_per_s / (term.rate_per_s * (gas.gamma - 1.0));
                let mode = &modes[at];
                at += 1;
                let temperature_coordinate = state.tilde[mode.start]
                    + (radial * volume / (gas.density * cp * gas.temperature)).sqrt() * pressure;
                let temperature =
                    temperature_coordinate * (gas.temperature / (volume * gas.density * cp)).sqrt();
                temperature_sum += radial.sqrt() * temperature;
                entropy.push(
                    (gas.density * cp * temperature - radial.sqrt() * pressure) / gas.temperature,
                );
                for (j, &capacity) in capacities.iter().enumerate() {
                    let temperature = state.tilde[mode.start + j + 1]
                        * (gas.temperature / (volume * gas.density * capacity)).sqrt();
                    entropy.push(gas.density * capacity * temperature / gas.temperature);
                }
            }
            let correction =
                pressure / gas.pressure - temperature_sum / gas.temperature - pressure / bulk;
            for (i, phi) in phi.iter().enumerate() {
                compression[e * local + i] += volume * phi * correction;
            }
            for _ in spectral.series_terms() {
                viscous.push(gas.density.sqrt() * state.tilde[viscous_at]);
                viscous_at += 1;
            }
        }
    }
    assert_eq!(at, modes.len());
    assert_eq!(viscous_at, modes[0].start);
    compression.extend(momentum);
    compression.extend(entropy);
    compression.extend(viscous);
    compression
}

#[test]
fn spatial_molecular_update_preserves_independent_weak_eos_entropy_and_momenta() {
    let old = field([280.0, 350.0], [0.85, 0.35], 100_000.0);
    let new = field([410.0, 300.0], [0.85, 0.35], 115_000.0);
    let mut duct = base(&old, true);
    let mut state = excited(&duct);
    let before = state.clone();
    let original = duct.clone();
    let invariants = physical(&duct, &state);
    let report = duct
        .rebind_spatial_molecular_gas_field(&mut state, &new, 0.02, 1.0, &mut || false)
        .unwrap();
    let after = physical(&duct, &state);
    for (a, b) in after.iter().zip(&invariants) {
        close(*a, *b, b.abs().max(1e-20));
    }
    assert_eq!(duct.gas_field(), Some(&new));
    assert_eq!(state.step_count(), before.step_count());
    assert_eq!(state.initial_energy_j(), before.initial_energy_j());
    assert_eq!(
        state.cumulative_dissipated_j(),
        before.cumulative_dissipated_j()
    );
    close(
        report.energy_after_j - report.energy_before_j,
        report.parameter_work_j,
        report.energy_before_j + report.energy_after_j,
    );
    // The spatial map is live: a scalar inlet rescaling cannot represent it.
    let a = &original.elements[0];
    let b = &duct.elements[0];
    let map = conserved_map(&a.mp_sqrt, &b.mp_sqrt, &b.mp_inv_sqrt, 3).unwrap();
    assert!(
        map.iter()
            .enumerate()
            .filter(|(i, _)| i / 3 != i % 3)
            .any(|(_, v)| v.abs() > 1e-5)
    );
    let back = duct
        .rebind_spatial_molecular_gas_field(&mut state, &old, 0.02, 1.0, &mut || false)
        .unwrap();
    for (a, b) in state.tilde.iter().zip(&before.tilde) {
        close(*a, *b, before.stored_energy_j().sqrt());
    }
    close(
        report.parameter_work_j + back.parameter_work_j,
        0.0,
        before.stored_energy_j(),
    );
}

#[test]
fn spatial_molecular_update_keeps_structural_histories_and_continues_native_loss() {
    let old = field([300.0; 2], [0.7; 2], 100_000.0);
    let new = field([390.0, 315.0], [0.7; 2], 100_000.0);
    let pin = WallPin {
        surface_density: 0.4,
        stiffness_per_area: 4e8,
        resistance: 30.0,
    };
    let material = GeneralizedMaxwell::new(2e6, vec![(8e5, 0.001), (3e5, 0.0002)]).unwrap();
    let mut duct = base(&old, true)
        .with_compliant_wall_relaxation(pin, &material)
        .unwrap();
    let mut state = excited(&duct);
    let before = state.clone();
    let motion = duct.compliant_wall_state(&state).unwrap();
    let mode = duct.molecular.as_ref().unwrap().modes.last().unwrap();
    let start = mode.start + mode.rates.len() + 1;
    let original_step = duct.propagator(2e-6, 1e-13).unwrap();
    let report = duct
        .rebind_spatial_molecular_gas_field(&mut state, &new, 0.02, 1.0, &mut || false)
        .unwrap();
    assert_eq!(&state.tilde[start..], &before.tilde[start..]);
    let after = duct.compliant_wall_state(&state).unwrap();
    assert_eq!(motion, after);
    assert!(
        after
            .relaxation_pressure_pa
            .iter()
            .flatten()
            .any(|x| x.abs() > 0.0)
    );
    let checkpoint = state.clone();
    assert!(original_step.step(&mut state, &[]).is_err());
    assert_eq!(state, checkpoint);
    let stepper = duct.propagator(2e-6, 1e-13).unwrap();
    let mut energy = state.stored_energy_j();
    for _ in 0..10 {
        let result = stepper.step(&mut state, &[]).unwrap();
        assert!(state.stored_energy_j() <= energy + 1e-12 * before.initial_energy_j());
        energy = state.stored_energy_j();
        assert_eq!(
            result.ledger.cumulative_parameter_work_j,
            report.parameter_work_j
        );
        close(
            result.ledger.energy_balance_residual_j,
            0.0,
            before.initial_energy_j(),
        );
        assert!(duct.molecular_dissipated_power_w(&state).unwrap() >= 0.0);
        assert!(duct.wall_dissipated_power_w(&state).unwrap() >= 0.0);
    }
}

#[test]
fn spatial_molecular_update_noop_empty_law_and_uniform_transitions() {
    let original = field([310.0, 390.0], [0.7; 2], 100_000.0);
    let uniform = field([350.0; 2], [0.7; 2], 100_000.0);
    let mut duct = base(&original, true);
    let mut state = excited(&duct);
    let before = state.clone();
    let operator = duct.l_tilde.clone();
    let report = duct
        .rebind_spatial_molecular_gas_field(&mut state, &original, 0.02, 0.0, &mut || false)
        .unwrap();
    assert_eq!(state.tilde, before.tilde);
    assert_eq!(report.parameter_work_j, 0.0);
    assert_eq!(duct.l_tilde, operator);
    let physical_before = physical(&duct, &state);
    duct.rebind_spatial_molecular_gas_field(&mut state, &uniform, 0.02, 1.0, &mut || false)
        .unwrap();
    assert!(duct.molecular_relaxation().is_some());
    for (a, b) in physical(&duct, &state).iter().zip(&physical_before) {
        close(*a, *b, b.abs().max(1e-20));
    }
    duct.rebind_spatial_molecular_gas_field(&mut state, &original, 0.02, 1.0, &mut || false)
        .unwrap();
    assert!(duct.molecular_relaxation_field().is_some());
    for (a, b) in state.tilde.iter().zip(&before.tilde) {
        close(*a, *b, before.stored_energy_j().sqrt());
    }

    let mut marked = base(&original, false);
    let mut plain = marked.clone();
    plain.molecular = None;
    let mut a = marked.project_initial_state(|_| 11.0, |_| 2e-7);
    let mut b = a.clone();
    let expected = plain
        .rebind_spatial_viscothermal_gas_field(&mut b, &uniform, 0.02, 1.0, &mut || false)
        .unwrap();
    let actual = marked
        .rebind_spatial_molecular_gas_field(&mut a, &uniform, 0.02, 1.0, &mut || false)
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(a, b);
    assert_eq!(marked.l_tilde, plain.l_tilde);
    assert_eq!(marked.molecular_relaxation_law(), Some(&declaration(false)));
}

#[test]
fn spatial_molecular_update_refusals_and_cancellation_are_atomic() {
    let original = field([300.0, 340.0], [0.8, 0.4], 100_000.0);
    let replacement = field([400.0, 310.0], [0.8, 0.4], 100_000.0);
    let duct = base(&original, true);
    let state = excited(&duct);
    let changed_stations = AxialGasField::new(
        vec![0.0, 0.2, original.length_m()],
        original.static_pressure_pa(),
        vec![
            original.primitive_states()[0].clone(),
            PrimitiveGasState {
                temperature_k: 320.0,
                composition: original.primitive_states()[0].composition.clone(),
            },
            original.primitive_states()[1].clone(),
        ],
        original.model().clone(),
        QueryPoint::new(),
        original.phase(),
    )
    .unwrap();
    for (field, accuracy, allowance) in [
        (field([300.0, 340.0], [0.8, 0.3], 100_000.0), 0.02, 1.0),
        (changed_stations, 0.02, 1.0),
        (replacement.clone(), 0.02, 0.0),
        (replacement.clone(), 0.0, 1.0),
    ] {
        let mut candidate = duct.clone();
        let mut next = state.clone();
        assert!(
            candidate
                .rebind_spatial_molecular_gas_field(
                    &mut next,
                    &field,
                    accuracy,
                    allowance,
                    &mut || false
                )
                .is_err()
        );
        assert_eq!(next, state);
        assert_eq!(candidate.operator_identity(), duct.operator_identity());
        assert_eq!(candidate.gas_field(), duct.gas_field());
    }
    let mut calls = 0;
    let mut candidate = duct.clone();
    let mut next = state.clone();
    assert!(
        candidate
            .rebind_spatial_molecular_gas_field(&mut next, &replacement, 0.02, 1.0, &mut || {
                calls += 1;
                calls > 12
            })
            .is_err()
    );
    assert!(calls > 12);
    assert_eq!(next, state);
    assert_eq!(candidate.operator_identity(), duct.operator_identity());
}
