//! Manufactured material data -> physical frame ports -> original coupled motion.
//! Continuous references use geometry, consistent mass and constitutive data;
//! none is derived from the projected damping matrix or propagated endpoint.

use fs_couple::frame_vibration::{FrameVibrationBasis, FrameVibrationBudget, FrameVibrationError};
use fs_couple::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalConnection, ModalCouplingConfig, ModalCouplingError,
    PreparedModalBudget,
};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::{Density, DynViscosity, Pressure};
use fs_solid::frame2d::{Frame2d, FrameMassFormulation, FrameViscousSection, Geometry};

const DENSITY: f64 = 3.0;
const YOUNG: f64 = 1200.0;
const RADIUS: f64 = 0.1;
const RATE: u32 = 256;
const BAND: (f64, f64) = (0.0, 4096.0);

fn area() -> f64 {
    core::f64::consts::PI * RADIUS * RADIUS
}

fn source(eta: f64, band: (f64, f64)) -> FrameViscousSection {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [
        ("density", DENSITY, Density::DIMS),
        ("young_modulus", YOUNG, Pressure::DIMS),
        ("axial_eta", eta, DynViscosity::DIMS),
        ("bending_eta", 0.0, DynViscosity::DIMS),
    ] {
        let mut validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
        if name.ends_with("_eta") {
            validity = validity.with("omega", band.0, band.1);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured coupled-frame material".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured solid".into(),
            phase: "solid".into(),
            process: "linear viscosity reference".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with("omega", 0.5 * (band.0 + band.1))
        .unwrap();
    let material = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let section = ConcentricElasticSection::try_new(
        0.0,
        vec![ConcentricMaterialLayer {
            outer_radius_m: RADIUS,
            material,
        }],
    )
    .unwrap();
    FrameViscousSection::try_new(section, "axial_eta", "bending_eta").unwrap()
}

fn basis(viscosities: &[f64], retained: usize, band: (f64, f64)) -> FrameVibrationBasis {
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    for node in 0..=viscosities.len() {
        frame.add_node(node as f64 / viscosities.len() as f64, 0.0);
        frame.fix(
            node,
            if node == 0 {
                [true; 3]
            } else {
                [false, true, true]
            },
        );
    }
    for (i, &eta) in viscosities.iter().enumerate() {
        frame
            .add_material_viscous_element(i, i + 1, 5, source(eta, band))
            .unwrap();
    }
    FrameVibrationBasis::from_frame(
        &frame,
        retained,
        FrameVibrationBudget {
            max_frame_dofs: 12,
            max_dense_entries: 144,
        },
        || false,
    )
    .unwrap()
}

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 4,
        max_connections: 8,
        max_setup_terms: 8192,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 10.0,
        maximum_abs_pressure_pa: 100.0,
        maximum_abs_connection_force_n: 1000.0,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-10,
    }
}

fn prepared() -> PreparedModalBudget {
    PreparedModalBudget {
        max_modes: 4,
        max_dense_bytes: 8 * 1024 * 1024,
        max_multiply_terms: 1_000_000_000,
    }
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-8 * scale,
        "actual={actual:.16e}, expected={expected:.16e}, scale={scale:.3e}"
    );
}

fn physical(
    basis: &FrameVibrationBasis,
    system: &CoupledModalSystem,
    component: usize,
) -> (Vec<f64>, Vec<f64>) {
    let (mut u, mut v) = (vec![0.0; basis.ndof()], vec![0.0; basis.ndof()]);
    basis
        .reconstruct_into(
            system.components()[component].states(),
            &mut u,
            &mut v,
            || false,
        )
        .unwrap();
    (u, v)
}

fn scalar_state(basis: &FrameVibrationBasis, q: f64, v: f64) -> [ModalAcousticState; 1] {
    let shape = basis.modes()[0].displacement[3];
    [ModalAcousticState {
        displacement_m_sqrt_kg: q / shape,
        velocity_m_sqrt_kg_per_s: v / shape,
    }]
}

fn free_oscillator(q0: f64, v0: f64, omega2: f64, decay: f64, t: f64) -> (f64, f64) {
    let wd = (omega2 - decay * decay).sqrt();
    let (sin, cos) = (wd * t).sin_cos();
    let envelope = (-decay * t).exp();
    (
        envelope * (q0 * cos + (v0 + decay * q0) / wd * sin),
        envelope * (v0 * cos - (decay * v0 + omega2 * q0) / wd * sin),
    )
}

#[test]
fn data_only_viscosity_changes_actual_ordinary_and_prepared_motion_and_loss() {
    let (mass, k) = (DENSITY * area() / 3.0, YOUNG * area());
    let omega = (k / mass).sqrt();
    let (q0, v0) = (0.004, 0.03);
    let authored_zeta = 0.02;
    let initial = 0.5 * k * q0 * q0 + 0.5 * mass * v0 * v0;
    let mut final_energy = Vec::new();
    for eta in [2.0, 8.0] {
        for preparation in [None, Some(prepared())] {
            let basis = basis(&[eta], 1, BAND);
            let mut component = basis
                .coupled_component(
                    RATE,
                    &[authored_zeta],
                    ModalAcousticTimeBudget::audible_reference(),
                    18,
                    || false,
                )
                .unwrap();
            let initial_state = scalar_state(&basis, q0, v0);
            component.restore_states(&initial_state).unwrap();
            assert_eq!(component.states(), &initial_state);
            assert_eq!(component.material_connections().len(), 3);
            let (basis, mut system) = component
                .into_system(config(), preparation, &CancelGate::new_clock_free())
                .unwrap();
            let cm = eta * area();
            let ca = 2.0 * authored_zeta * omega * mass;
            let decay = (cm + ca) / (2.0 * mass);
            let (mut source_loss, mut authored_loss) = (0.0, 0.0);
            for sample in 1..=64 {
                let report = system.step(&[0.0]).unwrap().clone();
                let (q, v) =
                    free_oscillator(q0, v0, k / mass, decay, f64::from(sample) / f64::from(RATE));
                let (u, velocity) = physical(&basis, &system, 0);
                close(u[3], q, 0.01);
                close(velocity[3], v, 0.4);
                source_loss += report.connection_dissipation_j;
                authored_loss += report.component_dissipation_j;
                let lost = initial - 0.5 * k * q * q - 0.5 * mass * v * v;
                close(source_loss, lost * cm / (cm + ca), initial);
                close(authored_loss, lost * ca / (cm + ca), initial);
                assert_eq!(report.external_work_j, 0.0);
                assert_eq!(report.material_dissipation_j, 0.0);
                assert_eq!(report.observer_pressure_pa, 0.0);
                assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);
            }
            close(
                system.total_energy_j().unwrap() + source_loss + authored_loss,
                initial,
                initial,
            );
            final_energy.push(system.total_energy_j().unwrap());
        }
    }
    close(final_energy[0], final_energy[1], initial);
    close(final_energy[2], final_energy[3], initial);
    assert!(final_energy[2] < final_energy[0]);
}

fn chain_matrices() -> ([[f64; 2]; 2], [[f64; 2]; 2], [[f64; 2]; 2]) {
    let member_mass = DENSITY * area() * 0.5;
    let k = YOUNG * area() / 0.5;
    let (c1, c2) = (area() / 0.5, 5.0 * area() / 0.5);
    (
        [
            [2.0 * member_mass / 3.0, member_mass / 6.0],
            [member_mass / 6.0, member_mass / 3.0],
        ],
        [[2.0 * k, -k], [-k, k]],
        [[c1 + c2, -c2], [-c2, c2]],
    )
}

fn quadratic(matrix: &[[f64; 2]; 2], x: [f64; 2]) -> f64 {
    x[0] * (matrix[0][0] * x[0] + matrix[0][1] * x[1])
        + x[1] * (matrix[1][0] * x[0] + matrix[1][1] * x[1])
}

fn chain_states(basis: &FrameVibrationBasis, q: [f64; 2], v: [f64; 2]) -> Vec<ModalAcousticState> {
    let (m, _, _) = chain_matrices();
    basis
        .modes()
        .iter()
        .map(|mode| {
            let project = |x: [f64; 2]| {
                mode.displacement[3] * (m[0][0] * x[0] + m[0][1] * x[1])
                    + mode.displacement[6] * (m[1][0] * x[0] + m[1][1] * x[1])
            };
            ModalAcousticState {
                displacement_m_sqrt_kg: project(q),
                velocity_m_sqrt_kg_per_s: project(v),
            }
        })
        .collect()
}

// Independent physical-coordinate ODE and its constitutive dissipation integral.
fn chain_derivative(y: [f64; 5]) -> [f64; 5] {
    let (m, k, c) = chain_matrices();
    let f: [f64; 2] =
        std::array::from_fn(|i| -k[i][0] * y[0] - k[i][1] * y[1] - c[i][0] * y[2] - c[i][1] * y[3]);
    let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    [
        y[2],
        y[3],
        (m[1][1] * f[0] - m[0][1] * f[1]) / det,
        (m[0][0] * f[1] - m[1][0] * f[0]) / det,
        quadratic(&c, [y[2], y[3]]),
    ]
}

fn chain_reference(q: [f64; 2], v: [f64; 2]) -> Vec<[f64; 5]> {
    let mut y = [q[0], q[1], v[0], v[1], 0.0];
    let h = 1.0 / (f64::from(RATE) * 128.0);
    let mut result = Vec::new();
    for _ in 0..64 {
        for _ in 0..128 {
            let k1 = chain_derivative(y);
            let k2 = chain_derivative(std::array::from_fn(|i| y[i] + 0.5 * h * k1[i]));
            let k3 = chain_derivative(std::array::from_fn(|i| y[i] + 0.5 * h * k2[i]));
            let k4 = chain_derivative(std::array::from_fn(|i| y[i] + h * k3[i]));
            y = std::array::from_fn(|i| {
                y[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
            });
        }
        result.push(y);
    }
    result
}

#[test]
fn heterogeneous_full_damping_reaches_motion_and_loss_without_discarding_modal_cross_terms() {
    let (q0, v0) = ([0.002, -0.001], [0.01, -0.015]);
    let reference = chain_reference(q0, v0);
    let (m, k, c) = chain_matrices();
    let initial = 0.5 * quadratic(&k, q0) + 0.5 * quadratic(&m, v0);
    for preparation in [None, Some(prepared())] {
        let basis = basis(&[1.0, 5.0], 2, BAND);
        let mut component = basis
            .coupled_component(
                RATE,
                &[0.0; 2],
                ModalAcousticTimeBudget::audible_reference(),
                72,
                || false,
            )
            .unwrap();
        component
            .restore_states(&chain_states(&basis, q0, v0))
            .unwrap();
        let projected_cross: f64 = component
            .material_connections()
            .iter()
            .map(|port| port.damping_n_s_m * port.left.shapes[0] * port.left.shapes[1])
            .sum();
        let (a, b) = (
            &basis.modes()[0].displacement,
            &basis.modes()[1].displacement,
        );
        let expected_cross =
            a[3] * (c[0][0] * b[3] + c[0][1] * b[6]) + a[6] * (c[1][0] * b[3] + c[1][1] * b[6]);
        assert!(expected_cross.abs() > 0.1);
        close(projected_cross, expected_cross, expected_cross.abs());
        let (basis, mut system) = component
            .into_system(config(), preparation, &CancelGate::new_clock_free())
            .unwrap();
        let mut loss = 0.0;
        for expected in &reference {
            let report = system.step(&[0.0; 2]).unwrap().clone();
            let (q, v) = physical(&basis, &system, 0);
            close(q[3], expected[0], 0.01);
            close(q[6], expected[1], 0.01);
            close(v[3], expected[2], 1.0);
            close(v[6], expected[3], 1.0);
            loss += report.connection_dissipation_j;
            close(loss, expected[4], initial);
            assert_eq!(report.component_dissipation_j, 0.0);
            assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);
        }
        close(system.total_energy_j().unwrap() + loss, initial, initial);
    }
}

#[test]
fn larger_reciprocal_system_retains_remapped_material_ports_initial_state_and_physical_work() {
    let gate = CancelGate::new_clock_free();
    let (mass, k, c) = (DENSITY * area() / 3.0, YOUNG * area(), 8.0 * area());
    let spring = 0.3 * k;
    let force = 0.02;
    let (ql, qr, vl, vr) = (0.004, -0.001, 0.01, -0.015);
    let left = basis(&[8.0], 1, BAND);
    let right = basis(&[8.0], 1, BAND);
    let initial_left = scalar_state(&left, ql, vl);
    let initial_right = scalar_state(&right, qr, vr);
    let mut lc = left
        .coupled_component(
            RATE,
            &[0.0],
            ModalAcousticTimeBudget::audible_reference(),
            18,
            || false,
        )
        .unwrap();
    let mut rc = right
        .coupled_component(
            RATE,
            &[0.0],
            ModalAcousticTimeBudget::audible_reference(),
            18,
            || false,
        )
        .unwrap();
    lc.restore_states(&initial_left).unwrap();
    rc.restore_states(&initial_right).unwrap();
    let (left, lmodel, mut links) = lc.into_parts(0);
    let (right, rmodel, rports) = rc.into_parts(1);
    assert_eq!(lmodel.states(), &initial_left);
    assert_eq!(rmodel.states(), &initial_right);
    assert!(
        rports
            .iter()
            .all(|p| p.left.component == 1 && p.right.component == 1)
    );
    links.extend(rports);
    links.push(ModalConnection {
        left: left.translation_attachment(0, 1, [1.0, 0.0]).unwrap().port,
        right: right.translation_attachment(1, 1, [1.0, 0.0]).unwrap().port,
        stiffness_n_m: spring,
        damping_n_s_m: 0.0,
        rest_extension_m: 0.0,
    });
    let mut system =
        CoupledModalSystem::new_prepared(vec![lmodel, rmodel], links, config(), prepared(), &gate)
            .unwrap();
    let initial = 0.5 * mass * (vl * vl + vr * vr)
        + 0.5 * k * (ql * ql + qr * qr)
        + 0.5 * spring * (ql - qr) * (ql - qr);
    close(system.total_energy_j().unwrap(), initial, initial);
    let mut generalized = [0.0];
    left.project_force_into(
        &[0.0, 0.0, 0.0, force, 0.0, 0.0],
        [0.0; 2],
        &mut generalized,
        || false,
    )
    .unwrap();
    let (mut work, mut loss) = (0.0, 0.0);
    for sample in 1..=64 {
        let report = system.step(&[generalized[0], 0.0]).unwrap().clone();
        let t = f64::from(sample) / f64::from(RATE);
        let (qs, vs) = free_oscillator(
            0.5 * (ql + qr) - force / (2.0 * k),
            0.5 * (vl + vr),
            k / mass,
            c / (2.0 * mass),
            t,
        );
        let (qa, va) = free_oscillator(
            0.5 * (ql - qr) - force / (2.0 * (k + 2.0 * spring)),
            0.5 * (vl - vr),
            (k + 2.0 * spring) / mass,
            c / (2.0 * mass),
            t,
        );
        let (qs, qa) = (
            qs + force / (2.0 * k),
            qa + force / (2.0 * (k + 2.0 * spring)),
        );
        let (ul, velocity_l) = physical(&left, &system, 0);
        let (ur, velocity_r) = physical(&right, &system, 1);
        close(ul[3], qs + qa, 0.01);
        close(ur[3], qs - qa, 0.01);
        close(velocity_l[3], vs + va, 0.5);
        close(velocity_r[3], vs - va, 0.5);
        work += report.external_work_j;
        loss += report.connection_dissipation_j;
        let expected_work = force * (qs + qa - ql);
        let energy = mass * (vs * vs + va * va) + k * qs * qs + (k + 2.0 * spring) * qa * qa;
        close(work, expected_work, initial);
        close(loss, initial + expected_work - energy, initial);
        assert_eq!(report.connection_forces_n.len(), 7);
        assert_eq!(report.component_dissipation_j, 0.0);
    }
    close(
        system.total_energy_j().unwrap() + loss,
        initial + work,
        initial,
    );
}

#[test]
fn source_full_spectrum_projection_and_cancellation_admission_preserve_retry() {
    let b = basis(&[8.0], 1, BAND);
    let modal_budget = ModalAcousticTimeBudget::audible_reference();
    assert!(matches!(
        b.time_model(RATE, &[0.0], modal_budget),
        Err(FrameVibrationError::Invalid(_))
    ));
    assert!(matches!(
        b.coupled_component(RATE, &[0.0], modal_budget, 17, || false),
        Err(FrameVibrationError::Invalid(_))
    ));
    for band in [(1.0, 4096.0), (0.0, 400.0)] {
        let bounded = basis(&[8.0], 1, band);
        assert!(matches!(
            bounded.coupled_component(RATE, &[0.0], modal_budget, 18, || false),
            Err(FrameVibrationError::Solid(_))
        ));
    }
    // Native chain frequencies are approximately 32 and 113 rad/s. Retaining
    // the first does not authorize a clock with Nyquist 24*pi = 75 rad/s.
    let truncated = basis(&[1.0, 5.0], 1, BAND);
    assert_eq!(truncated.available_mode_count(), 2);
    assert!(matches!(
        truncated.coupled_component(24, &[0.0], modal_budget, 36, || false),
        Err(FrameVibrationError::Solid(_))
    ));
    let mut polls = 0;
    let mut clean = b
        .coupled_component(RATE, &[0.0], modal_budget, 18, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        b.coupled_component(RATE, &[0.0], modal_budget, 18, || {
            seen += 1;
            seen == polls
        }),
        Err(FrameVibrationError::Cancelled)
    ));
    let mut retried = b
        .coupled_component(RATE, &[0.0], modal_budget, 18, || false)
        .unwrap();
    assert_eq!(clean.material_connections(), retried.material_connections());
    let state = scalar_state(&b, 0.003, 0.02);
    clean.restore_states(&state).unwrap();
    retried.restore_states(&state).unwrap();
    let gate = CancelGate::new_clock_free();
    let (_, mut expected) = clean
        .into_system(config(), Some(prepared()), &gate)
        .unwrap();
    let (_, mut actual) = retried
        .into_system(config(), Some(prepared()), &gate)
        .unwrap();
    assert_eq!(
        expected.step(&[0.01]).unwrap(),
        actual.step(&[0.01]).unwrap()
    );
    let accepted = actual.last_frame().unwrap().clone();
    let accepted_state = actual.components()[0].states().to_vec();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        actual.step_under_gate(&[0.01], &cancelled),
        Err(ModalCouplingError::Cancelled)
    ));
    assert_eq!(actual.last_frame(), Some(&accepted));
    assert_eq!(actual.components()[0].states(), accepted_state.as_slice());
    assert_eq!(
        expected.step(&[0.01]).unwrap(),
        actual.step_under_gate(&[0.01], &gate).unwrap()
    );
    assert_eq!(
        expected.components()[0].states(),
        actual.components()[0].states()
    );
}
