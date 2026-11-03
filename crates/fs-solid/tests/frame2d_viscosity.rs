//! Source viscosity reaches the real native frame history, not a modal fit.
//! Data are manufactured; the references are the continuous mechanical ODEs.

use fs_evidence::ValidityDomain;
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
use fs_solid::frame2d::{FrameMassFormulation, FrameViscousSection};
use fs_solid::{DynamicHistory, Frame2d, Geometry, Rayleigh, SectionModel};

const RHO: f64 = 1000.0;
const YOUNG: f64 = 2e9;
const RADIUS: f64 = 0.02;
const BAND: (f64, f64) = (0.0, 1e7);

fn source(eta_axial: f64, eta_bending: f64, band: (f64, f64)) -> FrameViscousSection {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [
        ("density", RHO, Density::DIMS),
        ("young_modulus", YOUNG, Pressure::DIMS),
        ("axial_eta", eta_axial, DynViscosity::DIMS),
        ("bending_eta", eta_bending, DynViscosity::DIMS),
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
                    source: "manufactured frame viscosity".into(),
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
            process: "linear viscosity fixture".into(),
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

fn guided(section: FrameViscousSection, axis: usize) -> Frame2d {
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0);
    frame.add_node(1.0, 0.0);
    frame.fix(0, [true; 3]);
    frame.fix(
        1,
        if axis == 0 {
            [false, true, true]
        } else {
            [true, false, true]
        },
    );
    frame
        .add_material_viscous_element(0, 1, 5, section)
        .unwrap();
    frame
}

fn near(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={absolute:.3e}"
    );
}

fn closed_energy(history: &DynamicHistory, initial_energy: f64) {
    // The native relative-input metric has a zero denominator in free decay.
    // Check the same actual work identity against the nonzero initial storage.
    for i in 0..history.t.len() {
        let residual = history.input_energy[i]
            - (history.kinetic_energy[i] - history.kinetic_energy[0])
            - history.damping_energy[i]
            - history.internal_work[i];
        near(residual, 0.0, 1e-7 * initial_energy);
    }
}

#[test]
fn source_axial_and_bending_swaps_reproduce_continuous_damped_motion_and_loss() {
    let area = core::f64::consts::PI * RADIUS.powi(2);
    let inertia = area * RADIUS.powi(2) / 4.0;
    for axis in [0, 1] {
        let (mass, stiffness) = if axis == 0 {
            (RHO * area / 3.0, YOUNG * area)
        } else {
            // Unit-length guided bending includes the section's slope-spin
            // kinetic energy as well as centerline translation.
            (
                RHO * (area * 13.0 / 35.0 + 6.0 * inertia / 5.0),
                12.0 * YOUNG * inertia,
            )
        };
        let omega = (stiffness / mass).sqrt();
        let dt = 0.01 / omega;
        let rayleigh = Rayleigh {
            a0: 0.05 * omega,
            a1: 0.02 / omega,
        };
        let cr = rayleigh.a0 * mass + rayleigh.a1 * stiffness;
        let initial = 1e-5;
        let initial_energy = 0.5 * stiffness * initial * initial;
        let mut final_energy = Vec::new();
        for scale in [0.0, 1.0, 4.0] {
            let mut frame = guided(source(scale * 2e5, scale * 4e6, BAND), axis);
            let bound = frame.material_viscous_section(0).unwrap();
            near(
                bound.axial().viscous_rigidity_n_s(),
                scale * 2e5 * area,
                1e-10,
            );
            near(
                bound.bending().viscous_stiffness_n_m2_s(),
                scale * 4e6 * inertia,
                1e-12,
            );
            let material_c = if axis == 0 {
                scale * 2e5 * area
            } else {
                12.0 * scale * 4e6 * inertia
            };
            near(frame.natural_frequencies().unwrap()[0], omega, 1e-9 * omega);
            let mut force = vec![0.0; frame.ndof()];
            force[3 + axis] = stiffness * initial;
            frame.static_load(&force, 1).unwrap();
            near(frame.displacements()[3 + axis], initial, 1e-12 * initial);
            let history = frame
                .newmark(&[0.0; 801], dt, (0.0, 0.0), rayleigh)
                .unwrap();
            let decay = (material_c + cr) / (2.0 * mass);
            let wd = (omega * omega - decay * decay).sqrt();
            for (time, displacement) in history.t.iter().zip(&history.u) {
                let expected = initial
                    * (-decay * time).exp()
                    * ((wd * time).cos() + decay / wd * (wd * time).sin());
                near(displacement[3 + axis], expected, 8e-5 * initial);
            }
            let time = *history.t.last().unwrap();
            let q = initial
                * (-decay * time).exp()
                * ((wd * time).cos() + decay / wd * (wd * time).sin());
            let v = -initial * (-decay * time).exp() * omega * omega / wd * (wd * time).sin();
            let remaining = 0.5 * stiffness * q * q + 0.5 * mass * v * v;
            let lost = initial_energy - remaining;
            let material = *history.material_damping_energy.last().unwrap();
            let authored = *history.rayleigh_damping_energy.last().unwrap();
            near(
                material,
                lost * material_c / (material_c + cr),
                1e-4 * initial_energy,
            );
            near(
                authored,
                lost * cr / (material_c + cr),
                1e-4 * initial_energy,
            );
            near(
                material + authored,
                *history.damping_energy.last().unwrap(),
                1e-12 * initial_energy,
            );
            closed_energy(&history, initial_energy);
            assert!(material >= 0.0 && authored > 0.0);
            if scale == 0.0 {
                assert_eq!(material, 0.0);
            }
            final_energy.push(remaining);
        }
        assert!(final_energy[2] < final_energy[1] && final_energy[1] < final_energy[0]);
    }
}

#[test]
fn physical_operator_keeps_rotational_coupling_without_damping_rigid_motion() {
    let frame = guided(source(2e5, 4e6, BAND), 1);
    let operator = frame.material_damping_operator().unwrap().unwrap();
    let matrix = operator.entries();
    let cb = 4e6 * core::f64::consts::PI * RADIUS.powi(4) / 4.0;
    // Hermite curvature-rate integrals on L=1: 12, -6 and 4.
    near(matrix[4 * 6 + 4], 12.0 * cb, 1e-12);
    near(matrix[4 * 6 + 5], -6.0 * cb, 1e-12);
    near(matrix[5 * 6 + 5], 4.0 * cb, 1e-12);
    for velocity in [
        [1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 1.0, 1.0],
    ] {
        for row in 0..6 {
            near(
                (0..6)
                    .map(|col| matrix[row * 6 + col] * velocity[col])
                    .sum(),
                0.0,
                1e-12,
            );
        }
    }

    // Non-unit length and oblique orientation expose the metre-valued
    // rotational weights that a reciprocal scalar force port must retain.
    let mut oblique =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    oblique.add_node(0.0, 0.0);
    oblique.add_node(1.2, 1.6);
    oblique.fix(0, [true; 3]);
    oblique
        .add_material_viscous_element(0, 1, 5, source(2e5, 4e6, BAND))
        .unwrap();
    let operator = oblique.material_damping_operator().unwrap().unwrap();
    let ports = operator.ports();
    assert_eq!(ports.len(), 3);
    let ca = 2e5 * core::f64::consts::PI * RADIUS.powi(2);
    for (port, coefficient, weights) in [
        (&ports[0], ca / 2.0, [-0.6, -0.8, 0.0, 0.6, 0.8, 0.0]),
        (&ports[1], 3.0 * cb / 8.0, [-1.6, 1.2, 2.0, 1.6, -1.2, 2.0]),
        (&ports[2], cb / 8.0, [0.0, 0.0, 2.0, 0.0, 0.0, -2.0]),
    ] {
        assert_eq!(port.element, 0);
        assert_eq!(port.dofs, [0, 1, 2, 3, 4, 5]);
        near(port.damping_n_s_m, coefficient, 1e-12);
        for (actual, expected) in port.displacement_weights.iter().zip(weights) {
            near(*actual, expected, 1e-14);
        }
    }
    let velocity = [0.3, -0.7, 0.2, -0.1, 0.5, -0.9];
    let port_power: f64 = ports
        .iter()
        .map(|port| {
            let rate: f64 = (0..6)
                .map(|i| port.displacement_weights[i] * velocity[port.dofs[i]])
                .sum();
            port.damping_n_s_m * rate * rate
        })
        .sum();
    let matrix_power: f64 = (0..6)
        .map(|i| {
            (0..6)
                .map(|j| velocity[i] * operator.entries()[i * 6 + j] * velocity[j])
                .sum::<f64>()
        })
        .sum();
    near(port_power, matrix_power, 1e-12);
}

fn axial_chain() -> Frame2d {
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    for node in 0..3 {
        frame.add_node(0.5 * node as f64, 0.0);
        frame.fix(
            node,
            if node == 0 {
                [true; 3]
            } else {
                [false, true, true]
            },
        );
    }
    frame
        .add_material_viscous_element(0, 1, 3, source(2e4, 0.0, BAND))
        .unwrap();
    frame
        .add_material_viscous_element(1, 2, 7, source(1.6e5, 0.0, BAND))
        .unwrap();
    let mut force = vec![0.0; 9];
    force[6] = 1.0;
    frame.static_load(&force, 1).unwrap();
    frame
}

// Independent continuous equation: M*q'' + C*q' + K*q = 0, and Q'=q'^T C q'.
// No frame matrix, native mode, Newmark recurrence or reported loss is used.
fn reference_derivative(y: [f64; 5]) -> [f64; 5] {
    let area = core::f64::consts::PI * RADIUS.powi(2);
    let m = RHO * area * 0.5;
    let (m00, m01, m11) = (2.0 * m / 3.0, m / 6.0, m / 3.0);
    let k = YOUNG * area / 0.5;
    let (c1, c2) = (2e4 * area / 0.5, 1.6e5 * area / 0.5);
    let viscous = [(c1 + c2) * y[2] - c2 * y[3], c2 * (y[3] - y[2])];
    let rhs = [
        -k * (2.0 * y[0] - y[1]) - viscous[0],
        -k * (y[1] - y[0]) - viscous[1],
    ];
    let det = m00 * m11 - m01 * m01;
    [
        y[2],
        y[3],
        (m11 * rhs[0] - m01 * rhs[1]) / det,
        (m00 * rhs[1] - m01 * rhs[0]) / det,
        y[2] * viscous[0] + y[3] * viscous[1],
    ]
}

fn continuous_reference(end: f64, steps: usize) -> Vec<[f64; 5]> {
    let compliance = 0.5 / (YOUNG * core::f64::consts::PI * RADIUS.powi(2));
    let mut states = vec![[compliance, 2.0 * compliance, 0.0, 0.0, 0.0]];
    let h = end / steps as f64;
    for _ in 0..steps {
        let y = *states.last().unwrap();
        let k1 = reference_derivative(y);
        let k2 = reference_derivative(std::array::from_fn(|i| y[i] + 0.5 * h * k1[i]));
        let k3 = reference_derivative(std::array::from_fn(|i| y[i] + 0.5 * h * k2[i]));
        let k4 = reference_derivative(std::array::from_fn(|i| y[i] + h * k3[i]));
        states.push(std::array::from_fn(|i| {
            y[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
        }));
    }
    states
}

#[test]
fn heterogeneous_nonproportional_loss_retains_cross_terms_and_refines_to_continuous_ode() {
    let duration = 0.002;
    let reference = continuous_reference(duration, 8192);
    let amplitude = reference[0][1];
    let operator = axial_chain().material_damping_operator().unwrap().unwrap();
    assert_eq!(operator.ndof(), 9);
    assert_eq!(operator.free_dofs(), &[3, 6]);
    assert_eq!(operator.sections().len(), 2);
    assert_eq!(operator.omega_band_rad_s(), BAND);
    operator.admit_timestep(duration / 512.0).unwrap();
    let area = core::f64::consts::PI * RADIUS.powi(2);
    let (c1, c2) = (2e4 * area / 0.5, 1.6e5 * area / 0.5);
    let entries = operator.entries();
    near(entries[3 * 9 + 3], c1 + c2, 1e-11);
    near(entries[3 * 9 + 6], -c2, 1e-11);
    near(entries[6 * 9 + 6], c2, 1e-11);
    assert_eq!(entries[3 * 9 + 6], entries[6 * 9 + 3]);
    let mut errors = Vec::new();
    for steps in [256, 512] {
        let mut frame = axial_chain();
        let history = frame
            .newmark(
                &vec![0.0; steps + 1],
                duration / steps as f64,
                (0.0, 0.0),
                Rayleigh { a0: 0.0, a1: 0.0 },
            )
            .unwrap();
        let mut error = 0.0_f64;
        for (sample, u) in history.u.iter().enumerate() {
            let expected = reference[sample * (8192 / steps)];
            error = error
                .max((u[3] - expected[0]).abs())
                .max((u[6] - expected[1]).abs());
        }
        errors.push(error / amplitude);
        near(
            *history.material_damping_energy.last().unwrap(),
            reference.last().unwrap()[4],
            0.005 * amplitude,
        );
        assert_eq!(*history.rayleigh_damping_energy.last().unwrap(), 0.0);
        closed_energy(&history, 0.5 * amplitude); // Initial tip load is exactly 1 N.
    }
    assert!(errors[0] < 0.01, "coarse relative error {}", errors[0]);
    assert!(
        errors[1] < 0.3 * errors[0],
        "second-order refinement: {errors:?}"
    );
}

#[test]
fn continuation_and_rejected_band_update_preserve_actual_velocity_and_source_history() {
    let mut full = axial_chain();
    let mut split = full.clone();
    let dt = 5e-6;
    let damping = Rayleigh { a0: 0.7, a1: 1e-7 };
    let whole = full.newmark(&[0.0; 401], dt, (0.0, 0.0), damping).unwrap();
    let first = split.newmark(&[0.0; 151], dt, (0.0, 0.0), damping).unwrap();
    let mut retry = split.clone();
    // This reaches the source admission after a nonzero velocity was accepted.
    assert!(split.newmark(&[0.0; 2], 1e-9, (0.0, 0.0), damping).is_err());
    assert_eq!(split.displacements(), retry.displacements());
    let second = split.newmark(&[0.0; 251], dt, (0.0, 0.0), damping).unwrap();
    let repeated = retry.newmark(&[0.0; 251], dt, (0.0, 0.0), damping).unwrap();
    assert_eq!(second.u, repeated.u);
    assert_eq!(second.kinetic_energy, repeated.kinetic_energy);
    assert_eq!(
        second.material_damping_energy,
        repeated.material_damping_energy
    );
    for (&a, &b) in full.displacements().iter().zip(split.displacements()) {
        near(a, b, 1e-13);
    }
    near(
        first.material_damping_energy.last().unwrap()
            + second.material_damping_energy.last().unwrap(),
        *whole.material_damping_energy.last().unwrap(),
        1e-12,
    );
    near(
        first.rayleigh_damping_energy.last().unwrap()
            + second.rayleigh_damping_energy.last().unwrap(),
        *whole.rayleigh_damping_energy.last().unwrap(),
        1e-12,
    );
}

#[test]
fn source_loss_never_bypasses_linear_inertial_or_frequency_admission() {
    let section = source(2e5, 4e6, BAND);
    let mut rotating = Frame2d::new(Geometry::Corotational);
    rotating.add_node(0.0, 0.0);
    rotating.add_node(1.0, 0.0);
    assert!(
        rotating
            .add_material_viscous_element(0, 1, 5, section.clone())
            .is_err()
    );
    assert!(rotating.material_section(0).is_none());
    let mut lumped = Frame2d::new(Geometry::Linear);
    lumped.add_node(0.0, 0.0);
    lumped.add_node(1.0, 0.0);
    lumped.fix(0, [true; 3]);
    lumped
        .add_material_viscous_element(0, 1, 5, section)
        .unwrap();
    let none = Rayleigh { a0: 0.0, a1: 0.0 };
    // Source-derived section spin gives lumped material rotations real inertia.
    assert!(lumped.newmark(&[0.0; 2], 1e-5, (0.0, 0.0), none).is_ok());
    // A genuinely massless raw member still cannot bypass viscous admission.
    let tip = lumped.add_node(2.0, 0.0);
    lumped
        .add_element(1, tip, 5, 0.0, &|| SectionModel::Elastic {
            ea: 1e6,
            ei: 1e3,
        })
        .unwrap();
    assert!(lumped.newmark(&[0.0; 2], 1e-5, (0.0, 0.0), none).is_err());
    for band in [(1.0, 1e7), (0.0, 100.0)] {
        let mut frame = guided(source(2e5, 4e6, band), 1);
        assert!(frame.newmark(&[0.0; 2], 1e-5, (0.0, 0.0), none).is_err());
    }
    let mut unresolved = guided(source(2e5, 4e6, BAND), 0);
    assert!(
        unresolved
            .newmark(&[0.0; 2], 0.01, (0.0, 0.0), none)
            .is_err()
    );
    let mut mixed = guided(source(2e5, 4e6, BAND), 1);
    mixed
        .add_element(0, 1, 5, 0.0, &|| SectionModel::Elastic {
            ea: -1.0,
            ei: 1.0,
        })
        .unwrap();
    assert!(mixed.newmark(&[0.0; 2], 1e-5, (0.0, 0.0), none).is_err());
}

fn with_harmonic_cx(f: impl FnOnce(&fs_exec::Cx<'_>)) {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 5,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

/// G1: the source's full heterogeneous loss matrix acts on the native forced
/// frame, including mutual inertia, phase and independently computed heat.
#[test]
fn harmonic_response_retains_nonproportional_source_loss_and_consistent_mass() {
    use fs_math::c64::C64;
    use fs_modal::harmonic::HarmonicResponseControls;
    with_harmonic_cx(|cx| {
        let frame = axial_chain();
        let accepted = frame.displacements().to_vec();
        let controls = HarmonicResponseControls {
            maximum_coordinates: 9,
            maximum_payload_bytes: 1 << 20,
            maximum_relative_residual: 1e-12,
        };
        let area = core::f64::consts::PI * RADIUS.powi(2);
        let member_mass = RHO * area * 0.5;
        let stiffness = YOUNG * area / 0.5;
        let (c1, c2) = (2e4 * area / 0.5, 1.6e5 * area / 0.5);
        let force = [C64::new(0.7, -0.2), C64::new(-0.3, 0.5)];
        for omega in [0.0, 1000.0, 2300.0, 5000.0] {
            let prepared = frame
                .prepare_harmonic_response(omega, Rayleigh { a0: 0.0, a1: 0.0 }, controls, cx)
                .unwrap();
            assert_eq!(prepared.free_dofs(), &[3, 6]);
            let entry = |k: f64, m: f64, c: f64| C64::new(k - omega * omega * m, -omega * c);
            let a = entry(2.0 * stiffness, 2.0 * member_mass / 3.0, c1 + c2);
            let b = entry(-stiffness, member_mass / 6.0, -c2);
            let d = entry(stiffness, member_mass / 3.0, c2);
            let determinant = a * d - b * b;
            let expected = [
                (d * force[0] - b * force[1]) / determinant,
                (a * force[1] - b * force[0]) / determinant,
            ];
            let response = prepared.solve(&force, cx).unwrap();
            for (actual, reference) in response.displacement.iter().zip(expected) {
                assert!((*actual - reference).abs() < 3e-12 * reference.abs());
            }
            let v = expected.map(|q| C64::new(q.im, -q.re).scale(omega));
            // Original axial strain-rate ports: eta*A/L * |v_j-v_i|^2 / 2.
            let heat = 0.5 * (c1 * v[0].norm_sq() + c2 * (v[1] - v[0]).norm_sq());
            near(
                response.average_viscous_power_w,
                heat,
                5e-12 * heat.max(1e-30),
            );
            near(
                response.average_input_power_w,
                heat,
                5e-12 * heat.max(1e-30),
            );
        }
        assert_eq!(frame.displacements(), accepted);
    });
}

/// G0/G1: a harmonic query admits its actual source frequency, not a guessed
/// time grid; the original source band and viscosity remain authoritative.
#[test]
fn harmonic_viscosity_checks_original_frequency_band_and_counts_rayleigh_once() {
    use fs_math::c64::C64;
    use fs_modal::harmonic::HarmonicResponseControls;
    with_harmonic_cx(|cx| {
        let controls = HarmonicResponseControls {
            maximum_coordinates: 6,
            maximum_payload_bytes: 1 << 20,
            maximum_relative_residual: 1e-12,
        };
        let frame = guided(source(2e5, 4e6, (100.0, 5000.0)), 0);
        let rayleigh = Rayleigh { a0: 0.7, a1: 1e-7 };
        for omega in [0.0, 99.0, 5001.0] {
            assert!(
                frame
                    .prepare_harmonic_response(omega, rayleigh, controls, cx)
                    .is_err()
            );
        }
        let area = core::f64::consts::PI * RADIUS.powi(2);
        let mass = RHO * area / 3.0;
        let stiffness = YOUNG * area;
        let damping = 2e5 * area + rayleigh.a0 * mass + rayleigh.a1 * stiffness;
        for omega in [100.0, 2500.0, 5000.0] {
            let response = frame
                .prepare_harmonic_response(omega, rayleigh, controls, cx)
                .unwrap()
                .solve(&[C64::ONE], cx)
                .unwrap();
            let q = C64::ONE / C64::new(stiffness - omega * omega * mass, -omega * damping);
            assert!((response.displacement[0] - q).abs() < 3e-12 * q.abs());
            let heat = 0.5 * damping * omega * omega * q.norm_sq();
            near(response.average_viscous_power_w, heat, 5e-12 * heat);
            near(response.average_input_power_w, heat, 5e-12 * heat);
        }
        assert_eq!(
            frame
                .material_viscous_section(0)
                .unwrap()
                .omega_band_rad_s(),
            (100.0, 5000.0)
        );
    });
}

/// G1: free translation and the elastic resonance coexist in the same source
/// frame. A positive harmonic frequency does not require removal of rigid modes.
#[test]
fn harmonic_source_frame_keeps_rigid_motion_and_its_elastic_resonance() {
    use fs_math::c64::C64;
    use fs_modal::harmonic::HarmonicResponseControls;
    with_harmonic_cx(|cx| {
        let mut frame = Frame2d::with_mass_formulation(
            Geometry::Linear,
            FrameMassFormulation::ConsistentReference,
        );
        frame.add_node(0.0, 0.0);
        frame.add_node(1.0, 0.0);
        frame.fix(0, [false, true, true]);
        frame.fix(1, [false, true, true]);
        frame
            .add_material_viscous_element(0, 1, 5, source(2e5, 4e6, BAND))
            .unwrap();
        // Preserve the existing public transient owner's unresolved-mode rule.
        assert!(frame.material_damping_operator().is_err());
        let controls = HarmonicResponseControls {
            maximum_coordinates: 6,
            maximum_payload_bytes: 1 << 20,
            maximum_relative_residual: 1e-12,
        };
        let area = core::f64::consts::PI * RADIUS.powi(2);
        let mass = RHO * area;
        let stiffness = YOUNG * area;
        let damping = 2e5 * area;
        for omega in [1000.0, (12.0 * stiffness / mass).sqrt()] {
            let prepared = frame
                .prepare_harmonic_response(omega, Rayleigh { a0: 0.0, a1: 0.0 }, controls, cx)
                .unwrap();
            assert_eq!(prepared.free_dofs(), &[0, 3]);
            let rigid = prepared.solve(&[C64::ONE, C64::ONE], cx).unwrap();
            let translation = -2.0 / (mass * omega * omega);
            for q in &rigid.displacement {
                assert!((*q - C64::from_re(translation)).abs() < 2e-12 * translation.abs());
            }
            assert!(rigid.average_viscous_power_w.abs() < 1e-22);
            let differential = prepared.solve(&[C64::ONE, C64::from_re(-1.0)], cx).unwrap();
            let extension = C64::ONE
                / C64::new(
                    2.0 * stiffness - omega * omega * mass / 6.0,
                    -2.0 * omega * damping,
                );
            assert!((differential.displacement[0] - extension).abs() < 3e-12 * extension.abs());
            assert!((differential.displacement[1] + extension).abs() < 3e-12 * extension.abs());
            let heat = 2.0 * damping * omega * omega * extension.norm_sq();
            near(differential.average_viscous_power_w, heat, 5e-12 * heat);
            near(differential.average_input_power_w, heat, 5e-12 * heat);
        }
    });
}

#[test]
fn impulse_right_limit_uses_source_viscosity_and_actual_continuation_clock() {
    let area = core::f64::consts::PI * RADIUS.powi(2);
    let mass = RHO * area / 3.0;
    let stiffness = YOUNG * area;
    let eta = 2e5;
    let material_c = eta * area;
    let damping = Rayleigh { a0: 0.25, a1: 1e-6 };
    let c = material_c + damping.a0 * mass + damping.a1 * stiffness;
    let mut original = guided(source(eta, 4e6, BAND), 0);
    let mut preload = vec![0.0; original.ndof()];
    preload[3] = stiffness * 2e-6;
    original.static_load(&preload, 1).unwrap();
    let mut impulse = vec![0.0; original.ndof()];
    impulse[3] = mass * 0.125;
    let mut post_load = preload.clone();
    post_load[3] *= 0.5;
    let expected = (post_load[3] - stiffness * original.displacements()[3] - c * 0.125) / mass;
    let mut first = original.clone();
    let mut second = original.clone();
    let a = first
        .apply_impulse(&impulse, &post_load, damping, 1e-6)
        .unwrap();
    let b = second
        .apply_impulse(&impulse, &post_load, damping, 2e-6)
        .unwrap();
    assert_eq!(
        a, b,
        "clock admission must not turn impulse into a duration-dependent force"
    );
    assert_eq!(first.velocities(), second.velocities());
    assert_eq!(first.accelerations(), second.accelerations());
    near(first.velocities()[3], 0.125, 1e-14);
    near(first.accelerations()[3], expected, 2e-10);
    near(a.impulse_work_j, 0.5 * mass * 0.125_f64.powi(2), 1e-14);
    assert_eq!(first.displacements(), original.displacements());
    assert_eq!(
        first.section_deformations(0),
        original.section_deformations(0)
    );
    let acceleration = first.accelerations().to_vec();
    let history = first
        .newmark_with_load(&[0.0], 1e-6, (0.0, 0.0), damping, &post_load)
        .unwrap();
    near(history.kinetic_energy[0], a.kinetic_energy_after_j, 1e-14);
    assert_eq!(first.accelerations(), acceleration);
    for dt in [1e-9, 1.0] {
        let mut refused = original.clone();
        let before = format!("{refused:?}");
        assert!(
            refused
                .apply_impulse(&impulse, &post_load, damping, dt)
                .is_err()
        );
        assert_eq!(
            format!("{refused:?}"),
            before,
            "a source-band or complete-spectrum refusal must not publish the proposed velocity"
        );
    }
}
