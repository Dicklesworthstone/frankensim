//! Independent continuous forced shell/attachment oracles and real source loss.
#![deny(unsafe_code)]

#[path = "coating_shell.rs"]
mod fixture;

use fs_couple::shear_shell::{
    CoatingSide, MaterialShellImage, ShellProjectionBudget, ShellReductionOptions,
    assembly::{
        MaterialShellAssembly, MaterialShellComponent, ShellAssemblyForcePort,
        ShellAssemblyReductionOptions,
    },
    coating::{
        poles::{
            ShellHarmonicPort, ShellResonanceBudget, ShellResonanceError, ShellResonanceLoad,
            response::{ShellForcedBudget, ShellForcedResponse},
        },
        relaxation::CoatingRelaxationOptions,
        strain::ShellCoatingStrainBudget,
    },
    constant_loss::relaxation::{SubstrateRelaxationOptions, SubstrateStorageReference},
    pressure::ShellPressureSurfaces,
};
use fs_exec::CancelGate;
use fs_material::coating::ResolvedCoating;
use fs_math::c64::C64;
use fs_phs::{
    PortHamiltonian, QuadraticStorage,
    harmonic::{HarmonicBudget, HarmonicError, HarmonicModel},
};
use fs_plate::shell::{
    junction::TieBudget,
    mitc::ShellMass,
    pressure::{OrientedShellFacet, ShellFace},
};

fn image(gate: &CancelGate) -> MaterialShellImage {
    image_with_modes(gate, 1)
}

fn image_with_modes(gate: &CancelGate, modes: usize) -> MaterialShellImage {
    let (shell, directors, constraints) = fixture::coupon(gate);
    let shell = shell
        .with_coatings(vec![fixture::layer(CoatingSide::Positive, 20e-6)], gate)
        .unwrap()
        .with_coating_relaxation(
            CoatingRelaxationOptions {
                terms: 16,
                storage_relative_tolerance: 0.001,
                loss_relative_tolerance: 0.001,
                max_memory_coordinates: 256,
            },
            gate,
        )
        .unwrap()
        .with_substrate_relaxation(
            SubstrateRelaxationOptions {
                storage_reference: SubstrateStorageReference::Relaxed,
                max_memory_coordinates: 256,
            },
            gate,
        )
        .unwrap();
    let facets = |side| {
        (0..shell.mesh().tris.len())
            .map(|triangle| OrientedShellFacet { triangle, side })
            .collect()
    };
    let surfaces = ShellPressureSurfaces {
        exterior: facets(ShellFace::Positive),
        cavity: facets(ShellFace::Negative),
    };
    let shell = shell.with_pressure_surfaces(surfaces, gate).unwrap();
    let ports: Vec<_> = [5, 6]
        .into_iter()
        .map(|node| {
            let mut row = vec![0.0; 6 * shell.mesh().nodes.len()];
            row[6 * node + 2] = 1.0;
            row
        })
        .collect();
    shell
        .reduce(
            &directors,
            &ShellReductionOptions {
                constraints,
                mass: ShellMass::Consistent,
                drilling_alpha: 1e-3,
                eigenvalue_window: (0.0, 1e9),
                n_modes: modes,
                assembly_budget: fixture::budget(&shell),
                projection_budget: ShellProjectionBudget {
                    max_modes: modes,
                    max_projection_terms: 100_000_000,
                    reference_resistance_n_s_m: 1.0,
                },
            },
            &ports,
            gate,
        )
        .unwrap()
}

fn budget() -> ShellResonanceBudget {
    ShellResonanceBudget {
        maximum_coordinates: 128,
        maximum_payload_bytes: 16_000_000,
        maximum_preparation_work: 100_000_000_000,
        physical_relative_tolerance: 1e-10,
    }
}

fn forced_budget() -> ShellForcedBudget {
    ShellForcedBudget {
        maximum_payload_bytes: 1_000_000,
        maximum_observation_work: 10_000_000,
    }
}

// Physical attachment state [spring displacement relative to support, momentum].
// Ground damping acts on momentum. Each signed map multiplies support velocity;
// the full mutual impedance is retained by the one original load state pair.
fn attachment(mass: f64, stiffness: f64, damping: f64, maps: &[f64]) -> PortHamiltonian {
    let mut g = vec![0.0; 2 * maps.len()];
    for (i, &a) in maps.iter().enumerate() {
        g[i] = -a;
    }
    PortHamiltonian::new(
        2,
        maps.len(),
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, damping],
        g,
        Box::new(QuadraticStorage::new(vec![stiffness, 0.0, 0.0, 1.0 / mass], 2).unwrap()),
    )
    .unwrap()
}

fn branches(image: &MaterialShellImage) -> Vec<(&[f64], f64)> {
    let substrate = image.substrate_relaxation().unwrap();
    let coating = image.coating_relaxation().unwrap();
    substrate
        .coupling_rows_per_s()
        .iter()
        .zip(substrate.relaxation_times_s())
        .chain(
            coating
                .coupling_rows_per_s()
                .iter()
                .zip(coating.relaxation_times_s()),
        )
        .map(|(row, &tau)| (row.as_slice(), tau))
        .collect()
}

fn body_operator(image: &MaterialShellImage, s: C64, w: f64) -> C64 {
    let mut value = s * s + s.scale(image.damping().matrix_per_s()[0]) + C64::from_re(w * w);
    for (row, tau) in branches(image) {
        value = value + s.scale(tau * row[0] * row[0]) / (C64::ONE + s.scale(tau));
    }
    value
}

// Independent two-body elimination with the load denominator multiplied out:
// this oracle is defined even at the exact standalone attachment resonance.
#[allow(clippy::too_many_arguments)]
fn analytic(
    image: &MaterialShellImage,
    omega: f64,
    w: f64,
    beta: f64,
    mass: f64,
    stiffness: f64,
    gamma: f64,
    force: C64,
) -> (C64, C64, C64) {
    let s = C64::new(0.0, -omega);
    let a = s * s + s.scale(gamma);
    let load = a + C64::from_re(stiffness / mass);
    let determinant = body_operator(image, s, w) * load + a.scale(beta * beta * stiffness);
    (
        force * load / determinant,
        force * a.scale(-beta) / determinant,
        force * s.scale(beta * stiffness) / determinant,
    )
}

fn near(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 3e-10 * expected.abs().max(1e-200),
        "actual {actual:?}, expected {expected:?}"
    );
}

fn near_real(actual: f64, expected: f64) {
    near(C64::from_re(actual), C64::from_re(expected));
}

fn source_channels(image: &MaterialShellImage, result: &ShellForcedResponse, q: C64, w: f64) {
    let omega = result.physical.omega_rad_s;
    let s = C64::new(0.0, -omega);
    let v = s * q;
    near(result.displacement_m_sqrt_kg[0], q);
    near(result.velocity_m_sqrt_kg_per_s[0], v);
    near_real(
        result.mean_mechanical_storage_j,
        0.25 * (v.norm_sq() + q.scale(w).norm_sq()),
    );
    near_real(
        result.background_dissipation_w,
        0.5 * image.damping().matrix_per_s()[0] * v.norm_sq(),
    );
    assert_eq!(result.interface_dissipation_w, 0.0);
    let all = branches(image);
    let mut total = result.mean_mechanical_storage_j;
    let mut heat = result.background_dissipation_w;
    for (i, ((row, tau), (&storage, &loss))) in all
        .iter()
        .zip(
            result
                .substrate_branch_mean_storage_j
                .iter()
                .chain(&result.coating_branch_mean_storage_j)
                .zip(
                    result
                        .substrate_branch_dissipation_w
                        .iter()
                        .chain(&result.coating_branch_dissipation_w),
                ),
        )
        .enumerate()
    {
        let z = q.scale(row[0]) / (C64::ONE + s.scale(*tau));
        let extension = q.scale(row[0]) - z;
        near(result.physical.state[2 + i], z);
        near_real(storage, 0.25 * extension.norm_sq());
        near_real(loss, 0.5 * extension.norm_sq() / tau);
        total += storage;
        heat += loss;
    }
    assert_eq!(
        result.physical.state.len(),
        2 + all.len() + result.load.as_ref().map_or(0, |l| l.state.len())
    );
    if let Some(load) = &result.load {
        total += load.mean_stored_energy_j;
        heat += load.dissipated_power_w;
    }
    near_real(result.physical.mean_stored_energy_j, total);
    near_real(result.physical.dissipated_power_w, heat);
    near_real(result.physical.supplied_power_w, heat);
    assert!(result.physical_relative_residual < 1e-10);
    assert!(result.original_equation_relative_residual < 1e-10);
    for (i, row) in image.force_ports().iter().enumerate() {
        near(result.force_port_displacement_m[i], q.scale(row[0]));
        near(result.force_port_velocity_m_s()[i], v.scale(row[0]));
    }
    near(
        result.exterior_volume_velocity_m3_s(),
        v.scale(image.volume_shapes_m2_per_sqrt_kg()[0]),
    );
    near(
        result.cavity_volume_velocity_m3_s(),
        v.scale(image.cavity_volume_shapes_m2_per_sqrt_kg()[0]),
    );
}

#[test]
fn g1_joint_forcing_matches_independent_causal_two_body_response_and_actual_heat() {
    let gate = CancelGate::new_clock_free();
    let image = image(&gate);
    let bare = image.prepare_coating_resonance(budget(), &gate).unwrap();
    let w = bare.equilibrium_angular_frequencies_rad_s().next().unwrap();
    let beta = image.force_ports()[0][0];
    let mass = 0.75 / (beta * beta);
    let stiffness = mass * 1_000_000.0;
    let gamma = 20.0;
    let prepared = image
        .prepare_coating_resonance_with_load(
            ShellResonanceLoad {
                system: attachment(mass, stiffness, mass * gamma, &[1.0]),
                ports: &[ShellHarmonicPort::Force(0)],
                band_rad_s: bare.source_band_rad_s(),
            },
            budget(),
            &gate,
        )
        .unwrap();
    assert_eq!(prepared.physical_port_count(), 4);
    let input = [C64::new(0.7, -0.2), C64::ZERO, C64::ZERO, C64::ZERO];
    for omega in [800.0, 1000.0, 1200.0, 5000.0] {
        let factor = prepared
            .harmonic_factor(omega, forced_budget(), &gate)
            .unwrap();
        assert_eq!(factor.omega_rad_s(), omega);
        let result = factor.response(&input, &gate).unwrap();
        let (q, r, p) = analytic(
            &image,
            omega,
            w,
            beta,
            mass,
            stiffness,
            gamma,
            input[0].scale(beta),
        );
        source_channels(&image, &result, q, w);
        let load = result.load.as_ref().unwrap();
        near(load.state[0], r);
        near(load.state[1], p);
        near(load.effort[0], r.scale(-stiffness));
        near(load.flow[0], C64::new(0.0, -omega) * q.scale(beta));
        near_real(
            load.mean_stored_energy_j,
            0.25 * (stiffness * r.norm_sq() + p.norm_sq() / mass),
        );
        near_real(load.dissipated_power_w, 0.5 * gamma * p.norm_sq() / mass);
        near_real(load.absorbed_port_power_w[0], load.dissipated_power_w);
        assert!(result.substrate_branch_dissipation_w.iter().sum::<f64>() > 0.0);
        assert!(result.coating_branch_dissipation_w.iter().sum::<f64>() > 0.0);
        assert!(load.dissipated_power_w > 0.0);
    }
}

#[test]
fn g1_exact_bare_load_pole_is_a_finite_joint_response_and_host_antiresonance() {
    let gate = CancelGate::new_clock_free();
    let image = image(&gate);
    let bare = image.prepare_coating_resonance(budget(), &gate).unwrap();
    let beta = image.force_ports()[0][0];
    let standalone = HarmonicModel::from_ode(
        &attachment(1.0, 1_000_000.0, 0.0, &[1.0]),
        HarmonicBudget {
            max_dimension: 2,
            max_dense_bytes: 10_000,
            max_dense_work: 10_000,
            relative_residual_tolerance: 1e-10,
        },
        || false,
    )
    .unwrap();
    assert!(matches!(
        standalone.factor(1000.0, || false),
        Err(HarmonicError::Singular { .. })
    ));
    let prepared = image
        .prepare_coating_resonance_with_load(
            ShellResonanceLoad {
                system: attachment(1.0, 1_000_000.0, 0.0, &[1.0]),
                ports: &[ShellHarmonicPort::Force(0)],
                band_rad_s: bare.source_band_rad_s(),
            },
            budget(),
            &gate,
        )
        .unwrap();
    let input = [C64::new(0.6, -0.8), C64::ZERO, C64::ZERO, C64::ZERO];
    let result = prepared
        .harmonic_factor(1000.0, forced_budget(), &gate)
        .unwrap()
        .response(&input, &gate)
        .unwrap();
    let load = result.load.as_ref().unwrap();
    let r = input[0].scale(-1e-6);
    let p = C64::new(0.0, -1000.0) * r;
    near(load.state[0], r);
    near(load.state[1], p);
    near(load.effort[0], input[0]);
    assert!(result.displacement_m_sqrt_kg[0].abs() < 1e-14 * r.abs() / beta.abs());
    assert!(result.velocity_m_sqrt_kg_per_s[0].abs() < 1e-14 * p.abs() / beta.abs());
    assert_eq!(load.dissipated_power_w, 0.0);
    near_real(
        load.mean_stored_energy_j,
        0.25 * (1_000_000.0 * r.norm_sq() + p.norm_sq()),
    );
    assert!(load.mean_stored_energy_j > 0.0);
    assert!(load.absorbed_port_power_w[0].abs() < 1e-14 * 1000.0 * load.mean_stored_energy_j);
    assert!(result.physical_relative_residual < 1e-10);
    assert!(result.original_equation_relative_residual < 1e-10);
}

#[test]
fn g3_mutual_signed_work_pressure_drives_and_multiple_inputs_share_the_full_factor() {
    let gate = CancelGate::new_clock_free();
    let image = image(&gate);
    let bare = image.prepare_coating_resonance(budget(), &gate).unwrap();
    let w = bare.equilibrium_angular_frequencies_rad_s().next().unwrap();
    let maps = [
        1.0,
        -0.35 * image.force_ports()[0][0] / image.force_ports()[1][0],
    ];
    let beta = maps[0] * image.force_ports()[0][0] + maps[1] * image.force_ports()[1][0];
    let mass = 0.5 / (beta * beta);
    let stiffness = mass * 1_000_000.0;
    let gamma = 15.0;
    let prepared = image
        .prepare_coating_resonance_with_load(
            ShellResonanceLoad {
                system: attachment(mass, stiffness, mass * gamma, &maps),
                ports: &[ShellHarmonicPort::Force(0), ShellHarmonicPort::Force(1)],
                band_rad_s: bare.source_band_rad_s(),
            },
            budget(),
            &gate,
        )
        .unwrap();
    let input = [
        C64::new(0.7, 0.3),
        C64::new(-0.1, 0.2),
        C64::new(200.0, 30.0),
        C64::new(50.0, -20.0),
    ];
    let omega = 1100.0;
    let factor = prepared
        .harmonic_factor(omega, forced_budget(), &gate)
        .unwrap();
    let result = factor.response(&input, &gate).unwrap();
    let force = input[0].scale(image.force_ports()[0][0])
        + input[1].scale(image.force_ports()[1][0])
        + input[2].scale(image.volume_shapes_m2_per_sqrt_kg()[0])
        + input[3].scale(image.cavity_volume_shapes_m2_per_sqrt_kg()[0]);
    let (q, r, p) = analytic(&image, omega, w, beta, mass, stiffness, gamma, force);
    source_channels(&image, &result, q, w);
    let load = result.load.as_ref().unwrap();
    near(load.state[0], r);
    near(load.state[1], p);
    assert!(load.absorbed_port_power_w[0] > 0.0);
    assert!(load.absorbed_port_power_w[1] < 0.0);
    near_real(
        load.absorbed_port_power_w.iter().sum(),
        load.dissipated_power_w,
    );
    for i in 0..2 {
        near(load.effort[i], r.scale(-maps[i] * stiffness));
        near(
            load.flow[i],
            C64::new(0.0, -omega) * q.scale(image.force_ports()[i][0]),
        );
        near_real(
            load.absorbed_port_power_w[i],
            0.5 * (load.flow[i].conj() * load.effort[i]).re,
        );
    }
    assert!(result.physical.port_power_w[2] * result.physical.port_power_w[3] < 0.0);
    for i in 0..4 {
        near_real(
            result.physical.port_power_w[i],
            0.5 * (input[i].conj() * result.physical.output[i]).re,
        );
    }
    let phase = C64::new(0.3, -0.4);
    let next_input: Vec<_> = input.iter().map(|&u| phase * u).collect();
    let next = factor.response(&next_input, &gate).unwrap();
    for (&actual, &before) in next.physical.state.iter().zip(&result.physical.state) {
        near(actual, phase * before);
    }
    near_real(
        next.physical.mean_stored_energy_j,
        phase.norm_sq() * result.physical.mean_stored_energy_j,
    );
    near_real(
        next.physical.dissipated_power_w,
        phase.norm_sq() * result.physical.dissipated_power_w,
    );
    assert_eq!(factor.response(&input, &gate).unwrap(), result);
}

#[test]
fn g0_forced_admission_cancel_and_replay_keep_the_selected_source_unchanged() {
    let gate = CancelGate::new_clock_free();
    let image = image(&gate);
    let bare = image.prepare_coating_resonance(budget(), &gate).unwrap();
    let prepared = image
        .prepare_coating_resonance_with_load(
            ShellResonanceLoad {
                system: attachment(1.0, 1_000_000.0, 10.0, &[1.0]),
                ports: &[ShellHarmonicPort::Force(0)],
                band_rad_s: (900.0, 1500.0),
            },
            budget(),
            &gate,
        )
        .unwrap();
    assert_eq!(prepared.source_band_rad_s(), (900.0, 1500.0));
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY, 899.0, 1501.0] {
        assert!(matches!(
            prepared.harmonic_factor(omega, forced_budget(), &gate),
            Err(ShellResonanceError::Invalid(_))
        ));
    }
    let mut small = forced_budget();
    small.maximum_payload_bytes = 0;
    assert!(matches!(
        prepared.harmonic_factor(1100.0, small, &gate),
        Err(ShellResonanceError::Budget {
            quantity: "forced physical observation bytes",
            ..
        })
    ));
    small = forced_budget();
    small.maximum_observation_work = 0;
    assert!(matches!(
        prepared.harmonic_factor(1100.0, small, &gate),
        Err(ShellResonanceError::Budget {
            quantity: "forced physical observation work",
            ..
        })
    ));
    let input = [C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO];
    let factor = prepared
        .harmonic_factor(1100.0, forced_budget(), &gate)
        .unwrap();
    let before = factor.response(&input, &gate).unwrap();
    assert!(factor.response(&input[..3], &gate).is_err());
    assert!(
        factor
            .response(&[C64::new(f64::NAN, 0.0); 4], &gate)
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        prepared
            .harmonic_factor(1100.0, forced_budget(), &cancelled)
            .is_err()
    );
    assert!(factor.response(&input, &cancelled).is_err());
    assert_eq!(factor.response(&input, &gate).unwrap(), before);
    let zero = factor.response(&[C64::ZERO; 4], &gate).unwrap();
    assert!(zero.physical.state.iter().all(|&x| x == C64::ZERO));
    assert_eq!(zero.physical.mean_stored_energy_j, 0.0);
    assert_eq!(zero.physical.dissipated_power_w, 0.0);
    assert_eq!(zero.physical_relative_residual, 0.0);
    let w = bare.equilibrium_angular_frequencies_rad_s().next().unwrap();
    let response = bare
        .harmonic_factor(1100.0, forced_budget(), &gate)
        .unwrap()
        .response(&input, &gate)
        .unwrap();
    let q =
        C64::from_re(image.force_ports()[0][0]) / body_operator(&image, C64::new(0.0, -1100.0), w);
    source_channels(&image, &response, q, w);
    assert!(response.load.is_none());
    assert_eq!(prepared.dimension(), bare.dimension() + 2);
}

#[test]
fn g1_two_shell_coordinates_keep_full_cross_memory_and_nonproportional_loss() {
    let gate = CancelGate::new_clock_free();
    let image = image_with_modes(&gate, 2);
    let bare = image.prepare_coating_resonance(budget(), &gate).unwrap();
    let w: Vec<_> = bare.equilibrium_angular_frequencies_rad_s().collect();
    let beta = &image.force_ports()[0];
    let mass = 0.5 / (beta[0] * beta[0] + beta[1] * beta[1]);
    let stiffness = mass * 1_000_000.0;
    let gamma = 20.0;
    let prepared = image
        .prepare_coating_resonance_with_load(
            ShellResonanceLoad {
                system: attachment(mass, stiffness, mass * gamma, &[1.0]),
                ports: &[ShellHarmonicPort::Force(0)],
                band_rad_s: bare.source_band_rad_s(),
            },
            budget(),
            &gate,
        )
        .unwrap();
    let omega = 1500.0;
    let s = C64::new(0.0, -omega);
    let input = [
        C64::new(0.7, -0.2),
        C64::new(-0.15, 0.3),
        C64::from_re(40.0),
        C64::from_re(-15.0),
    ];
    let result = prepared
        .harmonic_factor(omega, forced_budget(), &gate)
        .unwrap()
        .response(&input, &gate)
        .unwrap();
    let all = branches(&image);
    assert!(all.iter().any(|(row, _)| row[0] * row[1] != 0.0));
    // Source loss is carried by the retained memory here; the independent
    // damped attachment couples both host coordinates through the same state.
    assert_ne!(beta[0] * beta[1], 0.0);
    let a = s * s + s.scale(gamma);
    let load_operator = a + C64::from_re(stiffness / mass);
    let mut matrix = [C64::ZERO; 4];
    let mut force = [C64::ZERO; 2];
    for i in 0..2 {
        force[i] = input[0].scale(image.force_ports()[0][i])
            + input[1].scale(image.force_ports()[1][i])
            + input[2].scale(image.volume_shapes_m2_per_sqrt_kg()[i])
            + input[3].scale(image.cavity_volume_shapes_m2_per_sqrt_kg()[i]);
        for j in 0..2 {
            let mut body = s.scale(image.damping().matrix_per_s()[2 * i + j]);
            if i == j {
                body = body + s * s + C64::from_re(w[i] * w[i]);
            }
            for (row, tau) in &all {
                body = body + s.scale(tau * row[i] * row[j]) / (C64::ONE + s.scale(*tau));
            }
            matrix[2 * i + j] = body * load_operator + a.scale(stiffness * beta[i] * beta[j]);
        }
    }
    let determinant = matrix[0] * matrix[3] - matrix[1] * matrix[2];
    let q = [
        (matrix[3] * force[0] - matrix[1] * force[1]) * load_operator / determinant,
        (matrix[0] * force[1] - matrix[2] * force[0]) * load_operator / determinant,
    ];
    let v = [s * q[0], s * q[1]];
    let support = q[0].scale(beta[0]) + q[1].scale(beta[1]);
    let r = -a * support / load_operator;
    let p = s.scale(stiffness) * support / load_operator;
    let load = result.load.as_ref().unwrap();
    near(load.state[0], r);
    near(load.state[1], p);
    let mut energy = 0.25 * (stiffness * r.norm_sq() + p.norm_sq() / mass);
    let mut heat = 0.5 * gamma * p.norm_sq() / mass;
    near_real(load.mean_stored_energy_j, energy);
    near_real(load.dissipated_power_w, heat);
    let mut mechanical = 0.0;
    let mut background = 0.0;
    for i in 0..2 {
        near(result.displacement_m_sqrt_kg[i], q[i]);
        near(result.velocity_m_sqrt_kg_per_s[i], v[i]);
        mechanical += 0.25 * (v[i].norm_sq() + q[i].scale(w[i]).norm_sq());
        for j in 0..2 {
            background += 0.5 * (v[i].conj() * v[j]).re * image.damping().matrix_per_s()[2 * i + j];
        }
    }
    near_real(result.mean_mechanical_storage_j, mechanical);
    near_real(result.background_dissipation_w, background);
    energy += mechanical;
    heat += background;
    for (j, (row, tau)) in all.iter().enumerate() {
        let extension = q[0].scale(row[0]) + q[1].scale(row[1]);
        let z = extension / (C64::ONE + s.scale(*tau));
        near(result.physical.state[4 + j], z);
        energy += 0.25 * (extension - z).norm_sq();
        heat += 0.5 * (extension - z).norm_sq() / tau;
    }
    near_real(result.physical.mean_stored_energy_j, energy);
    near_real(result.physical.dissipated_power_w, heat);
    near_real(result.physical.supplied_power_w, heat);
    for (i, row) in image.force_ports().iter().enumerate() {
        near(
            result.force_port_displacement_m[i],
            q[0].scale(row[0]) + q[1].scale(row[1]),
        );
        near(
            result.force_port_velocity_m_s()[i],
            v[0].scale(row[0]) + v[1].scale(row[1]),
        );
    }
    assert!(result.original_equation_relative_residual < 1e-10);
    assert!(result.physical_relative_residual < 1e-10);
}


fn strain_budget(evaluations: usize) -> ShellCoatingStrainBudget {
    ShellCoatingStrainBudget {
        maximum_bytes: 8_000_000,
        maximum_work: 100_000_000,
        evaluations,
    }
}

#[test]
fn g1_actual_forcing_checks_the_whole_cycle_before_publishing_and_reuses_the_factor() {
    let gate = CancelGate::new_clock_free();
    let image = image(&gate);
    let prepared = image
        .prepare_coating_resonance_with_load(
            ShellResonanceLoad {
                system: attachment(1.0, 1_000_000.0, 20.0, &[1.0]),
                ports: &[ShellHarmonicPort::Force(0)],
                band_rad_s: (900.0, 1500.0),
            },
            budget(),
            &gate,
        )
        .unwrap();
    let linear = prepared
        .harmonic_factor(1100.0, forced_budget(), &gate)
        .unwrap();
    let unit_input = [C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO];
    let unit = linear.response(&unit_input, &gate).unwrap();
    let mut observer = image
        .prepare_coating_strain(strain_budget(5), &gate)
        .unwrap();
    let unit_strain = observer.evaluate(&[C64::ONE], &gate).unwrap();
    let cap = unit_strain.sufficient_common_amplitude_scale.unwrap();
    // Each quadrature snapshot is inside the source envelope, while the
    // intermediate phase reaches sqrt(2) times either snapshot's strain.
    // Determine the actual force from this retained physical coordinate;
    // never use a free-pole normalization as the forced amplitude.
    let target = C64::new(0.8 * cap, 0.8 * cap);
    let input = [
        target / unit.displacement_m_sqrt_kg[0],
        C64::ZERO,
        C64::ZERO,
        C64::ZERO,
    ];
    let raw = linear.response(&input, &gate).unwrap();
    let q = raw.displacement_m_sqrt_kg[0];
    near(q, target);
    for snapshot in [q.re, q.im] {
        let report = observer.evaluate_snapshot(&[snapshot], &gate).unwrap();
        assert!(report.within_authored_envelope);
        near_real(report.films[0].limit_ratio, 0.8);
    }
    let intermediate = observer
        .evaluate_snapshot(&[(q.re + q.im) / 2.0_f64.sqrt()], &gate)
        .unwrap();
    assert!(!intermediate.within_authored_envelope);
    let cycle = observer.evaluate(&[q], &gate).unwrap();
    assert!(!cycle.within_authored_envelope);
    near_real(
        cycle.films[0].frobenius_upper_bound,
        intermediate.films[0].frobenius_upper_bound,
    );

    assert!(
        prepared
            .harmonic_factor_with_coating_strain(
                1100.0,
                forced_budget(),
                ShellCoatingStrainBudget {
                    maximum_work: 0,
                    ..strain_budget(4)
                },
                &gate,
            )
            .is_err()
    );
    let mut factor = prepared
        .harmonic_factor_with_coating_strain(
            1100.0,
            forced_budget(),
            strain_budget(4),
            &gate,
        )
        .unwrap();
    assert_eq!(factor.omega_rad_s(), 1100.0);
    assert_eq!(factor.strain_resources().coated_elements, 18);
    assert!(factor.response(&input[..3], &gate).is_err());
    assert!(factor.response(&[C64::new(f64::NAN, 0.0); 4], &gate).is_err());
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(factor.response(&input, &stopped).is_err());
    assert_eq!(factor.completed_strain_evaluations(), 0);
    assert!(matches!(
        factor.response(&input, &gate),
        Err(ShellResonanceError::Invalid(
            "forced shell motion exceeds the sufficient original film strain envelope"
        ))
    ));
    assert_eq!(factor.completed_strain_evaluations(), 1);
    let safe_input = input.map(|u| u.scale(0.5));
    let expected = linear.response(&safe_input, &gate).unwrap();
    let accepted = factor.response(&safe_input, &gate).unwrap();
    assert_eq!(accepted.shell, expected);
    assert!(accepted.coating_strain.within_authored_envelope);
    let film = &accepted.coating_strain.films[0];
    assert!(std::ptr::eq(film.source, &image.source().coatings()[0]));
    assert_eq!(film.source.input().material.strain_limit(), 0.0003);
    near_real(film.limit_ratio, 0.5 * cycle.films[0].limit_ratio);
    let w = prepared.equilibrium_angular_frequencies_rad_s().next().unwrap();
    source_channels(&image, &accepted.shell, q.scale(0.5), w);
    assert!(accepted.shell.load.as_ref().unwrap().mean_stored_energy_j > 0.0);
    let rotated_input = safe_input.map(|u| C64::new(0.0, 1.0) * u);
    let rotated = factor.response(&rotated_input, &gate).unwrap();
    near_real(
        rotated.coating_strain.films[0].frobenius_upper_bound,
        film.frobenius_upper_bound,
    );
    assert_eq!(factor.response(&safe_input, &gate).unwrap().shell, expected);
    assert_eq!(factor.completed_strain_evaluations(), 4);
    assert!(matches!(
        factor.response(&safe_input, &gate),
        Err(ShellResonanceError::Invalid(
            "forced shell coating strain query allowance exhausted"
        ))
    ));
    assert_eq!(linear.response(&unit_input, &gate).unwrap(), unit);
}

#[test]
fn g3_assembled_forcing_keeps_every_original_component_film_and_its_own_limit() {
    let gate = CancelGate::new_clock_free();
    let limits = [0.0003, 0.000003];
    let mut components = Vec::new();
    let mut rows = Vec::new();
    for (i, &limit) in limits.iter().enumerate() {
        let (shell, directors, constraints) = fixture::coupon(&gate);
        let mut film = fixture::layer(CoatingSide::Positive, (20.0 + 5.0 * i as f64) * 1e-6);
        // Different caller-authored envelopes on the original physical films.
        // The source fit, its conditioning and its Poisson ratio stay unchanged.
        film.material = ResolvedCoating::from_material(
            film.material.source().clone(),
            film.material.poisson_ratio(),
            limit,
        )
        .unwrap();
        let shell = shell
            .with_coatings(vec![film], &gate)
            .unwrap()
            .with_coating_relaxation(
                CoatingRelaxationOptions {
                    terms: 16,
                    storage_relative_tolerance: 0.001,
                    loss_relative_tolerance: 0.001,
                    max_memory_coordinates: 256,
                },
                &gate,
            )
            .unwrap()
            .with_substrate_relaxation(
                SubstrateRelaxationOptions {
                    storage_reference: SubstrateStorageReference::Relaxed,
                    max_memory_coordinates: 256,
                },
                &gate,
            )
            .unwrap();
        let facets = |side| {
            (0..shell.mesh().tris.len())
                .map(|triangle| OrientedShellFacet { triangle, side })
                .collect()
        };
        let surfaces = ShellPressureSurfaces {
            exterior: facets(ShellFace::Positive),
            cavity: facets(ShellFace::Negative),
        };
        let shell = shell.with_pressure_surfaces(surfaces, &gate).unwrap();
        let mut row = vec![0.0; 6 * shell.mesh().nodes.len()];
        row[6 * 5 + 2] = 1.0;
        rows.push(row);
        components.push(MaterialShellComponent {
            name: format!("coupon-{i}"),
            assembly_budget: fixture::budget(&shell),
            specimen: shell,
            directors,
            constraints,
            mass: ShellMass::Consistent,
            drilling_alpha: 1e-3,
        });
    }
    let ports: Vec<_> = (0..2)
        .map(|driven| ShellAssemblyForcePort {
            components: rows
                .iter()
                .enumerate()
                .map(|(i, row)| {
                    if i == driven {
                        row.clone()
                    } else {
                        vec![0.0; row.len()]
                    }
                })
                .collect(),
        })
        .collect();
    // Two independently clamped components share the original assembled basis.
    // No seam is authored, so force port 1 excites only component 1.
    let source = MaterialShellAssembly {
        components,
        edge_ties: vec![],
        face_ties: vec![],
        tie_budget: TieBudget {
            max_segments: 32,
            max_slave_dofs: 24,
            gap_tolerance_m: 1e-10,
            angle_tolerance_rad: 1e-10,
        },
    };
    let image = source
        .reduce(
            &ShellAssemblyReductionOptions {
                eigenvalue_window: (0.0, 1e9),
                n_modes: 2,
                projection_budget: ShellProjectionBudget {
                    max_modes: 2,
                    max_projection_terms: 100_000_000,
                    reference_resistance_n_s_m: 1.0,
                },
                max_total_nodes: 32,
                max_total_triangles: 36,
            },
            &ports,
            &gate,
        )
        .unwrap();
    let prepared = image.prepare_coating_resonance(budget(), &gate).unwrap();
    let linear = prepared
        .harmonic_factor(1100.0, forced_budget(), &gate)
        .unwrap();
    let mut observer = image
        .prepare_coating_strain(strain_budget(2), &gate)
        .unwrap();
    let second_only = [C64::ZERO, C64::ONE, C64::ZERO, C64::ZERO];
    let raw = linear.response(&second_only, &gate).unwrap();
    let strain = observer.evaluate(&raw.displacement_m_sqrt_kg, &gate).unwrap();
    assert_eq!(strain.films.len(), 2);
    assert!(strain.films[1].limit_ratio > 0.0);
    let over_scale = 1.2 / strain.films[1].limit_ratio;
    assert!(over_scale * strain.films[0].limit_ratio < 1.0);
    // A default 0.0003 limit would incorrectly accept the second component.
    assert!(over_scale * strain.films[1].frobenius_upper_bound < limits[0]);
    let mut factor = prepared
        .harmonic_factor_with_coating_strain(
            1100.0,
            forced_budget(),
            strain_budget(3),
            &gate,
        )
        .unwrap();
    let excessive = second_only.map(|u| u.scale(over_scale));
    assert!(matches!(
        factor.response(&excessive, &gate),
        Err(ShellResonanceError::Invalid(
            "forced shell motion exceeds the sufficient original film strain envelope"
        ))
    ));
    let safe = excessive.map(|u| u.scale(0.5));
    let accepted = factor.response(&safe, &gate).unwrap();
    near_real(accepted.coating_strain.films[1].limit_ratio, 0.6);
    assert_eq!(accepted.shell, linear.response(&safe, &gate).unwrap());
    let both = [C64::ONE, C64::new(0.5, -0.5), C64::ZERO, C64::ZERO];
    let raw = linear.response(&both, &gate).unwrap();
    let strain = observer.evaluate(&raw.displacement_m_sqrt_kg, &gate).unwrap();
    let safe_scale = 0.5 * strain.sufficient_common_amplitude_scale.unwrap();
    let both = both.map(|u| u.scale(safe_scale));
    let accepted = factor.response(&both, &gate).unwrap();
    assert_eq!(accepted.shell, linear.response(&both, &gate).unwrap());
    assert_eq!(accepted.coating_strain.films.len(), 2);
    for (i, report) in accepted.coating_strain.films.iter().enumerate() {
        let part = &image.source().components[i];
        assert_eq!(report.component_index, Some(i));
        assert_eq!(report.component_name, Some(part.name.as_str()));
        assert_eq!(report.film_index, 0);
        assert_eq!(report.covered_elements, part.specimen.mesh().tris.len());
        assert!(std::ptr::eq(report.source, &part.specimen.coatings()[0]));
        assert_eq!(report.source.input().material.strain_limit(), limits[i]);
        assert!(report.frobenius_upper_bound > 0.0);
        near_real(report.limit_ratio, report.frobenius_upper_bound / limits[i]);
        assert!(report.limit_ratio <= 1.0);
    }
    assert_eq!(factor.completed_strain_evaluations(), 3);
    assert_eq!(factor.strain_resources().coated_elements, 36);
}

#[path = "coating_dynamic_response/descriptor.rs"]
mod descriptor;
