//! Joint physical shell/attachment poles, without isolated impedance inversion.
#![deny(unsafe_code)]

#[path = "coating_shell.rs"]
mod fixture;

use fs_couple::shear_shell::{
    CoatingSide, MaterialShellImage, ShellProjectionBudget, ShellReductionOptions,
    coating::{
        poles::{
            BeynOptions, PoleOptions, ResonanceRegion, ShellHarmonicPort, ShellResonanceBudget,
            ShellResonanceError, ShellResonanceLoad,
        },
        relaxation::CoatingRelaxationOptions,
    },
    constant_loss::relaxation::{SubstrateRelaxationOptions, SubstrateStorageReference},
    pressure::ShellPressureSurfaces,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_phs::{PortHamiltonian, QuadraticStorage};
use fs_plate::shell::{
    mitc::ShellMass,
    pressure::{OrientedShellFacet, ShellFace},
};

fn image(gate: &CancelGate) -> MaterialShellImage {
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
                n_modes: 1,
                assembly_budget: fixture::budget(&shell),
                projection_budget: ShellProjectionBudget {
                    max_modes: 1,
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
fn options(w: f64) -> PoleOptions {
    PoleOptions {
        spectral_scale: w,
        step_tolerance: 1e-13,
        residual_tolerance: 1e-14,
        max_iterations: 64,
        max_dimension: 128,
        max_dense_bytes: 16_000_000,
        max_work: 1_000_000_000,
    }
}
fn contour(w: f64) -> BeynOptions {
    BeynOptions {
        nodes: 128,
        probe_columns: 3,
        rank_relative_tolerance: 1e-10,
        count_tolerance: 1e-5,
        quadrature_tolerance: 1e-5,
        contour_min_reciprocal_condition: 1e-16,
        linear_residual_tolerance: 1e-12,
        svd_tolerance: 1e-13,
        svd_max_sweeps: 100,
        eigen_max_sweeps: 120,
        max_dimension: 128,
        max_dense_bytes: 16_000_000,
        max_work: 10_000_000_000,
        refinement: options(w),
    }
}

// A real mass on a spring driven by a moving spring support. The state is
// [relative spring displacement, mass momentum]; input is support velocity.
// The damper acts between the mass and ground. No fitted modal damping is used.
// Multiple support maps create an original full mutual port impedance.
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

// Independent scalar body/attachment dispersion, used only as a test oracle.
// It multiplies out the load denominator, so its bare poles are not divided by.
fn joint_dispersion(
    s: C64,
    image: &MaterialShellImage,
    w: f64,
    wl: f64,
    gamma: f64,
    coupling_stiffness: f64,
) -> (C64, C64) {
    let mut body = s * s + s.scale(image.damping().matrix_per_s()[0]) + C64::from_re(w * w);
    let mut derivative = s.scale(2.0) + C64::from_re(image.damping().matrix_per_s()[0]);
    for (row, tau) in branches(image) {
        let d = C64::ONE + s.scale(tau);
        body = body + s.scale(tau * row[0] * row[0]) / d;
        derivative = derivative + C64::from_re(tau * row[0] * row[0]) / (d * d);
    }
    let load = s * s + s.scale(gamma) + C64::from_re(wl * wl);
    (
        body * load + (s * s + s.scale(gamma)).scale(coupling_stiffness),
        derivative * load
            + body * (s.scale(2.0) + C64::from_re(gamma))
            + (s.scale(2.0) + C64::from_re(gamma)).scale(coupling_stiffness),
    )
}
fn root(
    mut s: C64,
    image: &MaterialShellImage,
    w: f64,
    wl: f64,
    gamma: f64,
    stiffness: f64,
) -> C64 {
    for _ in 0..64 {
        let (f, df) = joint_dispersion(s, image, w, wl, gamma, stiffness);
        let correction = f / df;
        s = s - correction;
        if correction.abs() < 1e-14 * w {
            return s;
        }
    }
    panic!("independent joint dispersion did not converge");
}
fn near(a: C64, b: C64, scale: f64) {
    assert!(
        (a - b).abs() < 2e-10 * scale,
        "{a:?} versus {b:?}, scale={scale}"
    );
}

#[test]
fn g1_dynamic_attachment_splits_hereditary_shell_resonance_and_keeps_both_power_accounts() {
    let gate = CancelGate::new_clock_free();
    let image = image(&gate);
    let bare = image.prepare_coating_resonance(budget(), &gate).unwrap();
    let w = bare.equilibrium_angular_frequencies_rad_s().next().unwrap();
    let beta = image.force_ports()[0][0];
    let mass = 0.5 / (beta * beta);
    let stiffness = mass * w * w;
    let gamma = 0.015 * w;
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
    assert_eq!(prepared.dimension(), bare.dimension() + 2);
    assert_eq!(prepared.load_state_count(), 2);
    assert_eq!(prepared.coating_memory_count(), bare.coating_memory_count());
    assert_eq!(
        prepared.substrate_memory_count(),
        bare.substrate_memory_count()
    );
    let mut rates = Vec::new();
    for fraction in [0.72, 1.43] {
        let initial = C64::new(-0.015 * w, -fraction * w);
        let expected = root(initial, &image, w, w, gamma, beta * beta * stiffness);
        let seed = prepared
            .displacement_seed(initial, &[C64::ONE], &gate)
            .unwrap();
        let pole = prepared
            .resonance_pole(initial, &seed, &options(w), &gate)
            .unwrap();
        near(pole.mode.laplace_pole_per_s, expected, w);
        let load = pole.load.as_ref().unwrap();
        let s = pole.mode.laplace_pole_per_s;
        let flow = pole.velocity_m_sqrt_kg_per_s[0].scale(beta);
        near(load.flow[0], flow, flow.abs());
        near(
            load.effort[0],
            load.state[0].scale(-stiffness),
            load.effort[0].abs(),
        );
        near(
            s * load.state[0],
            load.state[1].scale(1.0 / mass) - flow,
            flow.abs(),
        );
        near(
            s * load.state[1],
            load.state[0].scale(-stiffness) - load.state[1].scale(gamma),
            load.effort[0].abs(),
        );
        let h = 0.5 * (stiffness * load.state[0].norm_sq() + load.state[1].norm_sq() / mass);
        let d = gamma * load.state[1].norm_sq() / mass;
        assert!((load.stored_energy_j - h).abs() < 1e-12 * h);
        assert!((load.dissipated_power_w - d).abs() < 1e-12 * d);
        assert!(load.stored_energy_j > 0.0 && load.dissipated_power_w > 0.0);
        assert!(pole.substrate_branch_dissipation_w.iter().sum::<f64>() > 0.0);
        assert!(pole.coating_branch_dissipation_w.iter().sum::<f64>() > 0.0);
        let scale = pole.mode.dissipated_power_w;
        assert!(load.decay_power_residual_w.abs() < 1e-9 * scale);
        assert!(pole.shell_decay_power_residual_w.abs() < 1e-9 * scale);
        assert!(pole.decay_power_residual_w.abs() < 1e-9 * scale);
        assert!(load.relative_equation_residual < 1e-11);
        rates.push(s);
    }
    assert!((rates[0].im - rates[1].im).abs() > 0.4 * w);
    let spectrum = prepared
        .resonance_spectrum(
            ResonanceRegion {
                center: C64::new(-0.03 * w, -1.1 * w),
                real_radius: 0.2 * w,
                imaginary_radius: 0.7 * w,
            },
            &contour(w),
            &gate,
        )
        .unwrap();
    assert_eq!(spectrum.argument_count, 2);
    assert_eq!(spectrum.poles.len(), 2);
    for s in rates {
        assert!(
            spectrum
                .poles
                .iter()
                .any(|p| (p.mode.laplace_pole_per_s - s).abs() < 1e-9 * w)
        );
    }
}

#[test]
fn g3_mutual_load_uses_signed_original_force_projections_without_diagonalizing() {
    let gate = CancelGate::new_clock_free();
    let image = image(&gate);
    let bare = image.prepare_coating_resonance(budget(), &gate).unwrap();
    let w = bare.equilibrium_angular_frequencies_rad_s().next().unwrap();
    let maps = [0.8, -0.3];
    let beta = maps[0] * image.force_ports()[0][0] + maps[1] * image.force_ports()[1][0];
    let mass = 0.3 / (beta * beta);
    let k = mass * w * w;
    let gamma = 0.02 * w;
    let prepared = image
        .prepare_coating_resonance_with_load(
            ShellResonanceLoad {
                system: attachment(mass, k, mass * gamma, &maps),
                ports: &[ShellHarmonicPort::Force(0), ShellHarmonicPort::Force(1)],
                band_rad_s: bare.source_band_rad_s(),
            },
            budget(),
            &gate,
        )
        .unwrap();
    let start = C64::new(-0.02 * w, -0.8 * w);
    let expected = root(start, &image, w, w, gamma, beta * beta * k);
    let seed = prepared
        .displacement_seed(start, &[C64::ONE], &gate)
        .unwrap();
    let pole = prepared
        .resonance_pole(start, &seed, &options(w), &gate)
        .unwrap();
    near(pole.mode.laplace_pole_per_s, expected, w);
    let load = pole.load.unwrap();
    for i in 0..2 {
        near(
            load.effort[i],
            load.state[0].scale(-maps[i] * k),
            load.effort[i].abs(),
        );
        near(
            load.flow[i],
            pole.velocity_m_sqrt_kg_per_s[0].scale(image.force_ports()[i][0]),
            load.flow[i].abs(),
        );
    }
    assert!(load.absorbed_port_power_w.iter().all(|x| *x != 0.0));
}

#[test]
fn g1_zero_net_volume_coupling_retains_independent_load_resonance_at_its_bare_pole() {
    let gate = CancelGate::new_clock_free();
    let image = image(&gate);
    assert_eq!(
        image.volume_shapes_m2_per_sqrt_kg()[0],
        -image.cavity_volume_shapes_m2_per_sqrt_kg()[0]
    );
    let bare = image.prepare_coating_resonance(budget(), &gate).unwrap();
    let w = bare.equilibrium_angular_frequencies_rad_s().next().unwrap() * 1.5;
    let gamma = 0.02 * w;
    let prepared = image
        .prepare_coating_resonance_with_load(
            ShellResonanceLoad {
                system: attachment(1.0, w * w, gamma, &[1.0, 1.0]),
                ports: &[
                    ShellHarmonicPort::ExteriorVolume,
                    ShellHarmonicPort::CavityVolume,
                ],
                band_rad_s: bare.source_band_rad_s(),
            },
            budget(),
            &gate,
        )
        .unwrap();
    let expected = C64::new(-0.5 * gamma, -(w * w - 0.25 * gamma * gamma).sqrt());
    let mut state = vec![C64::ZERO; prepared.dimension()];
    state[bare.dimension()] = C64::ONE;
    state[bare.dimension() + 1] = expected;
    let pole = prepared
        .resonance_pole(expected, &state, &options(w), &gate)
        .unwrap();
    near(pole.mode.laplace_pole_per_s, expected, w);
    assert!(pole.displacement_m_sqrt_kg[0].abs() < 1e-14);
    assert!(pole.exterior_volume_velocity_m3_s().abs() < 1e-14);
    let load = pole.load.unwrap();
    assert_eq!(load.state.len(), 2);
    assert!(load.stored_energy_j > 0.0);
    assert!(load.dissipated_power_w > 0.0);
    assert!(load.absorbed_port_power_w.iter().all(|p| p.abs() < 1e-20));
}

#[test]
fn g0_joint_load_admission_preserves_source_on_capacity_band_port_and_cancellation_refusal() {
    let gate = CancelGate::new_clock_free();
    let image = image(&gate);
    let bare = image.prepare_coating_resonance(budget(), &gate).unwrap();
    let make = || attachment(1.0, 1.0, 0.1, &[1.0]);
    let mut small = budget();
    small.maximum_coordinates = bare.dimension() + 1;
    assert!(matches!(
        image.prepare_coating_resonance_with_load(
            ShellResonanceLoad {
                system: make(),
                ports: &[ShellHarmonicPort::Force(0)],
                band_rad_s: bare.source_band_rad_s(),
            },
            small,
            &gate
        ),
        Err(ShellResonanceError::Budget {
            quantity: "complete coordinates",
            ..
        })
    ));
    for (ports, band) in [
        (vec![ShellHarmonicPort::Force(8)], bare.source_band_rad_s()),
        (vec![], bare.source_band_rad_s()),
        (vec![ShellHarmonicPort::Force(0)], (0.0, 1.0)),
    ] {
        assert!(
            image
                .prepare_coating_resonance_with_load(
                    ShellResonanceLoad {
                        system: make(),
                        ports: &ports,
                        band_rad_s: band,
                    },
                    budget(),
                    &gate
                )
                .is_err()
        );
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        image
            .prepare_coating_resonance_with_load(
                ShellResonanceLoad {
                    system: make(),
                    ports: &[ShellHarmonicPort::Force(0)],
                    band_rad_s: bare.source_band_rad_s(),
                },
                budget(),
                &cancelled
            )
            .is_err()
    );
    let after = image.prepare_coating_resonance(budget(), &gate).unwrap();
    assert_eq!(bare.dimension(), after.dimension());
    assert_eq!(bare.source_band_rad_s(), after.source_band_rad_s());
}
