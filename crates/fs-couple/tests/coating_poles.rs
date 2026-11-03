//! Original shell material poles: independent dispersion and physical accounts.
#![deny(unsafe_code)]

#[path = "coating_shell.rs"]
mod fixture;

use fs_couple::shear_shell::{
    CoatingSide, MaterialShellImage, ResolvedShearShell, ShellProjectionBudget,
    ShellReductionOptions,
    assembly::{
        MaterialShellAssembly, MaterialShellComponent, ShellAssemblyForcePort,
        ShellAssemblyReductionOptions,
    },
    coating::{
        poles::{
            BeynOptions, PoleOptions, PreparedShellResonance, ResonanceRegion,
            ShellResonanceBudget, ShellResonanceError, ShellResonancePole,
        },
        relaxation::CoatingRelaxationOptions,
    },
    constant_loss::relaxation::{SubstrateRelaxationOptions, SubstrateStorageReference},
    pressure::ShellPressureSurfaces,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_plate::shell::{
    junction::TieBudget,
    mitc::{NodalDirectors, ShellConstraint, ShellMass},
    pressure::{OrientedShellFacet, ShellFace},
};

fn source(
    gate: &CancelGate,
    film_memory: bool,
    substrate_memory: bool,
) -> (ResolvedShearShell, NodalDirectors, Vec<ShellConstraint>) {
    let (shell, directors, constraints) = fixture::coupon(gate);
    let mut shell = shell
        .with_coatings(vec![fixture::layer(CoatingSide::Positive, 20e-6)], gate)
        .unwrap();
    if film_memory {
        shell = shell
            .with_coating_relaxation(
                CoatingRelaxationOptions {
                    terms: 16,
                    storage_relative_tolerance: 0.001,
                    loss_relative_tolerance: 0.001,
                    max_memory_coordinates: 256,
                },
                gate,
            )
            .unwrap();
    }
    if substrate_memory {
        shell = shell
            .with_substrate_relaxation(
                SubstrateRelaxationOptions {
                    storage_reference: SubstrateStorageReference::Relaxed,
                    max_memory_coordinates: 256,
                },
                gate,
            )
            .unwrap();
    }
    let facets = |side| {
        (0..shell.mesh().tris.len())
            .map(|triangle| OrientedShellFacet { triangle, side })
            .collect()
    };
    let surfaces = ShellPressureSurfaces {
        exterior: facets(ShellFace::Positive),
        cavity: facets(ShellFace::Negative),
    };
    (
        shell.with_pressure_surfaces(surfaces, gate).unwrap(),
        directors,
        constraints,
    )
}

fn ports(shell: &ResolvedShearShell) -> Vec<Vec<f64>> {
    [5, 6]
        .map(|node| {
            let mut row = vec![0.0; 6 * shell.mesh().nodes.len()];
            row[6 * node + 2] = 1.0;
            row
        })
        .to_vec()
}
fn projection_budget() -> ShellProjectionBudget {
    ShellProjectionBudget {
        max_modes: 3,
        max_projection_terms: 100_000_000,
        reference_resistance_n_s_m: 1.0,
    }
}
fn image(gate: &CancelGate, modes: usize, film: bool, substrate: bool) -> MaterialShellImage {
    let (shell, directors, constraints) = source(gate, film, substrate);
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
                projection_budget: projection_budget(),
            },
            &ports(&shell),
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
fn options(scale: f64) -> PoleOptions {
    PoleOptions {
        spectral_scale: scale,
        step_tolerance: 1e-13,
        residual_tolerance: 1e-13,
        max_iterations: 48,
        max_dimension: 128,
        max_dense_bytes: 16_000_000,
        max_work: 1_000_000_000,
    }
}
fn contour_options(scale: f64, nodes: usize) -> BeynOptions {
    BeynOptions {
        nodes,
        probe_columns: 4,
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
        refinement: options(scale),
    }
}
fn memory(image: &MaterialShellImage) -> Vec<(&[f64], f64)> {
    let substrate = image.substrate_relaxation().unwrap();
    let films = image.coating_relaxation().unwrap();
    substrate
        .coupling_rows_per_s()
        .iter()
        .zip(substrate.relaxation_times_s())
        .chain(
            films
                .coupling_rows_per_s()
                .iter()
                .zip(films.relaxation_times_s()),
        )
        .map(|(row, &tau)| (row.as_slice(), tau))
        .collect()
}
fn near(a: C64, b: C64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:?} versus {b:?}; scale={scale}"
    );
}

// Independent scalar Schur elimination of every ORIGINAL Maxwell arm. The
// production pencil retains z and never divides by these denominators.
fn dispersion(s: C64, omega: f64, damping: f64, branches: &[(&[f64], f64)]) -> (C64, C64) {
    let mut value = s * s + s.scale(damping) + C64::from_re(omega * omega);
    let mut derivative = s.scale(2.0) + C64::from_re(damping);
    for &(row, tau) in branches {
        let denominator = C64::ONE + s.scale(tau);
        value = value + s.scale(tau * row[0] * row[0]) / denominator;
        derivative = derivative + C64::from_re(tau * row[0] * row[0]) / (denominator * denominator);
    }
    (value, derivative)
}
fn scalar_root(start: C64, omega: f64, image: &MaterialShellImage) -> C64 {
    let branches = memory(image);
    let mut s = start;
    for _ in 0..32 {
        let (value, derivative) =
            dispersion(s, omega, image.damping().matrix_per_s()[0], &branches);
        let correction = value / derivative;
        s = s - correction;
        if correction.abs() < 2e-15 * omega {
            break;
        }
    }
    assert!(
        dispersion(s, omega, image.damping().matrix_per_s()[0], &branches)
            .0
            .abs()
            < 1e-10 * omega * omega
    );
    s
}
fn first_mode(prepared: &PreparedShellResonance<'_>, gate: &CancelGate) -> ShellResonancePole {
    let omega = prepared
        .equilibrium_angular_frequencies_rad_s()
        .next()
        .unwrap();
    let rate = C64::new(-0.01 * omega, -1.01 * omega);
    let mut displacement = vec![C64::ZERO; prepared.mode_count()];
    displacement[0] = C64::ONE;
    let seed = prepared
        .displacement_seed(rate, &displacement, gate)
        .unwrap();
    prepared
        .resonance_pole(rate, &seed, &options(omega), gate)
        .unwrap()
}

#[test]
fn g1_scalar_causal_shell_pole_matches_independent_dispersion_and_physical_decay() {
    let gate = CancelGate::new_clock_free();
    let image = image(&gate, 1, true, true);
    let prepared = image.prepare_coating_resonance(budget(), &gate).unwrap();
    let omega = prepared
        .equilibrium_angular_frequencies_rad_s()
        .next()
        .unwrap();
    let pole = first_mode(&prepared, &gate);
    let expected = scalar_root(C64::new(-0.01 * omega, -omega), omega, &image);
    near(pole.mode.laplace_pole_per_s, expected, omega, 1e-11);
    assert!(pole.mode.laplace_pole_per_s.re < 0.0 && pole.oscillation_within_source_band);
    assert!(prepared.substrate_memory_count() > 0 && prepared.coating_memory_count() > 1);
    assert_eq!(prepared.dimension(), 2 + memory(&image).len());
    let q = pole.displacement_m_sqrt_kg[0];
    let s = pole.mode.laplace_pole_per_s;
    near(pole.velocity_m_sqrt_kg_per_s[0], s * q, 1.0, 1e-12);
    for (index, &(row, tau)) in memory(&image).iter().enumerate() {
        let z = pole.mode.state[2 + index];
        near(z, q.scale(row[0]) / (C64::ONE + s.scale(tau)), 1.0, 1e-12);
    }
    let total_heat = pole
        .substrate_branch_dissipation_w
        .iter()
        .chain(&pole.coating_branch_dissipation_w)
        .sum::<f64>()
        + pole.background_dissipation_w
        + pole.interface_dissipation_w;
    assert!(pole.substrate_branch_dissipation_w.iter().sum::<f64>() > 0.0);
    assert!(pole.coating_branch_dissipation_w.iter().sum::<f64>() > 0.0);
    assert!((total_heat - pole.mode.dissipated_power_w).abs() < 1e-11 * total_heat);
    assert!(pole.decay_power_residual_w.abs() < 1e-10 * total_heat);
    near(
        pole.exterior_volume_velocity_m3_s(),
        -pole.cavity_volume_velocity_m3_s(),
        pole.exterior_volume_velocity_m3_s().abs(),
        1e-12,
    );
    for (i, row) in image.force_ports().iter().enumerate() {
        near(
            pole.force_port_displacement_m[i],
            q.scale(row[0]),
            q.abs() * row[0].abs(),
            1e-12,
        );
        near(
            pole.force_port_velocity_m_s()[i],
            s * pole.force_port_displacement_m[i],
            pole.force_port_velocity_m_s()[i].abs(),
            1e-11,
        );
    }
    assert!(pole.mode.relative_residual < 1e-11);
    assert!(pole.physical_relative_residual < 1e-11);
}

#[test]
fn g1_real_material_relaxation_is_retained_without_inventing_an_oscillation() {
    let gate = CancelGate::new_clock_free();
    let image = image(&gate, 1, true, true);
    let prepared = image.prepare_coating_resonance(budget(), &gate).unwrap();
    let omega = prepared
        .equilibrium_angular_frequencies_rad_s()
        .next()
        .unwrap();
    let branches = memory(&image);
    let slow_rate = branches
        .iter()
        .map(|(_, tau)| 1.0 / tau)
        .fold(f64::INFINITY, f64::min);
    let c = image.damping().matrix_per_s()[0];
    let (mut low, mut high) = (-slow_rate * (1.0 - 1e-12), 0.0);
    assert!(dispersion(C64::from_re(low), omega, c, &branches).0.re < 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (low + high);
        if dispersion(C64::from_re(mid), omega, c, &branches).0.re < 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    let expected = 0.5 * (low + high);
    let start = C64::from_re(expected * (1.0 - 1e-6));
    let seed = prepared
        .displacement_seed(start, &[C64::ONE], &gate)
        .unwrap();
    // The slow internal root has a much smaller mechanical participation than
    // the oscillatory mode. Refine its full state more tightly so the same
    // original componentwise physical gate can be met.
    let refinement = PoleOptions {
        residual_tolerance: 1e-15,
        step_tolerance: 1e-15,
        ..options(slow_rate)
    };
    let pole = prepared
        .resonance_pole(start, &seed, &refinement, &gate)
        .unwrap();
    near(
        pole.mode.laplace_pole_per_s,
        C64::from_re(expected),
        slow_rate,
        2e-10,
    );
    assert_eq!(pole.mode.laplace_pole_per_s.im, 0.0);
    assert!(!pole.oscillation_within_source_band);
    assert_eq!(pole.mode.state.len(), prepared.dimension());
    assert!(pole.mode.stored_energy_j > 0.0 && pole.mode.dissipated_power_w > 0.0);
    assert!(pole.decay_power_residual_w.abs() < 1e-10 * pole.mode.dissipated_power_w);
}

#[test]
fn g1_contour_discovers_the_same_full_material_mode_on_a_quadrature_ladder() {
    let gate = CancelGate::new_clock_free();
    let image = image(&gate, 1, true, true);
    let prepared = image.prepare_coating_resonance(budget(), &gate).unwrap();
    let omega = prepared
        .equilibrium_angular_frequencies_rad_s()
        .next()
        .unwrap();
    let expected = scalar_root(C64::new(-0.01 * omega, -omega), omega, &image);
    let region = ResonanceRegion {
        center: expected + C64::new(0.025 * omega, 0.03 * omega),
        real_radius: 0.15 * omega,
        imaginary_radius: 0.18 * omega,
    };
    let mut errors = Vec::new();
    // The nested eight-node count is unresolved at sixteen fine nodes. The
    // owner must refuse rather than publish a missing or guessed mode count.
    assert!(matches!(
        prepared.resonance_spectrum(region, &contour_options(omega, 16), &gate),
        Err(ShellResonanceError::Spectrum(
            fs_phs::harmonic::HarmonicSpectrumError::Contour(
                fs_la::nonlinear_eig::beyn::BeynError::CountUnresolved { .. }
            )
        ))
    ));
    for nodes in [32, 64, 128] {
        let report = prepared
            .resonance_spectrum(region, &contour_options(omega, nodes), &gate)
            .unwrap();
        assert_eq!(report.argument_count, 1);
        assert_eq!(report.resolved_rank, 1);
        assert_eq!(report.poles.len(), 1);
        assert_eq!(report.quadrature_nodes, nodes);
        assert_eq!(report.poles[0].mode.state.len(), prepared.dimension());
        near(
            report.poles[0].mode.laplace_pole_per_s,
            expected,
            omega,
            2e-10,
        );
        errors.push((report.unrefined_rates_per_s[0] - expected).abs() / omega);
    }
    assert!(errors[2] < 1e-10, "{errors:?}");
    assert!(
        errors[0] < 1e-10 || errors[2] < errors[0] * 0.1,
        "{errors:?}"
    );
}

#[test]
fn g3_full_cross_mode_memory_and_single_component_assembly_have_the_same_poles() {
    let gate = CancelGate::new_clock_free();
    let direct = image(&gate, 2, true, true);
    let prepared = direct.prepare_coating_resonance(budget(), &gate).unwrap();
    let pole = first_mode(&prepared, &gate);
    let frequencies: Vec<_> = prepared.equilibrium_angular_frequencies_rad_s().collect();
    let s = pole.mode.laplace_pole_per_s;
    let q = &pole.displacement_m_sqrt_kg;
    let branches = memory(&direct);
    assert!(branches.iter().any(|(row, _)| row[0] * row[1] != 0.0));
    for i in 0..2 {
        let mut residual = (s * s + C64::from_re(frequencies[i] * frequencies[i])) * q[i];
        let mut scale = (s.abs() * s.abs() + frequencies[i] * frequencies[i]) * q[i].abs();
        for j in 0..2 {
            let term = s * q[j].scale(direct.damping().matrix_per_s()[2 * i + j]);
            residual = residual + term;
            scale += term.abs();
        }
        for (arm, &(row, _)) in branches.iter().enumerate() {
            let strain = q[0].scale(row[0]) + q[1].scale(row[1]) - pole.mode.state[4 + arm];
            let term = strain.scale(row[i]);
            residual = residual + term;
            scale += term.abs();
        }
        assert!(residual.abs() < 1e-11 * scale, "{residual:?} scale={scale}");
    }
    let (shell, directors, constraints) = source(&gate, true, true);
    let assembly_ports: Vec<_> = ports(&shell)
        .into_iter()
        .map(|row| ShellAssemblyForcePort {
            components: vec![row],
        })
        .collect();
    let assembly = MaterialShellAssembly {
        components: vec![MaterialShellComponent {
            name: "coupon".into(),
            assembly_budget: fixture::budget(&shell),
            specimen: shell,
            directors,
            constraints,
            mass: ShellMass::Consistent,
            drilling_alpha: 1e-3,
        }],
        edge_ties: vec![],
        face_ties: vec![],
        tie_budget: TieBudget {
            max_segments: 32,
            max_slave_dofs: 24,
            gap_tolerance_m: 1e-10,
            angle_tolerance_rad: 1e-10,
        },
    };
    let assembled = assembly
        .reduce(
            &ShellAssemblyReductionOptions {
                eigenvalue_window: (0.0, 1e9),
                n_modes: 2,
                projection_budget: projection_budget(),
                max_total_nodes: 16,
                max_total_triangles: 18,
            },
            &assembly_ports,
            &gate,
        )
        .unwrap();
    let assembled_prepared = assembled
        .prepare_coating_resonance(budget(), &gate)
        .unwrap();
    let assembled_pole = first_mode(&assembled_prepared, &gate);
    near(assembled_pole.mode.laplace_pole_per_s, s, s.abs(), 1e-8);
    let first = pole.force_port_displacement_m[0];
    let other = assembled_pole.force_port_displacement_m[0];
    near(
        pole.force_port_displacement_m[1] / first,
        assembled_pole.force_port_displacement_m[1] / other,
        1.0,
        1e-8,
    );
    near(
        pole.exterior_volume_velocity_m3_s() / first,
        assembled_pole.exterior_volume_velocity_m3_s() / other,
        (pole.exterior_volume_velocity_m3_s() / first).abs(),
        1e-8,
    );
    assert_eq!(prepared.dimension(), assembled_prepared.dimension());
}

#[test]
fn g0_g4_frozen_sources_complete_state_resources_cancellation_and_replay() {
    let gate = CancelGate::new_clock_free();
    let image = image(&gate, 1, true, true);
    let prepared = image.prepare_coating_resonance(budget(), &gate).unwrap();
    let before = first_mode(&prepared, &gate);
    for bad in [
        ShellResonanceBudget {
            maximum_coordinates: prepared.dimension() - 1,
            ..budget()
        },
        ShellResonanceBudget {
            maximum_payload_bytes: 0,
            ..budget()
        },
        ShellResonanceBudget {
            maximum_preparation_work: 0,
            ..budget()
        },
    ] {
        assert!(matches!(
            image.prepare_coating_resonance(bad, &gate),
            Err(ShellResonanceError::Budget { .. })
        ));
    }
    let omega = prepared
        .equilibrium_angular_frequencies_rad_s()
        .next()
        .unwrap();
    let start = C64::new(-0.01 * omega, -omega);
    assert!(prepared.displacement_seed(start, &[], &gate).is_err());
    let tau = image.coating_relaxation().unwrap().relaxation_times_s()[0];
    assert!(
        prepared
            .displacement_seed(C64::from_re(-1.0 / tau), &[C64::ONE], &gate)
            .is_err()
    );
    let seed = prepared
        .displacement_seed(start, &[C64::ONE], &gate)
        .unwrap();
    assert!(
        prepared
            .resonance_pole(start, &seed[..seed.len() - 1], &options(omega), &gate)
            .is_err()
    );
    assert!(
        prepared
            .resonance_pole(
                start,
                &seed,
                &PoleOptions {
                    max_dimension: prepared.dimension() - 1,
                    ..options(omega)
                },
                &gate
            )
            .is_err()
    );
    assert!(
        prepared
            .resonance_pole(
                start,
                &seed,
                &PoleOptions {
                    max_dense_bytes: 0,
                    ..options(omega)
                },
                &gate
            )
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        image
            .prepare_coating_resonance(budget(), &cancelled)
            .is_err()
    );
    assert!(
        prepared
            .resonance_pole(start, &seed, &options(omega), &cancelled)
            .is_err()
    );
    assert_eq!(before, first_mode(&prepared, &gate));
    for (film, substrate) in [(false, true), (true, false)] {
        let frozen = self::image(&gate, 1, film, substrate);
        assert!(matches!(
            frozen.prepare_coating_resonance(budget(), &gate),
            Err(ShellResonanceError::Invalid(_))
        ));
    }
}
