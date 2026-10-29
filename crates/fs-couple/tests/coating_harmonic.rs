//! Source-backed coupon response checks; not empirical violin validation.
#![deny(unsafe_code)]

#[path = "coating_shell.rs"]
mod fixture;

use fs_couple::acoustic_realize::AcousticRealizeError;
use fs_couple::shear_shell::{
    CoatingSide, MaterialShellImage, ResolvedShearShell, ShellProjectionBudget,
    ShellReductionOptions,
    assembly::{
        MaterialShellAssembly, MaterialShellComponent, ShellAssemblyForcePort,
        ShellAssemblyReductionOptions,
    },
    coating::{harmonic::ShellCoatingHarmonicOptions, relaxation::CoatingRelaxationOptions},
    pressure::ShellPressureSurfaces,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_plate::shell::{
    junction::TieBudget,
    mitc::{NodalDirectors, ShellConstraint, ShellMass},
    pressure::{OrientedShellFacet, ShellFace},
};
use std::f64::consts::TAU;

fn source(
    gate: &CancelGate,
    memory: bool,
) -> (ResolvedShearShell, NodalDirectors, Vec<ShellConstraint>) {
    let (shell, directors, constraints) = fixture::coupon(gate);
    let shell = shell
        .with_coatings(vec![fixture::layer(CoatingSide::Positive, 20e-6)], gate)
        .unwrap();
    let shell = if memory {
        shell
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
    } else {
        shell
    };
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

fn physical_ports(shell: &ResolvedShearShell) -> Vec<Vec<f64>> {
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

fn make_image(gate: &CancelGate, memory: bool) -> MaterialShellImage {
    reduce_image(source(gate, memory), gate)
}

fn reduce_image(
    (shell, directors, constraints): (ResolvedShearShell, NodalDirectors, Vec<ShellConstraint>),
    gate: &CancelGate,
) -> MaterialShellImage {
    shell
        .reduce(
            &directors,
            &ShellReductionOptions {
                constraints,
                mass: ShellMass::Consistent,
                drilling_alpha: 1e-3,
                eigenvalue_window: (0.0, 1e9),
                n_modes: 3,
                assembly_budget: fixture::budget(&shell),
                projection_budget: projection_budget(),
            },
            &physical_ports(&shell),
            gate,
        )
        .unwrap()
}

fn options() -> ShellCoatingHarmonicOptions {
    ShellCoatingHarmonicOptions {
        maximum_work: 1_000_000,
        solve_relative_tolerance: 1e-11,
    }
}

fn relative_complex(a: &[C64], b: &[C64]) -> f64 {
    let numerator = a
        .iter()
        .zip(b)
        .map(|(&a, &b)| (a - b).norm_sq())
        .sum::<f64>()
        .sqrt();
    numerator / b.iter().map(|z| z.norm_sq()).sum::<f64>().sqrt().max(1e-30)
}

#[test]
fn g0_source_film_response_preserves_power_port_reciprocity_and_pressure_signs() {
    let gate = CancelGate::new_clock_free();
    let image = make_image(&gate, true);
    let force = |port: usize| {
        image.force_ports()[port]
            .iter()
            .copied()
            .map(C64::from_re)
            .collect::<Vec<_>>()
    };
    for hz in [120.0, 200.0, 1000.0, 4000.0] {
        let a = image
            .coating_harmonic_response(TAU * hz, &force(0), options(), &gate)
            .unwrap();
        let b = image
            .coating_harmonic_response(TAU * hz, &force(1), options(), &gate)
            .unwrap();
        assert!(a.solve_relative_residual < 1e-12);
        assert!(a.external_power_w > 0.0 && a.background_dissipation_w > 0.0);
        assert_eq!(a.interface_dissipation_w, 0.0);
        assert!(a.film_dissipation_w[0] > 0.0);
        assert!(a.power_residual_w.abs() < 1e-10 * a.total_dissipation_w());
        assert!(
            relative_complex(
                &[a.force_port_velocity_m_s[1]],
                &[b.force_port_velocity_m_s[0]]
            ) < 1e-10
        );
        assert!(
            relative_complex(
                &[a.exterior_volume_velocity_m3_s],
                &[-a.cavity_volume_velocity_m3_s]
            ) < 1e-12
        );
        assert_eq!(a.external_power_w, 0.5 * a.force_port_velocity_m_s[0].re);
    }
    let zero = image
        .coating_harmonic_response(TAU * 1000.0, &[C64::ZERO; 3], options(), &gate)
        .unwrap();
    assert_eq!(zero.displacement_m_sqrt_kg, [C64::ZERO; 3]);
    assert_eq!(zero.total_dissipation_w(), 0.0);
    assert_eq!(zero.power_residual_w, 0.0);
    let bare = reduce_image(fixture::coupon(&gate), &gate);
    let bare_force: Vec<_> = bare.force_ports()[0].iter().copied().map(C64::from_re).collect();
    let baseline = bare.coating_harmonic_response(TAU*1000.0, &bare_force, options(), &gate).unwrap();
    let coated = image.coating_harmonic_response(TAU*1000.0, &force(0), options(), &gate).unwrap();
    assert!(baseline.film_dissipation_w.is_empty());
    assert!(baseline.external_power_w > 0.0);
    assert!(baseline.power_residual_w.abs() < 1e-10*baseline.total_dissipation_w());
    assert!(relative_complex(&coated.force_port_velocity_m_s, &baseline.force_port_velocity_m_s) > 1e-7);
}

#[test]
fn g3_exact_film_operator_matches_complete_projected_prony_memory_in_band() {
    let gate = CancelGate::new_clock_free();
    let image = make_image(&gate, true);
    let memory = image.coating_relaxation().unwrap();
    let n = memory.mode_count();
    assert_eq!(memory.films().len(), 1);
    let film = &memory.films()[0];
    let h = film.unit_stiffness_per_pa_s2();
    assert!(
        h.iter()
            .enumerate()
            .any(|(ij, value)| ij / n != ij % n && value.abs() > 1e-12)
    );
    for hz in [120.0, 200.0, 500.0, 1000.0, 2500.0, 4000.0] {
        let omega = TAU * hz;
        let (storage, loss) = film.source().input().material.law().modulus(omega).unwrap();
        let e0 = film.source().reference_storage_pa();
        let mut retained: Vec<_> = h.iter().map(|h| C64::from_re(e0 * h)).collect();
        for (row, &tau) in memory
            .coupling_rows_per_s()
            .iter()
            .zip(memory.relaxation_times_s())
        {
            let x = omega * tau;
            let response = C64::new(x * x / (1.0 + x * x), -x / (1.0 + x * x));
            for i in 0..n {
                for j in 0..n {
                    retained[i * n + j] = retained[i * n + j] + response.scale(row[i] * row[j]);
                }
            }
        }
        // Storage and loss are checked independently; E0 occurs once.
        let exact: Vec<_> = h.iter().map(|h| C64::new(storage * h, -loss * h)).collect();
        for component in [false, true] {
            let select = |z: &C64| if component { z.im } else { z.re };
            let numerator = retained
                .iter()
                .zip(&exact)
                .map(|(a, b)| (select(a) - select(b)).powi(2))
                .sum::<f64>()
                .sqrt();
            let denominator = exact.iter().map(|z| select(z).powi(2)).sum::<f64>().sqrt();
            assert!(
                numerator / denominator < 0.001,
                "component={component}, Hz={hz}"
            );
        }
        let response = image
            .coating_harmonic_response(
                omega,
                &[C64::ONE, C64::new(0.0, 0.35), C64::from_re(-0.22)],
                options(),
                &gate,
            )
            .unwrap();
        let q = &response.displacement_m_sqrt_kg;
        // Actual transient energy rows, with no second scalar fit.
        let mut memory_power = 0.0;
        for (row, &tau) in memory
            .coupling_rows_per_s()
            .iter()
            .zip(memory.relaxation_times_s())
        {
            let strain = row
                .iter()
                .zip(q)
                .fold(C64::ZERO, |sum, (&c, &q)| sum + q.scale(c));
            let x = omega * tau;
            memory_power += 0.5 * omega * x / (1.0 + x * x) * strain.norm_sq();
        }
        let exact_power = response.film_dissipation_w[0];
        assert!((memory_power - exact_power).abs() / exact_power < 0.001);
    }
}

#[test]
fn g4_harmonic_source_support_complete_forcing_budget_and_cancellation_refuse() {
    let gate = CancelGate::new_clock_free();
    let image = make_image(&gate, true);
    for omega in [0.0, -1.0, f64::NAN, TAU * 20.0, TAU * 20000.0] {
        assert!(
            image
                .coating_harmonic_response(omega, &[C64::ONE; 3], options(), &gate)
                .is_err()
        );
    }
    assert!(
        image
            .coating_harmonic_response(TAU * 1000.0, &[C64::ONE; 2], options(), &gate)
            .is_err()
    );
    let mut insufficient = options();
    insufficient.maximum_work = 0;
    assert!(
        image
            .coating_harmonic_response(TAU * 1000.0, &[C64::ONE; 3], insufficient, &gate)
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        image.coating_harmonic_response(TAU * 1000.0, &[C64::ONE; 3], options(), &cancelled),
        Err(AcousticRealizeError::Cancelled)
    ));
    let frozen = make_image(&gate, false);
    assert!(
        frozen
            .coating_harmonic_response(TAU * 1000.0, &[C64::ONE; 3], options(), &gate)
            .is_err()
    );
}

#[test]
fn g3_single_component_assembly_preserves_source_response_and_signed_volume() {
    let gate = CancelGate::new_clock_free();
    let direct = make_image(&gate, true);
    let (shell, directors, constraints) = source(&gate, true);
    let ports: Vec<_> = physical_ports(&shell)
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
                n_modes: 3,
                projection_budget: projection_budget(),
                max_total_nodes: 16,
                max_total_triangles: 18,
            },
            &ports,
            &gate,
        )
        .unwrap();
    let force = |rows: &[Vec<f64>]| {
        rows[0]
            .iter()
            .copied()
            .map(C64::from_re)
            .collect::<Vec<_>>()
    };
    let a = direct
        .coating_harmonic_response(TAU * 1000.0, &force(direct.force_ports()), options(), &gate)
        .unwrap();
    let b = assembled
        .coating_harmonic_response(
            TAU*1000.0, &force(assembled.force_ports()), options(), &gate).unwrap();
    assert!(relative_complex(&a.force_port_velocity_m_s, &b.force_port_velocity_m_s) < 1e-8);
    assert!(relative_complex(&[a.exterior_volume_velocity_m3_s], &[b.exterior_volume_velocity_m3_s]) < 1e-8);
    assert!(b.power_residual_w.abs() < 1e-10*b.total_dissipation_w());
}
 
