//! One full material shell, original P2 interior and independently observed exterior.

use fs_bem::helmholtz::Medium as ExteriorMedium;
use fs_bem::near_field::Options;
use fs_bem::reciprocal_radiation::{Budget, CoherentPortFieldLimits};
use fs_couple::cavity_modal::interior::structure::shell::exterior::TensorShellInteriorExteriorControls;
use fs_couple::cavity_modal::interior::structure::shell::{
    ConformingShellFace, TensorShellInterior, TensorShellInteriorControls, TensorShellInteriorError,
};
use fs_couple::tensor_shell_radiation::{
    TensorShellRadiation, TensorShellRadiationAcousticControls, TensorShellRadiationControls,
    TensorShellRadiationFrequencyControls,
};
use fs_exec::{CancelGate, Cx};
use fs_feec::acoustics::helmholtz::heterogeneous::HeterogeneousHelmholtzFe;
use fs_feec::acoustics::helmholtz::heterogeneous::time::ResistiveBoundary;
use fs_feec::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;
use fs_modal::harmonic::{HarmonicPortControls, HarmonicResponseControls};
use fs_plate::shell::viscosity::ShellViscosityBudget;
use fs_plate::shell::viscosity::harmonic::{
    TensorPronyShell, TensorPronyShellSection, TensorShellHarmonicError,
};
use fs_plate::shell::{CstDktMass, ShellMesh, ShellSupport};
use fs_rep_mesh::TetComplex;

#[path = "tensor_shell_interior_exterior/sources.rs"]
mod sources;

const POINTS: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [0.4, 0.0, 0.0],
    [0.0, 0.3, 0.0],
    [0.0, 0.0, 0.5],
];
const FACETS: [[usize; 3]; 4] = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
const OMEGA: f64 = 100.0;
const WALL_Y: f64 = 2e-5;

fn interfaces() -> [ConformingShellFace; 4] {
    FACETS.map(|shell| ConformingShellFace {
        fluid: shell.map(|v| v as u32),
        shell,
    })
}
fn source(cx: &Cx<'_>) -> TensorPronyShell {
    // Explicit synthetic constitutive and fluid data, not measured validation.
    let c: Mandel6 = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            1e6 * if i < 3 && j < 3 {
                4.0 + if i == j { 6.0 } else { 0.0 }
            } else if i == j {
                6.0
            } else {
                0.0
            }
        })
    });
    let scaled = |s| c.map(|row| row.map(|v| v * s));
    let law = TensorProny::new_mandel(
        c,
        vec![
            PronyBranch::new(scaled(0.4), 0.001).unwrap(),
            PronyBranch::new(scaled(0.2), 0.02).unwrap(),
        ],
        0.02,
    )
    .unwrap()
    .with_kelvin_viscosity_mandel(scaled(1e-5))
    .unwrap()
    .with_band((0.0, 150.0))
    .unwrap();
    let section = TensorPronyShellSection::new(law, 0.01, 1200.0).unwrap();
    TensorPronyShell::assemble(
        ShellMesh::new(POINTS.to_vec(), FACETS.to_vec()).unwrap(),
        vec![section; 4],
        &[0, 1, 2],
        ShellSupport::Clamped,
        CstDktMass::Consistent,
        0.0,
        ShellViscosityBudget {
            max_elements: 4,
            max_full_dofs: 24,
            max_dofs: 6,
            max_modes: 6,
            max_work: 20_000_000,
        },
        cx,
    )
    .unwrap()
}
fn harmonic(coordinates: usize) -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: coordinates,
        maximum_payload_bytes: 16 << 20,
        maximum_relative_residual: 2e-11,
    }
}
fn interior_controls() -> TensorShellInteriorControls {
    TensorShellInteriorControls {
        material: harmonic(6),
        coupled: harmonic(16),
        maximum_interface_faces: 4,
        maximum_trace_bytes: 64 << 10,
        maximum_trace_work: 1_000_000,
        power_relative_tolerance: 2e-8,
    }
}
fn boundary_budget() -> Budget {
    Budget {
        maximum_unknowns: 128,
        maximum_cache_bytes: 16 << 20,
        maximum_work: 2_000_000_000,
        angular_refinement: 0,
    }
}
fn exterior_controls() -> TensorShellInteriorExteriorControls {
    TensorShellInteriorExteriorControls {
        acoustic: TensorShellRadiationAcousticControls {
            maximum_boundary_payload_bytes: 16 << 20,
            field: CoherentPortFieldLimits {
                maximum_payload_bytes: 1 << 20,
                maximum_work: 10_000_000,
            },
            maximum_radiated_power_work: 100_000_000,
        },
        ports: HarmonicPortControls {
            maximum_ports: 6,
            maximum_projection_terms: 1_000_000,
            maximum_payload_bytes: 4 << 20,
        },
        maximum_force_map_bytes: 8 * 16 * 6,
        power_relative_tolerance: 2e-8,
    }
}
fn medium() -> ExteriorMedium {
    ExteriorMedium {
        density: 1000.0,
        sound_speed: 1480.0,
    }
}
fn wall() -> [ResistiveBoundary; 1] {
    [ResistiveBoundary {
        faces: vec![[0, 2, 1]],
        admittance_m_pa_s: WALL_Y,
    }]
}
fn fixture(
    test: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &TensorPronyShell,
        &HeterogeneousHelmholtzFe<'_, '_>,
        &TensorShellRadiation<'_>,
    ),
) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 89,
                kernel_id: 31,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let shell = source(&cx);
        let complex = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
        let gas = Medium {
            density_kg_m3: 12.0,
            sound_speed_m_s: 200.0,
        };
        let geometry = HelmholtzFe::new(&complex, &POINTS, 2, gas).unwrap();
        let interior = HeterogeneousHelmholtzFe::new(
            &geometry,
            5,
            &mut |_, _| Ok::<_, AcousticFeError>(gas),
            &mut || false,
        )
        .unwrap();
        let exterior = TensorShellRadiation::new(
            &shell,
            medium(),
            TensorShellRadiationControls {
                acoustic_refinements: 1,
                maximum_panels: 16,
                maximum_coordinates: 6,
                maximum_setup_bytes: 1 << 20,
                maximum_setup_work: 10_000_000,
            },
            &cx,
        )
        .unwrap();
        test(&cx, &gate, &shell, &interior, &exterior);
    });
}
fn force() -> [C64; 6] {
    [
        C64::new(1e-4, -2e-5),
        C64::new(-3e-5, 1e-5),
        C64::new(2e-5, 4e-5),
        C64::new(1e-6, -2e-7),
        C64::new(2e-7, 1e-7),
        C64::new(-5e-7, 1e-7),
    ]
}
fn near(a: C64, b: C64) {
    assert!((a - b).abs() <= 3e-7 * b.abs().max(1e-22), "{a:?} != {b:?}");
}
fn work(force: &[C64], velocity: &[C64]) -> f64 {
    0.5 * force
        .iter()
        .zip(velocity)
        .map(|(&f, &v)| (f * v.conj()).re)
        .sum::<f64>()
}
fn product(matrix: &fs_sparse::Csr, x: &[C64], row: usize) -> C64 {
    let (columns, values) = matrix.row(row);
    columns
        .iter()
        .zip(values)
        .fold(C64::ZERO, |sum, (&j, &a)| sum + x[j].scale(a))
}

// Cubic-exact integral of original P2 pressure times the actual apex P1 shape.
// Original oriented area vectors are computed independently of adapter B/normal rows.
fn interior_force(interior: &HeterogeneousHelmholtzFe<'_, '_>, pressure: &[C64]) -> [C64; 6] {
    let rule = [
        ([1.0 / 3.0; 3], -27.0 / 48.0),
        ([0.6, 0.2, 0.2], 25.0 / 48.0),
        ([0.2, 0.6, 0.2], 25.0 / 48.0),
        ([0.2, 0.2, 0.6], 25.0 / 48.0),
    ];
    let mut force = [C64::ZERO; 6];
    for face in FACETS {
        let a: [f64; 3] = core::array::from_fn(|d| POINTS[face[1]][d] - POINTS[face[0]][d]);
        let b: [f64; 3] = core::array::from_fn(|d| POINTS[face[2]][d] - POINTS[face[0]][d]);
        let area = [
            0.5 * (a[1] * b[2] - a[2] * b[1]),
            0.5 * (a[2] * b[0] - a[0] * b[2]),
            0.5 * (a[0] * b[1] - a[1] * b[0]),
        ];
        for (l, weight) in rule {
            let mut barycentric = [0.0; 4];
            for j in 0..3 {
                barycentric[face[j]] = l[j];
            }
            let p = interior.pressure_at(pressure, 0, barycentric).unwrap();
            for d in 0..3 {
                force[d] = force[d] + p.scale(weight * barycentric[3] * area[d]);
            }
        }
    }
    force
}

#[test]
fn one_full_factor_balances_original_shell_interior_and_independently_observed_exterior() {
    fixture(|cx, _, shell, interior, geometry| {
        // Paired caller winding changes neither physical side nor work orientation.
        let reversed = interfaces().map(|mut face| {
            face.fluid.swap(0, 1);
            face.shell.swap(0, 1);
            face
        });
        let owner = TensorShellInterior::new(
            interior,
            shell,
            &reversed,
            &wall(),
            OMEGA,
            interior_controls(),
            cx,
        )
        .unwrap();
        let boundary = geometry.prepare(120.0, boundary_budget(), cx).unwrap();
        let prepared = owner
            .prepare_with_exterior(&boundary, exterior_controls(), cx)
            .unwrap();
        let flow = owner
            .project_inward_velocity(
                &[(
                    [0, 2, 1],
                    [
                        C64::new(1e-8, -2e-9),
                        C64::new(-3e-9, 4e-9),
                        C64::new(5e-9, 1e-9),
                    ],
                )],
                cx,
            )
            .unwrap();
        let response = prepared.solve(&force(), &flow, cx).unwrap();
        let field = &response.coupled.field;
        assert_eq!(field.structural_displacement.len(), 6);
        assert_eq!(field.pressure_pa.len(), 10); // All original P2 coordinates, not a cavity mode.
        let internal = interior_force(interior, &field.pressure_pa);
        for j in 0..6 {
            near(field.fluid_force[j], internal[j]);
        }

        // Independently integrate exterior traction from the returned physical
        // surface pressure and analytic apex trace z/0.5; do not use adapter rows.
        let mut external = [C64::ZERO; 6];
        let mut velocity = Vec::new();
        for ((&point, &normal), (&area, &pressure)) in geometry
            .surface()
            .centroids()
            .iter()
            .zip(geometry.surface().normals())
            .zip(
                geometry
                    .surface()
                    .areas()
                    .iter()
                    .zip(response.exterior.surface_pressure_pa()),
            )
        {
            let shape = point[2] / 0.5;
            let mut v = C64::ZERO;
            for d in 0..3 {
                external[d] = external[d] + pressure.scale(shape * normal[d] * area);
                v = v + field.structural_velocity[d].scale(shape * normal[d]);
            }
            velocity.push(v);
        }
        let load = field.structural_port_load.as_ref().unwrap();
        for j in 0..6 {
            near(load.flow[j], field.structural_velocity[j]);
            near(load.effort[j], external[j]);
            let material = owner.material_operator();
            let k = product(material.stiffness(), &field.structural_displacement, j);
            let m =
                product(material.mass(), &field.structural_displacement, j).scale(OMEGA * OMEGA);
            let c = product(material.damping(), &field.structural_velocity, j);
            let residual = k - m + c - internal[j] + external[j] - force()[j];
            assert!(
                residual.abs()
                    < 2e-8
                        * (k.abs()
                            + m.abs()
                            + c.abs()
                            + internal[j].abs()
                            + external[j].abs()
                            + force()[j].abs())
            );
        }
        let net: Vec<_> = (0..6)
            .map(|j| force()[j] + internal[j] - external[j])
            .collect();
        let native = shell
            .prepare_harmonic_response(OMEGA, harmonic(6), cx)
            .unwrap()
            .solve(&net, cx)
            .unwrap();
        for j in 0..6 {
            near(field.structural_displacement[j], native.displacement[j]);
        }
        let surface = interior.geometry().boundary_mass(&[[0, 2, 1]]).unwrap();
        for row in 0..interior.ndof() {
            let residual = product(interior.stiffness(), &field.pressure_potential_pa_s, row)
                - product(interior.mass(), &field.pressure_potential_pa_s, row)
                    .scale(OMEGA * OMEGA)
                + product(&surface, &field.pressure_pa, row).scale(WALL_Y)
                + field.displaced_flow_m3_s[row]
                - flow[row];
            assert!(residual.abs() < 1e-15);
        }
        let direct = boundary
            .boundary()
            .solve_batch(OMEGA / medium().sound_speed, medium(), &[&velocity], || {
                false
            })
            .unwrap();
        for (&a, &b) in response
            .exterior
            .surface_pressure_pa()
            .iter()
            .zip(&direct[0].pressure)
        {
            near(a, b);
        }
        let receivers = prepared
            .prepare_receivers(
                &[[0.1, -0.3, 0.8], [-0.3, 0.1, 0.2]],
                0.01,
                Options::default(),
                true,
                cx,
            )
            .unwrap();
        let measured = receivers
            .evaluate_coherent_velocity(&response.exterior)
            .unwrap();
        let reference = receivers.evaluate_velocity(&direct[0]).unwrap();
        for r in 0..2 {
            near(measured.scalar.pressure[r], reference.scalar.pressure[r]);
            for d in 0..3 {
                near(
                    measured.particle_velocity_m_s[r][d],
                    reference.particle_velocity_m_s[r][d],
                );
            }
        }
        let scale = response.coupled.apparent_power_scale_w;
        assert!(
            (work(&internal, &field.structural_velocity)
                - field.average_fluid_to_structure_power_w)
                .abs()
                < 2e-8 * scale
        );
        assert!(
            (work(&external, &field.structural_velocity) - response.radiation.outgoing_power_w)
                .abs()
                < 2e-8 * scale
        );
        assert!(field.average_structural_dissipation_w > 0.0);
        assert!(field.average_wall_dissipation_w > 0.0);
        assert!(response.radiation.outgoing_power_w > 0.0);
        assert!(field.structural_power_balance_residual_w().abs() < 2e-8 * scale);
        assert!(field.acoustic_power_balance_residual_w().abs() < 2e-8 * scale);
        assert!(response.power_balance_residual_w.abs() < 2e-8 * scale);
    });
}

#[test]
fn both_fluids_change_motion_detachment_recovers_exterior_and_coherent_sources_superpose() {
    fixture(|cx, _, shell, interior, geometry| {
        let owner = TensorShellInterior::new(
            interior,
            shell,
            &interfaces(),
            &wall(),
            OMEGA,
            interior_controls(),
            cx,
        )
        .unwrap();
        let boundary = geometry.prepare(120.0, boundary_budget(), cx).unwrap();
        let prepared = owner
            .prepare_with_exterior(&boundary, exterior_controls(), cx)
            .unwrap();
        let zero = vec![C64::ZERO; interior.ndof()];
        let both = prepared.solve(&force(), &zero, cx).unwrap();
        let interior_only = owner
            .prepare(cx)
            .unwrap()
            .solve(&force(), &zero, cx)
            .unwrap();
        let detached = TensorShellInterior::new(
            interior,
            shell,
            &[],
            &wall(),
            OMEGA,
            interior_controls(),
            cx,
        )
        .unwrap();
        let exterior_only = detached
            .prepare_with_exterior(&boundary, exterior_controls(), cx)
            .unwrap()
            .solve(&force(), &zero, cx)
            .unwrap();
        let limits = exterior_controls();
        let native = boundary
            .prepare_frequency(
                OMEGA,
                TensorShellRadiationFrequencyControls {
                    mechanical: harmonic(6),
                    ports: limits.ports,
                    maximum_boundary_payload_bytes: limits.acoustic.maximum_boundary_payload_bytes,
                    field: limits.acoustic.field,
                    maximum_radiated_power_work: limits.acoustic.maximum_radiated_power_work,
                    power_relative_tolerance: limits.power_relative_tolerance,
                },
                cx,
            )
            .unwrap()
            .solve(&force(), cx)
            .unwrap();
        let scale = both
            .coupled
            .field
            .structural_displacement
            .iter()
            .take(3)
            .map(|q| q.abs())
            .fold(0.0, f64::max);
        for j in 0..6 {
            near(
                exterior_only.coupled.field.structural_displacement[j],
                native.motion().displacement[j],
            );
        }
        assert!(
            exterior_only
                .coupled
                .field
                .pressure_pa
                .iter()
                .all(|p| p.abs() < 1e-20)
        );
        for reference in [&interior_only.field, &exterior_only.coupled.field] {
            assert!(
                reference
                    .structural_displacement
                    .iter()
                    .take(3)
                    .zip(&both.coupled.field.structural_displacement)
                    .any(|(&a, &b)| (a - b).abs() > 1e-4 * scale)
            );
        }
        let flow = owner
            .project_inward_velocity(&[([0, 2, 1], [C64::new(1e-8, 4e-9); 3])], cx)
            .unwrap();
        let acoustic = prepared.solve(&[C64::ZERO; 6], &flow, cx).unwrap();
        let phase = C64::new(-0.3, 0.8);
        let mixed_flow: Vec<_> = flow.iter().map(|&v| phase * v).collect();
        let mixed = prepared.solve(&force(), &mixed_flow, cx).unwrap();
        for j in 0..6 {
            near(
                mixed.coupled.field.structural_displacement[j],
                both.coupled.field.structural_displacement[j]
                    + phase * acoustic.coupled.field.structural_displacement[j],
            );
        }
        for j in 0..interior.ndof() {
            near(
                mixed.coupled.field.pressure_pa[j],
                both.coupled.field.pressure_pa[j] + phase * acoustic.coupled.field.pressure_pa[j],
            );
        }
        for ((&a, &b), &sum) in both
            .exterior
            .surface_pressure_pa()
            .iter()
            .zip(acoustic.exterior.surface_pressure_pa())
            .zip(mixed.exterior.surface_pressure_pa())
        {
            near(sum, a + phase * b);
        }
        let silence = prepared.solve(&[C64::ZERO; 6], &zero, cx).unwrap();
        assert_eq!(silence.radiation.outgoing_power_w, 0.0);
        assert_eq!(silence.power_balance_residual_w, 0.0);
    });
}

#[test]
fn wrong_source_wrong_fluid_side_bands_budgets_strain_and_cancellation_refuse() {
    fixture(|cx, gate, shell, interior, geometry| {
        let owner = TensorShellInterior::new(
            interior,
            shell,
            &interfaces(),
            &[],
            OMEGA,
            interior_controls(),
            cx,
        )
        .unwrap();
        let boundary = geometry.prepare(120.0, boundary_budget(), cx).unwrap();
        let other_shell = source(cx);
        let other = TensorShellInterior::new(
            interior,
            &other_shell,
            &interfaces(),
            &[],
            OMEGA,
            interior_controls(),
            cx,
        )
        .unwrap();
        assert!(matches!(
            other.prepare_with_exterior(&boundary, exterior_controls(), cx),
            Err(TensorShellInteriorError::Invalid(_))
        ));
        let mut outside = POINTS;
        outside[3][2] = -0.5;
        let complex = TetComplex::from_tets(4, vec![[0, 2, 1, 3]]);
        let gas = Medium {
            density_kg_m3: 12.0,
            sound_speed_m_s: 200.0,
        };
        let wrong_geometry = HelmholtzFe::new(&complex, &outside, 2, gas).unwrap();
        let wrong_fluid = HeterogeneousHelmholtzFe::new(
            &wrong_geometry,
            5,
            &mut |_, _| Ok::<_, AcousticFeError>(gas),
            &mut || false,
        )
        .unwrap();
        let wrong_side = TensorShellInterior::new(
            &wrong_fluid,
            shell,
            &interfaces()[..1],
            &[],
            OMEGA,
            interior_controls(),
            cx,
        )
        .unwrap();
        assert!(matches!(
            wrong_side.prepare_with_exterior(&boundary, exterior_controls(), cx),
            Err(TensorShellInteriorError::Invalid(_))
        ));
        let short_band = geometry.prepare(80.0, boundary_budget(), cx).unwrap();
        assert!(
            owner
                .prepare_with_exterior(&short_band, exterior_controls(), cx)
                .is_err()
        );
        for limits in [
            TensorShellInteriorExteriorControls {
                maximum_force_map_bytes: 8 * 16 * 6 - 1,
                ..exterior_controls()
            },
            TensorShellInteriorExteriorControls {
                ports: HarmonicPortControls {
                    maximum_ports: 5,
                    ..exterior_controls().ports
                },
                ..exterior_controls()
            },
            TensorShellInteriorExteriorControls {
                acoustic: TensorShellRadiationAcousticControls {
                    maximum_boundary_payload_bytes: 1,
                    ..exterior_controls().acoustic
                },
                ..exterior_controls()
            },
            TensorShellInteriorExteriorControls {
                power_relative_tolerance: f64::NAN,
                ..exterior_controls()
            },
        ] {
            assert!(owner.prepare_with_exterior(&boundary, limits, cx).is_err());
        }
        let prepared = owner
            .prepare_with_exterior(&boundary, exterior_controls(), cx)
            .unwrap();
        let zero = vec![C64::ZERO; interior.ndof()];
        assert!(prepared.solve(&[C64::ZERO; 5], &zero, cx).is_err());
        assert!(
            prepared
                .solve(&[C64::new(f64::NAN, 0.0); 6], &zero, cx)
                .is_err()
        );
        assert!(matches!(
            prepared.solve(&[C64::from_re(1e9); 6], &zero, cx),
            Err(TensorShellInteriorError::Shell(
                TensorShellHarmonicError::Material { .. }
            ))
        ));
        gate.request();
        assert!(prepared.solve(&force(), &zero, cx).is_err());
        assert!(
            owner
                .prepare_with_exterior(&boundary, exterior_controls(), cx)
                .is_err()
        );
    });
}
