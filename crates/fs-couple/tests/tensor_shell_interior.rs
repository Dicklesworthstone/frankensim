//! Original tensor-shell motion, complete P2 cavity fields and reciprocal work.

use fs_couple::cavity_modal::interior::structure::shell::{
    ConformingShellFace, TensorShellInterior, TensorShellInteriorControls, TensorShellInteriorError,
};
use fs_exec::{CancelGate, Cx};
use fs_feec::acoustics::helmholtz::heterogeneous::HeterogeneousHelmholtzFe;
use fs_feec::acoustics::helmholtz::heterogeneous::time::ResistiveBoundary;
use fs_feec::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;
use fs_modal::harmonic::HarmonicResponseControls;
use fs_plate::shell::viscosity::ShellViscosityBudget;
use fs_plate::shell::viscosity::harmonic::{
    TensorPronyShell, TensorPronyShellSection, TensorShellHarmonicError,
};
use fs_plate::shell::{CstDktMass, ShellMesh, ShellSupport};
use fs_rep_mesh::TetComplex;

fn controls() -> TensorShellInteriorControls {
    let harmonic = HarmonicResponseControls {
        maximum_coordinates: 64,
        maximum_payload_bytes: 4 * 1024 * 1024,
        maximum_relative_residual: 1e-10,
    };
    TensorShellInteriorControls {
        material: harmonic,
        coupled: harmonic,
        maximum_interface_faces: 3,
        maximum_trace_bytes: 64 * 1024,
        maximum_trace_work: 1_000_000,
        power_relative_tolerance: 2e-9,
    }
}

fn faces() -> [ConformingShellFace; 3] {
    [[0, 1], [1, 2], [2, 0]].map(|[a, b]| ConformingShellFace {
        fluid: [a, b, 4],
        shell: [a as usize, b as usize, 3],
    })
}

fn fixture(
    test: impl FnOnce(&Cx<'_>, &CancelGate, &TensorPronyShell, &HeterogeneousHelmholtzFe<'_, '_>),
) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 73,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        // Split the tetrahedron around its base-face center. Each original
        // boundary triangle has area 1/6 and an independently known normal -z.
        let points = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0 / 3.0, 1.0 / 3.0, 0.0],
        ];
        let complex = TetComplex::from_tets(5, vec![[0, 1, 4, 3], [1, 2, 4, 3], [2, 0, 4, 3]]);
        let medium = Medium {
            density_kg_m3: 1.2,
            sound_speed_m_s: 3.0,
        };
        let geometry = HelmholtzFe::new(&complex, &points, 2, medium).unwrap();
        let interior = HeterogeneousHelmholtzFe::new(
            &geometry,
            5,
            &mut |_, _| Ok::<_, AcousticFeError>(medium),
            &mut || false,
        )
        .unwrap();
        let c: Mandel6 = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                if i == j {
                    [12.0, 18.0, 24.0, 6.0, 9.0, 7.0][i]
                } else if i < 3 && j < 3 {
                    3.0
                } else {
                    0.0
                }
            })
        });
        let mut branch = c.map(|row| row.map(|v| 0.6 * v));
        branch[0][0] *= 1.7; // An anisotropic, nonproportional relaxation branch.
        let law = TensorProny::new_mandel(c, vec![PronyBranch::new(branch, 0.3).unwrap()], 0.02)
            .unwrap()
            .with_kelvin_viscosity_mandel(c.map(|row| row.map(|v| 0.02 * v)))
            .unwrap()
            .with_band((0.0, 40.0))
            .unwrap();
        let section = TensorPronyShellSection::new(law, 0.08, 4.0).unwrap();
        let shell = TensorPronyShell::assemble(
            ShellMesh::new(
                vec![points[0], points[1], points[2], points[4]],
                vec![[0, 1, 3], [1, 2, 3], [2, 0, 3]],
            )
            .unwrap(),
            vec![section; 3],
            &[0, 1, 2],
            ShellSupport::Clamped,
            CstDktMass::Consistent,
            0.0,
            ShellViscosityBudget {
                max_elements: 3,
                max_full_dofs: 24,
                max_dofs: 6,
                max_modes: 6,
                max_work: 20_000_000,
            },
            &cx,
        )
        .unwrap();
        assert_eq!(shell.model().free, 6);
        test(&cx, &gate, &shell, &interior);
    });
}

fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-13),
        "{actual:?} != {expected:?}"
    );
}
fn real_near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-22),
        "{actual:e} != {expected:e}"
    );
}
fn work(force: &[C64], velocity: &[C64]) -> f64 {
    0.5 * force
        .iter()
        .zip(velocity)
        .map(|(&f, &v)| (f * v.conj()).re)
        .sum::<f64>()
}
fn product(matrix: &fs_sparse::Csr, field: &[C64], row: usize) -> C64 {
    let (columns, values) = matrix.row(row);
    columns
        .iter()
        .zip(values)
        .fold(C64::ZERO, |sum, (&j, &a)| sum + field[j].scale(a))
}
fn force() -> [C64; 6] {
    [
        C64::new(3e-8, -2e-8),
        C64::new(-2e-8, 1e-8),
        C64::new(7e-8, 3e-8),
        C64::new(2e-9, -3e-9),
        C64::new(-1e-9, 2e-9),
        C64::new(4e-11, 1e-11),
    ]
}

// Cubic-exact integration of the actual P1 center shape times P2 pressure.
// This oracle uses neither the exported coupling matrix nor boundary projection.
fn pressure_force(interior: &HeterogeneousHelmholtzFe<'_, '_>, pressure: &[C64]) -> C64 {
    let quadrature = [
        ([1.0 / 3.0; 3], -27.0 / 48.0),
        ([0.6, 0.2, 0.2], 25.0 / 48.0),
        ([0.2, 0.6, 0.2], 25.0 / 48.0),
        ([0.2, 0.2, 0.6], 25.0 / 48.0),
    ];
    let mut traction = C64::ZERO;
    for element in 0..3 {
        for (l, weight) in quadrature {
            let p = interior
                .pressure_at(pressure, element, [l[0], l[1], l[2], 0.0])
                .unwrap();
            traction = traction - p.scale(weight * l[2] / 6.0);
        }
    }
    traction
}

#[test]
fn original_p1_pressure_work_recovers_all_native_material_shell_coordinates() {
    fixture(|cx, _, shell, interior| {
        let heterogeneous = HeterogeneousHelmholtzFe::new(
            interior.geometry(),
            5,
            &mut |_, point: [f64; 3]| {
                Ok::<_, AcousticFeError>(Medium {
                    density_kg_m3: 1.1 + 0.2 * point[0],
                    sound_speed_m_s: 3.0 + 0.4 * point[1],
                })
            },
            &mut || false,
        )
        .unwrap();
        let interior = &heterogeneous;
        let applied = force();
        for omega in [0.7, 3.1] {
            let model =
                TensorShellInterior::new(interior, shell, &faces(), &[], omega, controls(), cx)
                    .unwrap();
            let response = model
                .prepare(cx)
                .unwrap()
                .solve(&applied, &vec![C64::ZERO; interior.ndof()], cx)
                .unwrap();
            let field = &response.field;
            assert!(core::ptr::eq(model.material_operator().source(), shell));
            assert_eq!(response.angular_frequency_rad_s, omega);
            assert_eq!(field.structural_displacement.len(), 6);
            assert_eq!(field.pressure_pa.len(), interior.ndof());
            assert!(interior.ndof() > interior.geometry().positions().len());
            let actual_traction = pressure_force(interior, &field.pressure_pa);
            near(field.fluid_force[2], actual_traction, 2e-9);
            for coordinate in [0, 1, 3, 4, 5] {
                assert_eq!(field.fluid_force[coordinate], C64::ZERO);
                assert!(field.structural_displacement[coordinate].abs() > 1e-16);
            }
            let mut total = applied;
            total[2] = total[2] + actual_traction;
            let native = shell
                .prepare_harmonic_response(omega, controls().material, cx)
                .unwrap()
                .solve(&total, cx)
                .unwrap();
            for coordinate in 0..6 {
                near(
                    field.structural_displacement[coordinate],
                    native.displacement[coordinate],
                    2e-8,
                );
                near(
                    field.structural_velocity[coordinate],
                    native.velocity[coordinate],
                    2e-8,
                );
            }
            real_near(
                field.average_structural_dissipation_w,
                native.average_viscous_power_w,
                2e-8,
            );
            let interface = 0.5 * (actual_traction * field.structural_velocity[2].conj()).re;
            assert!(
                (interface - field.average_fluid_to_structure_power_w).abs()
                    < 2e-9 * response.apparent_power_scale_w
            );
            // The fluid rows use original inverse-density/compliance matrices,
            // with p=-i*omega*psi and the same independently checked boundary work.
            for row in 0..interior.ndof() {
                let residual = product(interior.stiffness(), &field.pressure_potential_pa_s, row)
                    - product(interior.mass(), &field.pressure_potential_pa_s, row)
                        .scale(omega * omega)
                    + field.displaced_flow_m3_s[row];
                assert!(residual.abs() < 2e-13);
                near(
                    field.pressure_pa[row],
                    C64::new(
                        omega * field.pressure_potential_pa_s[row].im,
                        -omega * field.pressure_potential_pa_s[row].re,
                    ),
                    1e-10,
                );
            }
            assert!(response.relative_power_balance < 2e-9);
        }
    });
}

#[test]
fn retained_compression_stiffens_the_shell_and_detachment_recovers_the_dry_source() {
    fixture(|cx, _, shell, interior| {
        let omega = 0.03;
        let mut applied = [C64::ZERO; 6];
        applied[2] = C64::new(1e-7, -2e-8);
        let zero = vec![C64::ZERO; interior.ndof()];
        let wet = TensorShellInterior::new(interior, shell, &faces(), &[], omega, controls(), cx)
            .unwrap();
        let detached =
            TensorShellInterior::new(interior, shell, &[], &[], omega, controls(), cx).unwrap();
        let response = wet
            .prepare(cx)
            .unwrap()
            .solve(&applied, &zero, cx)
            .unwrap()
            .field;
        let dry = detached
            .prepare(cx)
            .unwrap()
            .solve(&applied, &zero, cx)
            .unwrap()
            .field;
        let native = shell
            .prepare_harmonic_response(omega, controls().material, cx)
            .unwrap()
            .solve(&applied, cx)
            .unwrap();
        assert_eq!(wet.structural_dofs(), 6);
        assert_eq!(wet.pressure_dofs(), interior.ndof());
        for coordinate in 0..6 {
            near(
                dry.structural_displacement[coordinate],
                native.displacement[coordinate],
                2e-8,
            );
        }
        assert!(dry.pressure_pa.iter().all(|p| p.abs() < 1e-18));
        assert!(
            response.structural_displacement[2].abs() < 0.95 * dry.structural_displacement[2].abs()
        );
        // V=1/6, integral(center P1 trace)=1/6 and bulk modulus=10.8.
        // Thus the exact uniform-compression limit is Kc=bulk*Aeff^2/V=1.8.
        let reaction = response.fluid_force[2] / response.structural_displacement[2];
        near(reaction, C64::from_re(-1.8), 3e-4);
        for element in 0..3 {
            let p = interior
                .pressure_at(&response.pressure_pa, element, [0.25; 4])
                .unwrap();
            near(p, response.structural_displacement[2].scale(10.8), 3e-4);
        }
        assert!(response.energy.fluid_compression_j > 0.0);
        // Reversing paired vertex winding cannot reverse the fluid-owned normal.
        let reversed = faces().map(|mut face| {
            face.fluid.swap(0, 1);
            face.shell.swap(0, 1);
            face
        });
        let reverse =
            TensorShellInterior::new(interior, shell, &reversed, &[], omega, controls(), cx)
                .unwrap();
        for (a, b) in wet.coupling().iter().zip(reverse.coupling()) {
            assert!((a - b).abs() < 1e-13);
        }
    });
}

#[test]
fn coherent_force_and_flow_preserve_signed_interface_work_and_separate_physical_heat() {
    fixture(|cx, _, shell, interior| {
        let wall = ResistiveBoundary {
            faces: vec![[0, 1, 3]],
            admittance_m_pa_s: 0.08,
        };
        let model =
            TensorShellInterior::new(interior, shell, &faces(), &[wall], 1.7, controls(), cx)
                .unwrap();
        let prepared = model.prepare(cx).unwrap();
        let force = force();
        let flow = model
            .project_inward_velocity(
                &[(
                    [1, 2, 3],
                    [
                        C64::new(2e-7, -1e-7),
                        C64::new(-3e-8, 8e-8),
                        C64::new(1e-7, 4e-8),
                    ],
                )],
                cx,
            )
            .unwrap();
        let mechanical = prepared
            .solve(&force, &vec![C64::ZERO; interior.ndof()], cx)
            .unwrap();
        let acoustic = prepared.solve(&[C64::ZERO; 6], &flow, cx).unwrap();
        assert!(mechanical.field.average_fluid_to_structure_power_w < 0.0);
        assert!(acoustic.field.average_fluid_to_structure_power_w > 0.0);
        let phase = C64::new(0.4, -0.7);
        let phased_flow: Vec<_> = flow.iter().map(|&z| phase * z).collect();
        let coherent = prepared.solve(&force, &phased_flow, cx).unwrap();
        for (sum, (a, b)) in coherent.field.structural_displacement.iter().zip(
            mechanical
                .field
                .structural_displacement
                .iter()
                .zip(&acoustic.field.structural_displacement),
        ) {
            near(*sum, *a + phase * *b, 3e-8);
        }
        for (sum, (a, b)) in coherent.field.pressure_pa.iter().zip(
            mechanical
                .field
                .pressure_pa
                .iter()
                .zip(&acoustic.field.pressure_pa),
        ) {
            near(*sum, *a + phase * *b, 3e-8);
        }
        for response in [&mechanical, &acoustic, &coherent] {
            let field = &response.field;
            let traction = pressure_force(interior, &field.pressure_pa);
            real_near(
                field.average_fluid_to_structure_power_w,
                0.5 * (traction * field.structural_velocity[2].conj()).re,
                3e-8,
            );
            let viscous_force: Vec<_> = (0..6)
                .map(|row| {
                    product(
                        model.material_operator().damping(),
                        &field.structural_velocity,
                        row,
                    )
                })
                .collect();
            real_near(
                field.average_structural_dissipation_w,
                work(&viscous_force, &field.structural_velocity),
                3e-8,
            );
            // Positive degree-four triangle quadrature integrates |P2 pressure|^2
            // on the independently known area-1/2 wall (tet 0, center omitted).
            let mut wall_heat = 0.0;
            for (a, weight) in [
                (0.445948490915965, 0.223381589678011),
                (0.091576213509771, 0.109951743655322),
            ] {
                for l in [
                    [a, a, 1.0 - 2.0 * a],
                    [a, 1.0 - 2.0 * a, a],
                    [1.0 - 2.0 * a, a, a],
                ] {
                    let p = interior
                        .pressure_at(&field.pressure_pa, 0, [l[0], l[1], 0.0, l[2]])
                        .unwrap();
                    wall_heat += 0.5 * 0.08 * 0.5 * weight * (p * p.conj()).re;
                }
            }
            assert!(wall_heat > 0.0 && field.average_structural_dissipation_w > 0.0);
            real_near(field.average_wall_dissipation_w, wall_heat, 3e-8);
            let tolerance = 2e-9 * response.apparent_power_scale_w;
            assert!(field.structural_power_balance_residual_w().abs() < tolerance);
            assert!(field.acoustic_power_balance_residual_w().abs() < tolerance);
            assert!(field.power_balance_residual_w.abs() < tolerance);
        }
    });
}

#[test]
fn interface_domains_original_strain_and_complete_coordinate_budgets_refuse() {
    fixture(|cx, gate, shell, interior| {
        let build = |interfaces: &[ConformingShellFace], omega, limits| {
            TensorShellInterior::new(interior, shell, interfaces, &[], omega, limits, cx)
        };
        let valid = faces();
        let mut mismatched = valid;
        mismatched[0].fluid.swap(0, 1);
        let mut internal = valid;
        internal[0].fluid = [0, 4, 3];
        assert!(build(&mismatched, 1.7, controls()).is_err());
        assert!(build(&internal, 1.7, controls()).is_err());
        assert!(build(&[valid[0]; 2], 1.7, controls()).is_err());
        assert!(matches!(
            build(&valid, 41.0, controls()),
            Err(TensorShellInteriorError::Shell(
                TensorShellHarmonicError::Material { .. }
            ))
        ));
        let mut limits = controls();
        limits.coupled.maximum_coordinates = shell.model().free + interior.ndof() - 1;
        assert!(build(&valid, 1.7, limits).is_err());
        for limits in [
            TensorShellInteriorControls {
                maximum_interface_faces: 2,
                ..controls()
            },
            TensorShellInteriorControls {
                maximum_trace_bytes: 1,
                ..controls()
            },
            TensorShellInteriorControls {
                maximum_trace_work: 1,
                ..controls()
            },
        ] {
            assert!(matches!(
                build(&valid, 1.7, limits),
                Err(TensorShellInteriorError::Budget { .. })
            ));
        }
        let model = build(&valid, 1.7, controls()).unwrap();
        let prepared = model.prepare(cx).unwrap();
        let zero = vec![C64::ZERO; interior.ndof()];
        assert!(prepared.solve(&[C64::ZERO; 5], &zero, cx).is_err());
        assert!(
            prepared
                .solve(&[C64::new(f64::NAN, 0.0); 6], &zero, cx)
                .is_err()
        );
        assert!(matches!(
            prepared.solve(&[C64::from_re(1e6); 6], &zero, cx),
            Err(TensorShellInteriorError::Shell(
                TensorShellHarmonicError::Material { .. }
            ))
        ));
        gate.request();
        assert!(prepared.solve(&force(), &zero, cx).is_err());
        assert!(model.prepare(cx).is_err());
        assert!(build(&valid, 1.7, controls()).is_err());
    });
}
