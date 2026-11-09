//! G1/G3/G4: reciprocal lossy-volume/wall fields, independent work and refusals.

use fs_feec::acoustics::helmholtz::heterogeneous::time::{
    ResistiveBoundary,
    mechanical::{MechanicalBoundary, traction::WallTractionProfile},
};
use fs_feec::acoustics::helmholtz::heterogeneous::{
    BoundaryAdmittance, EquivalentFluidCoefficients, EquivalentFluidControls, HeterogeneousError,
    HeterogeneousHelmholtzFe, MechanicalHarmonicLimits, PreparedEquivalentFluid, VelocityTrace,
};
use fs_feec::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_feec::fixtures::{masked_cube_grid, single_tet};
use fs_math::c64::C64;
use fs_phs::{WallPin, acoustic_outlet::AcousticWallRelaxation};
use fs_rep_mesh::TetComplex;
use std::cell::Cell;

const MEDIUM: Medium = Medium {
    density_kg_m3: 1.3,
    sound_speed_m_s: 2.0,
};
const MATERIAL: EquivalentFluidCoefficients = EquivalentFluidCoefficients {
    dynamic_density_kg_m3: C64::new(1.3, -0.4),
    dynamic_bulk_modulus_pa: C64::new(5.2, 0.8),
};
const AREA: f64 = 0.06;
const OMEGA: f64 = 2.2;

fn brick(nx: usize) -> (TetComplex, Vec<[f64; 3]>) {
    let (complex, mut points) = masked_cube_grid(nx, 1, 1, &|_, _, _| true);
    for p in &mut points {
        p[0] /= nx as f64;
        p[1] *= 0.2;
        p[2] *= 0.3;
    }
    (complex, points)
}

fn plane(fe: &HelmholtzFe<'_>, x: f64) -> Vec<[u32; 3]> {
    fe.boundary_faces()
        .into_iter()
        .filter(|face| face.iter().all(|&v| fe.positions()[v as usize][0] == x))
        .collect()
}

fn velocity(fe: &HelmholtzFe<'_>, value: C64) -> Vec<VelocityTrace> {
    plane(fe, 0.0)
        .into_iter()
        .map(|face| (face, [value; 3]))
        .collect()
}

fn wall(faces: Vec<[u32; 3]>) -> MechanicalBoundary {
    MechanicalBoundary {
        faces,
        wall: WallPin {
            surface_density: 0.6,
            stiffness_per_area: 3.0,
            resistance: 0.4,
        },
        relaxation: vec![
            AcousticWallRelaxation {
                stiffness_per_area: 2.5,
                relaxation_time_s: 0.3,
            },
            AcousticWallRelaxation {
                stiffness_per_area: 0.7,
                relaxation_time_s: 0.08,
            },
        ],
    }
}

fn prepare<'f, 'c>(
    geometry: &'f HelmholtzFe<'c>,
    wall: &MechanicalBoundary,
) -> PreparedEquivalentFluid<'f, 'c> {
    PreparedEquivalentFluid::new_with_coupled_mechanical_walls(
        geometry,
        OMEGA / core::f64::consts::TAU,
        &[],
        core::slice::from_ref(wall),
        &[],
        EquivalentFluidControls::default(),
        &mut |_, _| Ok::<_, HeterogeneousError>(MATERIAL),
        &mut || false,
    )
    .unwrap()
}

fn sin(z: C64) -> C64 {
    C64::new(z.re.sin() * z.im.cosh(), z.re.cos() * z.im.sinh())
}
fn cos(z: C64) -> C64 {
    C64::new(z.re.cos() * z.im.cosh(), -z.re.sin() * z.im.sinh())
}

/// Independently solve the continuum layer transfer matrix and wall force law.
/// p(L) + traction = Z_wall v(L), positive v along the original x coordinate.
fn exact(wall: &MechanicalBoundary, inlet: C64, traction: C64, x: f64) -> (C64, C64) {
    let rho = MATERIAL.dynamic_density_kg_m3;
    let bulk = MATERIAL.dynamic_bulk_modulus_pa;
    let k = (rho / bulk).sqrt().scale(OMEGA);
    let z = (rho * bulk).sqrt();
    let mut zw = C64::new(
        wall.wall.resistance,
        OMEGA * wall.wall.surface_density - wall.wall.stiffness_per_area / OMEGA,
    );
    for arm in &wall.relaxation {
        zw = zw
            + C64::from_re(arm.stiffness_per_area) / C64::new(1.0 / arm.relaxation_time_s, OMEGA);
    }
    let imaginary = C64::new(0.0, 1.0);
    let end_velocity =
        (inlet + imaginary * sin(k) * traction / z) / (cos(k) + imaginary * sin(k) * zw / z);
    let end_pressure = zw * end_velocity - traction;
    let angle = k.scale(1.0 - x);
    (
        cos(angle) * end_pressure + imaginary * z * sin(angle) * end_velocity,
        imaginary * sin(angle) * end_pressure / z + cos(angle) * end_velocity,
    )
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance * expected.abs().max(1.0),
        "actual={actual:?}, expected={expected:?}, tolerance={tolerance:e}"
    );
}

fn scalar_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance * expected.abs().max(1e-12),
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn lossy_slab_and_all_wall_memories_converge_to_independent_continuum_work() {
    let mut errors = Vec::new();
    for nx in [2, 4, 8] {
        let (mesh, points) = brick(nx);
        let geometry = HelmholtzFe::new(&mesh, &points, 2, MEDIUM).unwrap();
        let boundary = wall(plane(&geometry, 1.0));
        let prepared = prepare(&geometry, &boundary);
        let response = prepared
            .solve(&velocity(&geometry, C64::ONE), &mut || false)
            .unwrap();
        assert_eq!(prepared.wall_count(), 1);
        assert_eq!(response.walls[0].relaxation_memory.len(), 2);
        assert_eq!(response.pressure_pa.len(), geometry.ndof());
        let error = points
            .iter()
            .zip(&response.pressure_pa)
            .map(|(point, &p)| (p - exact(&boundary, C64::ONE, C64::ZERO, point[0]).0).abs())
            .fold(0.0, f64::max);
        errors.push(error);
        assert!(response.power_balance_residual_w().abs() < 2e-11);
        assert!(response.fluid_power_balance_residual_w().abs() < 2e-11);
        assert!(response.wall_accounting.power_balance_residual_w().abs() < 2e-11);
        assert!(response.backward_error < 5e-13);
        if nx != 8 {
            continue;
        }
        let expected_v = exact(&boundary, C64::ONE, C64::ZERO, 1.0).1;
        for (local, &global) in prepared.wall_dofs(0).unwrap().iter().enumerate() {
            if global >= points.len() {
                continue;
            }
            close(
                response.walls[0].normal_velocity_m_s[local],
                expected_v,
                2e-4,
            );
            close(
                response.walls[0].displacement_m[local],
                expected_v / C64::new(0.0, OMEGA),
                2e-4,
            );
            for (arm, memory) in boundary
                .relaxation
                .iter()
                .zip(&response.walls[0].relaxation_memory)
            {
                let rate = 1.0 / arm.relaxation_time_s;
                let expected_z = expected_v.scale(arm.stiffness_per_area.sqrt() * rate)
                    / (C64::new(0.0, OMEGA) * C64::new(rate, OMEGA));
                close(memory[local], expected_z, 2e-4);
            }
        }
        let wall_heat = AREA * 0.5 * expected_v.norm_sq();
        scalar_close(
            response.wall_accounting.resistance_dissipation_w,
            wall_heat * boundary.wall.resistance,
            3e-4,
        );
        let relaxation_resistance: f64 = boundary
            .relaxation
            .iter()
            .map(|arm| {
                let rate = 1.0 / arm.relaxation_time_s;
                arm.stiffness_per_area * rate / (rate * rate + OMEGA * OMEGA)
            })
            .sum();
        scalar_close(
            response.wall_accounting.relaxation_dissipation_w,
            wall_heat * relaxation_resistance,
            3e-4,
        );
        // Simpson integration of the independent exact p/v solution. Density
        // heat uses -omega Im(rho)|v|²/2, not the FE pressure-gradient form.
        let mut density = 0.0;
        let mut compression = 0.0;
        for i in 0..=512 {
            let weight = if i == 0 || i == 512 {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            let (p, v) = exact(&boundary, C64::ONE, C64::ZERO, i as f64 / 512.0);
            density += weight * v.norm_sq();
            compression += weight * p.norm_sq();
        }
        let measure = AREA / (3.0 * 512.0) * 0.5 * OMEGA;
        scalar_close(
            response.density_dissipation_w,
            -measure * MATERIAL.dynamic_density_kg_m3.im * density,
            3e-4,
        );
        scalar_close(
            response.compressional_dissipation_w,
            -measure * MATERIAL.dynamic_bulk_modulus_pa.recip().im * compression,
            3e-4,
        );
        scalar_close(
            response.input_power_w,
            0.5 * AREA * exact(&boundary, C64::ONE, C64::ZERO, 0.0).0.re,
            3e-4,
        );
    }
    assert!(
        errors[1] < errors[0] / 4.0 && errors[2] < errors[1] / 4.0,
        "{errors:?}"
    );
    assert!(errors[2] < 2e-4, "{errors:?}");
}

#[test]
fn coherent_wall_traction_drives_the_lossy_fluid_and_keeps_owner_and_work() {
    let (mesh, points) = brick(4);
    let geometry = HelmholtzFe::new(&mesh, &points, 2, MEDIUM).unwrap();
    let boundary = wall(plane(&geometry, 1.0));
    let prepared = prepare(&geometry, &boundary);
    let profile = WallTractionProfile {
        wall: 0,
        coefficients: points
            .iter()
            .enumerate()
            .filter(|(_, p)| p[0] == 1.0)
            .map(|(i, _)| (i, 1.0))
            .collect(),
    };
    let source = prepared
        .prepare_wall_traction(&profile, &mut || false)
        .unwrap();
    let inlet = C64::new(0.7, -0.25);
    let traction = C64::new(0.4, 0.15);
    let fluid = prepared
        .solve(&velocity(&geometry, inlet), &mut || false)
        .unwrap();
    let mechanical = prepared
        .solve_with_wall_tractions(&[], &[(&source, traction)], &mut || false)
        .unwrap();
    let both = prepared
        .solve_with_wall_tractions(
            &velocity(&geometry, inlet),
            &[(&source, traction)],
            &mut || false,
        )
        .unwrap();
    for ((&p, &a), &b) in both
        .pressure_pa
        .iter()
        .zip(&fluid.pressure_pa)
        .zip(&mechanical.pressure_pa)
    {
        close(p, a + b, 2e-11);
    }
    for (point, &p) in points.iter().zip(&both.pressure_pa) {
        close(p, exact(&boundary, inlet, traction, point[0]).0, 2e-3);
    }
    let expected_v = exact(&boundary, inlet, traction, 1.0).1;
    scalar_close(
        both.wall_accounting.wall_traction_power_w,
        0.5 * AREA * (traction.conj() * expected_v).re,
        2e-3,
    );
    assert_eq!(mechanical.input_power_w, 0.0);
    assert!(mechanical.wall_accounting.fluid_to_wall_power_w < 0.0);
    for response in [&fluid, &mechanical, &both] {
        assert!(response.power_balance_residual_w().abs() < 2e-11);
        assert!(response.fluid_power_balance_residual_w().abs() < 2e-11);
        assert!(response.wall_accounting.power_balance_residual_w().abs() < 2e-11);
    }
    assert!(
        (both.density_dissipation_w
            - fluid.density_dissipation_w
            - mechanical.density_dissipation_w)
            .abs()
            > 1e-5
    );
    let phase = C64::new(0.6, 0.8);
    let rotated = prepared
        .solve_with_wall_tractions(
            &velocity(&geometry, inlet * phase),
            &[(&source, traction * phase)],
            &mut || false,
        )
        .unwrap();
    for (&actual, &original) in rotated.pressure_pa.iter().zip(&both.pressure_pa) {
        close(actual, original * phase, 2e-11);
    }
    scalar_close(
        rotated.density_dissipation_w,
        both.density_dissipation_w,
        2e-11,
    );
    assert!(matches!(
        prepared.solve_with_wall_tractions(&[], &[(&source, traction)], &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    let other = prepare(&geometry, &boundary);
    assert!(
        other
            .solve_with_wall_tractions(&[], &[(&source, traction)], &mut || false)
            .is_err()
    );
    let repeat = prepared
        .solve_with_wall_tractions(
            &velocity(&geometry, inlet),
            &[(&source, traction)],
            &mut || false,
        )
        .unwrap();
    assert_eq!(repeat.pressure_pa, both.pressure_pa);
}

#[test]
fn zero_local_wall_impedance_keeps_finite_motion_and_original_soft_pressure_limit() {
    let (mesh, points) = brick(2);
    let geometry = HelmholtzFe::new(&mesh, &points, 2, MEDIUM).unwrap();
    let boundary = MechanicalBoundary {
        faces: plane(&geometry, 1.0),
        wall: WallPin {
            surface_density: 1.0,
            stiffness_per_area: OMEGA * OMEGA,
            resistance: 0.0,
        },
        relaxation: vec![],
    };
    let prepared = prepare(&geometry, &boundary);
    let response = prepared
        .solve(&velocity(&geometry, C64::ONE), &mut || false)
        .unwrap();
    let soft = PreparedEquivalentFluid::new(
        &geometry,
        OMEGA / core::f64::consts::TAU,
        &[],
        &boundary.faces,
        EquivalentFluidControls::default(),
        &mut |_, _| Ok::<_, HeterogeneousError>(MATERIAL),
        &mut || false,
    )
    .unwrap()
    .solve(&velocity(&geometry, C64::ONE), &mut || false)
    .unwrap();
    for (&actual, &expected) in response.pressure_pa.iter().zip(&soft.pressure_pa) {
        close(actual, expected, 2e-10);
    }
    assert!(
        response.walls[0]
            .normal_velocity_m_s
            .iter()
            .any(|v| v.abs() > 0.1)
    );
    assert!(response.wall_accounting.kinetic_energy_j > 0.0);
    scalar_close(
        response.wall_accounting.kinetic_energy_j,
        response.wall_accounting.elastic_energy_j,
        2e-12,
    );
    assert_eq!(response.wall_accounting.dissipation_w(), 0.0);
    assert!(response.density_dissipation_w > 0.0);
    assert!(response.power_balance_residual_w().abs() < 2e-11);
}

#[test]
fn lossless_volume_agrees_with_original_mechanical_owner_and_adjacent_release() {
    let (mesh, points) = single_tet();
    let geometry = HelmholtzFe::new(&mesh, &points, 2, MEDIUM).unwrap();
    let boundary = wall(vec![[0, 1, 3]]);
    let release = [[0, 1, 2]];
    let robin = [BoundaryAdmittance {
        faces: vec![[0, 2, 3]],
        admittance_m_pa_s: C64::from_re(0.12),
    }];
    let prepared = PreparedEquivalentFluid::new_with_coupled_mechanical_walls(
        &geometry,
        OMEGA / core::f64::consts::TAU,
        &robin,
        core::slice::from_ref(&boundary),
        &release,
        EquivalentFluidControls::default(),
        &mut |_, _| {
            Ok::<_, HeterogeneousError>(EquivalentFluidCoefficients {
                dynamic_density_kg_m3: C64::from_re(MEDIUM.density_kg_m3),
                dynamic_bulk_modulus_pa: C64::from_re(
                    MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s.powi(2),
                ),
            })
        },
        &mut || false,
    )
    .unwrap();
    let original = HeterogeneousHelmholtzFe::new(
        &geometry,
        4,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap();
    let original = original
        .prepare_with_coupled_mechanical_walls(
            OMEGA / core::f64::consts::TAU,
            &[ResistiveBoundary {
                faces: robin[0].faces.clone(),
                admittance_m_pa_s: 0.12,
            }],
            core::slice::from_ref(&boundary),
            &release,
            MechanicalHarmonicLimits::default(),
            &mut || false,
        )
        .unwrap();
    let traces = [([1, 2, 3], [C64::new(0.2, -0.1); 3])];
    let actual = prepared.solve(&traces, &mut || false).unwrap();
    let expected = original.solve(&traces, &mut || false).unwrap();
    assert_eq!(prepared.wall_dofs(0), original.wall_dofs(0));
    assert!(prepared.wall_dofs(0).unwrap().contains(&0));
    assert_eq!(actual.pressure_pa[0], C64::ZERO);
    for (&a, &b) in actual
        .pressure_pa
        .iter()
        .zip(&expected.acoustic.pressure_pa)
    {
        close(a, b, 2e-11);
    }
    for (&a, &b) in actual.walls[0]
        .normal_velocity_m_s
        .iter()
        .zip(&expected.walls[0].normal_velocity_m_s)
    {
        close(a, b, 2e-11);
    }
    scalar_close(
        actual.wall_accounting.dissipation_w(),
        expected.wall_accounting.dissipation_w(),
        2e-11,
    );
    let original_heat = original
        .wall_dissipation(&expected.walls, &mut || false)
        .unwrap();
    let actual_heat = prepared
        .wall_dissipation(&actual.walls, &mut || false)
        .unwrap();
    scalar_close(actual_heat[0].area_m2, original_heat[0].area_m2, 2e-13);
    scalar_close(
        actual_heat[0].total_dissipation_w(),
        original_heat[0].total_dissipation_w(),
        2e-11,
    );
    assert_eq!(actual.density_dissipation_w, 0.0);
    assert_eq!(actual.compressional_dissipation_w, 0.0);
}

#[test]
fn complete_state_and_dense_limits_and_face_conflicts_precede_material_queries() {
    let (mesh, points) = single_tet();
    let geometry = HelmholtzFe::new(&mesh, &points, 2, MEDIUM).unwrap();
    let boundary = wall(vec![[0, 1, 3]]);
    let calls = Cell::new(0);
    let mut material = |_: usize, _: [f64; 3]| {
        calls.set(calls.get() + 1);
        Ok::<_, HeterogeneousError>(MATERIAL)
    };
    for dense in [false, true] {
        let mut controls = EquivalentFluidControls::default();
        if dense {
            controls.harmonic.max_dense_bytes = 32 * geometry.ndof().pow(2);
        } else {
            controls.harmonic.max_dofs = geometry.ndof();
        }
        let result = PreparedEquivalentFluid::new_with_coupled_mechanical_walls(
            &geometry,
            OMEGA / core::f64::consts::TAU,
            &[],
            core::slice::from_ref(&boundary),
            &[],
            controls,
            &mut material,
            &mut || false,
        );
        assert!(matches!(result, Err(HeterogeneousError::Limit { .. })));
        assert_eq!(calls.get(), 0);
    }
    let result = PreparedEquivalentFluid::new_with_coupled_mechanical_walls(
        &geometry,
        OMEGA / core::f64::consts::TAU,
        &[],
        core::slice::from_ref(&boundary),
        &boundary.faces,
        EquivalentFluidControls::default(),
        &mut material,
        &mut || false,
    );
    assert!(result.is_err());
    assert_eq!(calls.get(), 0);
}

#[test]
fn original_patch_area_mass_and_positive_heat_keep_hierarchical_and_local_identity() {
    let (mesh, points) = brick(2);
    let geometry = HelmholtzFe::new(&mesh, &points, 2, MEDIUM).unwrap();
    let side_faces = geometry
        .boundary_faces()
        .into_iter()
        .filter(|face| face.iter().all(|&v| points[v as usize][1] == 0.0))
        .collect();
    let mut side = wall(side_faces);
    side.wall.surface_density = 0.9;
    side.wall.resistance = 0.8;
    let boundaries = [wall(plane(&geometry, 1.0)), side];
    let prepared = PreparedEquivalentFluid::new_with_coupled_mechanical_walls(
        &geometry,
        OMEGA / core::f64::consts::TAU,
        &[],
        &boundaries,
        &[],
        EquivalentFluidControls::default(),
        &mut |_, _| Ok::<_, HeterogeneousError>(MATERIAL),
        &mut || false,
    )
    .unwrap();
    let response = prepared
        .solve(&velocity(&geometry, C64::ONE), &mut || false)
        .unwrap();
    let native = prepared
        .wall_dissipation(&response.walls, &mut || false)
        .unwrap();
    let areas = [AREA, 0.3];
    for (index, power) in native.iter().enumerate() {
        scalar_close(power.area_m2, areas[index], 2e-13);
        scalar_close(
            power.mass_kg,
            areas[index] * boundaries[index].wall.surface_density,
            2e-13,
        );
        assert!(power.resistance_dissipation_w > 0.0 && power.relaxation_dissipation_w > 0.0);
    }
    scalar_close(
        native.iter().map(|p| p.resistance_dissipation_w).sum(),
        response.wall_accounting.resistance_dissipation_w,
        2e-13,
    );
    scalar_close(
        native.iter().map(|p| p.relaxation_dissipation_w).sum(),
        response.wall_accounting.relaxation_dissipation_w,
        2e-13,
    );

    // Authored coherent constant motion has vertex coefficients only. Check
    // each wall's physical loss against its area and independent local law;
    // hierarchical bubble coordinates must not contribute a spurious area.
    let velocities = [C64::new(0.2, -0.3), C64::new(-0.4, 0.1)];
    let mut fields = response.walls.clone();
    for (index, field) in fields.iter_mut().enumerate() {
        for (local, &global) in prepared.wall_dofs(index).unwrap().iter().enumerate() {
            let v = if global < points.len() {
                velocities[index]
            } else {
                C64::ZERO
            };
            let displacement = v / C64::new(0.0, OMEGA);
            field.normal_velocity_m_s[local] = v;
            field.displacement_m[local] = displacement;
            for (arm, law) in boundaries[index].relaxation.iter().enumerate() {
                let coupling = law.stiffness_per_area.sqrt();
                let e = v.scale(coupling) / C64::new(1.0 / law.relaxation_time_s, OMEGA);
                field.relaxation_extension[arm][local] = e;
                field.relaxation_memory[arm][local] = displacement.scale(coupling) - e;
            }
        }
    }
    let powers = prepared.wall_dissipation(&fields, &mut || false).unwrap();
    for (index, power) in powers.iter().enumerate() {
        let mean_velocity_squared = 0.5 * areas[index] * velocities[index].norm_sq();
        let relaxation: f64 = boundaries[index]
            .relaxation
            .iter()
            .map(|arm| {
                let rate = 1.0 / arm.relaxation_time_s;
                arm.stiffness_per_area * rate / (rate * rate + OMEGA * OMEGA)
            })
            .sum();
        scalar_close(
            power.resistance_dissipation_w,
            mean_velocity_squared * boundaries[index].wall.resistance,
            2e-13,
        );
        scalar_close(
            power.relaxation_dissipation_w,
            mean_velocity_squared * relaxation,
            2e-13,
        );
    }
    let mut transformed = fields.clone();
    let field = &mut transformed[0];
    for value in field
        .displacement_m
        .iter_mut()
        .chain(&mut field.normal_velocity_m_s)
        .chain(field.relaxation_memory.iter_mut().flatten())
        .chain(field.relaxation_extension.iter_mut().flatten())
    {
        *value = *value * C64::new(0.0, 2.0);
    }
    let changed = prepared
        .wall_dissipation(&transformed, &mut || false)
        .unwrap();
    scalar_close(
        changed[0].total_dissipation_w(),
        4.0 * powers[0].total_dissipation_w(),
        2e-13,
    );
    assert_eq!(changed[1], powers[1]);
    assert!(
        prepared
            .wall_dissipation(&fields[..1], &mut || false)
            .is_err()
    );
    let mut broken = fields.clone();
    broken[0].normal_velocity_m_s.pop();
    assert!(prepared.wall_dissipation(&broken, &mut || false).is_err());
    assert!(matches!(
        prepared.wall_dissipation(&fields, &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    assert_eq!(
        prepared.wall_dissipation(&fields, &mut || false).unwrap(),
        powers
    );
}
