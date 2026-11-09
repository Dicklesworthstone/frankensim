//! Original spatial drives in the complete lossy pressure/wall solve.
//! Numerical materials are explicit fixtures, not calibrated specimens.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_feec::acoustics::helmholtz::heterogeneous::time::mechanical::{
    MechanicalBoundary, traction::WallTractionProfile,
};
use fs_feec::acoustics::helmholtz::heterogeneous::{
    BoundaryAdmittance, EquivalentFluidCoefficients, EquivalentFluidControls,
    EquivalentFluidPressureField, HarmonicLimits, HeterogeneousError, HeterogeneousHelmholtzFe,
    MAX_BOUNDARY_SOURCES, MAX_VOLUME_SOURCES, PreparedEquivalentFluid,
};
use fs_feec::acoustics::helmholtz::{
    AcousticFeError, HelmholtzFe, Medium, QuadraticGeometryControls,
};
use fs_feec::fixtures::single_tet;
use fs_math::c64::C64;
use fs_phs::{WallPin, acoustic_outlet::AcousticWallRelaxation};
use std::cell::Cell;

const MEDIUM: Medium = Medium {
    density_kg_m3: 1.3,
    sound_speed_m_s: 2.0,
};
const MATERIAL: EquivalentFluidCoefficients = EquivalentFluidCoefficients {
    dynamic_density_kg_m3: C64::new(1.3, -0.4),
    dynamic_bulk_modulus_pa: C64::new(5.2, 0.8),
};
const OMEGA: f64 = 2.2;
const VOLUME: f64 = 1.0 / 6.0;
const FACE: [u32; 3] = [0, 1, 2];
const MIDSIDE: [f64; 6] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];

fn prepare<'f, 'c>(
    geometry: &'f HelmholtzFe<'c>,
    material: EquivalentFluidCoefficients,
    walls: &[MechanicalBoundary],
    admittances: &[BoundaryAdmittance],
    release: &[[u32; 3]],
) -> PreparedEquivalentFluid<'f, 'c> {
    PreparedEquivalentFluid::new_with_coupled_mechanical_walls(
        geometry,
        OMEGA / core::f64::consts::TAU,
        admittances,
        walls,
        release,
        EquivalentFluidControls::default(),
        &mut |_, _| Ok::<_, HeterogeneousError>(material),
        &mut || false,
    )
    .unwrap()
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-20),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual={actual:?}, expected={expected:?}"
    );
}

fn constant(field: &EquivalentFluidPressureField, pressure: C64, vertices: usize) {
    for (i, &value) in field.pressure_pa.iter().enumerate() {
        close(
            value,
            if i < vertices { pressure } else { C64::ZERO },
            4e-12,
        );
    }
}

fn balance(field: &EquivalentFluidPressureField) {
    let scale = field.input_power_w.abs()
        + field.wall_accounting.wall_traction_power_w.abs()
        + field.boundary_dissipation_w
        + field.density_dissipation_w
        + field.compressional_dissipation_w
        + field.wall_accounting.resistance_dissipation_w
        + field.wall_accounting.relaxation_dissipation_w;
    assert!(field.power_balance_residual_w().abs() < 3e-11 * scale.max(1e-20));
    assert!(field.fluid_power_balance_residual_w().abs() < 3e-11 * scale.max(1e-20));
    assert!(field.wall_accounting.power_balance_residual_w().abs() < 3e-11 * scale.max(1e-20));
}

#[test]
fn original_volume_injection_recovers_complex_compliance_and_signed_boundary_work() {
    let (mesh, positions) = single_tet();
    let geometry = HelmholtzFe::new(&mesh, &positions, 2, MEDIUM).unwrap();
    let material_calls = Cell::new(0);
    let prepared = PreparedEquivalentFluid::new(
        &geometry,
        OMEGA / core::f64::consts::TAU,
        &[],
        &[],
        EquivalentFluidControls::default(),
        &mut |_, _| {
            material_calls.set(material_calls.get() + 1);
            Ok::<_, HeterogeneousError>(MATERIAL)
        },
        &mut || false,
    )
    .unwrap();
    let source_calls = Cell::new(0);
    let volume = prepared
        .volume_velocity_source(
            4,
            &mut |_, _| {
                source_calls.set(source_calls.get() + 1);
                Ok::<_, HeterogeneousError>(1.0 / VOLUME)
            },
            &mut || false,
        )
        .unwrap();
    let boundary = prepared
        .quadratic_velocity_source(&[(FACE, [1.0; 6])], &mut || false)
        .unwrap();
    relative(volume.integrated_shape(), 1.0, 3e-14);
    relative(boundary.integrated_shape_m2(), 0.5, 3e-14);
    let counts = (material_calls.get(), source_calls.get());
    let q = C64::new(0.02, -0.007);
    // Continuity gives p=K*Q/(i*omega*V) throughout this closed homogeneous cell.
    let expected = MATERIAL.dynamic_bulk_modulus_pa * q / C64::new(0.0, OMEGA * VOLUME);
    let work = MATERIAL.dynamic_bulk_modulus_pa.im * q.norm_sq() / (2.0 * OMEGA * VOLUME);
    let field = prepared
        .solve_with_sources(&[], &[], &[(&volume, q)], &[], &mut || false)
        .unwrap();
    constant(&field, expected, mesh.vertex_count);
    relative(field.volume_source_power_w, work, 2e-11);
    relative(field.input_power_w, work, 2e-11);
    relative(field.compressional_dissipation_w, work, 2e-11);
    assert!(field.density_dissipation_w < 3e-12 * work);
    assert_eq!(field.boundary_source_power_w, 0.0);
    balance(&field);

    // Equal and opposite raw/prepared motion leaves the physical field intact,
    // while the prepared source independently ABSORBS one quarter of the work.
    let inward = q.scale(-0.25 / 0.5);
    let coherent = prepared
        .solve_with_sources(
            &[(FACE, [inward.scale(-1.0); 3])],
            &[(&boundary, inward)],
            &[(&volume, q)],
            &[],
            &mut || false,
        )
        .unwrap();
    constant(&coherent, expected, mesh.vertex_count);
    relative(coherent.boundary_source_power_w, -0.25 * work, 2e-11);
    relative(coherent.volume_source_power_w, work, 2e-11);
    relative(coherent.input_power_w, work, 2e-11);
    balance(&coherent);
    assert_eq!(
        (material_calls.get(), source_calls.get()),
        counts,
        "warm solves cannot resample either material or the original source"
    );
}

#[test]
fn midside_motion_zero_net_injection_and_all_wall_memories_share_one_coherent_solve() {
    let (mesh, positions) = single_tet();
    let geometry = HelmholtzFe::new(&mesh, &positions, 2, MEDIUM).unwrap();
    let wall = MechanicalBoundary {
        faces: vec![[1, 2, 3]],
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
    };
    let prepared = prepare(&geometry, MATERIAL, &[wall], &[], &[]);
    let boundary = prepared
        .quadratic_velocity_source(&[(FACE, MIDSIDE)], &mut || false)
        .unwrap();
    relative(boundary.integrated_shape_m2(), 1.0 / 6.0, 3e-14);
    assert!(boundary.has_nonzero_projection());
    let volume = prepared
        .volume_velocity_source(
            4,
            &mut |_, x| Ok::<_, HeterogeneousError>(6.0 * (x[0] - 0.25)),
            &mut || false,
        )
        .unwrap();
    assert!(volume.integrated_shape().abs() < 3e-14);
    let traction = prepared
        .prepare_wall_traction(
            &WallTractionProfile {
                wall: 0,
                coefficients: vec![(1, 1.0), (2, 1.0), (3, 1.0)],
            },
            &mut || false,
        )
        .unwrap();
    let raw = [([0, 2, 3], [C64::new(0.02, 0.01); 3])];
    let u = C64::new(0.03, -0.01);
    let q = C64::new(0.04, 0.02);
    let force = C64::new(0.07, -0.03);
    let original = prepared
        .solve_with_wall_tractions(&raw, &[(&traction, force)], &mut || false)
        .unwrap();
    let delegated = prepared
        .solve_with_sources(&raw, &[], &[], &[(&traction, force)], &mut || false)
        .unwrap();
    assert_eq!(original.pressure_pa, delegated.pressure_pa);
    assert_eq!(original.input_power_w, delegated.input_power_w);
    assert_eq!(delegated.boundary_source_power_w, 0.0);
    assert_eq!(delegated.volume_source_power_w, 0.0);
    let boundary_only = prepared
        .solve_with_sources(&[], &[(&boundary, u)], &[], &[], &mut || false)
        .unwrap();
    let volume_only = prepared
        .solve_with_sources(&[], &[], &[(&volume, q)], &[], &mut || false)
        .unwrap();
    assert!(boundary_only.pressure_pa.iter().any(|p| p.abs() > 1e-6));
    assert!(
        volume_only.density_dissipation_w > 0.0,
        "a zero-net source still excites spatial pressure gradients"
    );
    let combined = prepared
        .solve_with_sources(
            &raw,
            &[(&boundary, u)],
            &[(&volume, q)],
            &[(&traction, force)],
            &mut || false,
        )
        .unwrap();
    for i in 0..geometry.ndof() {
        close(
            combined.pressure_pa[i],
            original.pressure_pa[i] + boundary_only.pressure_pa[i] + volume_only.pressure_pa[i],
            2e-12,
        );
    }
    assert_eq!(combined.walls[0].relaxation_extension.len(), 2);
    for i in 0..combined.walls[0].normal_velocity_m_s.len() {
        close(
            combined.walls[0].normal_velocity_m_s[i],
            original.walls[0].normal_velocity_m_s[i]
                + boundary_only.walls[0].normal_velocity_m_s[i]
                + volume_only.walls[0].normal_velocity_m_s[i],
            2e-12,
        );
    }
    assert!(combined.wall_accounting.resistance_dissipation_w > 0.0);
    assert!(combined.wall_accounting.relaxation_dissipation_w > 0.0);
    for field in [&original, &boundary_only, &volume_only, &combined] {
        balance(field);
    }
    let cancelled = prepared
        .solve_with_sources(
            &[],
            &[],
            &[(&volume, q), (&volume, q.scale(-1.0))],
            &[],
            &mut || false,
        )
        .unwrap();
    assert!(cancelled.pressure_pa.iter().all(|&p| p == C64::ZERO));
    assert_eq!(cancelled.input_power_w, 0.0);
}

#[test]
fn curved_geometry_retains_original_volume_and_quadratic_boundary_moments() {
    let (mesh, original) = single_tet();
    let a = 0.6;
    let map = |[x, y, z]: [f64; 3]| [x, y, z * (1.0 + a * x)];
    let positions = original.iter().copied().map(map).collect::<Vec<_>>();
    let midpoints = mesh
        .edges
        .iter()
        .map(|&[i, j]| {
            map(core::array::from_fn(|axis| {
                f64::midpoint(original[i as usize][axis], original[j as usize][axis])
            }))
        })
        .collect::<Vec<_>>();
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 19,
                kernel_id: 93,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let geometry = HelmholtzFe::new_quadratic(
            &mesh,
            &positions,
            &midpoints,
            2,
            7,
            MEDIUM,
            QuadraticGeometryControls::default(),
            &cx,
        )
        .unwrap();
        let prepared = prepare(&geometry, MATERIAL, &[], &[], &[]);
        let volume = prepared
            .volume_velocity_source(7, &mut |_, _| Ok::<_, HeterogeneousError>(1.0), &mut || {
                false
            })
            .unwrap();
        let boundary = prepared
            .quadratic_velocity_source_with_quadrature(&[([0, 1, 3], MIDSIDE)], 7, &mut || false)
            .unwrap();
        // det(J)=1+a*x. On y=0, the same factor multiplies physical area.
        // Exact simplex moments: integral 4*lambda0*lambda1 = 1/6,
        // integral 4*lambda0*lambda1^2 = 1/15 on the reference face.
        let physical_volume = 1.0 / 6.0 + a / 24.0;
        relative(volume.integrated_shape(), physical_volume, 3e-13);
        relative(boundary.integrated_shape_m2(), 1.0 / 6.0 + a / 15.0, 3e-13);
        let q = C64::new(0.02, -0.007);
        let field = prepared
            .solve_with_sources(&[], &[], &[(&volume, q)], &[], &mut || false)
            .unwrap();
        constant(
            &field,
            MATERIAL.dynamic_bulk_modulus_pa * q / C64::new(0.0, OMEGA),
            mesh.vertex_count,
        );
        relative(
            field.volume_source_power_w,
            physical_volume * MATERIAL.dynamic_bulk_modulus_pa.im * q.norm_sq() / (2.0 * OMEGA),
            3e-11,
        );
        balance(&field);
    });
}

#[test]
fn lossless_limit_retains_original_source_projection_and_full_pressure_response() {
    let (mesh, positions) = single_tet();
    let geometry = HelmholtzFe::new(&mesh, &positions, 2, MEDIUM).unwrap();
    let real = HeterogeneousHelmholtzFe::new(
        &geometry,
        4,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap();
    let admittance = [BoundaryAdmittance {
        faces: vec![[1, 2, 3]],
        admittance_m_pa_s: C64::new(0.4, 0.2),
    }];
    let numerical = EquivalentFluidCoefficients {
        dynamic_density_kg_m3: C64::from_re(MEDIUM.density_kg_m3),
        dynamic_bulk_modulus_pa: C64::from_re(
            MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s.powi(2),
        ),
    };
    let equivalent = prepare(&geometry, numerical, &[], &admittance, &[]);
    let original = real
        .prepare(
            OMEGA / core::f64::consts::TAU,
            &admittance,
            HarmonicLimits::default(),
            &mut || false,
        )
        .unwrap();
    let volume = equivalent
        .volume_velocity_source(
            4,
            &mut |_, p| Ok::<_, HeterogeneousError>(1.0 + p[0]),
            &mut || false,
        )
        .unwrap();
    let original_volume = real
        .volume_velocity_source(
            4,
            &mut |_, p| Ok::<_, HeterogeneousError>(1.0 + p[0]),
            &mut || false,
        )
        .unwrap();
    let boundary = equivalent
        .quadratic_velocity_source(&[(FACE, MIDSIDE)], &mut || false)
        .unwrap();
    let original_boundary = real
        .quadratic_velocity_source(&[(FACE, MIDSIDE)], &mut || false)
        .unwrap();
    assert_eq!(
        volume.integrated_shape(),
        original_volume.integrated_shape()
    );
    assert_eq!(
        boundary.integrated_shape_m2(),
        original_boundary.integrated_shape_m2()
    );
    let q = C64::new(0.02, -0.007);
    let u = C64::new(-0.03, 0.01);
    let field = equivalent
        .solve_with_sources(&[], &[(&boundary, u)], &[(&volume, q)], &[], &mut || false)
        .unwrap();
    let reference = original
        .solve_with_boundary_sources(
            &[],
            &[(&original_boundary, u)],
            &[(&original_volume, q)],
            &mut || false,
        )
        .unwrap();
    for (&a, &b) in field.pressure_pa.iter().zip(&reference.pressure_pa) {
        close(a, b, 3e-12);
    }
    relative(field.input_power_w, reference.input_power_w, 3e-11);
    relative(
        field.boundary_source_power_w,
        reference.boundary_source_power_w,
        3e-11,
    );
    relative(
        field.volume_source_power_w,
        reference.volume_source_power_w,
        3e-11,
    );
    assert_eq!(field.density_dissipation_w, 0.0);
    assert_eq!(field.compressional_dissipation_w, 0.0);
    balance(&field);
    // The shared type never weakens the original source identity contract.
    assert!(
        equivalent
            .solve_with_sources(&[], &[], &[(&original_volume, C64::ZERO)], &[], &mut || {
                false
            })
            .is_err()
    );
    assert!(
        equivalent
            .solve_with_sources(
                &[],
                &[(&original_boundary, C64::ZERO)],
                &[],
                &[],
                &mut || false
            )
            .is_err()
    );
    assert!(
        original
            .solve_with_volume_sources(&[], &[(&volume, C64::ZERO)], &mut || false)
            .is_err()
    );
    assert!(
        original
            .solve_with_boundary_sources(&[], &[(&boundary, C64::ZERO)], &[], &mut || false)
            .is_err()
    );
}

#[test]
fn source_identity_support_budgets_and_cancellation_preserve_the_original_owner() {
    let (mesh, positions) = single_tet();
    let geometry = HelmholtzFe::new(&mesh, &positions, 2, MEDIUM).unwrap();
    let prepared = prepare(&geometry, MATERIAL, &[], &[], &[]);
    let another = prepare(&geometry, MATERIAL, &[], &[], &[]);
    let volume = prepared
        .volume_velocity_source(
            4,
            &mut |_, _| Ok::<_, HeterogeneousError>(1.0 / VOLUME),
            &mut || false,
        )
        .unwrap();
    let boundary = prepared
        .quadratic_velocity_source(&[(FACE, MIDSIDE)], &mut || false)
        .unwrap();
    assert!(
        another
            .solve_with_sources(&[], &[], &[(&volume, C64::ZERO)], &[], &mut || false)
            .is_err()
    );
    assert!(
        another
            .solve_with_sources(&[], &[(&boundary, C64::ZERO)], &[], &[], &mut || false)
            .is_err()
    );
    assert!(
        prepared
            .solve_with_sources(
                &[],
                &[],
                &[(&volume, C64::new(f64::NAN, 0.0))],
                &[],
                &mut || false
            )
            .is_err()
    );
    assert!(
        prepared
            .quadratic_velocity_source(&[(FACE, [f64::INFINITY; 6])], &mut || false)
            .is_err()
    );
    let calls = Cell::new(0);
    assert!(
        prepared
            .volume_velocity_source(
                3,
                &mut |_, _| {
                    calls.set(calls.get() + 1);
                    Ok::<_, HeterogeneousError>(1.0)
                },
                &mut || false
            )
            .is_err()
    );
    assert_eq!(calls.get(), 0);
    assert!(matches!(
        prepared.volume_velocity_source(
            4,
            &mut |_, _| Err::<f64, _>(HeterogeneousError::Cancelled),
            &mut || false
        ),
        Err(HeterogeneousError::Cancelled)
    ));
    let release = prepare(&geometry, MATERIAL, &[], &[], &[FACE]);
    let released_source = release
        .quadratic_velocity_source(&[(FACE, [0.0; 6])], &mut || false)
        .unwrap();
    assert!(
        release
            .solve_with_sources(&[], &[(&released_source, C64::ZERO)], &[], &[], &mut || {
                false
            })
            .is_err()
    );
    let at_cap = vec![(&volume, C64::ZERO); MAX_VOLUME_SOURCES];
    assert!(
        prepared
            .solve_with_sources(&[], &[], &at_cap, &[], &mut || false)
            .is_ok()
    );
    assert!(matches!(
        prepared.solve_with_sources(
            &[],
            &[],
            &vec![(&volume, C64::ZERO); MAX_VOLUME_SOURCES + 1],
            &[],
            &mut || false
        ),
        Err(HeterogeneousError::Limit { .. })
    ));
    assert!(matches!(
        prepared.solve_with_sources(
            &[],
            &vec![(&boundary, C64::ZERO); MAX_BOUNDARY_SOURCES + 1],
            &[],
            &[],
            &mut || false
        ),
        Err(HeterogeneousError::Limit { .. })
    ));
    let q = C64::new(0.02, -0.007);
    let expected = prepared
        .solve_with_sources(&[], &[], &[(&volume, q)], &[], &mut || false)
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        prepared.solve_with_sources(&[], &[], &[(&volume, q)], &[], &mut || {
            polls += 1;
            polls == 12
        }),
        Err(HeterogeneousError::Cancelled)
    ));
    let again = prepared
        .solve_with_sources(&[], &[], &[(&volume, q)], &[], &mut || false)
        .unwrap();
    assert_eq!(again.pressure_pa, expected.pressure_pa);
    assert_eq!(again.volume_source_power_w, expected.volume_source_power_w);
}
