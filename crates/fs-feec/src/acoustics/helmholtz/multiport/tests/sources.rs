//! Original pressure sources: independent compression and full-field equations.

use super::*;
use crate::acoustics::helmholtz::QuadraticGeometryControls;
use crate::acoustics::helmholtz::heterogeneous::{
    HarmonicLimits, HeterogeneousError, HeterogeneousHelmholtzFe, MAX_VOLUME_SOURCES,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use std::cell::Cell;

fn complex(value: [f64; 2]) -> C64 {
    C64::new(value[0], value[1])
}

#[test]
fn g1_curved_full_pressure_injection_retains_exact_compression_and_work() {
    let (mesh, original) = single_tet();
    let a = 0.6;
    let map = |[x, y, z]: [f64; 3]| [x, y, z * (1.0 + a * x)];
    let positions = original.iter().copied().map(map).collect::<Vec<_>>();
    let midpoints = mesh.edges.iter().map(|&[i, j]| {
        map(core::array::from_fn(|axis| {
            f64::midpoint(original[i as usize][axis], original[j as usize][axis])
        }))
    }).collect::<Vec<_>>();
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate, arena,
            StreamKey { seed: 19, kernel_id: 93, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic,
        );
        let geometry = HelmholtzFe::new_quadratic(
            &mesh, &positions, &midpoints, 2, 7, AIR,
            QuadraticGeometryControls::default(), &cx,
        ).unwrap();
        let ports = [BoundaryMode::uniform(&[[0, 1, 2]])];
        let model = InteriorMultiport::new_full(&geometry, &ports, geometry.ndof(), &|| false)
            .unwrap();
        // The source identity survives movement; it is not an old stack address.
        let geometry = Box::new(geometry);
        let calls = Cell::new(0);
        let volume = model.volume_velocity_source(
            &geometry, 7, &mut |_, _| {
                calls.set(calls.get() + 1);
                Ok::<_, HeterogeneousError>(1.0)
            }, &mut || false,
        ).unwrap();
        let boundary = model.quadratic_velocity_source_with_quadrature(
            &geometry, &[([0, 1, 3], [0.0, 0.0, 0.0, 1.0, 0.0, 0.0])],
            7, &mut || false,
        ).unwrap();
        // det(J)=1+a*x; on y=0 the same factor multiplies physical area.
        let physical_volume = 1.0 / 6.0 + a / 24.0;
        assert!((volume.integrated_shape() - physical_volume).abs() < 3e-13);
        assert!((boundary.integrated_shape_m2() - (1.0 / 6.0 + a / 15.0)).abs() < 3e-13);
        let sampled = calls.get();
        let frequency = 137.0;
        let omega = core::f64::consts::TAU * frequency;
        let bulk = C64::new(AIR.density_kg_m3 * AIR.sound_speed_m_s.powi(2), 8000.0);
        let q = C64::new(0.02, -0.007);
        // Uniform continuity i*w*p/K=q, with every physical boundary rigid.
        // The external incident wave supplies exactly P/(2*sqrt(R)), so Q_port=0.
        let exact = bulk * q / C64::new(0.0, omega);
        for resistance in [1.1e3_f64, 7.3e3] {
            let prepared = model.prepare_scattering_with_bulk_modulus(
                frequency, bulk, &[resistance], &|| false,
            ).unwrap();
            let field = prepared.response_with_sources(
                &[exact.scale(0.5 / resistance.sqrt())], &[], &[(&volume, q)], &|| false,
            ).unwrap();
            for (i, &p) in field.amplitudes.iter().enumerate() {
                let expected = if i < mesh.vertex_count { exact } else { C64::ZERO };
                assert!((complex(p) - expected).abs() < 2e-9 * exact.abs());
            }
            assert!(field.inward_flow_m3_s[0].abs() < 2e-9 * q.abs() * physical_volume);
            let work = physical_volume * bulk.im * q.norm_sq() / (2.0 * omega);
            let loss = prepared.bulk_dissipation_w(&field.amplitudes, &|| false).unwrap();
            assert!((field.volume_source_power_w - work).abs() < 2e-9 * work);
            assert!((loss - work).abs() < 2e-9 * work);
            assert!((field.input_power_w() + field.volume_source_power_w - loss).abs()
                < 2e-9 * work);
            assert_eq!(field.boundary_source_power_w, 0.0);
            let zero = prepared.response_with_sources(
                &[C64::ZERO], &[], &[(&volume, q), (&volume, -q)], &|| false,
            ).unwrap();
            assert!(zero.amplitudes.iter().all(|p| *p == [0.0; 2]));
            assert_eq!(zero.volume_source_power_w, 0.0);
        }
        assert_eq!(calls.get(), sampled, "warm solves must retain the original source projection");
    });
}

#[test]
fn g2_signed_full_pressure_sources_match_original_heterogeneous_equations() {
    let (mesh, points) = single_tet();
    let medium = Medium { density_kg_m3: 1.3, sound_speed_m_s: 2.0 };
    let geometry = HelmholtzFe::new(&mesh, &points, 2, medium).unwrap();
    let ports = [BoundaryMode::uniform(&[[1, 2, 3]])];
    let model = InteriorMultiport::new_full(&geometry, &ports, geometry.ndof(), &|| false)
        .unwrap();
    let original = HeterogeneousHelmholtzFe::new(
        &geometry, 4, &mut |_, _| Ok::<_, AcousticFeError>(medium), &mut || false,
    ).unwrap();
    let shape = [0.0, 0.0, 0.0, 1.0, -1.0, 0.0];
    let boundary = model.quadratic_velocity_source(&geometry, &[([0, 1, 2], shape)], &mut || false)
        .unwrap();
    let volume = model.volume_velocity_source(
        &geometry, 4, &mut |_, p| Ok::<_, HeterogeneousError>(p[0] - 0.25), &mut || false,
    ).unwrap();
    let original_boundary = original.quadratic_velocity_source(
        &[([0, 1, 2], shape)], &mut || false,
    ).unwrap();
    let original_volume = original.volume_velocity_source(
        4, &mut |_, p| Ok::<_, HeterogeneousError>(p[0] - 0.25), &mut || false,
    ).unwrap();
    assert!(boundary.integrated_shape_m2().abs() < 1e-14);
    assert!(boundary.has_nonzero_projection());
    assert!(volume.integrated_shape().abs() < 1e-14);
    let frequency = 0.35;
    let q = C64::new(0.02, -0.007);
    let u = C64::new(-0.03, 0.01);
    let port_flow = C64::new(0.04, 0.02);
    let prepared = model.prepare(frequency, &|| false).unwrap();
    let field = prepared.response_with_sources(
        &[port_flow], &[(&boundary, u)], &[(&volume, q)], &|| false,
    ).unwrap();
    let reference = original.prepare(
        frequency, &[], HarmonicLimits::default(), &mut || false,
    ).unwrap().solve_with_boundary_sources(
        &[([1, 2, 3], [port_flow.scale(1.0 / model.port_areas_m2()[0]); 3])],
        &[(&original_boundary, u)], &[(&original_volume, q)], &mut || false,
    ).unwrap();
    let scale = reference.pressure_pa.iter().map(|p| p.abs()).fold(1e-30_f64, f64::max);
    assert!(field.amplitudes.iter().zip(&reference.pressure_pa)
        .all(|(&p, &expected)| (complex(p) - expected).abs() < 2e-10 * scale));
    let work_scale = scale * (port_flow.abs() + 0.5 * u.abs() + q.abs() / 6.0);
    assert!((field.boundary_source_power_w - reference.boundary_source_power_w).abs()
        < 2e-10 * work_scale);
    assert!((field.volume_source_power_w - reference.volume_source_power_w).abs()
        < 2e-10 * work_scale);
    assert!((field.input_power_w() + field.boundary_source_power_w + field.volume_source_power_w)
        .abs() < 2e-10 * work_scale);
    assert!(field.amplitudes.iter().skip(mesh.vertex_count).any(|p| complex(*p).abs() > 1e-6 * scale));

    // Exact owner and support admission holds even when the requested amplitude is zero.
    let other_geometry = HelmholtzFe::new(&mesh, &points, 2, medium).unwrap();
    assert!(model.volume_velocity_source(
        &other_geometry, 4, &mut |_, _| Ok::<_, HeterogeneousError>(1.0), &mut || false,
    ).is_err());
    let other_model = InteriorMultiport::new_full(&geometry, &ports, geometry.ndof(), &|| false)
        .unwrap();
    let foreign = other_model.volume_velocity_source(
        &geometry, 4, &mut |_, _| Ok::<_, HeterogeneousError>(1.0), &mut || false,
    ).unwrap();
    assert!(prepared.response_with_sources(&[C64::ZERO], &[], &[(&foreign, C64::ZERO)], &|| false)
        .is_err());
    assert!(prepared.response_with_sources(
        &[C64::ZERO], &[], &[(&original_volume, C64::ZERO)], &|| false,
    ).is_err());
    assert!(model.quadratic_velocity_source(
        &geometry, &[([1, 2, 3], [0.0; 6])], &mut || false,
    ).is_err());
    let reduced = InteriorMultiport::new(&geometry, &complete(&geometry), &ports, false, &|| false)
        .unwrap();
    assert!(reduced.volume_velocity_source(
        &geometry, 4, &mut |_, _| Ok::<_, HeterogeneousError>(1.0), &mut || false,
    ).is_err());
    let legacy = prepared.response(&[port_flow], &|| false).unwrap();
    let empty = prepared.response_with_sources(&[port_flow], &[], &[], &|| false).unwrap();
    assert_eq!(legacy.amplitudes, empty.amplitudes);
    assert!(prepared.response_with_sources(
        &[C64::ZERO], &[], &[(&volume, C64::new(f64::NAN, 0.0))], &|| false,
    ).is_err());
    let mut at_cap = vec![(&volume, C64::ZERO); MAX_VOLUME_SOURCES];
    assert!(prepared.response_with_sources(&[C64::ZERO], &[], &at_cap, &|| false).is_ok());
    at_cap.push((&volume, C64::ZERO));
    assert!(matches!(
        prepared.response_with_sources(&[C64::ZERO], &[], &at_cap, &|| false),
        Err(HeterogeneousError::Limit { maximum: MAX_VOLUME_SOURCES, .. })
    ));
    let polls = Cell::new(0);
    assert!(prepared.response_with_sources(
        &[port_flow], &[(&boundary, u)], &[(&volume, q)], &|| {
            polls.set(polls.get() + 1);
            polls.get() >= 8
        },
    ).is_err());
    let replay = prepared.response_with_sources(
        &[port_flow], &[(&boundary, u)], &[(&volume, q)], &|| false,
    ).unwrap();
    assert_eq!(replay.amplitudes, field.amplitudes);
    assert_eq!(replay.volume_source_power_w, field.volume_source_power_w);
}
