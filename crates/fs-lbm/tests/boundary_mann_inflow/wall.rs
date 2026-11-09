//! The actual D3Q19 profile consumes the original harmonic wall observation.

use super::*;
use fs_atmo::mann_field::wall::MannPlaneWall;

fn settings(height: f64) -> BoundaryMannConfig3 {
    let mut config = config(Face3::XMin);
    config.grid_origin_in_field_m[2] = height;
    config
}

fn frozen_wall(source: &HomogeneousMannField) -> BoundaryMannSource3<'_> {
    BoundaryMannSource3::Frozen {
        source,
        time_s: 0.03,
        reference_time_s: 0.01,
        convection_velocity_m_s: [0.4, 0.1, 0.0],
    }
}

fn location(node: [usize; 3], height: f64, elapsed: f64) -> [f64; 3] {
    [
        0.2 + node[1] as f64 * 0.25 - 0.1 * elapsed,
        -0.3 - node[0] as f64 * 0.25 + 0.4 * elapsed,
        height + node[2] as f64 * 0.25,
    ]
}

#[test]
fn points_retain_no_penetration_and_solenoidal_source_in_the_native_sgs_step() {
    let gate = CancelGate::new_clock_free();
    let source = snapshot(41);
    let unchanged = source.clone();
    let config = settings(0.0);
    let wall = MannPlaneWall::new(0.0).unwrap();
    let inlet = BoundaryMannInflow3::new_with_plane_wall(
        source.model(),
        config,
        BoundaryMannObservation3::Point,
        wall,
        0,
        &gate,
    )
    .unwrap();
    assert!(inlet.spatial_filter().is_none());
    assert_eq!(inlet.wall_spatial_filter().unwrap().wall(), wall);
    assert_eq!(inlet.workspace_bytes(), 0);
    let mut domain = grid(config.face, false);
    inlet
        .apply(
            &mut domain,
            frozen_wall(&source),
            InletFlux3::Unconstrained,
            &gate,
        )
        .unwrap();
    let profile = domain.velocity_profile(config.face).unwrap().to_vec();
    for (node, &actual) in nodes(config.face).zip(&profile) {
        let expected = source
            .sample_wall_blocked(location(node, 0.0, 0.02), wall)
            .unwrap();
        near(actual, velocity(config.mean_velocity_m_s, expected.u));
        let trace = expected.grad[0][0] + expected.grad[1][1] + expected.grad[2][2];
        assert!(trace.abs() < 2e-14);
        if node[2] == 0 {
            assert_eq!(actual[2], 0.0);
            assert_eq!(expected.grad[2][0], 0.0);
            assert_eq!(expected.grad[2][1], 0.0);
        }
    }
    let mut workspace = domain.sgs_workspace().unwrap();
    let closure = SgsConfig {
        model: SgsModel::Wale { cw: 0.5 },
        filter_width_lu: 1.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    };
    for _ in 0..2 {
        domain
            .step_with_sgs(&gate, closure, &mut workspace, None)
            .unwrap();
    }
    near(domain.velocity(0, 4, 3), profile[3 * 8 + 4]);
    assert_eq!(source, unchanged);
}

#[test]
fn finite_faces_integrate_the_evanescent_wall_correction_without_restoring_energy() {
    let gate = CancelGate::new_clock_free();
    let source = snapshot(41);
    let config = settings(0.125);
    let wall = MannPlaneWall::new(0.0).unwrap();
    let inlet = BoundaryMannInflow3::new_with_plane_wall(
        source.model(),
        config,
        BoundaryMannObservation3::FaceAverage,
        wall,
        16 * 32,
        &gate,
    )
    .unwrap();
    assert_eq!(inlet.workspace_bytes(), 16 * 32);
    assert_eq!(
        inlet.wall_spatial_filter().unwrap().vertical_radius_m(),
        0.125
    );
    let mut domain = grid(config.face, false);
    inlet
        .apply(
            &mut domain,
            frozen_wall(&source),
            InletFlux3::Unconstrained,
            &gate,
        )
        .unwrap();
    // Independent tensor-product spatial quadrature of point observations.
    // The lowest face touches the wall and retains its complete finite area.
    let rule = [
        (-0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
        (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
    ];
    let bulk = MannBoxFilter::prepare(
        source.model(),
        [[0.0; 3], [0.125, 0.0, 0.0], [0.0, 0.0, 0.125]],
        16 * 8,
        || false,
    )
    .unwrap();
    let mut blocking = 0.0_f64;
    for (node, &actual) in nodes(config.face).zip(domain.velocity_profile(config.face).unwrap()) {
        let center = location(node, 0.125, 0.02);
        let mut mean = [0.0; 3];
        for (a, wa) in rule {
            for (b, wb) in rule {
                let sample = source
                    .sample_wall_blocked(
                        [center[0] + 0.125 * a, center[1], center[2] + 0.125 * b],
                        wall,
                    )
                    .unwrap();
                for i in 0..3 {
                    mean[i] += 0.25 * wa * wb * sample.u[i];
                }
            }
        }
        near(actual, velocity(config.mean_velocity_m_s, mean));
        let homogeneous = velocity(
            config.mean_velocity_m_s,
            source.sample_filtered(center, &bulk).unwrap().u,
        );
        blocking = blocking.max((actual[2] - homogeneous[2]).abs());
    }
    assert!(
        blocking > 1e-6,
        "a homogeneous filter cannot replace the wall correction"
    );
}

#[test]
fn wall_blocked_endpoint_and_interval_mean_keep_the_original_ou_path_and_clock() {
    let gate = CancelGate::new_clock_free();
    let config = settings(0.125);
    let wall = MannPlaneWall::new(0.0).unwrap();
    let mut source =
        FixedStepMannField::from_snapshot(snapshot(41), 83, config.dt_s, 1 << 20).unwrap();
    let mut average = MannStepAverage::prepare(&source, 1 << 20, || false).unwrap();
    source
        .advance_one_with_average(&mut average, 16 * 8, || false)
        .unwrap();
    let accepted_source = source.checkpoint();
    let accepted_average = average.clone();
    let inlet = BoundaryMannInflow3::new_with_plane_wall(
        source.model(),
        config,
        BoundaryMannObservation3::FaceAverage,
        wall,
        16 * 32,
        &gate,
    )
    .unwrap();
    let filter = inlet.wall_spatial_filter().unwrap();
    let mut domain = grid(config.face, false);
    let endpoint = inlet
        .apply(
            &mut domain,
            BoundaryMannSource3::Evolving(&source),
            InletFlux3::Unconstrained,
            &gate,
        )
        .unwrap();
    assert_eq!(endpoint.sampled_interval_s, [source.time_s(); 2]);
    for (node, &actual) in nodes(config.face).zip(domain.velocity_profile(config.face).unwrap()) {
        near(
            actual,
            velocity(
                config.mean_velocity_m_s,
                source
                    .sample_wall_blocked_filtered(location(node, 0.125, 0.0), filter)
                    .unwrap()
                    .u,
            ),
        );
    }
    let endpoint_values = domain.velocity_profile(config.face).unwrap().to_vec();
    let report = inlet
        .apply(
            &mut domain,
            BoundaryMannSource3::StepAverage {
                source: &source,
                average: &average,
            },
            InletFlux3::Unconstrained,
            &gate,
        )
        .unwrap();
    assert_eq!(report.sampled_interval_s, average.interval().unwrap());
    for (node, &actual) in nodes(config.face).zip(domain.velocity_profile(config.face).unwrap()) {
        near(
            actual,
            velocity(
                config.mean_velocity_m_s,
                average
                    .sample_wall_blocked_filtered(location(node, 0.125, 0.0), filter)
                    .unwrap()
                    .u,
            ),
        );
    }
    assert_ne!(
        domain.velocity_profile(config.face).unwrap(),
        endpoint_values
    );
    assert_eq!(source.checkpoint(), accepted_source);
    assert_eq!(average, accepted_average);
}

#[test]
fn wall_domain_tangency_source_and_budget_refusals_preserve_the_accepted_profile() {
    let gate = CancelGate::new_clock_free();
    let source = snapshot(41);
    let wall = MannPlaneWall::new(0.0).unwrap();
    let config = settings(0.125);
    let inlet = BoundaryMannInflow3::new_with_plane_wall(
        source.model(),
        config,
        BoundaryMannObservation3::FaceAverage,
        wall,
        16 * 32,
        &gate,
    )
    .unwrap();
    let mut domain = grid(config.face, false);
    inlet
        .apply(
            &mut domain,
            frozen_wall(&source),
            InletFlux3::Unconstrained,
            &gate,
        )
        .unwrap();
    let accepted = domain.velocity_profile(config.face).unwrap().to_vec();
    let populations = bits(&domain);
    // Reversed vertical axis: some valid samples are staged before a later
    // footprint crosses the plane. The complete published profile must survive.
    let mut crossing = settings(1.0);
    crossing.field_to_grid = [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]];
    let crossing = BoundaryMannInflow3::new_with_plane_wall(
        source.model(),
        crossing,
        BoundaryMannObservation3::FaceAverage,
        wall,
        16 * 32,
        &gate,
    )
    .unwrap();
    assert!(matches!(
        crossing.apply(
            &mut domain,
            frozen_wall(&source),
            InletFlux3::Unconstrained,
            &gate
        ),
        Err(BoundaryMannError3::Atmosphere(_))
    ));
    assert!(matches!(
        inlet.apply(
            &mut domain,
            frozen(&source),
            InletFlux3::Unconstrained,
            &gate
        ),
        Err(BoundaryMannError3::Input(_))
    ));
    assert!(matches!(
        inlet.apply(
            &mut domain,
            frozen_wall(&snapshot(42)),
            InletFlux3::Unconstrained,
            &gate
        ),
        Err(BoundaryMannError3::Atmosphere(_))
    ));
    assert!(
        BoundaryMannInflow3::new_with_plane_wall(
            source.model(),
            config,
            BoundaryMannObservation3::FaceAverage,
            wall,
            16 * 32 - 1,
            &gate
        )
        .is_err()
    );
    let mut normal_mean = config;
    normal_mean.mean_velocity_m_s[2] = 0.01;
    assert!(matches!(
        BoundaryMannInflow3::new_with_plane_wall(
            source.model(),
            normal_mean,
            BoundaryMannObservation3::Point,
            wall,
            0,
            &gate
        ),
        Err(BoundaryMannError3::Input(_))
    ));
    let mut normal_face = settings(0.0);
    normal_face.face = Face3::ZMin;
    let normal_face = BoundaryMannInflow3::new_with_plane_wall(
        source.model(),
        normal_face,
        BoundaryMannObservation3::Point,
        wall,
        0,
        &gate,
    )
    .unwrap();
    assert!(matches!(
        normal_face.apply(
            &mut grid(Face3::ZMin, false),
            frozen_wall(&source),
            InletFlux3::MeanNormal(0.0),
            &gate
        ),
        Err(BoundaryMannError3::Input(_))
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        inlet.apply(
            &mut domain,
            frozen_wall(&source),
            InletFlux3::Unconstrained,
            &cancelled
        ),
        Err(BoundaryMannError3::Boundary(
            BoundaryProfileError3::Cancelled
        ))
    ));
    assert_eq!(domain.velocity_profile(config.face).unwrap(), accepted);
    assert_eq!(bits(&domain), populations);
}
