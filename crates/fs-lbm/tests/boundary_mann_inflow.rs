//! The original Mann field drives the actual D3Q19 profile and SGS step.
#![cfg(feature = "turbulent-inflow")]

use fs_atmo::mann::MannParams;
use fs_atmo::mann_field::filter::MannBoxFilter;
use fs_atmo::mann_field::fixed_step::{FixedStepMannField, MannStepAverage};
use fs_atmo::mann_field::{HomogeneousMannField, MannQuadrature, MannSpectralModel};
use fs_exec::{CancelGate, TilePool};
use fs_lbm::d3q19::{
    BoundaryGrid3, BoundaryMannConfig3, BoundaryMannError3, BoundaryMannInflow3,
    BoundaryMannObservation3, BoundaryMannSource3, BoundaryProfileError3, BoundarySpec3, Face3,
    FaceBoundary3, InletFlux3,
};
use fs_lbm::sgs3::{SgsConfig, SgsModel, SgsUse};

#[path = "boundary_mann_inflow/wall.rs"]
mod wall;

fn snapshot(seed: u64) -> HomogeneousMannField {
    let model = MannSpectralModel::build(
        seed,
        MannParams {
            alpha_eps23: 0.025,
            length_m: 1.0,
            gamma: 3.9,
        },
        MannQuadrature {
            n_modes: 16,
            min_kl: 0.03,
            max_kl: 2.0,
        },
        [0.3, -0.1, 0.0],
    )
    .unwrap();
    HomogeneousMannField::at_anchor(model, 79, vec![0.02; 16]).unwrap()
}

fn config(face: Face3) -> BoundaryMannConfig3 {
    let mut mean = [0.0; 3];
    mean[face as usize / 2] = if face as usize % 2 == 0 { 0.2 } else { -0.2 };
    BoundaryMannConfig3 {
        face,
        dx_m: 0.25,
        dt_s: 0.005,
        mean_velocity_m_s: mean,
        grid_origin_in_field_m: [0.2, -0.3, 0.4],
        field_to_grid: [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        maximum_mode_samples: 64 * 16,
    }
}

fn grid(face: Face3, periodic: bool) -> BoundaryGrid3 {
    let mut faces = [if periodic {
        FaceBoundary3::Periodic
    } else {
        FaceBoundary3::stationary_wall()
    }; 6];
    let settings = config(face);
    faces[face as usize] = FaceBoundary3::Velocity {
        velocity: settings.mean_velocity_m_s.map(|u| u * 0.02),
    };
    faces[(face as usize) ^ 1] = FaceBoundary3::Pressure { density: 1.0 };
    BoundaryGrid3::new(8, 8, 8, 0.8, [0.0; 3], BoundarySpec3::new(faces))
}

fn nodes(face: Face3) -> impl Iterator<Item = [usize; 3]> {
    (0..8)
        .flat_map(|z| (0..8).flat_map(move |y| (0..8).map(move |x| [x, y, z])))
        .filter(move |p| p[face as usize / 2] == if face as usize % 2 == 0 { 0 } else { 7 })
}

fn frozen(source: &HomogeneousMannField) -> BoundaryMannSource3<'_> {
    BoundaryMannSource3::Frozen {
        source,
        time_s: 0.03,
        reference_time_s: 0.01,
        convection_velocity_m_s: [0.4, 0.1, -0.2],
    }
}

// Independently expanded 90-degree frame map, at on-site nodes. The frozen
// translation is in grid axes; evolving/averaged sources already own advection.
fn point(node: [usize; 3], elapsed: f64) -> [f64; 3] {
    [
        0.2 + node[1] as f64 * 0.25 - 0.1 * elapsed,
        -0.3 - node[0] as f64 * 0.25 + 0.4 * elapsed,
        0.4 + node[2] as f64 * 0.25 + 0.2 * elapsed,
    ]
}

fn velocity(mean: [f64; 3], fluctuation: [f64; 3]) -> [f64; 3] {
    [
        (mean[0] - fluctuation[1]) * 0.02,
        (mean[1] + fluctuation[0]) * 0.02,
        (mean[2] + fluctuation[2]) * 0.02,
    ]
}

fn near(actual: [f64; 3], expected: [f64; 3]) {
    for a in 0..3 {
        assert!(
            (actual[a] - expected[a]).abs() < 2e-14,
            "{actual:?} != {expected:?}"
        );
    }
}

fn bits(grid: &BoundaryGrid3) -> Vec<u64> {
    let mut values = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                values.extend(grid.populations(x, y, z).map(f64::to_bits));
            }
        }
    }
    values
}

#[test]
fn frozen_face_averages_drive_all_six_native_sgs_boundaries_and_workers() {
    let gate = CancelGate::new_clock_free();
    let pool = TilePool::for_host(2, 0x4d414e4e);
    let source = snapshot(41);
    let unchanged = source.clone();
    let closure = SgsConfig {
        model: SgsModel::Wale { cw: 0.5 },
        filter_width_lu: 1.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    };
    for face in Face3::ALL {
        let settings = config(face);
        let inlet = BoundaryMannInflow3::new(
            source.model(),
            settings,
            BoundaryMannObservation3::FaceAverage,
            16 * 8,
            &gate,
        )
        .unwrap();
        let mut domain = grid(face, false);
        let report = inlet
            .apply(
                &mut domain,
                frozen(&source),
                InletFlux3::Unconstrained,
                &gate,
            )
            .unwrap();
        assert_eq!(report.profile.sampled_cells, 64);
        assert_eq!(report.mode_samples, 64 * 16);
        assert_eq!(report.sampled_interval_s, [0.03; 2]);
        assert_eq!(report.spatial_phase_bound_rad, 0.5);
        assert_eq!(report.profile.normal_velocity_correction_lu, 0.0);
        let mut edges = [[0.0, -0.125, 0.0], [0.125, 0.0, 0.0], [0.0, 0.0, 0.125]];
        edges[face as usize / 2] = [0.0; 3];
        let oracle = MannBoxFilter::prepare(source.model(), edges, 128, || false).unwrap();
        let profile = domain.velocity_profile(face).unwrap();
        let mut difference_from_point = 0.0_f64;
        for (p, &actual) in nodes(face).zip(profile) {
            let location = point(p, 0.02);
            near(
                actual,
                velocity(
                    settings.mean_velocity_m_s,
                    source.sample_filtered(location, &oracle).unwrap().u,
                ),
            );
            let unfiltered = velocity(
                settings.mean_velocity_m_s,
                source.sample(location).unwrap().u,
            );
            for a in 0..3 {
                difference_from_point =
                    difference_from_point.max((actual[a] - unfiltered[a]).abs());
            }
        }
        assert!(
            difference_from_point > 1e-9,
            "face observation must retain its spectral attenuation"
        );
        let mut parallel = domain.clone();
        let mut serial_workspace = domain.sgs_workspace().unwrap();
        let mut parallel_workspace = parallel.sgs_workspace().unwrap();
        for step in 0..2 {
            let serial = domain
                .step_with_sgs(&gate, closure, &mut serial_workspace, None)
                .unwrap();
            let pooled = parallel
                .step_pooled_with_sgs(&pool, &gate, closure, &mut parallel_workspace, None)
                .unwrap();
            assert_eq!(serial, pooled);
            assert_eq!(bits(&domain), bits(&parallel));
            if step == 1 {
                assert!(serial.sgs.maximum_eddy_viscosity_lu > 0.0);
            }
        }
        let mut p = [3, 4, 3];
        p[face as usize / 2] = if face as usize % 2 == 0 { 0 } else { 7 };
        near(
            domain.velocity(p[0], p[1], p[2]),
            velocity(
                settings.mean_velocity_m_s,
                source.sample_filtered(point(p, 0.02), &oracle).unwrap().u,
            ),
        );
    }
    assert_eq!(source, unchanged);
}

#[test]
fn evolving_endpoints_and_exact_interval_means_use_the_original_source_without_advancing() {
    let gate = CancelGate::new_clock_free();
    let settings = config(Face3::XMin);
    let frozen_field = snapshot(41);
    let inlet = BoundaryMannInflow3::new(
        frozen_field.model(),
        settings,
        BoundaryMannObservation3::Point,
        0,
        &gate,
    )
    .unwrap();
    assert_eq!(inlet.spatial_filter().unwrap().workspace_bytes(), 0);
    let mut source =
        FixedStepMannField::from_snapshot(frozen_field, 83, settings.dt_s, 1 << 20).unwrap();
    let mut average = MannStepAverage::prepare(&source, 1 << 20, || false).unwrap();
    source
        .advance_one_with_average(&mut average, 16 * 8, || false)
        .unwrap();
    let accepted_source = source.checkpoint();
    let accepted_average = average.clone();
    let mut domain = grid(settings.face, false);
    let endpoint = inlet
        .apply(
            &mut domain,
            BoundaryMannSource3::Evolving(&source),
            InletFlux3::Unconstrained,
            &gate,
        )
        .unwrap();
    assert_eq!(endpoint.sampled_interval_s, [source.time_s(); 2]);
    for (p, &actual) in nodes(settings.face).zip(domain.velocity_profile(settings.face).unwrap()) {
        near(
            actual,
            velocity(
                settings.mean_velocity_m_s,
                source.sample(point(p, 0.0)).unwrap().u,
            ),
        );
    }
    let endpoint_profile = domain.velocity_profile(settings.face).unwrap().to_vec();
    let mean = inlet
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
    assert_eq!(mean.sampled_interval_s, average.interval().unwrap());
    for (p, &actual) in nodes(settings.face).zip(domain.velocity_profile(settings.face).unwrap()) {
        near(
            actual,
            velocity(
                settings.mean_velocity_m_s,
                average.sample(point(p, 0.0)).unwrap().u,
            ),
        );
    }
    assert_ne!(
        domain.velocity_profile(settings.face).unwrap(),
        endpoint_profile.as_slice()
    );
    assert_eq!(source.checkpoint(), accepted_source);
    assert_eq!(average, accepted_average);
    let accepted_profile = domain.velocity_profile(settings.face).unwrap().to_vec();
    let different_clock =
        FixedStepMannField::from_snapshot(snapshot(41), 83, 2.0 * settings.dt_s, 1 << 20).unwrap();
    assert!(matches!(
        inlet.apply(
            &mut domain,
            BoundaryMannSource3::Evolving(&different_clock),
            InletFlux3::Unconstrained,
            &gate
        ),
        Err(BoundaryMannError3::Input(_))
    ));
    source.advance_to(2, 16 * 4, || false).unwrap();
    assert!(matches!(
        inlet.apply(
            &mut domain,
            BoundaryMannSource3::StepAverage {
                source: &source,
                average: &average
            },
            InletFlux3::Unconstrained,
            &gate
        ),
        Err(BoundaryMannError3::Input(_))
    ));
    assert_eq!(
        domain.velocity_profile(settings.face).unwrap(),
        accepted_profile.as_slice()
    );
}

#[test]
fn resolution_source_topology_work_and_cancellation_refusals_preserve_the_accepted_trace() {
    let gate = CancelGate::new_clock_free();
    let source = snapshot(41);
    let settings = config(Face3::XMin);
    let inlet = BoundaryMannInflow3::new(
        source.model(),
        settings,
        BoundaryMannObservation3::Point,
        0,
        &gate,
    )
    .unwrap();
    let mut domain = grid(settings.face, false);
    inlet
        .apply(
            &mut domain,
            frozen(&source),
            InletFlux3::Unconstrained,
            &gate,
        )
        .unwrap();
    let accepted = domain.velocity_profile(settings.face).unwrap().to_vec();
    let populations = bits(&domain);
    let too_fast = BoundaryMannSource3::Frozen {
        source: &source,
        time_s: 0.03,
        reference_time_s: 0.01,
        convection_velocity_m_s: [1e6, 0.0, 0.0],
    };
    assert!(matches!(
        inlet.apply(&mut domain, too_fast, InletFlux3::Unconstrained, &gate),
        Err(BoundaryMannError3::TemporalNyquist { .. })
    ));
    let mut budget = settings;
    budget.maximum_mode_samples -= 1;
    let bounded = BoundaryMannInflow3::new(
        source.model(),
        budget,
        BoundaryMannObservation3::Point,
        0,
        &gate,
    )
    .unwrap();
    assert!(matches!(
        bounded.apply(
            &mut domain,
            frozen(&source),
            InletFlux3::Unconstrained,
            &gate
        ),
        Err(BoundaryMannError3::Work {
            required: 1024,
            maximum: 1023
        })
    ));
    let different_source = snapshot(42);
    assert!(matches!(
        inlet.apply(
            &mut domain,
            frozen(&different_source),
            InletFlux3::Unconstrained,
            &gate
        ),
        Err(BoundaryMannError3::Atmosphere(_))
    ));
    assert!(matches!(
        inlet.apply(
            &mut domain,
            frozen(&source),
            InletFlux3::MeanNormal(-0.004),
            &gate
        ),
        Err(BoundaryMannError3::Boundary(
            BoundaryProfileError3::MixedWallSeam
        ))
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        inlet.apply(
            &mut domain,
            frozen(&source),
            InletFlux3::Unconstrained,
            &cancelled
        ),
        Err(BoundaryMannError3::Boundary(
            BoundaryProfileError3::Cancelled
        ))
    ));
    assert_eq!(
        domain.velocity_profile(settings.face).unwrap(),
        accepted.as_slice()
    );
    assert_eq!(bits(&domain), populations);
    let mut periodic = grid(settings.face, true);
    assert!(matches!(
        inlet.apply(
            &mut periodic,
            frozen(&source),
            InletFlux3::Unconstrained,
            &gate
        ),
        Err(BoundaryMannError3::Input(_))
    ));
    assert!(periodic.velocity_profile(settings.face).is_none());
    let mut aliased = settings;
    aliased.dx_m = 2.0;
    assert!(matches!(
        BoundaryMannInflow3::new(
            source.model(),
            aliased,
            BoundaryMannObservation3::Point,
            0,
            &gate
        ),
        Err(BoundaryMannError3::SpatialNyquist { .. })
    ));
    let mut reflected = settings;
    reflected.field_to_grid[2][2] = -1.0;
    assert!(matches!(
        BoundaryMannInflow3::new(
            source.model(),
            reflected,
            BoundaryMannObservation3::Point,
            0,
            &gate
        ),
        Err(BoundaryMannError3::Input(_))
    ));
}
