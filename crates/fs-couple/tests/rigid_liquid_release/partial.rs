//! Selected and retained matter must together reproduce the original body.

use super::*;
use fs_couple::liquid_transfer::rigid::RigidLiquidReleaseReport;
use fs_couple::liquid_transfer::rigid::partial::{
    PartialRigidLiquidReleaseControls, PartialRigidLiquidReleaseError, release_rigid_liquid_cells,
};
use fs_lbm::material::free_surface::enthalpy::MaterialFreeSurfaceEnthalpy3;
use fs_lbm::material::free_surface::transfer::LiquidTransferMoments;
use fs_mbd::{Gravity, RigidBodyIntegrator, Wrench};

fn partial_controls() -> PartialRigidLiquidReleaseControls {
    let mut release = controls();
    release.transfer.mechanical.max_tetrahedra = 24;
    PartialRigidLiquidReleaseControls {
        release,
        principal_relative_tolerance: 1e-12,
        minimum_retained_clearance_m: 0.001,
    }
}

fn parcel_center(i: usize) -> GaVec3 {
    GaVec3::new(0.024 * i as f64, 0.0, 0.0)
}
fn parcel_lengths(i: usize) -> GaVec3 {
    [
        GaVec3::new(0.004, 0.003, 0.002),
        GaVec3::new(0.003, 0.003, 0.002),
        GaVec3::new(0.002, 0.002, 0.002),
    ][i]
}
fn authored_rotation() -> So3 {
    rotate(GaVec3::new(0.3, -0.2, 0.1))
}
fn parcels(card: &MaterialCard, count: usize) -> (TetrahedralDensityBody, Vec<f64>) {
    let mut p = Vec::new();
    let mut topology = Vec::new();
    let mut h = Vec::new();
    for i in 0..count {
        let center = parcel_center(i);
        let vertices = points(center, parcel_lengths(i), authored_rotation());
        h.extend(
            vertices
                .iter()
                .map(|p| 180000.0 + 20000.0 * i as f64 + 2e5 * (p.x - center.x) - 1e5 * p.y),
        );
        p.extend(vertices);
        topology.extend(TOPOLOGY.map(|tet| tet.map(|v| v + 8 * i)));
    }
    (geometry(card, &p, &topology, &vec![400.0; p.len()]), h)
}

// Independent rectangular-volume integrals. No returned partition frame or
// inertia tensor is used to construct these expected momenta and energies.
fn analytic_parcel(
    i: usize,
    offset: GaVec3,
    rotation: So3,
    parent_center: GaVec3,
    parent_velocity: GaVec3,
    omega: GaVec3,
    origin: [f64; 3],
) -> LiquidTransferMoments {
    let lengths = parcel_lengths(i);
    let mass = RHO * lengths.x * lengths.y * lengths.z;
    let center = offset + rotation.rotate(parcel_center(i)).unwrap();
    let velocity = parent_velocity + omega.cross(center - parent_center);
    let moments = GaVec3::new(
        mass * (lengths.y.powi(2) + lengths.z.powi(2)) / 12.0,
        mass * (lengths.x.powi(2) + lengths.z.powi(2)) / 12.0,
        mass * (lengths.x.powi(2) + lengths.y.powi(2)) / 12.0,
    );
    let local_omega = authored_rotation()
        .inverse()
        .rotate(rotation.inverse().rotate(omega).unwrap())
        .unwrap();
    let spin = rotation
        .rotate(
            authored_rotation()
                .rotate(GaVec3::new(
                    moments.x * local_omega.x,
                    moments.y * local_omega.y,
                    moments.z * local_omega.z,
                ))
                .unwrap(),
        )
        .unwrap();
    let momentum = velocity.scale(mass);
    let relative = center - vector(origin);
    LiquidTransferMoments {
        mass_kg: mass,
        first_mass_moment_kg_m: array(relative.scale(mass)),
        momentum_kg_m_s: array(momentum),
        angular_momentum_kg_m2_s: array(relative.cross(momentum) + spin),
        kinetic_energy_j: 0.5 * mass * velocity.dot(velocity) + 0.5 * omega.dot(spin),
    }
}

fn compare_moments(actual: LiquidTransferMoments, expected: LiquidTransferMoments) {
    close(actual.mass_kg, expected.mass_kg, 2e-18);
    close_vector(
        actual.first_mass_moment_kg_m,
        vector(expected.first_mass_moment_kg_m),
        2e-19,
    );
    close_vector(
        actual.momentum_kg_m_s,
        vector(expected.momentum_kg_m_s),
        2e-18,
    );
    close_vector(
        actual.angular_momentum_kg_m2_s,
        vector(expected.angular_momentum_kg_m2_s),
        2e-20,
    );
    close(actual.kinetic_energy_j, expected.kinetic_energy_j, 2e-20);
}

fn live_liquid_moments(
    liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    report: &RigidLiquidReleaseReport,
) -> (LiquidTransferMoments, f64) {
    let mut totals = LiquidTransferMoments::default();
    let mut h = 0.0;
    for component in &report.transfer.mechanical.components {
        for &coordinates in &component.target_cells {
            let [x, y, z] = coordinates;
            let sample = liquid.flow().sample_at_si(x, y, z).unwrap().unwrap();
            let mass = sample.density_kg_m3 * sample.fill_fraction * DX.powi(3);
            let velocity = vector(sample.velocity_m_s);
            let momentum = velocity.scale(mass);
            let r = vector(coordinates.map(|v| (v as f64 + 0.5) * DX));
            totals.mass_kg += mass;
            totals.first_mass_moment_kg_m =
                array(vector(totals.first_mass_moment_kg_m) + r.scale(mass));
            totals.momentum_kg_m_s = array(vector(totals.momentum_kg_m_s) + momentum);
            totals.angular_momentum_kg_m2_s =
                array(vector(totals.angular_momentum_kg_m2_s) + r.cross(momentum));
            totals.kinetic_energy_j += 0.5 * mass * velocity.dot(velocity);
            h += mass * liquid.specific_enthalpy_at(coordinates).unwrap();
        }
    }
    (totals, h)
}

#[test]
fn g1_partial_release_retains_full_inertia_motion_heat_and_native_continuation() {
    let card = card("partial native rigid liquid", [RHO; 2]);
    let query = point(400.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &query,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let surface = resolve_liquid_surface_state(
        &card,
        &query,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let caloric = chart(&card, SolidLiquidPhase::Liquid);
    let (body, h) = parcels(&card, 2);
    let frame = body.principal_frame(1e-12, || false).unwrap();
    let untouched = (body.clone(), h.clone());
    let rotation = rotate(GaVec3::new(-0.2, 0.4, 0.1));
    let omega = GaVec3::new(0.3, -0.4, 0.5);
    let gate = CancelGate::new_clock_free();
    for (origin, boost) in [
        ([0.0; 3], [0.0; 3]),
        ([-0.02, 0.03, -0.01], [0.01, -0.005, 0.002]),
    ] {
        let offset = vector(origin) + GaVec3::new(0.008, 0.008, 0.008);
        let parent_center = offset + rotation.rotate(body.inertia().center_of_mass()).unwrap();
        let velocity = GaVec3::new(0.005, -0.004, 0.002) + vector(boost);
        let motion = motion(&frame, parent_center, rotation, velocity, omega);
        let source = state(&body, &frame, motion, &h);
        with_cx(&gate, |cx| {
            let scaling = plan_free_surface_from_material(&surface, DX, DT, 0.1).unwrap();
            let (mut liquid, mut retained, report) = release_rigid_liquid_cells(
                source,
                &[5, 2, 4, 0, 3, 1],
                scaling,
                domain([0.0; 3]),
                origin,
                &law,
                &caloric,
                partial_controls(),
                cx,
                &gate,
            )
            .unwrap();
            let selected_expected =
                analytic_parcel(0, offset, rotation, parent_center, velocity, omega, origin);
            let retained_expected =
                analytic_parcel(1, offset, rotation, parent_center, velocity, omega, origin);
            compare_moments(report.released.rigid, selected_expected);
            compare_moments(report.retained, retained_expected);
            assert_eq!(
                report.selected_lineage.original_cell_indices,
                (0..6).collect::<Vec<_>>()
            );
            assert_eq!(
                report.retained_lineage.original_cell_indices,
                (6..12).collect::<Vec<_>>()
            );
            assert_eq!(
                report.retained_lineage.original_vertex_indices,
                (8..16).collect::<Vec<_>>()
            );
            assert_eq!(retained.geometry.cells(), &body.cells()[6..]);
            assert_eq!(retained.frame.inertia(), retained.geometry.inertia());
            assert!(retained.frame.inertia().inertia_com().m[1].abs() > 1e-14);
            assert!(report.source_queries > report.released.transfer.source_queries);
            assert!(
                report.source_queries <= partial_controls().release.transfer.max_source_queries
            );
            for (local, &original) in report
                .retained_lineage
                .original_vertex_indices
                .iter()
                .enumerate()
            {
                assert_eq!(
                    retained.geometry.vertices()[local],
                    body.vertices()[original]
                );
                assert_eq!(
                    retained.specific_enthalpy_j_kg[local].to_bits(),
                    h[original].to_bits()
                );
                let before = motion
                    .point_kinematics(
                        frame.properties(),
                        frame
                            .arm_for_body_point(body.vertices()[original].position_m())
                            .unwrap(),
                    )
                    .unwrap();
                let after = retained
                    .motion
                    .point_kinematics(
                        retained.frame.properties(),
                        retained
                            .frame
                            .arm_for_body_point(retained.geometry.vertices()[local].position_m())
                            .unwrap(),
                    )
                    .unwrap();
                close_vector(
                    [
                        after.point_world.x,
                        after.point_world.y,
                        after.point_world.z,
                    ],
                    GaVec3::new(
                        before.point_world.x,
                        before.point_world.y,
                        before.point_world.z,
                    ),
                    1e-14,
                );
                close_vector(
                    [
                        after.point_velocity_world.x,
                        after.point_velocity_world.y,
                        after.point_velocity_world.z,
                    ],
                    GaVec3::new(
                        before.point_velocity_world.x,
                        before.point_velocity_world.y,
                        before.point_velocity_world.z,
                    ),
                    1e-14,
                );
            }
            let (live, live_h) = live_liquid_moments(&liquid, &report.released);
            compare_moments(live, report.released.transfer.mechanical.target);
            close(
                live.mass_kg + retained_expected.mass_kg,
                body.inertia().mass(),
                2e-18,
            );
            close_vector(
                array(vector(live.momentum_kg_m_s) + vector(retained_expected.momentum_kg_m_s)),
                vector(report.original.momentum_kg_m_s),
                2e-18,
            );
            close_vector(
                array(
                    vector(live.angular_momentum_kg_m2_s)
                        + vector(retained_expected.angular_momentum_kg_m2_s),
                ),
                vector(report.original.angular_momentum_kg_m2_s),
                2e-20,
            );
            close(
                report.retained_enthalpy_j,
                retained_expected.mass_kg * 200000.0,
                2e-12,
            );
            close(
                live_h + report.retained_enthalpy_j,
                report.original_enthalpy_j,
                2e-10,
            );
            close(report.rigid_partition_kinetic_error_j, 0.0, 2e-20);
            close(
                report.kinetic_projection_change_j,
                live.kinetic_energy_j + retained_expected.kinetic_energy_j
                    - report.original.kinetic_energy_j,
                2e-20,
            );
            // The unchanged remainder cancels through native retained
            // increments, without subtracting its absolute kinetic inventory.
            close(
                (report.retained_kinetic_projection_change_j()
                    - report.released.retained_kinetic_projection_change_j())
                .to_f64(),
                report.rigid_partition_kinetic_error_j,
                1e-30,
            );
            close(retained.time_s, source.time_s, 0.0);
            close(liquid.physical_time_s(), source.time_s, 0.0);
            let step = RigidBodyIntegrator::new(Gravity::ZERO)
                .step(
                    retained.motion,
                    retained.frame.properties(),
                    Wrench::ZERO,
                    DT,
                )
                .unwrap();
            retained.motion = step.state_after;
            retained.time_s += step.duration_seconds;
            close(
                step.diagnostics_after.mechanical_energy,
                step.diagnostics_before.mechanical_energy,
                2e-20,
            );
            let thermal = liquid.step(&gate).unwrap();
            close(
                thermal.enthalpy_before_j,
                selected_expected.mass_kg * 180000.0,
                2e-10,
            );
            close(liquid.physical_time_s(), retained.time_s, 0.0);
            assert_eq!(source.motion, motion);
        });
    }
    assert_eq!((body, h), untouched);
}

#[test]
fn g1_successive_partial_then_full_releases_keep_original_lineage_and_world_momenta() {
    let card = card("successive original parcel releases", [RHO; 2]);
    let query = point(400.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &query,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let surface = resolve_liquid_surface_state(
        &card,
        &query,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let caloric = chart(&card, SolidLiquidPhase::Liquid);
    let (body, h) = parcels(&card, 3);
    let frame = body.principal_frame(1e-12, || false).unwrap();
    let rotation = rotate(GaVec3::new(-0.2, 0.4, 0.1));
    let offset = GaVec3::new(0.008, 0.008, 0.008);
    let parent_center = offset + rotation.rotate(body.inertia().center_of_mass()).unwrap();
    let motion = motion(
        &frame,
        parent_center,
        rotation,
        GaVec3::new(0.005, -0.004, 0.002),
        GaVec3::new(0.3, -0.4, 0.5),
    );
    let source = state(&body, &frame, motion, &h);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, DX, DT, 0.1).unwrap();
        let (liquid_a, remaining_ab, a) = release_rigid_liquid_cells(
            source,
            &[0, 1, 2, 3, 4, 5],
            scaling,
            domain([0.0; 3]),
            [0.0; 3],
            &law,
            &caloric,
            partial_controls(),
            cx,
            &gate,
        )
        .unwrap();
        let origin_b = array(rotation.rotate(parcel_center(1)).unwrap());
        let (liquid_b, remaining_c, b) = release_rigid_liquid_cells(
            remaining_ab.as_source(),
            &[0, 1, 2, 3, 4, 5],
            scaling,
            domain([0.0; 3]),
            origin_b,
            &law,
            &caloric,
            partial_controls(),
            cx,
            &gate,
        )
        .unwrap();
        let origin_c = array(rotation.rotate(parcel_center(2)).unwrap());
        let (liquid_c, c) = release_rigid_liquid(
            remaining_c.as_source(),
            scaling,
            domain([0.0; 3]),
            origin_c,
            &law,
            &caloric,
            controls(),
            cx,
            &gate,
        )
        .unwrap();
        let ancestry: Vec<_> = b
            .retained_lineage
            .original_cell_indices
            .iter()
            .map(|&i| a.retained_lineage.original_cell_indices[i])
            .collect();
        assert_eq!(ancestry, (12..18).collect::<Vec<_>>());
        assert_eq!(remaining_c.geometry.cells(), &body.cells()[12..]);
        assert_eq!(remaining_c.specific_enthalpy_j_kg, h[16..]);
        let (mut mass, mut extensive_h, mut kinetic) = (0.0, 0.0, 0.0);
        let (mut momentum, mut angular) = (GaVec3::new(0.0, 0.0, 0.0), GaVec3::new(0.0, 0.0, 0.0));
        for (liquid, release, origin) in [
            (&liquid_a, &a.released, [0.0; 3]),
            (&liquid_b, &b.released, origin_b),
            (&liquid_c, &c, origin_c),
        ] {
            let (live, live_h) = live_liquid_moments(liquid, release);
            compare_moments(live, release.transfer.mechanical.target);
            mass += live.mass_kg;
            extensive_h += live_h;
            kinetic += live.kinetic_energy_j;
            momentum = momentum + vector(live.momentum_kg_m_s);
            angular = angular
                + vector(live.angular_momentum_kg_m2_s)
                + vector(origin).cross(vector(live.momentum_kg_m_s));
            close(liquid.physical_time_s(), source.time_s, 0.0);
            assert_eq!(liquid.accepted_steps(), 0);
        }
        close(mass, a.original.mass_kg, 3e-18);
        close_vector(array(momentum), vector(a.original.momentum_kg_m_s), 3e-18);
        close_vector(
            array(angular),
            vector(a.original.angular_momentum_kg_m2_s),
            3e-20,
        );
        close(extensive_h, a.original_enthalpy_j, 3e-10);
        close(
            kinetic - a.original.kinetic_energy_j,
            (a.retained_kinetic_projection_change_j()
                + b.retained_kinetic_projection_change_j()
                + c.retained_kinetic_projection_change_j())
            .to_f64(),
            3e-20,
        );
    });
}

#[test]
fn g3_partial_release_refuses_missing_occupancy_cold_remainder_and_failed_budgets_atomically() {
    let card = card("partial release refusal source", [RHO; 2]);
    let query = point(400.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &query,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let surface = resolve_liquid_surface_state(
        &card,
        &query,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let caloric = chart(&card, SolidLiquidPhase::Liquid);
    let (body, h) = parcels(&card, 2);
    let frame = body.principal_frame(1e-12, || false).unwrap();
    let motion = motion(
        &frame,
        GaVec3::new(0.008, 0.008, 0.008) + body.inertia().center_of_mass(),
        So3::identity(),
        GaVec3::new(0.01, 0.0, 0.0),
        GaVec3::new(0.0, 0.0, 0.5),
    );
    let source = state(&body, &frame, motion, &h);
    let untouched = (body.clone(), h.clone(), motion);
    let selection = [0, 1, 2, 3, 4, 5];
    let phase_chart = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (0.0, 300.0, 0.0),
            (100000.0, 350.0, 0.0),
            (200000.0, 400.0, 1.0),
            (400000.0, 500.0, 1.0),
        ]
        .map(
            |(specific_enthalpy_j_kg, temperature_k, liquid_mass_fraction)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg,
                temperature_k,
                liquid_mass_fraction,
                bulk_density_kg_m3: RHO,
            },
        )
        .to_vec(),
    )
    .unwrap();
    let mut cold_h = vec![300000.0; 16];
    cold_h[8..].fill(50000.0);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, DX, DT, 0.1).unwrap();
        let mut covering_domain = domain([0.0; 3]);
        covering_domain.cells[0] = 64;
        covering_domain.max_cells = 64 * 16 * 16;
        assert!(matches!(
            release_rigid_liquid_cells(
                source,
                &selection,
                scaling,
                covering_domain,
                [0.0; 3],
                &law,
                &caloric,
                partial_controls(),
                cx,
                &gate
            ),
            Err(PartialRigidLiquidReleaseError::RetainedOverlapsReceiver {
                original_cell: 6..=11
            })
        ));
        let mut clearance = partial_controls();
        clearance.minimum_retained_clearance_m = 0.025;
        assert!(matches!(
            release_rigid_liquid_cells(
                source,
                &selection,
                scaling,
                domain([0.0; 3]),
                [0.0; 3],
                &law,
                &caloric,
                clearance,
                cx,
                &gate
            ),
            Err(PartialRigidLiquidReleaseError::RetainedOverlapsReceiver { .. })
        ));
        for invalid in [vec![], vec![0, 0], vec![12], (0..12).collect()] {
            assert!(matches!(
                release_rigid_liquid_cells(
                    source,
                    &invalid,
                    scaling,
                    domain([0.0; 3]),
                    [0.0; 3],
                    &law,
                    &caloric,
                    partial_controls(),
                    cx,
                    &gate
                ),
                Err(PartialRigidLiquidReleaseError::Geometry(_))
            ));
        }
        assert!(matches!(
            release_rigid_liquid_cells(
                state(&body, &frame, motion, &cold_h),
                &selection,
                scaling,
                domain([0.0; 3]),
                [0.0; 3],
                &law,
                &phase_chart,
                partial_controls(),
                cx,
                &gate
            ),
            Err(PartialRigidLiquidReleaseError::Invalid(_))
        ));
        for bad in [
            {
                let mut c = partial_controls();
                c.release.transfer.mechanical.max_overlap_tests = 1;
                c
            },
            {
                let mut c = partial_controls();
                c.release.transfer.max_source_queries = 16;
                c
            },
            {
                let mut c = partial_controls();
                c.release
                    .transfer
                    .mechanical
                    .maximum_projection_energy_change_j = 0.0;
                c
            },
        ] {
            assert!(
                release_rigid_liquid_cells(
                    source,
                    &selection,
                    scaling,
                    domain([0.0; 3]),
                    [0.0; 3],
                    &law,
                    &caloric,
                    bad,
                    cx,
                    &gate
                )
                .is_err()
            );
        }
        let (_, a, report_a) = release_rigid_liquid_cells(
            source,
            &selection,
            scaling,
            domain([0.0; 3]),
            [0.0; 3],
            &law,
            &caloric,
            partial_controls(),
            cx,
            &gate,
        )
        .unwrap();
        let (_, b, report_b) = release_rigid_liquid_cells(
            source,
            &[5, 4, 3, 2, 1, 0],
            scaling,
            domain([0.0; 3]),
            [0.0; 3],
            &law,
            &caloric,
            partial_controls(),
            cx,
            &gate,
        )
        .unwrap();
        assert_eq!(a, b);
        assert_eq!(report_a, report_b);
    });
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    with_cx(&cancelled, |cx| {
        let scaling = plan_free_surface_from_material(&surface, DX, DT, 0.1).unwrap();
        assert!(matches!(
            release_rigid_liquid_cells(
                source,
                &selection,
                scaling,
                domain([0.0; 3]),
                [0.0; 3],
                &law,
                &caloric,
                partial_controls(),
                cx,
                &cancelled
            ),
            Err(PartialRigidLiquidReleaseError::Release(
                RigidLiquidReleaseError::Cancelled
            ))
        ));
    });
    assert_eq!((body, h, motion), untouched);
}
