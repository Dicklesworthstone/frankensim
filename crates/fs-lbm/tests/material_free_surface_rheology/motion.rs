//! Moving-wall work and source rheology share the original native collision.

use super::*;
use fs_lbm::d3q19::freesurface3::{PlanarWall3, PlanarWallFieldControls3};
use fs_lbm::material::free_surface::MaterialPlanarWallVelocityField3;

fn controls() -> PlanarWallFieldControls3 {
    PlanarWallFieldControls3 {
        maximum_domain_cells: 512,
        maximum_wall_faces: 64,
        maximum_relative_mass_defect: 1e-11,
        maximum_relative_face_mass_defect: 1e-12,
    }
}

fn resting_cells() -> Vec<MaterialFreeSurfaceInitialCell3> {
    initial(true)
        .into_iter()
        .map(|mut cell| {
            cell.velocity_m_s = [0.0; 3];
            cell
        })
        .collect()
}

fn velocities() -> Vec<[f64; 2]> {
    (0..64)
        .map(|i| {
            [
                0.02 + 0.002 * ((i % 8) as f64 - 3.5),
                -0.001 * ((i / 8) as f64 - 3.5),
            ]
        })
        .collect()
}

fn motion(velocities: &[[f64; 2]]) -> MaterialPlanarWallVelocityField3<'_> {
    MaterialPlanarWallVelocityField3 {
        wall: PlanarWall3::ZMin,
        moment_origin_m: [-0.003, 0.009, 0.002],
        tangential_velocity_m_s: velocities,
    }
}

fn flow<'s, 'm>(
    binding: &'s IsothermalRheologyLaw<'m>,
    cells: &[MaterialFreeSurfaceInitialCell3],
    gravity: [f64; 3],
    cx: &Cx<'_>,
    gate: &CancelGate,
) -> MaterialFreeSurfaceRheology3<'s, 'm> {
    MaterialFreeSurfaceRheology3::from_cells_with_pressure_reference(
        binding,
        DX,
        DT,
        0.1,
        domain(gravity),
        cells,
        solve(),
        100_000.0,
        cx,
        gate,
    )
    .unwrap()
}

fn same_state(
    actual: &MaterialFreeSurfaceRheology3<'_, '_>,
    expected: &MaterialFreeSurfaceRheology3<'_, '_>,
) {
    assert_eq!(
        actual.lattice().grid().state_bits(),
        expected.lattice().grid().state_bits()
    );
    assert_eq!(
        actual.ledger_mass_kg().unwrap(),
        expected.ledger_mass_kg().unwrap()
    );
    assert_eq!(actual.lattice().carry(), expected.lattice().carry());
    assert_eq!(
        actual.lattice().conversions(),
        expected.lattice().conversions()
    );
    for cell in resting_cells() {
        assert_eq!(
            actual
                .sample_at_si(
                    cell.coordinates[0],
                    cell.coordinates[1],
                    cell.coordinates[2]
                )
                .unwrap(),
            expected
                .sample_at_si(
                    cell.coordinates[0],
                    cell.coordinates[1],
                    cell.coordinates[2]
                )
                .unwrap()
        );
        assert_eq!(
            actual.material_at(cell.coordinates),
            expected.material_at(cell.coordinates)
        );
        assert_eq!(
            actual.rheology_at(cell.coordinates),
            expected.rheology_at(cell.coordinates)
        );
    }
}

#[test]
fn constant_source_moving_wall_keeps_native_population_and_work_bits_without_caloric_data() {
    let gate = CancelGate::new_clock_free();
    let source = card(0.0, 1000.0, 0.001);
    let binding = law(&source);
    let frozen = ResolvedLiquidSurfaceState::try_from_material(
        bundle(&source),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
    )
    .unwrap();
    let gravity = [0.03, -0.02, 0.01];
    let velocities = velocities();
    with_cx(&gate, |cx| {
        for filled in [false, true] {
            let mut cells = resting_cells();
            if !filled {
                cells.retain(|cell| cell.coordinates[2] < 6);
            }
            let mut actual = flow(&binding, &cells, gravity, cx, &gate);
            let mut expected = plan_free_surface_from_material(&frozen, DX, DT, 0.1)
                .unwrap()
                .with_pressure_reference(100_000.0)
                .unwrap()
                .free_surface3_from_cells(domain(gravity), &cells, cx)
                .unwrap();
            for _ in 0..4 {
                let report = actual
                    .step_with_planar_wall_velocity_field_si(motion(&velocities), controls(), &gate)
                    .unwrap();
                let native_report = expected
                    .step_with_planar_wall_velocity_field_si(motion(&velocities), controls(), &gate)
                    .unwrap();
                assert_eq!(report, native_report);
                assert_eq!(
                    actual.lattice().grid().state_bits(),
                    expected.lattice().grid().state_bits()
                );
                assert_eq!(
                    actual.ledger_mass_kg().unwrap(),
                    expected.ledger_mass_kg().unwrap()
                );
                let material = actual.material_at([3, 3, 0]).unwrap();
                assert!(material.property("specific_heat_capacity").is_none());
                assert!(material.property("thermal_conductivity").is_none());
                assert_eq!(material.card_identity(), source.content_hash());
            }
        }
    });
}

#[test]
fn wall_generated_pressure_and_shear_resolve_one_law_and_change_actual_work() {
    let gate = CancelGate::new_clock_free();
    let cells = resting_cells();
    let velocities = velocities();
    let gravity = [0.03, -0.02, 0.01];
    let pressure_slope = 0.002;
    let mut results = Vec::new();
    with_cx(&gate, |cx| {
        for rate_slope in [-0.00005, 0.0, 0.001] {
            let source = pressure::pressure_card(pressure_slope, rate_slope, 5.0);
            let binding = pressure::binding(&source, 100_000.0);
            let mut liquid = flow(&binding, &cells, gravity, cx, &gate);
            let initial_mass = liquid.ledger_mass_kg().unwrap();
            let mut work = 0.0;
            let mut pressure_change = 0.0_f64;
            for _ in 0..6 {
                let report = liquid
                    .step_with_planar_wall_velocity_field_si(motion(&velocities), controls(), &gate)
                    .unwrap();
                let independent_work: f64 = report
                    .faces
                    .iter()
                    .zip(&velocities)
                    .map(|(face, u)| {
                        let impulse = face.exchange.fluid_impulse_kg_m_s;
                        assert_eq!(face.exchange.wall_work_j, -face.exchange.fluid_work_j);
                        for d in 0..3 {
                            assert_eq!(face.exchange.wall_impulse_kg_m_s[d], -impulse[d]);
                        }
                        impulse[0] * u[0] + impulse[1] * u[1]
                    })
                    .sum();
                near(report.exchange.fluid_work_j, independent_work, 2e-18);
                assert!(report.exchange.fluid_work_j > 0.0);
                work += report.exchange.fluid_work_j;
                near(report.mass_after_kg, initial_mass, 2e-14);
                for coordinates in [[2, 3, 0], [5, 4, 0], [3, 3, 2]] {
                    let [x, y, z] = coordinates;
                    let native = liquid.lattice().sample_at(x, y, z).unwrap();
                    // Independent SI scaling, rho0=1000 and dx/dt=1.
                    let p = 100_000.0 + (native.density - 1.0) * 1000.0 / 3.0;
                    pressure_change = pressure_change.max((p - 100_000.0).abs());
                    let a = 0.5 + 3.0 * (0.1 + pressure_slope * (p - 100_000.0));
                    let rhs = numerator(&liquid, gravity, coordinates) / DT;
                    // Solve gamma * (a + 3*slope*gamma) = actual Pi-norm/dt.
                    let rate = 2.0 * rhs / (a + (a * a + 12.0 * rate_slope * rhs).sqrt());
                    let mu = 0.1 + pressure_slope * (p - 100_000.0) + rate_slope * rate;
                    let material = liquid.material_at(coordinates).unwrap();
                    let sample = liquid.rheology_at(coordinates).unwrap();
                    let source_pressure = material
                        .query_point()
                        .iter()
                        .find(|(name, _)| name == "pressure")
                        .unwrap()
                        .1;
                    near(source_pressure, p, 3e-11);
                    near(sample.shear_rate_per_s, rate, 1e-8);
                    near(
                        material.property("dynamic_viscosity").unwrap().value_si(),
                        mu,
                        1e-11,
                    );
                    near(sample.collision_tau, 0.5 + 3.0 * mu, 1e-10);
                    assert_eq!(material.card_identity(), source.content_hash());
                }
            }
            assert!(
                pressure_change > 1e-6,
                "actual wall-driven pressure must be sampled"
            );
            results.push((liquid.lattice().grid().state_bits(), work.to_bits()));
        }
    });
    for i in [0, 2] {
        assert_ne!(
            results[i].0, results[1].0,
            "source shear response must reach collision"
        );
        assert_ne!(
            results[i].1, results[1].1,
            "source shear response must change actuator work"
        );
    }
}

#[test]
fn late_wall_shear_refusal_preserves_all_material_state_and_a_safe_retry() {
    let gate = CancelGate::new_clock_free();
    let cells = resting_cells();
    let velocities = velocities();
    let source = card(0.003, 1e-6, 0.001);
    let binding = law(&source);
    with_cx(&gate, |cx| {
        let mut liquid = flow(&binding, &cells, [0.0; 3], cx, &gate);
        let mut expected = liquid.clone();
        // Construction admits the resting source. Actual wall reflection
        // creates out-of-domain strain on the still-private candidate.
        assert!(matches!(
            liquid.step_with_planar_wall_velocity_field_si(motion(&velocities), controls(), &gate),
            Err(FreeSurfaceRheologyError::Material(_))
        ));
        same_state(&liquid, &expected);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            liquid
                .step_with_planar_wall_velocity_field_si(
                    motion(&velocities),
                    controls(),
                    &cancelled
                )
                .is_err()
        );
        same_state(&liquid, &expected);
        let resting = [[0.0; 2]; 64];
        let retry = liquid
            .step_with_planar_wall_velocity_field_si(motion(&resting), controls(), &gate)
            .unwrap();
        let direct = expected
            .step_with_planar_wall_velocity_field_si(motion(&resting), controls(), &gate)
            .unwrap();
        assert_eq!(retry, direct);
        same_state(&liquid, &expected);
    });
}

#[test]
fn pressure_refusal_preserves_the_last_accepted_moving_wall_state() {
    let gate = CancelGate::new_clock_free();
    let cells = resting_cells();
    let velocities = velocities();
    let source = pressure::pressure_card(0.002, 0.001, 1e-6);
    let binding = pressure::binding(&source, 100_000.0);
    with_cx(&gate, |cx| {
        let mut liquid = flow(&binding, &cells, [0.0; 3], cx, &gate);
        let resting = liquid.lattice().grid().state_bits();
        // At rest, opposite tangential reflection corrections cancel density
        // at EACH receiving face. First accept the generated shear, then test
        // pressure created when that nonuniform state streams on the next tick.
        let first = liquid
            .step_with_planar_wall_velocity_field_si(motion(&velocities), controls(), &gate)
            .unwrap();
        assert!(first.exchange.fluid_work_j > 0.0);
        assert_ne!(liquid.lattice().grid().state_bits(), resting);
        for cell in &cells {
            let [x, y, z] = cell.coordinates;
            let native = liquid.lattice().sample_at(x, y, z).unwrap();
            assert!(((native.density - 1.0) * 1000.0 / 3.0).abs() < 1e-6);
        }
        let accepted = liquid.clone();
        assert!(matches!(
            liquid.step_with_planar_wall_velocity_field_si(motion(&velocities), controls(), &gate),
            Err(FreeSurfaceRheologyError::Material(_))
        ));
        same_state(&liquid, &accepted);
    });
}
