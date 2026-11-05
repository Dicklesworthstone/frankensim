//! Local surface heat must drive resolved convection without changing the
//! uniform-wall solver or weakening its conservative, transactional checks.

use fs_lbm::ThermalLbm;
use fs_lbm::d3q19::ThermalLbm3;
use fs_lbm::thermal::boundary::{ThermalWall, ThermalWallField};
use fs_lbm::thermal::sources::{ThermalSourceError, ThermalSourceLimits};

fn limits() -> ThermalSourceLimits {
    ThermalSourceLimits {
        max_cells: 4096, temperature_range: [-2.0, 3.0], balance_tolerance: 1e-11,
    }
}

fn bits<const Q: usize>(populations: &[[f64; Q]]) -> Vec<u64> {
    populations.iter().flatten().map(|x| x.to_bits()).collect()
}

fn same2(actual: &ThermalLbm, expected: &ThermalLbm) {
    assert_eq!(bits(&actual.gpop), bits(&expected.gpop));
    assert_eq!(bits(&actual.grid.f), bits(&expected.grid.f));
    assert_eq!(bits(&actual.grid.fext), bits(&expected.grid.fext));
    assert_eq!(actual.last_wall_exchange(), expected.last_wall_exchange());
}

fn same3(actual: &ThermalLbm3, expected: &ThermalLbm3) {
    assert_eq!(bits(&actual.gpop), bits(&expected.gpop));
    assert_eq!(bits(&actual.grid.f), bits(&expected.grid.f));
    assert_eq!(bits(&actual.grid.fext), bits(&expected.grid.fext));
    assert_eq!(actual.last_wall_exchange(), expected.last_wall_exchange());
}

#[test]
fn uniform_spatial_fields_preserve_native_bits_through_source_steps() {
    for walls in [
        [ThermalWall::Temperature(1.0), ThermalWall::Temperature(0.0)],
        [ThermalWall::HeatFlux(0.0004), ThermalWall::Adiabatic],
        [ThermalWall::Convective { transfer: 0.03, ambient: 1.1 },
         ThermalWall::Convective { transfer: 0.02, ambient: -0.1 }],
    ] {
        let mut actual2 = ThermalLbm::slab(4, 6, 0.9, 0.8, 1e-4);
        let mut expected2 = ThermalLbm::slab(4, 6, 0.9, 0.8, 1e-4);
        let mut actual3 = ThermalLbm3::slab(3, 2, 6, 0.9, 0.8, 1e-4);
        let mut expected3 = ThermalLbm3::slab(3, 2, 6, 0.9, 0.8, 1e-4);
        actual2.set_spatial_thermal_walls(walls.map(|wall| vec![wall; 4])).unwrap();
        actual3.set_spatial_thermal_walls(walls.map(|wall| vec![wall; 6])).unwrap();
        expected2.set_thermal_walls(walls).unwrap();
        expected3.set_thermal_walls(walls).unwrap();
        for slab in [&mut actual2, &mut expected2] { slab.perturb(0.001); }
        for slab in [&mut actual3, &mut expected3] { slab.perturb(0.001); }
        let (n2, n3) = (actual2.gpop.len(), actual3.gpop.len());
        let mut heat2 = vec![0.0; n2];
        let mut heat3 = vec![0.0; n3];
        heat2[4..n2 - 4].fill(0.0001);
        heat3[6..n3 - 6].fill(0.0001);
        let force2 = vec![[1e-5, -2e-5]; n2];
        let force3 = vec![[1e-5, -2e-5, 1e-5]; n3];
        for _ in 0..8 {
            assert_eq!(
                actual2.step_with_sources(&heat2, &force2, limits()).unwrap(),
                expected2.step_with_sources(&heat2, &force2, limits()).unwrap(),
            );
            assert_eq!(
                actual3.step_with_sources(&heat3, &force3, limits()).unwrap(),
                expected3.step_with_sources(&heat3, &force3, limits()).unwrap(),
            );
            same2(&actual2, &expected2);
            same3(&actual3, &expected3);
            assert_eq!(actual2.rayleigh().to_bits(), expected2.rayleigh().to_bits());
            assert_eq!(actual3.nusselt().to_bits(), expected3.nusselt().to_bits());
        }
    }
}

fn patterned() -> (ThermalLbm, ThermalLbm3) {
    let mut two = ThermalLbm::slab(4, 6, 0.9, 0.8, 1e-4);
    let mut three = ThermalLbm3::slab(3, 2, 6, 0.9, 0.8, 1e-4);
    let bottom = [ThermalWall::HeatFlux(0.0001), ThermalWall::Adiabatic,
        ThermalWall::Temperature(0.9), ThermalWall::Convective { transfer: 0.03, ambient: 1.0 }];
    let top = [ThermalWall::Adiabatic, ThermalWall::Temperature(0.1),
        ThermalWall::HeatFlux(-0.0001), ThermalWall::Convective { transfer: 0.02, ambient: 0.1 }];
    two.set_spatial_thermal_walls([bottom.to_vec(), top.to_vec()]).unwrap();
    three.set_spatial_thermal_walls([
        (0..6).map(|i| bottom[i % 4]).collect(),
        (0..6).map(|i| top[i % 4]).collect(),
    ]).unwrap();
    (two, three)
}

#[test]
fn invalid_fields_and_late_source_rejection_preserve_patterns_reports_and_history() {
    let (mut actual2, mut actual3) = patterned();
    let (mut expected2, mut expected3) = patterned();
    for _ in 0..3 {
        actual2.step(); expected2.step();
        actual3.step(); expected3.step();
    }
    for lengths in [[3, 4], [4, 5], [0, 0]] {
        assert!(actual2.set_spatial_thermal_walls(lengths.map(|n| vec![ThermalWall::Adiabatic; n])).is_err());
    }
    for lengths in [[5, 6], [6, 7], [0, 0]] {
        assert!(actual3.set_spatial_thermal_walls(lengths.map(|n| vec![ThermalWall::Adiabatic; n])).is_err());
    }
    for invalid in [ThermalWall::Temperature(f64::NAN), ThermalWall::HeatFlux(f64::INFINITY),
        ThermalWall::Convective { transfer: -0.1, ambient: 0.5 }]
    {
        let mut walls2 = [vec![ThermalWall::Adiabatic; 4], vec![ThermalWall::Adiabatic; 4]];
        let mut walls3 = [vec![ThermalWall::Adiabatic; 6], vec![ThermalWall::Adiabatic; 6]];
        walls2[1][3] = invalid;
        walls3[1][5] = invalid;
        assert!(actual2.set_spatial_thermal_walls(walls2).is_err());
        assert!(actual3.set_spatial_thermal_walls(walls3).is_err());
    }
    same2(&actual2, &expected2);
    same3(&actual3, &expected3);
    assert_eq!(actual2.thermal_wall_fields(), expected2.thermal_wall_fields());
    assert_eq!(actual3.thermal_wall_fields(), expected3.thermal_wall_fields());

    let (n2, n3) = (actual2.gpop.len(), actual3.gpop.len());
    let mut heat2 = vec![0.0; n2];
    let mut heat3 = vec![0.0; n3];
    let force2 = vec![[0.0; 2]; n2];
    let force3 = vec![[0.0; 3]; n3];
    // The first half-kick is admitted. The final half-kick exceeds the
    // interval after the candidate has already streamed its patterned walls.
    heat2[actual2.grid.idx(0, 3)] = 2.0;
    heat3[actual3.grid.idx(0, 0, 3)] = 2.0;
    let tight = ThermalSourceLimits { temperature_range: [-0.5, 1.6], ..limits() };
    assert!(matches!(actual2.step_with_sources(&heat2, &force2, tight), Err(ThermalSourceError::Temperature { .. })));
    assert!(matches!(actual3.step_with_sources(&heat3, &force3, tight), Err(ThermalSourceError::Temperature { .. })));
    same2(&actual2, &expected2);
    same3(&actual3, &expected3);
    assert_eq!(actual2.thermal_wall_fields(), expected2.thermal_wall_fields());
    assert_eq!(actual3.thermal_wall_fields(), expected3.thermal_wall_fields());
    heat2.fill(0.0); heat3.fill(0.0);
    actual2.step_with_sources(&heat2, &force2, limits()).unwrap();
    expected2.step_with_sources(&heat2, &force2, limits()).unwrap();
    actual3.step_with_sources(&heat3, &force3, limits()).unwrap();
    expected3.step_with_sources(&heat3, &force3, limits()).unwrap();
    same2(&actual2, &expected2);
    same3(&actual3, &expected3);

    // Existing setters select a uniform wall again and invalidate old reports.
    actual2.set_thermal_walls([ThermalWall::Adiabatic; 2]).unwrap();
    actual3.reset_thermal_walls();
    assert_eq!(actual2.thermal_wall_fields(), [ThermalWallField::Uniform(ThermalWall::Adiabatic); 2]);
    assert_eq!(actual3.thermal_wall_fields(), actual3.thermal_walls().map(ThermalWallField::Uniform));
    assert!(actual2.last_wall_exchange().is_none());
    assert!(actual3.last_wall_exchange().is_none());
}

fn storage<const Q: usize>(populations: &[[f64; Q]], plane: usize) -> f64 {
    populations[plane..populations.len() - plane].iter().flatten().sum()
}

#[test]
fn local_heating_and_cooling_drive_opposite_buoyant_plumes_and_conserve_heat() {
    let mut two = ThermalLbm::slab(12, 8, 0.9, 1.0, 0.001);
    let mut three = ThermalLbm3::slab(12, 8, 8, 0.9, 1.0, 0.001);
    // Begin isothermal and quiescent. No seeded flow or mean temperature
    // difference can manufacture the response to these local surface loads.
    let equilibrium2 = two.gpop[0].map(|g| 0.5 * g);
    let equilibrium3 = three.gpop[0].map(|g| 0.5 * g);
    two.gpop.fill(equilibrium2);
    three.gpop.fill(equilibrium3);
    let mut walls2 = [vec![ThermalWall::Adiabatic; 12], vec![ThermalWall::Adiabatic; 12]];
    let mut walls3 = [vec![ThermalWall::Adiabatic; 96], vec![ThermalWall::Adiabatic; 96]];
    for x in [2, 3] {
        walls2[0][x] = ThermalWall::HeatFlux(0.002);
        walls2[0][x + 6] = ThermalWall::HeatFlux(-0.002);
        for y in [2, 3] {
            walls3[0][x + 12 * y] = ThermalWall::HeatFlux(0.002);
            walls3[0][x + 6 + 12 * (y + 4)] = ThermalWall::HeatFlux(-0.002);
        }
    }
    two.set_spatial_thermal_walls(walls2).unwrap();
    three.set_spatial_thermal_walls(walls3).unwrap();
    let initial2 = storage(&two.gpop, 12);
    let initial3 = storage(&three.gpop, 96);
    for _ in 0..120 { two.step(); three.step(); }
    assert!(two.temperature(2, 2) > 0.5 && two.temperature(8, 2) < 0.5);
    assert!(three.temperature(2, 2, 2) > 0.5 && three.temperature(8, 6, 2) < 0.5);
    let hot2 = two.grid.moments(two.grid.idx(2, 2)).u[1];
    let cold2 = two.grid.moments(two.grid.idx(8, 2)).u[1];
    let hot3 = three.grid.moments3(three.grid.idx(2, 2, 2)).1[2];
    let cold3 = three.grid.moments3(three.grid.idx(8, 6, 2)).1[2];
    assert!(hot2 > 1e-5 && cold2 < -1e-5, "2D plume velocities: {hot2:e}, {cold2:e}");
    assert!(hot3 > 1e-5 && cold3 < -1e-5, "3D plume velocities: {hot3:e}, {cold3:e}");
    assert!((storage(&two.gpop, 12) - initial2).abs() < 1e-10);
    assert!((storage(&three.gpop, 96) - initial3).abs() < 1e-10);
    for report in [two.last_wall_exchange().unwrap(), three.last_wall_exchange().unwrap()] {
        assert!(report[0].heat_in.abs() < 1e-15);
        assert_eq!(report[1].heat_in, 0.0);
    }
}
