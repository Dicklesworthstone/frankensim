use fs_lbm::material::{
    MaterialSlab, MaterialSlab3, MaterialThermalSlab, MaterialThermalSlab3,
    MaterialThermalWall, plan_scaling_from_material,
};
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::fluid_transport::FluidTransportState;

const DX: f64 = 0.1;
const DT: f64 = 0.05;
const DEPTH: f64 = 0.25;
const RHO: f64 = 2.0;
const CP: f64 = 5.0;
const K: f64 = 0.2;
const NX: usize = 3;
const NZ: usize = 6;

fn slabs(reverse: bool, gravity: f64) -> (MaterialSlab, MaterialSlab3) {
    let material = FluidTransportState::try_new(RHO, CP, K, 0.2).unwrap();
    let scale = plan_scaling_from_material(&material, DX, DT, 0.1).unwrap();
    let (bottom, top) = if reverse { (300.0, 320.0) } else { (320.0, 300.0) };
    let two = scale.thermal_slab(MaterialThermalSlab {
        nx: NX, ny: NZ, bottom_temperature_k: bottom, top_temperature_k: top,
        gravity_m_s2: gravity, expansion_per_k: 1e-3, max_cells: 128,
    }).unwrap();
    let three = scale.thermal_slab3(MaterialThermalSlab3 {
        nx: NX, ny: NX, nz: NZ, bottom_temperature_k: bottom, top_temperature_k: top,
        gravity_m_s2: gravity, expansion_per_k: 1e-3, max_cells: 128,
    }).unwrap();
    (two, three)
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:.3e}");
}

fn limits() -> HeatSourceLimits {
    HeatSourceLimits {
        max_cells: 128, temperature_range_k: [275.0, 370.0], balance_tolerance_j: 1e-11,
    }
}

fn fluid_sum<const Q: usize>(populations: &[[f64; Q]], plane: usize) -> f64 {
    populations[plane..populations.len() - plane].iter().flatten().sum()
}

#[test]
fn g1_two_robin_walls_match_physical_series_resistances_in_both_dimensions() {
    let (mut two, mut three) = slabs(false, 0.0);
    let walls = [
        MaterialThermalWall::Convective { h_w_m2_k: 0.7, ambient_temperature_k: 320.0 },
        MaterialThermalWall::Convective { h_w_m2_k: 0.4, ambient_temperature_k: 300.0 },
    ];
    two.set_thermal_walls(walls).unwrap();
    three.set_thermal_walls(walls).unwrap();
    assert!(two.last_wall_exchange(DEPTH).unwrap().is_none());
    assert!(three.last_wall_exchange().unwrap().is_none());
    for _ in 0..6000 {
        two.lattice.step();
        three.lattice.step();
    }
    let length = NZ as f64 * DX;
    let heat_flux = 20.0 / (1.0 / 0.7 + length / K + 1.0 / 0.4);
    let surface = [320.0 - heat_flux / 0.7, 300.0 + heat_flux / 0.4];
    let two_report = two.last_wall_exchange(DEPTH).unwrap().unwrap();
    let three_report = three.last_wall_exchange().unwrap().unwrap();
    for side in 0..2 {
        let signed_flux = if side == 0 { heat_flux } else { -heat_flux };
        close(two_report[side].heat_in_j, signed_flux * NX as f64 * DX * DEPTH * DT, 1e-11);
        close(three_report[side].heat_in_j, signed_flux * (NX * NX) as f64 * DX * DX * DT, 1e-11);
        close(two_report[side].mean_temperature_k, surface[side], 2e-8);
        close(three_report[side].mean_temperature_k, surface[side], 2e-8);
        for wall in [two_report[side], three_report[side]] {
            close(wall.min_temperature_k, surface[side], 2e-8);
            close(wall.max_temperature_k, surface[side], 2e-8);
        }
    }
    for layer in 1..=NZ {
        let expected = surface[0] - heat_flux * (layer as f64 - 0.5) * DX / K;
        close(two.temperature_si(two.lattice.temperature(0, layer)).unwrap(), expected, 2e-8);
        close(three.temperature_si(three.lattice.temperature(0, 0, layer)).unwrap(), expected, 2e-8);
    }
    // The surface trace differs from its adjacent cell by the half-cell drop.
    close(two_report[0].mean_temperature_k
        - two.temperature_si(two.lattice.temperature(0, 1)).unwrap(), heat_flux * 0.5 * DX / K, 2e-8);
    let deeper = two.last_wall_exchange(2.0 * DEPTH).unwrap().unwrap();
    assert_eq!(deeper[0].heat_in_j, 2.0 * two_report[0].heat_in_j);
    assert_eq!(deeper[0].mean_temperature_k, two_report[0].mean_temperature_k);
}

#[test]
fn g1_g3_imposed_wall_heat_and_volume_sources_close_with_signed_temperature_coordinates() {
    for reverse in [false, true] {
        let (mut two, mut three) = slabs(reverse, 9.81);
        let walls = [MaterialThermalWall::HeatFlux(25.0), MaterialThermalWall::Adiabatic];
        two.set_thermal_walls(walls).unwrap();
        three.set_thermal_walls(walls).unwrap();
        two.lattice.perturb(0.001);
        three.lattice.perturb(0.001);
        let count2 = two.lattice.gpop.len();
        let count3 = three.lattice.gpop.len();
        let mut power2 = vec![0.0; count2];
        let mut power3 = vec![0.0; count3];
        power2[NX..count2 - NX].fill(4.0);
        power3[NX * NX..count3 - NX * NX].fill(4.0);
        let force2 = vec![[0.01, 0.0]; count2];
        let force3 = vec![[0.01, 0.0, 0.0]; count3];
        let energy2 = RHO * CP * DX * DX * DEPTH * two.temperature_difference_k();
        let energy3 = RHO * CP * DX * DX * DX * three.temperature_difference_k();
        for _ in 0..8 {
            let before2 = fluid_sum(&two.lattice.gpop, NX);
            let before3 = fluid_sum(&three.lattice.gpop, NX * NX);
            let report2 = two.step_with_power_density(&power2, &force2, DEPTH, limits()).unwrap();
            let report3 = three.step_with_power_density(&power3, &force3, limits()).unwrap();
            let wall2 = two.last_wall_exchange(DEPTH).unwrap().unwrap();
            let wall3 = three.last_wall_exchange().unwrap().unwrap();
            let raw2 = two.lattice.last_wall_exchange().unwrap();
            let raw3 = three.lattice.last_wall_exchange().unwrap();
            let (offset, scale) = if reverse { (320.0, -20.0) } else { (300.0, 20.0) };
            for (native, physical) in raw2.iter().zip(&wall2).chain(raw3.iter().zip(&wall3)) {
                let a = offset + scale * native.min_temperature;
                let b = offset + scale * native.max_temperature;
                close(physical.min_temperature_k, a.min(b), 1e-12);
                close(physical.max_temperature_k, a.max(b), 1e-12);
                assert!(physical.min_temperature_k <= physical.mean_temperature_k + 1e-12);
                assert!(physical.mean_temperature_k <= physical.max_temperature_k + 1e-12);
            }
            // The imposed roll creates distinct local traces, so the reversed
            // coordinate case cannot pass with accidentally swapped bounds.
            assert!(wall2[0].min_temperature_k < wall2[0].max_temperature_k);
            assert!(wall3[0].min_temperature_k < wall3[0].max_temperature_k);
            let incoming2 = 25.0 * NX as f64 * DX * DEPTH * DT;
            let incoming3 = 25.0 * (NX * NX) as f64 * DX * DX * DT;
            close(wall2[0].heat_in_j, incoming2, 1e-13);
            close(wall3[0].heat_in_j, incoming3, 1e-13);
            assert_eq!(wall2[1].heat_in_j, 0.0);
            assert_eq!(wall3[1].heat_in_j, 0.0);
            close(report2.bottom_wall_j, wall2[0].heat_in_j, 1e-13);
            close(report3.bottom_wall_j, wall3[0].heat_in_j, 1e-13);
            close((fluid_sum(&two.lattice.gpop, NX) - before2) * energy2,
                incoming2 + 4.0 * (NX * NZ) as f64 * DX * DX * DEPTH * DT, 2e-12);
            close((fluid_sum(&three.lattice.gpop, NX * NX) - before3) * energy3,
                incoming3 + 4.0 * (NX * NX * NZ) as f64 * DX * DX * DX * DT, 2e-12);
            assert!(report2.stored_change_j > 0.0 && report3.stored_change_j > 0.0);
        }
        let cell2 = two.lattice.grid.idx(0, 3);
        let cell3 = three.lattice.grid.idx(0, 0, 3);
        assert!(two.lattice.grid.moments(cell2).u[0] > 0.0);
        assert!(three.lattice.grid.moments3(cell3).1[0] > 0.0);
        assert_ne!(two.lattice.grid.fext[cell2][1], 0.0);
        assert_ne!(three.lattice.grid.fext[cell3][2], 0.0);

        // A refused next candidate retains the previously accepted report.
        let saved2 = (two.lattice.gpop.clone(), two.lattice.grid.f.clone(), two.lattice.last_wall_exchange());
        let saved3 = (three.lattice.gpop.clone(), three.lattice.grid.f.clone(), three.lattice.last_wall_exchange());
        let laws2 = two.lattice.thermal_walls();
        let laws3 = three.lattice.thermal_walls();
        power2[NX..count2 - NX].fill(1e9);
        power3[NX * NX..count3 - NX * NX].fill(1e9);
        assert!(two.step_with_power_density(&power2, &force2, DEPTH, limits()).is_err());
        assert!(three.step_with_power_density(&power3, &force3, limits()).is_err());
        assert_eq!((&two.lattice.gpop, &two.lattice.grid.f, two.lattice.last_wall_exchange()),
            (&saved2.0, &saved2.1, saved2.2));
        assert_eq!((&three.lattice.gpop, &three.lattice.grid.f, three.lattice.last_wall_exchange()),
            (&saved3.0, &saved3.1, saved3.2));
        assert_eq!(two.lattice.thermal_walls(), laws2);
        assert_eq!(three.lattice.thermal_walls(), laws3);
    }
}

#[test]
fn g0_insulated_material_slabs_preserve_nonuniform_thermal_inventory() {
    let (mut two, mut three) = slabs(false, 0.0);
    two.set_thermal_walls([MaterialThermalWall::Adiabatic; 2]).unwrap();
    three.set_thermal_walls([MaterialThermalWall::Adiabatic; 2]).unwrap();
    let before2 = fluid_sum(&two.lattice.gpop, NX);
    let before3 = fluid_sum(&three.lattice.gpop, NX * NX);
    let first = two.lattice.temperature(0, 1);
    for _ in 0..100 {
        two.lattice.step();
        three.lattice.step();
    }
    close(fluid_sum(&two.lattice.gpop, NX), before2, 2e-12);
    close(fluid_sum(&three.lattice.gpop, NX * NX), before3, 2e-12);
    assert!((two.lattice.temperature(0, 1) - first).abs() > 0.1);
    for wall in two.last_wall_exchange(DEPTH).unwrap().unwrap() {
        assert_eq!(wall.heat_in_j, 0.0);
    }
    for wall in three.last_wall_exchange().unwrap().unwrap() {
        assert_eq!(wall.heat_in_j, 0.0);
    }
}

#[test]
fn g0_invalid_physical_wall_data_leaves_both_original_laws_and_last_report_intact() {
    let (mut two, mut three) = slabs(false, 0.0);
    let original = [MaterialThermalWall::HeatFlux(1.0), MaterialThermalWall::Temperature(302.0)];
    two.set_thermal_walls(original).unwrap();
    three.set_thermal_walls(original).unwrap();
    two.lattice.step();
    three.lattice.step();
    let laws2 = two.lattice.thermal_walls();
    let laws3 = three.lattice.thermal_walls();
    let report2 = two.lattice.last_wall_exchange();
    let report3 = three.lattice.last_wall_exchange();
    for bad in [
        MaterialThermalWall::Temperature(0.0),
        MaterialThermalWall::Temperature(f64::NAN),
        MaterialThermalWall::HeatFlux(f64::NAN),
        MaterialThermalWall::HeatFlux(f64::from_bits(1)),
        MaterialThermalWall::Convective { h_w_m2_k: -1.0, ambient_temperature_k: 300.0 },
        MaterialThermalWall::Convective { h_w_m2_k: f64::INFINITY, ambient_temperature_k: 300.0 },
        MaterialThermalWall::Convective { h_w_m2_k: 0.0, ambient_temperature_k: 0.0 },
    ] {
        assert!(two.set_thermal_walls([MaterialThermalWall::Adiabatic, bad]).is_err());
        assert!(three.set_thermal_walls([MaterialThermalWall::Adiabatic, bad]).is_err());
        assert_eq!(two.lattice.thermal_walls(), laws2);
        assert_eq!(three.lattice.thermal_walls(), laws3);
        assert_eq!(two.lattice.last_wall_exchange(), report2);
        assert_eq!(three.lattice.last_wall_exchange(), report3);
    }
    assert!(two.last_wall_exchange(0.0).is_err());
    assert!(two.last_wall_exchange(f64::INFINITY).is_err());
    two.set_thermal_walls([MaterialThermalWall::Adiabatic; 2]).unwrap();
    three.set_thermal_walls([MaterialThermalWall::Adiabatic; 2]).unwrap();
    assert!(two.last_wall_exchange(DEPTH).unwrap().is_none());
    assert!(three.last_wall_exchange().unwrap().is_none());
}

#[test]
fn g0_material_queries_admit_active_temperatures_and_declared_walls_without_ghost_storage() {
    let gate = fs_exec::CancelGate::new();
    let tight = HeatSourceLimits { temperature_range_k: [309.0, 311.0], ..limits() };
    for reverse in [false, true] {
        let (mut two, mut three) = slabs(reverse, 0.0);
        let count2 = two.lattice.gpop.len();
        let count3 = three.lattice.gpop.len();
        let plane3 = NX * NX;
        // The constructor's bottom ghost is unit-temperature equilibrium.
        // Initialize only fluid cells at theta=0.5 (310 K in either mapping).
        // Its 300/320 K ghost planes must remain untouched and irrelevant.
        let equilibrium2 = two.lattice.gpop[0].map(|value| 0.5 * value);
        let equilibrium3 = three.lattice.gpop[0].map(|value| 0.5 * value);
        two.lattice.gpop[NX..count2 - NX].fill(equilibrium2);
        three.lattice.gpop[plane3..count3 - plane3].fill(equilibrium3);
        let ghosts2 = (two.lattice.gpop[..NX].to_vec(), two.lattice.gpop[count2 - NX..].to_vec());
        let ghosts3 = (three.lattice.gpop[..plane3].to_vec(), three.lattice.gpop[count3 - plane3..].to_vec());
        two.set_thermal_walls([MaterialThermalWall::Adiabatic; 2]).unwrap();
        three.set_thermal_walls([MaterialThermalWall::Adiabatic; 2]).unwrap();
        let power2 = vec![0.0; count2];
        let power3 = vec![0.0; count3];
        let force2 = vec![[0.0; 2]; count2];
        let force3 = vec![[0.0; 3]; count3];
        let mut queries2 = 0;
        let mut queries3 = 0;
        let report2 = two.step_with_material_transport(
            &power2, &force2, DEPTH, tight, &gate,
            |cell, kelvin, retained| {
                assert_eq!(cell, NX + queries2 % (count2 - 2 * NX));
                close(kelvin, 310.0, 1e-12);
                queries2 += 1;
                Ok(*retained)
            },
        ).unwrap();
        let report3 = three.step_with_material_transport(
            &power3, &force3, tight, &gate,
            |cell, kelvin, retained| {
                assert_eq!(cell, plane3 + queries3 % (count3 - 2 * plane3));
                close(kelvin, 310.0, 1e-12);
                queries3 += 1;
                Ok(*retained)
            },
        ).unwrap();
        assert_eq!(queries2, 2 * (count2 - 2 * NX));
        assert_eq!(queries3, 2 * (count3 - 2 * plane3));
        assert_eq!(report2.applied.len(), count2 - 2 * NX);
        assert_eq!(report2.endpoint.len(), count2 - 2 * NX);
        assert_eq!(report3.applied.len(), count3 - 2 * plane3);
        assert_eq!(report3.endpoint.len(), count3 - 2 * plane3);
        for heat in [report2.heat, report3.heat] {
            assert_eq!(heat.bottom_wall_j, 0.0);
            assert_eq!(heat.top_wall_j, 0.0);
            close(heat.stored_change_j, 0.0, 1e-12);
        }
        assert_eq!((&two.lattice.gpop[..NX], &two.lattice.gpop[count2 - NX..]),
            (ghosts2.0.as_slice(), ghosts2.1.as_slice()));
        assert_eq!((&three.lattice.gpop[..plane3], &three.lattice.gpop[count3 - plane3..]),
            (ghosts3.0.as_slice(), ghosts3.1.as_slice()));

        // A real fixed-temperature reservoir outside the material interval
        // must still refuse before querying the material or publishing state.
        let walls = [MaterialThermalWall::Temperature(320.0), MaterialThermalWall::Adiabatic];
        two.set_thermal_walls(walls).unwrap();
        three.set_thermal_walls(walls).unwrap();
        let saved2 = (two.lattice.gpop.clone(), two.lattice.grid.f.clone(),
            two.lattice.grid.fext.clone(), two.lattice.last_wall_exchange());
        let saved3 = (three.lattice.gpop.clone(), three.lattice.grid.f.clone(),
            three.lattice.grid.fext.clone(), three.lattice.last_wall_exchange());
        let laws2 = two.lattice.thermal_walls();
        let laws3 = three.lattice.thermal_walls();
        queries2 = 0;
        queries3 = 0;
        assert!(two.step_with_material_transport(
            &power2, &force2, DEPTH, tight, &gate,
            |_, _, retained| { queries2 += 1; Ok(*retained) },
        ).is_err());
        assert!(three.step_with_material_transport(
            &power3, &force3, tight, &gate,
            |_, _, retained| { queries3 += 1; Ok(*retained) },
        ).is_err());
        assert_eq!((queries2, queries3), (0, 0));
        assert_eq!((&two.lattice.gpop, &two.lattice.grid.f, &two.lattice.grid.fext,
            two.lattice.last_wall_exchange()), (&saved2.0, &saved2.1, &saved2.2, saved2.3));
        assert_eq!((&three.lattice.gpop, &three.lattice.grid.f, &three.lattice.grid.fext,
            three.lattice.last_wall_exchange()), (&saved3.0, &saved3.1, &saved3.2, saved3.3));
        assert_eq!(two.lattice.thermal_walls(), laws2);
        assert_eq!(three.lattice.thermal_walls(), laws3);
    }
}
