use super::*;

fn walls(links: usize) -> [Vec<MaterialThermalWall>; 2] {
    [
        (0..links).map(|i| MaterialThermalWall::HeatFlux(1.0 + 0.25 * i as f64)).collect(),
        (0..links).map(|i| match i % 3 {
            0 => MaterialThermalWall::Temperature(310.5),
            1 => MaterialThermalWall::Adiabatic,
            _ => MaterialThermalWall::Convective {
                h_w_m2_k: 0.7 + 0.1 * i as f64, ambient_temperature_k: 312.0,
            },
        }).collect(),
    ]
}

// Independent physical pair closure at an initially isothermal 310 K fluid.
// A lattice normal link has surface resistance dt/(2*w*rho*Cp*dx).
fn physical_pair(wall: MaterialThermalWall, weight: f64) -> (f64, f64) {
    let resistance = DT / (2.0 * weight * RHO * CP * DX);
    let flux = match wall {
        MaterialThermalWall::Temperature(kelvin) => (kelvin - 310.0) / resistance,
        MaterialThermalWall::Adiabatic => 0.0,
        MaterialThermalWall::HeatFlux(flux) => flux,
        MaterialThermalWall::Convective { h_w_m2_k, ambient_temperature_k } => {
            h_w_m2_k * (ambient_temperature_k - 310.0) / (1.0 + h_w_m2_k * resistance)
        }
    };
    (flux, 310.0 + resistance * flux)
}

#[test]
fn mixed_spatial_si_laws_preserve_local_heat_traces_and_signed_coordinates() {
    for reverse in [false, true] {
        let (mut two, mut three) = slabs(reverse, 0.0);
        let equilibrium2 = two.lattice.gpop[0].map(|g| 0.5 * g);
        let equilibrium3 = three.lattice.gpop[0].map(|g| 0.5 * g);
        two.lattice.gpop.fill(equilibrium2);
        three.lattice.gpop.fill(equilibrium3);
        let walls2 = walls(NX);
        let walls3 = walls(NX * NX);
        two.set_spatial_thermal_walls(walls2.clone()).unwrap();
        three.set_spatial_thermal_walls(walls3.clone()).unwrap();
        let n2 = two.lattice.gpop.len();
        let n3 = three.lattice.gpop.len();
        let heat2 = two.step_with_power_density(
            &vec![0.0; n2], &vec![[0.0; 2]; n2], DEPTH, limits(),
        ).unwrap();
        let heat3 = three.step_with_power_density(
            &vec![0.0; n3], &vec![[0.0; 3]; n3], limits(),
        ).unwrap();
        let report2 = two.last_wall_exchange(DEPTH).unwrap().unwrap();
        let report3 = three.last_wall_exchange().unwrap().unwrap();
        for side in 0..2 {
            let layer = if side == 0 { 1 } else { NZ };
            for (field, weight, area, report) in [
                (&walls2[side], 1.0 / 6.0, DX * DEPTH, report2[side]),
                (&walls3[side], 1.0 / 8.0, DX * DX, report3[side]),
            ] {
                let expected: Vec<_> = field.iter().map(|&law| physical_pair(law, weight)).collect();
                let incoming = expected.iter().map(|pair| pair.0).sum::<f64>() * area * DT;
                let mean = expected.iter().map(|pair| pair.1).sum::<f64>() / field.len() as f64;
                let minimum = expected.iter().map(|pair| pair.1).fold(f64::INFINITY, f64::min);
                let maximum = expected.iter().map(|pair| pair.1).fold(f64::NEG_INFINITY, f64::max);
                close(report.heat_in_j, incoming, 2e-13);
                close(report.mean_temperature_k, mean, 2e-12);
                close(report.min_temperature_k, minimum, 2e-12);
                close(report.max_temperature_k, maximum, 2e-12);
                for (link, (flux, _)) in expected.into_iter().enumerate() {
                    let actual = if field.len() == NX {
                        two.temperature_si(two.lattice.temperature(link, layer)).unwrap()
                    } else {
                        three.temperature_si(three.lattice.temperature(link % NX, link / NX, layer)).unwrap()
                    };
                    close(actual, 310.0 + flux * DT / (RHO * CP * DX), 2e-12);
                }
            }
        }
        for (heat, report) in [(heat2, report2), (heat3, report3)] {
            close(heat.bottom_wall_j, report[0].heat_in_j, 1e-13);
            close(heat.top_wall_j, report[1].heat_in_j, 1e-13);
            close(heat.stored_change_j, report[0].heat_in_j + report[1].heat_in_j, 1e-12);
            assert_eq!(heat.source_j, 0.0);
        }
    }
}

#[test]
fn invalid_si_patches_are_atomic_and_material_callbacks_admit_every_patch() {
    let (mut two, mut three) = slabs(true, 0.0);
    two.set_spatial_thermal_walls(walls(NX)).unwrap();
    three.set_spatial_thermal_walls(walls(NX * NX)).unwrap();
    two.lattice.step(); three.lattice.step();
    let before2 = (two.lattice.gpop.clone(), two.lattice.grid.f.clone(), two.lattice.last_wall_exchange());
    let before3 = (three.lattice.gpop.clone(), three.lattice.grid.f.clone(), three.lattice.last_wall_exchange());
    let fields2 = two.lattice.thermal_wall_fields().map(|field| field.laws().collect::<Vec<_>>());
    let fields3 = three.lattice.thermal_wall_fields().map(|field| field.laws().collect::<Vec<_>>());
    assert!(two.set_spatial_thermal_walls(walls(NX - 1)).is_err());
    assert!(three.set_spatial_thermal_walls(walls(NX * NX - 1)).is_err());
    for invalid in [MaterialThermalWall::Temperature(0.0), MaterialThermalWall::HeatFlux(f64::INFINITY),
        MaterialThermalWall::Convective { h_w_m2_k: -0.1, ambient_temperature_k: 310.0 }]
    {
        let mut bad2 = walls(NX);
        let mut bad3 = walls(NX * NX);
        bad2[1][NX - 1] = invalid;
        bad3[1][NX * NX - 1] = invalid;
        assert!(two.set_spatial_thermal_walls(bad2).is_err());
        assert!(three.set_spatial_thermal_walls(bad3).is_err());
    }
    assert_eq!((&two.lattice.gpop, &two.lattice.grid.f, two.lattice.last_wall_exchange()), (&before2.0, &before2.1, before2.2));
    assert_eq!((&three.lattice.gpop, &three.lattice.grid.f, three.lattice.last_wall_exchange()), (&before3.0, &before3.1, before3.2));
    assert_eq!(two.lattice.thermal_wall_fields().map(|field| field.laws().collect::<Vec<_>>()), fields2);
    assert_eq!(three.lattice.thermal_wall_fields().map(|field| field.laws().collect::<Vec<_>>()), fields3);

    // A finite declaration can still lie outside this step's material range.
    // Place it last on the upper plate so checking a default or only the first
    // patch cannot accidentally pass. No material query may run on refusal.
    let mut outside2 = walls(NX);
    let mut outside3 = walls(NX * NX);
    outside2[1][NX - 1] = MaterialThermalWall::Temperature(390.0);
    outside3[1][NX * NX - 1] = MaterialThermalWall::Temperature(390.0);
    two.set_spatial_thermal_walls(outside2).unwrap();
    three.set_spatial_thermal_walls(outside3).unwrap();
    let mut queries = 0;
    let gate = CancelGate::new_clock_free();
    let (n2, n3) = (two.lattice.gpop.len(), three.lattice.gpop.len());
    assert!(two.step_with_material_transport(
        &vec![0.0; n2], &vec![[0.0; 2]; n2], DEPTH, limits(), &gate,
        |_, _, retained| { queries += 1; Ok(*retained) },
    ).is_err());
    assert!(three.step_with_material_transport(
        &vec![0.0; n3], &vec![[0.0; 3]; n3], limits(), &gate,
        |_, _, retained| { queries += 1; Ok(*retained) },
    ).is_err());
    assert_eq!(queries, 0);
    assert_eq!((&two.lattice.gpop, &two.lattice.grid.f), (&before2.0, &before2.1));
    assert_eq!((&three.lattice.gpop, &three.lattice.grid.f), (&before3.0, &before3.1));
}

#[test]
fn patterned_wall_sgs_reconstructs_current_fluid_gradients_without_uniform_defaults() {
    let mut patterned = seeded_channel_slab(false);
    let mut reconstructed = seeded_channel_slab(false);
    let plane = patterned.lattice.grid.nx * patterned.lattice.grid.ny;
    patterned.set_spatial_thermal_walls([
        (0..plane).map(|i| MaterialThermalWall::Temperature(315.0 + (i % 3) as f64)).collect(),
        (0..plane).map(|i| MaterialThermalWall::Temperature(303.0 + (i % 5) as f64)).collect(),
    ]).unwrap();
    // Nonzero flux walls select the independently existing three-fluid-layer
    // reconstruction. At the same accepted state a patterned wall uses that
    // same gradient, while its actual heat remains with the local wall law.
    reconstructed.set_thermal_walls([
        MaterialThermalWall::HeatFlux(1.0), MaterialThermalWall::HeatFlux(-1.0),
    ]).unwrap();
    let bytes = patterned.lattice.required_thermal_channel_sgs_workspace_bytes().unwrap();
    let mut work = patterned.lattice.thermal_channel_sgs_workspace(bytes).unwrap();
    let mut reference = reconstructed.lattice.thermal_channel_sgs_workspace(bytes).unwrap();
    let gate = CancelGate::new_clock_free();
    let actual = patterned.lattice.prepare_thermal_channel_sgs(channel_config(), None, &mut work, &gate).unwrap();
    let expected = reconstructed.lattice.prepare_thermal_channel_sgs(channel_config(), None, &mut reference, &gate).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(work.relaxation().momentum, reference.relaxation().momentum);
    assert_eq!(work.relaxation().thermal, reference.relaxation().thermal);
}
