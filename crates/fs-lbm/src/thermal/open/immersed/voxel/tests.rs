use super::*;
use crate::thermal::open::{ThermalOpenBalance3, ThermalOpenBoundaries3, ThermalOpenLimits3};

mod orthotropic;

fn voxel(cell: usize, temperature: f64, capacity: f64, conductivity: f64) -> ImmersedThermalVoxel3 {
    ImmersedThermalVoxel3 {
        cell,
        temperature,
        heat_capacity: capacity,
        conductivity,
        contact_resistance: 0.0,
        heat_per_step: 0.0,
    }
}

fn field(dims: [usize; 3], voxels: &[ImmersedThermalVoxel3]) -> ThermalOpenField3 {
    let count = dims.into_iter().product();
    let mut solid = vec![false; count];
    for voxel in voxels {
        solid[voxel.cell] = true;
    }
    ThermalOpenField3::try_new::<ThermalOpenError3>(
        dims,
        &solid,
        &vec![0.2; count],
        0.85,
        ThermalOpenBoundaries3 {
            faces: [ThermalOpenFace3::Periodic; 6],
            internal_wall: ThermalWall::Adiabatic,
        },
        ThermalOpenLimits3 {
            maximum_cells: count,
            maximum_storage_bytes: 1024 * 1024,
            temperature_range: [-2.0, 3.0],
            balance_tolerance: 1e-10,
        },
        0,
        |_| Ok([0.0; 3]),
        || Ok(()),
    )
    .unwrap()
}

fn advance(field: &mut ThermalOpenField3) -> ThermalOpenBalance3 {
    let heat = vec![0.0; field.populations().len()];
    let trial = field
        .begin_step::<ThermalOpenError3>(&heat, None, || Ok(()))
        .unwrap();
    let balance = field
        .finish_step::<ThermalOpenError3>(&heat, trial, || Ok(()))
        .unwrap();
    field.publish(field.completed_steps() + 1, balance);
    balance
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

fn expected_pair(
    a: ImmersedThermalVoxel3,
    b: ImmersedThermalVoxel3,
    faces: usize,
    g: f64,
) -> [f64; 2] {
    let ka = faces as f64 / (4.0 + 0.5 / a.conductivity + a.contact_resistance);
    let kb = faces as f64 / (4.0 + 0.5 / b.conductivity + b.contact_resistance);
    let aa = a.heat_capacity + ka + g;
    let bb = b.heat_capacity + kb + g;
    let ra = a.heat_capacity * a.temperature + 0.2 * ka + a.heat_per_step;
    let rb = b.heat_capacity * b.temperature + 0.2 * kb + b.heat_per_step;
    let determinant = aa * bb - g * g;
    [
        (bb * ra + g * rb) / determinant,
        (aa * rb + g * ra) / determinant,
    ]
}

#[test]
fn g1_heterogeneous_half_cells_match_the_analytic_conjugate_pair() {
    let mut a = voxel(27, 1.4, 2.0, 2.0);
    let mut b = voxel(28, 0.4, 5.0, 0.5);
    a.heat_per_step = 0.03;
    b.heat_per_step = -0.01;
    b.contact_resistance = 0.25;
    let expected = expected_pair(a, b, 5, 0.8);
    let mut owner = field([5, 4, 3], &[a, b]);
    let mapping = owner
        .set_immersed_thermal_voxels(&[b, a], Default::default())
        .unwrap();
    assert_eq!(mapping.cells(), &[27, 28]);
    assert_eq!(mapping.conduction_links(), 1);
    assert_eq!(mapping.node_for_cell(28), Some(1));
    assert_eq!(mapping.node_for_cell(29), None);
    assert_eq!(owner.immersed_body(0).unwrap().solid.resistance, 0.25);
    assert_eq!(owner.immersed_body(1).unwrap().solid.resistance, 1.25);
    let balance = advance(&mut owner);
    for (node, temperature) in expected.into_iter().enumerate() {
        near(
            owner.immersed_body(node).unwrap().solid.temperature,
            temperature,
            2e-14,
        );
    }
    let left = owner.immersed_body(0).unwrap().last_exchange.unwrap();
    let right = owner.immersed_body(1).unwrap().last_exchange.unwrap();
    near(left.conducted_in, 0.8 * (expected[1] - expected[0]), 2e-14);
    assert_eq!(left.conducted_in, -right.conducted_in);
    near(balance.combined_defect, 0.0, 1e-13);
}

#[test]
fn g3_periodic_parallel_faces_and_input_order_keep_actual_solid_geometry() {
    // With two x cells both x faces join the same two solid nodes. Four y/z
    // faces per node remain exposed; both x conductances must be retained.
    let a = voxel(8, 1.1, 2.0, 0.8);
    let b = voxel(9, 0.3, 3.0, 0.8);
    let expected = expected_pair(a, b, 4, 1.6);
    let mut first = field([2, 3, 3], &[a, b]);
    let mut reversed = field([2, 3, 3], &[a, b]);
    let map = first
        .set_immersed_thermal_voxels(&[a, b], Default::default())
        .unwrap();
    let other = reversed
        .set_immersed_thermal_voxels(&[b, a], Default::default())
        .unwrap();
    assert_eq!(map, other);
    assert_eq!(map.conduction_links(), 1);
    advance(&mut first);
    advance(&mut reversed);
    for (node, temperature) in expected.into_iter().enumerate() {
        assert_eq!(first.immersed_body(node), reversed.immersed_body(node));
        assert_eq!(first.immersed_body(node).unwrap().interface_links, 4);
        near(
            first.immersed_body(node).unwrap().solid.temperature,
            temperature,
            2e-14,
        );
    }
    assert_eq!(first.populations(), reversed.populations());
}

#[test]
fn g1_a_heated_interior_voxel_conducts_to_the_fluid_with_one_energy_inventory() {
    let mut voxels = Vec::new();
    for z in 1..4 {
        for y in 1..4 {
            for x in 1..4 {
                let cell = x + 5 * y + 25 * z;
                let mut node = voxel(cell, if cell == 62 { 1.4 } else { 0.2 }, 2.0, 0.7);
                if cell == 62 {
                    node.heat_per_step = 0.05;
                }
                voxels.push(node);
            }
        }
    }
    let mut owner = field([5, 5, 5], &voxels);
    let initial = owner.populations().iter().flatten().sum::<f64>()
        + voxels
            .iter()
            .map(|v| v.heat_capacity * v.temperature)
            .sum::<f64>();
    let map = owner
        .set_immersed_thermal_voxels(&voxels, Default::default())
        .unwrap();
    assert_eq!(map.conduction_links(), 54);
    let core = map.node_for_cell(62).unwrap();
    assert_eq!(owner.immersed_body(core).unwrap().interface_links, 0);
    let balance = advance(&mut owner);
    assert!(owner.immersed_body(core).unwrap().solid.temperature < 1.4);
    assert!(
        owner
            .immersed_body(map.node_for_cell(61).unwrap())
            .unwrap()
            .solid
            .temperature
            > 0.2
    );
    let final_energy = owner.populations().iter().flatten().sum::<f64>()
        + (0..map.cells().len())
            .map(|i| {
                let solid = owner.immersed_body(i).unwrap().solid;
                solid.heat_capacity * solid.temperature
            })
            .sum::<f64>();
    near(final_energy - initial, 0.05, 2e-12);
    near(balance.combined_defect, 0.0, 2e-12);
    near(
        balance.internal_wall.unwrap().solid.unwrap().source,
        0.05,
        1e-15,
    );
}

#[test]
fn g0_invalid_voxel_topology_or_budget_preserves_the_installed_network() {
    let a = voxel(27, 1.1, 2.0, 0.8);
    let b = voxel(28, 0.3, 3.0, 0.8);
    let mut owner = field([5, 4, 3], &[a, b]);
    owner
        .set_immersed_thermal_voxels(&[a, b], Default::default())
        .unwrap();
    advance(&mut owner);
    let accepted = [owner.immersed_body(0), owner.immersed_body(1)];
    let populations = owner.populations().to_vec();
    for bad in [
        vec![a, a],
        vec![ImmersedThermalVoxel3 { cell: 0, ..a }],
        vec![ImmersedThermalVoxel3 {
            conductivity: 0.0,
            ..a
        }],
        vec![ImmersedThermalVoxel3 {
            contact_resistance: -0.1,
            ..a
        }],
        vec![],
    ] {
        assert!(
            owner
                .set_immersed_thermal_voxels(&bad, Default::default())
                .is_err()
        );
    }
    for controls in [
        ImmersedVoxelControls3 {
            maximum_geometry_bytes: 0,
            ..Default::default()
        },
        ImmersedVoxelControls3 {
            maximum_body_storage_bytes: 0,
            ..Default::default()
        },
        ImmersedVoxelControls3 {
            conduction: WallConductionControls {
                max_links: 0,
                ..Default::default()
            },
            ..Default::default()
        },
    ] {
        assert!(
            owner
                .set_immersed_thermal_voxels(&[a, b], controls)
                .is_err()
        );
    }
    assert_eq!([owner.immersed_body(0), owner.immersed_body(1)], accepted);
    assert_eq!(owner.populations(), populations);
}
