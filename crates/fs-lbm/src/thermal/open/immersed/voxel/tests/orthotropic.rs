use super::*;
use crate::thermal::advection::D3Q7;

fn orthotropic(cell: usize) -> OrthotropicThermalVoxel3 {
    OrthotropicThermalVoxel3 {
        cell,
        temperature: 1.4,
        heat_capacity: 2.5,
        conductivity: [0.2, 1.5, 5.0],
        contact_resistance: [0.1, 0.3, 1.0],
        heat_per_step: 0.02,
    }
}

fn carrier(dims: [usize; 3], voxels: &[OrthotropicThermalVoxel3]) -> ThermalOpenField3 {
    let mask: Vec<_> = voxels
        .iter()
        .map(|v| voxel(v.cell, v.temperature, v.heat_capacity, 1.0))
        .collect();
    field(dims, &mask)
}

#[test]
fn g1_each_normal_uses_its_own_half_cell_and_actual_outgoing_trace() {
    let voxel = orthotropic(13);
    let mut owner = carrier([3; 3], &[voxel]);
    let neighbors = [[12, 14], [10, 16], [4, 22]];
    let traces = [[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]];
    for axis in 0..3 {
        for side in 0..2 {
            owner.populations[neighbors[axis][side]] =
                D3Q7.equilibrium(traces[axis][side], [0.0; 3]);
        }
    }
    owner
        .set_immersed_orthotropic_thermal_voxels(&[voxel], Default::default())
        .unwrap();
    let k: [f64; 3] = std::array::from_fn(|axis| {
        1.0 / (4.0 + 0.5 / voxel.conductivity[axis] + voxel.contact_resistance[axis])
    });
    let rhs = voxel.heat_capacity * voxel.temperature
        + voxel.heat_per_step
        + (0..3)
            .map(|axis| k[axis] * (traces[axis][0] + traces[axis][1]))
            .sum::<f64>();
    let expected = rhs / (voxel.heat_capacity + 2.0 * k.iter().sum::<f64>());
    let balance = advance(&mut owner);
    let body = owner.immersed_body(0).unwrap();
    near(body.solid.temperature, expected, 3e-15);
    for axis in 0..3 {
        near(
            body.surface_resistance[axis],
            0.5 / voxel.conductivity[axis] + voxel.contact_resistance[axis],
            0.0,
        );
        for side in 0..2 {
            let incoming = 1 + 2 * axis + usize::from(side == 0);
            let reflected = owner.populations()[neighbors[axis][side]][incoming];
            near(
                reflected - 0.125 * traces[axis][side],
                k[axis] * (expected - traces[axis][side]),
                3e-16,
            );
        }
    }
    let exchange = body.last_exchange.unwrap();
    near(
        exchange.stored_change + exchange.heat_in,
        voxel.heat_per_step,
        3e-15,
    );
    near(balance.combined_defect, 0.0, 2e-14);
}

#[test]
fn g1_rotating_solid_tensor_and_geometry_preserves_the_independent_pair_balance() {
    for parallel_faces in [1usize, 2] {
        let mut reference: Option<[f64; 2]> = None;
        for normal in 0..3 {
            let mut dims = [5; 3];
            dims[normal] = if parallel_faces == 2 { 2 } else { 5 };
            let stride = [1, dims[0], dims[0] * dims[1]];
            let mut coordinate = [2; 3];
            coordinate[normal] = if parallel_faces == 2 { 0 } else { 2 };
            let first = (0..3).map(|axis| coordinate[axis] * stride[axis]).sum();
            let mut a = orthotropic(first);
            let mut b = orthotropic(first + stride[normal]);
            for local in 0..3 {
                let axis = (normal + local) % 3;
                a.conductivity[axis] = [0.4, 2.0, 4.0][local];
                b.conductivity[axis] = [1.6, 0.8, 3.0][local];
                a.contact_resistance[axis] = [0.2, 0.4, 0.7][local];
                b.contact_resistance[axis] = [0.3, 0.1, 0.5][local];
            }
            a.heat_per_step = 0.06;
            b.temperature = 0.4;
            b.heat_capacity = 4.0;
            b.heat_per_step = -0.01;
            let mut owner = carrier(dims, &[a, b]);
            let map = owner
                .set_immersed_orthotropic_thermal_voxels(&[b, a], Default::default())
                .unwrap();
            assert_eq!(map.conduction_links(), 1);
            let fluid_k = |v: OrthotropicThermalVoxel3| -> f64 {
                (0..3)
                    .map(|axis| {
                        let faces = 2 - if axis == normal { parallel_faces } else { 0 };
                        faces as f64
                            / (4.0 + 0.5 / v.conductivity[axis] + v.contact_resistance[axis])
                    })
                    .sum()
            };
            let ka = fluid_k(a);
            let kb = fluid_k(b);
            let g = parallel_faces as f64
                / (0.5 / a.conductivity[normal] + 0.5 / b.conductivity[normal]);
            let aa = a.heat_capacity + ka + g;
            let bb = b.heat_capacity + kb + g;
            let ra = a.heat_capacity * a.temperature + 0.2 * ka + a.heat_per_step;
            let rb = b.heat_capacity * b.temperature + 0.2 * kb + b.heat_per_step;
            let determinant = aa * bb - g * g;
            let expected = [
                (bb * ra + g * rb) / determinant,
                (aa * rb + g * ra) / determinant,
            ];
            let balance = advance(&mut owner);
            for node in 0..2 {
                let body = owner.immersed_body(node).unwrap();
                near(body.solid.temperature, expected[node], 5e-15);
                assert_eq!(body.interface_links, 6 - parallel_faces);
                near(body.last_exchange.unwrap().defect(), 0.0, 5e-15);
                if let Some(reference) = reference {
                    near(body.solid.temperature, reference[node], 3e-15);
                }
            }
            near(
                owner
                    .immersed_body(0)
                    .unwrap()
                    .last_exchange
                    .unwrap()
                    .conducted_in,
                g * (expected[1] - expected[0]),
                3e-15,
            );
            near(balance.combined_defect, 0.0, 3e-14);
            reference = Some(expected);
        }
    }
}

#[test]
fn g3_isotropic_declarations_match_exactly_and_uniform_orthotropic_equilibrium_is_retained() {
    let a = voxel(27, 1.4, 2.0, 0.7);
    let b = voxel(28, 0.4, 5.0, 1.2);
    let mut isotropic = field([5, 4, 3], &[a, b]);
    let mut directional = field([5, 4, 3], &[a, b]);
    let original = isotropic
        .set_immersed_thermal_voxels(&[b, a], Default::default())
        .unwrap();
    let extended = directional
        .set_immersed_orthotropic_thermal_voxels(&[a.into(), b.into()], Default::default())
        .unwrap();
    assert_eq!(original, extended);
    for _ in 0..20 {
        assert_eq!(advance(&mut isotropic), advance(&mut directional));
        assert_eq!(isotropic.populations(), directional.populations());
        for node in 0..2 {
            assert_eq!(
                isotropic.immersed_body(node),
                directional.immersed_body(node)
            );
        }
    }
    for capacity in [f64::from_bits(1), 1e-24, 1.0, 1e24, f64::MAX] {
        let mut voxel = orthotropic(13);
        voxel.temperature = 0.5;
        voxel.heat_capacity = capacity;
        voxel.heat_per_step = 0.0;
        let mut owner = carrier([3; 3], &[voxel]);
        for cell in 0..27 {
            if !owner.solid[cell] {
                owner.populations[cell] = D3Q7.equilibrium(0.5, [0.0; 3]);
            }
        }
        owner
            .set_immersed_orthotropic_thermal_voxels(&[voxel], Default::default())
            .unwrap();
        advance(&mut owner);
        let body = owner.immersed_body(0).unwrap();
        assert_eq!(body.solid.temperature, 0.5);
        assert_eq!(body.last_exchange.unwrap().heat_in, 0.0);
        assert_eq!(body.last_exchange.unwrap().stored_change, 0.0);
    }
}

#[test]
fn g0_bad_directional_data_and_late_cancel_preserve_actual_native_history() {
    let v = orthotropic(13);
    let mut actual = carrier([3; 3], &[v]);
    let mut control = carrier([3; 3], &[v]);
    for owner in [&mut actual, &mut control] {
        owner
            .set_immersed_orthotropic_thermal_voxels(&[v], Default::default())
            .unwrap();
        advance(owner);
    }
    let accepted = actual.immersed_body(0);
    let populations = actual.populations().to_vec();
    let history = format!("{:?}", actual.history);
    let receipt = actual.last_balance();
    for axis in 0..3 {
        for invalid in [0.0, -0.1, f64::INFINITY, f64::NAN, f64::from_bits(1)] {
            let mut bad = v;
            bad.conductivity[axis] = invalid;
            assert!(
                actual
                    .set_immersed_orthotropic_thermal_voxels(&[bad], Default::default())
                    .is_err()
            );
        }
        let mut bad = v;
        bad.contact_resistance[axis] = -1.0;
        assert!(
            actual
                .set_immersed_orthotropic_thermal_voxels(&[bad], Default::default())
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
                max_factor_bytes: 1,
                ..Default::default()
            },
            ..Default::default()
        },
    ] {
        assert!(
            actual
                .set_immersed_orthotropic_thermal_voxels(&[v], controls)
                .is_err()
        );
    }
    let heat = [0.0; 27];
    let mut polls = 0;
    actual
        .begin_step::<ThermalOpenError3>(&heat, None, || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        actual.begin_step::<ThermalOpenError3>(&heat, None, || {
            count += 1;
            if count == polls {
                Err(ThermalOpenError3::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(ThermalOpenError3::Cancelled)
    ));
    assert_eq!(actual.immersed_body(0), accepted);
    assert_eq!(actual.populations(), populations);
    assert_eq!(format!("{:?}", actual.history), history);
    assert_eq!(actual.last_balance(), receipt);
    assert_eq!(advance(&mut actual), advance(&mut control));
    assert_eq!(actual.populations(), control.populations());
    assert_eq!(actual.immersed_body(0), control.immersed_body(0));
}
