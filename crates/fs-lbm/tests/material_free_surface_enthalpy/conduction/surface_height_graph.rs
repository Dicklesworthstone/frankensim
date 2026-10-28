//! Shared curved geometry and solved-interface temperature reach native stress.
//! Synthetic sources check conservation and units, not experimental calibration.
use super::*;
use fs_lbm::material::free_surface::enthalpy::thermocapillary::HeightGraphThermocapillarity;

fn height(u: f64, v: f64) -> f64 {
    // A smooth convex cap with a mixed slope, crossing lattice layers. Its
    // positive Laplace loading stays inside the native graph fill envelope;
    // a concave depression can enter the VOF owner's overfill hysteresis.
    2.85 + 0.1 * u + 0.08 * v - 0.002 * u * u - 0.001 * v * v + 0.001 * u * v
}

fn with_graph<R>(
    card: &MaterialCard,
    axis: usize,
    sign: i8,
    run: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let source = liquid_source(card, false);
    let law = LiquidTemperatureLaw::try_from_material(
        card,
        source.material().clone(),
        "T",
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let chart = caloric(card);
    let tangent: Vec<_> = (0..3).filter(|&d| d != axis).collect();
    let mut cells = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                let p = [x, y, z];
                let h = height(p[tangent[0]] as f64, p[tangent[1]] as f64);
                let fill = if sign > 0 {
                    h - p[axis] as f64
                } else {
                    p[axis] as f64 + 1.0 - h
                };
                if fill > 0.0 {
                    cells.push(MaterialFreeSurfaceInitialCell3 {
                        coordinates: p,
                        fill_fraction: fill.min(1.0),
                        velocity_m_s: [0.0; 3],
                    });
                }
            }
        }
    }
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut thermal = plan_free_surface_from_material(&source, DX, DT, 0.1)
            .unwrap()
            .free_surface3_from_cells_with_height_graph_thermocapillarity(
                domain(),
                &cells,
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                HeightGraphThermocapillarity {
                    normal_axis: axis,
                    outward_sign: sign,
                    maximum_slope: 0.6,
                    minimum_surface_fill: 0.001,
                },
                None,
                None,
                |_| 200000.0,
                cx,
                &gate,
            )
            .unwrap();
        run(&mut thermal, &gate)
    })
}

fn surface_controls() -> FreeSurfaceConductionControls {
    FreeSurfaceConductionControls {
        maximum_faces: 80,
        ..controls()
    }
}

#[test]
fn g1_curved_heat_patch_geometry_matches_quadratic_graph_in_all_orientations() {
    let card = thermal_card(-1e-5, 500.0);
    for axis in 0..3 {
        for sign in [-1, 1] {
            with_graph(&card, axis, sign, |thermal, gate| {
                let tangent: Vec<_> = (0..3).filter(|&d| d != axis).collect();
                let patches = thermal.height_graph_interface_patches(gate).unwrap();
                assert_eq!(patches.len(), 64);
                for patch in patches {
                    let u = patch.liquid[tangent[0]] as f64;
                    let v = patch.liquid[tangent[1]] as f64;
                    let slopes = [0.1 - 0.004 * u + 0.001 * v, 0.08 - 0.002 * v + 0.001 * u];
                    let metric = (1.0 + slopes[0] * slopes[0] + slopes[1] * slopes[1]).sqrt();
                    near(patch.centroid_m[axis], DX * height(u, v), 3e-13);
                    near(patch.centroid_m[tangent[0]], DX * (u + 0.5), 3e-13);
                    near(patch.centroid_m[tangent[1]], DX * (v + 0.5), 3e-13);
                    near(patch.area_m2, DX * DX * metric, 3e-13);
                    near(patch.outward_normal[axis], sign as f64 / metric, 3e-13);
                    for d in 0..2 {
                        near(
                            patch.outward_normal[tangent[d]],
                            -sign as f64 * slopes[d] / metric,
                            3e-13,
                        );
                    }
                    assert_eq!(
                        patch.gas[axis] as isize - patch.liquid[axis] as isize,
                        sign as isize
                    );
                    assert!(
                        thermal
                            .flow()
                            .lattice()
                            .sample_at(patch.liquid[0], patch.liquid[1], patch.liquid[2],)
                            .unwrap()
                            .fill_fraction
                            > 0.0
                    );
                    assert!(
                        thermal
                            .flow()
                            .lattice()
                            .sample_at(patch.gas[0], patch.gas[1], patch.gas[2],)
                            .is_none_or(|sample| sample.fill_fraction == 0.0)
                    );
                }
            });
        }
    }
}

#[test]
fn g1_curved_solved_surface_shares_native_sigma_geometry_and_conservative_heat() {
    let slope = -1e-5;
    let card = thermal_card(slope, 500.0);
    let gas = gas_source::model(0.0, 250000.0, 1000.0);
    with_graph(&card, 2, 1, |thermal, gate| {
        let patches = thermal.height_graph_interface_patches(gate).unwrap();
        let surfaces: Vec<_> = patches
            .iter()
            .map(|patch| {
                assert!(patch.area_m2 > DX * DX);
                patch.with_exchange([0.25 * DX, 0.5 * DX], 0.0, 400.0)
            })
            .collect();
        let points: Vec<_> = patches
            .iter()
            .map(|patch| gas_source::point(350.0 + 10.0 * patch.liquid[0] as f64, 100000.0))
            .collect();
        let boundaries: Vec<_> = surfaces
            .iter()
            .zip(&points)
            .map(|(&surface, context)| FreeSurfaceInterfaceContact {
                surface,
                model: &gas,
                context,
            })
            .collect();
        // Five Cartesian contacts leave only dx² on this surface cell in the
        // old voxel budget. The shared graph area correction must expose its
        // true, larger curved patch without weakening unrelated face limits.
        let center = patches
            .iter()
            .find(|patch| patch.liquid[..2] == [3, 3])
            .unwrap()
            .liquid;
        let faces: Vec<_> = [[-1, 0, 0], [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1]]
            .into_iter()
            .map(|offset| FreeSurfaceConductionFace {
                cells: [
                    center,
                    std::array::from_fn(|d| (center[d] as isize + offset[d]) as usize),
                ],
                area_m2: DX * DX,
                distances_m: [0.5 * DX; 2],
            })
            .collect();
        let energy = thermal.total_enthalpy_j().unwrap();
        let mass = thermal.flow().ledger_mass_kg().unwrap();
        let report = thermal
            .step_with_surface_exchange(
                &faces,
                &[],
                surface_controls(),
                &mut [],
                &[],
                &[],
                &boundaries,
                FreeSurfaceInterfaceSolve::default(),
                false,
                gate,
            )
            .unwrap()
            .conduction;
        let exchange = report.surface.as_ref().unwrap();
        assert_eq!(exchange.native_endpoint.len(), 64);
        assert_eq!(exchange.heat_from_radiation_j, 0.0);
        near(report.flow.external_heat_j, exchange.heat_from_gas_j, 2e-9);
        near(
            thermal.total_enthalpy_j().unwrap() - energy,
            exchange.heat_from_gas_j,
            2e-8,
        );
        near(thermal.flow().ledger_mass_kg().unwrap(), mass, 2e-12);

        let live = thermal.height_graph_interface_patches(gate).unwrap();
        let sample = |x, y| {
            exchange
                .native_endpoint
                .iter()
                .find(|sample| sample.surface.face.liquid[..2] == [x, y])
                .unwrap()
        };
        let sigma = |x, y| sample(x, y).surface_tension.as_ref().unwrap().value_si();
        for patch in &live {
            let resolved = sample(patch.liquid[0], patch.liquid[1]);
            assert_eq!(resolved.surface.face.liquid, patch.liquid);
            assert_eq!(resolved.surface.face.gas, patch.gas);
            near(resolved.surface.face.area_m2, patch.area_m2, 2e-13);
            let gl = resolved.liquid.thermal_conductivity_w_m_k() * patch.area_m2 / (0.25 * DX);
            let gg = resolved.gas.state().thermal_conductivity * patch.area_m2 / (0.5 * DX);
            let expected_ts =
                (gl * resolved.liquid_temperature_k + gg * resolved.gas_temperature_k) / (gl + gg);
            near(resolved.surface_temperature_k, expected_ts, 2e-11);
            let physical_sigma = 0.01 + slope * (expected_ts - 400.0);
            near(
                resolved.surface_tension.as_ref().unwrap().value_si(),
                physical_sigma,
                2e-12,
            );
            near(
                thermal
                    .flow()
                    .lattice()
                    .surface_tension_at(patch.liquid[0], patch.liquid[1], patch.liquid[2])
                    .unwrap(),
                physical_sigma * DT * DT / (1000.0 * DX.powi(3)),
                2e-12,
            );
        }
        let patch = live
            .iter()
            .find(|patch| patch.liquid[..2] == [3, 3])
            .unwrap();
        let ambient = [
            (sigma(4, 3) - sigma(2, 3)) / (2.0 * DX),
            (sigma(3, 4) - sigma(3, 2)) / (2.0 * DX),
            0.0,
        ];
        let normal = patch.outward_normal;
        let normal_gradient = (0..3).map(|d| ambient[d] * normal[d]).sum::<f64>();
        let mut impulse = [0.0; 3];
        for z in 0..8 {
            if let (Some(cell), Some(force)) = (
                thermal.flow().lattice().sample_at(3, 3, z),
                thermal.flow().lattice().force_density_at(3, 3, z),
            ) {
                for d in 0..3 {
                    impulse[d] += cell.fill_fraction * force[d];
                }
            }
        }
        for d in 0..3 {
            let physical_force = patch.area_m2 * (ambient[d] - normal_gradient * normal[d]);
            let expected = physical_force * DT * DT / (1000.0 * DX.powi(4));
            assert!((impulse[d] - expected).abs() < 2e-12 * expected.abs().max(1e-6));
        }
        assert!(impulse[0] < 0.0);
        assert!((0..3).map(|d| impulse[d] * normal[d]).sum::<f64>().abs() < 1e-16);
        assert!(
            (sample(3, 3).liquid_temperature_k - sample(3, 3).surface_temperature_k).abs() > 0.1
        );
    });
}

#[test]
fn g4_curved_heat_geometry_refusals_leave_flow_and_enthalpy_unchanged() {
    let gas = gas_source::model(0.0, 250000.0, 1000.0);
    let point = gas_source::point(350.0, 100000.0);
    // Include a true scalar sigma claim with no source axes. Geometry remains
    // mandatory even when interface temperature cannot change that constant.
    for card in [thermal_card(-1e-5, 500.0), source_card(0.01, false, 500.0)] {
        with_graph(&card, 2, 1, |thermal, gate| {
            let boundaries: Vec<_> = thermal
                .height_graph_interface_patches(gate)
                .unwrap()
                .into_iter()
                .map(|patch| FreeSurfaceInterfaceContact {
                    surface: patch.with_exchange([0.25 * DX, 0.5 * DX], 0.0, 400.0),
                    model: &gas,
                    context: &point,
                })
                .collect();
            let before = thermal.flow().lattice().grid().state_bits();
            let energy = thermal.total_enthalpy_j().unwrap().to_bits();
            for kind in 0..4 {
                let mut invalid = boundaries.clone();
                match kind {
                    0 => {
                        invalid.pop();
                    }
                    1 => invalid[0] = invalid[1],
                    2 => invalid[0].surface.face.area_m2 = DX * DX,
                    _ => invalid[0].surface.face.distances_m[0] = 0.0,
                }
                assert!(
                    thermal
                        .step_with_surface_exchange(
                            &[],
                            &[],
                            surface_controls(),
                            &mut [],
                            &[],
                            &[],
                            &invalid,
                            FreeSurfaceInterfaceSolve::default(),
                            false,
                            gate,
                        )
                        .is_err()
                );
                assert_eq!(before, thermal.flow().lattice().grid().state_bits());
                assert_eq!(energy, thermal.total_enthalpy_j().unwrap().to_bits());
            }
        });
    }
}
