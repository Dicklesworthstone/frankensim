//! Original thermal adhesion and sigma(T) meet the actual graph/Guo owner.

use super::*;
use fs_lbm::material::free_surface::enthalpy::thermocapillary::HeightGraphThermocapillarity;

fn graph_sources(hot_work: f64) -> (MaterialCard, InterfaceSystemCard) {
    let base = source_card(0.004, false, 500.0);
    let plateau = |hot| PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        // Real constant source intervals, with room for the existing thermal
        // ghost interpolation and small advective mixing at each endpoint.
        knots: vec![(300.0, 0.004), (406.0, 0.004), (409.0, hot), (500.0, hot)],
        dims: SURFACE_TENSION_DIMS,
    };
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == SURFACE_TENSION_PROPERTY {
            claim.value = plateau(0.0044);
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    let liquid = MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap();
    let solid = body::wall_card(0.5, 500.0);
    let original = equilibrium_interface(
        &liquid,
        &solid,
        EquilibriumAngleOrigin::YoungDupre,
        0.0,
        500.0,
    );
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        claim.value = plateau(hot_work);
        claim.interpolation = InterpolationPolicy::LinearInside;
        claims.insert_claim(claim).unwrap();
    }
    let interface = InterfaceSystemCard::assemble(
        original.surface_a().clone(),
        original.surface_b().clone(),
        original.context().clone(),
        claims,
        vec![],
    )
    .unwrap();
    (liquid, interface)
}

fn graph_cells(slope: f64) -> Vec<MaterialFreeSurfaceInitialCell3> {
    let mut initial = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                // Integrate an affine plane through each cell, including a
                // genuine two-cell partial interface band when slope = 0.3.
                let fill = (0..8)
                    .map(|i| {
                        let u = x as f64 + (i as f64 + 0.5) / 8.0;
                        let h = if slope == 0.0 {
                            3.5
                        } else {
                            4.0 + slope * (u - 4.0)
                        };
                        (h - z as f64).clamp(0.0, 1.0)
                    })
                    .sum::<f64>()
                    / 8.0;
                if fill > 0.0 {
                    initial.push(MaterialFreeSurfaceInitialCell3 {
                        coordinates: [x, y, z],
                        fill_fraction: fill,
                        velocity_m_s: [0.0; 3],
                    });
                }
            }
        }
    }
    initial
}

fn warm(p: [usize; 3]) -> bool {
    p[0] == 0 && (2..=5).contains(&p[1])
}

fn with_graph<R>(
    liquid: &MaterialCard,
    interface: &InterfaceSystemCard,
    initial: &[MaterialFreeSurfaceInitialCell3],
    initially_hot: bool,
    run: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let source = liquid_source(liquid, false);
    let law = LiquidTemperatureLaw::try_from_material(
        liquid,
        source.material().clone(),
        "T",
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let chart = caloric(liquid);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut thermal = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
            .unwrap()
            .free_surface3_from_cells_with_height_graph_thermocapillarity(
                domain(),
                initial,
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                HeightGraphThermocapillarity {
                    normal_axis: 2,
                    outward_sign: 1,
                    maximum_slope: 0.6,
                    minimum_surface_fill: 0.01,
                },
                None,
                None,
                |p| {
                    if initially_hot && warm(p) {
                        240000.0
                    } else {
                        200000.0
                    }
                },
                cx,
                &gate,
            )
            .unwrap()
            .with_equilibrium_planar_wetting(
                PlanarWall3::XMin,
                interface,
                EquilibriumAngleOrigin::YoungDupre,
                MaterialPropertySelection::SingleClaimOnly,
                64,
                &gate,
            )
            .unwrap();
        run(&mut thermal, &gate)
    })
}

fn heat_cells(
    thermal: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    initial: &[MaterialFreeSurfaceInitialCell3],
    selected: impl Fn([usize; 3]) -> bool,
) -> Vec<FreeSurfaceCellHeat> {
    initial
        .iter()
        .filter(|cell| selected(cell.coordinates))
        .map(|cell| {
            let p = cell.coordinates;
            FreeSurfaceCellHeat {
                coordinates: p,
                energy_j: 40000.0
                    * thermal
                        .flow()
                        .region_mass_kg(|x, y, z| [x, y, z] == p)
                        .unwrap(),
            }
        })
        .collect()
}

#[test]
fn g1_graph_wetting_uses_current_work_and_tension_before_the_same_native_impulse() {
    let initial = graph_cells(0.0);
    let run = |hot_work, initially_hot| {
        let (liquid, interface) = graph_sources(hot_work);
        with_graph(
            &liquid,
            &interface,
            &initial,
            initially_hot,
            |thermal, gate| {
                let work_id = thermal
                    .wetting_source()
                    .unwrap()
                    .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected;
                let sigma_id = thermal
                    .material_at([0, 3, 3])
                    .unwrap()
                    .material()
                    .property(SURFACE_TENSION_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected;
                let before_h = thermal.total_enthalpy_j().unwrap();
                let heat = if initially_hot {
                    vec![]
                } else {
                    heat_cells(thermal, &initial, warm)
                };
                let supplied = heat.iter().map(|q| q.energy_j).sum::<f64>();
                let report = thermal.step_with_heat(&heat, gate).unwrap();
                near(report.external_heat_j, supplied, 2e-12);
                near(
                    thermal.total_enthalpy_j().unwrap() - before_h,
                    supplied,
                    2e-10,
                );
                near(report.mass_before_kg, report.mass_after_kg, 2e-12);
                assert_eq!(report.viscous_heat_j, 0.0);
                let material = thermal.material_at([0, 3, 3]).unwrap();
                let sigma = material
                    .material()
                    .property(SURFACE_TENSION_PROPERTY)
                    .unwrap();
                let sample = &thermal.wetting_samples().unwrap()[&[0, 3, 3]];
                let work = sample
                    .material
                    .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY)
                    .unwrap();
                assert_eq!(work.answer().receipt.selected, work_id);
                assert_eq!(sigma.answer().receipt.selected, sigma_id);
                assert_eq!(
                    sample.material.state_coordinate("T"),
                    Some(temperature(material))
                );
                assert!(temperature(material) > 409.0);
                near(sigma.value_si(), 0.0044, 2e-14);
                near(work.value_si(), hot_work, 2e-14);
                near(
                    sample.applied_angle_radians,
                    (hot_work / 0.0044 - 1.0).acos(),
                    2e-13,
                );
                // rho=1000 and dx=dt=0.001 give identical SI/native sigma values.
                near(
                    thermal
                        .flow()
                        .lattice()
                        .surface_tension_at(0, 3, 3)
                        .unwrap(),
                    0.0044,
                    2e-14,
                );
                for (&p, sample) in thermal.wetting_samples().unwrap() {
                    if p[1] == 0 || p[1] == 7 {
                        assert_eq!(sample.applied_angle_radians, core::f64::consts::FRAC_PI_2);
                    }
                }
                let patch = thermal
                    .height_graph_interface_patches(gate)
                    .unwrap()
                    .into_iter()
                    .find(|patch| patch.liquid[0] == 0 && patch.liquid[1] == 3)
                    .unwrap();
                near(
                    patch.outward_normal[0],
                    sample.applied_angle_radians.cos(),
                    2e-13,
                );
                let force = thermal.flow().lattice().force_density_at(0, 3, 3).unwrap();
                assert!(force.iter().any(|component| component.abs() > 1e-12));
                near(
                    (0..3).map(|i| force[i] * patch.outward_normal[i]).sum(),
                    0.0,
                    2e-13,
                );
                thermal.flow().lattice().grid().state_bits()
            },
        )
    };
    let heated = run(0.0048, false);
    let warm_start = run(0.0048, true);
    let largest_difference = |a: &[u64], b: &[u64]| {
        assert_eq!(a.len(), b.len());
        a.iter()
            .zip(b)
            .map(|(&x, &y)| (f64::from_bits(x) - f64::from_bits(y)).abs())
            .fold(0.0_f64, f64::max)
    };
    // Both candidates start their native interval at the SAME hot material,
    // geometry and physical rest. A lagged cold wall angle gives a different
    // graph metric and axial Marangoni force during this very interval.
    assert!(largest_difference(&heated, &warm_start) < 2e-12);
    assert!(largest_difference(&heated, &run(0.0044, false)) > 1e-12);
}

#[test]
fn g1_heated_neutral_transverse_corners_use_sourced_adhesion_in_native_momentum() {
    let initial = graph_cells(0.0);
    let run = |hot_work| {
        let (liquid, interface) = graph_sources(hot_work);
        with_graph(&liquid, &interface, &initial, false, |thermal, gate| {
            let source = thermal.wetting_source().unwrap()
                .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY).unwrap()
                .answer().receipt.selected;
            let before_h = thermal.total_enthalpy_j().unwrap();
            let before_mass = thermal.flow().ledger_mass_kg().unwrap();
            // Heat the complete sourced side wall, including BOTH junctions
            // with the independently neutral Cartesian transverse walls.
            let heat = heat_cells(thermal, &initial, |p| p[0] == 0);
            let supplied = heat.iter().map(|entry| entry.energy_j).sum::<f64>();
            let report = thermal.step_with_heat(&heat, gate).unwrap();
            near(report.external_heat_j, supplied, 2e-12);
            near(thermal.total_enthalpy_j().unwrap() - before_h, supplied, 2e-10);
            near(thermal.flow().ledger_mass_kg().unwrap(), before_mass, 2e-12);
            assert_eq!(report.viscous_heat_j, 0.0);
            let patches = thermal.height_graph_interface_patches(gate).unwrap();
            for edge in [0, 7] {
                let patch = patches.iter().find(|patch| {
                    patch.liquid[0] == 0 && patch.liquid[1] == edge
                }).unwrap();
                let sample = &thermal.wetting_samples().unwrap()[&patch.liquid];
                let material = thermal.material_at(patch.liquid).unwrap();
                let sigma = material.material().property(SURFACE_TENSION_PROPERTY).unwrap();
                let work = sample.material.property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY).unwrap();
                assert_eq!(work.answer().receipt.selected, source);
                assert_eq!(sample.material.state_coordinate("T"), Some(temperature(material)));
                assert!(temperature(material) > 409.0);
                near(work.value_si(), hot_work, 2e-14);
                near(sigma.value_si(), 0.0044, 2e-14);
                // Independent Young--Dupre and orthogonal-wall normal
                // constraints determine the metric, not the source sampler.
                let cosine = hot_work / 0.0044 - 1.0;
                let sine = (1.0 - cosine * cosine).sqrt();
                near(sample.applied_angle_radians, cosine.acos(), 2e-13);
                near(patch.outward_normal[0], cosine, 2e-13);
                assert_eq!(patch.outward_normal[1], 0.0);
                near(patch.outward_normal[2], sine, 2e-13);
                near(patch.area_m2, 1e-6 / sine, 2e-13);
                let [x, y, z] = patch.liquid;
                let force = thermal.flow().lattice().force_density_at(x, y, z).unwrap();
                assert!(force.iter().any(|component| component.abs() > 1e-12));
                assert!((0..3).map(|d| force[d] * patch.outward_normal[d])
                    .sum::<f64>().abs() < 2e-13);
            }
            thermal.flow().lattice().grid().state_bits()
        })
    };
    let neutral = run(0.0044);
    let adhesive = run(0.0048);
    assert_eq!(neutral.len(), adhesive.len());
    let difference = neutral.iter().zip(&adhesive)
        .map(|(&a, &b)| (f64::from_bits(a) - f64::from_bits(b)).abs())
        .fold(0.0_f64, f64::max);
    assert!(difference > 1e-12, "corner adhesion must reach the actual native populations");
}

#[test]
fn g3_graph_thermal_corner_source_band_and_slope_refusals_leave_every_owner_unchanged() {
    for case in 0..3 {
        let initial = graph_cells(if case == 1 { 0.3 } else { 0.0 });
        let (liquid, interface) = graph_sources(match case {
            // A supported corner must still refuse loss of finite-angle
            // partial wetting: W > 2*sigma cannot be replaced by a clamped angle.
            0 => 0.009,
            2 => 0.0075,
            _ => 0.0048,
        });
        with_graph(&liquid, &interface, &initial, false, |thermal, gate| {
            let before = thermal.flow().lattice().grid().state_bits();
            let samples = thermal.wetting_samples().unwrap().clone();
            let source = thermal.wetting_source().unwrap().clone();
            let time = thermal.physical_time_s();
            let energy = thermal.total_enthalpy_j().unwrap();
            let fields = initial
                .iter()
                .map(|cell| {
                    let [x, y, z] = cell.coordinates;
                    (
                        thermal.specific_enthalpy_at(cell.coordinates),
                        thermal.flow().lattice().surface_tension_at(x, y, z),
                        thermal.flow().lattice().force_density_at(x, y, z),
                    )
                })
                .collect::<Vec<_>>();
            let heat = heat_cells(thermal, &initial, |p| match case {
                0 => p[0] == 0 && p[1] == 0,
                1 => p == [0, 3, 3],
                _ => warm(p),
            });
            let error = thermal.step_with_heat(&heat, gate).unwrap_err().to_string();
            let reason = match case {
                0 => "partial-wetting equilibrium",
                1 => "occupied interface band",
                _ => "slope envelope",
            };
            assert!(error.contains(reason), "expected {reason}, got {error}");
            assert_eq!(thermal.flow().lattice().grid().state_bits(), before);
            assert_eq!(thermal.wetting_samples().unwrap(), &samples);
            assert_eq!(thermal.wetting_source().unwrap(), &source);
            assert_eq!(thermal.physical_time_s(), time);
            assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
            assert_eq!(
                initial
                    .iter()
                    .map(|cell| {
                        let [x, y, z] = cell.coordinates;
                        (
                            thermal.specific_enthalpy_at(cell.coordinates),
                            thermal.flow().lattice().surface_tension_at(x, y, z),
                            thermal.flow().lattice().force_density_at(x, y, z),
                        )
                    })
                    .collect::<Vec<_>>(),
                fields
            );
        })
    }
}
