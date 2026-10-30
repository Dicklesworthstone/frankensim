//! G1/G3/G4: reference-mass transfer retains actual source density and shear.
//! Synthetic material curves isolate native force and storage, not calibration.

use super::*;
use crate::material::free_surface::enthalpy::FreeSurfaceCellHeat;
use crate::material::free_surface::enthalpy::rheology::RheologySolve;
use crate::material::transport::boussinesq::BoussinesqConfig;
use fs_qty::Frequency;

fn card(density: Vec<(f64, f64)>, rheology: bool) -> MaterialCard {
    let original = if rheology {
        super::rheology::source_card(0.02, 100.0)
    } else {
        thermal_card(1000.0, 1e-4)
    };
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == "density" {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: density.clone(),
                dims: claim.key.dims(),
            };
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

fn reference_point(rheology: bool) -> QueryPoint {
    let point = point().with("T", 300.0).unwrap();
    if rheology {
        point
            .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 1.5)
            .unwrap()
    } else {
        point
    }
}

fn config(gravity: [f64; 3]) -> BoussinesqConfig<3> {
    BoussinesqConfig {
        gravity_m_s2: gravity,
        max_relative_density_departure: 0.03,
        max_relative_specific_heat_departure: 0.0,
    }
}

fn bounds() -> TetEnthalpyTransferLimits {
    let mut limits = thermal_limits();
    limits.mechanical.max_tetrahedra = 16;
    limits.max_source_queries = 512;
    limits.enthalpy_tolerance_j = 1e-9;
    limits
}

fn rheology() -> RheologySolve {
    RheologySolve {
        maximum_iterations: 128,
        relaxation: 0.5,
        absolute_relaxation_tolerance: 1e-13,
        relative_relaxation_tolerance: 1e-13,
    }
}

#[test]
fn source_density_changes_actual_gravity_after_conservative_fragment_transfer() {
    let run = |hot_density: f64, buoyant: bool| {
        let card = card(vec![(300.0, 1000.0), (500.0, hot_density)], false);
        let point = reference_point(false);
        let law = LiquidTemperatureLaw::try_new(
            &card,
            &point,
            "T",
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let surface = resolve_liquid_surface_state(
            &card,
            &point,
            LiquidViscosityBasis::Dynamic,
            LiquidCapillarity::Neglected,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let caloric = chart(&card, 1000.0);
        let (mut p, mut t, _) = cube([0.0; 3]);
        for point in &mut p {
            *point = point.map(|x| 2.0 * x);
        }
        let second: Vec<_> = p.iter().map(|p| [p[0] + 0.016, p[1], p[2]]).collect();
        p.extend(second);
        let second: Vec<_> = t.iter().map(|tet| tet.map(|i| i + 8)).collect();
        t.extend(second);
        let v: Vec<_> = (0..16)
            .map(|i| {
                if i < 8 {
                    [0.01, 0.0, 0.0]
                } else {
                    [-0.005, 0.0, 0.0]
                }
            })
            .collect();
        let h: Vec<_> = (0..16)
            .map(|i| if i < 8 { 200000.0 } else { 250000.0 })
            .collect();
        let before = (p.clone(), t.clone(), v.clone(), h.clone());
        let gate = CancelGate::new_clock_free();
        let result = with_cx(&gate, |cx| {
            let domain = MaterialFreeSurfaceDomain3 {
                cells: [32; 3],
                max_cells: 32768,
                gravity_m_s2: [0.0, 0.0, -1.0],
                ..domain()
            };
            let scaling = plan_free_surface_from_material(&surface, 0.001, 0.001, 0.1).unwrap();
            if hot_density != 1000.0 {
                assert!(
                    scaling
                        .free_surface3_from_tetrahedra_with_enthalpy(
                            domain,
                            [0.0; 3],
                            source(&p, &t, &v, &h),
                            &law,
                            &caloric,
                            bounds(),
                            Default::default(),
                            cx,
                            &gate,
                        )
                        .is_err()
                );
            }
            let (mut liquid, report) = if buoyant {
                scaling.free_surface3_from_tetrahedra_with_boussinesq_enthalpy(
                    domain,
                    [0.0; 3],
                    source(&p, &t, &v, &h),
                    &law,
                    &caloric,
                    bounds(),
                    Default::default(),
                    config(domain.gravity_m_s2),
                    None,
                    cx,
                    &gate,
                )
            } else {
                scaling.free_surface3_from_tetrahedra_with_enthalpy(
                    domain,
                    [0.0; 3],
                    source(&p, &t, &v, &h),
                    &law,
                    &caloric,
                    bounds(),
                    Default::default(),
                    cx,
                    &gate,
                )
            }
            .unwrap();
            let mass = 1000.0 * 0.008_f64.powi(3);
            let energy = mass * 450000.0;
            assert_eq!(report.mechanical.components.len(), 2);
            for (index, component) in report.mechanical.components.iter().enumerate() {
                close(component.target.mass_kg, mass, 1e-14);
                let expected = scale(
                    if index == 0 {
                        [0.01, 0.0, 0.0]
                    } else {
                        [-0.005, 0.0, 0.0]
                    },
                    mass,
                );
                for axis in 0..3 {
                    close(
                        component.target.momentum_kg_m_s[axis],
                        expected[axis],
                        1e-16,
                    );
                    close(
                        component.target.angular_momentum_kg_m2_s[axis],
                        component.source.angular_momentum_kg_m2_s[axis],
                        1e-18,
                    );
                }
            }
            close(report.target_enthalpy_j, energy, 1e-9);
            close(liquid.flow().ledger_mass_kg().unwrap(), 2.0 * mass, 1e-14);
            // Attaching the actual density force changes its Guo frame only;
            // it must not supply the extra half impulse often hidden here.
            close(
                liquid
                    .flow()
                    .sample_at_si(7, 7, 7)
                    .unwrap()
                    .unwrap()
                    .velocity_m_s[2],
                0.0,
                1e-14,
            );
            let rho_400 = 0.5 * (1000.0 + hot_density);
            close(
                liquid.material_at([7, 7, 7]).unwrap().density_kg_m3(),
                rho_400,
                1e-10,
            );
            liquid.step(&gate).unwrap();
            let first = liquid
                .flow()
                .sample_at_si(7, 7, 7)
                .unwrap()
                .unwrap()
                .velocity_m_s[2];
            close(first, -0.001 * rho_400 / 1000.0, 2e-11);
            let mut heat = Vec::new();
            for z in 0..32 {
                for y in 0..32 {
                    for x in 0..32 {
                        if let Some(cell) = liquid.flow().sample_at_si(x, y, z).unwrap() {
                            if cell.fill_fraction > 0.0 {
                                heat.push(FreeSurfaceCellHeat {
                                    coordinates: [x, y, z],
                                    energy_j: cell.density_kg_m3
                                        * cell.fill_fraction
                                        * 1e-9
                                        * 100000.0,
                                });
                            }
                        }
                    }
                }
            }
            let step = liquid.step_with_heat(&heat, &gate).unwrap();
            let second = liquid
                .flow()
                .sample_at_si(7, 7, 7)
                .unwrap()
                .unwrap()
                .velocity_m_s[2];
            let rho_450 = 1000.0 + 0.75 * (hot_density - 1000.0);
            close(second - first, -0.001 * rho_450 / 1000.0, 2e-11);
            close(
                step.enthalpy_after_j,
                energy + heat.iter().map(|q| q.energy_j).sum::<f64>(),
                1e-8,
            );
            close(step.mass_after_kg, 2.0 * mass, 1e-14);
            close(step.viscous_heat_j, 0.0, 0.0);
            (liquid.flow().lattice().grid().state_bits(), first, second)
        });
        assert_eq!((p, t, v, h), before);
        result
    };
    let frozen = run(1000.0, false);
    let reference = run(1000.0, true);
    assert_eq!(frozen, reference);
    let expanding = run(980.0, true);
    assert!(expanding.1 > reference.1 && expanding.2 > reference.2);
}

#[test]
fn buoyancy_and_rheology_share_source_state_without_replacing_reference_inertia() {
    let run = |hot_density: f64, gravity: [f64; 3]| {
        let card = card(vec![(300.0, 1000.0), (500.0, hot_density)], true);
        let law = LiquidTemperatureLaw::try_new(
            &card,
            &reference_point(true),
            "T",
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let surface = ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
            law.reference_state().material().clone(),
            LiquidViscosityBasis::Dynamic,
            "T",
            "rate",
        )
        .unwrap();
        let caloric = chart(&card, 1000.0);
        let (p, t, _) = cube([0.0; 3]);
        let v: Vec<_> = p.iter().map(|p| [2.0 * (p[1] - 0.004), 0.0, 0.0]).collect();
        let h = vec![200000.0; p.len()];
        let gate = CancelGate::new_clock_free();
        with_cx(&gate, |cx| {
            let scaling = plan_free_surface_from_material(&surface, 0.001, 0.0005, 0.1).unwrap();
            let (mut liquid, report) = scaling
                .free_surface3_from_tetrahedra_with_boussinesq_enthalpy(
                    MaterialFreeSurfaceDomain3 {
                        gravity_m_s2: gravity,
                        ..domain()
                    },
                    [0.0; 3],
                    source(&p, &t, &v, &h),
                    &law,
                    &caloric,
                    bounds(),
                    Default::default(),
                    config(gravity),
                    Some(rheology()),
                    cx,
                    &gate,
                )
                .unwrap();
            for rate in report.source_shear_rates_per_s {
                close(rate, 2.0, 1e-13);
            }
            let energy = liquid.total_enthalpy_j().unwrap();
            let mut heat = 0.0;
            for _ in 0..3 {
                let step = liquid.step_with_viscous_heating(&[], &gate).unwrap();
                for sample in &step.applied {
                    let point = sample.material.material().query();
                    let mu =
                        0.1 - 0.0004 * (point.axes()["T"] - 300.0) + 0.02 * point.axes()["rate"];
                    close(sample.material.dynamic_viscosity_pa_s(), mu, 1e-11);
                    close(
                        sample.heat_j,
                        mu * sample.shear_rate_per_s.powi(2) * sample.liquid_volume_m3 * 0.0005,
                        1e-14,
                    );
                }
                heat += step.flow.viscous_heat_j;
            }
            assert!(heat > 0.0);
            close(liquid.total_enthalpy_j().unwrap(), energy + heat, 1e-8);
            let rho = liquid.material_at([3, 3, 3]).unwrap().density_kg_m3();
            close(
                liquid.flow().lattice().force_density_at(3, 3, 3).unwrap()[2],
                gravity[2] * 0.00025 * rho / 1000.0,
                1e-13,
            );
            (
                liquid.flow().lattice().grid().state_bits(),
                liquid.total_enthalpy_j().unwrap(),
            )
        })
    };
    // The unprojected source kinematic viscosities differ; actual collision
    // must use the same dynamic mu/rho0, including after molecular heat.
    assert_eq!(run(1000.0, [0.0; 3]), run(980.0, [0.0; 3]));
    assert_ne!(run(980.0, [0.0; 3]).0, run(980.0, [0.0, 0.0, -1.0]).0);
}

#[test]
fn original_source_interval_bounds_and_receiving_gravity_refuse_atomically() {
    let run = |density: Vec<(f64, f64)>, gravity: [f64; 3], cancelled: bool| {
        let card = card(density, false);
        let point = reference_point(false);
        let law = LiquidTemperatureLaw::try_new(
            &card,
            &point,
            "T",
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let surface = resolve_liquid_surface_state(
            &card,
            &point,
            LiquidViscosityBasis::Dynamic,
            LiquidCapillarity::Neglected,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let caloric = chart(&card, 1000.0);
        let (p, t, v) = cube([0.0; 3]);
        let h: Vec<_> = p
            .iter()
            .map(|p| if p[0] < 0.004 { 0.0 } else { 400000.0 })
            .collect();
        let unchanged = (p.clone(), t.clone(), v.clone(), h.clone());
        let gate = CancelGate::new_clock_free();
        if cancelled {
            gate.request();
        }
        let refused = with_cx(&gate, |cx| {
            plan_free_surface_from_material(&surface, 0.001, 0.001, 0.1)
                .unwrap()
                .free_surface3_from_tetrahedra_with_boussinesq_enthalpy(
                    MaterialFreeSurfaceDomain3 {
                        gravity_m_s2: [0.0, 0.0, -1.0],
                        ..domain()
                    },
                    [0.0; 3],
                    source(&p, &t, &v, &h),
                    &law,
                    &caloric,
                    bounds(),
                    Default::default(),
                    config(gravity),
                    None,
                    cx,
                    &gate,
                )
                .is_err()
        });
        assert_eq!((p, t, v, h), unchanged);
        refused
    };
    // Every nodal density equals rho0, but the original density claim has an
    // interior 5% departure. The original pinned source-line owner must see it.
    assert!(run(
        vec![(300.0, 1000.0), (400.0, 950.0), (500.0, 1000.0)],
        [0.0, 0.0, -1.0],
        false
    ));
    // Source and geometric transfer pass; native gravity attachment then
    // rejects a mismatched reference rather than adding an extra force owner.
    assert!(run(
        vec![(300.0, 1000.0), (500.0, 980.0)],
        [0.0, 0.0, -2.0],
        false
    ));
    assert!(run(
        vec![(300.0, 1000.0), (500.0, 980.0)],
        [0.0, 0.0, -1.0],
        true
    ));
    assert!(!run(
        vec![(300.0, 1000.0), (500.0, 980.0)],
        [0.0, 0.0, -1.0],
        false
    ));
}
