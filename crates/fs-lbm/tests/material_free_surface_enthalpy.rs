//! G1/G3/G4: original VOF routes carry extensive heat and source viscosity.
//! Synthetic constitutive data test coupling and conservation, not experiments.

use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_lbm::d3q19::freesurface3::ContactModel3;
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceCellHeat, FreeSurfaceEnthalpyControls, FreeSurfaceEnthalpyError,
};
use fs_lbm::material::free_surface::{
    MaterialFreeSurfaceDomain3, MaterialFreeSurfaceInitialCell3, plan_free_surface_from_material,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::surface::{LiquidCapillarity, resolve_liquid_surface_state};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};

#[cfg(feature = "enthalpy-coupling")]
#[path = "material_free_surface_enthalpy/conduction.rs"]
mod conduction;

#[path = "material_free_surface_enthalpy/buoyancy.rs"]
mod buoyancy;

#[path = "material_free_surface_enthalpy/viscous.rs"]
mod viscous;

#[path = "material_free_surface_enthalpy/rheology.rs"]
mod rheology;

#[path = "material_free_surface_enthalpy/height_graph.rs"]
mod height_graph;

#[path = "material_free_surface_enthalpy/mass_exchange.rs"]
mod mass_exchange;

fn card(hot_viscosity: f64, hot_density: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        ("density", Dims([-3, 1, 0, 0, 0, 0]), [1000.0, hot_density]),
        (
            "specific_heat_capacity",
            Dims([2, 0, -2, -1, 0, 0]),
            [2000.0; 2],
        ),
        ("thermal_conductivity", Dims([1, 1, -3, -1, 0, 0]), [0.5; 2]),
        (
            "dynamic_viscosity",
            Dims([-1, 1, -1, 0, 0, 0]),
            [0.1, hot_viscosity],
        ),
        ("surface_tension", Dims([0, 1, -2, 0, 0, 0]), [0.001; 2]),
    ] {
        let quantity = if name == "thermal_conductivity" {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (500.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 300.0, 500.0)
                    .with("pressure", 100000.0, 100000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic moving-liquid enthalpy fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", t)
        .unwrap()
        .with("pressure", 100000.0)
        .unwrap()
}

fn caloric(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Liquid,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 300.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 600000.0,
                temperature_k: 600.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
        ],
    )
    .unwrap()
}

fn domain() -> MaterialFreeSurfaceDomain3<'static> {
    MaterialFreeSurfaceDomain3 {
        cells: [8; 3],
        gravity_m_s2: [0.0; 3],
        contact: ContactModel3::Neutral,
        solid_tiles: &[],
        max_cells: 512,
    }
}

fn initial(full: bool) -> Vec<MaterialFreeSurfaceInitialCell3> {
    let mut result = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                if full || ((1..=3).contains(&x) && (1..=6).contains(&y) && (1..=6).contains(&z)) {
                    result.push(MaterialFreeSurfaceInitialCell3 {
                        coordinates: [x, y, z],
                        fill_fraction: 1.0,
                        velocity_m_s: [0.02, 0.0, 0.0],
                    });
                }
            }
        }
    }
    result
}

fn with_cx<T>(gate: &CancelGate, f: impl FnOnce(&fs_exec::Cx<'_>) -> T) -> T {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, Cx, ExecMode, StreamKey};
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn near(a: f64, b: f64, rel: f64) {
    assert!(
        (a - b).abs() <= rel * a.abs().max(b.abs()).max(1e-20),
        "{a} != {b}"
    );
}

#[test]
fn g1_heterogeneous_enthalpy_follows_actual_interface_conversion_mass() {
    let card = card(0.02, 1000.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point(400.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        &card,
        &point(400.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = caloric(&card);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let flow = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain(), &initial(false), cx)
            .unwrap();
        let mut thermal = flow
            .with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                |p| 100000.0 + 50000.0 * p[0] as f64,
                &gate,
            )
            .unwrap();
        let before = thermal.total_enthalpy_j().unwrap();
        let mass = thermal.flow().ledger_mass_kg().unwrap();
        let report = thermal.step(&gate).unwrap();
        let conversions = thermal.flow().lattice().conversions();
        assert!(conversions.to_fluid > 0 && conversions.gas_to_interface > 0);
        assert!(report.mass_routes > 0);
        near(report.enthalpy_after_j, before, 2e-11);
        near(report.mass_after_kg, mass, 1e-12);
        assert!(report.balance_defect_j.abs() < 1e-9);
        let advanced = thermal.specific_enthalpy_at([4, 3, 3]).unwrap();
        assert!((150000.0..=250000.0).contains(&advanced));
        let endpoint = thermal.material_at([4, 3, 3]).unwrap();
        let t = chart
            .state_at_specific_enthalpy(advanced)
            .unwrap()
            .temperature_k();
        assert_eq!(endpoint, &law.at_temperature(Temperature::new(t)).unwrap());
    });
}

fn heated_trajectory(hot_viscosity: f64, apply_heat: bool) -> (Vec<u64>, f64, f64) {
    let card = card(hot_viscosity, 1000.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point(400.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        &card,
        &point(400.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = caloric(&card);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let initial = initial(true);
        let flow = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain(), &initial, cx)
            .unwrap();
        let mut thermal = flow
            .with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                |_| 200000.0,
                &gate,
            )
            .unwrap();
        let before = thermal.total_enthalpy_j().unwrap();
        let heat = if apply_heat {
            initial
                .iter()
                .map(|cell| FreeSurfaceCellHeat {
                    coordinates: cell.coordinates,
                    energy_j: 0.1,
                })
                .collect::<Vec<_>>()
        } else {
            Vec::new()
        };
        let report = thermal.step_with_heat(&heat, &gate).unwrap();
        near(
            report.enthalpy_after_j,
            before + 0.1 * heat.len() as f64,
            2e-11,
        );
        for _ in 0..2 {
            thermal.step(&gate).unwrap();
        }
        let h = thermal.specific_enthalpy_at([3, 3, 3]).unwrap();
        near(h, if apply_heat { 300000.0 } else { 200000.0 }, 2e-10);
        let mu = thermal
            .material_at([3, 3, 3])
            .unwrap()
            .dynamic_viscosity_pa_s();
        (
            thermal.flow().lattice().grid().state_bits(),
            thermal.total_enthalpy_j().unwrap(),
            mu,
        )
    })
}

#[test]
fn g3_heat_and_data_only_viscosity_changes_reach_native_collision() {
    let cold = heated_trajectory(0.02, false);
    let hot = heated_trajectory(0.02, true);
    let substituted = heated_trajectory(0.08, true);
    assert_ne!(
        cold.0, hot.0,
        "caloric heating must change the actual momentum trajectory"
    );
    assert_ne!(
        hot.0, substituted.0,
        "a viscosity-curve data change must reach native collision"
    );
    assert!(hot.2 < cold.2 && hot.2 < substituted.2);
    near(hot.1, substituted.1, 2e-11);
}

#[test]
fn g4_heat_domain_budget_and_cancellation_refusals_preserve_every_accepted_owner() {
    let card = card(0.02, 1000.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point(400.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        &card,
        &point(400.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = caloric(&card);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let flow = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain(), &initial(false), cx)
            .unwrap();
        let mut thermal = flow
            .with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls {
                    max_mass_routes: 1,
                    ..FreeSurfaceEnthalpyControls::default()
                },
                |_| 200000.0,
                &gate,
            )
            .unwrap();
        let bits = thermal.flow().lattice().grid().state_bits();
        let energy = thermal.total_enthalpy_j().unwrap();
        let state = thermal.material_at([2, 3, 3]).unwrap().clone();
        // Caloric T=550K is supported by this chart, but not by the pinned source.
        let heat = [FreeSurfaceCellHeat {
            coordinates: [2, 3, 3],
            energy_j: 0.3,
        }];
        assert!(matches!(
            thermal.step_with_heat(&heat, &gate),
            Err(FreeSurfaceEnthalpyError::Material(_))
        ));
        assert!(
            thermal.step(&gate).is_err(),
            "the actual route budget must refuse after candidate work"
        );
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(thermal.step(&cancelled).is_err());
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
        assert_eq!(thermal.material_at([2, 3, 3]).unwrap(), &state);
        assert_eq!(thermal.flow().lattice().conversions().to_fluid, 0);
    });
}

#[test]
fn g0_density_change_and_sourced_capillarity_require_their_missing_work_owners() {
    let card = card(0.02, 900.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point(300.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = caloric(&card);
    let gate = CancelGate::new_clock_free();
    for capillarity in [LiquidCapillarity::Neglected, LiquidCapillarity::Sourced] {
        let source = resolve_liquid_surface_state(
            &card,
            &point(300.0),
            LiquidViscosityBasis::Dynamic,
            capillarity,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        with_cx(&gate, |cx| {
            let flow = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
                .unwrap()
                .free_surface3_from_cells(domain(), &initial(false), cx)
                .unwrap();
            let result = flow.with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                |_| 0.0,
                &gate,
            );
            if capillarity == LiquidCapillarity::Sourced {
                assert!(result.is_err());
            } else {
                let mut thermal = result.unwrap();
                let before = thermal.flow().lattice().grid().state_bits();
                assert!(
                    thermal
                        .step_with_heat(
                            &[FreeSurfaceCellHeat {
                                coordinates: [2, 3, 3],
                                energy_j: 0.1,
                            }],
                            &gate
                        )
                        .is_err()
                );
                assert_eq!(thermal.flow().lattice().grid().state_bits(), before);
                assert_eq!(thermal.total_enthalpy_j().unwrap(), 0.0);
            }
        });
    }
}

#[test]
fn g0_freezing_on_a_latent_plateau_refuses_without_erasing_motion_or_heat() {
    let card = card(0.02, 1000.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point(400.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        &card,
        &point(400.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (0.0, 300.0, 0.0),
            (100000.0, 350.0, 0.0),
            (200000.0, 350.0, 1.0),
            (500000.0, 500.0, 1.0),
        ]
        .into_iter()
        .map(|(h, t, liquid)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: liquid,
            bulk_density_kg_m3: 1000.0,
        })
        .collect(),
    )
    .unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let flow = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain(), &initial(false), cx)
            .unwrap();
        let mut thermal = flow
            .with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                |_| 300000.0,
                &gate,
            )
            .unwrap();
        let bits = thermal.flow().lattice().grid().state_bits();
        let enthalpy = thermal.total_enthalpy_j().unwrap();
        assert!(
            thermal
                .step_with_heat(
                    &[FreeSurfaceCellHeat {
                        coordinates: [2, 3, 3],
                        energy_j: -0.15,
                    }],
                    &gate
                )
                .is_err()
        );
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(thermal.total_enthalpy_j().unwrap(), enthalpy);
        assert_eq!(thermal.specific_enthalpy_at([2, 3, 3]), Some(300000.0));
    });
}
