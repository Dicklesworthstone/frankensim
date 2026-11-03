//! Native phase, conduction, finite walls and electrical heat share one clock.

use super::*;
use fs_couple::liquid_gas::evaporation::{
    conduction::{ConductingLiquidGasCandidate, advance_liquid_surfaces_with_gas_and_conduction},
    multiple::{DistributedLiquidGasSolve, advance_liquid_surfaces_with_gas},
};
use fs_evidence::ValidityDomain;
use fs_lbm::conjugate::{ElectricSetup, Electrode, ElectrodeDrive};
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceCellHeat,
    conduction::{
        FreeSurfaceBodyContact, FreeSurfaceBodyGasFace, FreeSurfaceBodyRadiationFace,
        FreeSurfaceConductionControls, FreeSurfaceConductionFace, FreeSurfaceElectricalInput,
        FreeSurfaceGasFace, FreeSurfaceThermalBody,
    },
};
use fs_matdb::*;
use fs_material::evaporation::gas::chamber::GasLiquidChamberSolve;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::{
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, ThermalMaterialPoint,
};
use fs_thermal::enthalpy::CellGeometry;

pub(super) fn joint_controls() -> DistributedLiquidGasSolve {
    let single = controls();
    DistributedLiquidGasSolve {
        chamber: GasLiquidChamberSolve {
            surface: single.phase,
            maximum_surfaces: 2,
        },
        native: single.native,
        mass_tolerance_kg: single.mass_tolerance_kg,
        energy_tolerance_j: single.energy_tolerance_j,
    }
}

pub(super) fn thermal_controls() -> FreeSurfaceConductionControls {
    FreeSurfaceConductionControls {
        maximum_faces: 8,
        maximum_substeps: 1000,
        maximum_dt_s: DT / 4.0,
        energy_tolerance_j: 1e-9,
    }
}

pub(super) fn film(liquid: [usize; 3]) -> FreeSurfaceGasFace {
    FreeSurfaceGasFace {
        liquid,
        gas: [liquid[0] + 1, liquid[1], liquid[2]],
        area_m2: 1e-6,
        distances_m: [0.0005; 2],
    }
}

pub(super) fn with_wall<R>(
    run: impl FnOnce(FreeSurfaceThermalBody<'_>, FreeSurfaceBodyContact) -> R,
) -> R {
    let mut claims = ClaimSet::new();
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    for (key, value, dims) in [
        (conductivity.clone(), 0.5, THERMAL_CONDUCTIVITY_DIMS),
        (
            PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            0.6,
            EMISSIVITY_DIMS,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 340.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured finite phase/conduction wall".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured solid".into(),
            phase: "solid".into(),
            process: "authored".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        [280.0, 340.0]
            .map(|temperature_k| EnthalpyPhaseKnot {
                temperature_k,
                specific_enthalpy_j_kg: 500.0 * (temperature_k - 280.0),
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 2000.0,
            })
            .to_vec(),
    )
    .unwrap();
    let point = QueryPoint::new().with("T", 310.0).unwrap();
    let body = FreeSurfaceThermalBody::try_new(
        7,
        ThermalMaterialPoint {
            material: &card,
            context: &point,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        &caloric,
        CellGeometry {
            start_m: [-0.001, 0.0035, 0.0035],
            end_m: [0.0, 0.0035, 0.0035],
            volume_m3: 1e-9,
            surface_area_m2: 6e-6,
        },
        15000.0,
        &conductivity,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
    .with_gas_faces(vec![FreeSurfaceBodyGasFace {
        wall: [-1, 0, 3],
        gas: [0, 0, 3],
        area_m2: 1e-6,
        distances_m: [0.0005; 2],
    }])
    .with_radiation_faces(
        vec![FreeSurfaceBodyRadiationFace {
            wall: [-1, 0, 3],
            gas: [0, 0, 3],
            area_m2: 1e-6,
            surroundings_k: 450.0,
        }],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    run(
        body,
        FreeSurfaceBodyContact {
            body_id: 7,
            liquid: [0, 3, 3],
            wall: [-1, 3, 3],
            area_m2: 1e-6,
            distances_m: [0.0005; 2],
        },
    )
}

fn inventory(
    liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    gas: &ClosedGasParcel,
    binding: &GasLiquidBinding,
    bodies: &[FreeSurfaceThermalBody<'_>],
    next: &ConductingLiquidGasCandidate<'_, '_, '_, '_>,
    external_j: f64,
) {
    let dm =
        next.liquid().flow().ledger_mass_kg().unwrap() - liquid.flow().ledger_mass_kg().unwrap();
    let dh = next.liquid().total_enthalpy_j().unwrap() - liquid.total_enthalpy_j().unwrap();
    let du = next.gas().mass_kg() * next.gas().state().specific_internal_energy_j_kg()
        - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
    assert_eq!(next.bodies().len(), bodies.len());
    let dw = next
        .bodies()
        .iter()
        .zip(bodies)
        .map(|(next, old)| next.enthalpy_j() - old.enthalpy_j())
        .sum::<f64>();
    near(
        dm + next.gas().mass_kg() - gas.mass_kg(),
        0.0,
        controls().mass_tolerance_kg,
    );
    near(
        dm + next.gas().species_mass_kg(&water()).unwrap() - gas.species_mass_kg(&water()).unwrap(),
        0.0,
        controls().mass_tolerance_kg,
    );
    near(
        next.gas().species_mass_kg(&carrier()).unwrap(),
        gas.species_mass_kg(&carrier()).unwrap(),
        1e-22,
    );
    near(
        dh + binding.liquid_enthalpy_offset_j_kg * dm + du + dw,
        external_j,
        controls().energy_tolerance_j,
    );
    near(
        dh,
        next.report().phase.native_liquid_enthalpy_change_j,
        1e-12,
    );
    assert_eq!(next.gas().volume_m3(), gas.volume_m3());
    assert_eq!(next.gas().point().axes().get("specimen-lot"), Some(&27.0));
    assert!(next.report().absolute_energy_defect_j <= controls().energy_tolerance_j);
    assert!(next.report().phase.absolute_mass_defect_kg <= controls().mass_tolerance_kg);
}

#[test]
fn g1_opposing_phase_streams_share_original_conduction_wall_and_finite_gas_storage() {
    let cold = [3, 2, 5];
    let warm = [3, 4, 2];
    with_wall(|body, contact| {
        with_liquid_interfaces(
            0.0,
            &[(cold, 0.4, 295.0), (warm, 0.4, 305.0)],
            |liquid, gate| {
                let source = source(-13.4e6);
                let gas = gas(&source, 300.0, 3500.0);
                let original_gas = gas.clone();
                let original_bits = liquid.flow().lattice().grid().state_bits();
                let initial_wall = body.enthalpy_j();
                let binding = binding(&source, 0.0);
                let surfaces = [cold, warm].map(|liquid| LiquidGasEvaporationSurface {
                    liquid,
                    gas: film(liquid).gas,
                    ..surface(&Phase, &binding)
                });
                let bodies = [body];
                let next = advance_liquid_surfaces_with_gas_and_conduction(
                    liquid,
                    &gas,
                    &source,
                    &surfaces,
                    &bodies,
                    &[FreeSurfaceConductionFace {
                        cells: [cold, [2, 2, 5]],
                        area_m2: 1e-6,
                        distances_m: [0.0005; 2],
                    }],
                    &[FreeSurfaceCellHeat {
                        coordinates: [1, 3, 3],
                        energy_j: 2e-6,
                    }],
                    &[contact],
                    &[],
                    &[film(cold), film(warm)],
                    joint_controls(),
                    thermal_controls(),
                    controls().phase.gas.gas,
                    None,
                    true,
                    gate,
                )
                .unwrap();
                let report = next.report();
                let thermal = &report.thermal.conduction;
                assert!(report.phase.surfaces[0].phase.transferred_mass_kg < 0.0);
                assert!(report.phase.surfaces[1].phase.transferred_mass_kg > 0.0);
                assert_eq!(report.phase.flow, thermal.flow);
                assert_eq!(report.phase.mass_exchange.physical_time_s, 0.0);
                assert_eq!(next.liquid().accepted_steps(), 1);
                assert_eq!(next.liquid().physical_time_s(), DT);
                assert!(thermal.thermal.substeps >= 4);
                assert!(!thermal.last_applied.is_empty());
                assert_eq!(thermal.environment.as_ref().unwrap().last_applied.len(), 2);
                assert!(thermal.bodies.as_ref().unwrap().heat_to_liquid_j > 0.0);
                assert!(thermal.body_environment.as_ref().unwrap().heat_to_solids_j < 0.0);
                let radiant = thermal.body_radiation.as_ref().unwrap().heat_to_solids_j;
                assert!(radiant > 0.0);
                // Initial populations are equilibrium: the original stress
                // estimator may return zero despite the spatial velocity gradient.
                // Check that its actual molecular receipts enter the flow once.
                assert!(!report.thermal.applied.is_empty());
                assert!(
                    report
                        .thermal
                        .applied
                        .iter()
                        .all(|sample| sample.heat_j.is_finite() && sample.heat_j >= 0.0)
                );
                let molecular_j = report
                    .thermal
                    .applied
                    .iter()
                    .map(|sample| sample.heat_j)
                    .sum::<f64>();
                near(thermal.flow.viscous_heat_j, molecular_j, 1e-24);
                assert_eq!(
                    report.thermal_gas_temperature_k[0],
                    report.phase.exchange.gas_temperature_k[1]
                );
                assert_eq!(
                    report.thermal_gas_temperature_k[1],
                    next.gas().state().state().temperature
                );
                inventory(
                    liquid,
                    &gas,
                    &binding,
                    &bodies,
                    &next,
                    2e-6 + thermal.flow.viscous_heat_j + radiant,
                );
                assert_eq!(gas, original_gas);
                assert_eq!(liquid.flow().lattice().grid().state_bits(), original_bits);
                assert_eq!(bodies[0].enthalpy_j(), initial_wall);
                assert!(std::ptr::eq(next.catalog(), &source));
            },
        );
    });
}

#[test]
fn g1_pressure_uses_the_same_post_thermal_gas_and_receipt_begins_after_phase() {
    with_liquid_pressure(0.0, 0.4, Some(100000.0), |liquid, gate| {
        let source = source(-13.4e6);
        let initial = gas(&source, 320.0, 6000.0);
        let gas = ClosedGasParcel::new(initial.model().clone(), initial.point().clone(), 3.686e-7)
            .unwrap();
        let binding = binding(&source, 0.0);
        let surfaces = [surface(&Phase, &binding)];
        let phase = advance_liquid_surfaces_with_gas(
            liquid,
            &gas,
            &source,
            &surfaces,
            joint_controls(),
            gate,
        )
        .unwrap();
        let advance = |pressure| {
            advance_liquid_surfaces_with_gas_and_conduction(
                liquid,
                &gas,
                &source,
                &surfaces,
                &[],
                &[],
                &[],
                &[],
                &[],
                &[film(INTERFACE)],
                joint_controls(),
                thermal_controls(),
                controls().phase.gas.gas,
                pressure,
                false,
                gate,
            )
            .unwrap()
        };
        let passive = advance(None);
        let loaded = advance(Some(0.05));
        let thermal = &loaded.report().thermal.conduction;
        let receipt = thermal.finite_gas.as_ref().unwrap().pressure.unwrap();
        assert_eq!(loaded.report().phase.phase, phase.report().phase);
        assert_eq!(loaded.gas(), passive.gas());
        assert_eq!(loaded.gas().model(), phase.gas().model());
        assert_eq!(loaded.gas().mass_kg(), phase.gas().mass_kg());
        assert_eq!(
            receipt.initial_pressure_pa,
            phase.gas().state().state().pressure
        );
        assert_eq!(
            receipt.applied_pressure_pa,
            loaded.gas().state().state().pressure
        );
        assert_eq!(loaded.report().phase.pressure, Some(receipt));
        assert_eq!(receipt.volume_m3, gas.volume_m3());
        assert_eq!(receipt.maximum_density_excursion, 0.05);
        assert!(
            loaded.report().thermal_gas_temperature_k[1]
                < loaded.report().thermal_gas_temperature_k[0]
        );
        assert_ne!(
            loaded.liquid().flow().lattice().grid().state_bits(),
            passive.liquid().flow().lattice().grid().state_bits()
        );
        assert_eq!(loaded.liquid().accepted_steps(), 1);
        assert_eq!(loaded.liquid().physical_time_s(), DT);
        inventory(liquid, &gas, &binding, &[], &loaded, 0.0);
    });
}

#[test]
fn g1_existing_electrical_core_supplies_phase_and_thermal_transaction_once() {
    with_electrical_liquid(0.0, 0.4, |liquid, gate| {
        let liquid = liquid
            .clone()
            .with_electrical_core(
                FreeSurfaceElectricalInput {
                    origin: [0, 3, 3],
                    dimensions: [3, 1, 1],
                    setup: ElectricSetup {
                        conductors: vec![],
                        electrodes: vec![
                            Electrode {
                                lo: [0; 3],
                                hi: [1; 3],
                                drive: ElectrodeDrive::Voltage(0.0),
                            },
                            Electrode {
                                lo: [2, 0, 0],
                                hi: [3, 1, 1],
                                drive: ElectrodeDrive::Voltage(0.1),
                            },
                        ],
                    },
                    contacts: vec![],
                    maximum_cells: 3,
                    maximum_solves: 32,
                    power_tolerance_w: 1e-9,
                },
                gate,
            )
            .unwrap();
        let source = source(-13.4e6);
        let gas = gas(&source, 320.0, 6000.0);
        let binding = binding(&source, 0.0);
        let surfaces = [surface(&Phase, &binding)];
        // The old direct-flow entry must retain its electrical-owner refusal.
        assert!(
            advance_liquid_surfaces_with_gas(
                &liquid,
                &gas,
                &source,
                &surfaces,
                joint_controls(),
                gate,
            )
            .is_err()
        );
        let next = advance_liquid_surfaces_with_gas_and_conduction(
            &liquid,
            &gas,
            &source,
            &surfaces,
            &[],
            &[],
            &[],
            &[],
            &[],
            &[film(INTERFACE)],
            joint_controls(),
            thermal_controls(),
            controls().phase.gas.gas,
            None,
            false,
            gate,
        )
        .unwrap();
        let electric = next
            .report()
            .thermal
            .conduction
            .electrical
            .as_ref()
            .unwrap();
        // One 1 mm bulk cell, sigma=1000 S/m, area=1 mm²: R=1 ohm.
        let supplied = 0.1_f64.powi(2) * DT;
        near(electric.joule_heat_j, supplied, 1e-14);
        near(electric.electrode_energy_j, supplied, 1e-14);
        assert_eq!(electric.substeps, 4);
        assert_eq!(electric.current_solves, 8);
        assert!(next.report().phase.phase.transferred_mass_kg < 0.0);
        assert_eq!(next.liquid().accepted_steps(), 1);
        assert_eq!(liquid.accepted_steps(), 0);
        inventory(&liquid, &gas, &binding, &[], &next, supplied);
    });
}

#[test]
fn g3_late_thermal_failure_discards_phase_gas_native_liquid_and_finite_solids() {
    with_wall(|body, contact| {
        with_liquid(0.0, 0.4, |liquid, gate| {
            let source = source(-13.4e6);
            let gas = gas(&source, 320.0, 6000.0);
            let original_gas = gas.clone();
            let original_bits = liquid.flow().lattice().grid().state_bits();
            let original_h = liquid.total_enthalpy_j().unwrap();
            let binding = binding(&source, 0.0);
            let initial_wall = body.enthalpy_j();
            let bodies = [body];
            let mut exhausted = thermal_controls();
            exhausted.maximum_substeps = 1;
            for (limits, heat) in [
                (exhausted, vec![]),
                (
                    thermal_controls(),
                    vec![FreeSurfaceCellHeat {
                        coordinates: [1, 3, 3],
                        energy_j: 1.0,
                    }],
                ),
            ] {
                let result = advance_liquid_surfaces_with_gas_and_conduction(
                    liquid,
                    &gas,
                    &source,
                    &[surface(&Phase, &binding)],
                    &bodies,
                    &[],
                    &heat,
                    &[contact],
                    &[],
                    &[film(INTERFACE)],
                    joint_controls(),
                    limits,
                    controls().phase.gas.gas,
                    None,
                    false,
                    gate,
                );
                assert!(matches!(result, Err(EvaporatingLiquidGasError::Liquid(_))));
                assert_eq!(gas, original_gas);
                assert_eq!(liquid.flow().lattice().grid().state_bits(), original_bits);
                assert_eq!(liquid.total_enthalpy_j().unwrap(), original_h);
                assert_eq!(liquid.accepted_steps(), 0);
                assert_eq!(bodies[0].enthalpy_j(), initial_wall);
            }
            // A refused candidate leaves the same owners usable for a valid step.
            let next = advance_liquid_surfaces_with_gas_and_conduction(
                liquid,
                &gas,
                &source,
                &[surface(&Phase, &binding)],
                &bodies,
                &[],
                &[],
                &[contact],
                &[],
                &[film(INTERFACE)],
                joint_controls(),
                thermal_controls(),
                controls().phase.gas.gas,
                None,
                false,
                gate,
            )
            .unwrap();
            assert_eq!(next.liquid().accepted_steps(), 1);
            assert_ne!(next.bodies()[0].enthalpy_j(), initial_wall);
            assert_ne!(next.gas(), &gas);
        });
    });
}
