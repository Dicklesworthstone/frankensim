//! G0 atomic bounded work and G1 independent original-NASA9 endpoint balances.
use super::*;
use crate::evaporation::gas::tests::{
    Phase, binding, catalog, controls as phase_controls, gas, liquid, temperature,
};
use crate::evaporation::gas::{GasEvaporationError, GasLiquidBinding};
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{GasEnergyError, parcel::ClosedGasParcel},
};
use fs_thermal::enthalpy::backward_euler::{
    BoundaryExchangeStepCoefficients, BoundaryStepCoefficients, StepCoefficients,
    TransportStepCoefficients,
};
use fs_thermal::enthalpy::evaporation::ThermalEvaporationError;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};
use fs_thermal::evaporation::LiquidSurfaceState;
use fs_thermochem::SpeciesId;
use std::cell::Cell;
use std::collections::BTreeMap;

const MOLAR_MASS: f64 = 0.018015;
const R: f64 = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / MOLAR_MASS;
const CARRIER_R: f64 = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
const LATENT: f64 = 200_000.0;
const LIQUID_MASS: f64 = 0.01;

// Adapt the existing gas-test liquid inversion to the native thermal owner.
// The gas remains the original variable-capacity NASA9 mixture, and the phase
// source and liquid caloric inversion are reused unchanged from gas/tests.rs.
struct OriginalLiquid;
impl ThermalLaw for OriginalLiquid {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = temperature(h).map_err(|_| DepositionError::MaterialDomain)?;
        if (200.0..=900.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / (4.2 * R)
    }
}
static ORIGINAL_LIQUID: OriginalLiquid = OriginalLiquid;

fn vapor() -> SpeciesId {
    SpeciesId::new("B").unwrap()
}
fn carrier() -> SpeciesId {
    SpeciesId::new("A").unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}
fn controls() -> ImplicitGasControls {
    ImplicitGasControls::new(phase_controls())
}
fn surface<'a>(
    id: u64,
    binding: &'a GasLiquidBinding,
    mass_coefficient: f64,
    heat_coefficient: f64,
) -> ImplicitGasSurface<'a> {
    ImplicitGasSurface {
        cell_id: CellId(id),
        law: &Phase,
        binding,
        pure_liquid_enthalpy_range_j_kg: [
            liquid(200.0).specific_enthalpy_j_kg,
            liquid(900.0).specific_enthalpy_j_kg,
        ],
        mass_transfer_coefficient_m_s: mass_coefficient,
        heat_transfer_coefficient_w_m2_k: heat_coefficient,
    }
}
fn cell(id: u64, t: f64, source_w: f64, area_m2: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id + 10,
        body_id: id,
        material_id: MaterialId(1),
        mass_kg: LIQUID_MASS,
        specific_h_j_kg: liquid(t).specific_enthalpy_j_kg,
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64, 0.01, 0.0],
            volume_m3: 1e-5,
            surface_area_m2: area_m2,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w,
        },
    }
}
fn thermal(
    cells: Vec<DepositedCell>,
    contact: Option<(f64, f64)>,
) -> DepositionSimulation<'static> {
    let mut events: Vec<_> = cells
        .into_iter()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    if let Some((area_m2, conductance_w_k)) = contact {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 9,
                a: CellId(1),
                b: CellId(2),
                area_m2,
                conductance_w_k,
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), &ORIGINAL_LIQUID)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 8,
            max_substeps: 32,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-6,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn coefficients<T>(payload: T) -> BoundaryExchangeStepCoefficients<T> {
    TransportStepCoefficients::from(BoundaryStepCoefficients::from(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload,
    }))
    .into()
}
fn empty_heat(
    _: &DepositionSimulation<'_>,
    _: &BTreeMap<CellId, LiquidSurfaceState>,
    _: &ClosedGasParcel,
) -> Result<BoundaryExchangeStepCoefficients<()>, ImplicitGasError> {
    Ok(coefficients(()))
}
fn unchanged(
    thermal: &DepositionSimulation<'_>,
    original_thermal: &DepositionSimulation<'_>,
    gas: &ClosedGasParcel,
    original_gas: &ClosedGasParcel,
) {
    assert_eq!(thermal.time_s(), original_thermal.time_s());
    assert_eq!(thermal.cells(), original_thermal.cells());
    assert_eq!(thermal.links(), original_thermal.links());
    assert_eq!(thermal.work_budget(), original_thermal.work_budget());
    assert_eq!(gas, original_gas);
}

// Independent integrals of the existing fixture's NASA9 coefficients. These
// supply reference arithmetic only; production evaluation uses the real source.
fn vapor_h(t: f64) -> f64 {
    R * (4.0 * t + 0.0005 * t * t)
}
fn vapor_u(t: f64) -> f64 {
    R * (3.0 * t + 0.0005 * t * t)
}
fn delta_h(a: f64, b: f64) -> f64 {
    R * (b - a) * (4.0 + 0.0005 * (a + b))
}
fn saturation(t: f64) -> f64 {
    1500.0 * fs_math::det::exp(LATENT / R * (1.0 / 300.0 - 1.0 / t))
}
fn unit_film_rate(surface_t: f64, gas_t: f64, pressure: f64, vapor_pressure: f64) -> f64 {
    MOLAR_MASS * pressure / (fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * gas_t)
        * ((pressure - vapor_pressure) / (pressure - saturation(surface_t))).ln()
}
fn gas_capacity_polynomial(gas: &ClosedGasParcel) -> (f64, f64) {
    let carrier_mass = gas.species_mass_kg(&carrier()).unwrap();
    let vapor_mass = gas.species_mass_kg(&vapor()).unwrap();
    (
        carrier_mass * CARRIER_R * 2.5 + vapor_mass * R * 3.0,
        vapor_mass * R * 0.001,
    )
}
fn quadratic_increment(linear: f64, quadratic: f64, rhs: f64) -> f64 {
    2.0 * rhs / (linear + (linear * linear + 4.0 * quadratic * rhs).sqrt())
}

struct Opposing {
    thermal: DepositionSimulation<'static>,
    gas: ClosedGasParcel,
    mass_coefficient: [f64; 2],
    gas_temperature: f64,
    transferred_mass: f64,
    supplied_heat: [f64; 2],
    sensible_heat: [f64; 2],
    external_heat: f64,
    gas_energy: f64,
}
impl Opposing {
    const DT: f64 = 20.0;
    const AREA: f64 = 0.001;
    const HTC: [f64; 2] = [20.0, 10.0];
}
fn opposing(catalog: &ThermochemicalGasModel) -> Opposing {
    let gas = gas(catalog, 320.0, 1600.0, 0.001);
    let d = 2.0f64.powi(-20);
    let (a, b) = gas_capacity_polynomial(&gas);
    let k = Opposing::HTC.map(|htc| htc * Opposing::AREA * Opposing::DT);
    // Equal opposite vapor masses retain composition. Gas warming satisfies
    // ΔU = d*(h_v(360)-h_v(Tg)) - Σ K_i*(Tg-T_i). This is a quadratic
    // obtained from the actual NASA9 source coefficients, not a solver result.
    let linear = a + b * 320.0 + d * R * 4.32 + k[0] + k[1];
    let quadratic = 0.5 * (b + d * R * 0.001);
    let rhs = d * delta_h(320.0, 360.0) + k[0] * 40.0 - k[1] * 40.0;
    let increment = quadratic_increment(linear, quadratic, rhs);
    let gas_temperature = 320.0 + increment;
    let carrier_mass = gas.species_mass_kg(&carrier()).unwrap();
    let vapor_mass = gas.species_mass_kg(&vapor()).unwrap();
    let pressure = (carrier_mass * CARRIER_R + vapor_mass * R) * gas_temperature / gas.volume_m3();
    let vapor_pressure = vapor_mass * R * gas_temperature / gas.volume_m3();
    let mass_coefficient = [360.0, 280.0].map(|t| {
        let signed = if t == 360.0 { d } else { -d };
        signed
            / (Opposing::AREA
                * Opposing::DT
                * unit_film_rate(t, gas_temperature, pressure, vapor_pressure))
    });
    let supplied_heat = [
        LIQUID_MASS * delta_h(355.0, 360.0) + d * LATENT,
        LIQUID_MASS * delta_h(285.0, 280.0) - d * (LATENT + delta_h(280.0, gas_temperature)),
    ];
    let sensible_heat = [
        k[0] * (gas_temperature - 360.0),
        k[1] * (gas_temperature - 280.0),
    ];
    let contact_w = 2.0 * (280.0 - 360.0);
    let source_w = [
        (supplied_heat[0] - sensible_heat[0]) / Opposing::DT - contact_w,
        (supplied_heat[1] - sensible_heat[1]) / Opposing::DT + contact_w,
    ];
    let thermal = thermal(
        vec![
            cell(1, 355.0, source_w[0], 0.0012),
            cell(2, 285.0, source_w[1], 0.0012),
        ],
        Some((0.0002, 2.0)),
    );
    Opposing {
        thermal,
        gas,
        mass_coefficient,
        gas_temperature,
        transferred_mass: d,
        supplied_heat,
        sensible_heat,
        external_heat: Opposing::DT * (source_w[0] + source_w[1]),
        gas_energy: increment * (a + 0.5 * b * (320.0 + gas_temperature)),
    }
}

#[test]
fn opposite_surfaces_share_finite_gas_heat_inventory_and_one_independent_nasa9_endpoint() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let fixture = opposing(&source);
    let original_thermal = fixture.thermal.clone();
    let original_gas = fixture.gas.clone();
    let config = controls();
    assert!(
        fixture.transferred_mass
            < config.phase.maximum_inventory_fraction
                * fixture.gas.species_mass_kg(&vapor()).unwrap()
    );
    let mut outcomes = Vec::new();
    for order in [[1, 2], [2, 1]] {
        let surfaces = order.map(|id| {
            let index = id as usize - 1;
            surface(
                id,
                &binding,
                fixture.mass_coefficient[index],
                Opposing::HTC[index],
            )
        });
        let calls = Cell::new(0usize);
        let candidate = solve_deposition_with_gas(
            &fixture.thermal,
            &surfaces,
            &fixture.gas,
            &source,
            Opposing::DT,
            config,
            |thermal, liquids, gas| {
                calls.set(calls.get() + 1);
                assert_eq!(thermal.cells()[&CellId(1)].mass_kg, LIQUID_MASS);
                assert_eq!(thermal.cells()[&CellId(2)].mass_kg, LIQUID_MASS);
                Ok(coefficients((
                    calls.get(),
                    gas.state().state().temperature,
                    liquids.clone(),
                )))
            },
            || Ok(()),
        )
        .unwrap();
        let report = candidate.report();
        let gas = candidate.gas();
        near(
            gas.state().state().temperature,
            fixture.gas_temperature,
            2e-8,
        );
        near(candidate.payload().1, gas.state().state().temperature, 2e-8);
        assert_eq!(candidate.payload().0, calls.get());
        assert!(calls.get() <= report.residual_evaluations);
        for (index, target) in [360.0, 280.0].into_iter().enumerate() {
            let id = CellId(index as u64 + 1);
            let signed_dm = if index == 0 {
                fixture.transferred_mass
            } else {
                -fixture.transferred_mass
            };
            let receipt = &report.thermal.thermal.surfaces[&id];
            near(candidate.thermal().temperature_k(id).unwrap(), target, 2e-8);
            near(receipt.transferred_mass_kg, signed_dm, 2e-14);
            near(receipt.liquid.mass_kg, LIQUID_MASS - signed_dm, 2e-14);
            near(
                candidate.payload().2[&id].mass_kg,
                receipt.liquid.mass_kg,
                2e-14,
            );
            near(receipt.liquid.area_m2, Opposing::AREA, 1e-17);
            near(receipt.supplied_heat_j, fixture.supplied_heat[index], 1e-7);
            near(
                report.sensible_heat_into_liquid_j[&id],
                fixture.sensible_heat[index],
                1e-7,
            );
            near(receipt.latent_energy_j, signed_dm * LATENT, 1e-8);
            let upstream = if index == 0 {
                target
            } else {
                gas.state().state().temperature
            };
            near(receipt.vapor_energy_j, signed_dm * vapor_h(upstream), 1e-8);
        }
        let hot = &report.thermal.thermal.surfaces[&CellId(1)];
        let cold = &report.thermal.thermal.surfaces[&CellId(2)];
        near(
            hot.transferred_mass_kg + cold.transferred_mass_kg,
            0.0,
            2e-14,
        );
        near(hot.latent_energy_j + cold.latent_energy_j, 0.0, 1e-8);
        assert!(hot.vapor_energy_j + cold.vapor_energy_j > 0.05);
        near(
            report.thermal.vapor_enthalpy_correction_j[&CellId(2)],
            -fixture.transferred_mass * delta_h(280.0, fixture.gas_temperature),
            1e-8,
        );
        near(report.external_heat_j, fixture.external_heat, 1e-7);
        near(
            candidate.gas_exchange().step().stored_energy_change_j,
            fixture.gas_energy,
            1e-7,
        );
        let actual_gas_energy = gas.mass_kg() * gas.state().specific_internal_energy_j_kg()
            - fixture.gas.mass_kg() * fixture.gas.state().specific_internal_energy_j_kg();
        near(actual_gas_energy, fixture.gas_energy, 1e-7);
        near(
            candidate.thermal().energy_j() - fixture.thermal.energy_j() + actual_gas_energy,
            fixture.external_heat,
            1e-7,
        );
        near(
            gas.species_mass_kg(&carrier()).unwrap(),
            fixture.gas.species_mass_kg(&carrier()).unwrap(),
            1e-14,
        );
        near(
            gas.species_mass_kg(&vapor()).unwrap(),
            fixture.gas.species_mass_kg(&vapor()).unwrap(),
            2e-14,
        );
        near(
            candidate.thermal().mass_kg() + gas.mass_kg(),
            fixture.thermal.mass_kg() + fixture.gas.mass_kg(),
            2e-14,
        );
        assert_eq!(gas.volume_m3(), fixture.gas.volume_m3());
        assert_eq!(gas.point().axes().get("local-lot"), Some(&11.0));
        assert!(std::ptr::eq(candidate.catalog(), &source));
        assert_eq!(candidate.thermal().links(), fixture.thermal.links());
        assert_eq!(candidate.thermal().time_s(), Opposing::DT);
        near(report.energy_residual_j, 0.0, 1e-7);
        near(
            report.gas_energy_residual_j,
            0.0,
            config.phase.gas.energy_tolerance_j,
        );
        near(
            report.gas_mass_residual_kg,
            0.0,
            config.phase.gas.mass_tolerance_kg,
        );
        assert!(report.newton_iterations > 0);
        assert!(
            report.gas_iterations > 0
                && report.gas_iterations <= config.phase.gas.gas.maximum_iterations
        );
        outcomes.push((
            candidate.thermal().cells().clone(),
            gas.clone(),
            report.clone(),
        ));
    }
    assert_eq!(outcomes[0], outcomes[1]);
    unchanged(
        &fixture.thermal,
        &original_thermal,
        &fixture.gas,
        &original_gas,
    );
}

#[test]
fn condensation_carries_endpoint_gas_enthalpy_and_cancels_large_formation_references() {
    let mut physical = Vec::new();
    for reference in [0.0, 1e12] {
        let source = catalog(reference);
        let binding = binding(reference);
        let initial = gas(&source, 340.0, 5000.0, 0.001);
        let d = 2.0f64.powi(-22);
        let (a, b) = gas_capacity_polynomial(&initial);
        // For condensation ΔU_g = -d*h_v(Tg). Cancelling the lost vapor's
        // original internal energy leaves ΔU_at_initial_masses = -d*R*Tg.
        let increment = quadratic_increment(a + b * 340.0 + d * R, 0.5 * b, -d * R * 340.0);
        let expected_gas_t = 340.0 + increment;
        let carrier_mass = initial.species_mass_kg(&carrier()).unwrap();
        let vapor_mass = initial.species_mass_kg(&vapor()).unwrap();
        let pressure = (carrier_mass * CARRIER_R + (vapor_mass - d) * R) * expected_gas_t
            / initial.volume_m3();
        let vapor_pressure = (vapor_mass - d) * R * expected_gas_t / initial.volume_m3();
        let coefficient =
            -d / (0.001 * unit_film_rate(300.0, expected_gas_t, pressure, vapor_pressure));
        let supplied = -d * (LATENT + delta_h(300.0, expected_gas_t));
        let old = thermal(vec![cell(1, 300.0, supplied, 0.001)], None);
        let snapshot = old.clone();
        let original_gas = initial.clone();
        let mut config = controls();
        config.phase.latent_compatibility_tolerance_j_kg = 0.01;
        let candidate = solve_deposition_with_gas(
            &old,
            &[surface(1, &binding, coefficient, 0.0)],
            &initial,
            &source,
            1.0,
            config,
            empty_heat,
            || Ok(()),
        )
        .unwrap();
        let report = candidate.report();
        let receipt = &report.thermal.thermal.surfaces[&CellId(1)];
        let tg = candidate.gas().state().state().temperature;
        let tl = candidate.thermal().temperature_k(CellId(1)).unwrap();
        near(tg, expected_gas_t, 2e-8);
        near(tl, 300.0, 2e-8);
        near(receipt.transferred_mass_kg, -d, 2e-14);
        near(receipt.liquid.mass_kg, LIQUID_MASS + d, 2e-14);
        assert!(
            candidate
                .gas()
                .species_partial_pressure_pa(&vapor())
                .unwrap()
                < 5000.0
        );
        near(
            candidate.gas().species_mass_kg(&vapor()).unwrap(),
            vapor_mass - d,
            2e-14,
        );
        near(
            candidate.gas().species_mass_kg(&carrier()).unwrap(),
            carrier_mass,
            1e-14,
        );
        near(
            receipt.vapor_energy_j,
            receipt.transferred_mass_kg * vapor_h(tg),
            1e-8,
        );
        near(
            receipt.latent_energy_j,
            receipt.transferred_mass_kg * LATENT,
            1e-10,
        );
        let correction = report.thermal.vapor_enthalpy_correction_j[&CellId(1)];
        near(
            correction,
            receipt.transferred_mass_kg * delta_h(tl, tg),
            1e-8,
        );
        assert!(
            correction < -0.01,
            "upstream sensible vapor energy must remain explicit"
        );
        near(
            candidate.gas_exchange().step().energy_into_gas_j,
            receipt.transferred_mass_kg * (vapor_h(tg) + reference),
            1e-6,
        );
        // Evaluate the independent physical balance after cancelling formation
        // references analytically, without subtracting huge absolute reports.
        let final_vapor_mass = candidate.gas().species_mass_kg(&vapor()).unwrap();
        let gas_change = carrier_mass * CARRIER_R * 2.5 * (tg - 340.0)
            + final_vapor_mass * R * (tg - 340.0) * (3.0 + 0.0005 * (tg + 340.0))
            + (final_vapor_mass - vapor_mass) * vapor_u(340.0);
        let liquid_change = receipt.liquid.mass_kg * delta_h(300.0, tl)
            - receipt.transferred_mass_kg * (vapor_h(300.0) - LATENT);
        near(gas_change + liquid_change, supplied, 1e-8);
        near(report.external_heat_j, supplied, 1e-10);
        near(report.energy_residual_j, 0.0, 1e-7);
        near(
            report.gas_energy_residual_j,
            0.0,
            config.phase.gas.energy_tolerance_j,
        );
        near(
            report.gas_mass_residual_kg,
            0.0,
            config.phase.gas.mass_tolerance_kg,
        );
        assert_eq!(report.sensible_heat_into_liquid_j[&CellId(1)], 0.0);
        physical.push((tg, tl, receipt.liquid.mass_kg, correction));
        unchanged(&old, &snapshot, &initial, &original_gas);
    }
    near(physical[0].0, physical[1].0, 2e-8);
    near(physical[0].1, physical[1].1, 2e-8);
    near(physical[0].2, physical[1].2, 2e-14);
    near(physical[0].3, physical[1].3, 1e-8);
}

#[test]
fn gross_vapor_donor_and_original_caloric_identity_cannot_be_hidden_by_net_transfer() {
    let source = catalog(0.0);
    let fixture = opposing(&source);
    let original_thermal = fixture.thermal.clone();
    let original_gas = fixture.gas.clone();
    let good = binding(0.0);
    let surfaces = [
        surface(1, &good, fixture.mass_coefficient[0], Opposing::HTC[0]),
        surface(2, &good, fixture.mass_coefficient[1], Opposing::HTC[1]),
    ];
    let mut restricted = controls();
    restricted.phase.maximum_inventory_fraction = 0.05;
    assert!(fixture.transferred_mass > 0.05 * fixture.gas.species_mass_kg(&vapor()).unwrap());
    let result = solve_deposition_with_gas(
        &fixture.thermal,
        &surfaces,
        &fixture.gas,
        &source,
        Opposing::DT,
        restricted,
        empty_heat,
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(ImplicitGasError::Phase(
            GasEvaporationError::StepSizeLimit { .. }
        ))
    ));
    unchanged(
        &fixture.thermal,
        &original_thermal,
        &fixture.gas,
        &original_gas,
    );

    let bad_datum = binding(100_000.0);
    let result = solve_deposition_with_gas(
        &fixture.thermal,
        &[surface(1, &bad_datum, fixture.mass_coefficient[0], 0.0)],
        &fixture.gas,
        &source,
        1.0,
        controls(),
        empty_heat,
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(ImplicitGasError::Phase(
            GasEvaporationError::CaloricMismatch { .. }
        ))
    ));
    let wrong_species = GasLiquidBinding {
        species: carrier(),
        liquid_enthalpy_offset_j_kg: 0.0,
    };
    let result = solve_deposition_with_gas(
        &fixture.thermal,
        &[surface(1, &wrong_species, fixture.mass_coefficient[0], 0.0)],
        &fixture.gas,
        &source,
        1.0,
        controls(),
        empty_heat,
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(ImplicitGasError::Phase(GasEvaporationError::Invalid(_)))
    ));
    unchanged(
        &fixture.thermal,
        &original_thermal,
        &fixture.gas,
        &original_gas,
    );
}

#[test]
fn nonlinear_residual_and_final_native_gas_work_are_separately_bounded() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let fixture = opposing(&source);
    let original_thermal = fixture.thermal.clone();
    let original_gas = fixture.gas.clone();
    let surfaces = [
        surface(1, &binding, fixture.mass_coefficient[0], Opposing::HTC[0]),
        surface(2, &binding, fixture.mass_coefficient[1], Opposing::HTC[1]),
    ];
    let mut limited = controls();
    limited.thermal.maximum_evaluations = 1;
    let calls = Cell::new(0usize);
    let result = solve_deposition_with_gas(
        &fixture.thermal,
        &surfaces,
        &fixture.gas,
        &source,
        Opposing::DT,
        limited,
        |_, _, _| {
            calls.set(calls.get() + 1);
            Ok(coefficients(()))
        },
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
            DepositionError::Budget("implicit finite-gas residual evaluations")
        )))
    ));
    assert_eq!(calls.get(), 1);
    unchanged(
        &fixture.thermal,
        &original_thermal,
        &fixture.gas,
        &original_gas,
    );

    limited = controls();
    limited.thermal.maximum_iterations = 1;
    assert!(matches!(
        solve_deposition_with_gas(
            &fixture.thermal,
            &surfaces,
            &fixture.gas,
            &source,
            Opposing::DT,
            limited,
            empty_heat,
            || Ok(()),
        ),
        Err(ImplicitGasError::NonConverged(_))
    ));
    unchanged(
        &fixture.thermal,
        &original_thermal,
        &fixture.gas,
        &original_gas,
    );

    limited = controls();
    limited.phase.gas.gas.maximum_iterations = 1;
    assert!(matches!(
        solve_deposition_with_gas(
            &fixture.thermal,
            &surfaces,
            &fixture.gas,
            &source,
            Opposing::DT,
            limited,
            empty_heat,
            || Ok(()),
        ),
        Err(ImplicitGasError::Phase(GasEvaporationError::Gas(
            GasEnergyError::NotConverged { iterations: 1, .. }
        )))
    ));
    unchanged(
        &fixture.thermal,
        &original_thermal,
        &fixture.gas,
        &original_gas,
    );
}

#[test]
fn final_shared_sources_are_fresh_funded_and_atomic_through_the_last_checkpoint() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let fixture = opposing(&source);
    let original_thermal = fixture.thermal.clone();
    let original_gas = fixture.gas.clone();
    let surfaces = [
        surface(1, &binding, fixture.mass_coefficient[0], Opposing::HTC[0]),
        surface(2, &binding, fixture.mass_coefficient[1], Opposing::HTC[1]),
    ];
    let calls = Cell::new(0usize);
    let mut polls = 0usize;
    let accepted = solve_deposition_with_gas(
        &fixture.thermal,
        &surfaces,
        &fixture.gas,
        &source,
        Opposing::DT,
        controls(),
        |_, _, gas| {
            calls.set(calls.get() + 1);
            Ok(coefficients((calls.get(), gas.clone())))
        },
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(accepted.payload().0, calls.get());
    assert_eq!(&accepted.payload().1, accepted.gas());
    let evaluations = accepted.report().residual_evaluations;
    assert!(evaluations > 2);
    assert!(calls.get() >= 2 && calls.get() <= evaluations);

    let final_call = calls.get();
    calls.set(0);
    let result = solve_deposition_with_gas(
        &fixture.thermal,
        &surfaces,
        &fixture.gas,
        &source,
        Opposing::DT,
        controls(),
        |_, _, _| {
            calls.set(calls.get() + 1);
            if calls.get() == final_call {
                Err(DepositionError::InvalidInput("final original gas/heat source refused").into())
            } else {
                Ok(coefficients(()))
            }
        },
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput("final original gas/heat source refused")
        )))
    ));
    assert_eq!(calls.get(), final_call);
    unchanged(
        &fixture.thermal,
        &original_thermal,
        &fixture.gas,
        &original_gas,
    );

    let mut limited = controls();
    limited.thermal.maximum_evaluations = evaluations - 1;
    assert!(matches!(
        solve_deposition_with_gas(
            &fixture.thermal,
            &surfaces,
            &fixture.gas,
            &source,
            Opposing::DT,
            limited,
            empty_heat,
            || Ok(()),
        ),
        Err(ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
            DepositionError::Budget("implicit finite-gas residual evaluations")
        )))
    ));
    unchanged(
        &fixture.thermal,
        &original_thermal,
        &fixture.gas,
        &original_gas,
    );

    for cancel_at in [1, polls / 2, polls] {
        let mut count = 0usize;
        let result = solve_deposition_with_gas(
            &fixture.thermal,
            &surfaces,
            &fixture.gas,
            &source,
            Opposing::DT,
            controls(),
            empty_heat,
            || {
                count += 1;
                if count == cancel_at {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        );
        assert!(matches!(
            result,
            Err(ImplicitGasError::Thermal(ThermalEvaporationError::Thermal(
                DepositionError::Cancelled
            )))
        ));
        assert_eq!(count, cancel_at);
        unchanged(
            &fixture.thermal,
            &original_thermal,
            &fixture.gas,
            &original_gas,
        );
    }
}
