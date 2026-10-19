//! Data-only material substitutions through actual flow, heat and motion owners.
use super::*;
use crate::cavity_phs::piston::material::MaterialGasPiston;
use crate::cavity_phs::piston::tests::{controls, drive, gas};
use fs_evidence::ValidityDomain;
use fs_ga::Vec3;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve,
    HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferSolve;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_mbd_material::{RegionBinding, VolumeMoments, resolve_body};
use fs_ports::{ModuleWindowAudit, PortedModule};
use fs_qty::{Density, Temperature};
use std::convert::Infallible;

fn flow() -> GasNozzleStepSolve {
    GasNozzleStepSolve { nozzle: controls().gas, maximum_donor_fraction: 0.1,
        transfer: GasMassTransferSolve { donor: controls().gas, receiver: controls().gas,
            mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 } }
}
fn exchange() -> HeatExchangeSolve {
    HeatExchangeSolve { heat_bracket_j: [-1.0, 1.0], source_tolerance_j: 1e-9,
        energy_tolerance_j: 1e-7, maximum_iterations: 128 }
}
fn nozzle() -> GasNozzle { GasNozzle { area_m2: 1e-5, discharge_coefficient: 0.8 } }
fn reservoir(x: f64) -> ClosedGasParcel {
    let mut gas = gas(x, 0.001);
    gas.add_heat(100.0, controls().gas, || false).unwrap(); gas
}
fn point(t: f64) -> QueryPoint { QueryPoint::new().with("T", t).unwrap().with("lot", 7.0).unwrap() }
fn material(k: f64) -> (MaterialCard, PropertyKey, ResolvedHeatCapacityCurve) {
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cp = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let mut claims = ClaimSet::new();
    for (key, lower, upper) in [(cp.clone(), 500.0, 500.0), (density.clone(), 4000.0, 4000.0),
        (conductivity.clone(), k, 2.0*k)] {
        let dims = key.quantity().dims();
        claims.insert_claim(PropertyClaim { key,
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, lower), (800.0, upper)], dims },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 800.0).with("lot", 7.0, 7.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            observations: vec![], provenance: Provenance { source: "synthetic thermal piston specimen".into(),
                license: "CC0-1.0".into(), artifact: None } }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic specimen".into(),
        phase: "solid".into(), process: "authored fixture".into(), revision: 0 }, claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &cp, density: &density, temperature_axis: "T",
        lower: &point(200.0), upper: &point(800.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-5).unwrap();
    (card, conductivity, caloric)
}
fn body(curve: &EquilibriumEnthalpyPhaseCurve) -> UniformEnthalpyStepInput<'_> {
    // A finite hot lining, separate from the translating inertial specimen.
    UniformEnthalpyStepInput { curve, initial: curve.state_at_specific_enthalpy(150000.0).unwrap(),
        mass_kg: 0.001, volume_m3: 0.001/4000.0, surface_area_m2: 0.02,
        internal_heat_j: 0.0, duration_s: drive(0.0).duration_s }
}
fn layers<'a>(card: &'a MaterialCard, key: &'a PropertyKey, context: &'a QueryPoint)
    -> MaterialConductionLayers<'a>
{
    MaterialConductionLayers { material: card, context, temperature_axis: "T", liquid_fraction_axis: None,
        conductivity: key, selection: MaterialPropertySelection::SingleClaimOnly,
        area_m2: 0.02, gas_distance_m: 0.001, body_distance_m: 0.02 }
}
fn piston() -> GasPiston { GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap() }
fn link(g: f64) -> impl FnMut(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
    -> Result<ConductanceSample<EquilibriumPhaseState>, Infallible>
{
    move |_, source| Ok(ConductanceSample { conductance_w_k: g, source })
}

#[test]
fn material_conductivity_and_inlet_composition_reach_the_same_heat_and_motion() {
    let run = |k, inlet_fraction| {
        let (card, key, caloric) = material(k); let context = point(300.0);
        let inertial = resolve_body(&[RegionBinding { key: "piston", material: &card, point: &context,
            geometry: VolumeMoments::cuboid(Vec3::new(0.1, 0.1, 0.025)).unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly }], 1, || false).unwrap();
        let mut piston = MaterialGasPiston::new(gas(0.75, 0.001), inertial.clone(), 0.01, 0.0).unwrap();
        let mut reservoir = reservoir(inlet_fraction); let mut body = body(caloric.curve());
        let mut heat = 0.0;
        for _ in 0..24 {
            let candidate = piston.nozzle_thermal_layers_candidate(&reservoir, nozzle(), drive(0.0), body,
                controls(), flow(), exchange(), layers(&card, &key, &context), || false).unwrap();
            let (next, next_reservoir, thermal, report) = candidate.into_parts();
            assert_eq!(next.body(), &inertial);
            assert_eq!(next.piston().mass_kg(), inertial.inertia().mass());
            assert_eq!(report.direction, Some(NozzleDirection::IntoPiston));
            assert!(report.total_energy_residual_j.abs() <= exchange().energy_tolerance_j);
            assert!(report.total_mass_residual_kg.abs() <= flow().transfer.mass_tolerance_kg);
            assert_eq!(thermal.state.material_card_identity(), card.content_hash());
            let actual = layers(&card, &key, &context).resolve(&report.thermal.piston.midpoint,
                report.thermal.thermal_midpoint).unwrap();
            assert_eq!(thermal.report, actual.source);
            assert_eq!(report.thermal.conductance_w_k, actual.conductance_w_k);
            assert_eq!(thermal.external_heat_j, -report.thermal.exchanged_heat_j);
            heat += report.thermal.exchanged_heat_j;
            body.initial = thermal.state; piston = next; reservoir = next_reservoir;
        }
        assert!(matches!(piston.piston().audit_last_window().unwrap(), ModuleWindowAudit::Unavailable { .. }));
        (piston, heat, body.initial.temperature_k())
    };
    let (low, low_q, low_t) = run(0.2, 0.01);
    let (high, high_q, high_t) = run(2.0, 0.01);
    assert!(high_q > low_q && low_q > 0.0);
    assert!(high_t < low_t);
    assert!(high.piston().velocity_m_s() > low.piston().velocity_m_s());
    assert_ne!(high.gas().state().state().sound_speed, low.gas().state().state().sound_speed);
    let (other_gas, other_q, _) = run(0.2, 0.9);
    assert_ne!(other_gas.gas().state().state().specific_gas_constant, low.gas().state().state().specific_gas_constant);
    assert_ne!(other_gas.piston().velocity_m_s(), low.piston().velocity_m_s());
    assert_ne!(other_q, low_q);
}

#[test]
fn closed_nozzle_and_insulated_wall_reduce_to_the_original_owners() {
    let (card, key, caloric) = material(0.2); let context = point(300.0);
    let old = piston(); let reservoir = reservoir(0.25); let body = body(caloric.curve());
    let closed = GasNozzle { area_m2: 0.0, ..nozzle() };
    let actual = old.nozzle_thermal_layers_candidate(&reservoir, closed, drive(0.02), body,
        controls(), flow(), exchange(), layers(&card, &key, &context), || false).unwrap();
    let expected = old.thermal_candidate(drive(0.02), body, controls(), exchange(),
        |gas, body| layers(&card, &key, &context).resolve(gas, body), || false).unwrap();
    assert_eq!(actual.piston(), expected.piston());
    assert_eq!(actual.piston().ports, expected.piston().ports);
    assert_eq!(actual.thermal(), expected.thermal());
    assert_eq!(&actual.step().thermal, expected.step());
    assert_eq!(actual.reservoir(), &reservoir);
    assert_eq!(actual.step().direction, None);

    let actual = old.nozzle_thermal_layers_candidate(&reservoir, nozzle(), drive(0.02), body,
        controls(), flow(), exchange(), MaterialConductionLayers { area_m2: 0.0,
            ..layers(&card, &key, &context) }, || false).unwrap();
    let mut expected = old.clone(); let mut expected_reservoir = reservoir.clone();
    let report = expected.step_with_nozzle(&mut expected_reservoir, nozzle(), drive(0.02),
        flow(), controls(), exchange().energy_tolerance_j, || false).unwrap();
    assert_eq!(actual.piston(), &expected);
    assert_eq!(actual.piston().ports, expected.ports);
    assert_eq!(actual.reservoir(), &expected_reservoir);
    assert_eq!(actual.step().thermal.piston, report.piston);
    assert_eq!(actual.step().transfer, report.transfer);
    assert_eq!(actual.thermal().state, body.initial);
    assert_eq!(actual.thermal().external_heat_j, 0.0);
}

#[test]
fn latent_wall_storage_and_both_sources_close_the_whole_trajectory() {
    let knot = |h, t, x| EnthalpyPhaseKnot { specific_enthalpy_j_kg: h, temperature_k: t,
        liquid_mass_fraction: x, bulk_density_kg_m3: 4000.0 };
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(fs_blake3::hash_bytes(b"synthetic latent lining"),
        vec![knot(0.0, 200.0, 0.0), knot(100000.0, 350.0, 0.0),
            knot(160000.0, 350.0, 1.0), knot(400000.0, 800.0, 1.0)]).unwrap();
    let mut body = UniformEnthalpyStepInput { initial: curve.state_at_specific_enthalpy(130000.0).unwrap(),
        internal_heat_j: 0.02, ..body(&curve) };
    let mut piston = piston(); let mut reservoir = reservoir(0.25);
    let initial = piston.clone(); let initial_reservoir = reservoir.clone(); let initial_body = body.initial;
    let input = PistonInput { heat_j: 0.03, external_force_n: 2.0, ..drive(0.0) };
    for count in 1..=20 {
        let candidate = piston.nozzle_thermal_candidate(&reservoir, nozzle(), input, body,
            controls(), flow(), exchange(), link(1.0), || false).unwrap();
        let (next, next_reservoir, thermal, report) = candidate.into_parts();
        assert!(report.thermal.exchanged_heat_j < 0.0);
        assert_eq!(thermal.state.temperature_k(), 350.0);
        assert!(thermal.state.liquid_mass_fraction() > body.initial.liquid_mass_fraction());
        let q = report.thermal.exchanged_heat_j;
        let advected = report.transfer.transferred_enthalpy_j;
        assert!((report.reservoir_energy_change_j + advected).abs() <= 1e-8);
        assert!((report.chamber_energy_change_j + report.thermal.piston.gas.work_by_gas_j
            - advected - input.heat_j - q).abs() <= 1e-7);
        body.initial = thermal.state; piston = next; reservoir = next_reservoir;
        // Independent extensive endpoint energies; neither transfer nor wall Q
        // appears as an external input to this enclosing physical boundary.
        let energy = |parcel: &ClosedGasParcel| parcel.mass_kg()*parcel.state().specific_internal_energy_j_kg();
        let change = (energy(piston.gas()) - energy(initial.gas()))
            + (energy(&reservoir) - energy(&initial_reservoir))
            + body.mass_kg*(body.initial.specific_enthalpy_j_kg()-initial_body.specific_enthalpy_j_kg())
            + 0.5*piston.mass_kg()*piston.velocity_m_s().powi(2)
            + input.exterior_pressure_pa*(piston.gas().volume_m3()-initial.gas().volume_m3())
            - input.external_force_n*(piston.position_m()-initial.position_m())
            - f64::from(count)*(input.heat_j+body.internal_heat_j);
        assert!(change.abs() <= f64::from(count)*exchange().energy_tolerance_j);
    }
}

#[test]
fn refusal_after_species_transfer_and_late_cancellation_preserve_all_inputs() {
    let (_, _, caloric) = material(0.2); let body = body(caloric.curve());
    let mut old = piston(); old.step(drive(0.0), controls(), || false).unwrap();
    let reservoir = reservoir(0.01); let before = old.clone(); let before_reservoir = reservoir.clone();
    let initial_r = old.gas().state().state().specific_gas_constant; let mut saw_new_composition = false;
    let failure = old.nozzle_thermal_candidate(&reservoir, nozzle(), drive(0.0), body,
        controls(), flow(), exchange(), |gas, _| {
            saw_new_composition = gas.state().specific_gas_constant != initial_r;
            Err::<ConductanceSample<()>, _>("wall material unavailable")
        }, || false);
    assert!(saw_new_composition);
    assert!(matches!(failure, Err(ThermalPistonError::Interface("wall material unavailable"))));
    let mut polls = 0;
    let complete = old.nozzle_thermal_candidate(&reservoir, nozzle(), drive(0.0), body,
        controls(), flow(), exchange(), link(1.0), || { polls += 1; false }).unwrap();
    let mut called = 0;
    let cancelled = old.nozzle_thermal_candidate(&reservoir, nozzle(), drive(0.0), body,
        controls(), flow(), exchange(), link(1.0), || { called += 1; called == polls });
    assert!(matches!(cancelled, Err(ThermalPistonError::Cancelled)));
    assert_eq!(old, before); assert_eq!(old.ports, before.ports); assert_eq!(reservoir, before_reservoir);
    assert_eq!(body.initial, caloric.curve().state_at_specific_enthalpy(150000.0).unwrap());
    let retry = old.nozzle_thermal_candidate(&reservoir, nozzle(), drive(0.0), body,
        controls(), flow(), exchange(), link(1.0), || false).unwrap();
    assert_eq!(retry, complete);
}

#[test]
fn material_heat_and_nozzle_motion_refine_with_the_declared_first_order_split() {
    let (card, key, caloric) = material(0.2); let context = point(300.0);
    let mut velocities = Vec::new();
    for count in [16, 32, 64] {
        let mut piston = piston(); let mut reservoir = reservoir(0.25);
        let input = PistonInput { duration_s: 0.002/f64::from(count), ..drive(0.0) };
        let mut body = UniformEnthalpyStepInput { duration_s: input.duration_s, ..body(caloric.curve()) };
        let solve = PistonSolve { momentum_tolerance_n_s: 2e-7*input.duration_s, ..controls() };
        for _ in 0..count {
            let candidate = piston.nozzle_thermal_layers_candidate(&reservoir, nozzle(), input, body,
                solve, flow(), exchange(), layers(&card, &key, &context), || false).unwrap();
            let (next, next_reservoir, thermal, report) = candidate.into_parts();
            assert!(report.total_energy_residual_j.abs() <= exchange().energy_tolerance_j);
            piston = next; reservoir = next_reservoir; body.initial = thermal.state;
        }
        velocities.push(piston.velocity_m_s());
    }
    let coarse = (velocities[0]-velocities[1]).abs(); let fine = (velocities[1]-velocities[2]).abs();
    assert!(fine > 0.0 && coarse/fine > 1.6 && coarse/fine < 2.6,
        "thermal/nozzle velocity refinement {velocities:?}; ratio {}", coarse/fine);
}
