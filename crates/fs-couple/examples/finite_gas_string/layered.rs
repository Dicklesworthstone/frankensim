//! A bounded executable adapter over the native layered string and thermal owners.
use super::*;
use fs_couple::acoustic_realize::concentric::{
    ConcentricStringThermalCell, ThermalConcentricStringRuntime,
    finite_gas::ConcentricFiniteGasSolve,
};
use fs_material::conductor::ElectricalPropertyBasis;
use fs_couple::acoustic_realize::finite_gas::flow::FiniteGasStringFlowDirection;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::section_temperature::{
    ConcentricSectionTemperatureLaw, thermoelastic::ConcentricThermoelasticLaw,
};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    THERMAL_CONDUCTIVITY_DIMS,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation,
    Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
    backward_euler::StepCoefficients, series_conductance_w_k,
};
use std::collections::BTreeMap;

fn conductances(
    bonded: &input::Bonded, section: &ConcentricElasticSection,
) -> Result<BTreeMap<u64, f64>> {
    bonded.sleeves.iter().enumerate().map(|(i, sleeve)| {
        let conductivity = |layer: usize| section.layers()[layer].material
            .property("thermal_conductivity").map(|p| p.value_si())
            .ok_or("layer has no sourced thermal conductivity");
        let shape = sleeve.resistance_geometry_per_m;
        Ok((i as u64, series_conductance_w_k(sleeve.interface_area_m2,
            [shape[0]*sleeve.interface_area_m2, shape[1]*sleeve.interface_area_m2],
            [conductivity(i)?, conductivity(i+1)?])?))
    }).collect()
}

fn coefficients(
    endpoint: &DepositionSimulation<'_>, law: &ConcentricSectionTemperatureLaw<'_>,
    bonded: &input::Bonded, power_w: f64,
) -> Result<StepCoefficients<ConcentricElasticSection>> {
    let temperatures = (0..bonded.sleeves.len()+1)
        .map(|i| endpoint.temperature_k(CellId(i as u64)))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let section = law.at_temperatures_k(&temperatures)?;
    Ok(StepCoefficients {
        sources_w: BTreeMap::from([(CellId(bonded.sleeves.len() as u64), power_w)]),
        conductances_w_k: conductances(bonded, &section)?,
        payload: section,
    })
}

pub(super) fn simulate(
    cx: &Cx<'_>, spec: &input::Spec, cards: &[MaterialCard], mut gas: IsochoricGasParcel,
    mut reservoir: Option<ClosedGasParcel>,
) -> Result<History> {
    spec.validate()?;
    let valid_reservoir = match (&spec.nozzle, &reservoir) {
        (None, None) => true,
        (Some(nozzle), Some(parcel)) => parcel.volume_m3().to_bits() == nozzle.reservoir_volume_m3.to_bits()
            && (spec.gas_solve.minimum_temperature_k..=spec.gas_solve.maximum_temperature_k)
                .contains(&parcel.state().state().temperature),
        _ => false,
    };
    if !valid_reservoir {
        return Err("nozzle data requires its original finite reservoir and admitted initial temperature".into());
    }
    let bonded = spec.bonded.as_ref().ok_or("bonded string geometry is required")?;
    let count = bonded.sleeves.len()+1;
    if cards.len() != count
        || !(spec.gas_solve.minimum_temperature_k..=spec.gas_solve.maximum_temperature_k)
            .contains(&gas.state().state().temperature)
    { return Err("layer/card count or initial gas temperature is unsupported".into()); }
    let reference_point = context(bonded.reference_temperature_k, spec.query_omega_rad_s)?;
    let mut material_layers = Vec::with_capacity(count);
    for (i, card) in cards.iter().enumerate() {
        let mut requirements = Vec::new();
        for (name, dims, admission) in [
            ("density", Density::DIMS, ScalarAdmissibility::StrictlyPositive),
            ("young_modulus", Pressure::DIMS, ScalarAdmissibility::StrictlyPositive),
            (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS,
                ScalarAdmissibility::Finite),
            (KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY, DynViscosity::DIMS,
                ScalarAdmissibility::NonNegative),
            ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS,
                ScalarAdmissibility::StrictlyPositive),
            ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS,
                ScalarAdmissibility::NonNegative),
        ] { requirements.push(ScalarPropertyRequirement::try_new(name, dims, admission)?); }
        if spec.drive.is_some() {
            // The native conductor refuses absent or ambiguous bases; neither
            // a material name nor a missing resistivity implies insulation.
            for basis in [ElectricalPropertyBasis::Resistivity, ElectricalPropertyBasis::Conductivity] {
                if card.claims().registered_quantity(basis.property_name()).is_some() {
                    requirements.push(basis.requirement()?);
                }
            }
        }
        material_layers.push(ConcentricMaterialLayer {
            outer_radius_m: if i == 0 { bonded.core_radius_m } else { bonded.sleeves[i-1].outer_radius_m },
            material: resolve_material_state_point(card, &reference_point, &requirements,
                MaterialPropertySelection::SingleClaimOnly)?,
        });
    }
    let reference = ConcentricElasticSection::try_new(0.0, material_layers)?;
    let card_refs: Vec<_> = cards.iter().collect();
    let temperature_law = ConcentricSectionTemperatureLaw::try_new(
        reference.clone(), &card_refs, &vec!["T"; count])?;
    let law = ConcentricThermoelasticLaw::try_new(temperature_law.clone())?;
    let caloric = cards.iter().map(|card| caloric_curve(spec, card)).collect::<Result<Vec<_>>>()?;
    let initial_section = temperature_law.at_temperatures_k(&vec![spec.solid_temperature_k; count])?;
    let initial_links = conductances(bonded, &initial_section)?;
    let mut events = Vec::with_capacity(2*count-1);
    let mut materials: Vec<(MaterialId, &dyn ThermalLaw)> = Vec::with_capacity(count);
    let mut inner = 0.0;
    for (i, ((geometry, layer), chart)) in reference.layer_integrals().iter()
        .zip(reference.layers()).zip(&caloric).enumerate()
    {
        materials.push((MaterialId(i as u64), chart.curve()));
        events.push(TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
            id: CellId(i as u64), segment_id: 0, body_id: 0, material_id: MaterialId(i as u64),
            mass_kg: geometry.line_density_kg_m()*spec.length_m,
            specific_h_j_kg: chart.state_at_temperature(spec.solid_temperature_k)?.specific_enthalpy_j_kg(),
            geometry: CellGeometry { start_m: [0.0; 3], end_m: [0.0, 0.0, spec.length_m],
                volume_m3: geometry.area_m2()*spec.length_m,
                surface_area_m2: 2.0*std::f64::consts::PI*(inner+layer.outer_radius_m)*spec.length_m
                    +2.0*geometry.area_m2() },
            boundary: Boundary { ambient_k: spec.solid_temperature_k, htc_w_m2_k: 0.0,
                reservoir: None, source_w: 0.0 },
        }) });
        inner = layer.outer_radius_m;
    }
    for (i, sleeve) in bonded.sleeves.iter().enumerate() {
        events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
            id: i as u64, a: CellId(i as u64), b: CellId((i+1) as u64),
            area_m2: sleeve.interface_area_m2, conductance_w_k: initial_links[&(i as u64)],
        }) });
    }
    let mut thermal = DepositionSimulation::new(events, &materials, DepositionBudget {
        max_cells: count, max_links: count-1, max_events: 2*count-1, max_substeps: 1,
    })?.advance_to(0.0, StepConfig { max_dt_s: 1.0, energy_tolerance_j: 0.0 },
        || cx.checkpoint().map_err(|_| fs_thermal::enthalpy::DepositionError::Cancelled))?.0;
    let bindings = caloric.iter().enumerate().map(|(i, caloric)| ConcentricStringThermalCell {
        cell: CellId(i as u64), caloric,
    }).collect();
    let mut live = ThermalConcentricStringRuntime::try_new(cx, string_template(spec), law,
        bindings, &thermal, bonded.reference_length_m, bonded.strain_limit, None,
        *gas.state().state(), spec.listener_m, spec.sample_rate_hz)?;
    live.initialize_point_load_equilibrium(cx,
        &[StringPointForce { station_frac: spec.station, force_n: spec.force_n }])?;
    let initial_mechanical = live.mechanical_internal_energy_j(cx)?;
    let initial_thermal = thermal.energy_j();
    let initial_gas = gas.state().specific_internal_energy_j_kg();
    let gas_mass = gas.mass_kg();
    let initial_context = context(spec.solid_temperature_k, spec.query_omega_rad_s)?;
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let outer = count-1;
    let transfer = MaterialConductionLayers { material: &cards[outer], context: &initial_context,
        temperature_axis: "T", liquid_fraction_axis: None, conductivity: &conductivity,
        selection: MaterialPropertySelection::SingleClaimOnly, area_m2: spec.transfer_area_m2,
        gas_distance_m: spec.gas_distance_m, body_distance_m: spec.solid_distance_m };
    let controls = ConcentricFiniteGasSolve {
        thermal: bonded.solve,
        gas: FiniteGasStringSolve { gas: spec.gas_solve,
            maximum_gas_temperature_increment_k: spec.maximum_gas_increment_k },
        exchange: spec.exchange,
    };
    let dt = 1.0/f64::from(spec.sample_rate_hz);
    let mut history = History {
        description: format!(
            "bonded_layers={count}, material_cards={:?}, caloric_curves={:?}, gas_volume_m3={:e}; fixed circular annuli, scalar alpha, affine E(T), perfect mechanical bond, authored radial resistance geometry; heat drives outer layer; initial/endpoint electrical columns retain their times, Joule quadrature uses accepted endpoint; solid_temperature_k is mass-weighted; finite gas is a thermal inventory, compact pressure has no radiation load; source identity is not measured validation{}",
            cards.iter().map(|c| c.content_hash().to_hex()).collect::<Vec<_>>(),
            caloric.iter().map(|c| c.curve().identity().to_hex()).collect::<Vec<_>>(),
            gas.volume_m3(), thermochemical::provenance(&[GasRecord::Thermochemical(gas.state().clone())])),
        pressure_pa: Vec::with_capacity(spec.samples), records: Vec::with_capacity(spec.samples),
        layer_temperatures: Vec::with_capacity(spec.samples),
        gas_flow: reservoir.as_ref().map(|parcel| flow::History::new(&gas, parcel, spec.samples)),
        electrical_feedback: None,
    };
    if let (Some(nozzle), Some(parcel)) = (&spec.nozzle, &reservoir) {
        history.description.push_str(&format!(
            "; nozzle={:?}, reservoir_volume_m3={:e}; first-order heat-then-flow, finite well-mixed volumes, no jet momentum or chamber load; reservoir{}",
            nozzle.geometry, parcel.volume_m3(),
            thermochemical::provenance(&[GasRecord::Thermochemical(parcel.state().clone())])));
    }
    let (mut material_heat, mut drag_heat, mut exchange_heat, mut joule_heat) = (0.0, 0.0, 0.0, 0.0);
    let (mut solid_heat, mut gas_heat, mut work, mut other_loss, mut allowance) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut reversible, mut axial, mut entropy) = (0.0, 0.0, 0.0);
    let mut energy_guess = vec![0.0; count];
    for sample in 0..spec.samples {
        let forces = live.project_point_forces(cx, &[StringPointForce { station_frac: spec.station,
            force_n: if sample < spec.release_sample { spec.force_n } else { 0.0 } }])?;
        let before_energy: Vec<_> = thermal.cells().values().map(|c| c.mass_kg*c.specific_h_j_kg).collect();
        let (current, resistance) = if let Some(drive) = spec.drive {
            let conductor = live.section().elastic().thermal_transport()?.with_ohmic_transport(spec.length_m)?;
            (conductor.evaluate(drive, 0.0)?.current_a, 1.0/conductor.conductance_s())
        } else { (0.0, 0.0) };
        let (frame, flux) = if let (Some(nozzle), Some(parcel)) = (&spec.nozzle, reservoir.as_mut()) {
            let report = live.step_with_isochoric_gas_nozzle_conduction(cx, &mut thermal, &mut gas,
                parcel, controls, &energy_guess, &forces, spec.gas_power_w*dt, spec.drive,
                &transfer, nozzle.geometry, nozzle.controls(spec.gas_solve),
                |endpoint| coefficients(endpoint, &temperature_law, bonded, spec.solid_power_w))?;
            let sign = match report.direction {
                FiniteGasStringFlowDirection::IntoSurroundings => 1.0,
                FiniteGasStringFlowDirection::OutOfSurroundings => -1.0,
            };
            let flux = (sign*report.flow.mass_flow_kg_s, sign*report.transfer.transferred_mass_kg,
                sign*report.transfer.transferred_enthalpy_j, report.energy_tolerance_j);
            (report.thermal, Some(flux))
        } else {
            (live.step_with_isochoric_gas_conduction(cx, &mut thermal, &mut gas,
                controls, &energy_guess, &forces, spec.gas_power_w*dt, spec.drive, &transfer,
                |endpoint| coefficients(endpoint, &temperature_law, bonded, spec.solid_power_w))?, None)
        };
        for ((guess, before), cell) in energy_guess.iter_mut().zip(before_energy).zip(thermal.cells().values()) {
            *guess = cell.mass_kg*cell.specific_h_j_kg-before;
        }
        let vibration = &frame.string.vibration;
        material_heat += vibration.dissipation.material_heat_j;
        drag_heat += frame.drag_heat_j;
        exchange_heat += frame.solid_to_gas_heat_j;
        let (endpoint_current, endpoint_resistance) = if let Some(electrical) = &frame.string.source.electrical {
            joule_heat += electrical.joule_heat_j;
            (electrical.electrical.current_a,
                1.0/live.section().elastic().thermal_transport()?.with_ohmic_transport(spec.length_m)?.conductance_s())
        } else { (0.0, 0.0) };
        solid_heat += spec.solid_power_w*dt;
        gas_heat += spec.gas_power_w*dt;
        work += vibration.acoustic.input_work_j;
        other_loss += vibration.dissipation.authored_loss_j+vibration.dissipation.roundoff_residual_j;
        allowance += flux.as_ref().map_or(frame.energy_tolerance_j, |flux| flux.3);
        for layer in &frame.string.layers {
            reversible += layer.reversible_heat_j;
            axial += layer.axial_free_energy_change_j;
            entropy += layer.mechanical_entropy_energy_change_j;
        }
        let mut stored = live.mechanical_internal_energy_j(cx)?-initial_mechanical
            +thermal.energy_j()-initial_thermal
            +gas_mass*(gas.state().specific_internal_energy_j_kg()-initial_gas);
        let mut gas_reference_scale = 0.0;
        if let (Some(gas_flow), Some(parcel)) = (&history.gas_flow, &reservoir) {
            let changes = gas_flow.energy_changes(&gas, parcel);
            stored = live.mechanical_internal_energy_j(cx)?-initial_mechanical
                +thermal.energy_j()-initial_thermal+changes[0]+changes[1];
            gas_reference_scale = gas_flow.reference_energy_scale(&gas, parcel);
        }
        let defect = stored+other_loss-work-joule_heat-solid_heat-gas_heat;
        let mut tolerance = allowance+64.0*f64::EPSILON*(sample+1) as f64
            *(stored.abs()+other_loss.abs()+work.abs()+joule_heat.abs()+solid_heat.abs()
                +gas_heat.abs()+initial_mechanical.abs()+initial_thermal.abs());
        if history.gas_flow.is_some() {
            tolerance += 64.0*f64::EPSILON*(sample+1) as f64*gas_reference_scale;
        }
        if !defect.is_finite() || !tolerance.is_finite() || defect.abs() > tolerance
            || vibration.epoch != live.accepted_samples()
            || live.accepted_samples() != (sample+1) as u64 || live.gas() != gas.state().state()
        { return Err("layered string/gas trajectory violates shared clock or cumulative energy".into()); }
        let temperatures = (0..count).map(|i| thermal.temperature_k(CellId(i as u64)))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        let mean = thermal.cells().values().zip(&temperatures)
            .map(|(cell, t)| cell.mass_kg*t).sum::<f64>()/thermal.mass_kg();
        let motion = live.point_motion(cx, spec.station)?;
        let pressure = vibration.acoustic.observer_pressure_pa;
        let row = [(sample+1) as f64*dt, pressure, motion.displacement_m, motion.velocity_m_s,
            mean, gas.state().state().temperature, gas.state().state().pressure,
            material_heat, drag_heat, exchange_heat, joule_heat, solid_heat, gas_heat,
            work, stored, defect, tolerance, current, endpoint_current, resistance, endpoint_resistance,
            live.modes()[0].damping_ratio, live.string().tension_n,
            live.modes()[0].angular_frequency_rad_s/(2.0*std::f64::consts::PI),
            reversible, axial, entropy];
        if row.iter().chain(&temperatures).any(|v| !v.is_finite()) {
            return Err("nonfinite layered trajectory observation".into());
        }
        if let (Some(gas_flow), Some(parcel), Some(flux)) =
            (&mut history.gas_flow, &reservoir, flux)
        { gas_flow.record(&gas, parcel, flux.0, flux.1, flux.2)?; }
        history.records.push(row);
        history.pressure_pa.push(pressure);
        history.layer_temperatures.push(temperatures);
    }
    cx.checkpoint().map_err(|_| "layered string cancelled")?;
    Ok(history)
}
