//! The authored finite inventories and powers, recorded before playback.
use super::json_string;
use fs_couple::bernoulli_aperture::dg::thermal_wall::FiniteWallThermalSolve;
use fs_couple::bernoulli_aperture::performance::AperturePerformance;
use fs_couple::bernoulli_aperture::performance::gas_wall::GasWallConduction;
use fs_couple::cavity_phs::piston::thermal::HeatExchangeSolve;
use fs_couple::dg_port::finite_gas::{FiniteGasDgLoss, FiniteGasDgSolve};
use fs_material::phase::UniformEnthalpyStepInput;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermochemical_gas::energy::GasEnergySolve;

pub(super) fn provenance(renderer: &AperturePerformance, clock_hz: u32) -> String {
    let Some(program) = renderer.finite_gas_heat_program() else { return String::new(); };
    let gas = renderer.finite_gas().expect("heat program retains its original inventory");
    let state = gas.state().state();
    let controls = program.controls.iter().map(|control| format!(
        "{{\"sample\":{},\"power_w\":{:e}}}", control.sample, control.power_w))
        .collect::<Vec<_>>().join(",");
    let budgets = gas_controls_json(program.caloric, program.coupling);
    format!(",\"finite_gas_heat\":{{\"clock_hz\":{clock_hz},\"phase\":\"frozen\",\"initial_mass_kg\":{:e},\"volume_m3\":{:e},\"initial_temperature_k\":{:e},\"initial_pressure_pa\":{:e},\"controls\":[{controls}],{budgets},\"scope\":\"signed external power held on the original physical clock; fixed-volume closed material inventory pays acoustic parameter work; all gas snapshots describe the initial state; no mean mass transport or computed exterior pressure jump\"}}",
        gas.mass_kg(), gas.volume_m3(), state.temperature, state.pressure)
}

/// Read the actual admitted owners before the renderer enters playback.
pub(super) fn wall_provenance(renderer: &AperturePerformance, clock_hz: u32) -> String {
    if let Some(program) = renderer.finite_wall_heat_program() {
        let body = body_json(renderer.finite_wall_body().expect("wall program retains its body"));
        let controls = program.controls.iter().map(|control| format!(
            "{{\"sample\":{},\"power_w\":{:e}}}", control.sample, control.power_w))
            .collect::<Vec<_>>().join(",");
        let coupling = wall_controls_json(program.coupling);
        return format!(",\"finite_wall_heat\":{{\"clock_hz\":{clock_hz},\"snapshot\":\"initial-before-playback\",\"body\":{body},\"controls\":[{controls}],\"coupling\":{coupling},\"scope\":\"signed external power held on the original physical clock; one all-solid body with fixed inventory, density and geometry contains the complete main-bore wall; body enthalpy pays native wall parameter work; optional native solid dissipation is deposited once; gas stays fixed and branch wall inventories remain independent\"}}");
    }
    let Some(program) = renderer.finite_gas_wall_heat_program() else { return String::new(); };
    let body = body_json(renderer.finite_gas_wall_body().expect("joint program retains its body"));
    let gas = renderer.finite_gas_wall_gas().expect("joint program retains its gas");
    let state = super::super::duct_metadata::gas_state(gas.state().state());
    let budgets = gas_controls_json(program.coupling.caloric, program.coupling.gas);
    let wall = wall_controls_json(program.coupling.wall);
    let exchange = exchange_json(program.coupling.exchange);
    let conduction = conduction_json(renderer.finite_gas_wall_conduction()
        .expect("joint program retains its original conduction source"));
    let controls = program.controls.iter().map(|control| format!(
        "{{\"sample\":{},\"gas_power_w\":{:e},\"wall_power_w\":{:e}}}",
        control.sample, control.gas_power_w, control.wall_power_w))
        .collect::<Vec<_>>().join(",");
    format!(",\"finite_gas_wall_heat\":{{\"clock_hz\":{clock_hz},\"snapshot\":\"initial-before-playback\",\"body\":{body},\"gas\":{{\"phase\":\"frozen\",\"mass_kg\":{:e},\"volume_m3\":{:e},\"specific_internal_energy_j_kg\":{:e},\"state\":{state},{budgets}}},\"controls\":[{controls}],\"wall_coupling\":{wall},\"interface\":{conduction},\"exchange\":{exchange},\"maximum_material_trials\":{},\"scope\":\"signed external powers held on the original physical clock; closed fixed-composition fixed-volume gas and one all-solid wall body with fixed inventory, density and geometry; original material and gas conductivities at their own energy midpoints determine reciprocal series conduction; interface heat is credited to gas and debited from body once; both inventories pay native parameter work; optional native solid heat is deposited once; no phase transfer, expansion, advection or resolved convection\"}}",
        gas.mass_kg(), gas.volume_m3(), gas.state().specific_internal_energy_j_kg(),
        program.coupling.maximum_material_trials)
}

fn body_json(body: UniformEnthalpyStepInput<'_>) -> String {
    let receipts = body.curve.source_receipts().iter().map(|receipt| {
        let sources = receipt.source_hashes.iter().map(|hash| format!("\"{}\"", hash.to_hex()))
            .collect::<Vec<_>>().join(",");
        let point = receipt.query_point.iter().map(|(axis, value)| {
            let quantity = receipt.axis_quantities.get(axis).map_or_else(|| "null".to_owned(),
                |quantity| format!("{:?}", quantity.canonical_bytes()));
            format!("{{\"axis\":{},\"value_si\":{value:e},\"quantity\":{quantity}}}", json_string(axis))
        }).collect::<Vec<_>>().join(",");
        format!("{{\"blake3\":\"{}\",\"property\":{},\"selected_claim\":\"{}\",\"policy\":{},\"observation_backed\":{},\"query\":[{point}],\"source_hashes\":[{sources}]}}",
            receipt.content_hash().to_hex(), json_string(&receipt.property),
            receipt.selected.0.to_hex(), json_string(receipt.policy), receipt.observation_backed)
    }).collect::<Vec<_>>().join(",");
    let knots = body.curve.knots();
    let first = knots.first().expect("admitted enthalpy chart has knots");
    let last = knots.last().expect("admitted enthalpy chart has knots");
    format!("{{\"mass_kg\":{:e},\"volume_m3\":{:e},\"surface_area_m2\":{:e},\"duration_s\":{:e},\"initial_temperature_k\":{:e},\"initial_specific_enthalpy_j_kg\":{:e},\"initial_density_kg_m3\":{:e},\"initial_liquid_mass_fraction\":{:e},\"material_card_blake3\":\"{}\",\"chart_blake3\":\"{}\",\"chart_knot_count\":{},\"chart_enthalpy_support_j_kg\":[{:e},{:e}],\"chart_source_receipts\":[{receipts}]}}",
        body.mass_kg, body.volume_m3, body.surface_area_m2, body.duration_s,
        body.initial.temperature_k(), body.initial.specific_enthalpy_j_kg(),
        body.initial.bulk_density_kg_m3(), body.initial.liquid_mass_fraction(),
        body.curve.material_card_identity().to_hex(), body.curve.identity().to_hex(),
        knots.len(), first.specific_enthalpy_j_kg, last.specific_enthalpy_j_kg)
}

fn conduction_json(conduction: &GasWallConduction) -> String {
    let context = conduction.context.axes().iter().map(|(axis, value)| {
        let quantity = conduction.context.axis_quantities().get(axis)
            .map_or_else(|| "null".to_owned(), |quantity| format!("{:?}", quantity.canonical_bytes()));
        format!("{{\"axis\":{},\"value_si\":{value:e},\"quantity\":{quantity}}}", json_string(axis))
    }).collect::<Vec<_>>().join(",");
    let selection = match &conduction.selection {
        MaterialPropertySelection::SingleClaimOnly => "\"single-claim-only\"".to_owned(),
        MaterialPropertySelection::PreferObservationBacked => "\"prefer-observation-backed\"".to_owned(),
        MaterialPropertySelection::PinnedByProperty(pins) => {
            let pins = pins.iter().map(|(property, claim)| format!(
                "{{\"property\":{},\"claim\":\"{}\"}}", json_string(property), claim.0.to_hex()))
                .collect::<Vec<_>>().join(",");
            format!("{{\"pinned_by_property\":[{pins}]}}")
        }
    };
    let liquid_axis = conduction.liquid_fraction_axis.as_deref()
        .map_or_else(|| "null".to_owned(), json_string);
    format!("{{\"model\":\"series-material-conduction\",\"area_m2\":{:e},\"gas_distance_m\":{:e},\"body_distance_m\":{:e},\"conductivity_property\":{},\"conductivity_quantity\":{:?},\"selection\":{selection},\"liquid_fraction_axis\":{liquid_axis},\"context\":[{context}]}}",
        conduction.area_m2, conduction.gas_distance_m, conduction.body_distance_m,
        json_string(conduction.conductivity.name()), conduction.conductivity.quantity().canonical_bytes())
}

fn exchange_json(controls: HeatExchangeSolve) -> String {
    format!("{{\"heat_bracket_j\":[{:e},{:e}],\"source_tolerance_j\":{:e},\"energy_tolerance_j\":{:e},\"maximum_iterations\":{}}}",
        controls.heat_bracket_j[0], controls.heat_bracket_j[1],
        controls.source_tolerance_j, controls.energy_tolerance_j, controls.maximum_iterations)
}

fn wall_controls_json(controls: FiniteWallThermalSolve) -> String {
    let exchange = exchange_json(controls.exchange);
    let heat = controls.dissipation_heat.map_or_else(|| "null".to_owned(), |budget| format!(
        "{{\"max_dimension\":{},\"max_dense_bytes\":{},\"max_multiply_terms\":{}}}",
        budget.max_dimension, budget.max_dense_bytes, budget.max_multiply_terms));
    format!("{{\"work_solve\":{exchange},\"maximum_absolute_parameter_work_j\":{:e},\"dissipation_heat\":{heat}}}",
        controls.maximum_absolute_parameter_work_j)
}

fn gas_controls_json(caloric: GasEnergySolve, coupling: FiniteGasDgSolve) -> String {
    let (loss, error) = match coupling.loss {
        FiniteGasDgLoss::Lossless => ("lossless", String::new()),
        FiniteGasDgLoss::UniformViscothermal { maximum_relative_error } =>
            ("uniform-viscothermal", format!(",\"maximum_relative_error\":{maximum_relative_error:e}")),
        FiniteGasDgLoss::SpatialViscothermal { maximum_relative_error } =>
            ("spatial-viscothermal", format!(",\"maximum_relative_error\":{maximum_relative_error:e}")),
        FiniteGasDgLoss::UniformMolecularDiffusion { maximum_relative_error } =>
            ("molecular-diffusion", format!(",\"maximum_relative_error\":{maximum_relative_error:e}")),
    };
    format!("\"caloric\":{{\"temperature_bracket_k\":[{:e},{:e}],\"energy_tolerance_j_kg\":{:e},\"maximum_iterations\":{}}},\"coupling\":{{\"work_bracket_j\":[{:e},{:e}],\"work_tolerance_j\":{:e},\"energy_tolerance_j\":{:e},\"maximum_iterations\":{},\"maximum_temperature_increment_k\":{:e},\"loss\":\"{loss}\"{error}}}",
        caloric.minimum_temperature_k, caloric.maximum_temperature_k,
        caloric.energy_tolerance_j_kg, caloric.maximum_iterations,
        coupling.work_bracket_j[0], coupling.work_bracket_j[1], coupling.work_tolerance_j,
        coupling.energy_tolerance_j, coupling.maximum_iterations,
        coupling.maximum_temperature_increment_k)
}
