//! Describe the actual coupled graph, not a guessed single inlet-to-outlet path.
use fs_couple::bernoulli_aperture::network::NetworkNode;
use fs_couple::bernoulli_aperture::performance::{ApertureObservation,CoupledAperture};
use fs_couple::bernoulli_aperture::performance::file::PlateValvePerformance;
use fs_vfit::impedance::SeriesImpedanceSpec;

pub(super) fn computed_baffled_provenance(
    record: &fs_couple::bernoulli_aperture::radiation::environment::multiport::source::ComputedBaffledRecord,
) -> String {
    let options = record.options();
    let limits = record.limits();
    let band = record.band_hz();
    let report = record.report();
    let apertures = record.apertures().iter().enumerate().map(|(port, aperture)| format!(
        "{{\"port\":{port},\"center_m\":[{:e},{:e}],\"radius_m\":{:e}}}",
        aperture.center_m[0], aperture.center_m[1], aperture.radius_m
    )).collect::<Vec<_>>().join(",");
    let number = |value: Option<f64>| value.filter(|v| v.is_finite())
        .map_or_else(|| "null".into(), |value| format!("{value:e}"));
    let rungs = report.rungs.iter().map(|rung| format!(
        "{{\"order\":{},\"fit_relative_error\":{},\"held_out_relative_error\":{},\"continuous_passive\":{},\"sampled_passive\":{},\"refusal\":{}}}",
        rung.order, number(rung.fit_relative_error), number(rung.held_out_relative_error),
        rung.continuous_passive, rung.sampled_passive,
        rung.refusal.as_ref().map_or_else(|| "null".into(), |text| super::json_string(text))
    )).collect::<Vec<_>>().join(",");
    let gas_sources = super::thermochemical::provenance(core::slice::from_ref(record.gas_record()));
    format!("{{\"source\":\"computed-baffled\",\"declaration_identity\":\"{}\",\"realization_identity\":\"{}\",\"exterior_gas\":{},\"band_hz\":[{:e},{:e}],\"ordered_apertures\":[{apertures}],\"frequency_samples\":{},\"orders\":[{},{},{}],\"fit_iterations\":{},\"maximum_relative_error\":{:e},\"max_fit_memory_bytes\":{},\"causality_subtractions\":{},\"causality_tail_bound_ratio\":{:e},\"causality_data_tolerance\":{:e},\"max_states\":{},\"max_memory_bytes\":{},\"sampled_passivity_points\":{},\"selected_order\":{},\"held_out_samples\":{},\"fit_relative_error\":{:e},\"held_out_relative_error\":{:e},\"maximum_causality_ratio\":{:e},\"continuous_passivity_class\":\"{:?}\",\"sampled_passivity_class\":\"{:?}\",\"rungs\":[{rungs}],\"scope\":\"actual coplanar circular openings and resolved exterior; shared closed-form piston and Rayleigh mutual operator; proper undelayed rational approximation over full mechanical Nyquist band; no passivity repair, delayed realization or measured-instrument qualification; sampled GridOnly is not a storage certificate\"{gas_sources}}}",
        fs_blake3::ContentHash(record.declaration_identity()).to_hex(),
        fs_blake3::ContentHash(record.realization_identity()).to_hex(),
        gas_state(record.gas_record().state()), band[0], band[1],
        options.frequency_samples, options.orders[0], options.orders[1], options.orders[2],
        options.fit_iterations, options.maximum_relative_error, options.max_fit_memory_bytes,
        options.causality.subtractions, options.causality.tail_bound_ratio,
        options.causality.relative_data_tolerance, limits.max_states, limits.max_memory_bytes,
        limits.sampled_passivity_points, report.selected_order, report.held_out_samples,
        report.maximum_fit_relative_error, report.maximum_held_out_relative_error,
        report.maximum_causality_ratio, report.continuous.class, report.sampled.class)
}

// Both native wind loaders describe the same original source artifact. Its
// byte identity is deliberately separate from caller-declared bake identities.
pub(super) fn mimo_radiation_provenance(
    record: &fs_couple::bernoulli_aperture::radiation::environment::multiport::input::MimoRadiationRecord,
) -> String {
    let limits = record.limits();
    let band = record.band_hz();
    let apertures = record.apertures().iter().enumerate().map(|(port, aperture)| format!(
        "{{\"port\":{port},\"center_m\":[{:e},{:e}],\"radius_m\":{:e}}}",
        aperture.center_m[0], aperture.center_m[1], aperture.radius_m
    )).collect::<Vec<_>>().join(",");
    let gas_sources = super::thermochemical::provenance(core::slice::from_ref(record.gas_record()));
    format!("{{\"schema\":\"frankensim-mimo-radiation-v1\",\"path\":{},\"artifact_blake3\":\"{}\",\"declaration_identity\":\"{}\",\"realization_identity\":\"{}\",\"exterior_gas\":{},\"band_hz\":[{:e},{:e}],\"ordered_apertures\":[{apertures}],\"max_states\":{},\"max_memory_bytes\":{},\"sampled_passivity_points\":{},\"source_convention\":\"proper undelayed reciprocal impedance in Pa s/m^3; s=+i omega\",\"work_scope\":\"signed terminal exchange, not dissipation or certified rational stored energy\",\"scope\":\"original externally produced operator, shared mean-recovery ZOH terminal host and coplanar infinite-baffle observer; declared bake identities are retained, not attested; continuous and sampled GridOnly checks do not qualify a physical instrument\"{gas_sources}}}",
        super::json_string(record.path()), record.artifact_hash().to_hex(),
        fs_blake3::ContentHash(record.declaration_identity()).to_hex(),
        fs_blake3::ContentHash(record.realization_identity()).to_hex(),
        gas_state(record.gas_record().state()), band[0], band[1],
        limits.max_states, limits.max_memory_bytes, limits.sampled_passivity_points)
}

fn impedance(load:SeriesImpedanceSpec)->String {
    let compliance=load.compliance_m3_pa.map_or_else(||"null".into(),|v|format!("{v:e}"));
    format!("\"resistance_pa_s_m3\":{:e},\"inertance_pa_s2_m3\":{:e},\"compliance_m3_pa\":{compliance}",
        load.resistance_pa_s_m3,load.inertance_pa_s2_m3)
}

pub(super) fn provenance(p:&PlateValvePerformance)->String {
    // Report the resolved default even on an unsegmented tube. Source identity
    // already binds the exact authored species, basis and fractions.
    let exterior = p.exterior_gas().map_or_else(String::new, |gas|
        format!(",\"exterior_gas\":{}", gas_state(gas)));
    let slit = p.renderer().system().aperture().viscous_slit().map_or_else(String::new, |law|
        format!(",\"slit_flow\":{{\"model\":\"spatial-laminar-bernoulli-v1\",\"length_m\":{:e},\"dynamic_viscosity_pa_s\":{:e},\"viscosity_source\":\"inlet-gas\",\"scope\":\"quasi-steady Newtonian slit with separated outflow; not resolved channel flow or pressure recovery\"}}",
            law.length_m(), law.dynamic_viscosity_pa_s()));
    let force_json = format!(",\"ambient_gas\":{}{exterior}{slit}{}", gas_state(p.ambient_gas()), force_provenance(p));
    let CoupledAperture::Network(n)=p.renderer().system() else {return force_json};
    let i=p.info();
    // Name the node even in a one-section graph: node zero need not be the inlet.
    // The numeric-reflection legacy tube never enters this graph path.
    let selected=match p.renderer().observation() {
        ApertureObservation::NetworkNode(node)|ApertureObservation::NetworkBaffled{node,..}=>node,
        ApertureObservation::Inlet=>n.spec().nodes.iter().position(|k|matches!(k,NetworkNode::Inlet))
            .expect("unique admitted inlet"),
        _=>unreachable!("network input rejects ambiguous tube observers"),
    };
    let nodes:Vec<_>=n.spec().nodes.iter().enumerate().map(|(index,node)| {
        let (kind,extra)=match *node {
            NetworkNode::Inlet=>("inlet",String::new()),
            NetworkNode::Driven=>("driven",String::new()),
            NetworkNode::Junction=>("junction",String::new()),
            NetworkNode::Termination{reflection}=>("reflection",format!(",\"reflection\":{reflection:e}")),
            NetworkNode::Impedance{load}=>("impedance",format!(",{}",impedance(load))),
            NetworkNode::ShuntAdmittance{load}=> {
                let terms:Vec<_>=load.terms().map(|t|format!(
                    "{{\"conductance_m3_pa_s\":{:e},\"rate_per_s\":{:e}}}",t.conductance_m3_pa_s,t.rate_per_s)).collect();
                ("shunt-admittance",format!(",\"conductance_m3_pa_s\":{:e},\"compliance_m3_pa\":{:e},\"relaxation_terms\":[{}]",
                    load.conductance_m3_pa_s(),load.compliance_m3_pa(),terms.join(",")))
            }
            NetworkNode::Relaxation{load}|NetworkNode::Series{load}|NetworkNode::Shunt{load}=> {
                let kind=match node {NetworkNode::Series{..}=>"series",NetworkNode::Shunt{..}=>"shunt",_=>"relaxation"};
                let terms:Vec<_>=load.terms().iter().map(|t|format!(
                    "{{\"resistance_pa_s_m3\":{:e},\"rate_per_s\":{:e}}}",t.resistance_pa_s_m3,t.rate_per_s)).collect();
                (kind,format!(",{},\"relaxation_terms\":[{}]",impedance(load.base()),terms.join(",")))
            }
        };
        let radiation=p.radiation_loads().iter().find(|(node,_)|*node==index)
            .map_or_else(String::new,|(_,load)|format!(",\"baffled_radiation_band_hz\":{:e}",load.maximum_frequency_hz()));
        format!("{{\"node\":{index},\"kind\":\"{kind}\"{extra}{radiation}}}")
    }).collect();
    let sections:Vec<_>=n.spec().sections.iter().zip(n.represented_sections()).enumerate().map(|(index,(s,r))| {
        let gas=n.section_gases().map_or_else(String::new, |gases|
            format!(",\"gas\":{}", gas_state(&gases[index])));
        format!("{{\"nodes\":[{},{}],\"radius_m\":{:e},\"requested_length_m\":{:e},\"represented_length_m\":{:e},\"one_way_mechanical_samples\":{},\"impedance_pa_s_m3\":{:e}{gas}}}",
            s.nodes[0],s.nodes[1],s.radius_m,s.length_m,r.represented_length_m,r.one_way_samples,r.impedance_pa_s_m3)
    }).collect();
    let loss_sections:Vec<_>=p.viscothermal_losses().iter().map(|r| {
        format!("{{\"source_section\":{},\"source_nodes\":[{},{}],\"source_length_m\":{:e},\"source_radius_m\":{:e},\"minimum_frequency_hz\":{:e},\"maximum_frequency_hz\":{:e},\"cells\":{},\"arms_per_load\":8,\"original_one_way_samples\":{},\"represented_length_m\":{:e},\"checked_max_complex_relative_error\":{:e},\"checked_max_loss_relative_error\":{:e},\"checked_max_scattering_error\":{:e},\"loss_node_range\":[{},{}],\"propagation_section_range\":[{},{}]}}",
            r.source.section,r.original_nodes[0],r.original_nodes[1],r.requested_length_m,r.loss.radius_m(),
            r.source.minimum_frequency_hz,r.source.maximum_frequency_hz,r.source.cells,
            r.one_way_samples,r.represented_length_m,r.loss.max_complex_error(),r.loss.max_real_loss_error(),
            r.max_scattering_error,r.node_range[0],r.node_range[1],r.section_range[0],r.section_range[1])
    }).collect();
    let loss_json=if loss_sections.is_empty() {String::new()} else {format!(
        ",\"viscothermal\":{{\"model\":\"wide-tube-zk-passive-rl-rc-v1\",\"source_sections\":[{}],\"extra_inviscid_inertia_compliance\":false,\"scope\":\"sampled finite-band first-order boundary layers; not DC, Poiseuille or thermal evolution\"}}",loss_sections.join(","))};
    format!("{force_json},\"duct_network\":{{\"node_count\":{},\"section_count\":{},\"radiation_terminal_count\":{},\"observed_node\":{selected},\"total_represented_section_length_m\":{:e},\"nodes\":[{}],\"sections\":[{}]{loss_json},\"scope\":\"one coupled graph; lossless propagating sections and explicit local loads; exterior output selects one outlet, not a sum or mutual exterior radiation model\"}}",
        i.duct_nodes,i.duct_sections,i.radiation_terminals,i.represented_tube_length_m,nodes.join(","),sections.join(","))
}

// Both standalone and ensemble output call this source-description seam.
fn force_provenance(p: &PlateValvePerformance) -> String {
    use fs_couple::bernoulli_aperture::dynamic::force::PlateForceFootprint;
    if p.info().force_ports == 0 { return String::new(); }
    let ports: Vec<_> = p.renderer().force_ports().iter().map(|port| {
        let site = match port.footprint() {
            PlateForceFootprint::Node(node) => format!("\"node\":{node}"),
            PlateForceFootprint::Patch(triangles) => format!("\"triangles\":{:?}",triangles),
        };
        format!("{{{site},\"force_velocity_projection\":{:e}}}",port.coefficient())
    }).collect();
    let internal = p.info().radiation_load.is_none() && !matches!(p.renderer().observation(),
        ApertureObservation::TubeBaffled(_) | ApertureObservation::NetworkBaffled {..});
    let scope = if internal { ",\"observation_scope\":\"internal coupled tube pressure; not an exterior microphone\"" } else { "" };
    format!("{scope},\"mechanical_forces\":{{\"ports\":[{}],\"events\":{},\"clock\":\"mechanical-sample-before-step\",\"initial_force_n\":0,\"scope\":\"prescribed total transverse forces on original plate; not lip-body or moving-lay dynamics\"}}",
        ports.join(","),p.info().force_controls)
}

// One description for standalone wind output and every reed ensemble source.
pub(super) fn reed_provenance(p: &fs_couple::render::schedule::reed::ReedPerformance) -> String {
    let (point, scope, observation_scope, receiver) = match p.outlet_receiver() {
        None => ("bore-jet-proxy",
            "authored reed primitives and static duct; bore-pressure plus compact-jet proxy; not a calibrated exterior microphone or material-identification claim",
            "bore-pressure-plus-compact-jet-proxy; not exterior microphone", String::new()),
        Some(info) => {
            let r = info.receiver;
            ("baffled-outlet",
                "terminal-only Rayleigh pressure; same regional TMM and flanged feedback load; finite-history sampled approximation, not a certified broadband radiation match; no side-hole, jet, structural or room radiation; not measured calibration",
                "terminal-only Rayleigh pressure in the bound exterior gas",
                format!(",\"outlet_receiver\":{{\"position_m\":[{:e},{:e},{:e}],\"radius_m\":{:e},\"gas\":{},\"radial_rings\":{},\"angular_points\":{},\"maximum_frequency_hz\":{:e},\"fft_samples\":{},\"axial_delay_mechanical_samples\":{},\"exterior_delay_mechanical_samples\":[{},{}],\"priming_source_sample\":-1,\"delay_compensated\":false}}",
                    r.position_m[0],r.position_m[1],r.position_m[2],info.radius_m,
                    gas_state(&info.gas),r.radial_rings,r.angular_points,r.maximum_frequency_hz,
                    info.fft_samples,info.axial_delay_samples,info.exterior_delay_samples.0,
                    info.exterior_delay_samples.1))
        }
    };
    let regions: Vec<_> = p.segment_gases().iter().map(gas_state).collect();
    let field = reed_field_provenance(p);
    format!(",\"pressure_point\":\"{point}\",\"model_scope\":\"{scope}\",\"observation_scope\":\"{observation_scope}\",\"ambient_gas\":{},\"segment_gases\":[{}]{receiver}{field}",
        gas_state(p.ambient_gas()),regions.join(","))
}

// Serialize actual admitted states only; never re-resolve names or reconstruct
// air from density/speed. The same snapshot vocabulary covers ambient and regions.
fn gas_state(g: &fs_material::gas::GasState) -> String {
    format!("{{\"temperature_k\":{:e},\"static_pressure_pa\":{:e},\"density_kg_m3\":{:e},\"sound_speed_m_s\":{:e},\"dynamic_viscosity_pa_s\":{:e},\"thermal_conductivity_w_m_k\":{:e},\"gamma\":{:e},\"prandtl\":{:e},\"water_mole_fraction\":{:e},\"specific_heat_cp_j_kg_k\":{:e},\"specific_gas_constant_j_kg_k\":{:e},\"characteristic_impedance_pa_s_m\":{:e}}}",
        g.temperature, g.pressure, g.density, g.sound_speed, g.dynamic_viscosity,
        g.thermal_conductivity, g.gamma, g.prandtl, g.water_mole_fraction,
        g.specific_heat_cp, g.specific_gas_constant, g.characteristic_impedance)
}

fn reed_field_provenance(p: &fs_couple::render::schedule::reed::ReedPerformance) -> String {
    let Some((field, _)) = p.gas_field() else { return String::new(); };
    use fs_couple::render::schedule::reed::GasFieldResolution;
    let (resolution, scope) = match p.gas_field_resolution().expect("field retains its resolution") {
        GasFieldResolution::Midpoint(cells) => (
            format!(",\"propagation\":\"midpoint-constant-v1\",\"cells_per_segment\":{cells},\"resolved_cells\":{}", p.segment_gases().len()),
            "frozen primitive temperature/composition field; EOS and transport reevaluated per midpoint; exact inlet and outlet states; same regional feedback and microphone; no mean-flow acoustic correction",
        ),
        GasFieldResolution::Continuous(steps) => (
            format!(",\"propagation\":\"continuous-magnus-v1\",\"steps_per_smooth_span\":{},\"cells_per_segment\":null,\"resolved_cells\":null,\"authored_segment_samples\":{}", steps.0, p.segment_gases().len()),
            "frozen primitive temperature/composition field; EOS and transport reevaluated at local Magnus quadrature points; exact inlet and outlet states; same continuous feedback and microphone; no mean-flow acoustic correction",
        ),
    };
    let species: Vec<_> = field.model().thermodynamics().components().iter()
        .map(|component| super::json_string(component.species().as_str())).collect();
    let stations: Vec<_> = field.stations_m().iter().zip(field.primitive_states())
        .map(|(x, state)| format!(
            "{{\"x_m\":{x:e},\"temperature_k\":{:e},\"mole_fractions\":{:?}}}",
            state.temperature_k, state.composition.fractions())).collect();
    let packs: Vec<_> = field.model().source_packs().iter().map(|receipt| format!(
        "{{\"species\":{},\"species_pack\":\"{}\",\"nasa9_pack\":\"{}\",\"transport_pack\":\"{}\"}}",
        super::json_string(receipt.species.as_str()), receipt.species_pack.to_hex(),
        receipt.nasa9_pack.to_hex(), receipt.transport_pack.to_hex())).collect();
    let steady = reed_steady_provenance(p);
    format!(",\"axial_gas_field\":{{\"static_pressure_pa\":{:e}{resolution},\"canonical_species\":[{}],\"stations\":[{}],\"source_packs\":[{}]{steady},\"scope\":\"{scope}\"}}",
        field.static_pressure_pa(), species.join(","),
        stations.join(","), packs.join(","))
}

fn reed_steady_provenance(p: &fs_couple::render::schedule::reed::ReedPerformance) -> String {
    let Some(generated) = p.steady_gas() else { return String::new(); };
    let problem = generated.problem();
    let solution = generated.solution();
    let controls = generated.controls();
    let cells: Vec<_> = problem.cells.iter().map(|cell| format!(
        "{{\"end_m\":{:e},\"volume_m3\":{:e},\"wall_area_m2\":{:e},\"wall_conductance_w_k\":{:e},\"wall_temperature_k\":{:e},\"initial_temperature_k\":{:e}}}",
        cell.end_m, cell.volume_m3, cell.wall_area_m2, cell.wall_conductance_w_k,
        cell.wall_temperature_k, cell.initial_temperature_k)).collect();
    let injections: Vec<_> = problem.injections.iter().map(|inlet| format!(
        "{{\"cell\":{},\"mass_flow_kg_s\":{:e},\"temperature_k\":{:e},\"mole_fractions\":{:?}}}",
        inlet.cell, inlet.mass_flow_kg_s, inlet.state.temperature_k,
        inlet.state.composition.fractions())).collect();
    format!(",\"steady_generation\":{{\"model\":\"shared-advected-enthalpy-species-v1\",\"inlet_mass_flow_kg_s\":{:e},\"inlet_temperature_k\":{:e},\"inlet_mole_fractions\":{:?},\"cells\":[{}],\"injections\":[{}],\"mass_flow_kg_s\":{:?},\"cell_power_residuals_w\":{:?},\"maximum_power_residual_bound_w\":{:e},\"power_tolerance_w\":{:e},\"maximum_energy_residual_j\":{:e},\"energy_tolerance_j\":{:e},\"caloric_residual_tolerance_j_kg\":{:e},\"temperature_bracket_k\":[{:e},{:e}],\"substeps\":{},\"maximum_substeps\":{},\"relaxation_time_s\":{:e},\"maximum_mean_mach\":{:e},\"mean_mach_limit\":{:e},\"scope\":\"source-bound stationary donor-cell heat and species balance; geometry from the same unbranched bore; prescribed total wall conductance; frozen acoustic field without mean-flow correction\"}}",
        problem.inlet_mass_flow_kg_s, problem.inlet.temperature_k,
        problem.inlet.composition.fractions(), cells.join(","), injections.join(","),
        solution.mass_flow_kg_s(), solution.cell_power_residuals_w(),
        solution.maximum_power_residual_bound_w(), controls.steady_power_tolerance_w,
        solution.maximum_energy_residual_j(), controls.energy_tolerance_j,
        controls.caloric.energy_tolerance_j_kg, controls.caloric.minimum_temperature_k,
        controls.caloric.maximum_temperature_k, solution.substeps(), controls.maximum_substeps,
        solution.relaxation_time_s(), solution.maximum_mean_mach(), controls.maximum_mean_mach)
}
