//! Describe the actual coupled graph, not a guessed single inlet-to-outlet path.
use fs_couple::bernoulli_aperture::network::NetworkNode;
use fs_couple::bernoulli_aperture::performance::{ApertureObservation,CoupledAperture};
use fs_couple::bernoulli_aperture::performance::file::PlateValvePerformance;
use fs_vfit::impedance::SeriesImpedanceSpec;

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
    let force_json = format!(",\"ambient_gas\":{}{exterior}{}", gas_state(p.ambient_gas()), force_provenance(p));
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
    format!(",\"pressure_point\":\"{point}\",\"model_scope\":\"{scope}\",\"observation_scope\":\"{observation_scope}\",\"ambient_gas\":{},\"segment_gases\":[{}]{receiver}",
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
