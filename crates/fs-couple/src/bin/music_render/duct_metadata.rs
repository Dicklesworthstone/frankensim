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
    let CoupledAperture::Network(n)=p.renderer().system() else {return String::new()};
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
    let sections:Vec<_>=n.spec().sections.iter().zip(n.represented_sections()).map(|(s,r)|format!(
        "{{\"nodes\":[{},{}],\"radius_m\":{:e},\"requested_length_m\":{:e},\"represented_length_m\":{:e},\"one_way_mechanical_samples\":{},\"impedance_pa_s_m3\":{:e}}}",
        s.nodes[0],s.nodes[1],s.radius_m,s.length_m,r.represented_length_m,r.one_way_samples,r.impedance_pa_s_m3)).collect();
    let loss_sections:Vec<_>=p.viscothermal_losses().iter().map(|r| {
        format!("{{\"source_section\":{},\"source_nodes\":[{},{}],\"source_length_m\":{:e},\"source_radius_m\":{:e},\"minimum_frequency_hz\":{:e},\"maximum_frequency_hz\":{:e},\"cells\":{},\"arms_per_load\":8,\"original_one_way_samples\":{},\"represented_length_m\":{:e},\"checked_max_complex_relative_error\":{:e},\"checked_max_loss_relative_error\":{:e},\"checked_max_scattering_error\":{:e},\"loss_node_range\":[{},{}],\"propagation_section_range\":[{},{}]}}",
            r.source.section,r.original_nodes[0],r.original_nodes[1],r.requested_length_m,r.loss.radius_m(),
            r.source.minimum_frequency_hz,r.source.maximum_frequency_hz,r.source.cells,
            r.one_way_samples,r.represented_length_m,r.loss.max_complex_error(),r.loss.max_real_loss_error(),
            r.max_scattering_error,r.node_range[0],r.node_range[1],r.section_range[0],r.section_range[1])
    }).collect();
    let loss_json=if loss_sections.is_empty() {String::new()} else {format!(
        ",\"viscothermal\":{{\"model\":\"wide-tube-zk-passive-rl-rc-v1\",\"source_sections\":[{}],\"extra_inviscid_inertia_compliance\":false,\"scope\":\"sampled finite-band first-order boundary layers; not DC, Poiseuille or thermal evolution\"}}",loss_sections.join(","))};
    format!(",\"duct_network\":{{\"node_count\":{},\"section_count\":{},\"radiation_terminal_count\":{},\"observed_node\":{selected},\"total_represented_section_length_m\":{:e},\"nodes\":[{}],\"sections\":[{}]{loss_json},\"scope\":\"one coupled graph; lossless propagating sections and explicit local loads; exterior output selects one outlet, not a sum or mutual exterior radiation model\"}}",
        i.duct_nodes,i.duct_sections,i.radiation_terminals,i.represented_tube_length_m,nodes.join(","),sections.join(","))
}
