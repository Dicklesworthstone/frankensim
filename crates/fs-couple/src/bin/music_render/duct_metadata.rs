//! Describe the actual coupled graph, not a guessed single inlet-to-outlet path.
use fs_couple::bernoulli_aperture::network::NetworkNode;
use fs_couple::bernoulli_aperture::performance::{ApertureObservation,CoupledAperture};
use fs_couple::bernoulli_aperture::performance::file::PlateValvePerformance;
use fs_vfit::impedance::SeriesImpedanceSpec;

/// Applied physical gas law and its actual bounded per-surface preparation.
/// A missing report keeps the legacy lossless observation explicit.
pub(super) fn observer_gas_provenance(
    report: Option<&fs_couple::pcm_wav::baffled::gas::GasPropagationReport>,
) -> String {
    use fs_couple::air_path::{GasAbsorptionModel, ObserverFilterReceipt};
    let Some(report) = report else { return "null".into(); };
    let law = match report.law {
        GasAbsorptionModel::StokesKirchhoff => "stokes-kirchhoff",
        GasAbsorptionModel::Iso9613Air { .. } => "iso9613",
        GasAbsorptionModel::MolecularRelaxation => "molecular-relaxation",
        GasAbsorptionModel::Lossless => "lossless",
    };
    // The finite native sidecar is a summary, not an unbounded second copy of
    // every path receipt. Complete per-path evidence remains in the report API.
    let bounds = |value: fn(&ObserverFilterReceipt) -> f64| {
        if report.filters.is_empty() { return "null".to_owned(); }
        let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
        for filter in &report.filters {
            let x = value(filter);
            low = low.min(x);
            high = high.max(x);
        }
        format!("[{low:e},{high:e}]")
    };
    let total_taps: usize = report.filters.iter().map(|filter| filter.taps).sum();
    let maximum_taps = report.filters.iter().map(|filter| filter.taps).max().unwrap_or(0);
    let maximum_truncated_energy = report.filters.iter().map(|filter| filter.truncated_energy)
        .fold(0.0_f64, f64::max);
    let maximum_budget_fraction = report.filters.iter().filter(|filter| filter.budget_energy > 0.0)
        .map(|filter| filter.truncated_energy / filter.budget_energy).fold(0.0_f64, f64::max);
    let molecular_relaxation_omitted = report.filters.iter()
        .any(|filter| filter.absorption.relaxation_omitted.is_some());
    let validity = report.filters.iter().filter_map(|filter| filter.classical_validity_exceeded_from_hz)
        .reduce(f64::min).map_or_else(|| "null".to_owned(), |frequency| format!("{frequency:e}"));
    let filters = format!("{{\"total_taps\":{total_taps},\"maximum_taps\":{maximum_taps},\"range_m\":{},\"maximum_truncated_energy\":{maximum_truncated_energy:e},\"budget_energy\":{},\"maximum_truncation_budget_fraction\":{maximum_budget_fraction:e},\"low_frequency_delay_s\":{},\"receiver_delay_s\":{},\"filter_latency_s\":{},\"classical_validity_exceeded_from_hz\":{validity},\"molecular_relaxation_omitted\":{molecular_relaxation_omitted}}}",
        bounds(|filter| filter.range_m), bounds(|filter| filter.budget_energy),
        bounds(|filter| filter.low_frequency_delay_s), bounds(|filter| filter.receiver_delay_s),
        bounds(|filter| filter.latency_s));
    format!("{{\"law\":\"{law}\",\"gas\":{},\"surface_paths\":{},\"unique_paths\":{},\"additional_latency_samples\":{},\"kernel_terms\":{},\"filter_summary\":{filters},\"scope\":\"physical exterior propagation on original quadrature paths; unchanged source and reaction; no bore bulk relaxation or free-field body scattering\"}}",
        gas_state(&report.gas), report.surface_filter_indices.len(), report.filters.len(),
        report.additional_latency_samples, report.kernel_terms)
}

pub(super) fn computed_baffled_provenance(
    record: &fs_couple::bernoulli_aperture::radiation::environment::multiport::source::ComputedBaffledRecord,
) -> String {
    let options = record.options();
    let limits = record.limits();
    let band = record.band_hz();
    let report = record.report();
    let driven_mouth = record.driven_mouth();
    let offset = usize::from(driven_mouth.is_some());
    let terminals = record.apertures().iter().enumerate().map(|(port, aperture)| format!(
        "{{\"port\":{},\"center_m\":[{:e},{:e}],\"radius_m\":{:e}}}",
        port + offset, aperture.center_m[0], aperture.center_m[1], aperture.radius_m
    )).collect::<Vec<_>>().join(",");
    let apertures = driven_mouth.map_or_else(|| terminals.clone(), |mouth| format!(
        "{{\"port\":0,\"shape\":\"rectangle\",\"center_m\":[{:e},{:e}],\"size_m\":[{:e},{:e}],\"cells\":[{},{}],\"outward_flow\":\"negative driven inlet flow\",\"terminal_radial_rings\":{},\"terminal_angular_points\":{}}},{terminals}",
        mouth.center_m[0], mouth.center_m[1], mouth.size_m[0], mouth.size_m[1],
        mouth.cells[0], mouth.cells[1], mouth.terminal_radial_rings, mouth.terminal_angular_points
    ));
    let number = |value: Option<f64>| value.filter(|v| v.is_finite())
        .map_or_else(|| "null".into(), |value| format!("{value:e}"));
    let rungs = report.rungs.iter().map(|rung| format!(
        "{{\"order\":{},\"fit_relative_error\":{},\"held_out_relative_error\":{},\"continuous_passive\":{},\"repair_relative_change\":{},\"sampled_passive\":{},\"refusal\":{}}}",
        rung.order, number(rung.fit_relative_error), number(rung.held_out_relative_error),
        rung.continuous_passive, number(rung.repair_relative_change), rung.sampled_passive,
        rung.refusal.as_ref().map_or_else(|| "null".into(), |text| super::json_string(text))
    )).collect::<Vec<_>>().join(",");
    let gas_sources = super::thermochemical::provenance(core::slice::from_ref(record.gas_record()));
    let source = if driven_mouth.is_some() { "computed-mixed-baffled" } else { "computed-baffled" };
    // An unrepaired fit keeps its original scope text byte for byte.
    let repair = report.rungs.iter().find(|rung| rung.order == report.selected_order)
        .and_then(|rung| rung.repair_relative_change)
        .map(|change| format!("residue passivity repair of the selected fit (in-band change {change:e}); no"));
    let scope = if driven_mouth.is_some() {
        format!("actual coplanar rectangular mouth and circular terminals in the bound exterior; complete shared power-consistent Rayleigh surface matrix and mutual source feedback; proper rational approximation over full mechanical Nyquist band; {} measured-instrument qualification; finite surface/frequency checks and sampled GridOnly are not continuum or storage certificates",
            repair.as_deref().unwrap_or("no passivity repair or"))
    } else {
        format!("actual coplanar circular openings and resolved exterior; shared closed-form piston and Rayleigh mutual operator; proper undelayed rational approximation over full mechanical Nyquist band; {} delayed realization or measured-instrument qualification; sampled GridOnly is not a storage certificate",
            repair.as_deref().unwrap_or("no passivity repair,"))
    };
    format!("{{\"source\":\"{source}\",\"declaration_identity\":\"{}\",\"realization_identity\":\"{}\",\"exterior_gas\":{},\"band_hz\":[{:e},{:e}],\"ordered_apertures\":[{apertures}],\"frequency_samples\":{},\"orders\":[{},{},{}],\"fit_iterations\":{},\"maximum_relative_error\":{:e},\"max_fit_memory_bytes\":{},\"causality_subtractions\":{},\"causality_tail_bound_ratio\":{:e},\"causality_data_tolerance\":{:e},\"max_states\":{},\"max_memory_bytes\":{},\"sampled_passivity_points\":{},\"selected_order\":{},\"held_out_samples\":{},\"fit_relative_error\":{:e},\"held_out_relative_error\":{:e},\"maximum_causality_ratio\":{:e},\"continuous_passivity_class\":\"{:?}\",\"sampled_passivity_class\":\"{:?}\",\"rungs\":[{rungs}],\"scope\":\"{scope}\"{gas_sources}}}",
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

// Read the same bound exterior used by an actual receiver. A DG reaction's
// retained record takes precedence over a separate outer-file declaration.
pub(super) fn aperture_exterior_gas(
    p: &PlateValvePerformance,
) -> Option<&fs_material::gas::GasState> {
    let bound = match p.renderer().system() {
        CoupledAperture::Dg(runtime) => runtime.load().radiation_record()
            .map(|record| record.gas_record().state()),
        _ => None,
    };
    bound.or_else(|| p.renderer().exterior_gas()).or_else(|| p.exterior_gas())
}

// The source record is meaningful even when the selected output is internal.
pub(super) fn dg_radiation_provenance(load: &fs_couple::dg_port::DgPortLoad) -> String {
    let mut output = load.radiation_record().map_or_else(String::new, |record|
        format!(",\"computed_radiation\":{}", computed_baffled_provenance(record)));
    output.push_str(&dg_molecular_provenance(load.duct()));
    if !load.terminal_channels().is_empty() {
        use fs_couple::dg_port::DgTerminalChannelKind;
        let channels = load.terminal_channels().iter().map(|bound| {
            let input = bound.input();
            let law = match input.kind {
                DgTerminalChannelKind::Closed => "\"kind\":\"closed\"".to_owned(),
                DgTerminalChannelKind::Slit {
                    width_m, height_m, length_m, transverse_cells, maximum_mach,
                } => format!("\"kind\":\"unsteady-fixed-slit\",\"width_m\":{width_m:e},\"height_m\":{height_m:e},\"length_m\":{length_m:e},\"transverse_cells\":{transverse_cells},\"maximum_mach\":{maximum_mach:e}"),
                DgTerminalChannelKind::Circular {
                    radius_m, length_m, radial_cells, maximum_mach,
                } => format!("\"kind\":\"unsteady-fixed-circular\",\"radius_m\":{radius_m:e},\"length_m\":{length_m:e},\"radial_cells\":{radial_cells},\"maximum_mach\":{maximum_mach:e}"),
            };
            format!("{{\"port\":{},{law},\"gas\":{}}}", input.port.0, gas_state(bound.gas()))
        }).collect::<Vec<_>>().join(",");
        output.push_str(&format!(",\"terminal_channels\":[{channels}],\"terminal_channel_scope\":\"same-step original transverse-momentum slit/circular states or exact zero-flow closures in the complete exterior matrix; additional fixed passage geometry and laminar port reduction; no moving-pad work, annular-gap resolution, edge shedding or turbulent-loss calibration\""));
    }
    if !load.duct().side_ports().is_empty() {
        let ports = load.terminal_ports().iter().map(|port| port.0.to_string())
            .collect::<Vec<_>>().join(",");
        let sides = load.duct().side_ports().iter().map(|port| {
            let source = load.duct().side_port_source(port.id)
                .expect("computed source admitted every original side-port gas");
            format!("{{\"port\":{},\"downstream_element\":{},\"station_m\":{:e},\"bore_area_m2\":{:e},\"gas\":{}}}",
                port.id.0, port.downstream_element, source.station_m,
                source.bore_area_m2, gas_state(&source.gas))
        }).collect::<Vec<_>>().join(",");
        let chimneys = load.branches().iter().map(|branch| {
            let duct = branch.duct();
            let mesh = duct.mesh();
            let band = duct.viscothermal_band_rad_s().map_or_else(|| "null".into(), |band|
                format!("[{:e},{:e}]", band[0], band[1]));
            let molecular = dg_molecular_provenance(duct);
            format!("{{\"port\":{},\"operator_identity\":\"{}\",\"length_m\":{:e},\"inlet_area_m2\":{:e},\"outlet_area_m2\":{:e},\"inlet_gas\":{},\"outlet_gas\":{},\"elements\":{},\"order\":{},\"total_dofs\":{},\"wall_band_rad_s\":{band}{molecular}}}",
                branch.port().0, duct.operator_identity().to_hex(), duct.total_length(),
                duct.inlet_area_m2(), duct.outlet_area_m2(), gas_state(duct.inlet_gas()),
                gas_state(duct.outlet_gas()), mesh.elements, mesh.order, duct.total_dofs())
        }).collect::<Vec<_>>().join(",");
        output.push_str(&format!(",\"terminal_port_order\":[{ports}],\"side_ports\":[{sides}],\"finite_chimneys\":[{chimneys}],\"side_opening_scope\":\"original DG face junctions and explicitly declared finite DG chimneys; branch outlets own exterior geometry and accepted radiated flow; complete mutual terminal radiation; no moving pad or three-dimensional junction correction\""));
    }
    output
}

fn dg_molecular_provenance(duct: &fs_duct::dg::DgDuct) -> String {
    if duct.molecular_relaxation_field().is_some() {
        let field = duct.gas_field().expect("spatial molecular owner retains its source field");
        let stations = field.stations_m().iter().map(|position| format!("{position:e}"))
            .collect::<Vec<_>>().join(",");
        return format!(",\"molecular_gas\":{{\"model\":\"shared-heat-capacity-reservoirs-v1\",\"medium\":\"original-axial-thermochemical-field\",\"stations_m\":[{stations}],\"scope\":\"same DG exponential; original local composition and temperature determine every molecular and radial thermal reservoir; integrated frozen pressure mass and power-preserving port coordinates; original Bessel error is not a combined-law certificate; no mean flow, axial conduction or live material transfer\"}}");
    }
    let Some(state) = duct.molecular_relaxation() else { return String::new(); };
    let reservoirs = molecular_reservoirs(state);
    format!(",\"molecular_gas\":{{\"model\":\"shared-heat-capacity-reservoirs-v1\",\"frozen_cv_j_kg_k\":{:e},\"frozen_sound_speed_m_s\":{:e},\"reservoirs\":[{reservoirs}],\"scope\":\"same DG exponential and original uniform equilibrium gas; coupled molecular and local radial thermal states at original nodal or spatial quadrature coordinates; original Bessel error is not a combined-law certificate; no axial conduction, additional bulk SK term or live chemistry\"}}",
        state.frozen_cv_j_kg_k(), state.frozen_sound_speed_m_s())
}

fn molecular_reservoirs(state: &fs_material::relaxation::RelaxationState) -> String {
    state.processes().iter().map(|process| format!(
        "{{\"species\":{},\"reservoir\":{},\"heat_capacity_j_kg_k\":{:e},\"relaxation_time_s\":{:e}}}",
        super::json_string(process.species.as_str()), process.index,
        process.heat_capacity_j_kg_k, process.relaxation_time_s,
    )).collect::<Vec<_>>().join(",")
}

// Both existing wind owners observe this same accepted DG terminal and record.
pub(super) fn dg_baffled_observation_provenance(
    load: &fs_couple::dg_port::DgPortLoad,
    receiver: fs_couple::pcm_wav::baffled::CircularOutletReceiver,
    microphone: &fs_couple::pcm_wav::baffled::BaffledPressure,
) -> String {
    dg_baffled_observation_on_clock(load, receiver, microphone, None)
}

pub(super) fn dg_baffled_projected_observation_provenance(
    load: &fs_couple::dg_port::DgPortLoad,
    receiver: fs_couple::pcm_wav::baffled::CircularOutletReceiver,
    microphone: &fs_couple::pcm_wav::baffled::BaffledPressure,
    observation_rate_hz: u32,
) -> String {
    dg_baffled_observation_on_clock(load, receiver, microphone, Some(observation_rate_hz))
}

fn dg_baffled_observation_on_clock(
    load: &fs_couple::dg_port::DgPortLoad,
    receiver: fs_couple::pcm_wav::baffled::CircularOutletReceiver,
    microphone: &fs_couple::pcm_wav::baffled::BaffledPressure,
    observation_rate_hz: Option<u32>,
) -> String {
    let source = load.radiation_record().expect("admitted bound DG radiation source");
    let exterior = gas_state(source.gas_record().state());
    let (clock, reconstruction, flow_scope) = if observation_rate_hz.is_some() {
        ("observation", "overlap-projected-terminal-flow-bandlimited-differentiating-Rayleigh",
            "interval-mean DG terminal flows projected onto the declared observation clock")
    } else {
        ("mechanical", "same-clock-bandlimited-differentiating-Rayleigh", "every accepted held DG terminal flow")
    };
    let clock_record = observation_rate_hz.map_or_else(String::new,
        |rate| format!(",\"sample_rate_hz\":{rate}"));
    if !load.duct().side_ports().is_empty() {
        let apertures = load.terminal_ports().iter().zip(source.apertures()).map(|(port, aperture)|
            format!("{{\"port\":{},\"center_m\":[{:e},{:e}],\"radius_m\":{:e}}}",
                port.0, aperture.center_m[0], aperture.center_m[1], aperture.radius_m)
        ).collect::<Vec<_>>().join(",");
        return format!(",\"observation_scope\":\"exterior Rayleigh pressure from {flow_scope} in original matrix order and bound source gas; uniform circular openings in one infinite rigid baffle; full mutual terminal feedback through any declared finite DG chimney; no mouth, free-field body or room radiation\",\"terminal_receiver\":{{\"position_m\":[{:e},{:e},{:e}],\"apertures\":[{apertures}],\"exterior_gas\":{exterior},\"radial_rings\":{},\"angular_points\":{},\"maximum_frequency_hz\":{:e},\"reconstruction\":\"{reconstruction}\",\"kernel_support_{clock}_samples\":[{},{}],\"filter_delay_{clock}_samples\":{:e}{clock_record},\"delay_compensated\":false,\"radiation_feedback_added\":true,\"delay_scope\":\"kernel support includes fractional travel-time interpolation; filter delay excludes physical propagation\"}}",
            receiver.position_m[0], receiver.position_m[1], receiver.position_m[2],
            receiver.radial_rings, receiver.angular_points, receiver.maximum_frequency_hz,
            microphone.delay_samples.0, microphone.delay_samples.1,
            microphone.filter_delay_samples());
    }
    let aperture = source.apertures().first().expect("admitted bound DG outlet");
    let flow_scope = if observation_rate_hz.is_none() { "each accepted held DG outlet flow" } else { flow_scope };
    format!(",\"observation_scope\":\"exterior Rayleigh pressure from {flow_scope} in the bound computed-source gas; one uniform circular opening in an infinite rigid baffle; no mouth, free-field body or room radiation\",\"outlet_receiver\":{{\"position_m\":[{:e},{:e},{:e}],\"radius_m\":{:e},\"exterior_gas\":{exterior},\"radial_rings\":{},\"angular_points\":{},\"maximum_frequency_hz\":{:e},\"reconstruction\":\"{reconstruction}\",\"kernel_support_{clock}_samples\":[{},{}],\"filter_delay_{clock}_samples\":{:e}{clock_record},\"delay_compensated\":false,\"radiation_feedback_added\":true,\"delay_scope\":\"kernel support includes fractional travel-time interpolation; filter delay excludes physical propagation\"}}",
        receiver.position_m[0], receiver.position_m[1], receiver.position_m[2],
        aperture.radius_m, receiver.radial_rings, receiver.angular_points,
        receiver.maximum_frequency_hz, microphone.delay_samples.0,
        microphone.delay_samples.1, microphone.filter_delay_samples())
}

fn impedance(load:SeriesImpedanceSpec)->String {
    let compliance=load.compliance_m3_pa.map_or_else(||"null".into(),|v|format!("{v:e}"));
    format!("\"resistance_pa_s_m3\":{:e},\"inertance_pa_s2_m3\":{:e},\"compliance_m3_pa\":{compliance}",
        load.resistance_pa_s_m3,load.inertance_pa_s2_m3)
}

pub(super) fn provenance(p:&PlateValvePerformance)->String {
    // Report the resolved default even on an unsegmented tube. Source identity
    // already binds the exact authored species, basis and fractions.
    let exterior = aperture_exterior_gas(p).map_or_else(String::new, |gas|
        format!(",\"exterior_gas\":{}", gas_state(gas)));
    let slit = p.renderer().system().aperture().viscous_slit().map_or_else(String::new, |law|
        format!(",\"slit_flow\":{{\"model\":\"spatial-laminar-bernoulli-v1\",\"length_m\":{:e},\"dynamic_viscosity_pa_s\":{:e},\"viscosity_source\":\"inlet-gas\",\"scope\":\"quasi-steady Newtonian slit with separated outflow; not resolved channel flow or pressure recovery\"}}",
            law.length_m(), law.dynamic_viscosity_pa_s()));
    let force_json = format!(",\"ambient_gas\":{}{exterior}{slit}{}", gas_state(p.ambient_gas()), force_provenance(p));
    if let CoupledAperture::Dg(runtime) = p.renderer().system() {
        let load = runtime.load();
        let duct = load.duct();
        let mesh = duct.mesh();
        let stepper = load.propagator();
        let band = duct.viscothermal_band_rad_s()
            .expect("authored DG selects admitted Bessel memory")
            .map(|w| w / core::f64::consts::TAU);
        let (wall_owner, series_poles, thermal_poles, quadrature, complex_error, loss_error) =
            if let Some(wall) = duct.viscothermal_model() {
                ("uniform-consistent-mass", wall.series_terms().len(), wall.shunt_terms().len(),
                    "null".to_owned(), wall.max_complex_error(), wall.max_real_loss_error())
            } else {
                let wall = duct.spatial_viscothermal_admission()
                    .expect("spatial authored DG retains its admitted wall rule");
                ("local-cylinder-positive-quadrature", wall.resolution.poles_per_field,
                    wall.resolution.poles_per_field, wall.resolution.quadrature_points.to_string(),
                    wall.max_complex_error, wall.max_real_loss_error)
            };
        let terminal = match duct.outlet() {
            fs_duct::dg::DgBoundary::Closed => "closed",
            fs_duct::dg::DgBoundary::IdealOpen => "ideal-open",
            fs_duct::dg::DgBoundary::Matched => "matched",
            fs_duct::dg::DgBoundary::Port(_) => {
                if load.radiation_record().is_some() { "computed-baffled" }
                else { "causal-impedance-port" }
            }
        };
        let field = duct.gas_field().map_or_else(String::new, dg_field_provenance);
        let radiation = dg_radiation_provenance(load);
        return format!("{force_json},\"dg_duct\":{{\"model\":\"exact-zoh-bessel-dg-v1\",\"operator_identity\":\"{}\",\"length_m\":{:e},\"inlet_area_m2\":{:e},\"outlet_area_m2\":{:e},\"inlet_gas\":{},\"outlet_gas\":{},\"elements\":{},\"order\":{},\"total_dofs\":{},\"wall_memory_dofs\":{},\"outlet\":\"{terminal}\",\"wall_owner\":\"{wall_owner}\",\"series_poles\":{series_poles},\"thermal_poles\":{thermal_poles},\"wall_quadrature_points\":{quadrature},\"band_hz\":[{:e},{:e}],\"checked_max_complex_relative_error\":{:e},\"checked_max_real_loss_relative_error\":{:e},\"exponential_tolerance\":{:e},\"exponential_certificate_residual\":{:e},\"mean_port_impedance_pa_s_m3\":{:e}{field}{radiation},\"scope\":\"original geometry and material-gas DG exponential with ExactZoh plate/Maxwell aperture; implicit same-step mean pressure and total flow; pressure observation declared separately; finite-band Bessel pole approximation; taper uses a local-cylinder Webster approximation; spatial and nonlinear time convergence separate; loss ledger inferred, not an independent heat integral; authored solid sections and spectra, not measured-material qualification\"}}",
            duct.operator_identity().to_hex(), duct.total_length(), duct.inlet_area_m2(),
            duct.outlet_area_m2(), gas_state(duct.inlet_gas()), gas_state(duct.outlet_gas()),
            mesh.elements, mesh.order, duct.total_dofs(),
            duct.total_dofs() - 2 * duct.field_dofs(), band[0], band[1], complex_error,
            loss_error, stepper.tolerance(), stepper.certificate_residual(),
            load.impedance_pa_s_m3());
    }
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
    let molecular = reed_molecular_provenance(p);
    let (point, scope, observation_scope) = if p.molecular_driving_point()
        .is_some_and(|source| !source.open_hole_sites().is_empty()) && p.outlet_receiver().is_some() {
        ("baffled-openings",
            "terminal and open-chimney Rayleigh pressure; same original molecular duct, finite chimneys and independent flanged exterior loads in feedback and observation; finite-history sampled approximation with parallel translated baffles and individually propagated mouth pressures; no mutual exterior loading, nonlinear vortex, jet, structural or room radiation; not measured calibration",
            "sum of terminal and all open or vented chimney pressures in the bound exterior gas; each mouth uses its own range and observer filter")
    } else if p.molecular_driving_point().is_some() && p.outlet_receiver().is_some() {
        (point,
            "terminal-only Rayleigh pressure; same original molecular duct and flanged exterior load in feedback and observation; finite-history sampled approximation; no tone-hole, jet, structural or room radiation; not measured calibration",
            observation_scope)
    } else { (point, scope, observation_scope) };
    format!(",\"pressure_point\":\"{point}\",\"model_scope\":\"{scope}\",\"observation_scope\":\"{observation_scope}\",\"ambient_gas\":{},\"segment_gases\":[{}]{receiver}{field}{molecular}",
        gas_state(p.ambient_gas()),regions.join(","))
}

fn reed_molecular_provenance(p: &fs_couple::render::schedule::reed::ReedPerformance) -> String {
    let Some(point) = p.molecular_driving_point() else { return String::new(); };
    let state = point.material().state();
    let reservoirs = molecular_reservoirs(state);
    let branches = point.geometry().segments.iter().any(|segment| matches!(segment, fs_duct::Segment::ToneHole { .. }));
    let scope = if branches {
        "original homogeneous gas reservoirs and declared wall diffusion in characteristic feedback, finite chimneys and terminal or mouth flows; shared converged FIR and Rayleigh owners; explicit independent flanged mouth loads; fixed vent-admittance mixture, not resolved pad geometry; observer propagation remains a separate declared law at each mouth range; no live chemistry, regional molecular field, mutual exterior radiation or measured source-data claim"
    } else {
        "original homogeneous gas reservoirs and declared wall diffusion in characteristic feedback and terminal flow; shared converged FIR and Rayleigh owners; observer propagation remains a separate declared law; no live chemistry, regional molecular field, compact tone-hole branch or measured source-data claim"
    };
    let openings = if point.open_hole_sites().is_empty() { String::new() } else {
        let sites = point.open_hole_sites().iter().enumerate().map(|(mouth, site)| format!(
            "{{\"mouth\":{mouth},\"axial_position_m\":{:e},\"hole_radius_m\":{:e},\"bore_radius_m\":{:e},\"chimney_height_m\":{:e},\"open_chimney_length_m\":{:e},\"open_fraction\":{:e}}}",
            site.axial_position_m, site.hole_radius_m, site.bore_radius_m,
            site.chimney_height_m, site.open_chimney_length_m, site.open_fraction,
        )).collect::<Vec<_>>().join(",");
        let observers = p.molecular_open_hole_observers().enumerate()
            .map(|(mouth, (info, receipt))| {
                use fs_couple::air_path::GasAbsorptionModel;
                let filter = receipt.filter;
                let law = match filter.absorption.law {
                    GasAbsorptionModel::Iso9613Air { .. } => "iso9613",
                    GasAbsorptionModel::MolecularRelaxation => "relaxation",
                    GasAbsorptionModel::StokesKirchhoff => "stokes-kirchhoff",
                    GasAbsorptionModel::Lossless => "none",
                };
                let optional = |value: Option<f64>| value.map_or_else(|| "null".to_owned(), |value| format!("{value:e}"));
                let r = info.receiver;
                format!("{{\"mouth\":{mouth},\"position_m\":[{:e},{:e},{:e}],\"radius_m\":{:e},\"gas\":{},\"radial_rings\":{},\"angular_points\":{},\"maximum_frequency_hz\":{:e},\"sample_rate_hz\":{},\"fft_samples\":{},\"axial_delay_mechanical_samples\":{},\"exterior_delay_mechanical_samples\":[{},{}],\"priming_source_sample\":-1,\"observer\":{{\"law\":\"{law}\",\"declared\":{},\"range_m\":{:e},\"range_spread_m\":{},\"taps\":{},\"latency_s\":{:e},\"low_frequency_delay_s\":{:e},\"alignment_samples\":{:e},\"padding_samples\":{},\"truncated_energy\":{:e},\"budget_energy\":{:e},\"receiver_delay_s\":{:e},\"delay_speed_m_s\":{},\"consumer_delay_s\":{},\"internal_part_unpropagated\":{},\"classical_validity_exceeded_from_hz\":{}}}}}",
                    r.position_m[0], r.position_m[1], r.position_m[2], info.radius_m,
                    gas_state(&info.gas), r.radial_rings, r.angular_points, r.maximum_frequency_hz,
                    info.sample_rate_hz, info.fft_samples, info.axial_delay_samples,
                    info.exterior_delay_samples.0, info.exterior_delay_samples.1,
                    p.observer_law_declared(), filter.range_m, optional(receipt.range_spread_m),
                    filter.taps, filter.latency_s, filter.low_frequency_delay_s,
                    filter.alignment_samples, filter.padding_samples, filter.truncated_energy,
                    filter.budget_energy, filter.receiver_delay_s, optional(filter.delay_speed_m_s),
                    optional(receipt.consumer_delay_s), receipt.internal_part_unpropagated,
                    optional(filter.classical_validity_exceeded_from_hz))
            }).collect::<Vec<_>>().join(",");
        format!(",\"open_mouths\":[{sites}],\"open_mouth_observers\":[{observers}],\"mouth_flow_scope\":\"original mouth volume flow includes the vent fraction once; full authored radius in both independent flanged load and Rayleigh disk; observer filters act before pressure summation\"")
    };
    format!(",\"molecular_gas\":{{\"model\":\"shared-heat-capacity-reservoirs-v1\",\"frozen_cv_j_kg_k\":{:e},\"frozen_sound_speed_m_s\":{:e},\"one_way_wavefront_s\":{:e},\"propagation_blocks\":{},\"material_points\":{},\"exterior_gas\":{},\"reservoirs\":[{reservoirs}]{openings},\"scope\":\"{scope}\"}}",
        state.frozen_cv_j_kg_k(), state.frozen_sound_speed_m_s(), point.one_way_wavefront_s(),
        point.material().blocks(), point.material().material_points(), gas_state(point.exterior_gas()))
}

// Serialize actual admitted states only; never re-resolve names or reconstruct
// air from density/speed. The same snapshot vocabulary covers ambient and regions.
/// The reed performance's observer receipt for its sidecar: the applied law
/// and its causal realization, the listener gas's whole figure, and what the
/// law leaves out (rulings 6510, 6532). An undeclared law that leaves out more
/// than the tolerance, or an unbounded or unrecorded amount, refuses here,
/// before any output exists.
pub(super) fn reed_observer_json(
    p: &fs_couple::render::schedule::reed::ReedPerformance,
) -> Result<String, fs_couple::acoustic_realize::AcousticRealizeError> {
    use fs_couple::air_path::GasAbsorptionModel;
    let omission = p.observer_admission()?;
    let receipt = p.observer();
    let filter = receipt.filter;
    let (law, omitted) = match filter.absorption.law {
        GasAbsorptionModel::Iso9613Air { .. } => ("iso9613", Some(0.0)),
        GasAbsorptionModel::MolecularRelaxation => ("relaxation", Some(0.0)),
        GasAbsorptionModel::StokesKirchhoff => ("stokes-kirchhoff", omission.omitted_db(true)),
        GasAbsorptionModel::Lossless => ("none", omission.omitted_db(false)),
    };
    let omitted = omitted.map_or_else(|| "null".to_owned(), |db| format!("{db:e}"));
    let observer = if receipt.internal_part_unpropagated { "proxy" } else { "baffled-outlet" };
    Ok(format!(
        ",\"observer\":{{\"observer\":\"{observer}\",\"law\":\"{law}\",\"declared\":{},\"omitted_db\":{omitted},\
         \"taps\":{},\"latency_s\":{:e},\"low_frequency_delay_s\":{:e},\"alignment_samples\":{:e},\
         \"padding_samples\":{},\"truncated_energy\":{:e},\"budget_energy\":{:e},\"receiver_delay_s\":{:e},\
         \"internal_part_unpropagated\":{},\"figure\":{}}}",
        p.observer_law_declared(),
        filter.taps,
        filter.latency_s,
        filter.low_frequency_delay_s,
        filter.alignment_samples,
        filter.padding_samples,
        filter.truncated_energy,
        filter.budget_energy,
        filter.receiver_delay_s,
        receipt.internal_part_unpropagated,
        omission.to_json(),
    ))
}

pub(super) fn gas_state(g: &fs_material::gas::GasState) -> String {
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

// Keep the original field primitives and source receipts with the DG observation.
pub(super) fn dg_field_provenance(field: &fs_material::gas_field::AxialGasField) -> String {
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
    format!(",\"axial_gas_field\":{{\"static_pressure_pa\":{:e},\"canonical_species\":[{}],\"stations\":[{}],\"source_packs\":[{}],\"scope\":\"frozen primitive temperature/composition; shared EOS and transport at original spatial quadrature; no mean-flow correction\"}}",
        field.static_pressure_pa(), species.join(","), stations.join(","), packs.join(","))
}

