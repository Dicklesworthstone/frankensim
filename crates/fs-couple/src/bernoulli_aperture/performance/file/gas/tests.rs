use super::*;
use super::super::{PlateValvePerformance, AperturePerformance, CoupledAperture};
use crate::bernoulli_aperture::cavity::HelmholtzLoadSpec;
use crate::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture};
use crate::bernoulli_aperture::dynamic::relaxation::InitialApertureMemory;
use crate::bernoulli_aperture::network::{ApertureNetwork, NetworkNode, TubeNetworkSpec, TubeSection};
use crate::bernoulli_aperture::radiation::BaffledRadiationLoad;
use crate::bernoulli_aperture::tube::TubeDrive;
use crate::bernoulli_aperture::viscothermal::{ViscothermalSectionSpec, with_regional_viscothermal_sections};
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use crate::pcm_wav::observation::{DecimatedRenderer, PressureRenderer};
use crate::pcm_wav::stream::{Pcm16WavStream, ScheduledWavProgress, render_pressure_pcm16};
use fs_exec::CancelGate;
use std::io::Cursor;

const INPUT: &str = include_str!("../../../../../examples/plate-valve.performance");
const REGIONAL: &str = include_str!("../../../../../examples/plate-valve-regional-gas.performance");
const BASE: &str = "mixture 293.15 100000 mole 2 N2 0.7 O2 0.3";
const ENRICHED: &str = "mixture 293.15 100000 mole 3 N2 0.4 O2 0.2 CO2 0.4";

fn declaration(text: &str, inherited: Option<f64>) -> Result<GasState, PlateValveInputError> {
    let mut row = Row { fields: text.split_ascii_whitespace(), line: 7 };
    let gas = read(&mut row, inherited)?;
    row.finish()?;
    Ok(gas)
}
fn load(text: &str, block: usize) -> PlateValvePerformance {
    let resolved = text.replace("audio 96000 9602", "audio 192000 19204");
    PlateValvePerformance::from_bytes(resolved.as_bytes(), block, &CancelGate::new()).unwrap()
}
fn network(renderer: &AperturePerformance) -> &ApertureNetwork {
    let CoupledAperture::Network(n) = renderer.system() else { panic!("expected physical graph") };
    n
}
fn uniform(mixture: &str) -> String {
    REGIONAL.replace("ambient 293.15 101325 0", &format!("ambient {mixture}"))
        .replace(" gas 313.15 0.5", "").replace(" gas 283.15 0", "")
}
fn render(text: &str, block: usize) -> (Vec<f64>, ApertureState) {
    let p = load(text, block);
    let mut values = vec![0.0; p.info().samples as usize];
    let mut renderer = p.into_renderer();
    for out in values.chunks_mut(block) { renderer.block(out).unwrap(); }
    (values, renderer.system().aperture().state())
}
fn assert_changed(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    assert!(a.iter().chain(b).all(|p| p.is_finite()));
    let energy: f64 = a.iter().map(|p| p*p).sum();
    let difference: f64 = a.iter().zip(b).map(|(p,q)| (p-q)*(p-q)).sum();
    assert!(energy > 0.0 && difference > 1e-8*energy, "physical change must exceed floating-point noise");
}

#[test]
fn legacy_declarations_keep_exact_dry_and_humid_material_states() {
    for (t, p, rh) in [(293.15, 101325.0, 0.0), (313.15, 90000.0, 0.5), (283.15, 101325.0, 1.0)] {
        let expected = GasState::try_new_moist_air(t, p, rh).unwrap();
        assert_eq!(declaration(&format!("{t} {p} {rh}"), None).unwrap(), expected);
        assert_eq!(declaration(&format!("{t} {rh}"), Some(p)).unwrap(), expected);
    }
    let explicit = REGIONAL.replace("gas 313.15 0.5", "gas 293.15 0")
        .replace("gas 283.15 0", "gas 293.15 0");
    let implicit = REGIONAL.replace(" gas 313.15 0.5", "").replace(" gas 283.15 0", "");
    assert_eq!(render(&explicit, 37), render(&implicit, 113));
    assert_eq!(*load(INPUT, 37).ambient_gas(), GasState::try_new_moist_air(293.15, 101325.0, 0.0).unwrap());
}

#[test]
fn mass_basis_and_species_aliases_reach_the_same_material_owner() {
    let composition = Composition::new(CompositionBasis::MassFraction, vec![0.75, 0.20, 0.05]).unwrap();
    let expected = resolve_gas_mixture(&[GasSpecies::Nitrogen, GasSpecies::Oxygen, GasSpecies::Argon],
        &composition, 300.0, 100000.0).unwrap();
    let parsed = declaration("mixture 300 100000 mass 3 nitrogen 0.75 oxygen 0.20 argon 0.05", None).unwrap();
    assert_eq!(parsed, expected);
    let input = INPUT.replace("ambient 293.15 101325 0", "ambient mixture 300 100000 mass 3 N2 0.75 O2 0.20 Ar 0.05");
    let loaded = load(&input, 37);
    assert_eq!(*loaded.ambient_gas(), expected);
    assert!(matches!(loaded.renderer().system(), CoupledAperture::Tube(_)), "mixtures must not force a second propagation implementation");
    assert_eq!(loaded.renderer().system().aperture().spec().density_kg_m3, expected.density);
    assert_changed(&render(INPUT, 37).0, &render(&input, 37).0);
}

#[test]
fn regional_mixtures_match_independent_aperture_cavity_loss_radiation_and_receiver_construction() {
    let inlet = "mixture 300 101325 mole 2 N2 0.7 O2 0.3";
    let outlet = "mixture 303.15 101325 mole 3 N2 0.6 O2 0.2 Ar 0.2";
    let cavity = "mixture 293.15 101325 mole 2 N2 0.5 CO2 0.5";
    let text = REGIONAL.replace("gas 313.15 0.5", &format!("gas {outlet}"))
        .replace("gas 283.15 0", &format!("gas {cavity}"))
        .replace("duct_section 0 1 0.125 0.007 0.002 viscothermal",
            &format!("duct_section 0 1 0.125 0.007 0.002 gas {inlet} viscothermal"));
    // The reference states bypass the input adapter entirely.
    let make = |species: &[GasSpecies], x, t| resolve_gas_mixture(species,
        &Composition::new(CompositionBasis::MoleFraction, x).unwrap(), t, 101325.0).unwrap();
    let gases = vec![make(&[GasSpecies::Nitrogen, GasSpecies::Oxygen], vec![0.7,0.3], 300.0),
        make(&[GasSpecies::Nitrogen, GasSpecies::Oxygen, GasSpecies::Argon], vec![0.6,0.2,0.2], 303.15),
        make(&[GasSpecies::Nitrogen, GasSpecies::CarbonDioxide], vec![0.5,0.5], 293.15)];
    let dt = 1.0 / 192000.0;
    let gate = CancelGate::new();
    let radiation = BaffledRadiationLoad::new(0.007, gases[1].density, gases[1].sound_speed, dt, 1000.0, &gate).unwrap();
    let chamber = HelmholtzLoadSpec { volume_m3: 1e-5, neck_radius_m: 0.0015,
        effective_neck_length_m: 0.02, resistance_pa_s_m3: 200000.0 };
    let original = TubeNetworkSpec {
        nodes: vec![NetworkNode::Inlet, NetworkNode::Junction, radiation.termination(),
            chamber.termination(gases[2].density, gases[2].sound_speed).unwrap()],
        sections: vec![TubeSection { nodes: [0,1], length_m: 0.125, radius_m: 0.007, max_length_error_m: 0.002 },
            TubeSection { nodes: [1,2], length_m: 0.125, radius_m: 0.007, max_length_error_m: 0.002 },
            TubeSection { nodes: [1,3], length_m: 0.043, radius_m: 0.002, max_length_error_m: 0.002 }],
        sound_speed_m_s: gases[0].sound_speed, max_wave_memory_bytes: 1<<20,
    };
    let select = |section| ViscothermalSectionSpec { section, cells: 6, minimum_frequency_hz: 100.0, maximum_frequency_hz: 1000.0 };
    let (graph, mapped, reports) = with_regional_viscothermal_sections(original, gases.clone(), dt, &[select(0), select(1)], &gate).unwrap();
    let loaded = load(&text, 37);
    assert_eq!(network(loaded.renderer()).spec(), &graph);
    assert_eq!(network(loaded.renderer()).section_gases(), Some(mapped.as_slice()));
    assert_eq!(loaded.renderer().system().aperture().spec().density_kg_m3, gases[0].density);
    assert_ne!(loaded.ambient_gas().density, gases[0].density, "explicit inlet must override ambient in the Bernoulli law");
    for (i, report) in reports.iter().enumerate() {
        assert_eq!(loaded.viscothermal_losses()[i].source, report.source);
        assert_eq!(loaded.viscothermal_losses()[i].loss.dynamic_viscosity_pa_s(), gases[i].dynamic_viscosity);
        assert_eq!(loaded.viscothermal_losses()[i].loss.prandtl(), gases[i].prandtl);
        assert_eq!(report.one_way_samples, (0.125/(gases[i].sound_speed*dt)).round() as usize);
        for index in report.section_range[0]..report.section_range[1] { assert_eq!(mapped[index], gases[i]); }
    }
    let seed = load(INPUT, 37);
    let a = seed.renderer().system().aperture();
    let valve = DynamicAperture::from_plate_with_closure_and_propagator(a.plate_reduction().unwrap().clone(),
        a.plate_closure().unwrap().spec().clone(), gases[0].density,
        graph.inlet_impedance_with_gases(&mapped).unwrap(), dt, a.spec().max_steps, a.state(),
        a.blade_propagator()).unwrap()
        .with_plate_relaxation(a.relaxation().unwrap().spec().clone(), InitialApertureMemory::Relaxed).unwrap();
    let mut direct = ApertureNetwork::with_section_gases(valve, graph, mapped).unwrap();
    let mut mic = BaffledPressure::circular_outlet(0.007, 192000, CircularOutletReceiver {
        position_m: [0.0,0.0,0.2], radial_rings: 8, angular_points: 32, maximum_frequency_hz: 1000.0,
    }, RayleighMedium { density: gases[1].density, sound_speed: gases[1].sound_speed }).unwrap();
    let mut actual = loaded.into_renderer();
    let mut baseline = load(REGIONAL, 37).into_renderer();
    let mut changed_motion = false;
    for n in 0..1024 {
        let before = direct.stored_energy_j();
        let f = direct.step(TubeDrive { upstream_pressure_pa: seed.renderer().schedule().sample("mouth", n*700/192000).unwrap(), body_flow_m3_s: 0.0 }).unwrap();
        let scale = before + f.stored_energy_j + f.dissipated_energy_j + f.upstream_work_j.abs();
        assert!(f.balance_residual_j().abs() <= 3e-10*scale.max(1e-30));
        let expected = mic.step(&[direct.node_frame(2).unwrap().net_flow_into_node_m3_s]).unwrap();
        let mut out = [0.0]; actual.block(&mut out).unwrap();
        assert_eq!(out[0].to_bits(), expected.to_bits());
        assert_eq!(actual.system().aperture().state(), direct.aperture().state());
        baseline.block(&mut [0.0]).unwrap();
        changed_motion |= actual.system().aperture().state() != baseline.system().aperture().state();
    }
    assert!(changed_motion, "composition must reach reciprocal mechanical feedback, not only the observer");
}

#[test]
fn composition_and_pressure_change_physical_parameters_and_rendered_sound() {
    let base = uniform(BASE);
    let enriched = uniform(ENRICHED);
    let low_pressure = uniform(&BASE.replace("100000", "50000"));
    let a = load(&base, 37); let b = load(&enriched, 37); let c = load(&low_pressure, 37);
    let ga = a.ambient_gas(); let gb = b.ambient_gas(); let gc = c.ambient_gas();
    assert!(gb.sound_speed < ga.sound_speed && gb.density > ga.density);
    assert_ne!(gb.specific_heat_cp, ga.specific_heat_cp);
    assert_ne!(gb.dynamic_viscosity, ga.dynamic_viscosity);
    assert_ne!(gb.thermal_conductivity, ga.thermal_conductivity);
    assert_eq!(ga.sound_speed, gc.sound_speed, "pressure must not become an invented detuning control");
    assert_eq!(ga.density, 2.0*gc.density);
    assert_eq!(ga.characteristic_impedance, 2.0*gc.characteristic_impedance);
    assert_eq!(ga.dynamic_viscosity, gc.dynamic_viscosity);
    assert_ne!(network(a.renderer()).spec().nodes[2], network(c.renderer()).spec().nodes[2]);
    assert_ne!(network(a.renderer()).spec().nodes[3], network(c.renderer()).spec().nodes[3]);
    for (x,y) in network(a.renderer()).represented_sections().iter().zip(network(c.renderer()).represented_sections()) {
        assert_eq!(x.one_way_samples, y.one_way_samples);
    }
    assert_eq!(a.renderer().system().aperture().spec().mass_kg, b.renderer().system().aperture().spec().mass_kg);
    let (pa, _) = render(&base, 37); let (pb, sb) = render(&enriched, 37); let (pc, _) = render(&low_pressure, 37);
    assert_changed(&pa, &pb); assert_changed(&pa, &pc);
    assert_eq!((pb.clone(), sb), render(&enriched, 113));
    for (name, gas, pressure) in [("baseline", ga, &pa), ("CO2-enriched", gb, &pb), ("lower pressure", gc, &pc)] {
        let rms = (pressure.iter().map(|p| p*p).sum::<f64>()/pressure.len() as f64).sqrt();
        eprintln!("{name}: rho={} kg/m3 c={} m/s cp={} J/kg/K mu={} Pa s k={} W/m/K rendered_rms={} Pa",
            gas.density, gas.sound_speed, gas.specific_heat_cp, gas.dynamic_viscosity, gas.thermal_conductivity, rms);
    }
}

#[test]
fn mixture_pcm_cancel_resume_preserves_loaded_and_observer_history() {
    let text = uniform(ENRICHED);
    let build =
        || DecimatedRenderer::new(load(&text, 37).into_renderer(), 192000, 48000, 37).unwrap();
    let stream = || Pcm16WavStream::new(Cursor::new(Vec::new()), 48000, 0.01, 37).unwrap();
    let mut direct = build(); let mut expected = stream(); let mut scratch = [0.0;37];
    render_pressure_pcm16(&mut direct, &mut expected, &CancelGate::new(), &mut scratch, 4801).unwrap();
    let mut resumed = build(); let mut actual = stream();
    render_pressure_pcm16(&mut resumed, &mut actual, &CancelGate::new(), &mut scratch, 1001).unwrap();
    let before = network(resumed.source()).stored_energy_j().to_bits();
    let state = resumed.source().system().aperture().state();
    let gate = CancelGate::new(); gate.request(); scratch.fill(-123.0);
    assert_eq!(render_pressure_pcm16(&mut resumed, &mut actual, &gate, &mut scratch, 3800).unwrap(), ScheduledWavProgress::Cancelled { samples: 0 });
    assert_eq!(scratch, [-123.0;37]);
    assert_eq!(before, network(resumed.source()).stored_energy_j().to_bits());
    assert_eq!(state, resumed.source().system().aperture().state());
    render_pressure_pcm16(&mut resumed, &mut actual, &CancelGate::new(), &mut scratch, 3800).unwrap();
    let bytes = actual.finish().unwrap().0.into_inner();
    assert_eq!(&bytes[..4], b"RIFF");
    assert_eq!(bytes, expected.finish().unwrap().0.into_inner());
}

#[test]
fn invalid_mixtures_and_static_pressure_jumps_never_fall_back() {
    for input in ["mixture", "mixture 300 100000 mole 0", "mixture 300 100000 mole 7",
        "mixture 300 100000 volume 1 N2 1", "mixture 300 100000 mole 1 Xe 1",
        "mixture 300 100000 mole 2 N2 0.5 N2 0.5", "mixture 300 100000 mole 2 N2 1 O2 0",
        "mixture 300 100000 mole 2 N2 0.7 O2 0.2", "mixture 300 100000 mole 2 N2 1.1 O2 -0.1",
        "mixture 300 100000 mole 1 N2 NaN", "mixture 300 100000 mole 2 N2 1",
        "mixture 300 100000 mole 1 N2 1 O2 0.2", "mixture 400 100000 mole 1 N2 1",
        "mixture 300 300000 mole 1 N2 1", "mixture 293.15 100000 mole 2 N2 0.9 H2O 0.1"] {
        let source = INPUT.replace("ambient 293.15 101325 0", &format!("ambient {input}"));
        let err = PlateValvePerformance::from_bytes(source.as_bytes(), 37, &CancelGate::new()).err().expect("invalid mixture accepted");
        assert!(err.to_string().contains("line 3"), "lost source line: {err}");
    }
    let text = REGIONAL.replace("gas 313.15 0.5", "gas mixture 300 90000 mole 1 N2 1");
    let err = PlateValvePerformance::from_bytes(text.as_bytes(), 37, &CancelGate::new()).err().unwrap();
    assert!(err.to_string().contains("static pressure"));
    // Even all-explicit, mutually equal section pressures cannot ignore ambient.
    let text = uniform(BASE).replace("0.002 viscothermal", "0.002 gas mixture 300 90000 mole 1 N2 1 viscothermal")
        .replace("0.043 0.002 0.002\n", "0.043 0.002 0.002 gas mixture 300 90000 mole 1 N2 1\n");
    assert!(PlateValvePerformance::from_bytes(text.as_bytes(), 37, &CancelGate::new()).is_err());
    let error = declaration("mixture 293.15 100000 mole 1 H2O 1", None).unwrap_err();
    assert!(error.to_string().contains("saturation"));
    assert!(std::error::Error::source(&error).is_some());
}

#[test]
fn valid_pure_vapor_is_not_rejected_as_an_air_only_network() {
    let steam = declaration("mixture 323.15 10000 mole 1 H2O 1", None).unwrap();
    let graph = TubeNetworkSpec { nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection: 0.0 }],
        sections: vec![TubeSection { nodes: [0,1], length_m: 0.3, radius_m: 0.006, max_length_error_m: 0.005 }],
        sound_speed_m_s: steam.sound_speed, max_wave_memory_bytes: 100000 };
    graph.realize_waveguide(&[steam], 1.0/48000.0, &CancelGate::new()).unwrap();
    let mut invalid = steam; invalid.water_mole_fraction = 1.01;
    assert!(graph.validate_section_gases(&[invalid]).is_err());
}
