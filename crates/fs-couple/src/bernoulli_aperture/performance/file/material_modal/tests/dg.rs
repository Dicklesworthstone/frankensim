//! The original regional material source retains every mode in the DG host.
use super::{pack, BladePropagator, CancelGate, LoadedAperture, ModalAperture,
    NormalizedMaterialCardPack, PlateValveInputError};
use crate::bernoulli_aperture::dg::loaded::performance::LoadedDgPerformance;
use crate::bernoulli_aperture::performance::file::material_modal::dg::MaterialModalDgValvePerformance;
use crate::pcm_wav::observation::PressureRenderer;

const INPUT: &str = include_str!("../../../../../../examples/material_modal_valve/paired-dg.performance");
const TENSOR: &str = "loss tensor eta11 eta22 eta66 eta12 eta16 eta26 1 1000000";

fn input() -> String {
    INPUT.replace("audio 96000 96000 4800 20", "audio 96000 96000 512 20")
}

fn build(text: &str, block: usize) -> MaterialModalDgValvePerformance {
    MaterialModalDgValvePerformance::from_bytes_with_sources(
        text.as_bytes(), block, &CancelGate::new_clock_free(), None,
        |path| Ok(pack(if path == "altered.matpack" { 1150.0 } else { 900.0 }, false)),
    ).unwrap()
}

fn valve(source: &MaterialModalDgValvePerformance) -> &ModalAperture {
    let LoadedAperture::Modal(valve) = source.renderer().system().aperture() else {
        panic!("the authored material source must retain its modal owner");
    };
    valve
}

fn render(source: MaterialModalDgValvePerformance) -> (Vec<f64>, LoadedDgPerformance) {
    let samples = source.info().samples as usize;
    let block = source.renderer().config().max_block;
    let mut renderer = source.into_renderer();
    let mut pressure = vec![0.0; samples];
    for chunk in pressure.chunks_mut(block) { renderer.block(chunk).unwrap(); }
    (pressure, renderer)
}

#[test]
fn original_material_modes_gas_and_dg_source_are_bound_and_block_invariant() {
    let text = input();
    let mut loads = 0;
    let source = MaterialModalDgValvePerformance::from_bytes_with_sources(
        text.as_bytes(), 37, &CancelGate::new_clock_free(), None,
        |path| { assert_eq!(path, "blade.matpack"); loads += 1; Ok(pack(900.0, false)) },
    ).unwrap();
    assert_eq!(loads, 1);
    assert_eq!((source.info().blades, source.info().modes, source.info().nodes,
        source.info().triangles), (2, 6, 42, 48));
    assert_eq!(source.material_records().len(), 4);
    let material = valve(&source);
    let load = source.renderer().system().load();
    assert_eq!(material.blade_propagator(), BladePropagator::ExactZoh);
    assert_eq!(material.gas(), Some(load.duct().inlet_gas()));
    assert_eq!(material.spec().impedance_pa_s_m3.to_bits(),
        load.impedance_pa_s_m3().to_bits());
    assert_eq!(material.viscosity_matrix_n_s_m().unwrap().len(), 36);
    assert_eq!(material.material_plates()[0].mode_range(), 0..3);
    assert_eq!(material.material_plates()[1].mode_range(), 3..6);
    assert_eq!(load.state().step_count(), 0);
    assert!(load.duct().has_viscothermal_memory());

    let (reference, rendered) = render(source);
    assert!(reference.iter().all(|p| p.is_finite()));
    assert!(reference.iter().any(|p| p.abs() > 1e-15));
    assert_eq!(rendered.samples_rendered(), 512);
    assert_eq!(*reference.last().unwrap(), rendered.last_frame().unwrap().aperture.bore_pressure_pa);
    for block in [1, 512] {
        let (other, source) = render(build(&text, block));
        assert_eq!(reference, other);
        assert_eq!(rendered.system().load().state(), source.system().load().state());
    }
    for changed in [
        text.replacen("section blade.matpack", "section altered.matpack", 1),
        text.replace("101325", "81060"),
    ] {
        let (other, _) = render(build(&changed, 37));
        assert!(reference.iter().zip(other).any(|(a, b)| (a - b).abs() > 1e-15));
    }
}

#[test]
fn original_spatial_pressure_force_ports_and_maxwell_memory_reach_the_same_dg_history() {
    let pressure_nodes = "pressure_node 1\n".repeat(21);
    let additions = format!(
        "pressure_shapes 1\npressure_shape 0 1000000 authored-dg-pressure\npressure_nodes 21\n{pressure_nodes}\
source_forces 2\nsource_force 3 1 -0.0001\nsource_force 97 1 0\n\
source_pressures 2\nsource_pressure 5 0 1\nsource_pressure 101 0 0\n"
    );
    let text = input()
        .replace(TENSOR, "loss regional-maxwell relaxed 1 boundary_modulus boundary_time")
        .replace("pressure-pa\t5e0", "pressure-pa\t0e0")
        .replace("compile_limits", &format!("{additions}compile_limits"));
    let source = build(&text, 37);
    assert!(valve(&source).material_relaxation().unwrap().is_regional());
    assert_eq!(valve(&source).material_pressure_ports()[0].as_ref().unwrap()
        .shape().provenance(), "authored-dg-pressure");
    let mut renderer = source.into_renderer();
    let mut pressure_work = 0.0;
    let mut mechanical_work = 0.0;
    for _ in 0..512 {
        renderer.block(&mut [0.0]).unwrap();
        let frame = renderer.last_frame().unwrap();
        pressure_work += frame.spatial_pressure[0].work_j.abs();
        mechanical_work += frame.aperture.mechanical_work_j.abs();
    }
    assert!(pressure_work > 0.0 && mechanical_work > 0.0);
    let LoadedAperture::Modal(material) = renderer.system().aperture() else {
        panic!("original modal source");
    };
    assert!(material.material_relaxation().unwrap().memory_sqrt_j()
        .iter().flatten().any(|memory| *memory != 0.0));
    assert_eq!(renderer.system().load().state().step_count(), 512);
    assert!(renderer.last_frame().unwrap().aperture.bore_pressure_pa.is_finite());
}

#[test]
fn invalid_clock_structure_and_unbound_observation_refuse_without_substitution() {
    let text = input();
    for early in [
        text.replace("audio 96000 96000", "audio 96000 48000"),
        text.replace("propagator exact-zoh", "propagator midpoint"),
        text.replace("propagator exact-zoh\n", ""),
    ] {
        assert!(MaterialModalDgValvePerformance::from_bytes_with_sources(
            early.as_bytes(), 37, &CancelGate::new_clock_free(), None,
            |_| -> Result<NormalizedMaterialCardPack, String> {
                panic!("clock and structural admission precede material loading")
            },
        ).is_err());
    }
    for refused in [
        text.replace("observation inlet-pressure",
            "observation baffled-outlet 0 0 0.25 8 32 1000"),
        text.replace("compile_limits",
            "source_forces 1\nsource_force 512 0 0.0001\ncompile_limits"),
        text.replace("compile_limits",
            "source_pressures 1\nsource_pressure 0 0 1\ncompile_limits"),
    ] {
        assert!(MaterialModalDgValvePerformance::from_bytes_with_sources(
            refused.as_bytes(), 37, &CancelGate::new_clock_free(), None,
            |_| Ok(pack(900.0, false)),
        ).is_err());
    }
    let error = MaterialModalDgValvePerformance::from_bytes_with_sources(
        text.as_bytes(), 37, &CancelGate::new_clock_free(), None,
        |_| Ok(pack(900.0, true)),
    ).err().unwrap();
    assert!(matches!(error, PlateValveInputError::Material { .. }));
}

#[test]
fn the_same_material_source_drives_reed_plates_and_the_original_dg_wall() {
    let declaration = "dg_wall isotropic blade.matpack 293.15 100 0.00002 stress-free 50 2000 maxwell 1 1 0.5 boundary_modulus boundary_time\n";
    let text = input().replace("slit_length", &format!("{declaration}slit_length"));
    let mut loads = 0;
    let source = MaterialModalDgValvePerformance::from_bytes_with_sources(
        text.as_bytes(), 37, &CancelGate::new_clock_free(), None, |path| {
            assert_eq!(path, "blade.matpack"); loads += 1; Ok(pack(900.0, false))
        },
    ).unwrap();
    assert_eq!(loads, 1, "reed and wall must share the original immutable artifact");
    assert_eq!(source.wall_record().unwrap().pack(), source.material_records()[0].pack());
    let duct = source.renderer().system().load().duct();
    let wall = duct.compliant_wall().unwrap();
    assert!((wall.surface_density - 0.018).abs() < 1e-15);
    assert!((wall.stiffness_per_area - 8e10).abs() < 1e-4);
    assert!(duct.has_viscothermal_memory());
    assert_eq!(duct.total_dofs(), 63, "gas, diffusion, radial motion and one physical Maxwell arm");
    let reference = render(source);
    assert_eq!(reference.0, render(build(&text, 1)).0);
    let load = reference.1.system().load();
    let motion = load.duct().compliant_wall_state(load.state()).unwrap();
    assert!(motion.stored_energy_j > 0.0);
    assert!(motion.relaxation_stored_energy_j > 0.0);
    assert_eq!(motion.relaxation_pressure_pa.len(), 1);
    let changed = text.replacen("isotropic blade.matpack", "isotropic altered.matpack", 1);
    let changed = build(&changed, 37);
    assert_eq!(changed.material_records()[0].pack(), &pack(900.0, false),
        "only the bore material changed; reed plates keep their original source");
    assert_ne!(reference.0, render(changed).0);
}

#[test]
fn physical_wall_regions_share_reed_artifacts_and_retain_local_material_memory() {
    let declaration = concat!(
        "dg_wall regions 2\n",
        "dg_wall_region 0 0.025 isotropic blade.matpack 293.15 100 0.00002 stress-free 50 2000 maxwell 1 1 0.5 boundary_modulus boundary_time\n",
        "dg_wall_region 0.025 0.05 isotropic blade.matpack 293.15 100 0.00001 stress-free 50 2000 maxwell 1 1 0.5 boundary_modulus boundary_time\n",
    );
    let text = input().replace("open 1 2 8", "open 2 2 8")
        .replace("slit_length", &format!("{declaration}slit_length"));
    let mut loads = 0;
    let source = MaterialModalDgValvePerformance::from_bytes_with_sources(
        text.as_bytes(), 37, &CancelGate::new_clock_free(), None, |path| {
            assert_eq!(path, "blade.matpack");
            loads += 1;
            Ok(pack(900.0, false))
        },
    ).unwrap();
    assert_eq!(loads, 1, "all reed and wall regions retain the one immutable source");
    let regions = source.wall_layout().regions();
    assert_eq!(regions.len(), 2);
    assert!(std::ptr::eq(regions[0].source().pack(), source.material_records()[0].pack()));
    let duct = source.renderer().system().load().duct();
    assert!(duct.has_viscothermal_memory());
    assert!(duct.compliant_wall().is_none());
    let first = duct.compliant_wall_at_element(0).unwrap();
    let second = duct.compliant_wall_at_element(1).unwrap();
    assert!((first.surface_density - 0.018).abs() < 1e-15);
    assert_eq!(second.surface_density, first.surface_density * 0.5);
    assert_eq!(second.stiffness_per_area, first.stiffness_per_area * 0.5);
    let reference = render(source);
    assert_eq!(reference.0, render(build(&text, 1)).0);
    let load = reference.1.system().load();
    let state = load.duct().compliant_wall_state(load.state()).unwrap();
    assert_eq!(state.relaxation_pressure_pa.len(), 2);
    assert!(state.relaxation_stored_energy_j > 0.0);
    let changed = text.replace("0.00001 stress-free", "0.00002 stress-free");
    assert_ne!(reference.0, render(build(&changed, 37)).0);
}

#[test]
fn finite_chimney_material_shares_reed_sources_and_changes_the_original_modal_feedback() {
    let declaration = concat!(
        "dg_wall assembly 2\n",
        "dg_wall_main rigid\n",
        "dg_wall_chimney 7 isotropic blade.matpack 293.15 100 0.00002 stress-free 50 2000 maxwell 1 1 0.5 boundary_modulus boundary_time\n",
    );
    let text = include_str!("../../../../../../examples/material_modal_valve/paired-dg-chimney.performance")
        .replace("audio 96000 96000 4800 20", "audio 96000 96000 512 20")
        .replace("slit_length", &format!("{declaration}slit_length"));
    let mut loads = 0;
    let source = MaterialModalDgValvePerformance::from_bytes_with_sources(
        text.as_bytes(), 37, &CancelGate::new_clock_free(), None, |path| {
            assert_eq!(path, "blade.matpack"); loads += 1; Ok(pack(900.0, false))
        },
    ).unwrap();
    assert_eq!(loads, 1);
    assert!(source.wall_record().is_none());
    let chimney = &source.wall_layout().chimneys()[0];
    assert_eq!(chimney.port(), fs_duct::dg::PortId(7));
    assert!(std::ptr::eq(chimney.layout().uniform_record().unwrap().pack(),
        source.material_records()[0].pack()));
    let load = source.renderer().system().load();
    assert!(load.duct().compliant_wall().is_none());
    let branch = &load.branches()[0];
    assert_ne!(branch.duct().inlet_gas(), load.duct().inlet_gas());
    assert!((branch.duct().compliant_wall().unwrap().surface_density - 0.018).abs() < 1e-15);
    let reference = render(source);
    assert_eq!(reference.0, render(build(&text, 1)).0);
    let branch = &reference.1.system().load().branches()[0];
    let state = branch.duct().compliant_wall_state(branch.state()).unwrap();
    assert!(state.stored_energy_j > 0.0 && state.relaxation_stored_energy_j > 0.0);
    let changed = text.replace("dg_wall_chimney 7 isotropic blade.matpack",
        "dg_wall_chimney 7 isotropic altered.matpack");
    let changed = build(&changed, 37);
    assert_eq!(changed.material_records()[0].pack(), &pack(900.0, false));
    let changed = render(changed);
    assert_ne!(reference.0, changed.0, "only the finite branch material changed");
    let states = |source: &LoadedDgPerformance| match source.system().aperture() {
        LoadedAperture::Modal(valve) => valve.modal_states().to_vec(),
        _ => panic!("the original material reed must remain active"),
    };
    assert_ne!(states(&reference.1), states(&changed.1));
}
