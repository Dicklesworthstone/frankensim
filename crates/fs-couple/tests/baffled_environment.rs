//! Exterior material data determines the same physical load and pressure medium.
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::radiation::{BaffledRadiationLoad,
    environment::bind_baffled_environment};
use fs_couple::pcm_wav::baffled::CircularOutletReceiver;
use fs_couple::pcm_wav::baffled::openings::duct::{
    DuctOpening, DuctRadiationSpec, DuctReconstruction,
};
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};

const RATE: u32 = 48_000;
fn gas(t: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), t, 101_325.0).unwrap()
}
fn graph(g: GasState, radius: f64) -> TubeNetworkSpec {
    TubeNetworkSpec { nodes: vec![NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Termination { reflection: -0.7 }],
        sections: vec![
            TubeSection { nodes: [0,1], length_m: 0.12, radius_m: radius, max_length_error_m: 0.005 },
            TubeSection { nodes: [1,2], length_m: 0.18, radius_m: radius, max_length_error_m: 0.005 }],
        sound_speed_m_s: g.sound_speed, max_wave_memory_bytes: 1<<20 }
}
fn radiation(exterior: Option<GasState>, count: usize) -> DuctRadiationSpec {
    DuctRadiationSpec {
        openings: (0..count).map(|network| DuctOpening { network, node: 2,
            center_m: [0.5*network as f64, 0.0] }).collect(), exterior,
        receiver: CircularOutletReceiver { position_m: [0.0,0.0,0.5], radial_rings: 16,
            angular_points: 64, maximum_frequency_hz: 400.0 },
        reconstruction: DuctReconstruction::Differentiating {
            mechanics_rate_hz: RATE, output_rate_hz: RATE/2 },
    }
}

#[test]
fn actual_radius_and_one_exterior_bind_identical_owner_loads_without_changing_the_bore() {
    let local = [gas(293.15), gas(315.15)]; let exterior = gas(333.15);
    let mut graphs = [graph(local[0],0.007), graph(local[1],0.005)];
    let before = graphs.clone(); let gases = [[local[0];2], [local[1];2]];
    let input = radiation(Some(exterior),2);
    let bound = bind_baffled_environment(&mut graphs,&[&gases[0],&gases[1]],input.clone(),
        &CancelGate::new_clock_free()).unwrap();
    assert_eq!(bound.radiation(),&input);
    for (i,binding) in bound.bindings().iter().enumerate() {
        let expected = BaffledRadiationLoad::new(before[i].sections[1].radius_m,
            exterior.density,exterior.sound_speed,1.0/f64::from(RATE),400.0,
            &CancelGate::new_clock_free()).unwrap();
        assert_eq!(binding.load(),&expected);
        assert_eq!(graphs[i].nodes[2],expected.termination());
        assert_eq!(binding.previous_reflection(),-0.7);
        assert_eq!(*binding.interior_gas(),local[i]);
        assert_eq!(binding.opening(),input.openings[i]);
        assert_eq!(graphs[i].nodes[..2],before[i].nodes[..2]);
        assert_eq!(format!("{:?}",graphs[i].sections),format!("{:?}",before[i].sections));
        assert_eq!(graphs[i].sound_speed_m_s,before[i].sound_speed_m_s);
        assert_eq!(graphs[i].max_wave_memory_bytes,before[i].max_wave_memory_bytes);
    }
}

#[test]
fn implicit_single_port_exterior_is_resolved_from_its_local_gas_not_the_inlet() {
    let inlet = gas(293.15); let outlet = gas(315.15);
    let mut graphs = [graph(inlet,0.007)];
    let bound = bind_baffled_environment(&mut graphs,&[&[inlet,outlet]],radiation(None,1),
        &CancelGate::new_clock_free()).unwrap();
    assert_eq!(bound.radiation().exterior,Some(outlet));
    assert_eq!(bound.bindings()[0].load().density_kg_m3(),outlet.density);
    assert_eq!(bound.bindings()[0].load().sound_speed_m_s(),outlet.sound_speed);
}

#[test]
fn later_boundary_failure_reactive_replacement_and_cancel_leave_all_graphs_unchanged() {
    let g=gas(293.15); let gases=[[g;2];2];
    let mut graphs=[graph(g,0.007),graph(g,0.08)];
    let before=format!("{graphs:?}");
    // The receiver spatial rule admits both disks, but the second load violates
    // ka<=0.5. The first successful load must not already have changed its node.
    assert!(bind_baffled_environment(&mut graphs,&[&gases[0],&gases[1]],radiation(Some(g),2),
        &CancelGate::new_clock_free()).is_err());
    assert_eq!(format!("{graphs:?}"),before);
    graphs[1]=graph(g,0.005);
    let existing=BaffledRadiationLoad::new(0.005,g.density,g.sound_speed,
        1.0/f64::from(RATE),400.0,&CancelGate::new_clock_free()).unwrap();
    graphs[1].nodes[2]=existing.termination(); let before=format!("{graphs:?}");
    assert!(bind_baffled_environment(&mut graphs,&[&gases[0],&gases[1]],radiation(Some(g),2),
        &CancelGate::new_clock_free()).is_err());
    assert_eq!(format!("{graphs:?}"),before);
    let gate=CancelGate::new_clock_free();gate.request();
    assert!(bind_baffled_environment(&mut graphs,&[&gases[0],&gases[1]],radiation(Some(g),2),&gate).is_err());
    assert_eq!(format!("{graphs:?}"),before);
}

#[test]
fn absent_mappings_duplicate_ports_and_pressure_mismatch_do_not_become_a_boundary_model() {
    let g=gas(293.15);let gases=[g;2];let mut graphs=[graph(g,0.007)];
    let original=format!("{graphs:?}");
    let gate=CancelGate::new_clock_free();
    assert!(bind_baffled_environment(&mut graphs,&[],radiation(None,1),&gate).is_err());
    let mut duplicate=radiation(Some(g),1);duplicate.openings.push(duplicate.openings[0]);
    assert!(bind_baffled_environment(&mut graphs,&[&gases],duplicate,&gate).is_err());
    let wrong=GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,90000.0).unwrap();
    assert!(bind_baffled_environment(&mut graphs,&[&gases],radiation(Some(wrong),1),&gate).is_err());
    assert_eq!(format!("{graphs:?}"),original);
}
