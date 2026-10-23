//! G0/G3/G4: material keys share the actual two-sided valve/radiation transaction.
use super::*;
use fs_couple::bernoulli_aperture::dynamic::relaxation::{
    InitialApertureMemory, PlateRelaxationRegion, PlateRelaxationSpec,
};
use fs_couple::bernoulli_aperture::network::boundaries::{BoundaryAperture, BoundaryDrive};
use fs_couple::bernoulli_aperture::plate::{PlateApertureOptions, PlateApertureReduction};
use fs_couple::bernoulli_aperture::slit_channel::{SlitChannel, SlitChannelSpec, SlitExit};
use fs_couple::bernoulli_aperture::two_sided::supply::combined_inlet_state_with_supply;
use fs_couple::bernoulli_aperture::viscous::ViscousSlit;
use fs_material::visco::GeneralizedMaxwell;
use fs_plate::{AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateSection, SliceOptions};

fn plate(density: f64) -> PlateApertureReduction {
    let chart = PlateChart::with_boundary_and_regions(PlateMesh::rectangle(0.025, 0.01, 6, 2),
        PlateSection::isotropic(4e9, 0.3, 0.0003, density).unwrap(),
        vec![0, 7, 14], vec![]).unwrap();
    PlateApertureReduction::from_chart(chart, PlateApertureOptions {
        assembly: AssemblyOptions { pretension: 0.0, support: EdgeSupport::Clamped },
        eigenvalue_window: (1.0, (core::f64::consts::TAU * 450.0).powi(2)),
        mode_index: 0, eigensolver: SliceOptions::default(), max_nodes: 100, max_triangles: 100,
        slit_edges: vec![[6, 13], [13, 20]], rest_opening_m: 0.0002,
        damping_ratio: 0.0, max_slit_mode_variation: 0.25, max_slope: 0.1,
    }, &CancelGate::new_clock_free()).unwrap()
}

fn ducts() -> [ApertureDuctSpec; 2] {
    let mut ducts = branched_sides();
    for (duct, temperature) in ducts.iter_mut().zip([283.15, 313.15]) {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, gas().pressure).unwrap();
        duct.graph.nodes.push(NetworkNode::Driven);
        duct.graph.sections.push(TubeSection { nodes: [1, 4],
            length_m: 5.0 * gas.sound_speed * DT, radius_m: 0.003, max_length_error_m: 1e-12 });
        duct.gases.push(gas);
    }
    ducts
}

fn key(plate: &PlateApertureReduction, gas: GasState, z: f64, steps: u64) -> BoundaryAperture {
    let relaxation = PlateRelaxationSpec { regions: vec![PlateRelaxationRegion {
        triangles: (0..plate.chart().mesh.tris.len()).collect(),
        material: GeneralizedMaxwell::new(4e9, vec![(2e9, 0.001), (1e9, 0.01)]).unwrap(),
        poisson_ratio: 0.3, band_hz: (0.0, 1000.0),
        provenance: "synthetic material key, not measured woodwind data".into(),
    }], max_branches: 2, max_dt_over_tau: 0.1, max_angular_step: 0.1 };
    let contact = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e10, 2.0,
        "synthetic key contact".into()).unwrap().with_internal_loss(0.5).unwrap();
    let mut aperture = DynamicAperture::from_plate_with_propagator(plate.clone(), gas.density,
        z, DT, steps, ApertureState { opening_m: plate.options().rest_opening_m,
            opening_velocity_m_s: 0.0 }, contact, BladePropagator::ExactZoh).unwrap()
        .with_plate_relaxation(relaxation, InitialApertureMemory::Relaxed).unwrap();
    aperture.set_viscous_slit(ViscousSlit::from_gas(&gas, 0.0003).unwrap()).unwrap();
    BoundaryAperture { node: 4, aperture }
}

fn add_keys(mut system: TwoSidedAperture, plate: &PlateApertureReduction) -> TwoSidedAperture {
    for side in SIDES {
        let section = system.graph(side).sections.iter().position(|s| s.nodes.contains(&4)).unwrap();
        let gas = system.section_gases(side)[section];
        let z = system.represented_sections(side)[section].impedance_pa_s_m3;
        let key = key(plate, gas, z, system.valve().max_steps());
        system = system.with_boundary_apertures(side, vec![key], 1,
            &CancelGate::new_clock_free()).unwrap();
    }
    system
}

fn make(plate: &PlateApertureReduction, kind: Kind, steps: u64, budget: TwoSidedLimits)
    -> TwoSidedAperture
{
    let ducts = ducts();
    let mut valve = valve(kind, &ducts);
    if steps > valve.max_steps() {
        match &mut valve {
            LoadedAperture::Single(a) => a.extend_step_budget(steps),
            LoadedAperture::Paired(a) => a.extend_step_budget(steps),
            LoadedAperture::Modal(a) => a.extend_step_budget(steps),
        }.unwrap();
    }
    let system = TwoSidedAperture::new(valve, ducts, budget, &CancelGate::new_clock_free()).unwrap();
    add_keys(attach(system, 0.04), plate)
}

fn drives(force: f64, exterior: f64) -> [BoundaryDrive; 1] {
    [BoundaryDrive { node: 4, exterior_pressure_pa: exterior, opening_force_n: force }]
}

fn same_state(actual: &TwoSidedAperture, expected: &TwoSidedAperture) {
    assert_eq!(actual.valve().accepted_steps(), expected.valve().accepted_steps());
    assert_eq!(actual.stored_energy_j(), expected.stored_energy_j());
    for side in SIDES {
        assert_eq!(actual.mimo_pole_states(side), expected.mimo_pole_states(side));
        assert_eq!(actual.mimo_port_frames(side), expected.mimo_port_frames(side));
        for node in 0..actual.graph(side).nodes.len() {
            assert_eq!(actual.node_frame(side, node), expected.node_frame(side, node));
        }
        let a = actual.boundary_apertures(side).unwrap();
        let b = expected.boundary_apertures(side).unwrap();
        assert_eq!(a.last_frame(), b.last_frame());
        let (a, b) = (a.aperture(4).unwrap(), b.aperture(4).unwrap());
        assert_eq!(a.state(), b.state());
        assert_eq!(a.accepted_steps(), b.accepted_steps());
        assert_eq!(a.relaxation().unwrap().memory_sqrt_j(), b.relaxation().unwrap().memory_sqrt_j());
    }
    assert_eq!(actual.supply_channel().map(|s| s.last_step()),
        expected.supply_channel().map(|s| s.last_step()));
}

#[test]
fn material_keys_drive_all_three_valves_and_complete_external_work_once() {
    let light = plate(900.0);
    let dense = plate(1800.0);
    for kind in KINDS {
        let mut system = make(&light, kind, 1024, limits());
        let mut changed = make(&dense, kind, 1024, limits());
        let (mut mechanics, mut upstream, mut downstream) = (false, false, false);
        let mut balance = -system.stored_energy_j();
        let mut supplied_key_work = 0.0;
        let mut exterior_work = 0.0;
        for sample in 0..512 {
            let up = drives(if sample < 256 { -0.01 } else { 0.0 }, 0.3);
            let down = drives(if sample < 256 { -0.02 } else { 0.0 }, 0.7);
            let input = drive(if sample < 128 { 20.0 } else { 0.0 });
            let before = system.stored_energy_j();
            let frame = system.step_with_boundaries(input, [&up, &down]).unwrap();
            let other = changed.step_with_boundaries(input, [&up, &down]).unwrap();
            assert!(frame.balance_residual_j().abs() < 1e-12);
            assert!((frame.storage_change_j - (frame.stored_energy_j - before)).abs() < 1e-12);
            let mut boundary_work = 0.0;
            for (side, network, boundary) in SIDES.into_iter()
                .zip([frame.upstream, frame.downstream]).zip(frame.boundaries)
                .map(|((s, n), b)| (s, n, b.unwrap())) {
                assert_eq!(boundary.network, network);
                assert!(boundary.balance_residual_j().abs() < 1e-12);
                let keys = system.boundary_apertures(side).unwrap();
                let key = keys.aperture_frame(4).unwrap();
                assert_eq!(key.step, frame.aperture.step);
                assert_eq!(system.node_frame(side, 4).unwrap().pressure_pa, key.bore_pressure_pa);
                assert_eq!(keys.aperture(4).unwrap().viscous_slit().unwrap().gas(), keys.gas(4));
                assert_ne!(keys.gas(4), Some(&gas()));
                supplied_key_work += boundary.mechanical_work_j.abs();
                exterior_work += boundary.exterior_work_j.abs();
                boundary_work += boundary.mechanical_work_j + boundary.exterior_work_j;
            }
            balance += frame.dissipated_energy_j + frame.coupled_terminal_work_j
                - frame.source_work_j - frame.aperture.mechanical_work_j - boundary_work;
            mechanics |= frame.aperture.state != other.aperture.state;
            upstream |= frame.upstream.inlet_pressure_pa != other.upstream.inlet_pressure_pa;
            downstream |= frame.downstream.inlet_pressure_pa != other.downstream.inlet_pressure_pa;
        }
        assert!(mechanics && upstream && downstream && supplied_key_work > 0.0 && exterior_work > 0.0);
        assert!((system.stored_energy_j() + balance).abs() < 1e-10);
        for side in SIDES {
            assert!(system.boundary_apertures(side).unwrap().aperture(4).unwrap()
                .relaxation().unwrap().dissipated_energy_j() > 0.0);
        }
    }
}

#[test]
fn zero_breath_closure_and_reopening_retain_chimney_pressure_and_material_memory() {
    let plate = plate(900.0);
    let mut system = make(&plate, Kind::Paired, 6144, limits());
    let (mut closed, mut reopened, mut trapped, mut reed_motion) = (false, false, false, false);
    for sample in 0..6144 {
        let up = drives(0.0, 0.0);
        let down = drives(if sample < 2048 { -0.03 } else { 0.0 }, 0.0);
        if sample == 2048 {
            assert!(system.stored_energy_j() > 0.0);
            assert!(system.boundary_apertures(ValveSide::Downstream).unwrap().aperture(4).unwrap()
                .relaxation().unwrap().memory_sqrt_j().iter().any(|v| *v != 0.0));
        }
        let frame = system.step_with_boundaries(drive(0.0), [&up, &down]).unwrap();
        assert_eq!(frame.source_work_j, 0.0);
        assert!(frame.balance_residual_j().abs() < 1e-12);
        let key = system.boundary_apertures(ValveSide::Downstream).unwrap().aperture_frame(4).unwrap();
        if sample < 2048 && key.midpoint_opening_m == 0.0 {
            closed = true;
            assert_eq!(key.jet_flow_m3_s, 0.0);
            trapped |= key.bore_pressure_pa.abs() > 0.01;
        }
        reopened |= sample >= 2048 && key.midpoint_opening_m > 1e-4;
        reed_motion |= frame.aperture.state.opening_velocity_m_s != 0.0;
    }
    assert!(closed && reopened && trapped && reed_motion);
}

fn supplied(plate: &PlateApertureReduction, steps: u64) -> TwoSidedAperture {
    let mut ducts = ducts();
    for duct in &mut ducts {
        for (section, gas) in duct.graph.sections.iter_mut().zip(&duct.gases) {
            section.length_m += 0.25 * gas.sound_speed * DT;
        }
        duct.viscothermal = vec![ViscothermalSectionSpec { section: 2, cells: 2,
            minimum_frequency_hz: 200.0, maximum_frequency_hz: 400.0 }];
    }
    let channel = SlitChannel::from_gas(&gas(), SlitChannelSpec { length_m: 0.015,
        width_m: 0.01, height_m: 0.001, transverse_cells: 8, sample_period_s: DT,
        exit: SlitExit::StaticPressure }).unwrap();
    let (gas, z) = combined_inlet_state_with_supply(&ducts, &channel).unwrap();
    let policy = FractionalTubePolicy { maximum_frequency_hz: 1000.0,
        maximum_summed_phase_error_rad: 0.1 };
    let system = TwoSidedAperture::with_supply_channel(paired(z, gas, steps, 1.5).into(),
        ducts, channel, limits(), Some([policy; 2]), &CancelGate::new_clock_free()).unwrap();
    add_keys(attach(system, 0.04), plate)
}

#[test]
fn second_side_refusal_cancellation_and_horizon_extension_preserve_every_history() {
    let plate = plate(900.0);
    let mut actual = supplied(&plate, 64);
    let mut reference = supplied(&plate, 96);
    let up = drives(-0.005, 0.3);
    let down = drives(-0.01, 0.7);
    for _ in 0..32 {
        assert_eq!(actual.step_with_boundaries(drive(20.0), [&up, &down]).unwrap(),
            reference.step_with_boundaries(drive(20.0), [&up, &down]).unwrap());
    }
    for bad in [vec![], vec![BoundaryDrive { node: 3, ..down[0] }],
        vec![BoundaryDrive { opening_force_n: f64::MAX, ..down[0] }]] {
        assert!(actual.step_with_boundaries(drive(20.0), [&up, &bad]).is_err());
        same_state(&actual, &reference);
    }
    let gate = CancelGate::new_clock_free(); gate.request();
    assert!(actual.step_with_boundaries_under_gate(drive(20.0), [&up, &down], &gate).is_err());
    same_state(&actual, &reference);
    assert!(actual.extend_step_budget(64).is_err());
    assert!(actual.extend_step_budget(u64::MAX).is_err());
    same_state(&actual, &reference);
    actual.extend_step_budget(96).unwrap();
    for side in SIDES {
        assert_eq!(actual.boundary_apertures(side).unwrap().aperture(4).unwrap().spec().max_steps, 96);
        assert!(actual.fractional_sections(side).is_some());
        assert_eq!(actual.viscothermal_sections(side).len(), 1);
    }
    for _ in 32..96 {
        assert_eq!(actual.step_with_boundaries(drive(0.0), [&up, &down]).unwrap(),
            reference.step_with_boundaries(drive(0.0), [&up, &down]).unwrap());
    }
    same_state(&actual, &reference);
}

#[test]
fn aggregate_energy_refusal_cannot_publish_either_key_or_radiation_candidate() {
    let plate = plate(900.0);
    let mut reference = make(&plate, Kind::Modal, 1024, limits());
    let mut probe = make(&plate, Kind::Modal, 1024, limits());
    let quiet = drives(0.0, 0.0);
    let mut peak = 0.0_f64;
    for _ in 0..64 {
        let frame = reference.step_with_boundaries(drive(0.1), [&quiet, &quiet]).unwrap();
        assert_eq!(frame, probe.step_with_boundaries(drive(0.1), [&quiet, &quiet]).unwrap());
        peak = peak.max(frame.stored_energy_j);
    }
    let force = drives(-0.03, 1.0);
    let next = probe.step_with_boundaries(drive(100.0), [&force, &force]).unwrap();
    assert!(next.stored_energy_j > peak);
    let budget = TwoSidedLimits { maximum_energy_j: 0.5 * (peak + next.stored_energy_j), ..limits() };
    let mut limited = make(&plate, Kind::Modal, 1024, budget);
    for _ in 0..64 { limited.step_with_boundaries(drive(0.1), [&quiet, &quiet]).unwrap(); }
    let error = limited.step_with_boundaries(drive(100.0), [&force, &force]).unwrap_err();
    assert!(matches!(error, fs_couple::acoustic_realize::AcousticRealizeError::InvalidDescription {
        what: "two-sided valve refused the complete storage, work or flow balance",
    }));
    same_state(&limited, &reference);
    assert_eq!(limited.step_with_boundaries(drive(0.0), [&quiet, &quiet]).unwrap(),
        reference.step_with_boundaries(drive(0.0), [&quiet, &quiet]).unwrap());
}

#[path = "keyed/schedule.rs"]
mod schedule;
