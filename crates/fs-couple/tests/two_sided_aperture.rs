//! Tests exercise the actual valve/network owners, not a stand-in output filter.
use fs_couple::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture, DynamicApertureSpec};
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::paired::{ApertureBlade, BladePropagator, BladeState,
    PairedAperture, PairedApertureLimits, PairedApertureSpec};
use fs_couple::bernoulli_aperture::two_sided::{ApertureDuctSpec, LoadedAperture,
    TwoSidedAperture, TwoSidedDrive, TwoSidedFrame, TwoSidedLimits, ValveSide, combined_inlet_state};
use fs_couple::bernoulli_aperture::{BernoulliAperture, dynamic::ApertureTerminal};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};

const DT: f64 = 1.0 / 48_000.0;
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn duct(g: GasState, delay: usize, reflection: f64) -> ApertureDuctSpec {
    ApertureDuctSpec {
        graph: TubeNetworkSpec {
            nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection }],
            sections: vec![TubeSection { nodes: [0, 1], length_m: delay as f64 * g.sound_speed * DT,
                radius_m: 0.007, max_length_error_m: 1e-12 }],
            sound_speed_m_s: g.sound_speed, max_wave_memory_bytes: 1 << 20,
        }, gases: vec![g], viscothermal: vec![],
    }
}
fn contact() -> Obstacle {
    Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
        "numerical two-sided valve contact fixture".into()).unwrap()
}
fn limits() -> TwoSidedLimits {
    TwoSidedLimits { maximum_energy_j: 1.0,
        energy_absolute_tolerance_j: 1e-16, energy_relative_tolerance: 1e-7,
        flow_absolute_tolerance_m3_s: 1e-16, flow_relative_tolerance: 1e-8 }
}
fn blade() -> ApertureBlade {
    ApertureBlade { mass_kg: 1e-5, stiffness_n_m: 500.0,
        damping_n_s_m: 0.03, pressure_area_m2: 1.25e-4 }
}
fn paired(z: f64, g: GasState, steps: u64, second_mass_ratio: f64) -> PairedAperture {
    let mut b = [blade(); 2];
    b[1].mass_kg *= second_mass_ratio;
    // Outside the midpoint one-cent band at 48 kHz (DEF-04): exact blade blocks.
    PairedAperture::new_with_propagator(PairedApertureSpec {
        rest_opening_m: 4e-4, width_m: 0.013, blades: b,
        density_kg_m3: g.density, impedance_pa_s_m3: z, time_step_s: DT, max_steps: steps,
        limits: PairedApertureLimits {
            maximum_abs_displacement_m: 0.01, maximum_abs_velocity_m_s: 100.0,
            maximum_energy_j: 1.0, max_evaluations: 256,
            velocity_tolerance_m_s: 1e-13, relative_tolerance: 1e-11,
        },
    }, [BladeState::default(); 2], contact(), BladePropagator::ExactZoh).unwrap()
}
fn system(up_delay: usize, mass_ratio: f64, steps: u64) -> TwoSidedAperture {
    let sides = [duct(gas(), up_delay, 0.7), duct(gas(), 32, -0.8)];
    let (g, z) = combined_inlet_state(&sides).unwrap();
    TwoSidedAperture::new(paired(z, g, steps, mass_ratio).into(), sides,
        limits(), &CancelGate::new_clock_free()).unwrap()
}
fn drive(p: f64) -> TwoSidedDrive {
    TwoSidedDrive { source_pressure_pa: p, blade_forces_n: [0.0; 2] }
}

#[test]
fn physical_faces_share_the_pressure_drop_flow_and_complete_work_balance() {
    let mut m = system(12, 1.7, 512);
    let gate = CancelGate::new_clock_free();
    let mut energy = 0.0;
    let mut total_source_work = 0.0;
    let mut total_loss = 0.0;
    let mut reflected = false;
    for k in 0..512 {
        let f = m.step_under_gate(drive(if k < 128 { 100.0 } else { 0.0 }), &gate).unwrap();
        let q = f.aperture.bore_flow_m3_s;
        assert!((f.upstream.inlet_flow_m3_s + q).abs() < 1e-14);
        assert!((f.downstream.inlet_flow_m3_s - q).abs() < 1e-14);
        assert_eq!(f.aperture.bore_pressure_pa, f.downstream.inlet_pressure_pa);
        assert_eq!(f.pressure_drop_pa, f.upstream_pressure_pa - f.aperture.bore_pressure_pa);
        assert!((f.aperture.pressure_work_j - DT*f.pressure_drop_pa*q).abs() < 1e-13);
        assert!((f.stored_energy_j - energy - f.storage_change_j).abs() < 1e-13);
        assert!(f.balance_residual_j().abs() < 1e-12);
        reflected |= f.upstream.incoming_pressure_pa != 0.0;
        total_source_work += f.source_work_j;
        total_loss += f.dissipated_energy_j;
        energy = f.stored_energy_j;
    }
    assert!(reflected, "the upstream network must actually react back");
    assert!(energy > 0.0, "source shutoff must not erase stored vibration");
    assert!((energy + total_loss - total_source_work).abs() < 1e-11);
}

#[test]
fn upstream_geometry_changes_the_reed_only_after_the_physical_return_time() {
    let mut short = system(12, 1.0, 128);
    let mut long = system(24, 1.0, 128);
    let gate = CancelGate::new_clock_free();
    let mut changed = false;
    for k in 0..128 {
        let a = short.step_under_gate(drive(100.0), &gate).unwrap();
        let b = long.step_under_gate(drive(100.0), &gate).unwrap();
        if k < 24 {
            assert_eq!(a.aperture.state, b.aperture.state);
            assert_eq!(a.aperture.bore_flow_m3_s, b.aperture.bore_flow_m3_s);
        } else {
            changed |= (a.aperture.state.opening_m - b.aperture.state.opening_m).abs() > 1e-12;
        }
    }
    assert!(changed, "upstream geometry must change mechanics, not just playback");
}

#[test]
fn both_original_mechanical_images_see_the_same_symmetric_load() {
    let sides = [duct(gas(), 12, 0.7), duct(gas(), 32, -0.8)];
    let (g, z) = combined_inlet_state(&sides).unwrap();
    let b = blade();
    // The symmetric image of the exact blades is the exact opening coordinate.
    let single = DynamicAperture::new_with_propagator(DynamicApertureSpec {
        aperture: BernoulliAperture { rest_opening_m: 4e-4, width_m: 0.013,
            closing_pressure_pa: (0.5*b.stiffness_n_m)*4e-4/b.pressure_area_m2 },
        mass_kg: 0.5*b.mass_kg, stiffness_n_m: 0.5*b.stiffness_n_m,
        damping_ratio: b.damping_n_s_m/(2.0*(b.mass_kg*b.stiffness_n_m).sqrt()),
        density_kg_m3: g.density, impedance_pa_s_m3: z, time_step_s: DT, max_steps: 128,
    }, ApertureState { opening_m: 4e-4, opening_velocity_m_s: 0.0 }, contact(),
        BladePropagator::ExactZoh).unwrap();
    let gate = CancelGate::new_clock_free();
    let mut a = TwoSidedAperture::new(single.into(), sides, limits(), &gate).unwrap();
    let mut b = system(12, 1.0, 128);
    for _ in 0..128 {
        let fa = a.step_under_gate(drive(20.0), &gate).unwrap();
        let fb = b.step_under_gate(drive(20.0), &gate).unwrap();
        assert!((fa.aperture.state.opening_m-fb.aperture.state.opening_m).abs() < 1e-9);
        assert!((fa.aperture.bore_pressure_pa-fb.aperture.bore_pressure_pa).abs() < 1e-3);
    }
}

#[test]
fn invalid_samples_and_cancellation_preserve_both_networks_and_both_blades() {
    let mut m = system(12, 1.5, 8);
    let gate = CancelGate::new_clock_free();
    m.step_under_gate(drive(50.0), &gate).unwrap();
    let stored = m.stored_energy_j();
    let up = *m.node_frame(ValveSide::Upstream, 0).unwrap();
    let down = *m.node_frame(ValveSide::Downstream, 0).unwrap();
    let states = match m.valve() { LoadedAperture::Paired(a) => *a.blade_states(), _ => unreachable!() };
    for bad in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(m.step_under_gate(drive(bad), &gate).is_err());
        assert_eq!(m.stored_energy_j(), stored);
        assert_eq!(*m.node_frame(ValveSide::Upstream, 0).unwrap(), up);
        assert_eq!(*m.node_frame(ValveSide::Downstream, 0).unwrap(), down);
        assert_eq!(m.valve().accepted_steps(), 1);
        match m.valve() { LoadedAperture::Paired(a) => assert_eq!(*a.blade_states(), states), _ => unreachable!() }
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(m.step_under_gate(drive(50.0), &cancelled).is_err());
    assert_eq!(m.valve().accepted_steps(), 1);
    m.step_under_gate(drive(50.0), &gate).unwrap();
    assert_eq!(m.valve().accepted_steps(), 2);
}

#[test]
fn lifetime_budget_and_block_partition_do_not_reset_any_participant() {
    let gate = CancelGate::new_clock_free();
    let mut blocks = system(12, 1.0, 3);
    let mut scalar = system(12, 1.0, 8);
    let drives = [drive(50.0); 8];
    let sentinel = TwoSidedFrame { upstream_pressure_pa: -999.0, ..TwoSidedFrame::default() };
    let mut out = [sentinel; 8];
    let progress = blocks.advance_block(&drives, &mut out, &gate).unwrap();
    assert_eq!(progress.completed, 3);
    assert_eq!(progress.terminal, ApertureTerminal::BudgetExhausted);
    assert!(out[3..].iter().all(|x| *x == sentinel));
    blocks.extend_step_budget(8).unwrap();
    blocks.advance_block(&drives[3..], &mut out[3..], &gate).unwrap();
    for frame in out { assert_eq!(frame, scalar.step_under_gate(drive(50.0), &gate).unwrap()); }
}

#[test]
fn material_mass_changes_both_loads_and_the_complete_coupled_trajectory() {
    let mut a = system(12, 1.0, 96);
    let mut b = system(12, 2.5, 96);
    let gate = CancelGate::new_clock_free();
    let mut mechanical = false;
    let mut upstream = false;
    let mut downstream = false;
    for _ in 0..96 {
        let fa = a.step_under_gate(drive(70.0), &gate).unwrap();
        let fb = b.step_under_gate(drive(70.0), &gate).unwrap();
        mechanical |= fa.aperture.state != fb.aperture.state;
        upstream |= fa.upstream.inlet_pressure_pa != fb.upstream.inlet_pressure_pa;
        downstream |= fa.downstream.inlet_pressure_pa != fb.downstream.inlet_pressure_pa;
    }
    assert!(mechanical && upstream && downstream);
}

#[test]
fn mismatched_inlet_gases_and_a_one_sided_valve_load_are_rejected() {
    let gate = CancelGate::new_clock_free();
    let mut sides = [duct(gas(), 12, 0.7), duct(gas(), 32, -0.8)];
    let (g, z) = combined_inlet_state(&sides).unwrap();
    let one_sided_z = sides[1].graph.inlet_impedance_with_gases(&sides[1].gases).unwrap();
    assert!(TwoSidedAperture::new(paired(one_sided_z,g,8,1.0).into(), sides.clone(), limits(), &gate).is_err());
    let warm = GasState::try_new(&GasSpec::dry_air_ussa1976(), 303.15, 101_325.0).unwrap();
    sides[0] = duct(warm, 12, 0.7);
    assert!(combined_inlet_state(&sides).is_err());
    assert!(TwoSidedAperture::new(paired(z,g,8,1.0).into(), sides, limits(), &gate).is_err());
}

#[path = "two_sided_aperture/mimo.rs"]
mod mimo;

