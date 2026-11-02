//! Material single/independent reeds share fixed-channel momentum and both ducts.
use super::*;
use fs_couple::bernoulli_aperture::{
    dynamic::{ApertureState, DynamicAperture},
    network::{NetworkNode, TubeNetworkSpec, TubeSection},
    paired::{BladeState, PairedAperture, PairedApertureLimits},
    slit_channel::{SlitChannel, SlitChannelSpec, SlitExit},
    two_sided::{ApertureDuctSpec, LoadedAperture, TwoSidedAperture, TwoSidedDrive,
        TwoSidedLimits, ValveSide, combined_inlet_state, observer::TwoSidedPressure,
        supply::combined_inlet_state_with_supply},
    viscous::ViscousSlit,
};
use fs_couple::pcm_wav::{baffled::CircularOutletReceiver, observation::PressureRenderer};
use fs_dcontact::Obstacle;
use fs_material::gas::{GasSpec, GasState};

const RATE: u32 = 192_000;
const DT: f64 = 1.0 / RATE as f64;
fn gas(viscosity: f64) -> GasState {
    let mut spec = GasSpec::dry_air_ussa1976();
    // Transport-only intervention, explicitly synthetic rather than a named gas.
    spec.sutherland_beta *= viscosity;
    GasState::try_new(&spec, 293.15, 101_325.0).unwrap()
}
fn ducts(g: GasState) -> [ApertureDuctSpec; 2] {
    [48, 96].map(|delay| ApertureDuctSpec { graph: TubeNetworkSpec {
        nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection: -0.8 }],
        sections: vec![TubeSection { nodes: [0,1],
            length_m: f64::from(delay) * g.sound_speed * DT,
            radius_m: 0.007, max_length_error_m: 1e-12 }],
        sound_speed_m_s: g.sound_speed, max_wave_memory_bytes: 1 << 20,
    }, gases: vec![g], viscothermal: vec![] })
}
fn channel(g: GasState, cells: usize) -> SlitChannel {
    SlitChannel::from_gas(&g, SlitChannelSpec { length_m: 0.015,
        width_m: 0.01, height_m: 0.001, transverse_cells: cells,
        sample_period_s: DT, exit: SlitExit::StaticPressure }).unwrap()
}
fn limits() -> TwoSidedLimits {
    TwoSidedLimits { maximum_energy_j: 1.0, energy_absolute_tolerance_j: 1e-16,
        energy_relative_tolerance: 1e-7, flow_absolute_tolerance_m3_s: 1e-16,
        flow_relative_tolerance: 1e-8 }
}
fn valve(plates: &[PlateApertureReduction;2], paired: bool, g: GasState, z: f64) -> LoadedAperture {
    let lay = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
        "synthetic supply/reed integration contact".into()).unwrap();
    // This is the MOVING reed gap's separate wall law. The fixed supply channel
    // is additional physical geometry, not a second copy of this slit.
    let slit = ViscousSlit::from_gas(&g, 0.001).unwrap();
    if paired {
        // Exact ZOH: this pair spends 87 % of the one-cent midpoint band on its two ducts alone (rate 6680 of
        // 7680 rad/s), so the channel's own friction (about 1e6 Pa s/m^3 at the reed's frequency against the
        // 9.6e5 left) cannot fit the midpoint band. These tests are about the shared flow and the work account,
        // not the midpoint band; the supply admission (parasitic root) still applies on this propagator.
        let mut a = PairedAperture::from_plates_with_propagator((*plates).clone(), g, z, DT, 4096,
            PairedApertureLimits { maximum_abs_displacement_m: 0.01,
                maximum_abs_velocity_m_s: 100.0, maximum_energy_j: 1.0,
                max_evaluations: 256, velocity_tolerance_m_s: 1e-13,
                relative_tolerance: 1e-11 }, [BladeState::default();2], lay,
            fs_couple::bernoulli_aperture::paired::BladePropagator::ExactZoh).unwrap();
        a.set_viscous_slit(slit).unwrap();
        a.into()
    } else {
        let mut a = DynamicAperture::from_plate(plates[0].clone(), g.density, z, DT, 4096,
            ApertureState { opening_m: plates[0].options().rest_opening_m,
                opening_velocity_m_s: 0.0 }, lay).unwrap();
        a.set_viscous_slit(slit).unwrap();
        a.into()
    }
}
fn model(plates: &[PlateApertureReduction;2], paired: bool, g: GasState,
    cells: usize, limits: TwoSidedLimits) -> TwoSidedAperture {
    let sides = ducts(g);
    let supply = channel(g, cells);
    let (inlet, z) = combined_inlet_state_with_supply(&sides, &supply).unwrap();
    TwoSidedAperture::with_supply_channel(valve(plates, paired, inlet, z),
        sides, supply, limits, None, &CancelGate::new_clock_free()).unwrap()
}
fn drive(p: f64) -> TwoSidedDrive {
    TwoSidedDrive { source_pressure_pa: p, blade_forces_n: [0.0;2] }
}
fn plates() -> [PlateApertureReduction;2] {
    [reduction(4e9, 900.0), reduction(4e9, 1200.0)]
}

#[test]
fn single_and_independent_material_reeds_close_one_flow_and_complete_channel_work() {
    let plates = plates();
    let gate = CancelGate::new_clock_free();
    for paired in [false, true] {
        let g = gas(1.0);
        let sides = ducts(g);
        let (_, z) = combined_inlet_state(&sides).unwrap();
        let mut plain = TwoSidedAperture::new(valve(&plates, paired, g, z),
            sides, limits(), &gate).unwrap();
        let first_without_passage = plain.step_under_gate(drive(3.0), &gate).unwrap();
        let mut m = model(&plates, paired, g, 32, limits());
        let mut energy = m.stored_energy_j();
        let initial_energy = energy;
        let mut work = 0.0;
        let mut loss = 0.0;
        let mut moving_release = false;
        for n in 0..512 {
            let pressure = if n < 128 { 3.0 } else { 0.0 };
            let f = m.step_under_gate(drive(pressure), &gate).unwrap();
            let supply = f.supply.unwrap();
            let q = f.aperture.bore_flow_m3_s;
            if n == 0 {
                assert_ne!(q, first_without_passage.aperture.bore_flow_m3_s);
                assert_ne!(f.aperture.state, first_without_passage.aperture.state);
            }
            assert!((supply.channel.flow_m3_s - q).abs() < 1e-14);
            assert!((f.upstream.inlet_flow_m3_s + q).abs() < 1e-14);
            assert!((f.downstream.inlet_flow_m3_s - q).abs() < 1e-14);
            assert_eq!(f.upstream_pressure_pa,
                pressure + f.upstream.inlet_pressure_pa - supply.pressure_drop_pa);
            assert_eq!(f.aperture.bore_pressure_pa, f.downstream.inlet_pressure_pa);
            assert!((f.aperture.pressure_work_j - DT * f.pressure_drop_pa * q).abs() < 1e-13);
            assert!((f.stored_energy_j - energy - f.storage_change_j).abs() < 1e-13);
            assert!(f.balance_residual_j().abs() < 1e-12);
            assert_eq!(supply.channel.exit_work_j, 0.0);
            assert!(supply.channel.viscous_dissipation_j >= 0.0);
            assert!(supply.channel.numerical_dissipation_j >= 0.0);
            assert_eq!(m.supply_channel().unwrap().last_step(), supply.channel);
            if n == 128 {
                assert_eq!(f.source_work_j, 0.0);
                moving_release = supply.channel.kinetic_energy_j > 0.0;
                assert!(m.supply_channel().unwrap().pressure_load().unwrap().history_pressure_pa != 0.0);
            }
            work += f.source_work_j;
            loss += f.dissipated_energy_j;
            energy = f.stored_energy_j;
        }
        assert!(moving_release, "source release must not erase the gas profile");
        assert!((energy - initial_energy + loss - work).abs() < 1e-11);
        match m.valve() {
            LoadedAperture::Single(a) =>
                assert_eq!(a.plate_reduction().unwrap().mass_kg(), plates[0].mass_kg()),
            LoadedAperture::Paired(a) => {
                assert_eq!(a.plates().unwrap()[0].mass_kg(), plates[0].mass_kg());
                assert_eq!(a.plates().unwrap()[1].mass_kg(), plates[1].mass_kg());
            }
            LoadedAperture::Modal(_) => unreachable!(),
        }
    }
}

#[test]
fn material_transport_reaches_real_reed_motion_and_the_existing_microphone() {
    let plates = plates();
    let light = gas(1.0);
    let viscous = gas(2.0);
    assert_eq!(light.density, viscous.density);
    assert_eq!(light.sound_speed, viscous.sound_speed);
    let receiver = CircularOutletReceiver { position_m: [0.0,0.0,0.1],
        radial_rings: 4, angular_points: 16, maximum_frequency_hz: 1000.0 };
    let render = |g, partition: usize| {
        let mut p = TwoSidedPressure::new(model(&plates, true, g, 32, limits()),
            ValveSide::Downstream, 1, receiver, RATE, 1024).unwrap();
        p.set_drive(drive(3.0)).unwrap();
        let mut signal = vec![0.0;1024];
        for block in signal.chunks_mut(partition) { p.block(block).unwrap(); }
        let states = match p.system().valve() {
            LoadedAperture::Paired(a) => *a.blade_states(), _ => unreachable!(),
        };
        let profile = p.system().supply_channel().unwrap().velocities_m_s().to_vec();
        (signal, states, profile)
    };
    let reference = render(light, 1024);
    let partitioned = render(light, 37);
    assert_eq!(reference, partitioned);
    let changed = render(viscous, 37);
    assert_ne!(reference.0, changed.0);
    assert_ne!(reference.1, changed.1);
    assert_ne!(reference.2, changed.2);
    assert!(reference.0.iter().any(|p| p.abs() > 1e-12));
}

#[test]
fn rejected_or_cancelled_trials_publish_neither_channel_nor_valve_nor_network() {
    let plates = plates();
    let gate = CancelGate::new_clock_free();
    let mut m = model(&plates, true, gas(1.0), 32, limits());
    m.step_under_gate(drive(3.0), &gate).unwrap();
    let profile = m.supply_channel().unwrap().velocities_m_s().to_vec();
    let last = m.supply_channel().unwrap().last_step();
    let energy = m.stored_energy_j();
    let up = *m.node_frame(ValveSide::Upstream, 0).unwrap();
    let down = *m.node_frame(ValveSide::Downstream, 0).unwrap();
    let blades = match m.valve() {
        LoadedAperture::Paired(a) => *a.blade_states(), _ => unreachable!(),
    };
    for p in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(m.step_under_gate(drive(p), &gate).is_err());
        assert_eq!(m.supply_channel().unwrap().velocities_m_s(), profile.as_slice());
        assert_eq!(m.supply_channel().unwrap().last_step(), last);
        assert_eq!(m.stored_energy_j(), energy);
        assert_eq!(*m.node_frame(ValveSide::Upstream, 0).unwrap(), up);
        assert_eq!(*m.node_frame(ValveSide::Downstream, 0).unwrap(), down);
        assert_eq!(m.valve().accepted_steps(), 1);
        match m.valve() {
            LoadedAperture::Paired(a) => assert_eq!(*a.blade_states(), blades), _ => unreachable!(),
        }
    }
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(m.step_under_gate(drive(3.0), &cancelled).is_err());
    assert_eq!(m.supply_channel().unwrap().last_step(), last);
    assert_eq!(*m.node_frame(ValveSide::Upstream, 0).unwrap(), up);
    assert_eq!(*m.node_frame(ValveSide::Downstream, 0).unwrap(), down);
    match m.valve() {
        LoadedAperture::Paired(a) => assert_eq!(*a.blade_states(), blades), _ => unreachable!(),
    }
    // Refuse AFTER all candidate physics exists, at the complete energy gate.
    let mut tiny = limits(); tiny.maximum_energy_j = 1e-24;
    let mut bounded = model(&plates, true, gas(1.0), 32, tiny);
    assert!(bounded.step_under_gate(drive(3.0), &gate).is_err());
    assert_eq!(bounded.valve().accepted_steps(), 0);
    assert!(bounded.supply_channel().unwrap().velocities_m_s().iter().all(|v| *v == 0.0));
    m.step_under_gate(drive(3.0), &gate).unwrap();
    assert_eq!(m.valve().accepted_steps(), 2);
}

// The paired fixture spends 87 % of the one-cent midpoint band on its two ducts alone (rate 6680 of 7680 rad/s at
// dt = 1/192 kHz). The channel's own friction at the plates' frequency (about 1.04e6 Pa s/m^3; 1.32e6 at the band
// edge) exceeds the 9.6e5 left, so the MIDPOINT propagator is correctly refused. `valve()` therefore builds this pair
// on exact ZOH; this test pins the refusal on the very same case so the band is covered, not walked around.
#[test]
fn the_paired_fixture_on_the_midpoint_propagator_is_refused_by_the_band_with_its_supply() {
    use fs_couple::acoustic_realize::AcousticRealizeError;
    let plates = plates();
    let g = gas(1.0);
    let sides = ducts(g);
    let lay = || Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
        "synthetic supply/reed integration contact".into()).unwrap();
    let limits = || PairedApertureLimits { maximum_abs_displacement_m: 0.01,
        maximum_abs_velocity_m_s: 100.0, maximum_energy_j: 1.0, max_evaluations: 256,
        velocity_tolerance_m_s: 1e-13, relative_tolerance: 1e-11 };
    // Ducts alone: admitted, inside the band.
    let (_, duct_z) = combined_inlet_state(&sides).unwrap();
    assert!(PairedAperture::from_plates(plates.clone(), g, duct_z, DT, 4096, limits(),
        [BladeState::default(); 2], lay()).is_ok());
    // With the passage's wall friction on top: refused, with the rate above 0.04 / dt.
    let (_, z) = combined_inlet_state_with_supply(&sides, &channel(g, 32)).unwrap();
    match PairedAperture::from_plates(plates, g, z, DT, 4096, limits(), [BladeState::default(); 2], lay()) {
        Err(AcousticRealizeError::TimeResolution(fs_phs::linear_exact::ResolvedBandError::Outside {
            rate_rad_s, dt_s, maximum_angular_step, .. })) => {
            assert_eq!(dt_s, DT);
            assert_eq!(maximum_angular_step, 0.04);
            assert!(rate_rad_s > maximum_angular_step / dt_s, "rate {rate_rad_s}");
        }
        other => panic!("the midpoint pair with its supply must leave the band, got {:?}", other.map(|_| ())),
    }
}

#[test]
fn supply_selection_refuses_stale_load_gas_clock_and_nonlinear_terminal() {
    let plates = plates();
    let gate = CancelGate::new_clock_free();
    let g = gas(1.0);
    let sides = ducts(g);
    let (_, old_z) = combined_inlet_state(&sides).unwrap();
    assert!(TwoSidedAperture::with_supply_channel(valve(&plates, true, g, old_z),
        sides.clone(), channel(g, 32), limits(), None, &gate).is_err());
    assert!(combined_inlet_state_with_supply(&sides, &channel(gas(2.0), 32)).is_err());
    let mut spec = *channel(g, 32).spec();
    spec.exit = SlitExit::UniformKineticHead;
    let nonlinear = SlitChannel::from_gas(&g, spec).unwrap();
    assert!(combined_inlet_state_with_supply(&sides, &nonlinear).is_err());
    spec.exit = SlitExit::StaticPressure;
    spec.sample_period_s *= 2.0;
    let wrong_clock = SlitChannel::from_gas(&g, spec).unwrap();
    let (_, z) = combined_inlet_state_with_supply(&sides, &wrong_clock).unwrap();
    assert!(TwoSidedAperture::with_supply_channel(valve(&plates, false, g, z),
        sides, wrong_clock, limits(), None, &gate).is_err());
}

#[test]
fn supplied_initial_fluid_momentum_is_retained_and_can_drive_an_unforced_material_reed() {
    let plates = plates();
    let gate = CancelGate::new_clock_free();
    let g = gas(1.0);
    let sides = ducts(g);
    let mut supply = channel(g, 32);
    for _ in 0..8 { supply.step(1.0).unwrap(); }
    let initial_profile = supply.velocities_m_s().to_vec();
    let initial_energy = supply.last_step().kinetic_energy_j;
    assert!(initial_energy > 0.0);
    let (_, z) = combined_inlet_state_with_supply(&sides, &supply).unwrap();
    let mut m = TwoSidedAperture::with_supply_channel(valve(&plates, true, g, z),
        sides, supply, limits(), None, &gate).unwrap();
    assert_eq!(m.supply_channel().unwrap().velocities_m_s(), initial_profile.as_slice());
    assert_eq!(m.stored_energy_j(), initial_energy);
    let mut loss = 0.0;
    let mut flow = false;
    for _ in 0..32 {
        let f = m.step_under_gate(drive(0.0), &gate).unwrap();
        assert_eq!(f.source_work_j, 0.0);
        loss += f.dissipated_energy_j;
        flow |= f.aperture.bore_flow_m3_s != 0.0;
    }
    assert!(flow, "retained channel momentum must react on the actual valve");
    assert!((m.stored_energy_j() + loss - initial_energy).abs() < 1e-13);
}

// ---- The passage's added mass against the reed: the parasitic root of the three-level coupling ----
//
// The passage steps on the reed's step-mean flow, so the coupled recursion has one real root besides the physical
// pair. For a free mass m with face area A against a passage of inertance I = rho L / A_t, mu = I A^2 / m, that root
// is -mu/(mu+2) per step (a rank-1 coupling gives mu = I sum(A_i^2 / m_i)). It is inside the unit circle for every
// finite mu, but tends to -1: a heavy passage rings at the Nyquist rate. The admission bounds it.
//
// Known answer: this root is -mu/(mu+2); here the jet is shut (a 2 um gap with a laminar slit law has admittance
// ~1e-4 of I/dt), which is the worst case. The reed is exact-ZOH, so the band admission does not apply.
mod added_mass {
    use super::*;
    use fs_couple::acoustic_realize::AcousticRealizeError;
    use fs_couple::bernoulli_aperture::{BernoulliAperture, dynamic::DynamicApertureSpec,
        paired::BladePropagator, two_sided::supply::maximum_parasitic_root};

    const M_REED: f64 = 1.0e-5;
    const K_REED: f64 = 500.0;
    const A_FACE: f64 = 1.25e-4;
    const H_REST: f64 = 2.0e-6;
    const V0: f64 = 2.0e-4;

    fn free_reed(inlet: GasState, z: f64) -> LoadedAperture {
        let lay = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
            "synthetic supply added-mass contact".into()).unwrap();
        let mut a = DynamicAperture::new_with_propagator(DynamicApertureSpec {
            aperture: BernoulliAperture { rest_opening_m: H_REST, width_m: 0.013,
                closing_pressure_pa: K_REED * H_REST / A_FACE },
            mass_kg: M_REED, stiffness_n_m: K_REED, damping_ratio: 0.02,
            density_kg_m3: inlet.density, impedance_pa_s_m3: z, time_step_s: DT, max_steps: 4096,
        }, ApertureState { opening_m: H_REST, opening_velocity_m_s: V0 }, lay, BladePropagator::ExactZoh).unwrap();
        a.set_viscous_slit(ViscousSlit::from_gas(&inlet, 0.001).unwrap()).unwrap();
        a.into()
    }
    fn long_ducts(g: GasState) -> [ApertureDuctSpec; 2] {
        // Reflections return after 8192 steps: outside every window below.
        [4096, 4096].map(|delay| ApertureDuctSpec { graph: TubeNetworkSpec {
            nodes: vec![NetworkNode::Inlet, NetworkNode::Termination { reflection: -0.8 }],
            sections: vec![TubeSection { nodes: [0,1],
                length_m: f64::from(delay) * g.sound_speed * DT,
                radius_m: 0.007, max_length_error_m: 1e-12 }],
            sound_speed_m_s: g.sound_speed, max_wave_memory_bytes: 1 << 20,
        }, gases: vec![g], viscothermal: vec![] })
    }
    /// A passage whose added-mass ratio against the reed is `mu`: L = mu m A_t / (rho A^2).
    fn passage(g: GasState, mu: f64) -> SlitChannel {
        let (width, height) = (0.01, 0.001);
        SlitChannel::from_gas(&g, SlitChannelSpec { length_m: mu * M_REED * (width * height) / (g.density * A_FACE * A_FACE),
            width_m: width, height_m: height, transverse_cells: 8, sample_period_s: DT,
            exit: SlitExit::StaticPressure }).unwrap()
    }
    fn system(mu: f64) -> Result<TwoSidedAperture, AcousticRealizeError> {
        let g = gas(1.0);
        let sides = long_ducts(g);
        let supply = passage(g, mu);
        let (inlet, z) = combined_inlet_state_with_supply(&sides, &supply).unwrap();
        TwoSidedAperture::with_supply_channel(free_reed(inlet, z), sides, supply, limits(), None,
            &CancelGate::new_clock_free())
    }
    /// Reed velocity and shared flow of `n` free steps after the velocity kick.
    fn free_run(mu: f64, n: usize) -> (Vec<f64>, Vec<f64>) {
        let mut m = system(mu).unwrap();
        let gate = CancelGate::new_clock_free();
        let (mut v, mut q) = (Vec::new(), Vec::new());
        for _ in 0..n {
            let f = m.step_under_gate(drive(0.0), &gate).unwrap();
            v.push(f.aperture.state.opening_velocity_m_s);
            q.push(f.aperture.bore_flow_m3_s);
        }
        (v, q)
    }
    /// Second difference: a geometric alternating component a r^n has D2_n = a r^(n-2) (r-1)^2, so D2_n / D2_(n-1)
    /// is r, while the smooth reed motion contributes only O((omega dt)^2).
    fn d2(x: &[f64], n: usize) -> f64 { x[n] - 2.0 * x[n - 1] + x[n - 2] }
    fn alternating_ratio(x: &[f64], from: usize, to: usize) -> f64 {
        let (mut num, mut den) = (0.0, 0.0);
        for n in from..to { num += d2(x, n) * d2(x, n - 1); den += d2(x, n - 1) * d2(x, n - 1); }
        num / den
    }
    fn maximum_mass_ratio() -> f64 {
        // Independent of the crate: rho_max = (1e-3)^(0.04 / tau), mu_max = 2 rho / (1 - rho).
        let rho = 1.0e-3_f64.powf(0.04 / core::f64::consts::TAU);
        2.0 * rho / (1.0 - rho)
    }

    #[test]
    fn the_parasitic_root_is_minus_mu_over_mu_plus_two_in_reed_velocity_and_flow() {
        // mu = 60 would be refused (mu_max = 44.5), which is the next test's point.
        for (mu, window) in [(10.0, 3..13), (40.0, 3..33)] {
            let (v, q) = free_run(mu, 64);
            let expected = -mu / (mu + 2.0);
            for (name, series) in [("velocity", &v), ("flow", &q)] {
                let r = alternating_ratio(series, window.start, window.end);
                assert!((r / expected - 1.0).abs() < 0.01,
                    "mu {mu}: measured {name} root {r:.5}, expected {expected:.5}");
            }
        }
    }

    #[test]
    fn a_supply_too_heavy_for_its_parasitic_root_is_refused_on_the_exact_propagator_too() {
        let maximum = maximum_parasitic_root();
        assert!((maximum - 1.0e-3_f64.powf(0.04 / core::f64::consts::TAU)).abs() < 1e-14);
        // Both sides of the bound: 1.02 mu_max is refused, 0.98 mu_max is admitted (below). mu = 60, the value the
        // first design note carried before the decay was fixed at -60 dB, is above mu_max = 44.5 and refused too.
        for mu in [1.02 * maximum_mass_ratio(), 60.0, 100.0, 1000.0] {
            match system(mu) {
                Err(AcousticRealizeError::SupplyParasiticRoot { mass_ratio, root, maximum_root }) => {
                    assert!((mass_ratio / mu - 1.0).abs() < 1e-9, "mu {mu}: reported {mass_ratio}");
                    assert!((root - mu / (mu + 2.0)).abs() < 1e-9);
                    assert_eq!(maximum_root, maximum);
                    assert!(root > maximum_root);
                }
                other => panic!("mu {mu} must be refused by its parasitic root, got {:?}", other.map(|_| ())),
            }
        }
        // Just inside the same bound the very same build is admitted.
        let inside = 0.98 * maximum_mass_ratio();
        assert!(inside / (inside + 2.0) < maximum);
        assert!(system(inside).is_ok());
    }

    #[test]
    fn at_the_bound_the_parasitic_mode_has_decayed_to_the_declared_decay_in_one_fast_period() {
        let steps = (core::f64::consts::TAU / 0.04).ceil() as usize; // one period of the fastest admitted mode
        let (v, q) = free_run(0.98 * maximum_mass_ratio(), 2 + steps + 8);
        for (name, series) in [("velocity", &v), ("flow", &q)] {
            let start = d2(series, 2).abs();
            let end = d2(series, 2 + steps).abs();
            assert!(end <= 1.0e-3 * start,
                "{name}: the alternating component kept {:.3e} of its start after {steps} steps", end / start);
        }
    }
}
