//! Synthetic physical wall inputs; these are model tests, not identified materials.
use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::dynamic::{ApertureState, ApertureTerminal, DynamicAperture, DynamicApertureSpec};
use fs_couple::bernoulli_aperture::network::{ApertureNetwork, ApertureNetworkFrame};
use fs_couple::bernoulli_aperture::tube::{TubeDrive, UniformTubeSpec};
use fs_couple::bernoulli_aperture::wall::{WallPatch, WallPin, WallState, lined_tube};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

const DT: f64 = 1e-5;
fn tube() -> UniformTubeSpec {
    UniformTubeSpec { length_m: 64.0 * 343.0 * DT, radius_m: 0.007, sound_speed_m_s: 343.0,
        terminal_reflection: -0.8, max_length_error_m: 1e-12, max_wave_memory_bytes: 1 << 20 }
}
fn wall(stiffness: f64, resistance: f64) -> WallPin {
    WallPin { surface_density: 0.2, stiffness_per_area: stiffness, resistance }
}
fn model(stiffness: f64, resistance: f64, budget: u64) -> (ApertureNetwork, Vec<WallPatch>) {
    let lined = lined_tube(tube(), wall(stiffness, resistance), 4, 1 << 20).unwrap();
    let spec = DynamicApertureSpec {
        aperture: BernoulliAperture { rest_opening_m: 4e-4, width_m: 0.013, closing_pressure_pa: 6000.0 },
        mass_kg: 1e-5, stiffness_n_m: 500.0, damping_ratio: 0.35, density_kg_m3: 1.2,
        impedance_pa_s_m3: lined.network.inlet_impedance(1.2).unwrap(), time_step_s: DT, max_steps: budget,
    };
    let lay = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
        "synthetic wall-coupling test contact".into()).unwrap().with_internal_loss(5.0).unwrap();
    let aperture = DynamicAperture::new(spec, ApertureState {
        opening_m: 4e-4, opening_velocity_m_s: 0.0,
    }, lay).unwrap();
    (ApertureNetwork::new(aperture, lined.network).unwrap(), lined.patches)
}
fn observe(model: &ApertureNetwork, patches: &[WallPatch]) -> Vec<WallState> {
    patches.iter().enumerate().map(|(i, p)| p.observe(model, i + 1).unwrap()).collect()
}
fn bits(f: ApertureNetworkFrame) -> Vec<u64> {
    let mut b = vec![f.aperture.step];
    b.extend([f.aperture.time_s, f.aperture.state.opening_m, f.aperture.state.opening_velocity_m_s,
        f.aperture.outgoing_pressure_pa, f.aperture.bore_pressure_pa, f.aperture.jet_flow_m3_s,
        f.aperture.swept_flow_m3_s, f.network.wave_stored_energy_j, f.network.load_stored_energy_j,
        f.network.interior_loss_j, f.network.terminal_loss_j, f.network.junction_residual_j,
        f.stored_energy_j, f.storage_change_j, f.dissipated_energy_j, f.upstream_work_j,
        f.body_work_j, f.balance_residual_j()].map(f64::to_bits));
    b
}

#[test]
fn centered_wall_geometry_covers_the_requested_area_once_and_reuses_the_owner_law() {
    let t = tube();
    let w = wall(2e7, 300.0);
    let lined = lined_tube(t, w, 4, 1 << 20).unwrap();
    assert_eq!(lined.network.nodes.len(), 6);
    assert_eq!(lined.network.sections.len(), 5);
    let area: f64 = lined.patches.iter().map(|p| p.area_m2).sum();
    assert!((area - 2.0 * core::f64::consts::PI * t.radius_m * t.length_m).abs() < 1e-16);
    let length: f64 = lined.network.sections.iter().map(|s| s.length_m).sum();
    assert!((length - t.length_m).abs() < 1e-15);
    assert_eq!(lined.network.sections[0].length_m * 2.0, lined.network.sections[1].length_m);
    for (i, section) in lined.network.sections.iter().enumerate() {
        assert_eq!(section.nodes, [i, i + 1]);
        assert_eq!(section.max_length_error_m, t.max_length_error_m / 5.0);
    }
    let omega = 1700.0;
    let patch = lined.patches[0];
    let load = patch.impedance().unwrap();
    let physical = fs_phs::wall_specific_impedance(&w, omega).unwrap().scale(1.0 / patch.area_m2);
    let acoustic = C64::new(load.resistance_pa_s_m3,
        -omega * load.inertance_pa_s2_m3 + 1.0 / (omega * load.compliance_m3_pa.unwrap()));
    assert!((physical - acoustic).abs() < 1e-12 * physical.abs());
    let (runtime, _) = model(2e7, 300.0, 1);
    let delays: Vec<_> = runtime.represented_sections().iter().map(|s| s.one_way_samples).collect();
    assert_eq!(delays, [8, 16, 16, 16, 8]);
}

#[test]
fn solved_wall_motion_satisfies_independent_mechanics_and_coupled_energy() {
    let mut omitted_storage = 0.0_f64;
    let mut omitted_loss = 0.0_f64;
    let mut max_motion = 0.0_f64;
    for resistance in [0.0, 300.0] {
        for stiffness in [2e7, 8e7] {
            let (mut model, patches) = model(stiffness, resistance, 1024);
            let mut old = observe(&model, &patches);
            for n in 0..1024 {
                let before = model.stored_energy_j();
                let f = model.step(TubeDrive {
                    upstream_pressure_pa: if n < 256 { 1200.0 } else { 0.0 },
                    body_flow_m3_s: if n < 512 { 2e-7 * (0.11 * f64::from(n)).sin() } else { 0.0 },
                }).unwrap();
                let next = observe(&model, &patches);
                let scale = (before + f.stored_energy_j + f.dissipated_energy_j
                    + f.upstream_work_j.abs() + f.body_work_j.abs()).max(f64::MIN_POSITIVE);
                let mut wall_energy = 0.0;
                let mut wall_change = 0.0;
                let mut wall_loss = 0.0;
                for (i, patch) in patches.iter().enumerate() {
                    let (a, b) = (old[i], next[i]);
                    let vm = 0.5 * (a.velocity_m_s + b.velocity_m_s);
                    let xm = 0.5 * (a.displacement_m + b.displacement_m);
                    let inertial = patch.wall.surface_density * (b.velocity_m_s - a.velocity_m_s) / DT;
                    let viscous = resistance * vm;
                    let spring = stiffness * xm;
                    let node = model.node_frame(i + 1).unwrap();
                    let force_scale = (inertial.abs() + viscous.abs() + spring.abs()
                        + node.pressure_pa.abs()).max(f64::MIN_POSITIVE);
                    assert!((inertial + viscous + spring - node.pressure_pa).abs() <= 2e-10 * force_scale);
                    let dx_scale = (a.displacement_m.abs() + b.displacement_m.abs()
                        + (DT * vm).abs()).max(f64::MIN_POSITIVE);
                    assert!((b.displacement_m - a.displacement_m - DT * vm).abs() <= 2e-12 * dx_scale);
                    let energy = 0.5 * patch.area_m2 * (patch.wall.surface_density * b.velocity_m_s.powi(2)
                        + stiffness * b.displacement_m.powi(2));
                    let loss = patch.area_m2 * resistance * vm * vm * DT;
                    assert!((energy - node.stored_energy_j).abs() <= 2e-12 * scale);
                    assert!((loss - node.absorbed_energy_j).abs() <= 2e-12 * scale);
                    assert!((node.load_flow_m3_s - patch.area_m2 * vm).abs()
                        <= 2e-12 * node.load_flow_m3_s.abs().max(f64::MIN_POSITIVE));
                    wall_energy += energy;
                    wall_change += energy - a.stored_energy_j;
                    wall_loss += loss;
                    max_motion = max_motion.max(b.displacement_m.abs());
                }
                assert!((wall_energy - f.network.load_stored_energy_j).abs() <= 2e-12 * scale);
                assert!((wall_loss - f.network.interior_loss_j).abs() <= 2e-12 * scale);
                assert!(f.balance_residual_j().abs() <= 3e-10 * scale);
                assert!(f.dissipated_energy_j >= 0.0);
                omitted_storage = omitted_storage.max((f.balance_residual_j() - wall_change).abs());
                omitted_loss = omitted_loss.max((f.balance_residual_j() - wall_loss).abs());
                if n >= 512 { assert!(f.stored_energy_j <= before + 3e-10 * scale); }
                old = next;
            }
        }
    }
    assert!(max_motion > 1e-6, "wall must actually move");
    assert!(omitted_storage > 1e-8 && omitted_loss > 1e-9, "uncredited wall work must fail");
}

#[test]
fn wall_stiffness_changes_valve_motion_only_after_the_physical_return_path() {
    let (mut soft, _) = model(2e7, 300.0, 512);
    let (mut stiff, _) = model(8e7, 300.0, 512);
    let mut first = None;
    let mut pressure_difference = 0.0_f64;
    let mut opening_difference = 0.0_f64;
    for n in 0..512 {
        let drive = TubeDrive { upstream_pressure_pa: 1200.0, body_flow_m3_s: 0.0 };
        let a = soft.step(drive).unwrap().aperture;
        let b = stiff.step(drive).unwrap().aperture;
        if a.outgoing_pressure_pa.to_bits() != b.outgoing_pressure_pa.to_bits() && first.is_none() { first = Some(n); }
        pressure_difference = pressure_difference.max((a.outgoing_pressure_pa - b.outgoing_pressure_pa).abs());
        opening_difference = opening_difference.max((a.state.opening_m - b.state.opening_m).abs());
    }
    assert_eq!(first, Some(16));
    assert!(pressure_difference > 1.0 && opening_difference > 1e-8);
}

#[test]
fn wall_and_wave_histories_survive_cancel_budget_resume_and_refused_samples() {
    let inputs: Vec<_> = (0..512).map(|n| TubeDrive {
        upstream_pressure_pa: if n < 128 { 800.0 } else { 0.0 },
        body_flow_m3_s: 2e-7 * (0.1 * f64::from(n)).sin(),
    }).collect();
    let gate = CancelGate::new_clock_free();
    let (mut complete, patches) = model(2e7, 300.0, 516);
    let expected: Vec<_> = inputs.iter().map(|&d| complete.step(d).unwrap()).collect();
    let (mut paused, _) = model(2e7, 300.0, 123);
    let sentinel = ApertureNetworkFrame::default();
    let mut actual = vec![sentinel; 512];
    paused.advance_block(&inputs[..80], &mut actual[..80], &gate).unwrap();
    let saved = observe(&paused, &patches);
    let cancel = CancelGate::new_clock_free(); cancel.request();
    let p = paused.advance_block(&inputs[80..], &mut actual[80..], &cancel).unwrap();
    assert_eq!((p.completed, p.terminal), (0, ApertureTerminal::Cancelled));
    assert_eq!(observe(&paused, &patches), saved);
    assert!(actual[80..].iter().all(|f| *f == sentinel));
    let p = paused.advance_block(&inputs[80..], &mut actual[80..], &gate).unwrap();
    assert_eq!((p.completed, p.terminal), (43, ApertureTerminal::BudgetExhausted));
    assert!(actual[123..].iter().all(|f| *f == sentinel));
    paused.extend_step_budget(516).unwrap();
    paused.advance_block(&inputs[123..], &mut actual[123..], &gate).unwrap();
    for (a, b) in actual.into_iter().zip(expected) { assert_eq!(bits(a), bits(b)); }
    assert_eq!(observe(&paused, &patches), observe(&complete, &patches));
    let old = observe(&paused, &patches);
    let energy = paused.stored_energy_j();
    for bad in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(paused.step(TubeDrive { upstream_pressure_pa: 500.0, body_flow_m3_s: bad }).is_err());
        assert_eq!(observe(&paused, &patches), old);
        assert_eq!(paused.stored_energy_j(), energy);
        assert_eq!(paused.aperture().accepted_steps(), 512);
    }
    assert!(paused.set_terminal_reflection(1, 0.0).is_err());
    let drive = TubeDrive { upstream_pressure_pa: 500.0, body_flow_m3_s: 0.0 };
    assert_eq!(bits(paused.step(drive).unwrap()), bits(complete.step(drive).unwrap()));
}

#[test]
fn geometry_and_observation_refuse_unfunded_or_incompatible_wall_models() {
    let t = tube(); let w = wall(2e7, 300.0);
    assert!(lined_tube(t, w, 0, 1 << 20).is_err());
    assert!(lined_tube(t, w, usize::MAX, usize::MAX).is_err());
    assert!(lined_tube(t, w, 4, 1).is_err());
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(WallPatch { area_m2: bad, wall: w }.impedance().is_err());
        assert!(lined_tube(t, WallPin { stiffness_per_area: bad, ..w }, 4, 1 << 20).is_err());
        assert!(lined_tube(t, WallPin { surface_density: bad, ..w }, 4, 1 << 20).is_err());
    }
    assert!(lined_tube(t, WallPin { resistance: -1.0, ..w }, 4, 1 << 20).is_err());
    let (network, patches) = model(2e7, 300.0, 1);
    assert!(patches[0].observe(&network, 0).is_err());
    assert!(patches[0].observe(&network, 5).is_err());
    assert!(patches[0].observe(&network, usize::MAX).is_err());
    let other = WallPatch { area_m2: 2.0 * patches[0].area_m2, ..patches[0] };
    assert!(other.observe(&network, 1).is_err());
    assert_eq!(patches[0].observe(&network, 1).unwrap(), WallState::default());
}
