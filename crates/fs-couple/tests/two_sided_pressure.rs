//! Tests exercise the actual valve/network owners, not a stand-in output filter.
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::paired::{ApertureBlade, BladePropagator, BladeState,
    PairedAperture, PairedApertureLimits, PairedApertureSpec};
use fs_couple::bernoulli_aperture::two_sided::{ApertureDuctSpec,
    TwoSidedAperture, TwoSidedDrive, TwoSidedLimits, ValveSide, combined_inlet_state};
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
    // These ~1.1 kHz, strongly loaded blades lie outside the midpoint one-cent
    // band at 48 kHz (DEF-04), so their linear blocks are propagated exactly.
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

fn pressure_renderer() -> fs_couple::bernoulli_aperture::two_sided::observer::TwoSidedPressure {
    use fs_couple::bernoulli_aperture::two_sided::observer::TwoSidedPressure;
    use fs_couple::pcm_wav::baffled::CircularOutletReceiver;
    TwoSidedPressure::new(system(12,1.3,256), ValveSide::Downstream, 1,
        CircularOutletReceiver { position_m: [0.0,0.0,0.25], radial_rings: 4,
            angular_points: 16, maximum_frequency_hz: 1000.0 }, 48_000, 32).unwrap()
}

#[test]
fn exterior_pressure_uses_retarded_outlet_flow_and_is_block_partition_invariant() {
    use fs_couple::pcm_wav::observation::PressureRenderer;
    let mut a = pressure_renderer();
    let mut b = pressure_renderer();
    a.set_drive(drive(50.0)).unwrap(); b.set_drive(drive(50.0)).unwrap();
    let mut block = [0.0; 256];
    for piece in block.chunks_mut(32) { a.block(piece).unwrap(); }
    let mut singles = [0.0; 256];
    for piece in singles.chunks_mut(1) { b.block(piece).unwrap(); }
    assert_eq!(block, singles);
    // The 32-sample tube transit precedes the explicitly retained receiver delay.
    assert!(block[..32].iter().all(|p| *p == 0.0));
    assert!(block[64..].iter().any(|p| p.abs() > 0.0));
    assert_eq!(a.samples_rendered(),256);
    assert!(a.validate_sample_count(1).is_err());
}

#[test]
fn invalid_pressure_callbacks_do_not_advance_or_overwrite_output() {
    use fs_couple::pcm_wav::observation::PressureRenderer;
    let mut a = pressure_renderer();
    assert!(a.validate_sample_rate(24_000).is_err());
    assert!(a.set_drive(drive(f64::NAN)).is_err());
    let mut too_large = [987.0; 33];
    assert!(a.block(&mut too_large).is_err());
    assert_eq!(too_large,[987.0;33]);
    assert_eq!(a.samples_rendered(),0);
    assert_eq!(a.system().valve().accepted_steps(),0);
    a.set_drive(drive(50.0)).unwrap();
    a.block(&mut [0.0;32]).unwrap();
    assert_eq!(a.samples_rendered(),32);
}

mod scheduled {
    use super::*;
    use fs_couple::bernoulli_aperture::performance::AperturePerformanceConfig;
    use fs_couple::bernoulli_aperture::two_sided::observer::schedule::ScheduledTwoSidedPressure;
    use fs_couple::pcm_wav::observation::PressureRenderer;
    use fs_couple::render::{ControlDelta, RenderError};
    use fs_scenario::gesture::GestureSchedule;

    fn phrase() -> GestureSchedule {
        GestureSchedule::from_canonical_bytes(b"frankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t1\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t0e0\nevents\t2\nevent\t0e0\t1e-3\tpressure-pa\t5e1\nevent\t3e-3\t0e0\tpressure-pa\t0e0\n").unwrap()
    }
    fn config() -> AperturePerformanceConfig {
        AperturePerformanceConfig { sample_rate_hz: 48_000, samples: 256, max_block: 32,
            max_compile_work: 100_000, max_controls: 16 }
    }
    fn scheduled() -> ScheduledTwoSidedPressure {
        ScheduledTwoSidedPressure::new(pressure_renderer(), phrase(), config()).unwrap()
    }

    #[test]
    fn canonical_controls_match_direct_loaded_physics_and_retain_independent_forces() {
        let forces = [1e-5, -2e-5];
        let mut source = pressure_renderer();
        source.set_drive(TwoSidedDrive { source_pressure_pa: 0.0, blade_forces_n: forces }).unwrap();
        let mut p = ScheduledTwoSidedPressure::new(source, phrase(), config()).unwrap();
        let controls = p.pending_controls().to_vec();
        assert!(controls.iter().any(|e| e.sample % 32 != 0));
        let mut direct = pressure_renderer();
        let mut held = 0.0;
        for sample in 0..256 {
            for event in controls.iter().filter(|e| e.sample == sample) {
                if let ControlDelta::SetBlowingPressure { pressure_pa, .. } = event.delta { held = pressure_pa; }
            }
            direct.set_drive(TwoSidedDrive { source_pressure_pa: held, blade_forces_n: forces }).unwrap();
            let (mut actual, mut expected) = ([0.0], [0.0]);
            p.block(&mut actual).unwrap(); direct.block(&mut expected).unwrap();
            assert_eq!(actual, expected);
            assert_eq!(p.source().last_frame(), direct.last_frame());
        }
        assert!(p.pending_controls().is_empty());
        assert_eq!(p.applied_controls().len(), controls.len());
    }

    #[test]
    fn pressure_phrase_is_callback_independent_and_release_preserves_both_wave_histories() {
        let mut reference = scheduled();
        let mut expected = [0.0; 256];
        for slot in &mut expected { reference.block(core::slice::from_mut(slot)).unwrap(); }
        for size in [7, 31, 32] {
            let mut p = scheduled();
            let mut actual = [0.0; 256];
            for block in actual.chunks_mut(size) { p.block(block).unwrap(); }
            assert_eq!(actual, expected);
            assert_eq!(p.source().last_frame(), reference.source().last_frame());
            assert_eq!(p.remaining_samples(), 0);
            assert!(p.source().system().stored_energy_j() > 0.0);
            assert_eq!(p.source().last_frame().unwrap().source_work_j, 0.0);
            let mut untouched = [123.0];
            assert!(p.block(&mut untouched).is_err());
            assert_eq!(untouched, [123.0]);
        }
    }

    #[test]
    fn schedule_admission_and_failed_physics_do_not_publish_future_controls() {
        let mut zero = config(); zero.samples = 0;
        assert!(ScheduledTwoSidedPressure::new(pressure_renderer(), phrase(), zero).is_err());
        let mut limited = config(); limited.max_controls = 1;
        assert!(ScheduledTwoSidedPressure::new(pressure_renderer(), phrase(), limited).is_err());
        let mut already_set = pressure_renderer(); already_set.set_drive(drive(1.0)).unwrap();
        assert!(ScheduledTwoSidedPressure::new(already_set, phrase(), config()).is_err());
        let mut late = pressure_renderer(); late.block(&mut [0.0]).unwrap();
        assert!(ScheduledTwoSidedPressure::new(late, phrase(), config()).is_err());
        let mut p = scheduled();
        assert!(p.set_blade_forces([f64::NAN, 0.0]).is_err());
        assert_eq!(p.samples_rendered(), 0);
        let controls = p.pending_controls().len();
        p.set_blade_forces([f64::MAX, 0.0]).unwrap();
        let mut output = [123.0; 8];
        assert!(p.block(&mut output).is_err());
        assert_eq!(output, [123.0; 8]);
        assert_eq!(p.samples_rendered(), 0);
        assert_eq!(p.pending_controls().len(), controls);
        assert!(p.applied_controls().is_empty());
        assert!(matches!(p.validate_sample_count(0), Err(RenderError::Poisoned)));
    }
}
