//! Tests exercise the actual valve/network owners, not a stand-in output filter.
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::paired::{ApertureBlade, BladeState, PairedAperture,
    PairedApertureLimits, PairedApertureSpec};
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
    PairedAperture::new(PairedApertureSpec {
        rest_opening_m: 4e-4, width_m: 0.013, blades: b,
        density_kg_m3: g.density, impedance_pa_s_m3: z, time_step_s: DT, max_steps: steps,
        limits: PairedApertureLimits {
            maximum_abs_displacement_m: 0.01, maximum_abs_velocity_m_s: 100.0,
            maximum_energy_j: 1.0, max_evaluations: 256,
            velocity_tolerance_m_s: 1e-13, relative_tolerance: 1e-11,
        },
    }, [BladeState::default(); 2], contact()).unwrap()
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
