//! Source -> shared nonlinear/contact mechanics -> physical receiver pressure.
//! The receiver is an explicit synthetic breathing-sphere fixture, not a
//! measured violin body or a claim about instrument timbre.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::TransverseDirections;
use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::runtime::MaterialStringRun;
use fs_couple::broadband_radiation::{
    BroadbandRadiationControls, ComplexShTrainingSample, DirectFarFieldHeldOutSample,
    HarmonicTimeConvention, RadiationSampleDiagnostics, SampledRadiationData, WeightPreset,
};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::render::plate::impact::{ImpactBody, ImpactConfig};
use fs_couple::render::plate::impact::audio::{
    ImpactListener, ImpactPressureRenderer, ImpactRadiation, ImpactSource, VelocityProjection,
};
use fs_dcontact::Obstacle;
use fs_exec::{CancelGate, Cx};
use fs_math::c64::C64;

const RATE: u32 = 48_000;
const SAMPLES: usize = 512;

fn radiation() -> ImpactRadiation {
    let medium = gas();
    let radius = 0.08;
    // Exact outgoing-wave solution for a uniformly pulsating sphere under
    // exp(-i*omega*t), in the SAME gas as the material string's drag.
    let transfer = |hz: f64| {
        let omega = core::f64::consts::TAU * hz;
        let ka = omega * radius / medium.sound_speed;
        C64::new(0.0, -omega * medium.density * radius * radius)
            * C64::new(1.0, -ka).recip()
            * C64::new(ka.cos(), -ka.sin())
    };
    let diagnostics = RadiationSampleDiagnostics {
        captured_fraction: 1.0, panels_per_wavelength: 100.0, condition_lower_bound: 1.0,
    };
    let root = (4.0 * core::f64::consts::PI).sqrt();
    let samples = SampledRadiationData {
        source_id: "analytic breathing sphere in the original string atmosphere".into(),
        harmonic_time_convention: HarmonicTimeConvention::ExpNegativeIOmegaT,
        l_max: 0,
        input_ids: vec!["radial_surface_velocity_m_s".into()],
        training: (0..24).map(|i| {
            let hz = 30.0 + f64::from(i) * 80.0;
            ComplexShTrainingSample {
                omega_rad_s: core::f64::consts::TAU * hz,
                coefficients_by_input: vec![vec![transfer(hz).scale(root)]],
                diagnostics,
            }
        }).collect(),
        held_out: (0..23).map(|i| {
            let hz = 70.0 + f64::from(i) * 80.0;
            DirectFarFieldHeldOutSample {
                omega_rad_s: core::f64::consts::TAU * hz,
                directions: vec![[0.,0.,1.], [0.,0.,-1.], [1.,0.,0.], [-1.,0.,0.],
                    [0.,1.,0.], [0.,-1.,0.], [1.,1.,1.], [-1.,-1.,-1.]],
                far_field_by_input: vec![vec![transfer(hz);8]],
                diagnostics,
            }
        }).collect(),
    };
    ImpactRadiation::from_velocity_samples(&samples, BroadbandRadiationControls {
        sample_rate_hz: f64::from(RATE), minimum_captured_fraction: 0.99,
        fit_order: 1, fit_iterations: 8, fit_weights: WeightPreset::Uniform,
        fit_d: true, far_field_signal_floor: 1e-10,
        maximum_normalized_error: 0.03, rms_normalized_error: 0.02,
    }, radius, medium.sound_speed).unwrap()
}

fn renderer<'a>(
    cx: &Cx<'_>, gate: &CancelGate, radiation: &'a ImpactRadiation,
    young: f64, bow: [f64;2], connected: bool,
) -> ImpactPressureRenderer<'a, MaterialStringRun> {
    let direction = TransverseDirections::try_new(bow, [1.0,0.0]).unwrap();
    let norm = (7800.0 * core::f64::consts::PI * 0.0003_f64.powi(2) * 0.5 * 0.5).sqrt();
    let mut state = vec![ModalAcousticState::default();4];
    // Same physical initial speed/kinetic energy; only its direction changes.
    state[0].velocity_m_sqrt_kg_per_s = 0.8 * bow[0] * norm;
    state[2].velocity_m_sqrt_kg_per_s = 0.8 * bow[1] * norm;
    let image = source(7800.0, young, 2).with_transverse_polarizations(direction).unwrap()
        .into_nonlinear(cx, gas(), RATE, 0.2, state, gate).unwrap();
    // Explicit radial effective mass and stiffness, independent of the string.
    // Its physical surface displacement is q/sqrt(m); initial acoustics are at rest.
    let receiver_mass = 0.04_f64;
    let receiver_shape = receiver_mass.sqrt().recip();
    let receiver = ImpactBody {
        potential: BodyPotential::Linear(vec![core::f64::consts::TAU * 300.0]),
        initial: vec![ModalAcousticState::default()], damping_per_s: vec![3.0],
    };
    let mut row = image.bridge_shape(0.25).unwrap();
    row.push(-receiver_shape);
    let contact = Obstacle::new(row, 1, 5, vec![0.0001], vec![1.0], 50_000.0, 1.5,
        "explicit synthetic filament/sphere surface contact".into()).unwrap()
        .with_internal_loss(0.05).unwrap();
    let mechanics = image.prepare(vec![receiver],
        if connected { vec![contact] } else { vec![] }, vec![], vec![],
        ImpactConfig { max_steps: SAMPLES as u64, ..config() }, gate).unwrap();
    ImpactPressureRenderer::new(mechanics, radiation,
        vec![VelocityProjection { input_id: "radial_surface_velocity_m_s".into(),
            weights: vec![0.0,0.0,0.0,0.0,receiver_shape] }],
        vec![0.0;5], ImpactListener {
            position_m: [0.0,0.0,1.0], maximum_delay_error_s: 1.01/f64::from(RATE),
            maximum_delay_samples: 4096, maximum_abs_pressure_pa: 100.0,
        }, RATE, 64).unwrap()
}

#[test]
fn material_and_force_direction_changes_reach_the_same_nonlinear_receiver_pressure_path() {
    with_cx(|cx, gate| {
        let radiation = radiation();
        let render = |young, bow, connected, chunk| {
            let mut r = renderer(cx, gate, &radiation, young, bow, connected);
            let initial = r.mechanics().system().stored_energy_j();
            let mut pressure = vec![0.0;SAMPLES];
            for block in pressure.chunks_mut(chunk) { r.block(block).unwrap(); }
            assert!(pressure.iter().all(|p| p.is_finite()));
            assert_eq!(r.mechanics().source().mode_count(), 4);
            assert_eq!(r.mechanics().source().transverse_mode_count(), 4);
            assert_eq!(r.mechanics().mode_count(), 5);
            assert_eq!(r.mechanics().samples_rendered(), SAMPLES as u64);
            assert!(r.mechanics().system().stored_energy_j() <= initial + 1e-8);
            assert!(r.last_mechanical_frame().unwrap().balance_residual_j.abs() < 1e-8);
            pressure
        };
        let oblique = render(200e9, [0.6,0.8], true, 64);
        let partitioned = render(200e9, [0.6,0.8], true, 7);
        assert_eq!(oblique.iter().map(|p| p.to_bits()).collect::<Vec<_>>(),
            partitioned.iter().map(|p| p.to_bits()).collect::<Vec<_>>());
        let disconnected = render(200e9, [0.6,0.8], false, 64);
        let orthogonal = render(200e9, [0.0,1.0], true, 64);
        assert!(disconnected.iter().all(|p| *p == 0.0));
        // Shared tension modulates an excited plane; it does not fabricate
        // displacement in an exactly unexcited orthogonal plane.
        assert!(orthogonal.iter().all(|p| *p == 0.0));
        let aligned = render(200e9, [1.0,0.0], true, 64);
        let softer = render(50e9, [0.6,0.8], true, 64);
        let norm = |p: &[f64]| p.iter().map(|p| p*p).sum::<f64>().sqrt();
        let difference = |a: &[f64], b: &[f64]| a.iter().zip(b)
            .map(|(a,b)| (a-b)*(a-b)).sum::<f64>().sqrt();
        assert!(norm(&oblique) > 1e-8);
        assert!(difference(&oblique, &aligned) > 1e-4*norm(&oblique));
        assert!(difference(&oblique, &softer) > 1e-5*norm(&oblique));
    });
}
