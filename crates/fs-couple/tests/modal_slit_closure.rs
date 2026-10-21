//! Bead frankensim-0ll17: the modal slit's exact full-C closure scenario (modal_aperture propagator test) refused a
//! step with "uniform aperture has unresolved gap momentum". The scalar gap solve bisects on a residual that contains
//! the zero-gap slit lay's discrete contact force. A reed held shut by a steady pressure takes tiny gap steps, so that
//! force is evaluated as a potential secant over a tiny increment in deep contact.
//!
//! The first test pins the property the solve needs: at every step size, including the 1e-14 |q| threshold where the
//! discrete gradient drops its correction, the secant stays the force of the midpoint penetration. The second replays
//! the closure scenario and reports the state, gap velocity and contact force of any step the solve refuses.
use fs_couple::bernoulli_aperture::dynamic::ApertureDrive;
use fs_couple::bernoulli_aperture::modal::{ApertureMode, ModalAperture, ModalApertureSpec};
use fs_couple::bernoulli_aperture::paired::{BladePropagator, BladeState, PairedApertureLimits};
use fs_couple::bernoulli_aperture::viscous::ViscousSlit;
use fs_couple::unilateral_contact::slit_contact_discrete_coefficients;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

const K: f64 = 1e8;
const ALPHA: f64 = 2.0;

fn contact() -> Obstacle {
    Obstacle::new(
        vec![-1.],
        1,
        1,
        vec![0.],
        vec![1.],
        K,
        ALPHA,
        "explicit synthetic modal slit contact".into(),
    )
    .unwrap()
}

/// Largest |secant - midpoint force| / midpoint force over the step sizes a resting closed reed takes.
fn worst_secant_error(lay: &Obstacle, opening: f64) -> (f64, f64) {
    let mut worst = (0.0_f64, 0.0_f64);
    for exponent in 8..=20 {
        for sign in [1.0, -1.0] {
            for mantissa in [1.0, 1.5, 3.0] {
                let dx = sign * mantissa * 10f64.powi(-exponent);
                let (elastic, _) = slit_contact_discrete_coefficients(lay, opening, opening + dx)
                    .expect("contact coefficients");
                let midpoint = K * (-(opening + 0.5 * dx)).powf(ALPHA);
                let error = (elastic - midpoint).abs() / midpoint;
                if error > worst.0 {
                    worst = (error, dx);
                }
            }
        }
    }
    worst
}

#[test]
fn zero_gap_slit_lay_secant_is_the_midpoint_force_at_tiny_steps() {
    // Penetrations of a reed held shut by 0.1-1 N on the fixture's K = 1e8 N/m^2, alpha = 2 lay.
    let lay = contact();
    for opening in [-1.0e-5, -3.2e-5, -5.0e-5, -1.0e-4] {
        let (error, dx) = worst_secant_error(&lay, opening);
        eprintln!("opening {opening:e}: worst relative secant error {error:e} at dx {dx:e}");
        assert!(
            error <= 1e-8,
            "opening {opening:e}: the secant leaves the midpoint force by {error:e} at dx {dx:e}"
        );
    }
}

fn limits() -> PairedApertureLimits {
    PairedApertureLimits {
        maximum_abs_displacement_m: 0.01,
        maximum_abs_velocity_m_s: 100.,
        maximum_energy_j: 1.,
        max_evaluations: 256,
        velocity_tolerance_m_s: 1e-13,
        relative_tolerance: 1e-11,
    }
}

/// The modal_aperture fixture's three-mode slit, at the closure scenario's clock and rest opening.
fn closure_spec() -> ModalApertureSpec {
    let mode =
        |mass_kg, stiffness_n_m, pressure_area_m2, gap_projection, force_projection| ApertureMode {
            mass_kg,
            stiffness_n_m,
            damping_n_s_m: 0.0,
            pressure_area_m2,
            gap_projection,
            force_projection,
        };
    ModalApertureSpec {
        rest_opening_m: 4e-5,
        width_m: 0.013,
        modes: vec![
            mode(1e-5, 500., 1.25e-4, 1., [1., 0.]),
            mode(8e-6, 2500., -4e-5, -0.2, [-0.3, 0.7]),
            mode(6e-6, 6000., 7e-5, 0., [0.2, -0.4]),
        ],
        density_kg_m3: 1.2,
        impedance_pa_s_m3: 1e6,
        time_step_s: 1.0 / 48_000.0,
        max_steps: 4096 * 16,
        limits: limits(),
    }
}

fn closure_model(state: Vec<BladeState>) -> ModalAperture {
    let mut a = ModalAperture::new_with_propagator(
        closure_spec(),
        state,
        contact(),
        BladePropagator::ExactZoh,
    )
    .unwrap()
    .with_viscosity(
        vec![0.03, 0.005, -0.002, 0.005, 0.02, 0.003, -0.002, 0.003, 0.01],
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    a.set_viscous_slit(ViscousSlit::new(1.81e-5, 0.001).unwrap())
        .unwrap();
    a
}

fn drive(p: f64) -> ApertureDrive {
    ApertureDrive {
        upstream_pressure_pa: p,
        incoming_pressure_pa: 0.,
        body_flow_m3_s: 0.,
    }
}

#[test]
fn closure_scenario_reports_any_refused_step() {
    let mut a = closure_model(vec![BladeState::default(); 3]);
    let gap_velocity = |modes: &[BladeState]| {
        closure_spec()
            .modes
            .iter()
            .zip(modes)
            .map(|(m, s)| m.gap_projection * s.velocity_m_s)
            .sum::<f64>()
    };
    let mut last = (f64::NAN, f64::NAN, f64::NAN, f64::NAN);
    for sample in 0..1_024 {
        let before = a.modal_states().to_vec();
        let pressure = if sample < 512 { 2_000.0 } else { -100.0 };
        match a.step(drive(pressure), [0.001, -0.0003]) {
            Ok(f) => {
                last = (
                    f.aperture.state.opening_m,
                    gap_velocity(&f.modes),
                    f.contact_force_n,
                    f.pressure_drop_pa,
                );
                assert!(f.work_certificate.residual <= f.work_certificate.effective_budget);
                assert!(f.flow_certificate.residual <= f.flow_certificate.effective_budget);
            }
            Err(e) => {
                let (opening, velocity, force, dp) = last;
                let (error, dx) = worst_secant_error(&contact(), opening.min(-1e-12));
                panic!(
                    "sample {sample} (drive {pressure} Pa) refused: {e:?}\n  state before: {before:?}\n  \
                     last accepted frame: opening {opening:e} m, gap velocity {velocity:e} m/s, contact force \
                     {force:e} N, pressure drop {dp:e} Pa\n  secant at that opening: worst relative error \
                     {error:e} at dx {dx:e}"
                );
            }
        }
    }
}
