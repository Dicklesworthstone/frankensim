//! G1/G5 analytic disclosure of the current linear owners, not IRP acceptance.
//! Existing midpoint paths remain unmigrated. This audit must not turn a
//! measured >1-cent defect into a fidelity claim. The network-piston and piano
//! rows require their separate owner-local tests; they are explicitly unmeasured
//! here rather than being represented by a surrogate scalar recurrence.

use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::dynamic::{
    ApertureDrive, ApertureState, DynamicAperture, DynamicApertureSpec,
};
use fs_couple::bernoulli_aperture::modal::{ApertureMode, ModalAperture, ModalApertureSpec};
use fs_couple::bernoulli_aperture::paired::{
    ApertureBlade, BladeState, PairedAperture, PairedApertureLimits, PairedApertureSpec,
};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::plate::impact::{BodyPotential, ImpactBody, ImpactConfig, ImpactSystem};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig,
};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};

const AMPLITUDE: f64 = 1e-6;
const REST: f64 = 1e-3;

#[derive(Clone, Copy, Debug)]
enum Integrator {
    ExactZoh,
    GonzalezMidpoint,
    Trapezoidal,
    FirstOrder,
}
#[derive(Clone, Copy)]
enum Owner {
    Phs,
    Impact,
    Dynamic,
    Paired,
    Modal,
    Coupled,
    Exact,
    Speaker,
    SeparateTest,
}
struct AuditRow {
    path: &'static str,
    integrator: Integrator,
    bound: Option<f64>,
    owner: Owner,
}

// Inventory re-derived from HEAD call sites, including prepared Gonzalez users.
// `None` is an absent admission guard, never an invented infinite validity band.
const PATHS: &[AuditRow] = &[
    AuditRow {
        path: "acoustic_realize::realize_coupled_kc / fs_phs::step",
        integrator: Integrator::GonzalezMidpoint,
        bound: None,
        owner: Owner::Phs,
    },
    AuditRow {
        path: "glottis / fs_phs::step",
        integrator: Integrator::GonzalezMidpoint,
        bound: None,
        owner: Owner::Phs,
    },
    AuditRow {
        path: "speaker::TsDriver / exact full quadratic step",
        integrator: Integrator::ExactZoh,
        bound: None,
        owner: Owner::Speaker,
    },
    AuditRow {
        path: "thin_plate::VkBody and PlateBank cavity / fs_phs::step",
        integrator: Integrator::GonzalezMidpoint,
        bound: None,
        owner: Owner::Phs,
    },
    AuditRow {
        path: "acoustic_realize::realize_coupled_ode / step_chain_inlet / step_plate_chain",
        integrator: Integrator::GonzalezMidpoint,
        bound: None,
        owner: Owner::Phs,
    },
    AuditRow {
        path: "acoustic_realize::realize_dirac_join / step_descriptor",
        integrator: Integrator::GonzalezMidpoint,
        bound: None,
        owner: Owner::SeparateTest,
    },
    AuditRow {
        path: "acoustic_realize::realize_blown_ode and realize_reed_ode",
        integrator: Integrator::GonzalezMidpoint,
        bound: None,
        owner: Owner::Phs,
    },
    AuditRow {
        path: "brass_loop::lip and cavity_phs::realize_plate_cavity / fs_phs::step",
        integrator: Integrator::GonzalezMidpoint,
        bound: None,
        owner: Owner::Phs,
    },
    AuditRow {
        path: "render::plate::ImpactSystem",
        integrator: Integrator::GonzalezMidpoint,
        bound: Some(0.9 * core::f64::consts::PI),
        owner: Owner::Impact,
    },
    AuditRow {
        path: "bernoulli_aperture::DynamicAperture",
        integrator: Integrator::GonzalezMidpoint,
        bound: None,
        owner: Owner::Dynamic,
    },
    AuditRow {
        path: "bernoulli_aperture::PairedAperture",
        integrator: Integrator::GonzalezMidpoint,
        bound: None,
        owner: Owner::Paired,
    },
    AuditRow {
        path: "bernoulli_aperture::ModalAperture",
        integrator: Integrator::GonzalezMidpoint,
        bound: None,
        owner: Owner::Modal,
    },
    AuditRow {
        path: "bernoulli_aperture::uniform (through ModalAperture)",
        integrator: Integrator::GonzalezMidpoint,
        bound: None,
        owner: Owner::Modal,
    },
    AuditRow {
        path: "render::schedule::force::CoupledModalSystem connections",
        integrator: Integrator::Trapezoidal,
        bound: Some(0.9 * core::f64::consts::PI),
        owner: Owner::Coupled,
    },
    AuditRow {
        path: "modal_acoustic_time::ModalAcousticTimeModel",
        integrator: Integrator::ExactZoh,
        bound: Some(0.9 * core::f64::consts::PI),
        owner: Owner::Exact,
    },
    AuditRow {
        path: "cavity_phs::piston::GasPiston",
        integrator: Integrator::GonzalezMidpoint,
        bound: None,
        owner: Owner::SeparateTest,
    },
    AuditRow {
        path: "cavity_phs::piston::thermal::boundary::network",
        integrator: Integrator::FirstOrder,
        bound: None,
        owner: Owner::SeparateTest,
    },
    AuditRow {
        path: "engine-piano::PianoAudio 1x",
        integrator: Integrator::Trapezoidal,
        bound: None,
        owner: Owner::SeparateTest,
    },
    AuditRow {
        path: "engine-piano::PianoAudio 2x",
        integrator: Integrator::Trapezoidal,
        bound: None,
        owner: Owner::SeparateTest,
    },
    AuditRow {
        path: "engine-piano::PianoAudio 4x",
        integrator: Integrator::Trapezoidal,
        bound: None,
        owner: Owner::SeparateTest,
    },
];

fn limits() -> PairedApertureLimits {
    PairedApertureLimits {
        maximum_abs_displacement_m: 1.0,
        maximum_abs_velocity_m_s: 1e4,
        maximum_energy_j: 1e6,
        max_evaluations: 256,
        velocity_tolerance_m_s: 1e-13,
        relative_tolerance: 1e-12,
    }
}
fn obstacle() -> Obstacle {
    Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        1e5,
        2.0,
        "authored analytic warping fixture; contact remains inactive".into(),
    )
    .unwrap()
}
fn impact(
    omega: f64,
    dt: f64,
    z: [f64; 2],
) -> Result<ImpactSystem, fs_couple::render::plate::impact::ImpactError> {
    ImpactSystem::new(
        vec![ImpactBody {
            potential: BodyPotential::Linear(vec![omega]),
            initial: vec![ModalAcousticState {
                displacement_m_sqrt_kg: AMPLITUDE * z[0],
                velocity_m_sqrt_kg_per_s: AMPLITUDE * omega * z[1],
            }],
            damping_per_s: vec![0.0],
        }],
        vec![],
        vec![],
        vec![],
        ImpactConfig {
            dt_s: dt,
            max_steps: 1,
            maximum_energy_j: 1e6,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-10,
            maximum_generalized_force: 1e6,
        },
    )
}
fn modal_model(omega: f64, rate: u32, z: [f64; 2]) -> ModalAcousticTimeModel {
    let mut model = ModalAcousticTimeModel::try_new(
        rate,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::new(0.0, 0.0),
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    model
        .restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: AMPLITUDE * z[0],
            velocity_m_sqrt_kg_per_s: AMPLITUDE * omega * z[1],
        }])
        .unwrap();
    model
}
fn normalized(x: &[f64], omega: f64) -> [f64; 2] {
    [x[0] / AMPLITUDE, x[1] / (AMPLITUDE * omega)]
}

fn advance(owner: Owner, omega: f64, rate: u32, z: [f64; 2]) -> [f64; 2] {
    let dt = 1.0 / f64::from(rate);
    let (q, v) = (AMPLITUDE * z[0], AMPLITUDE * omega * z[1]);
    match owner {
        Owner::Phs => {
            let model = fs_phs::mass_spring_damper(1.0, omega * omega, 0.0).unwrap();
            let x = fs_phs::step(&model, &[q, v], &[0.0], dt).unwrap().x;
            normalized(&x, omega)
        }
        Owner::Impact => {
            let mut model = impact(omega, dt, z).unwrap();
            model.step(&[0.0], &CancelGate::new()).unwrap();
            normalized(model.state(), omega)
        }
        Owner::Exact => {
            let mut model = modal_model(omega, rate, z);
            model.step(&[0.0]).unwrap();
            let x = model.states()[0];
            normalized(
                &[x.displacement_m_sqrt_kg, x.velocity_m_sqrt_kg_per_s],
                omega,
            )
        }
        Owner::Modal => {
            // Zero signed pressure area is an explicit uncoupled physical port;
            // gap motion still runs the actual uniform contact/flow owner.
            let spec = ModalApertureSpec {
                rest_opening_m: REST,
                width_m: 0.01,
                modes: vec![ApertureMode {
                    mass_kg: 1.0,
                    stiffness_n_m: omega * omega,
                    damping_n_s_m: 0.0,
                    pressure_area_m2: 0.0,
                    gap_projection: 1.0,
                    force_projection: [0.0; 2],
                }],
                density_kg_m3: 1.2,
                impedance_pa_s_m3: 1e6,
                time_step_s: dt,
                max_steps: 1,
                limits: limits(),
            };
            let mut model = ModalAperture::new(
                spec,
                vec![BladeState {
                    displacement_m: q,
                    velocity_m_s: v,
                }],
                obstacle(),
            )
            .unwrap();
            let frame = model.step(ApertureDrive::default(), [0.0; 2]).unwrap();
            normalized(
                &[frame.modes[0].displacement_m, frame.modes[0].velocity_m_s],
                omega,
            )
        }
        Owner::Dynamic | Owner::Paired => {
            let mass = 1e-5;
            let stiffness = mass * omega * omega;
            let area = 1e-6;
            // Manufacture zero differential pressure by cancelling the swept
            // volume with an externally prescribed body flow. This isolates
            // the actual mechanical map without approximating the jet law.
            let h = 0.5 * dt;
            let vm = (v - h * omega * omega * q) / (1.0 + h * h * omega * omega);
            let drive = ApertureDrive {
                body_flow_m3_s: area * vm,
                ..ApertureDrive::default()
            };
            if matches!(owner, Owner::Dynamic) {
                let mut model = DynamicAperture::new(
                    DynamicApertureSpec {
                        aperture: BernoulliAperture {
                            rest_opening_m: REST,
                            width_m: 0.01,
                            closing_pressure_pa: stiffness * REST / area,
                        },
                        mass_kg: mass,
                        stiffness_n_m: stiffness,
                        damping_ratio: 0.0,
                        density_kg_m3: 1.2,
                        impedance_pa_s_m3: 1e6,
                        time_step_s: dt,
                        max_steps: 1,
                    },
                    ApertureState {
                        opening_m: REST + q,
                        opening_velocity_m_s: v,
                    },
                    obstacle(),
                )
                .unwrap();
                let frame = model.step(drive).unwrap();
                assert!(
                    frame.bore_pressure_pa.abs() < 1e-7,
                    "manufactured pressure isolation failed"
                );
                normalized(
                    &[
                        frame.state.opening_m - REST,
                        frame.state.opening_velocity_m_s,
                    ],
                    omega,
                )
            } else {
                let blade = ApertureBlade {
                    mass_kg: mass,
                    stiffness_n_m: stiffness,
                    damping_n_s_m: 0.0,
                    pressure_area_m2: area,
                };
                let mut model = PairedAperture::new(
                    PairedApertureSpec {
                        rest_opening_m: REST,
                        width_m: 0.01,
                        blades: [blade; 2],
                        density_kg_m3: 1.2,
                        impedance_pa_s_m3: 1e6,
                        time_step_s: dt,
                        max_steps: 1,
                        limits: limits(),
                    },
                    [
                        BladeState {
                            displacement_m: q,
                            velocity_m_s: v,
                        },
                        BladeState::default(),
                    ],
                    obstacle(),
                )
                .unwrap();
                let frame = model.step(drive).unwrap();
                assert!(
                    frame.pressure_drop_pa.abs() < 1e-7,
                    "manufactured pressure isolation failed"
                );
                normalized(
                    &[frame.blades[0].displacement_m, frame.blades[0].velocity_m_s],
                    omega,
                )
            }
        }
        Owner::Coupled => {
            let config = ModalCouplingConfig {
                max_modes: 2,
                max_connections: 1,
                max_setup_terms: 64,
                nyquist_guard_fraction: 0.9,
                maximum_total_energy_j: 1e6,
                maximum_abs_pressure_pa: 1e6,
                maximum_abs_connection_force_n: 1e6,
                solve_relative_tolerance: 1e-11,
                energy_absolute_tolerance_j: 1e-12,
                energy_relative_tolerance: 1e-10,
            };
            let connection = ModalConnection {
                left: ModalAttachment {
                    component: 0,
                    shapes: vec![1.0],
                },
                right: ModalAttachment {
                    component: 1,
                    shapes: vec![1.0],
                },
                stiffness_n_m: 0.1 * omega * omega,
                damping_n_s_m: 0.0,
                rest_extension_m: 0.0,
            };
            let mut model = CoupledModalSystem::new(
                vec![
                    modal_model(omega, rate, z),
                    modal_model(omega, rate, [-z[0], -z[1]]),
                ],
                vec![connection],
                config,
                &CancelGate::new(),
            )
            .unwrap();
            model.step(&[0.0; 2]).unwrap();
            let x = model.components()[0].states()[0];
            normalized(
                &[x.displacement_m_sqrt_kg, x.velocity_m_sqrt_kg_per_s],
                omega,
            )
        }
        Owner::SeparateTest | Owner::Speaker => panic!("separate observation seam required"),
    }
}

fn speaker_phase(omega: f64, rate: u32) -> f64 {
    // For a lossless coil/cone with L=M=Bl=1, the nonzero continuous
    // eigenfrequency is sqrt(k+1). Positive but negligible losses/radiating
    // area isolate that pencil without changing the production propagator.
    let card = fs_couple::speaker::TsCard {
        re_ohm: 1e-12,
        le_h: 1.0,
        bl: 1.0,
        mms_kg: 1.0,
        cms_m_per_n: 1.0 / (omega * omega - 1.0),
        rms_n_s_m: 1e-12,
        sd_m2: 1e-12,
        x_max_m: 1.0,
    };
    let mut driver = fs_couple::speaker::TsDriver::new(card, None, None).unwrap();
    let dt = 1.0 / f64::from(rate);
    driver.step(0.01, dt).unwrap();
    let mut samples = Vec::new();
    for _ in 0..64 {
        let (record, _) = driver.step(0.0, dt).unwrap();
        samples.push(record.x[1]);
    }
    // Differentiate away the constant-flux equilibrium. The remaining real
    // sequence is a conjugate-pole pair; its recurrence measures cos(theta).
    let differences: Vec<_> = samples.windows(2).map(|x| x[1] - x[0]).collect();
    let numerator: f64 = differences.windows(3).map(|x| x[1] * (x[0] + x[2])).sum();
    let denominator: f64 = differences.windows(3).map(|x| 2.0 * x[1] * x[1]).sum();
    assert!(denominator > 0.0);
    let cosine = numerator / denominator;
    assert!((-1.0..1.0).contains(&cosine));
    det::acos(cosine)
}

fn phase(owner: Owner, omega: f64, rate: u32) -> f64 {
    if matches!(owner, Owner::Speaker) {
        return speaker_phase(omega, rate);
    }
    let a = advance(owner, omega, rate, [1.0, 0.0]);
    let b = advance(owner, omega, rate, [0.0, 1.0]);
    let re = 0.5 * (a[0] + b[1]);
    let determinant = a[0] * b[1] - b[0] * a[1];
    assert!(
        determinant > re * re,
        "linear map must have an oscillatory pair"
    );
    det::atan2(det::sqrt(determinant - re * re), re).abs()
}
fn cents(ratio: f64) -> f64 {
    1200.0 * det::ln(ratio) / core::f64::consts::LN_2
}
fn midpoint_cents(z: f64) -> f64 {
    cents(2.0 * det::atan(0.5 * z) / z)
}
fn metadata(output: std::io::Result<std::process::Output>) -> String {
    let output = output.expect("audit metadata command");
    assert!(output.status.success(), "audit metadata command failed");
    let value = String::from_utf8(output.stdout)
        .expect("UTF-8 audit metadata")
        .trim()
        .to_owned();
    assert!(
        !value.is_empty()
            && value
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b".-_".contains(&b)),
        "JSON-safe metadata"
    );
    value
}

#[test]
fn midpoint_warping_matches_atan_formula() {
    for rate in [48_000, 96_000, 192_000] {
        for frequency in [1_000.0, 5_000.0, 10_000.0] {
            let omega = core::f64::consts::TAU * frequency;
            let measured = cents(phase(Owner::Phs, omega, rate) * f64::from(rate) / omega);
            assert!((measured - midpoint_cents(omega / f64::from(rate))).abs() <= 1e-6);
        }
    }
    assert!((midpoint_cents(core::f64::consts::TAU * 5_000.0 / 48_000.0) + 59.11).abs() < 0.01);
    assert!((midpoint_cents(1.0) + 130.73).abs() < 0.1);
}

#[test]
fn available_render_linear_owners_are_measured_without_promoting_unmigrated_paths() {
    let git_sha = metadata(
        std::process::Command::new("git")
            .args(["rev-parse", "HEAD"])
            .output(),
    );
    let host = metadata(std::process::Command::new("uname").arg("-n").output());
    for row in PATHS {
        if matches!(row.owner, Owner::SeparateTest) {
            println!(
                "{{\"path\":\"{}\",\"integrator\":\"{:?}\",\"measured_cents\":null,\"verdict\":\"requires-owner-local-audit\",\"git_sha\":\"{}\",\"host\":\"{}\"}}",
                row.path, row.integrator, git_sha, host
            );
            continue;
        }
        for rate in [48_000, 96_000, 192_000] {
            for frequency in [1_000.0, 5_000.0, 10_000.0] {
                let omega = core::f64::consts::TAU * frequency;
                let z = omega / f64::from(rate);
                if row.bound.is_some_and(|bound| z >= bound) {
                    continue;
                }
                let exact_omega = if matches!(row.owner, Owner::Coupled) {
                    omega * det::sqrt(1.2)
                } else {
                    omega
                };
                let measured = cents(phase(row.owner, omega, rate) * f64::from(rate) / exact_omega);
                let predicted = if matches!(row.owner, Owner::Exact | Owner::Speaker) {
                    Some(0.0)
                } else if matches!(row.owner, Owner::Coupled) {
                    None
                } else {
                    Some(midpoint_cents(z))
                };
                if let Some(predicted) = predicted {
                    assert!(
                        (measured - predicted).abs() <= 1e-6,
                        "{}: {measured} vs {predicted}",
                        row.path
                    );
                }
                let verdict = if measured.abs() <= 0.01 {
                    "exact-linear-control"
                } else {
                    "requires-migration"
                };
                println!(
                    "{{\"path\":\"{}\",\"integrator\":\"{:?}\",\"omega_hz\":{},\"sample_rate_hz\":{},\"omega_dt\":{},\"measured_cents\":{},\"predicted_cents\":{},\"tolerance_cents\":0.000001,\"verdict\":\"{}\",\"git_sha\":\"{}\",\"host\":\"{}\"}}",
                    row.path,
                    row.integrator,
                    exact_omega / core::f64::consts::TAU,
                    rate,
                    exact_omega / f64::from(rate),
                    measured,
                    predicted.map_or("null".to_owned(), |x| x.to_string()),
                    verdict,
                    git_sha,
                    host
                );
            }
        }
    }
}

#[test]
fn impact_guard_cap_and_cap_plus_one() {
    let cap = 0.9 * core::f64::consts::PI;
    assert!(impact(1.0, cap * (1.0 - 1e-12), [0.0; 2]).is_ok());
    assert!(impact(1.0, cap, [0.0; 2]).is_err());
    assert!(
        midpoint_cents(cap) < -600.0,
        "existing guard admits a large disclosed phase defect"
    );
}

#[test]
fn coupled_modal_connection_warping_vs_exact_eigenproblem() {
    // Two equal unit masses: K=[[w²+k,-k],[-k,w²+k]], antisymmetric
    // eigenvalue w²+2k = 1.2w². The fixture excites exactly that eigenspace.
    // Exact component propagation plus trapezoidal connections is NOT a pure
    // trapezoid map of K, so the atan formula is not a valid oracle here.
    let omega = core::f64::consts::TAU * 5_000.0;
    let errors: Vec<_> = [48_000, 96_000, 192_000]
        .into_iter()
        .map(|rate| {
            cents(phase(Owner::Coupled, omega, rate) * f64::from(rate) / (omega * det::sqrt(1.2)))
                .abs()
        })
        .collect();
    assert!(errors[0] > 0.01, "connection defect must be live");
    assert!(errors[0] > errors[1] && errors[1] > errors[2]);
    for pair in errors.windows(2) {
        let order = det::ln(pair[0] / pair[1]) / core::f64::consts::LN_2;
        assert!(
            (order - 2.0).abs() < 0.15,
            "observed coupled phase order {order}"
        );
    }
}
