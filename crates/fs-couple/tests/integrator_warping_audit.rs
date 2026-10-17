//! G1/G5 linear propagation audit. The closure test fails on unmeasured owners
//! and on undeclared non-exact bands. A measured integration tolerance alone
//! is not instrument fidelity. Network/piano companions keep their actual
//! owner-local observables, joined from successful native runs at this HEAD
//! and host through FS_IRP_NETWORK_AUDIT_LOG and FS_IRP_PIANO_AUDIT_LOG.

use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::dynamic::{
    ApertureDrive, ApertureState, DynamicAperture, DynamicApertureSpec,
};
use fs_couple::bernoulli_aperture::modal::{ApertureMode, ModalAperture, ModalApertureSpec};
use fs_couple::bernoulli_aperture::paired::{
    ApertureBlade, BladePropagator, BladeState, PairedAperture, PairedApertureLimits,
    PairedApertureSpec,
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
    PhsExact,
    Descriptor,
    Piston,
    Impact,
    Dynamic,
    Paired,
    PairedExact,
    Modal,
    Coupled,
    Uncoupled,
    Exact,
    Speaker,
    Network,
    Piano(u8),
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
        path: "acoustic_realize::realize_coupled_kc / resolved nonlinear step_render",
        integrator: Integrator::GonzalezMidpoint,
        bound: Some(fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP),
        owner: Owner::Phs,
    },
    AuditRow {
        path: "glottis / exact complete fold island",
        integrator: Integrator::ExactZoh,
        bound: None,
        owner: Owner::PhsExact,
    },
    AuditRow {
        path: "speaker::TsDriver / exact full quadratic step",
        integrator: Integrator::ExactZoh,
        bound: None,
        owner: Owner::Speaker,
    },
    AuditRow {
        path: "thin_plate::VkBody / resolved nonlinear step_render",
        integrator: Integrator::GonzalezMidpoint,
        bound: Some(fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP),
        owner: Owner::Phs,
    },
    AuditRow {
        path: "thin_plate::PlateBank cavity / exact held-flow island",
        integrator: Integrator::ExactZoh,
        bound: None,
        owner: Owner::PhsExact,
    },
    AuditRow {
        path: "acoustic_realize::realize_coupled_ode / step_chain_inlet / step_plate_chain",
        integrator: Integrator::ExactZoh,
        bound: None,
        owner: Owner::PhsExact,
    },
    AuditRow {
        path: "acoustic_realize::realize_dirac_join / step_descriptor",
        integrator: Integrator::GonzalezMidpoint,
        bound: Some(fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP),
        owner: Owner::Descriptor,
    },
    AuditRow {
        path: "acoustic_realize::realize_blown_ode and realize_reed_ode",
        integrator: Integrator::ExactZoh,
        bound: None,
        owner: Owner::PhsExact,
    },
    AuditRow {
        path: "brass_loop::lip / exact held-force island",
        integrator: Integrator::ExactZoh,
        bound: None,
        owner: Owner::PhsExact,
    },
    AuditRow {
        path: "cavity_phs::realize_plate_cavity / exact complete transformer",
        integrator: Integrator::ExactZoh,
        bound: None,
        owner: Owner::PhsExact,
    },
    AuditRow {
        path: "render::plate::ImpactSystem",
        integrator: Integrator::GonzalezMidpoint,
        bound: Some(fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP),
        owner: Owner::Impact,
    },
    AuditRow {
        path: "bernoulli_aperture::DynamicAperture",
        integrator: Integrator::GonzalezMidpoint,
        bound: Some(fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP),
        owner: Owner::Dynamic,
    },
    AuditRow {
        path: "bernoulli_aperture::PairedAperture",
        integrator: Integrator::GonzalezMidpoint,
        bound: Some(fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP),
        owner: Owner::Paired,
    },
    AuditRow {
        path: "bernoulli_aperture::PairedAperture (ExactZoh blades)",
        integrator: Integrator::ExactZoh,
        bound: None,
        owner: Owner::PairedExact,
    },
    AuditRow {
        path: "bernoulli_aperture::ModalAperture",
        integrator: Integrator::GonzalezMidpoint,
        bound: Some(fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP),
        owner: Owner::Modal,
    },
    AuditRow {
        path: "bernoulli_aperture::uniform (through ModalAperture)",
        integrator: Integrator::GonzalezMidpoint,
        bound: Some(fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP),
        owner: Owner::Modal,
    },
    AuditRow {
        path: "render::schedule::force::CoupledModalSystem connections",
        integrator: Integrator::ExactZoh,
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
        bound: Some(fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP),
        owner: Owner::Piston,
    },
    AuditRow {
        path: "cavity_phs::piston::thermal::boundary::network",
        integrator: Integrator::FirstOrder,
        bound: None,
        owner: Owner::Network,
    },
    AuditRow {
        path: "engine-piano::PianoAudio 1x",
        integrator: Integrator::Trapezoidal,
        bound: None,
        owner: Owner::Piano(1),
    },
    AuditRow {
        path: "engine-piano::PianoAudio 2x",
        integrator: Integrator::Trapezoidal,
        bound: None,
        owner: Owner::Piano(2),
    },
    AuditRow {
        path: "engine-piano::PianoAudio 4x",
        integrator: Integrator::Trapezoidal,
        bound: None,
        owner: Owner::Piano(4),
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
            // m=1/omega, k=omega gives the same physical pole, with q and
            // momentum p=m*v on equal scales. This prevents the reference
            // Newton solve's absolute residual from dominating a microcent
            // phase measurement through an ill-scaled (q,p) coordinate pair.
            let model = fs_phs::mass_spring_damper(1.0 / omega, omega, 0.0).unwrap();
            // Homogeneity permits unit-sized coordinates for this linear
            // control, so an absolute Newton budget resolves the phase too.
            let x = fs_phs::step(&model, &z, &[0.0], dt).unwrap().x;
            [x[0], x[1]]
        }
        Owner::PhsExact => {
            let model = fs_phs::mass_spring_damper(1.0, omega * omega, 0.0).unwrap();
            let split = fs_phs::linear_exact::split_quadratic(
                &model,
                fs_phs::linear_exact::LinearSplitBudget {
                    max_modes: 2,
                    max_omega_dt: core::f64::consts::PI,
                },
            )
            .unwrap();
            let x = fs_phs::linear_exact::step_render(&model, &split, &[q, v], &[0.0], dt)
                .unwrap()
                .x;
            normalized(&x, omega)
        }
        Owner::Descriptor => {
            let a = fs_phs::mass_spring_damper(1.0, omega * omega, 0.0).unwrap();
            let b = fs_phs::mass_spring_damper(1.0, omega * omega, 0.0).unwrap();
            let sys = fs_phs::join_port(a, b, 0, 0).unwrap();
            let rate = fs_phs::linear_exact::descriptor_pole_rate_bound(&sys).unwrap();
            fs_phs::linear_exact::admit_midpoint_band(rate, dt).unwrap();
            let record = fs_phs::step_descriptor(&sys, &[q, v, q, v, 0.0], &[], dt).unwrap();
            normalized(&record.x, omega)
        }
        Owner::Piston => piston_map::column(omega, dt, z),
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
        Owner::Dynamic | Owner::Paired | Owner::PairedExact => {
            let mass = 1e-5;
            let stiffness = mass * omega * omega;
            let area = 1e-6;
            // Manufacture zero differential pressure by cancelling the swept
            // volume with an externally prescribed body flow. This isolates
            // the actual mechanical map without approximating the jet law.
            // The exact blade sweeps its exact mean velocity dx/dt.
            let h = 0.5 * dt;
            let vm = if matches!(owner, Owner::PairedExact) {
                (q * (det::cos(omega * dt) - 1.0) + v * det::sin(omega * dt) / omega) / dt
            } else {
                (v - h * omega * omega * q) / (1.0 + h * h * omega * omega)
            };
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
                let mut model = PairedAperture::new_with_propagator(
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
                    if matches!(owner, Owner::PairedExact) {
                        BladePropagator::ExactZoh
                    } else {
                        BladePropagator::Midpoint
                    },
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
        Owner::Coupled | Owner::Uncoupled => {
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
                stiffness_n_m: if matches!(owner, Owner::Coupled) {
                    0.1 * omega * omega
                } else {
                    0.0
                },
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
        Owner::Network | Owner::Piano(_) | Owner::Speaker => {
            panic!("separate observation seam required")
        }
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

// Only the canonical flat JSONL emitted by the two named companion tests is
// consumed here. This is not a general JSON parser or a replacement test owner.
fn companion_field<'a>(line: &'a str, name: &str) -> &'a str {
    let key = format!("\"{name}\":");
    let mut fields = line.match_indices(&key);
    let (offset, _) = fields.next().expect("required companion field");
    assert!(fields.next().is_none(), "duplicate companion field {name}");
    line[offset + key.len()..]
        .split([',', '}'])
        .next()
        .unwrap()
        .trim()
}
fn companion_text<'a>(line: &'a str, name: &str) -> &'a str {
    let value = companion_field(line, name)
        .strip_prefix('"')
        .and_then(|v| v.strip_suffix('"'))
        .expect("quoted companion field");
    assert!(!value.contains(['"', '\\']), "canonical unescaped field");
    value
}
fn companion_number(line: &str, name: &str) -> f64 {
    let value: f64 = companion_field(line, name)
        .parse()
        .expect("numeric companion field");
    assert!(value.is_finite(), "finite companion field {name}");
    value
}
fn companion_rows(variable: &str, path: &str, git: &str, host: &str) -> Vec<String> {
    let file = std::env::var_os(variable).unwrap_or_else(|| {
        panic!("run the real companion test with --nocapture and set {variable}")
    });
    let log = std::fs::read_to_string(file).expect("read native companion test output");
    let outcomes: Vec<_> = log
        .lines()
        .filter(|l| l.starts_with("test result:"))
        .collect();
    assert!(
        !outcomes.is_empty() && outcomes.iter().all(|l| l.starts_with("test result: ok.")),
        "{variable} must come from a successful test run"
    );
    let rows: Vec<_> = log
        .lines()
        .filter_map(|line| {
            // libtest can prefix the first --nocapture line with the test name.
            let start = line.find("{\"path\":")?;
            let row = &line[start..];
            (companion_text(row, "path") == path).then(|| {
                assert_eq!(companion_text(row, "git_sha"), git, "stale companion HEAD");
                assert_eq!(
                    companion_text(row, "host"),
                    host,
                    "different companion host"
                );
                row.to_owned()
            })
        })
        .collect();
    assert_eq!(
        rows.len(),
        3,
        "all three companion rungs/factors required exactly once"
    );
    rows
}

fn measured_companion(row: &AuditRow, git: &str, host: &str) {
    match row.owner {
        Owner::Network => {
            let mut rows = companion_rows("FS_IRP_NETWORK_AUDIT_LOG", row.path, git, host);
            rows.sort_by(|a, b| {
                companion_number(a, "steps").total_cmp(&companion_number(b, "steps"))
            });
            let mut errors = Vec::new();
            for (i, (record, steps)) in rows.iter().zip([160.0, 320.0, 640.0]).enumerate() {
                assert_eq!(companion_text(record, "integrator"), "FirstOrder");
                assert_eq!(companion_number(record, "steps"), steps);
                assert_eq!(companion_number(record, "next_steps"), 2.0 * steps);
                assert!((companion_number(record, "dt_s") * steps - 0.02).abs() < 1e-15);
                let temperature = companion_number(record, "temperature_k");
                let next = companion_number(record, "next_temperature_k");
                if i != 0 {
                    assert_eq!(
                        companion_number(&rows[i - 1], "next_temperature_k"),
                        temperature
                    );
                }
                let error = (temperature - next).abs();
                assert!(error > 1e-6, "refinement must exceed solver noise");
                assert!((error - companion_number(record, "next_rung_difference_k")).abs() < 1e-12);
                errors.push(error);
                println!("{record}");
            }
            for pair in errors.windows(2) {
                let order = det::ln(pair[0] / pair[1]) / core::f64::consts::LN_2;
                assert!(
                    (order - 1.0).abs() <= 0.1,
                    "actual first-order companion: {order}"
                );
            }
        }
        Owner::Piano(factor) => {
            let rows = companion_rows(
                "FS_IRP_PIANO_AUDIT_LOG",
                "engine-piano::PianoAudio",
                git,
                host,
            );
            let factors: std::collections::BTreeSet<_> =
                rows.iter().map(|r| companion_field(r, "factor")).collect();
            assert_eq!(
                factors,
                ["1", "2", "4"]
                    .into_iter()
                    .collect::<std::collections::BTreeSet<_>>()
            );
            let record = rows
                .iter()
                .find(|r| companion_number(r, "factor") == f64::from(factor))
                .unwrap();
            assert_eq!(
                companion_text(record, "integrator"),
                "exact local ZOH plus finite-step coupled solve"
            );
            let rate = 48_000.0 * f64::from(factor);
            assert_eq!(companion_number(record, "sample_rate_hz"), rate);
            assert_eq!(companion_number(record, "tolerance_cents"), 1.0);
            let omega = core::f64::consts::TAU * companion_number(record, "omega_hz");
            assert!(omega > 0.0);
            assert!((companion_number(record, "omega_dt") - omega / rate).abs() < 1e-12);
            let a = companion_number(record, "window_a_omega_rad_s");
            let b = companion_number(record, "window_b_omega_rad_s");
            assert!(a > 0.0 && b > 0.0);
            let measured = cents(0.5 * (a + b) / omega);
            let discrepancy = cents(a / b).abs();
            assert!((measured - companion_number(record, "measured_cents")).abs() < 1e-6);
            assert!(
                (discrepancy - companion_number(record, "window_difference_cents")).abs() < 1e-6
            );
            assert!(
                discrepancy < 1.0 && measured.abs() + 0.5 * discrepancy <= 1.0,
                "{}: {measured} cents, window discrepancy {discrepancy}; migration required",
                row.path
            );
            // This is the bead's actual PianoAudio coupled-eigenproblem
            // observable, not an exact-map claim or a fabricated global guard.
            println!("{record}");
        }
        _ => panic!("not a companion owner"),
    }
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
        if matches!(row.owner, Owner::Network | Owner::Piano(_)) {
            measured_companion(row, &git_sha, &host);
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
                let predicted = if matches!(
                    row.owner,
                    Owner::Exact
                        | Owner::Speaker
                        | Owner::Coupled
                        | Owner::PhsExact
                        | Owner::PairedExact
                ) {
                    Some(0.0)
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
    let cap = fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP;
    assert!(impact(1.0, cap, [0.0; 2]).is_ok());
    assert!(matches!(
        impact(1.0, cap + 1e-12, [0.0; 2]),
        Err(fs_couple::render::plate::impact::ImpactError::TimeResolution(_))
    ));
    assert!(midpoint_cents(cap).abs() <= 1.0);
    assert!(matches!(
        impact(1.0, 0.9 * core::f64::consts::PI * (1.0 - 1e-12), [0.0; 2]),
        Err(fs_couple::render::plate::impact::ImpactError::TimeResolution(_))
    ));
}

#[test]
fn every_render_path_is_listed_with_integrator_and_measured_warping() {
    available_render_linear_owners_are_measured_without_promoting_unmigrated_paths();
    impact_guard_cap_and_cap_plus_one();
    valve_angular_step_cap_and_cap_plus_one();
    dynamic_aperture_refuses_beyond_one_cent_resolved_band();
    paired_aperture_refuses_beyond_one_cent_resolved_band();
    modal_aperture_includes_full_viscosity_in_the_one_cent_band();
    let mut names = std::collections::BTreeSet::new();
    let mut unresolved = Vec::new();
    for row in PATHS {
        assert!(
            !row.path.is_empty() && names.insert(row.path),
            "render paths must be named once"
        );
        let exact = matches!(row.integrator, Integrator::ExactZoh);
        // The FirstOrder row is measured by a four-solution order ladder;
        // PianoAudio by its public coupled-board observations at all factors.
        // Both have already executed mandatory companion checks above. Neither
        // may substitute a scalar midpoint control for its actual owner.
        let resolved = match row.owner {
            Owner::Network => matches!(row.integrator, Integrator::FirstOrder),
            Owner::Piano(_) => matches!(row.integrator, Integrator::Trapezoidal),
            _ => {
                exact
                    || row.bound.is_some_and(|z| {
                        z > 0.0 && z <= fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP
                    })
            }
        };
        if !resolved {
            unresolved.push(row.path);
        }
    }
    // Disclosure rows cannot silently become acceptance rows merely because
    // all the already-available propagators pass their numerical controls.
    assert!(
        unresolved.is_empty(),
        "unmeasured or unresolved render paths: {unresolved:?}"
    );
}

#[test]
fn valve_angular_step_cap_and_cap_plus_one() {
    use fs_couple::bernoulli_aperture::modal::material::admit_modal_angular_step;
    assert!(admit_modal_angular_step(1.0, 1.0 / 192_000.0).is_ok());
    assert!(matches!(
        admit_modal_angular_step(1.0 + 1e-12, 1.0 / 192_000.0),
        Err(fs_couple::acoustic_realize::AcousticRealizeError::InvalidDescription { .. })
    ));
    for (cap, dt) in [
        (0.0, 1e-6),
        (f64::NAN, 1e-6),
        (1.0, 0.0),
        (1.0, f64::INFINITY),
    ] {
        assert!(admit_modal_angular_step(cap, dt).is_err());
    }
}

#[test]
fn dynamic_aperture_refuses_beyond_one_cent_resolved_band() {
    let cap = fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP;
    let make = |dt| {
        DynamicAperture::new(
            DynamicApertureSpec {
                aperture: BernoulliAperture {
                    rest_opening_m: REST,
                    width_m: 0.01,
                    closing_pressure_pa: 1e6,
                },
                mass_kg: 1.0,
                stiffness_n_m: 1.0,
                damping_ratio: 0.0,
                density_kg_m3: 1.2,
                impedance_pa_s_m3: 1e6,
                time_step_s: dt,
                max_steps: 1,
            },
            ApertureState {
                opening_m: REST,
                opening_velocity_m_s: 0.0,
            },
            obstacle(),
        )
    };
    assert!(make(cap / (1.0 + 1e-12)).is_ok());
    assert!(matches!(
        make(cap + 1e-12),
        Err(
            fs_couple::acoustic_realize::AcousticRealizeError::TimeResolution(
                fs_phs::linear_exact::ResolvedBandError::Outside { .. }
            )
        )
    ));
}

#[test]
fn paired_aperture_refuses_beyond_one_cent_resolved_band() {
    let cap = fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP;
    let blade = ApertureBlade {
        mass_kg: 1.0,
        stiffness_n_m: 1.0,
        damping_n_s_m: 0.0,
        pressure_area_m2: 1e-9,
    };
    let make = |dt| {
        PairedAperture::new(
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
            [BladeState::default(); 2],
            obstacle(),
        )
    };
    assert!(make(cap / (1.0 + 2e-12)).is_ok());
    assert!(matches!(
        make(cap + 1e-12),
        Err(
            fs_couple::acoustic_realize::AcousticRealizeError::TimeResolution(
                fs_phs::linear_exact::ResolvedBandError::Outside { .. }
            )
        )
    ));
}

#[test]
fn modal_aperture_includes_full_viscosity_in_the_one_cent_band() {
    use fs_couple::acoustic_realize::AcousticRealizeError;
    let cap = fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP;
    let make = |dt| {
        ModalAperture::new(
            ModalApertureSpec {
                rest_opening_m: REST,
                width_m: 0.01,
                modes: vec![
                    ApertureMode {
                        mass_kg: 1.0,
                        stiffness_n_m: 1.0,
                        damping_n_s_m: 0.0,
                        pressure_area_m2: 0.0,
                        gap_projection: 1.0,
                        force_projection: [0.0; 2]
                    };
                    2
                ],
                density_kg_m3: 1.2,
                impedance_pa_s_m3: 1e6,
                time_step_s: dt,
                max_steps: 1,
                limits: limits(),
            },
            vec![BladeState::default(); 2],
            obstacle(),
        )
    };
    assert!(make(cap).is_ok());
    assert!(matches!(
        make(cap + 1e-12),
        Err(AcousticRealizeError::TimeResolution(_))
    ));
    // A full off-diagonal C changes the complete pole norm although each
    // bare mechanical frequency is still one radian per second.
    assert!(
        make(cap / 2.5)
            .unwrap()
            .with_viscosity(vec![1.0, 0.5, 0.5, 1.0], &CancelGate::new())
            .is_ok()
    );
    assert!(matches!(
        make(cap / 2.5 + 1e-12)
            .unwrap()
            .with_viscosity(vec![1.0, 0.5, 0.5, 1.0], &CancelGate::new()),
        Err(AcousticRealizeError::TimeResolution(_))
    ));
}

#[test]
fn coupled_modal_connection_warping_vs_exact_eigenproblem() {
    // Two equal unit masses: K=[[w²+k,-k],[-k,w²+k]], antisymmetric
    // eigenvalue w²+2k = 1.2w². The fixture excites exactly that eigenspace.
    // The full pencil exponential includes this eigenvalue; exact component
    // propagation alone cannot reproduce it.
    let omega = core::f64::consts::TAU * 5_000.0;
    let errors: Vec<_> = [48_000, 96_000, 192_000]
        .into_iter()
        .map(|rate| {
            cents(phase(Owner::Coupled, omega, rate) * f64::from(rate) / (omega * det::sqrt(1.2)))
                .abs()
        })
        .collect();
    assert!(
        errors.iter().all(|error| *error <= 0.01),
        "coupled phase errors in cents: {errors:?}"
    );
    // Run the real network again with the connection stiffness deleted. Its
    // pole must revert to omega: relative to the coupled oracle, the signed
    // error is -600*log2(1.2) cents, far beyond 100 times the accepted band.
    let dropped =
        cents(phase(Owner::Uncoupled, omega, 48_000) * 48_000.0 / (omega * det::sqrt(1.2)));
    let missing_connection = cents(1.0 / det::sqrt(1.2));
    assert!(
        (dropped - missing_connection).abs() <= 0.01,
        "component-only phase {dropped} cents vs {missing_connection} cents"
    );
    assert!(
        dropped.abs() > 100.0 * 0.01,
        "deleting the coupling must fail the 0.01-cent band: {dropped} cents"
    );
}

mod piston_map {
    use fs_couple::cavity_phs::piston::{GasPiston, PistonInput, PistonSolve};
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId,
        LawParameter, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
        Provenance, QueryPoint, UncertaintyModel,
    };
    use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
    use fs_material::thermochemical_gas::energy::GasEnergySolve;
    use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
    use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
    use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    use fs_qty::{DynViscosity, MolarMass, Pressure, QuantitySpec, Temperature};
    use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
    use fs_thermochem::{
        Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1,
        NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1,
        Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId,
        StandardStateConventionV1, StandardStatePhaseV1, required_nasa9_parameter_dims,
    };

    fn quantity(kind: QuantityKind) -> QuantitySpec {
        QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
    }
    fn provenance() -> Provenance {
        Provenance {
            source: "authored constant-cp analytic phase fixture".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        }
    }
    fn parcel() -> ClosedGasParcel {
        let mut coefficients = [0.0; 9];
        coefficients[2] = 3.5;
        let mut parameters = std::collections::BTreeMap::new();
        for ((key, dims), value) in required_nasa9_parameter_dims()
            .into_iter()
            .take(9)
            .zip(coefficients)
        {
            parameters.insert(key.to_owned(), LawParameter { value, dims });
        }
        parameters.insert(
            "reference_pressure".into(),
            LawParameter {
                value: 100_000.0,
                dims: Pressure::DIMS,
            },
        );
        let card = ConstitutiveModelCard {
            law: LawId(NASA9_LAW_ID_V1.into()),
            law_version: NASA9_LAW_VERSION_V1,
            parameters,
            state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
            sources: vec![],
            provenance: provenance(),
        };
        let species = Nasa9StandardStateModelV1::new(
            SpeciesId::new("phase-gas").unwrap(),
            MolarMass::new(0.028),
            StandardStateConventionV1::new(
                StandardStatePhaseV1::Gas,
                ReferenceEquationOfStateV1::IdealGas,
                Pressure::new(100_000.0),
                ElementalReferenceIdV1::new("analytic").unwrap(),
            )
            .unwrap(),
            vec![Nasa9RegionV1::from_card(card).unwrap()],
        )
        .unwrap();
        let mut claims = ClaimSet::new();
        for (key, dims, value) in [
            ("dynamic_viscosity", DynViscosity::DIMS, 2e-5),
            ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.03),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(key, dims),
                    value: PropertyValue::Curve {
                        abscissa: "temperature".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(200.0, value), (2000.0, value)],
                        dims,
                    },
                    validity: ValidityDomain::unconstrained()
                        .with_quantity(
                            "temperature",
                            quantity(QuantityKind::AbsoluteTemperature),
                            200.0,
                            2000.0,
                        )
                        .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0, 1e7),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::LinearInside,
                    observations: vec![],
                    provenance: provenance(),
                })
                .unwrap();
        }
        let transport = SpeciesTransportMaterial {
            species: SpeciesId::new("phase-gas").unwrap(),
            material: MaterialCard::assemble(
                MaterialStateId {
                    chemistry: "phase-gas".into(),
                    phase: "gas".into(),
                    process: "analytic".into(),
                    revision: 0,
                },
                claims,
                vec![],
            )
            .unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        };
        let model = ThermochemicalGasModel::new(
            FrozenIdealGasMixtureModelV1::new(
                vec![species],
                Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
            )
            .unwrap(),
            vec![transport],
            None,
        )
        .unwrap();
        let point = QueryPoint::new()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                400.0,
            )
            .unwrap()
            .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0)
            .unwrap();
        ClosedGasParcel::new(model, point, 0.001).unwrap()
    }
    fn advance(omega: f64, dt: f64, z: [f64; 2], amplitude: f64) -> [f64; 2] {
        let mut gas = parcel();
        let volume = gas.volume_m3();
        let area = 0.01;
        let medium = gas.state().state();
        let pressure = medium.pressure;
        let mass =
            medium.density * medium.sound_speed.powi(2) * area * area / (volume * omega * omega);
        if z[0] != 0.0 {
            gas.set_volume_isentropic(
                volume + area * amplitude * z[0],
                GasEntropySolve {
                    minimum_temperature_k: 200.0,
                    maximum_temperature_k: 1800.0,
                    entropy_tolerance_j_kg_k: 1e-11,
                    maximum_iterations: 128,
                },
                || false,
            )
            .unwrap();
        }
        let mut model = GasPiston::new(gas, mass, area, omega * amplitude * z[1]).unwrap();
        // Allocate material error in J, then convert to the source's J/kg.
        // Midpoint reserves another factor four for its two inversions.
        let gas_mass = model.gas().mass_kg();
        let energy_budget_j = 1e-10;
        let record = model
            .step(
                PistonInput {
                    duration_s: dt,
                    heat_j: 0.0,
                    exterior_pressure_pa: pressure,
                    external_force_n: 0.0,
                },
                PistonSolve {
                    volume_bracket_m3: [0.0007, 0.0015],
                    gas: GasEnergySolve {
                        minimum_temperature_k: 200.0,
                        maximum_temperature_k: 1800.0,
                        energy_tolerance_j_kg: energy_budget_j / (8.0 * gas_mass),
                        maximum_iterations: 128,
                    },
                    momentum_tolerance_n_s: 1e-14,
                    energy_tolerance_j: energy_budget_j,
                    maximum_iterations: 128,
                },
                || false,
            )
            .unwrap();
        assert!(record.total_energy_residual_j.abs() <= energy_budget_j);
        assert!(record.gas.energy_residual_j.abs() <= energy_budget_j / 8.0);
        [
            (model.gas().volume_m3() - volume) / area,
            model.velocity_m_s() / omega,
        ]
    }
    pub(super) fn column(omega: f64, dt: f64, z: [f64; 2]) -> [f64; 2] {
        // Symmetric perturbations plus Richardson remove the leading A²
        // source nonlinearity; the production gas/momentum owner still solves
        // every trajectory. No scalar midpoint surrogate replaces the piston.
        let maps: Vec<_> = [1e-4, 5e-5, 2.5e-5]
            .into_iter()
            .map(|a| {
                let plus = advance(omega, dt, z, a);
                let minus = advance(omega, dt, z, -a);
                [
                    (plus[0] - minus[0]) / (2.0 * a),
                    (plus[1] - minus[1]) / (2.0 * a),
                ]
            })
            .collect();
        let coarse = [
            (4.0 * maps[1][0] - maps[0][0]) / 3.0,
            (4.0 * maps[1][1] - maps[0][1]) / 3.0,
        ];
        let fine = [
            (4.0 * maps[2][0] - maps[1][0]) / 3.0,
            (4.0 * maps[2][1] - maps[1][1]) / 3.0,
        ];
        for (a, b) in coarse.iter().zip(fine) {
            assert!(
                (a - b).abs() < 1e-8,
                "piston linearization unresolved: {a} {b}"
            );
        }
        fine
    }
}
