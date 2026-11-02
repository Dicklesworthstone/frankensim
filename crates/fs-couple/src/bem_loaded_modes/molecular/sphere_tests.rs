//! G1/G2: physical breathing mass, Maxwell storage and relaxing spherical sound.
use crate::bem_loaded_modes::{BemHarmonicOptions, BemLoadedModes, BemModesOptions};
use fs_bem::helmholtz::{
    Formulation, Medium, laplace::complex::ComplexOptions, molecular::MolecularExterior,
};
use fs_bem::observation_transfers::{MolecularResponseLimits, ReceiverSpec};
use fs_bem::panel3d::SpherePanels;
use fs_material::{
    gas::{GasSpec, GasState},
    relaxation::{RelaxingMedium, input::read_relaxation},
};
use fs_math::{c64::C64, det};
use fs_phs::RelaxationBranch;
use std::f64::consts::PI;

const OMEGA: f64 = 1000.0;
const RADIUS: f64 = 0.2;
const MASS: f64 = 0.1;
const STIFFNESS: f64 = 180_000.0;
const RESISTANCE: f64 = 6.0;
const EXCESS_STIFFNESS: f64 = 50_000.0;
const MEMORY_TIME: f64 = 0.0007;

fn exterior(relaxing: bool) -> MolecularExterior {
    // This intentionally slow caloric reservoir is synthetic, not a measured
    // atmospheric source. It distinguishes molecular loss from structural heat.
    let data = if relaxing {
        "declared 1 synthetic 1 vibration 1 100 1 power 101.325 293.15 0 200 500"
    } else {
        "declared 1 synthetic 0"
    };
    let law = read_relaxation(&mut data.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    MolecularExterior::new(law.evaluate(&medium).unwrap(), [800.0, 1200.0]).unwrap()
}

fn breathing(surface: SpherePanels, gas: &GasState) -> BemLoadedModes {
    let panels = surface.areas().len();
    let dimension = panels + 2;
    // q=sqrt(M)*x, hence v_skin=q_dot/sqrt(M) and f_q=F/sqrt(M).
    // A unit skin with an independently rescaled force would miss this chart.
    let weight = 1.0 / det::sqrt(MASS);
    BemLoadedModes::build(
        surface,
        vec![vec![weight; panels]],
        vec![STIFFNESS / MASS],
        vec![RESISTANCE / MASS],
        Medium::from(gas),
        BemModesOptions {
            formulation: Formulation::PlainCbie,
            boundary: ComplexOptions {
                max_panels: panels,
                max_bytes: 64 * panels * panels,
                max_work: 512 * (panels * panels) as u64,
            },
            panel_refinements: 0,
            reference_angular_frequency_rad_s: OMEGA,
            physical_residual_tolerance: 1e-10,
            max_dimension: dimension,
            max_dense_bytes: 256 * dimension * dimension + 1024 * dimension,
            max_operator_evaluations: 4,
            max_projection_work: 32 * (dimension * dimension) as u64,
        },
        &mut || false,
    )
    .unwrap()
    .with_relaxation_branches(
        &[RelaxationBranch {
            projection: vec![weight],
            stiffness: EXCESS_STIFFNESS,
            relaxation_time_s: MEMORY_TIME,
        }],
        1,
        || false,
    )
    .unwrap()
}

fn solve_options(panels: usize) -> BemHarmonicOptions {
    let n = (panels + 2) as u64;
    BemHarmonicOptions {
        max_load_cases: 2,
        max_response_bytes: 131_072,
        max_work: 64 * n.pow(3) + 256 * n.pow(2),
        power_relative_tolerance: 1e-10,
    }
}

fn response_limits(panels: usize) -> MolecularResponseLimits {
    MolecularResponseLimits {
        maximum_retained_entries: 2 * panels + 32,
        maximum_panel_products: 3 * panels as u64,
    }
}

fn receivers() -> [ReceiverSpec; 2] {
    [
        ReceiverSpec {
            id: "near radial".into(),
            position_m: [0.0, 0.0, 0.6],
        },
        ReceiverSpec {
            id: "far radial".into(),
            position_m: [0.0, 0.9, 0.0],
        },
    ]
}

fn exp_i(z: C64) -> C64 {
    C64::new(det::cos(z.re), det::sin(z.re)).scale(det::exp(-z.im))
}

fn close(a: C64, b: C64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * b.abs().max(1e-30),
        "{a:?} != {b:?}"
    );
}

#[test]
fn loaded_relaxing_sphere_converges_with_original_mass_memory_and_receiver_work() {
    let gas = exterior(true);
    let k = gas.source().relaxing_wavenumber(OMEGA).unwrap();
    let rho = gas.source().gas().density;
    let force = C64::new(0.3, -0.1);
    let s = C64::new(0.0, -OMEGA);
    let area = 4.0 * PI * RADIUS * RADIUS;
    // Independent spherical solution p(r)=p(a)*a/r*exp(ik(r-a)), with
    // dp/dr=i*rho*omega*v. No discrete BEM pressure or fit enters this load.
    let impedance = C64::new(0.0, rho * OMEGA) / C64::new(-1.0 / RADIUS - k.im, k.re);
    let exterior_load = impedance.scale(area);
    let memory_impedance =
        C64::from_re(EXCESS_STIFFNESS * MEMORY_TIME) / (C64::ONE + s.scale(MEMORY_TIME));
    let mechanical =
        C64::from_re(RESISTANCE) + s.scale(MASS) + C64::from_re(STIFFNESS) / s + memory_impedance;
    let velocity = force / (mechanical + exterior_load);
    let displacement = velocity / s;
    let elastic_memory = velocity.scale(det::sqrt(EXCESS_STIFFNESS) * MEMORY_TIME)
        / (C64::ONE + s.scale(MEMORY_TIME));
    let expected_power = [
        0.5 * (force * velocity.conj()).re,
        0.5 * RESISTANCE * velocity.norm_sq(),
        0.5 * elastic_memory.norm_sq() / MEMORY_TIME,
        0.5 * exterior_load.re * velocity.norm_sq(),
    ];
    assert!(expected_power.iter().all(|p| *p > 0.0));
    assert!(
        (expected_power[0] - expected_power[1..].iter().sum::<f64>()).abs()
            < 1e-14 * expected_power[0]
    );
    assert!(k.im > 0.01 * k.re);
    let microphones = receivers();
    let expected_fields = [0.6, 0.9].map(|r| {
        let pressure = impedance * velocity * exp_i(k.scale(r - RADIUS)).scale(RADIUS / r);
        let radial_velocity =
            pressure * C64::new(-1.0 / r - k.im, k.re) / C64::new(0.0, rho * OMEGA);
        (pressure, radial_velocity)
    });
    let mut errors = Vec::new();
    for level in 0..=2 {
        let surface = SpherePanels::icosphere(RADIUS, level).unwrap();
        let panels = surface.areas().len();
        let mechanics = breathing(surface, gas.source().gas());
        let prepared = mechanics
            .prepare_molecular_harmonic(
                &gas,
                OMEGA,
                &microphones,
                0.05,
                None,
                response_limits(panels),
                solve_options(panels),
                || false,
            )
            .unwrap();
        let modal_force = [force.scale(1.0 / det::sqrt(MASS))];
        let solution = prepared
            .harmonic_response_batch(&[&modal_force], || false)
            .unwrap()
            .remove(0);
        let actual_velocity = solution.modal_velocity[0].scale(1.0 / det::sqrt(MASS));
        let actual_displacement = solution.modal_displacement[0].scale(1.0 / det::sqrt(MASS));
        close(actual_velocity, s * actual_displacement, 2e-13);
        // Material history must follow the actual loaded motion, without any
        // fitted Q or replacement of the Maxwell arm by Kelvin resistance.
        let actual_memory = actual_velocity.scale(det::sqrt(EXCESS_STIFFNESS) * MEMORY_TIME)
            / (C64::ONE + s.scale(MEMORY_TIME));
        close(solution.elastic_memory_sqrt_j[0], actual_memory, 2e-12);
        let memory_heat = 0.5 * actual_memory.norm_sq() / MEMORY_TIME;
        assert!((solution.power.relaxation_w - memory_heat).abs() < 2e-12 * memory_heat);
        assert!(solution.mechanical_residual < 1e-10 && solution.memory_residual < 1e-10);
        assert!(solution.power.relative_balance < 1e-10);
        assert_eq!(solution.power.wall_w, 0.0);
        let actual_power = [
            solution.power.input_w,
            solution.power.kelvin_w,
            solution.power.relaxation_w,
            solution.power.exterior_w,
        ];
        let power_error = (actual_power
            .iter()
            .zip(expected_power)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            / expected_power.iter().map(|p| p * p).sum::<f64>())
        .sqrt();
        let (mut pressure_error, mut pressure_norm, mut velocity_error, mut velocity_norm) =
            (0.0, 0.0, 0.0, 0.0);
        for (i, &(p, u)) in expected_fields.iter().enumerate() {
            let field = solution.receivers[i];
            let radial = field.particle_velocity_m_s[if i == 0 { 2 } else { 1 }];
            pressure_error += (field.pressure_pa - p).norm_sq();
            pressure_norm += p.norm_sq();
            velocity_error += (radial - u).norm_sq();
            velocity_norm += u.norm_sq();
            if level == 2 {
                let r: f64 = if i == 0 { 0.6 } else { 0.9 };
                let flux = 2.0 * PI * r * r * (field.pressure_pa * radial.conj()).re;
                let exact_flux = expected_power[3] * det::exp(-2.0 * k.im * (r - RADIUS));
                assert!(
                    (flux / exact_flux - 1.0).abs() < 0.12,
                    "receiver flux includes molecular power lost before r={r}: {flux} vs {exact_flux}"
                );
            }
        }
        errors.push([
            (actual_velocity - velocity).abs() / velocity.abs(),
            (pressure_error / pressure_norm).sqrt(),
            (velocity_error / velocity_norm).sqrt(),
            power_error,
        ]);
        if level == 2 {
            close(actual_displacement, displacement, 0.08);
            close(solution.elastic_memory_sqrt_j[0], elastic_memory, 0.08);
            close(
                solution.fluid_reaction[0].scale(det::sqrt(MASS)),
                -exterior_load * velocity,
                0.08,
            );
        }
    }
    eprintln!("loaded molecular sphere [motion, pressure, velocity, work] errors: {errors:?}");
    for quantity in 0..4 {
        assert!(
            errors[1][quantity] < 0.8 * errors[0][quantity],
            "{errors:?}"
        );
        assert!(
            errors[2][quantity] < 0.8 * errors[1][quantity],
            "{errors:?}"
        );
        assert!(errors[2][quantity] < 0.08, "{errors:?}");
    }
}

#[test]
fn empty_molecular_reservoirs_recover_full_mechanical_pressure_memory_solve() {
    let gas = exterior(false);
    let surface = SpherePanels::icosphere(RADIUS, 0).unwrap();
    let panels = surface.areas().len();
    let mechanics = breathing(surface, gas.source().gas());
    let force = [C64::new(0.3, -0.1).scale(1.0 / det::sqrt(MASS))];
    let options = solve_options(panels);
    let original = mechanics
        .harmonic_response_batch(OMEGA, &[&force], &options, || false)
        .unwrap()
        .remove(0);
    let microphones = receivers();
    let prepared = mechanics
        .prepare_molecular_harmonic(
            &gas,
            OMEGA,
            &microphones,
            0.05,
            None,
            response_limits(panels),
            options,
            || false,
        )
        .unwrap();
    let molecular = prepared
        .harmonic_response_batch(&[&force], || false)
        .unwrap()
        .remove(0);
    close(
        molecular.modal_displacement[0],
        original.modal_displacement[0],
        2e-10,
    );
    close(
        molecular.modal_velocity[0],
        original.modal_velocity[0],
        2e-10,
    );
    close(
        molecular.elastic_memory_sqrt_j[0],
        original.elastic_memory_sqrt_j[0],
        2e-10,
    );
    for (a, b) in [
        (molecular.power.input_w, original.power.input_w),
        (molecular.power.kelvin_w, original.power.kelvin_w),
        (molecular.power.relaxation_w, original.power.relaxation_w),
        (molecular.power.exterior_w, original.power.radiated_w),
    ] {
        assert!((a - b).abs() < 2e-10 * b.abs(), "{a} != {b}");
    }
}
