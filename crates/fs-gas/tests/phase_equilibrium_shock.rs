//! Original phase normal shocks: independent PR and ideal-gas jump references.
#![cfg(feature = "material")]

#[allow(dead_code)]
#[path = "../../fs-material/tests/common/nonideal_flash.rs"]
mod common;
use common::{MIXTURE_MASS, model, source};
use fs_gas::material::phase_shock::{
    MixtureNormalShock, MixtureShockControl as Control, MixtureShockError as Error,
    equilibrium_mixture_normal_shock,
};
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_material::thermochemical_gas::nonideal::flash::ThermochemicalFlashError;
use fs_qty::{Pressure, Temperature, Velocity};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}; difference {}",
        actual - expected
    );
}
fn controls() -> Control {
    Control::new(
        [Temperature::new(200.0), Temperature::new(300.0)],
        [Pressure::new(5e6), Pressure::new(5.6e6)],
    )
}
fn shock(
    source: &ThermochemicalPengRobinsonModel,
    speed: f64,
    control: Control,
    cancelled: impl FnMut() -> bool,
) -> Result<MixtureNormalShock, Error> {
    equilibrium_mixture_normal_shock(
        source,
        Temperature::new(230.0),
        Pressure::new(3e6),
        Velocity::new(speed),
        control,
        cancelled,
    )
}
fn conservation(jump: &MixtureNormalShock, control: Control) {
    let a = &jump.upstream;
    let b = jump.downstream.state();
    let u1 = jump.upstream_speed.value();
    let u2 = jump.downstream_speed.value();
    let dp = b.pressure_pa() - a.pressure_pa();
    let dh = jump.downstream.accepted_increment_j_kg().unwrap();
    let mass = b.density_kg_m3() * u2 - a.density_kg_m3() * u1;
    let momentum = dp + b.density_kg_m3() * u2 * u2 - a.density_kg_m3() * u1 * u1;
    let energy = dh + 0.5 * (u2 * u2 - u1 * u1);
    let hugoniot = dh - 0.5 * dp * (1.0 / a.density_kg_m3() + 1.0 / b.density_kg_m3());
    let mass_tol = control.absolute_mass_flux_tolerance_kg_m2_s
        + control.relative_tolerance * jump.mass_flux_kg_m2_s;
    let momentum_tol = control.absolute_momentum_tolerance_pa + control.relative_tolerance * dp;
    let energy_tol = control.absolute_energy_tolerance_j_kg + control.relative_tolerance * dh.abs();
    close(mass, 0.0, mass_tol);
    close(momentum, 0.0, momentum_tol + 1e-8);
    close(energy, 0.0, energy_tol + 1e-9);
    close(hugoniot, 0.0, energy_tol);
    close(jump.mass_flux_residual_kg_m2_s, mass, 1e-10);
    close(jump.momentum_flux_residual_pa, momentum, 1e-8);
    close(jump.stagnation_enthalpy_residual_j_kg, energy, 1e-9);
    close(jump.hugoniot_residual_j_kg, hugoniot, 1e-9);
    assert!(jump.upstream_mach > 1.0 && jump.downstream_mach < 1.0);
    assert!(jump.entropy_change_j_kg_k >= -control.entropy_tolerance_j_kg_k);
    assert_eq!(jump.upstream_acoustics.source(), a);
    assert_eq!(jump.downstream_acoustics.source(), b);
    assert!(jump.pt_evaluations <= control.maximum_pt_evaluations);
    assert!(jump.pressure_evaluations <= control.maximum_pressure_evaluations);
}

#[test]
fn g1_equilibrium_shock_redistributes_latent_inventory_and_closes_original_energy() {
    let original = source();
    let c = controls();
    let jump = shock(&original, 225.0, c, || false).unwrap();
    // Independent direct PR fugacity + Rayleigh + Hugoniot solve; sound speed
    // comes from a separately solved isentropic density derivative.
    let b = jump.downstream.state();
    close(b.pressure_pa(), 5_335_710.658_109_467, 0.02);
    close(b.temperature_k(), 249.715_050_208_738_66, 2e-6);
    close(b.density_kg_m3(), 184.671_169_201_781_7, 2e-5);
    close(jump.downstream_speed.value(), 109.618_035_319_891_62, 2e-5);
    close(jump.enthalpy_increase_j_kg, 19_304.443_166_304_97, 1e-3);
    close(jump.entropy_change_j_kg_k, 5.164_964_948_157_831, 2e-5);
    close(
        jump.upstream_acoustics.state().sound_speed_m_s(),
        169.158_150_07,
        0.01,
    );
    close(
        jump.downstream_acoustics.state().sound_speed_m_s(),
        144.772_478_39,
        0.01,
    );
    close(
        jump.upstream.vapor_mass_fraction().unwrap(),
        0.408_352_234_634_033_25,
        2e-7,
    );
    close(
        b.vapor_mass_fraction().unwrap(),
        0.311_245_469_569_273_1,
        2e-7,
    );
    assert!(b.phase().two_phase().is_some());
    conservation(&jump, c);
}

#[test]
fn g2_equilibrium_shock_can_cross_from_two_phases_to_a_homogeneous_dense_outlet() {
    let original = source();
    let mut c = controls();
    c.pressure_bracket = [Pressure::new(8e6), Pressure::new(8.5e6)];
    let jump = shock(&original, 300.0, c, || false).unwrap();
    let b = jump.downstream.state();
    assert!(jump.upstream.phase().two_phase().is_some());
    assert!(b.phase().two_phase().is_none());
    close(b.pressure_pa(), 8_251_041.743_996_359, 0.03);
    close(b.temperature_k(), 265.501_521_362_098_6, 2e-6);
    close(b.density_kg_m3(), 255.954_581_944_758_3, 2e-5);
    close(jump.downstream_speed.value(), 105.452_514_761_626_05, 2e-5);
    close(jump.enthalpy_increase_j_kg, 39_439.883_565_224_45, 1e-3);
    close(jump.entropy_change_j_kg_k, 34.878_320_839_234_22, 2e-5);
    close(
        jump.downstream_acoustics.state().sound_speed_m_s(),
        410.917_605_87,
        0.01,
    );
    conservation(&jump, c);
}

#[test]
fn g3_dilute_limit_recovers_normal_shock_relations_without_a_gamma_input() {
    let original = model([0.0; 2], [0.0; 2], [100.0, 1000.0]);
    let gamma = 4.0 / 3.0;
    let gas_constant = R / MIXTURE_MASS;
    let mach_squared = 4.0;
    let input_speed = 2.0 * (gamma * gas_constant * 400.0).sqrt();
    let pressure_ratio = (2.0 * gamma * mach_squared - (gamma - 1.0)) / (gamma + 1.0);
    let density_ratio = (gamma + 1.0) * mach_squared / ((gamma - 1.0) * mach_squared + 2.0);
    let c = Control::new(
        [Temperature::new(300.0), Temperature::new(800.0)],
        [Pressure::new(400.0), Pressure::new(480.0)],
    );
    let jump = equilibrium_mixture_normal_shock(
        &original,
        Temperature::new(400.0),
        Pressure::new(100.0),
        Velocity::new(input_speed),
        c,
        || false,
    )
    .unwrap();
    let b = jump.downstream.state();
    // Finite-pressure PR corrections are independently about .00047 Pa and
    // .00049 K; the runtime uses the native EOS throughout.
    close(b.pressure_pa(), 100.0 * pressure_ratio, 0.001);
    close(
        b.temperature_k(),
        400.0 * pressure_ratio / density_ratio,
        0.001,
    );
    close(
        jump.downstream_speed.value(),
        input_speed / density_ratio,
        0.001,
    );
    close(jump.enthalpy_increase_j_kg, 335_601.255_184_979_9, 2e-4);
    close(jump.entropy_change_j_kg_k, 124.687_026_885_058, 2e-6);
    conservation(&jump, c);
}

#[test]
fn g4_shock_ignores_arbitrary_shared_species_energy_and_entropy_datums() {
    let plain = source();
    let shifted = model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]);
    let a = shock(&plain, 225.0, controls(), || false).unwrap();
    let b = shock(&shifted, 225.0, controls(), || false).unwrap();
    assert_eq!(
        a.downstream.state().temperature_k(),
        b.downstream.state().temperature_k()
    );
    assert_eq!(
        a.downstream.state().pressure_pa(),
        b.downstream.state().pressure_pa()
    );
    assert_eq!(
        a.downstream.state().density_kg_m3(),
        b.downstream.state().density_kg_m3()
    );
    assert_eq!(a.downstream_speed, b.downstream_speed);
    assert_eq!(a.enthalpy_increase_j_kg, b.enthalpy_increase_j_kg);
    assert_eq!(a.entropy_change_j_kg_k, b.entropy_change_j_kg_k);
    assert_eq!(a.upstream_mach, b.upstream_mach);
    assert_eq!(a.downstream_mach, b.downstream_mach);
    assert_eq!(a.pt_evaluations, b.pt_evaluations);
    assert_eq!(a.pressure_evaluations, b.pressure_evaluations);
    conservation(&b, controls());
}

#[test]
fn g5_identity_subsonic_and_unphysical_rayleigh_brackets_refuse() {
    let original = source();
    for speed in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            shock(&original, speed, controls(), || false),
            Err(Error::Invalid(_))
        ));
    }
    assert!(matches!(
        shock(&original, 100.0, controls(), || false),
        Err(Error::Admissibility {
            downstream_mach: None,
            ..
        }),
    ));
    let mut c = controls();
    c.pressure_bracket = [Pressure::new(3e6), Pressure::new(5.6e6)];
    assert!(matches!(
        shock(&original, 225.0, c, || false),
        Err(Error::Invalid(_))
    ));
    c = controls();
    c.pressure_bracket[1] = Pressure::new(9e6);
    assert!(matches!(
        shock(&original, 225.0, c, || false),
        Err(Error::Invalid(_))
    ));
    c = controls();
    c.acoustics.phase.stability.tangent_plane_tolerance *= 10.0;
    assert!(matches!(
        shock(&original, 225.0, c, || false),
        Err(Error::Invalid(_))
    ));
}

#[test]
fn g6_pressure_brackets_and_original_temperature_support_are_not_extrapolated() {
    let original = source();
    let mut c = controls();
    c.pressure_bracket = [Pressure::new(4e6), Pressure::new(4.2e6)];
    assert!(matches!(
        shock(&original, 225.0, c, || false),
        Err(Error::OutsideBracket { .. })
    ));
    let limited = model([0.0; 2], [0.0; 2], [200.0, 260.0]);
    assert!(matches!(
        shock(&limited, 225.0, controls(), || false),
        Err(Error::Equilibrium {
            source: ThermochemicalFlashError::Material(_),
            ..
        }),
    ));
}

#[test]
fn g7_pressure_and_pt_budgets_cover_preparation_ph_and_both_acoustic_stencils() {
    let original = source();
    let accepted = shock(&original, 225.0, controls(), || false).unwrap();
    let mut c = controls();
    c.maximum_pressure_evaluations = 2;
    assert!(matches!(
        shock(&original, 225.0, c, || false),
        Err(Error::WorkLimit {
            pressure_evaluations: 2,
            ..
        }),
    ));
    c = controls();
    c.maximum_pt_evaluations = 18;
    assert_eq!(
        shock(&original, 225.0, c, || false).unwrap_err(),
        Error::WorkLimit {
            pressure_evaluations: 0,
            pt_evaluations: 0
        },
    );
    c.maximum_pt_evaluations = 19;
    assert_eq!(
        shock(&original, 225.0, c, || false).unwrap_err(),
        Error::WorkLimit {
            pressure_evaluations: 0,
            pt_evaluations: 10
        },
    );
    c.maximum_pt_evaluations = accepted.pt_evaluations;
    assert_eq!(shock(&original, 225.0, c, || false).unwrap(), accepted);
    c.maximum_pt_evaluations -= 1;
    assert!(matches!(
        shock(&original, 225.0, c, || false),
        Err(Error::WorkLimit { .. })
            | Err(Error::Equilibrium {
                source: ThermochemicalFlashError::WorkLimit { .. },
                ..
            }),
    ));
}

#[test]
fn g8_final_cancellation_publishes_no_jump_and_does_not_mutate_the_original_source() {
    let original = source();
    let mut polls = 0;
    let accepted = shock(&original, 225.0, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    let mut seen = 0;
    assert_eq!(
        shock(&original, 225.0, controls(), || {
            seen += 1;
            seen == polls
        })
        .unwrap_err(),
        Error::Cancelled
    );
    assert_eq!(seen, polls);
    assert_eq!(
        shock(&original, 225.0, controls(), || false).unwrap(),
        accepted
    );
    assert_eq!(
        shock(&original, 225.0, controls(), || true).unwrap_err(),
        Error::Cancelled
    );
}
