//! Equilibrium acoustics with independent phase/caloric and constrained-state oracles.
#[path = "common/nonideal_flash.rs"]
mod common;

use common::{MIXTURE_MASS, model, source};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics, FrozenHomogeneousMixture};
use fs_material::thermochemical_gas::nonideal::{
    ThermochemicalPengRobinsonModel,
    flash::{
        ThermochemicalFlashError as E,
        equilibrium::{
            PhaseEnthalpyControls,
            acoustics::{PhaseEquilibriumAcousticControls, PhaseEquilibriumAcoustics},
            entropy::PhaseEntropyControls,
            volume::PhaseEnergyVolumeControls,
        },
    },
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseRootV1;

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "{actual} != {expected}, difference {}",
        actual - expected,
    );
}
fn response(source: &ThermochemicalPengRobinsonModel, t: f64) -> PhaseEquilibriumAcoustics {
    source
        .phase_equilibrium_acoustics(
            Temperature::new(t),
            Pressure::new(3e6),
            PhaseEquilibriumAcousticControls::default(),
            || false,
        )
        .unwrap()
}

#[test]
fn binary_coexistence_matches_independent_cubic_fugacity_entropy_derivatives() {
    let source = source();
    let actual = response(&source, 220.0);
    let d = actual.derivatives();
    // Independent 65-digit simultaneous Newton solve of the two PR cubics
    // and both fugacity equalities, followed by entropy/volume differentiation
    // with relative step 1e-12. The ideal calorics are the fixture's original
    // 3.5R/4.5R sources; phase compositions and amounts both re-equilibrate.
    close(actual.state().density_kg_m3(), 130.846_120_827_461_9, 2e-9);
    close(
        actual.state().sound_speed_m_s(),
        129.899_258_105_219_43,
        2e-7,
    );
    for (found, expected) in [
        (actual.isobaric_heat_capacity_j_kg_k(), 7285.853_312_349_442),
        (d.isochoric_heat_capacity_j_kg_k, 3639.875_747_764_385_4),
        (d.pressure_temperature_at_density_pa_k, 48906.461_943_654_81),
        (
            d.pressure_density_at_temperature_m2_s2,
            8429.842_815_976_168,
        ),
        (d.pressure_density_at_entropy_m2_s2, 16873.817_256_286_417),
        (d.pressure_energy_at_density_kg_m3, 13.436_299_844_491_45),
        (d.pressure_density_at_energy_m2_s2, 14519.422_736_051_248),
    ] {
        close(found, expected, 2e-7);
    }
    assert!(actual.source().phase().two_phase().is_some());
    assert_eq!(actual.source().constitutive_law(), source.mechanics());
    assert_eq!(
        actual.state().density_kg_m3(),
        actual.source().density_kg_m3()
    );
    assert_eq!(actual.pt_evaluations(), 9);
    let controls = PhaseEquilibriumAcousticControls::default();
    assert_eq!(actual.source().phase().options(), controls.phase);
    assert!(
        actual.eos_evaluation_upper_bound() <= controls.maximum_eos_evaluation_bound().unwrap()
    );
    assert!(actual.relative_refinement_error() <= controls.relative_tolerance);
    assert!(actual.relative_maxwell_residual() <= controls.relative_tolerance);
    let consumer: &dyn FluidAcoustics = &actual;
    assert_eq!(consumer.acoustic_state().unwrap(), actual.state());

    // A separately assembled frozen-phase Wood mixture has the same density
    // but excludes the equilibrium mass/heat transfer responsible for softening.
    let split = actual.source().phase().two_phase().unwrap();
    let frozen_phase = |root: &PengRobinsonPhaseRootV1, fractions: &[f64]| {
        let mut cp = 0.0;
        let mut mass = 0.0;
        for (species, &fraction) in source
            .source()
            .thermodynamics()
            .components()
            .iter()
            .zip(fractions)
        {
            cp += fraction
                * species
                    .evaluate_properties(Temperature::new(220.0))
                    .unwrap()
                    .cp()
                    .value();
            mass += fraction * species.molar_mass().value();
        }
        let mechanics = root.caloric_mechanics(cp - R, mass).unwrap();
        FluidAcousticState::try_new(root.mass_density_kg_per_m3(), mechanics.sound_speed_m_s)
            .unwrap()
    };
    let wood = FrozenHomogeneousMixture::from_mass_fraction(
        frozen_phase(split.liquid(), split.liquid_mole_fractions()),
        frozen_phase(split.vapor(), split.vapor_mole_fractions()),
        actual.source().vapor_mass_fraction().unwrap(),
    )
    .unwrap();
    close(
        wood.state().density_kg_m3(),
        actual.state().density_kg_m3(),
        2e-10,
    );
    close(wood.state().sound_speed_m_s(), 203.416_997_714_477_97, 2e-9);
    assert!(actual.state().sound_speed_m_s() < 0.7 * wood.state().sound_speed_m_s());
}

#[test]
fn homogeneous_liquid_and_vapor_recover_original_analytic_caloric_mechanics() {
    let source = source();
    for temperature in [180.0, 280.0] {
        let actual = response(&source, temperature);
        let state = actual.source();
        let root = state.phase().homogeneous().unwrap();
        let exact = root
            .caloric_mechanics(state.reference().properties().cv().value(), MIXTURE_MASS)
            .unwrap();
        let d = actual.derivatives();
        let vm = root.molar_volume_m3_per_mol();
        close(
            actual.state().sound_speed_m_s(),
            exact.sound_speed_m_s,
            2e-7,
        );
        close(
            d.isochoric_heat_capacity_j_kg_k,
            exact.molar_cv_j_mol_k / MIXTURE_MASS,
            2e-7,
        );
        close(
            d.pressure_temperature_at_density_pa_k,
            exact.dp_dt_v_pa_k,
            2e-7,
        );
        close(
            d.pressure_density_at_temperature_m2_s2,
            -exact.dp_dv_t_pa_mol_m3 * vm * vm / MIXTURE_MASS,
            2e-7,
        );
        close(
            actual.isobaric_heat_capacity_j_kg_k(),
            (exact.molar_cv_j_mol_k
                - temperature * exact.dp_dt_v_pa_k * exact.dp_dt_v_pa_k / exact.dp_dv_t_pa_mol_m3)
                / MIXTURE_MASS,
            2e-7,
        );
    }
}

#[test]
fn pressure_jacobian_matches_independent_ps_and_conserved_uv_inversions() {
    let source = source();
    let actual = response(&source, 220.0);
    let d = actual.derivatives();
    let rho = actual.source().density_kg_m3();
    let bracket = [Temperature::new(210.0), Temperature::new(230.0)];
    let entropy_controls = PhaseEntropyControls::new(bracket);
    let isentrope = |p| {
        source
            .phase_equilibrium_pressure_entropy_increment(
                Temperature::new(220.0),
                Pressure::new(3e6),
                Pressure::new(p),
                0.0,
                entropy_controls,
                || false,
            )
            .unwrap()
    };
    let dp = 300.0;
    let low = isentrope(3e6 - dp);
    let high = isentrope(3e6 + dp);
    let dr = high.state().density_kg_m3() - low.state().density_kg_m3();
    close(d.pressure_density_at_entropy_m2_s2, 2.0 * dp / dr, 3e-6);
    close(
        d.temperature_density_at_entropy_k_m3_kg,
        (high.state().temperature_k() - low.state().temperature_k()) / dr,
        3e-6,
    );

    let controls = PhaseEnergyVolumeControls::new(
        [Pressure::new(2.9e6), Pressure::new(3.1e6)],
        PhaseEnthalpyControls::new(bracket),
    );
    let uv = |density, heat| {
        source
            .phase_equilibrium_internal_energy_volume_increment(
                Temperature::new(220.0),
                Pressure::new(3e6),
                1.0 / density,
                heat,
                controls,
                || false,
            )
            .unwrap()
    };
    let du = 10.0;
    let low = uv(rho, -du);
    let high = uv(rho, du);
    close(
        d.pressure_energy_at_density_kg_m3,
        (high.state().pressure_pa() - low.state().pressure_pa()) / (2.0 * du),
        2e-5,
    );
    close(
        d.isochoric_heat_capacity_j_kg_k,
        2.0 * du / (high.state().temperature_k() - low.state().temperature_k()),
        2e-5,
    );
    let dr = rho * 1e-4;
    let low = uv(rho - dr, 0.0);
    let high = uv(rho + dr, 0.0);
    close(
        d.pressure_density_at_energy_m2_s2,
        (high.state().pressure_pa() - low.state().pressure_pa()) / (2.0 * dr),
        2e-5,
    );
    // The conservative Euler closure and the entropy-coordinate closure agree.
    close(
        d.pressure_density_at_entropy_m2_s2,
        d.pressure_density_at_energy_m2_s2
            + actual.source().pressure_pa() / rho * d.pressure_energy_at_density_kg_m3 / rho,
        5e-15,
    );
}

#[test]
fn source_datums_cannot_change_equilibrium_response_and_step_refinement_converges() {
    let ordinary = source();
    let shifted = model([1e20, -2e20], [-3e20, 4e20], [100.0, 450.0]);
    let original = response(&ordinary, 220.0);
    let other = response(&shifted, 220.0);
    assert_ne!(
        original.source().specific_entropy_j_kg_k(),
        other.source().specific_entropy_j_kg_k()
    );
    assert_ne!(
        original.source().specific_internal_energy_j_kg(),
        other.source().specific_internal_energy_j_kg()
    );
    assert_eq!(original.state(), other.state());
    assert_eq!(original.derivatives(), other.derivatives());
    assert_eq!(
        original.isobaric_heat_capacity_j_kg_k(),
        other.isobaric_heat_capacity_j_kg_k()
    );
    let fine = ordinary
        .phase_equilibrium_acoustics(
            Temperature::new(220.0),
            Pressure::new(3e6),
            PhaseEquilibriumAcousticControls {
                relative_temperature_step: 5e-4,
                relative_pressure_step: 5e-4,
                ..Default::default()
            },
            || false,
        )
        .unwrap();
    close(
        fine.state().sound_speed_m_s(),
        original.state().sound_speed_m_s(),
        2e-7,
    );
    assert!(fine.relative_refinement_error() < 0.5 * original.relative_refinement_error());
}

#[test]
fn phase_boundaries_source_limits_unresolved_accuracy_and_cancellation_refuse() {
    let source = source();
    let controls = PhaseEquilibriumAcousticControls::default();
    // Independent bubble point: 205.7464092585643 K at 3 MPa.
    // The center is separated but the coarse low-T probe is liquid.
    assert!(matches!(
        source.phase_equilibrium_acoustics(
            Temperature::new(205.8),
            Pressure::new(3e6),
            controls,
            || false
        ),
        Err(E::Invalid(
            "equilibrium acoustic stencil crosses a phase boundary"
        ))
    ));
    let unresolved = PhaseEquilibriumAcousticControls {
        relative_tolerance: 1e-12,
        ..controls
    };
    assert!(matches!(
        source.phase_equilibrium_acoustics(
            Temperature::new(220.0),
            Pressure::new(3e6),
            unresolved,
            || false
        ),
        Err(E::Invalid(
            "equilibrium acoustic derivatives did not converge"
        ))
    ));
    let tiny = PhaseEquilibriumAcousticControls {
        relative_temperature_step: f64::MIN_POSITIVE,
        ..controls
    };
    assert!(matches!(
        source.phase_equilibrium_acoustics(
            Temperature::new(220.0),
            Pressure::new(3e6),
            tiny,
            || false
        ),
        Err(E::Invalid(
            "equilibrium acoustic stencil is not representable"
        ))
    ));
    let limited = model([0.0; 2], [0.0; 2], [219.95, 220.05]);
    assert!(matches!(
        limited.phase_equilibrium_acoustics(
            Temperature::new(220.0),
            Pressure::new(3e6),
            controls,
            || false
        ),
        Err(E::Material(_))
    ));
    let mut capped = controls;
    capped.phase.stability.maximum_iterations = 1;
    assert!(matches!(
        source.phase_equilibrium_acoustics(
            Temperature::new(220.0),
            Pressure::new(3e6),
            capped,
            || false
        ),
        Err(E::Phase(_))
    ));
    assert!(matches!(
        source.phase_equilibrium_acoustics(
            Temperature::new(220.0),
            Pressure::new(3e6),
            controls,
            || true
        ),
        Err(E::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        source.phase_equilibrium_acoustics(
            Temperature::new(220.0),
            Pressure::new(3e6),
            controls,
            || {
                polls += 1;
                polls >= 25
            }
        ),
        Err(E::Cancelled)
    ));
}
