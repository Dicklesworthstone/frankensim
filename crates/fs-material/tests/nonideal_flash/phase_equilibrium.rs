//! Independent physical regressions for phase-boundary-crossing caloric closure.
use super::{close, common};
use common::{METHANE, MIXTURE_MASS, model, source};
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_material::thermochemical_gas::nonideal::flash::{
    ThermochemicalFlashError as E,
    equilibrium::{PhaseEnthalpyControls, PhaseEnthalpySolution, ResolvedThermochemicalPhaseState},
};
use fs_qty::{MolarMass, Pressure, Temperature};
use fs_thermochem::eos::peng_robinson::{
    PengRobinsonComponentV1, PengRobinsonMixtureV1, flash::PengRobinsonFlashErrorV1 as PhaseError,
};

fn controls() -> PhaseEnthalpyControls {
    PhaseEnthalpyControls::new([Temperature::new(175.0), Temperature::new(285.0)])
}
fn pt(source: &ThermochemicalPengRobinsonModel, t: f64) -> ResolvedThermochemicalPhaseState {
    source
        .evaluate_phase_equilibrium(
            Temperature::new(t),
            Pressure::new(3e6),
            controls().phase_options(),
            || false,
        )
        .unwrap()
}
fn heat(
    source: &ThermochemicalPengRobinsonModel,
    initial: f64,
    delta: f64,
) -> PhaseEnthalpySolution {
    source
        .phase_equilibrium_pressure_enthalpy_increment(
            Temperature::new(initial),
            Pressure::new(3e6),
            Pressure::new(3e6),
            delta,
            controls(),
            || false,
        )
        .unwrap()
}

#[test]
fn g7_screened_pt_matches_independent_homogeneous_and_split_oracles() {
    let source = source();
    // 50-digit independent Decimal cubic/fugacity calculations, not the
    // production stationary iteration. Original component NASA Cp/R = 3.5, 4.5.
    for (t, h, s, rho, u) in [
        (
            180.0,
            -203_472.419_617_408_38,
            4589.992_780_575_417,
            496.760_385_342_204_13,
            -209_511.548_518_424_48,
        ),
        (
            280.0,
            346_713.375_192_810_4,
            7011.990_134_658_653,
            36.471_975_375_322_785,
            264_458.439_200_782_34,
        ),
    ] {
        let state = pt(&source, t);
        close(state.specific_enthalpy_j_kg(), h, 1e-7);
        close(state.specific_internal_energy_j_kg(), u, 1e-7);
        close(state.specific_entropy_j_kg_k(), s, 1e-8);
        close(state.density_kg_m3(), rho, 1e-9);
        close(
            state.specific_enthalpy_j_kg() - state.specific_internal_energy_j_kg(),
            state.pressure_pa() / state.density_kg_m3(),
            1e-8,
        );
        assert!(state.phase().two_phase().is_none());
        assert_eq!(state.vapor_mass_fraction(), None);
        assert_eq!(state.vapor_volume_fraction(), None);
        let screen = state.phase().stability();
        assert!(screen.instability_witness().is_none());
        assert_eq!(
            screen.eos_evaluations(),
            1 + screen.vapor_seed().iterations() + screen.liquid_seed().iterations()
        );
        assert_eq!(
            state.phase().eos_evaluation_upper_bound(),
            screen.eos_evaluations()
        );
        let root = state.phase().homogeneous().unwrap();
        assert!(root.pressure_derivatives().1 < 0.0);
        assert!(
            root.caloric_mechanics(state.reference().properties().cv().value(), MIXTURE_MASS)
                .is_some()
        );
        assert_eq!(state.constitutive_law(), source.mechanics());
        assert_eq!(state.phase().options(), controls().phase_options());
    }
    let state = pt(&source, 220.0);
    let split = state.phase().two_phase().unwrap();
    assert!(state.phase().stability().instability_witness().is_some());
    assert!(state.phase().homogeneous().is_none());
    close(split.vapor_fraction_mol(), 0.302_585_915_611_944_35, 2e-10);
    close(
        split.liquid_mole_fractions()[METHANE],
        0.364_782_477_343_904_27,
        2e-10,
    );
    close(state.density_kg_m3(), 130.846_120_827_461_9, 2e-8);
    close(
        state.specific_enthalpy_j_kg(),
        -16_172.250_165_027_297,
        2e-6,
    );
    close(
        state.vapor_mass_fraction().unwrap(),
        0.245_214_412_632_239_24,
        2e-10,
    );
}

#[test]
fn g8_enthalpy_heating_and_cooling_cross_both_independent_phase_boundaries() {
    let source = source();
    // Independent bubble/dew roots at 3 MPa are 205.7464092585643 K and
    // 250.68921833557063 K. Endpoints therefore lie on opposite homogeneous
    // sides and 220 K lies strictly inside the coexistence region.
    for (initial, increment, expected_t, split) in [
        (180.0, 187_300.169_452_381_1, 220.0, true),
        (220.0, 362_885.625_357_837_7, 280.0, false),
        (180.0, 550_185.794_810_218_8, 280.0, false),
        (280.0, -550_185.794_810_218_8, 180.0, false),
    ] {
        let solution = heat(&source, initial, increment);
        close(solution.state().temperature_k(), expected_t, 3e-7);
        close(solution.accepted_increment_j_kg().unwrap(), increment, 1e-5);
        assert!(solution.residual_j_kg().abs() <= 1e-5);
        assert_eq!(solution.state().phase().two_phase().is_some(), split);
        assert!(solution.pt_evaluations() >= 3);
        assert!(solution.eos_evaluation_upper_bound() > solution.pt_evaluations());
        assert!(
            solution.eos_evaluation_upper_bound()
                <= controls().maximum_eos_evaluation_bound().unwrap()
        );
        if initial == 180.0 && expected_t == 280.0 {
            close(
                solution.entropy_change_j_kg_k().unwrap(),
                2421.997_354_083_236,
                2e-6,
            );
        }
        if initial == 280.0 {
            close(
                solution.entropy_change_j_kg_k().unwrap(),
                -2421.997_354_083_236,
                2e-6,
            );
        }
    }
    let absolute = source
        .phase_equilibrium_pressure_enthalpy(
            Pressure::new(3e6),
            -16_172.250_165_027_297,
            controls(),
            || false,
        )
        .unwrap();
    close(absolute.state().temperature_k(), 220.0, 3e-7);
    assert!(absolute.state().phase().two_phase().is_some());
    assert_eq!(absolute.accepted_increment_j_kg(), None);
}

#[test]
fn g9_crossing_increments_ignore_species_reference_gauges_and_preserve_fixed_points() {
    let plain = source();
    let shifted = model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]);
    let a = heat(&plain, 180.0, 187_300.169_452_381_1);
    let b = heat(&shifted, 180.0, 187_300.169_452_381_1);
    assert_eq!(
        a.state().temperature_k().to_bits(),
        b.state().temperature_k().to_bits()
    );
    assert_eq!(
        a.state().density_kg_m3().to_bits(),
        b.state().density_kg_m3().to_bits()
    );
    assert_eq!(a.accepted_increment_j_kg(), b.accepted_increment_j_kg());
    assert_eq!(a.entropy_change_j_kg_k(), b.entropy_change_j_kg_k());
    let mut fixed_controls = controls();
    fixed_controls.enthalpy.maximum_flash_evaluations = 1;
    for t in [180.0, 220.0, 280.0] {
        let initial = pt(&shifted, t);
        let fixed = shifted
            .phase_equilibrium_pressure_enthalpy_increment(
                Temperature::new(t),
                Pressure::new(3e6),
                Pressure::new(3e6),
                0.0,
                fixed_controls,
                || false,
            )
            .unwrap();
        assert_eq!(fixed.state(), &initial);
        assert_eq!(fixed.pt_evaluations(), 1);
        assert_eq!(fixed.accepted_increment_j_kg(), Some(0.0));
        assert_eq!(fixed.entropy_change_j_kg_k(), Some(0.0));
        assert_eq!(
            fixed.eos_evaluation_upper_bound(),
            initial.phase().eos_evaluation_upper_bound()
        );
    }
    assert!(matches!(
        shifted.phase_equilibrium_pressure_enthalpy(
            Pressure::new(3e6),
            pt(&shifted, 220.0).specific_enthalpy_j_kg(),
            controls(),
            || false,
        ),
        Err(E::Invalid(_))
    ));
}

#[test]
fn g10_screen_never_falls_back_after_unresolved_instability_or_amount_ambiguity() {
    let source = source();
    let mut options = controls().phase_options();
    options.stability.maximum_iterations = 1;
    assert!(matches!(
        source.evaluate_phase_equilibrium(
            Temperature::new(220.0),
            Pressure::new(3e6),
            options,
            || false,
        ),
        Err(E::Phase(PhaseError::StabilityIterationStalled { .. }))
    ));
    options = controls().phase_options();
    options.flash.maximum_iterations = 1;
    // The screen at this same state is independently negative in g7. Failing
    // its split must be a refusal, never a homogeneous caloric result.
    assert!(matches!(
        source.evaluate_phase_equilibrium(
            Temperature::new(220.0),
            Pressure::new(3e6),
            options,
            || false,
        ),
        Err(E::Phase(_))
    ));
    let pure = PengRobinsonMixtureV1::new(
        vec![PengRobinsonComponentV1 {
            critical_temperature: Temperature::new(190.564),
            critical_pressure: Pressure::new(4.5992e6),
            acentric_factor: 0.01142,
            molar_mass: MolarMass::new(0.016_042_46),
        }],
        vec![1.0],
        vec![],
    )
    .unwrap();
    // Independent Maxwell equal-area saturation oracle already used by the
    // original saturation owner's tests: two stable roots, equal Gibbs energy.
    assert!(matches!(
        pure.phase_screen_pt_cancellable(
            Temperature::new(150.0),
            Pressure::new(1_047_062.532_316_561_1),
            controls().phase_options(),
            || false,
        ),
        Err(PhaseError::Unrepresentable(
            "equal-Gibbs phase amount is undetermined"
        ))
    ));
    // The old public contract stays strictly two-phase.
    assert!(matches!(
        source.two_phase_flash_pressure_enthalpy(
            Pressure::new(3e6),
            -16_172.250_165_027_297,
            controls().enthalpy,
            || false,
        ),
        Err(E::Phase(_))
    ));
}

#[test]
fn g11_crossing_respects_total_work_source_support_and_final_cancellation() {
    let source = source();
    let mut capped = controls();
    capped.enthalpy.maximum_flash_evaluations = 2;
    assert!(matches!(
        source.phase_equilibrium_pressure_enthalpy_increment(
            Temperature::new(180.0),
            Pressure::new(3e6),
            Pressure::new(3e6),
            187_300.169_452_381_1,
            capped,
            || false,
        ),
        Err(E::WorkLimit {
            flash_evaluations: 2
        })
    ));
    let limited = model([0.0; 2], [0.0; 2], [200.0, 260.0]);
    assert!(matches!(
        limited.phase_equilibrium_pressure_enthalpy(Pressure::new(3e6), 0.0, controls(), || false,),
        Err(E::Material(_))
    ));
    let mut checkpoints = 0;
    source
        .phase_equilibrium_pressure_enthalpy(
            Pressure::new(3e6),
            -16_172.250_165_027_297,
            controls(),
            || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    assert_eq!(
        source
            .phase_equilibrium_pressure_enthalpy(
                Pressure::new(3e6),
                -16_172.250_165_027_297,
                controls(),
                || {
                    seen += 1;
                    seen == checkpoints
                },
            )
            .unwrap_err(),
        E::Cancelled
    );
    assert_eq!(seen, checkpoints);
}
