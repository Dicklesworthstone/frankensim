//! G0/G1/G3/G4 conserved U,V closure through the original phase and caloric laws.
use super::{close, common};
use common::{MIXTURE_MASS, model, source};
use fs_material::thermochemical_gas::nonideal::{
    ThermochemicalPengRobinsonModel,
    flash::{
        ThermochemicalFlashError as E,
        equilibrium::{
            PhaseEnthalpyControls, ResolvedThermochemicalPhaseState,
            volume::{PhaseEnergyVolumeControls, PhaseEnergyVolumeSolution},
        },
    },
};
use fs_qty::{Pressure, Temperature};

fn controls() -> PhaseEnergyVolumeControls {
    PhaseEnergyVolumeControls::new(
        [Pressure::new(2e6), Pressure::new(4e6)],
        PhaseEnthalpyControls::new([Temperature::new(170.0), Temperature::new(330.0)]),
    )
}
fn pt(
    source: &ThermochemicalPengRobinsonModel,
    t: f64,
    p: f64,
) -> ResolvedThermochemicalPhaseState {
    source
        .evaluate_phase_equilibrium(
            Temperature::new(t),
            Pressure::new(p),
            controls().phase.phase_options(),
            || false,
        )
        .unwrap()
}
fn volume(state: &ResolvedThermochemicalPhaseState) -> f64 {
    state.molar_volume_m3_mol() / MIXTURE_MASS
}
fn heat(
    source: &ThermochemicalPengRobinsonModel,
    initial: &ResolvedThermochemicalPhaseState,
    target_volume: f64,
    delta: f64,
) -> PhaseEnergyVolumeSolution {
    source
        .phase_equilibrium_internal_energy_volume_increment(
            Temperature::new(initial.temperature_k()),
            Pressure::new(initial.pressure_pa()),
            target_volume,
            delta,
            controls(),
            || false,
        )
        .unwrap()
}

#[test]
fn g12_uv_recovers_independent_homogeneous_and_coexistence_oracles() {
    let source = source();
    // Independent 50-digit cubic/fugacity solutions retained in g7. The
    // coexistence U follows the independent caloric H and original U=H-pV.
    for (t, u, rho, split) in [
        (
            180.0,
            -209_511.548_518_424_48,
            496.760_385_342_204_13,
            false,
        ),
        (
            220.0,
            -16_172.250_165_027_297 - 3e6 / 130.846_120_827_461_9,
            130.846_120_827_461_9,
            true,
        ),
        (280.0, 264_458.439_200_782_34, 36.471_975_375_322_785, false),
    ] {
        let target_volume = 1.0 / rho;
        let solved = source
            .phase_equilibrium_internal_energy_volume(target_volume, u, controls(), || false)
            .unwrap();
        close(solved.state().temperature_k(), t, 3e-7);
        close(solved.state().pressure_pa(), 3e6, 0.02);
        close(volume(solved.state()), target_volume, 1e-12);
        close(solved.state().specific_internal_energy_j_kg(), u, 1e-5);
        close(
            solved.energy_residual_j_kg(),
            solved.state().specific_internal_energy_j_kg() - u,
            0.0,
        );
        close(
            solved.volume_residual_m3_kg(),
            volume(solved.state()) - target_volume,
            0.0,
        );
        assert_eq!(solved.state().phase().two_phase().is_some(), split);
        assert_eq!(solved.accepted_increment_j_kg(), None);
        assert_eq!(solved.entropy_change_j_kg_k(), None);
        assert!(solved.pressure_evaluations() >= 3);
        assert!(solved.pt_evaluations() > solved.pressure_evaluations());
        assert!(solved.pt_evaluations() <= controls().maximum_pt_evaluations);
        assert!(
            solved.eos_evaluation_upper_bound()
                <= controls().maximum_eos_evaluation_bound().unwrap()
        );
        assert_eq!(solved.state().constitutive_law(), source.mechanics());
    }
}

#[test]
fn g13_uv_increment_crosses_phase_boundaries_without_species_datum_dependence() {
    let plain = source();
    let shifted = model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]);
    let start = pt(&plain, 180.0, 3e6);
    let shifted_start = pt(&shifted, 180.0, 3e6);
    let target = pt(&plain, 280.0, 3e6);
    let delta = 473_969.987_719_206_8;
    let a = heat(&plain, &start, volume(&target), delta);
    let b = heat(&shifted, &shifted_start, volume(&target), delta);
    close(a.state().temperature_k(), 280.0, 3e-7);
    close(a.state().pressure_pa(), 3e6, 0.02);
    close(a.accepted_increment_j_kg().unwrap(), delta, 1e-5);
    close(
        a.entropy_change_j_kg_k().unwrap(),
        2421.997_354_083_236,
        2e-6,
    );
    assert_eq!(
        a.state().temperature_k().to_bits(),
        b.state().temperature_k().to_bits()
    );
    assert_eq!(
        a.state().pressure_pa().to_bits(),
        b.state().pressure_pa().to_bits()
    );
    assert_eq!(
        a.state().density_kg_m3().to_bits(),
        b.state().density_kg_m3().to_bits()
    );
    assert_eq!(a.accepted_increment_j_kg(), b.accepted_increment_j_kg());
    assert_eq!(a.entropy_change_j_kg_k(), b.entropy_change_j_kg_k());
    assert_eq!(a.pt_evaluations(), b.pt_evaluations());
    let back = heat(&plain, a.state(), volume(&start), -delta);
    close(back.state().temperature_k(), 180.0, 5e-7);
    close(back.state().pressure_pa(), 3e6, 0.04);
    assert!(matches!(
        shifted.phase_equilibrium_internal_energy_volume(
            volume(&target),
            b.state().specific_internal_energy_j_kg(),
            controls(),
            || false
        ),
        Err(E::Invalid(_))
    ));
}

#[test]
fn g14_closed_volume_heat_changes_pressure_phase_amount_and_entropy() {
    let source = source();
    let start = pt(&source, 220.0, 3e6);
    let v = volume(&start);
    let warmed = heat(&source, &start, v, 1000.0);
    assert!(warmed.state().temperature_k() > start.temperature_k());
    assert!(warmed.state().pressure_pa() > start.pressure_pa());
    assert!(warmed.entropy_change_j_kg_k().unwrap() > 0.0);
    assert_ne!(
        warmed.state().vapor_mass_fraction(),
        start.vapor_mass_fraction()
    );
    close(volume(warmed.state()), v, 1e-12);
    close(warmed.accepted_increment_j_kg().unwrap(), 1000.0, 1e-5);
    close(
        warmed.state().specific_internal_energy_j_kg() - start.specific_internal_energy_j_kg(),
        1000.0,
        1e-5,
    );
    // No separate latent source is added: the original phase internal energy
    // already includes the accepted change in phase amount.
    close(
        warmed.state().specific_enthalpy_j_kg() - start.specific_enthalpy_j_kg(),
        1000.0 + (warmed.state().pressure_pa() - start.pressure_pa()) * v,
        2e-5,
    );
    let cooled = heat(&source, warmed.state(), v, -1000.0);
    close(cooled.state().temperature_k(), start.temperature_k(), 4e-7);
    close(cooled.state().pressure_pa(), start.pressure_pa(), 0.03);
    close(
        cooled.entropy_change_j_kg_k().unwrap() + warmed.entropy_change_j_kg_k().unwrap(),
        0.0,
        2e-6,
    );
    // A small nonzero increment must not disappear under a loose user tolerance.
    let mut loose = controls();
    loose.internal_energy_tolerance_j_kg = 1.0;
    let small = source
        .phase_equilibrium_internal_energy_volume_increment(
            Temperature::new(220.0),
            Pressure::new(3e6),
            v,
            1e-3,
            loose,
            || false,
        )
        .unwrap();
    assert!(small.accepted_increment_j_kg().unwrap() > 0.0);
    close(small.accepted_increment_j_kg().unwrap(), 1e-3, 2.5e-4);
}

#[test]
fn g15_uv_zero_heat_is_exact_at_one_original_pt_call() {
    let source = model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]);
    let mut cap = controls();
    cap.maximum_pt_evaluations = 1;
    cap.maximum_pressure_evaluations = 1;
    cap.phase.enthalpy.maximum_flash_evaluations = 1;
    for t in [180.0, 220.0, 280.0] {
        let start = pt(&source, t, 3e6);
        let solved = source
            .phase_equilibrium_internal_energy_volume_increment(
                Temperature::new(t),
                Pressure::new(3e6),
                volume(&start),
                0.0,
                cap,
                || false,
            )
            .unwrap();
        assert_eq!(solved.state(), &start);
        assert_eq!(solved.pt_evaluations(), 1);
        assert_eq!(solved.pressure_evaluations(), 0);
        assert_eq!(solved.accepted_increment_j_kg(), Some(0.0));
        assert_eq!(solved.entropy_change_j_kg_k(), Some(0.0));
        assert_eq!(solved.energy_residual_j_kg(), 0.0);
        assert_eq!(solved.volume_residual_m3_kg(), 0.0);
        assert_eq!(
            solved.eos_evaluation_upper_bound(),
            start.phase().eos_evaluation_upper_bound()
        );
    }
}

#[test]
fn g16_uv_bounds_total_nested_work_and_preserves_phase_refusals() {
    let source = source();
    let start = pt(&source, 220.0, 3e6);
    let run = |cap| {
        source.phase_equilibrium_internal_energy_volume_increment(
            Temperature::new(220.0),
            Pressure::new(3e6),
            volume(&start),
            1000.0,
            cap,
            || false,
        )
    };
    let mut cap = controls();
    cap.maximum_pt_evaluations = 3;
    assert!(matches!(
        run(cap),
        Err(E::WorkLimit {
            flash_evaluations: 3
        })
    ));
    cap = controls();
    cap.maximum_pressure_evaluations = 1;
    assert!(matches!(
        run(cap),
        Err(E::PressureWorkLimit {
            pressure_evaluations: 1
        })
    ));
    cap = controls();
    cap.phase.enthalpy.maximum_flash_evaluations = 2;
    assert!(matches!(
        run(cap),
        Err(E::WorkLimit {
            flash_evaluations: 3
        })
    ));
    cap = controls();
    cap.phase.stability.maximum_iterations = 1;
    assert!(matches!(run(cap), Err(E::Phase(_))));
    cap = controls();
    cap.pressure_bracket = [Pressure::new(1.7e6), Pressure::new(2e6)];
    assert!(matches!(
        source.phase_equilibrium_internal_energy_volume(
            volume(&start),
            start.specific_internal_energy_j_kg(),
            cap,
            || false
        ),
        Err(E::OutsideVolumeBracket)
    ));
    let limited = model([0.0; 2], [0.0; 2], [200.0, 260.0]);
    assert!(matches!(
        limited.phase_equilibrium_internal_energy_volume(
            volume(&start),
            start.specific_internal_energy_j_kg(),
            controls(),
            || false
        ),
        Err(E::Material(_))
    ));
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            source.phase_equilibrium_internal_energy_volume(invalid, 0.0, controls(), || false),
            Err(E::Invalid(_))
        ));
    }
    cap = controls();
    cap.maximum_pt_evaluations = 16_385;
    assert!(matches!(run(cap), Err(E::Invalid(_))));
}

#[test]
fn g17_uv_final_publication_cancellation_is_atomic_and_replayable() {
    let source = source();
    let start = pt(&source, 220.0, 3e6);
    let mut count = 0;
    let solved = source
        .phase_equilibrium_internal_energy_volume(
            volume(&start),
            start.specific_internal_energy_j_kg(),
            controls(),
            || {
                count += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    assert_eq!(
        source
            .phase_equilibrium_internal_energy_volume(
                volume(&start),
                start.specific_internal_energy_j_kg(),
                controls(),
                || {
                    seen += 1;
                    seen == count
                }
            )
            .unwrap_err(),
        E::Cancelled
    );
    assert_eq!(seen, count);
    let replay = source
        .phase_equilibrium_internal_energy_volume(
            volume(&start),
            start.specific_internal_energy_j_kg(),
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(solved, replay);
}
