//! Finite-time accuracy against an independently integrated physical rate law.

use super::*;

#[test]
fn atomic_transient_refinement_converges_to_independent_beb_trajectory() {
    // Independent oracle: NIST generalized BEB Eq.1 with H 1s B=U=13.6057 eV,
    // N=1, Q=0.5668 and a0=0.52918e-10 m, integrated against the Maxwellian
    // energy density by SciPy QUADPACK. No production Rust rate or fitted rate
    // surrogate enters the oracle. Sources:
    // https://physics.nist.gov/PhysRefData/Ionization/Eqs/latex.html
    // https://physics.nist.gov/cgi-bin/Ionization/ion_data.php?differential=Y&id=HI&initial=&ision=I
    //
    // The conserved-energy trajectory uses
    // Te(x)=[(3R/2)*0.05*30000-chi*(x-0.05)]/[(3R/2)*x], and
    // dx/dt=NA*k(Te)*x*[(1-x)-x*x/Kc(Te)] for one mole nuclei in 1 m³.
    // An adaptive DOP853 solve (rtol=2e-13) and a separate scalar solve of
    // t=integral(dx/(dx/dt)) agree at t=1e-7 s to 8e-17 in x. The latter gives
    // x=0.056015150225581782 and Te=15481.401668882725 K. This checks the
    // numerical transient of the declared ground-state model, not experiment.
    let expected_fraction = 0.056_015_150_225_581_78;
    let expected_temperature = 15_481.401_668_882_725;
    let (state, channel) = initial(0.05, 30_000.0, 4000.0, 0.0, 0.0);
    let mut fraction_errors = [0.0; 3];
    let mut temperature_errors = [0.0; 3];
    let mut accepted_steps = [0; 3];
    for (index, relative) in [1e-3, 1e-5, 1e-7].into_iter().enumerate() {
        let mut refined = controls();
        refined.relative_tolerance = relative;
        let result = advance_atomic_ionization(
            &state,
            &channel,
            Power::new(0.0),
            Power::new(0.0),
            Time::new(1e-7),
            refined,
            || false,
        )
        .unwrap();
        fraction_errors[index] =
            (result.state().amounts()[2].value() / expected_fraction - 1.0).abs();
        temperature_errors[index] =
            (result.state().electron_temperature().value() / expected_temperature - 1.0).abs();
        accepted_steps[index] = result.balance().accepted_steps;
        assert_eq!(result.elapsed_time(), Time::new(1e-7));
        assert_eq!(
            result.state().heavy_temperature(),
            state.heavy_temperature()
        );
        assert!(result.balance().maximum_normalized_local_error <= 1.0);
        assert!(result.balance().energy_residual_j.abs() <= refined.energy_tolerance_j);
    }
    eprintln!(
        "atomic transient relative fraction errors={fraction_errors:?}, relative Te errors={temperature_errors:?}, accepted steps={accepted_steps:?}"
    );
    // Backward Euler with full/two-half local error control has first-order
    // global accuracy. A 100-fold local tolerance change should substantially
    // reduce the global error; the local tolerance is not a global guarantee.
    for errors in [fraction_errors, temperature_errors] {
        assert!(errors[1] < 0.2 * errors[0], "transient errors={errors:?}");
        assert!(errors[2] < 0.2 * errors[1], "transient errors={errors:?}");
    }
    assert!(
        fraction_errors[2] < 1e-5,
        "fraction errors={fraction_errors:?}"
    );
    assert!(
        temperature_errors[2] < 2e-4,
        "Te errors={temperature_errors:?}"
    );
}
