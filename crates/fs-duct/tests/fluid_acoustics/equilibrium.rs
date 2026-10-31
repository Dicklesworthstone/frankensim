//! G1/G3: the explicit fully relaxed phase response drives actual pipe ports.
//! These are synthetic calorics and PR critical data, not measured acoustics.

use super::*;

#[path = "../../../fs-material/tests/common/nonideal_equilibrium.rs"]
mod source;

#[test]
fn phase_equilibrium_sound_speed_reaches_wave_and_closed_pipe_impedance() {
    let radius: f64 = 0.001;
    let length: f64 = 0.2;
    let area = core::f64::consts::PI * radius.powi(2);
    let duct = cylinder(radius, length);
    // Independent Maxwell-area coexistence and entropy-constrained compression
    // oracles at 150 K, with the shared fixture's original Cp/R=3.5. Rounded
    // values retain considerably more precision than this propagation check.
    let cases = [
        (0.01, 317.158_519_69, 29.335_537_60),
        (0.5, 31.693_091_062_654_844, 178.869_227_66),
        (0.9, 18.269_517_85, 252.240_558_19),
    ];
    let mut original_waves = Vec::new();
    for datum in [0.0, 1e20] {
        let model = source::model(datum, [120.0, 185.0]);
        for (quality, density, sound_speed) in cases {
            let state = model
                .evaluate_pure_equilibrium_density_temperature(
                    density,
                    Temperature::new(150.0),
                    source::controls(),
                    || false,
                )
                .unwrap();
            assert!((state.vapor_mass_fraction() - quality).abs() < 1e-9);
            let acoustic = state.equilibrium_acoustics().unwrap();
            assert!(std::ptr::eq(acoustic.source(), &state));
            for omega in [20.0, 200.0] {
                let wave = fluid::segment_wave(&acoustic, radius, omega).unwrap();
                let k = omega / sound_speed;
                let zc = density * sound_speed / area;
                near(wave.wavenumber, C64::from_re(k), k, 2e-8);
                near(wave.characteristic_impedance, C64::from_re(zc), zc, 2e-8);
                let closed = solve(&duct, &acoustic, omega, None, WallModel::Rigid);
                near(
                    closed.z_in,
                    C64::new(0.0, zc / (k * length).tan()),
                    zc,
                    5e-8,
                );
                let values = (wave.wavenumber, wave.characteristic_impedance, closed.z_in);
                if datum == 0.0 {
                    original_waves.push(values);
                } else {
                    let index = cases.iter().position(|case| case.0 == quality).unwrap() * 2
                        + usize::from(omega == 200.0);
                    assert_eq!(values, original_waves[index]);
                }
            }
        }
    }
}
