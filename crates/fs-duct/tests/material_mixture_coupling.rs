//! Cross-domain regression: the same material gas mixture reaches duct
//! impedance, pitch and viscothermal propagation without a gas-specific
//! branch in the acoustic solver. All species except USSA air are synthetic
//! ideal-gas transport fixtures, NOT an experimental gas calibration.
use fs_duct::{Duct, LossModel, Segment, Termination, input_impedance};
use fs_material::gas::{ConductivityModel, GasSpec, GasState};
use fs_qty::semantic::{Composition, CompositionBasis};

const TAU: f64 = std::f64::consts::TAU;

fn medium(helium_fraction: f64, pressure: f64) -> GasState {
    let air = GasSpec::dry_air_ussa1976();
    let helium_like = GasSpec {
        molar_mass: 0.004, gamma: 5.0 / 3.0, sutherland_beta: 1.0e-6,
        sutherland_s: 80.0, conductivity: ConductivityModel::Eucken,
    };
    let x = Composition::new(CompositionBasis::MoleFraction,
        vec![1.0 - helium_fraction, helium_fraction]).unwrap();
    GasState::try_new_mixture(&[air, helium_like], &x, 300.0, pressure, None).unwrap()
}

fn duct() -> Duct {
    Duct { segments: vec![Segment::Cylinder { radius: 0.008, length: 0.45 }] }
}
fn close(a: f64, b: f64, rel: f64) {
    assert!((a - b).abs() <= rel * a.abs().max(b.abs()), "{a} != {b}");
}

/// The first |Z_in| peak of the closed tube, SOLVED from the duct's own wide-tube
/// input impedance by golden-section search on (0.8, 1.2) c/2L.
fn solved_resonance(gas: &GasState) -> f64 {
    let model = duct();
    let magnitude = |f: f64| input_impedance(&model, gas, TAU * f, LossModel::WideTube,
        Termination::Closed).unwrap().impedance.abs();
    let half_wave = gas.sound_speed / (2.0 * 0.45);
    let (mut a, mut b) = (0.8 * half_wave, 1.2 * half_wave);
    let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
    while b - a > 1e-11 * half_wave {
        let (c, d) = (b - ratio * (b - a), a + ratio * (b - a));
        if magnitude(c) > magnitude(d) { b = d } else { a = c }
    }
    0.5 * (a + b)
}

/// Independent prediction from the record's state alone (Kirchhoff's wide-tube
/// dispersion): the phase constant is (omega/c)(1 + s), with
/// s = (1 + (gamma-1)/sqrt(Pr)) / (sqrt(2) r sqrt(rho omega/mu)), and the peak
/// sits where it times L is pi. Returns (f, s) at that frequency.
fn predicted_resonance(gas: &GasState) -> (f64, f64) {
    let dispersion = |f: f64| (1.0 + (gas.gamma - 1.0) / gas.prandtl.sqrt())
        / (2.0_f64.sqrt() * 0.008 * (gas.density * TAU * f / gas.dynamic_viscosity).sqrt());
    let mut f = gas.sound_speed / (2.0 * 0.45);
    for _ in 0..64 { f = gas.sound_speed / (2.0 * 0.45 * (1.0 + dispersion(f))); }
    (f, dispersion(f))
}

#[test]
fn composition_enters_physical_pitch_and_impedance_without_dedicated_solver_code() {
    let air_rich = medium(0.1, 101325.0);
    let helium_rich = medium(0.8, 101325.0);
    let model = duct();
    assert!(helium_rich.sound_speed > air_rich.sound_speed);
    assert!(helium_rich.density < air_rich.density);

    // The half-wave |Z_in| peak of the SAME geometry, solved, tracks the changed
    // gas. The prediction is first order in s. Damping that grows as sqrt(f)
    // moves the |Z| peak below the phase condition by s^2/2 (from |cot kL| with
    // alpha L = pi s), so each gas's solved peak must lie within s^2 of it.
    let mut solved = Vec::new();
    for gas in [&air_rich, &helium_rich] {
        let (f_solved, (f_predicted, s)) = (solved_resonance(gas), predicted_resonance(gas));
        let error = (f_solved / f_predicted - 1.0).abs();
        println!("{{\"suite\":\"fs-duct\",\"case\":\"mixture-solved-resonance\",\"c\":{},\"f_solved\":{f_solved:.6},\"f_predicted\":{f_predicted:.6},\"s\":{s:.4e},\"error\":{error:.3e},\"bound_s2\":{:.3e}}}",
            gas.sound_speed, s * s);
        assert!(error <= s * s, "solved {f_solved} Hz against predicted {f_predicted} Hz: {error:e} > s^2 {:e}", s * s);
        solved.push((f_solved, f_predicted, s));
    }
    let [(f_air, p_air, s_air), (f_helium, p_helium, s_helium)] = [solved[0], solved[1]];
    // Direction and amount: the helium-rich peak is higher, by the predicted ratio.
    assert!(f_helium > f_air);
    assert!(((f_helium / f_air) / (p_helium / p_air) - 1.0).abs() <= s_air * s_air + s_helium * s_helium);
    // Evaluate the same dimensionless acoustic frequency kL=pi/3:
    // lossless impedance is Z0 times a geometry-only function, so the
    // ratio must equal the material-derived characteristic impedances.
    let frequency = |gas: &GasState| gas.sound_speed / (6.0 * 0.45);
    let za = input_impedance(&model, &air_rich, TAU * frequency(&air_rich),
        LossModel::Lossless, Termination::Closed).unwrap().impedance;
    let zb = input_impedance(&model, &helium_rich, TAU * frequency(&helium_rich),
        LossModel::Lossless, Termination::Closed).unwrap().impedance;
    close(zb.im / za.im,
        helium_rich.characteristic_impedance / air_rich.characteristic_impedance, 1e-10);
    assert!(za.re.abs() < 1e-9 * za.im.abs());
    assert!(zb.re.abs() < 1e-9 * zb.im.abs());
}

#[test]
fn pressure_changes_acoustic_loading_without_false_ideal_gas_detuning() {
    let low = medium(0.3, 50000.0);
    let high = medium(0.3, 100000.0);
    close(low.sound_speed, high.sound_speed, 2e-12);
    close(high.density, 2.0 * low.density, 2e-12);
    let model = duct();
    let a = input_impedance(&model, &low, TAU * 300.0,
        LossModel::Lossless, Termination::Closed).unwrap().impedance;
    let b = input_impedance(&model, &high, TAU * 300.0,
        LossModel::Lossless, Termination::Closed).unwrap().impedance;
    close(b.im, 2.0 * a.im, 1e-10);
    let loss_a = input_impedance(&model, &low, TAU * 300.0,
        LossModel::WideTube, Termination::Closed).unwrap().impedance;
    let loss_b = input_impedance(&model, &high, TAU * 300.0,
        LossModel::WideTube, Termination::Closed).unwrap().impedance;
    // Viscothermal loss scales with kinematic viscosity and so does
    // NOT simply double with acoustic impedance under pressure changes.
    assert!((loss_b.re - 2.0 * loss_a.re).abs() > 1e-5 * loss_b.re.abs().max(1.0));
}
