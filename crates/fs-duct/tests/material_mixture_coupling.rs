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

#[test]
fn composition_enters_physical_pitch_and_impedance_without_dedicated_solver_code() {
    let air_rich = medium(0.1, 101325.0);
    let helium_rich = medium(0.8, 101325.0);
    let model = duct();
    assert!(helium_rich.sound_speed > air_rich.sound_speed);
    assert!(helium_rich.density < air_rich.density);

    // Quarter-wave resonance of the SAME geometry tracks the changed gas.
    let f_air = air_rich.sound_speed / (4.0 * 0.45);
    let f_helium = helium_rich.sound_speed / (4.0 * 0.45);
    assert!(f_helium > f_air);
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
