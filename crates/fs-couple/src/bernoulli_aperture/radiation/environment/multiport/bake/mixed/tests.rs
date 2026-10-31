use super::*;
use crate::bernoulli_aperture::radiation::environment::multiport::source::ComputedBaffledInput;
use fs_material::gas::GasSpec;

fn gas(temperature: f64, pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, pressure).unwrap()
}

fn fixture() -> (DrivenMouthSurface, [CircularOpening; 1], BaffledMimoBakeOptions) {
    let input = ComputedBaffledInput::parse(
        "1 33 2 4 6 4 0.1 67108864 1 2 0.001 32 1048576 65").unwrap();
    (DrivenMouthSurface { center_m: [-0.01, 0.0], size_m: [0.0006, 0.0004],
        cells: [8; 2], terminal_radial_rings: 8, terminal_angular_points: 32 },
        [CircularOpening { center_m: [0.0; 2], radius_m: 0.0003 }], input.options())
}

#[test]
fn mixed_openings_retain_mutual_phase_work_sign_and_physical_flow_normalization() {
    let (mouth, terminals, options) = fixture();
    let gas = gas(293.15, 101_325.0);
    let gate = CancelGate::new_clock_free();
    let surface = MixedSurface::new(mouth, &terminals, &gas, [1.0, 4000.0], options, &gate).unwrap();
    let frequencies = [100.0, 1000.0, 4000.0];
    let rows = surface.rows(&gas, &frequencies, &gate).unwrap();
    for (&frequency, row) in frequencies.iter().zip(&rows) {
        assert_eq!(row[1], row[2]);
        // Independent compact-aperture Rayleigh limit in s=+iw convention:
        // Zmt = i*rho*omega*exp(-i*k*d)/(2*pi*d). These apertures are small
        // against both wavelength and their 10 mm physical separation.
        let distance = 0.01;
        let phase = core::f64::consts::TAU * frequency * distance / gas.sound_speed;
        let expected = C64::new(det::sin(phase), det::cos(phase))
            .scale(gas.density * frequency / distance);
        let relative = det::sqrt((row[1] - expected).norm_sq() / expected.norm_sq());
        assert!(relative < 0.002, "f={frequency}: relative mutual error {relative}");
        let same_direction = row[0].re + row[3].re + 2.0 * row[1].re;
        let opposite_direction = row[0].re + row[3].re - 2.0 * row[1].re;
        assert!(opposite_direction >= 0.0);
        assert!(opposite_direction < same_direction);
    }
    // Both unit-Q self fields share the low-frequency monopole coefficient,
    // despite unequal shape and area. A missing 1/A on either side fails this.
    let coefficient = gas.density * core::f64::consts::TAU * frequencies[0].powi(2)
        / gas.sound_speed;
    for diagonal in [0, 3] {
        assert!((rows[0][diagonal].re / coefficient - 1.0).abs() < 1e-5);
    }
}

#[test]
fn material_pressure_and_temperature_change_the_complete_mixed_image() {
    let (mouth, terminals, options) = fixture();
    let ordinary = gas(293.15, 101_325.0);
    let doubled = gas(293.15, 202_650.0);
    let warmer = gas(333.15, 101_325.0);
    let gate = CancelGate::new_clock_free();
    let surface = MixedSurface::new(mouth, &terminals, &ordinary,
        [1.0, 4000.0], options, &gate).unwrap();
    let base = surface.rows(&ordinary, &[4000.0], &gate).unwrap();
    let dense = surface.rows(&doubled, &[4000.0], &gate).unwrap();
    let warm = surface.rows(&warmer, &[4000.0], &gate).unwrap();
    assert_eq!(ordinary.sound_speed, doubled.sound_speed);
    for (actual, original) in dense[0].iter().zip(&base[0]) {
        let expected = original.scale(2.0);
        assert!(det::sqrt((*actual - expected).norm_sq() / expected.norm_sq()) < 1e-14);
    }
    assert_ne!(base, warm);
    // Lower propagation phase in the warmer material, independently of density.
    assert!(warm[0][1].re / warm[0][1].im < base[0][1].re / base[0][1].im);
}

#[test]
fn original_surface_and_work_admissions_refuse_overlaps_unresolved_shapes_and_excess_memory() {
    let (mouth, terminals, options) = fixture();
    let gas = gas(293.15, 101_325.0);
    let gate = CancelGate::new_clock_free();
    let overlap = DrivenMouthSurface { center_m: [0.0; 2], ..mouth };
    assert!(MixedSurface::new(overlap, &terminals, &gas,
        [1.0, 4000.0], options, &gate).is_err());
    let coarse = DrivenMouthSurface { cells: [1; 2], ..mouth };
    assert!(matches!(MixedSurface::new(coarse, &terminals, &gas,
        [1.0, 4000.0], options, &gate),
        Err(BaffledMimoBakeError::Invalid(what)) if what.contains("static")));
    assert!(matches!(MixedSurface::new(coarse, &terminals, &gas,
        [1.0, 96_000.0], options, &gate), Err(BaffledMimoBakeError::Surface(_))));
    let mut tiny_memory = options;
    tiny_memory.max_fit_memory_bytes = 1;
    assert!(matches!(MixedSurface::new(mouth, &terminals, &gas,
        [1.0, 4000.0], tiny_memory, &gate), Err(BaffledMimoBakeError::Memory { .. })));
    let expensive = DrivenMouthSurface { cells: [128; 2], ..mouth };
    assert!(MixedSurface::new(expensive, &terminals, &gas,
        [1.0, 4000.0], options, &gate).is_err());
}

mod multiple;
