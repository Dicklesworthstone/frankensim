use super::*;
use fs_duct::compose::ComposeOptions;
use fs_duct::molecular::MolecularWallResolution;
use fs_duct::{ConeSlices, Duct, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};
use fs_thermochem::SpeciesId;

fn gas(temperature_k: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature_k, 101325.0).unwrap()
}

fn point(rate_scale: f64, wall: bool, exterior: GasState) -> MolecularDrivingPoint {
    // Authored synthetic data exercise the original input and constitutive
    // owners; these are deliberately not asserted to be measured air constants.
    let data = format!(
        "declared 1 synthetic 2 vibration 1 1300 1 power {} 293.15 0 200 500 \
         vibration 1 1900 1 power {} 293.15 0 200 500",
        101325.0 * 2e-5 * rate_scale,
        101325.0 * 3.2e-4 * rate_scale,
    );
    let state = read_relaxation(&mut data.split_ascii_whitespace())
        .unwrap()
        .evaluate(
            &RelaxingMedium::new(
                gas(293.15),
                vec![(SpeciesId::new("synthetic").unwrap(), 1.0)],
            )
            .unwrap(),
        )
        .unwrap();
    MolecularDrivingPoint::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.001,
                length: 0.15,
            }],
        },
        state,
        wall.then_some(MolecularWallResolution {
            poles_per_field: 16,
            band_rad_s: [1e-4, 30000.0],
            maximum_relative_error: 0.05,
        }),
        ConeSlices::Fixed(8),
        ComposeOptions::default(),
        exterior,
    )
    .unwrap()
}

#[test]
fn molecular_outlet_matches_independent_loaded_cylinder_solution() {
    let model = point(1.0, false, gas(320.0));
    let radius = model.inlet_radius_m();
    let area = core::f64::consts::PI * radius * radius;
    let gas = model.inlet_gas();
    let zc = gas.density * gas.sound_speed / area;
    for omega in [100.0, 1000.0, 6000.0, 24000.0] {
        let k = model.material().state().relaxing_wavenumber(omega).unwrap();
        let phase = k.scale(0.15);
        let cosine = C64::new(
            phase.re.cos() * phase.im.cosh(),
            -phase.re.sin() * phase.im.sinh(),
        );
        let sine = C64::new(
            phase.re.sin() * phase.im.cosh(),
            phase.re.cos() * phase.im.sinh(),
        );
        let wave_z = C64::from_re(omega * gas.density / area) / k;
        let load = fs_duct::compact_radiation_impedance(
            Termination::FlangedOpen,
            model.exterior_gas(),
            radius,
            omega,
        )
        .unwrap();
        // Direct uniform-cylinder solution of [Pin,Uin] = ABCD [Zload,1] Uout
        // and p+ = (Pin + Zc Uin)/2, independent of chain/scattering composition.
        let mixed = wave_z + load.scale(zc) / wave_z;
        let denominator = (load + C64::from_re(zc)) * cosine + C64::new(mixed.im, -mixed.re) * sine;
        let expected = C64::from_re(2.0) / denominator;
        let actual = model.flanged_outlet_flow(omega).unwrap();
        assert!(
            (actual - expected).abs() < 2e-11 * expected.abs(),
            "omega={omega}: {actual:?} versus {expected:?}"
        );
    }
}

#[test]
fn changing_only_reservoir_data_changes_physical_pressure_after_the_frozen_wavefront() {
    let baseline = point(1.0, true, gas(320.0));
    let changed = point(0.1, true, gas(320.0));
    let receiver = CircularOutletReceiver {
        position_m: [0.0, 0.0, 1.0],
        radial_rings: 2,
        angular_points: 8,
        maximum_frequency_hz: 800.0,
    };
    let mut a = baseline.baffled_outlet(8000, 1024, receiver).unwrap();
    let mut b = changed.baffled_outlet(8000, 1024, receiver).unwrap();
    assert_eq!(a.info().gas, gas(320.0));
    assert_eq!(
        a.info().axial_delay_samples,
        (baseline.one_way_wavefront_s() * 8000.0).floor() as usize
    );
    assert_eq!(a.info().axial_delay_samples, b.info().axial_delay_samples);
    let first = a.info().axial_delay_samples + a.info().exterior_delay_samples.0;
    let (mut energy, mut difference) = (0.0, 0.0);
    for n in 0..2000 {
        let input = if n == 0 { 100.0 } else { 0.0 };
        let pa = a.step(input).unwrap();
        let pb = b.step(input).unwrap();
        if n < first {
            assert_eq!(pa, 0.0);
            assert_eq!(pb, 0.0);
        }
        energy += pa * pa;
        difference += (pa - pb).powi(2);
    }
    assert!(
        energy > 1e-12 && difference > 1e-6 * energy,
        "energy={energy}, difference={difference}"
    );
}
