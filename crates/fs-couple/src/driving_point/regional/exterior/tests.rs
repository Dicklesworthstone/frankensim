use super::*;
use fs_duct::{Duct, HoleState, Segment, input_impedance_load};
use fs_material::gas::GasSpec;
use fs_material::gas_species::{GasSpecies, resolve_gas_mixture};
use fs_qty::semantic::{Composition, CompositionBasis};
use crate::pcm_wav::baffled::CircularOutletReceiver;

fn air(pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap()
}
fn enriched() -> GasState {
    resolve_gas_mixture(&[GasSpecies::Nitrogen, GasSpecies::CarbonDioxide],
        &Composition::new(CompositionBasis::MoleFraction, vec![0.3, 0.7]).unwrap(),
        293.15, 100000.0).unwrap()
}
fn duct() -> Duct {
    Duct { segments: vec![Segment::Cylinder { radius: 0.001, length: 0.2 },
        Segment::Cylinder { radius: 0.001, length: 0.3 }] }
}
fn bits(r: DuctResponse) -> [u64; 9] {
    [r.omega, r.impedance.re, r.impedance.im, r.min_shear_number, r.mouth_ka,
        r.u_mouth_over_u_in.re, r.u_mouth_over_u_in.im,
        r.p_mouth_over_u_in.re, r.p_mouth_over_u_in.im].map(f64::to_bits)
}

#[test]
fn explicit_identical_exterior_preserves_uniform_and_regional_response_and_history() {
    let a = air(100000.0); let b = enriched();
    for gases in [[a, a], [a, b]] {
        let original = RegionalDuct::new(&duct(), &gases).unwrap();
        let explicit = original.clone().with_exterior_gas(gases[1]).unwrap();
        for termination in [Termination::Closed, Termination::IdealOpen,
            Termination::UnflangedOpen, Termination::FlangedOpen] {
            for frequency in [173.0, 719.0, 2031.0] {
                let omega = core::f64::consts::TAU * frequency;
                assert_eq!(bits(original.response(omega, LossModel::Bessel, termination, None).unwrap()),
                    bits(explicit.response(omega, LossModel::Bessel, termination, None).unwrap()));
            }
        }
        let mut x = original.characteristic_line(Termination::UnflangedOpen, 48000, 4801, None, None).unwrap();
        let mut y = explicit.characteristic_line(Termination::UnflangedOpen, 48000, 4801, None, None).unwrap();
        for n in 0..1600 {
            let input = if n % 17 == 0 { 1.0 } else { 0.0 };
            assert_eq!(x.push(input).unwrap().to_bits(), y.push(input).unwrap().to_bits());
        }
    }
}

#[test]
fn exterior_load_reaches_the_inlet_without_replacing_internal_propagation() {
    let geometry = duct(); let a = air(100000.0); let b = enriched();
    let baseline = RegionalDuct::new(&geometry, &[a, b]).unwrap();
    let actual = baseline.clone().with_exterior_gas(a).unwrap();
    assert_eq!(actual.segment_gases(), baseline.segment_gases());
    assert_eq!(actual.inlet_impedance(), baseline.inlet_impedance());
    assert_eq!(actual.one_way_time_s(), baseline.one_way_time_s());
    let omega = core::f64::consts::TAU * 731.0;
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        // Independent assembly: external radiation -> heavy-gas tail -> air head.
        let load = fs_duct::compact_radiation_impedance(Termination::UnflangedOpen, &a, 0.001, omega).unwrap();
        let tail = input_impedance_load(&Duct { segments: vec![geometry.segments[1].clone()] }, &b,
            omega, loss, load, None).unwrap();
        let head = input_impedance_load(&Duct { segments: vec![geometry.segments[0].clone()] }, &a,
            omega, loss, tail.impedance, None).unwrap();
        let r = actual.response(omega, loss, Termination::UnflangedOpen, None).unwrap();
        assert_eq!(r.impedance, head.impedance);
        assert_eq!(r.u_mouth_over_u_in, tail.u_mouth_over_u_in * head.u_mouth_over_u_in);
        assert_eq!(r.p_mouth_over_u_in, tail.p_mouth_over_u_in * head.u_mouth_over_u_in);
        assert_eq!(r.mouth_ka, omega / a.sound_speed * 0.001);
        assert_ne!(r.impedance, baseline.response(omega, loss, Termination::UnflangedOpen, None).unwrap().impedance);
        for condition in [Termination::Closed, Termination::IdealOpen] {
            assert_eq!(bits(actual.response(omega, loss, condition, None).unwrap()),
                bits(baseline.response(omega, loss, condition, None).unwrap()));
        }
    }
}

#[test]
fn unflanged_band_admission_uses_exterior_sound_speed_not_bore_sound_speed() {
    let a = air(100000.0); let b = enriched();
    let omega = 0.99 * a.sound_speed / 0.001;
    let slow_bore = RegionalDuct::new(&duct(), &[b, b]).unwrap();
    assert!(matches!(slow_bore.response(omega, LossModel::Bessel, Termination::UnflangedOpen, None),
        Err(DuctError::RadiationKaTooLarge { .. })));
    assert!(slow_bore.with_exterior_gas(a).unwrap()
        .response(omega, LossModel::Bessel, Termination::UnflangedOpen, None).is_ok());
    let fast_bore = RegionalDuct::new(&duct(), &[a, a]).unwrap().with_exterior_gas(b).unwrap();
    assert!(matches!(fast_bore.response(omega, LossModel::Bessel, Termination::UnflangedOpen, None),
        Err(DuctError::RadiationKaTooLarge { .. })));
}

#[test]
fn independent_exterior_changes_receiver_travel_and_pressure_but_not_axial_travel() {
    let a = air(100000.0); let b = enriched();
    let original = RegionalDuct::new(&duct(), &[a, b]).unwrap();
    let external = original.clone().with_exterior_gas(a).unwrap();
    let receiver = CircularOutletReceiver { position_m: [0.0, 0.0, 1.0],
        radial_rings: 2, angular_points: 8, maximum_frequency_hz: 2000.0 };
    let mut x = original.baffled_outlet(48000, 1024, receiver, None).unwrap();
    let mut y = external.baffled_outlet(48000, 1024, receiver, None).unwrap();
    assert_eq!(x.info().gas, b); assert_eq!(y.info().gas, a);
    assert_eq!(x.info().axial_delay_samples, y.info().axial_delay_samples);
    assert!(y.info().exterior_delay_samples.0 < x.info().exterior_delay_samples.0);
    let first = y.info().axial_delay_samples + y.info().exterior_delay_samples.0;
    let (mut energy, mut difference) = (0.0, 0.0);
    for n in 0..2000 {
        let input = if n == 0 { 100.0 } else { 0.0 };
        let px = x.step(input).unwrap(); let py = y.step(input).unwrap();
        if n < first { assert_eq!(py, 0.0); }
        assert!(px.is_finite() && py.is_finite());
        energy += py * py; difference += (px - py) * (px - py);
    }
    assert!(energy > 0.0 && difference > 1e-8 * energy);
}

#[test]
fn incoherent_and_pressure_mismatched_media_refuse() {
    let a = air(100000.0); let b = enriched();
    let baseline = RegionalDuct::new(&duct(), &[a, b]).unwrap();
    assert!(baseline.clone().with_exterior_gas(air(50000.0)).is_err());
    let mut bad = a; bad.density *= 1.5;
    assert!(baseline.clone().with_exterior_gas(bad).is_err());
    assert_eq!(*baseline.exterior_gas(), b);
    for state in [HoleState::Closed, HoleState::Open, HoleState::Vent(0.3)] {
        let geometry = Duct { segments: vec![duct().segments[0].clone(),
            Segment::ToneHole { hole_radius: 0.0004, chimney_height: 0.001,
                bore_radius: 0.001, state }, duct().segments[1].clone()] };
        let region = RegionalDuct::new(&geometry, &[a, a, b]).unwrap();
        assert!(region.clone().with_exterior_gas(a).is_ok());
        assert!(region.with_exterior_gas(b).is_ok());
    }
}

#[test]
fn g3_exterior_reaches_upstream_holes_with_a_closed_terminal_and_a_baked_load() {
    let a = air(100000.0); let b = enriched();
    let geometry = Duct { segments: vec![duct().segments[0].clone(),
        Segment::ToneHole { hole_radius: 0.0004, chimney_height: 0.003,
            bore_radius: 0.001, state: HoleState::Vent(0.4) },
        duct().segments[1].clone()] };
    let default = RegionalDuct::new(&geometry, &[a, a, b]).unwrap();
    let original = default.clone().with_exterior_gas(a).unwrap();
    let external = default.clone().with_exterior_gas(b).unwrap();
    let omega = core::f64::consts::TAU * 731.0;
    for termination in [Termination::Closed, Termination::FlangedOpen] {
        let first = original.response_with_wall_and_holes(omega, LossModel::Bessel,
            termination, WallModel::Rigid).unwrap();
        let second = external.response_with_wall_and_holes(omega, LossModel::Bessel,
            termination, WallModel::Rigid).unwrap();
        let default_response = default.response_with_wall_and_holes(omega, LossModel::Bessel,
            termination, WallModel::Rigid).unwrap();
        assert_eq!(bits(default_response.0), bits(second.0));
        assert_eq!(default_response.1, second.1);
        assert_eq!(first.1.len(), 1); assert_eq!(second.1.len(), 1);
        assert_ne!(first.0.impedance, second.0.impedance);
        assert_ne!(first.1[0].u_hole_over_u_in, second.1[0].u_hole_over_u_in);
        assert_eq!(bits(second.0), bits(external.response(omega, LossModel::Bessel,
            termination, None).unwrap()));
        assert_eq!(original.one_way_time_s(), external.one_way_time_s());
    }
    // A baked terminal load does not replace the explicit side-hole medium.
    let z_load = fs_math::c64::C64::from_re(300000.0);
    let loaded = |region: &RegionalDuct| region.loaded_reflection(
        omega, z_load, WallModel::Rigid).unwrap();
    assert_ne!(loaded(&original), loaded(&external));
    assert_eq!(loaded(&default), loaded(&external));
}

#[test]
fn independent_exterior_reaches_uniform_holes_with_closed_or_ideal_open_terminal() {
    let a = air(100000.0); let b = enriched();
    let geometry = Duct { segments: vec![duct().segments[0].clone(),
        Segment::ToneHole { hole_radius: 0.0004, chimney_height: 0.003,
            bore_radius: 0.001, state: HoleState::Open },
        duct().segments[1].clone()] };
    let original = RegionalDuct::new(&geometry, &[a, a, a]).unwrap();
    let external = original.clone().with_exterior_gas(b).unwrap();
    for termination in [Termination::Closed, Termination::IdealOpen] {
        let omega = core::f64::consts::TAU * 731.0;
        assert_ne!(original.response(omega, LossModel::Bessel, termination, None).unwrap().impedance,
            external.response(omega, LossModel::Bessel, termination, None).unwrap().impedance);
        let mut first = original.characteristic_line(termination, 8000, 2048, None, Some(512)).unwrap();
        let mut second = external.characteristic_line(termination, 8000, 2048, None, Some(512)).unwrap();
        let mut difference = 0.0;
        for sample in 0..2048 {
            let input = if sample == 0 { 1.0 } else { 0.0 };
            let a = first.push(input).unwrap();
            let b = second.push(input).unwrap();
            assert!(a.is_finite() && b.is_finite());
            difference += (a - b).powi(2);
        }
        assert!(difference > 0.0, "terminal condition must not bypass side-hole loading");
    }
}
