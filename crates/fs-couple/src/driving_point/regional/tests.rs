use super::*;
use fs_material::gas::{GasSpec, GasState};
use fs_material::gas_species::{GasSpecies, resolve_gas_mixture};
use fs_math::c64::C64;
use fs_qty::semantic::{Composition, CompositionBasis};

fn air() -> GasState { GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 100000.0).unwrap() }
fn co2() -> GasState {
    resolve_gas_mixture(&[GasSpecies::Nitrogen, GasSpecies::CarbonDioxide],
        &Composition::new(CompositionBasis::MoleFraction, vec![0.3, 0.7]).unwrap(),
        293.15, 100000.0).unwrap()
}
fn cylinder(length: f64) -> Segment { Segment::Cylinder { radius: 0.001, length } }
fn hole(state: HoleState) -> Segment {
    Segment::ToneHole { hole_radius: 0.0004, chimney_height: 0.001, bore_radius: 0.001, state }
}
fn close(a: C64, b: C64) { assert!((a-b).abs() < 2e-10*b.abs().max(1.0), "{a:?} != {b:?}"); }
fn bits(r: DuctResponse) -> [u64; 9] {
    [r.omega, r.impedance.re, r.impedance.im, r.min_shear_number, r.mouth_ka,
        r.u_mouth_over_u_in.re, r.u_mouth_over_u_in.im,
        r.p_mouth_over_u_in.re, r.p_mouth_over_u_in.im].map(f64::to_bits)
}

#[test]
fn uniform_binding_keeps_complete_tmm_and_fir_bit_identical() {
    let duct = Duct { segments: vec![cylinder(0.12), hole(HoleState::Open),
        cylinder(0.1), hole(HoleState::Vent(0.3)),
        Segment::Cone { inlet_radius: 0.001, outlet_radius: 0.0008, length: 0.15 }] };
    let gas = air();
    let regional = RegionalDuct::new(&duct, &vec![gas; duct.segments.len()]).unwrap();
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        for termination in [Termination::Closed, Termination::IdealOpen, Termination::UnflangedOpen] {
            for f in [173.0, 389.0, 997.0] {
                let omega = core::f64::consts::TAU*f;
                assert_eq!(bits(regional.response(omega, loss, termination, None).unwrap()),
                    bits(input_impedance_wall(&duct, &gas, omega, loss, termination, None).unwrap()));
            }
        }
    }
    let mut a = super::super::characteristic_line(&duct, &gas, Termination::UnflangedOpen,
        48000, 4801, regional.inlet_impedance(), None).unwrap();
    let mut b = regional.characteristic_line(Termination::UnflangedOpen, 48000, 4801, None, None).unwrap();
    for n in 0..8192 {
        let input = if n%37 == 0 { 1.0 } else { 0.0 };
        assert_eq!(a.incoming().to_bits(), b.incoming().to_bits());
        assert_eq!(a.push(input).unwrap().to_bits(), b.push(input).unwrap().to_bits());
    }
}

// Independent closed-form lossless cylinder ABCD oracle; production delegates
// every region's propagation and load to fs-duct's numerical basis construction.
fn multiply(a: [C64;4], b: [C64;4]) -> [C64;4] {
    [a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
        a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3]]
}
#[test]
fn interfaces_match_independent_three_medium_transfer_matrix() {
    let lengths = [0.12, 0.17, 0.09];
    let gases = [air(), co2(), GasState::try_new(&GasSpec::dry_air_ussa1976(), 320.0, 100000.0).unwrap()];
    let duct = Duct { segments: lengths.map(cylinder).to_vec() };
    let regional = RegionalDuct::new(&duct, &gases).unwrap();
    let time: f64 = lengths.iter().zip(gases).map(|(l,g)| l/g.sound_speed).sum();
    assert_eq!(regional.one_way_time_s(), time);
    for i in 0..80 {
        let omega = core::f64::consts::TAU*(123.0+13.1*i as f64);
        let mut m = [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE];
        for (l,g) in lengths.iter().zip(gases) {
            let phase = omega*l/g.sound_speed;
            let z = g.density*g.sound_speed/(core::f64::consts::PI*1e-6);
            m = multiply(m, [C64::from_re(phase.cos()), C64::new(0.0,-z*phase.sin()),
                C64::new(0.0,-phase.sin()/z), C64::from_re(phase.cos())]);
        }
        let open = regional.response(omega, LossModel::Lossless, Termination::IdealOpen, None).unwrap();
        close(open.impedance, m[1]/m[3]);
        close(open.u_mouth_over_u_in, m[3].recip());
        close(open.p_mouth_over_u_in, C64::ZERO);
        let closed = regional.response(omega, LossModel::Lossless, Termination::Closed, None).unwrap();
        close(closed.impedance, m[0]/m[2]);
        close(closed.u_mouth_over_u_in, C64::ZERO);
    }
}

#[test]
fn terminal_radiation_and_viscothermal_losses_use_their_actual_media() {
    let duct = Duct { segments: vec![cylinder(0.13), cylinder(0.21)] };
    let gases = [air(), co2()];
    let regional = RegionalDuct::new(&duct, &gases).unwrap();
    let omega = core::f64::consts::TAU*731.0;
    let tail = input_impedance_wall(&Duct { segments: vec![duct.segments[1].clone()] }, &gases[1],
        omega, LossModel::Bessel, Termination::UnflangedOpen, None).unwrap();
    let head = input_impedance_load(&Duct { segments: vec![duct.segments[0].clone()] }, &gases[0],
        omega, LossModel::Bessel, tail.impedance, None).unwrap();
    let result = regional.response(omega, LossModel::Bessel, Termination::UnflangedOpen, None).unwrap();
    assert_eq!(result.mouth_ka.to_bits(), tail.mouth_ka.to_bits());
    close(result.impedance, head.impedance);
    close(result.u_mouth_over_u_in, tail.u_mouth_over_u_in*head.u_mouth_over_u_in);
    close(result.p_mouth_over_u_in, tail.p_mouth_over_u_in*head.u_mouth_over_u_in);
    assert_eq!(result.min_shear_number, head.min_shear_number.min(tail.min_shear_number));
    let mut actual = regional.characteristic_line(Termination::UnflangedOpen, 48000, 4801, None, None).unwrap();
    let mut wrong = super::super::characteristic_line(&duct, &gases[0], Termination::UnflangedOpen,
        48000, 4801, regional.inlet_impedance(), None).unwrap();
    let mut difference = 0.0;
    for n in 0..4096 {
        let input = if n == 0 { 1.0 } else { 0.0 };
        let a = actual.push(input).unwrap(); let b = wrong.push(input).unwrap();
        assert!(a.is_finite() && b.is_finite());
        difference += (a-b)*(a-b);
    }
    assert!(difference > 1e-8, "a downstream material change must reach the inlet feedback");
}

#[test]
fn incomplete_binding_pressure_jumps_and_unmodeled_mutual_holes_refuse() {
    let duct = Duct { segments: vec![cylinder(0.1), cylinder(0.2)] };
    let a = air(); let b = co2();
    assert!(RegionalDuct::new(&duct, &[a]).is_err());
    let mut bad = b; bad.pressure *= 0.5;
    assert!(RegionalDuct::new(&duct, &[a, bad]).is_err());
    let mut bad = b; bad.density *= 1.5;
    assert!(RegionalDuct::new(&duct, &[a, bad]).is_err());
    let low = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 50000.0).unwrap();
    assert!(RegionalDuct::new(&duct, &[a, low]).is_err());
    let holes = Duct { segments: vec![cylinder(0.1), hole(HoleState::Open), cylinder(0.1),
        hole(HoleState::Open), cylinder(0.1)] };
    assert!(RegionalDuct::new(&holes, &[a,a,b,b,b]).is_err());
    assert!(RegionalDuct::new(&holes, &[a,a,a,a,b]).is_ok());
    let too_many = Duct { segments: vec![cylinder(0.01); MAX_REGIONAL_SEGMENTS+1] };
    assert!(RegionalDuct::new(&too_many, &vec![a; MAX_REGIONAL_SEGMENTS+1]).is_err());
    let regional = RegionalDuct::new(&duct, &[a,b]).unwrap();
    assert!(regional.characteristic_line(Termination::Closed, 48000, 4, None, None).is_err());
    assert!(regional.response(f64::NAN, LossModel::Lossless, Termination::Closed, None).is_err());
    let omega = 1.01*b.sound_speed/0.001;
    assert!(matches!(regional.response(omega, LossModel::Bessel, Termination::UnflangedOpen, None),
        Err(DuctError::RadiationKaTooLarge { .. })));
}
