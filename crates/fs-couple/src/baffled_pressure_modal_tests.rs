//! Analytic radial-mode directivity and independent finite-distance quadrature.
use super::*;
use fs_material::gas::GasSpec;
use fs_math::c64::C64;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn spec(position_m: [f64; 3]) -> CircularOutletReceiver {
    CircularOutletReceiver { position_m, radial_rings: 128, angular_points: 256,
        maximum_frequency_hz: 6000.0 }
}
fn options() -> ModalOutletOptions {
    ModalOutletOptions { relative_kernel_tolerance: 3e-5, maximum_kernel_terms: 32_000_000 }
}

// Independent entire-function series; no production Bessel or modal shape call.
fn bessel(x: f64, order: usize) -> f64 {
    let mut term = if order == 0 {1.0} else {0.5*x};
    let mut sum = term;
    for m in 1..64 {
        term *= -0.25*x*x/(m*(m+order)) as f64;
        sum += term;
    }
    sum
}
const ROOTS: [f64; 3] = [0.0, 3.831_705_970_207_512_5, 7.015_586_669_815_619];

// Composite Simpson radial integration and a separate periodic angular grid.
// Integrate exp(ikR)/R directly; no production surface, FIR, delay or harmonic
// observer helper is used. Peak phasors follow exp(-i omega t).
fn reference(radius: f64, position: [f64; 3], frequency: f64, gas: &GasState) -> [C64; 3] {
    let (radial, angular) = (1024, 512);
    let omega = TAU*frequency;
    let mut pressure = [C64::ZERO; 3];
    for r in 0..=radial {
        let eta = r as f64/radial as f64;
        let weight = if r == 0 || r == radial {1.0} else if r % 2 == 0 {2.0} else {4.0};
        let mut green = C64::ZERO;
        for a in 0..angular {
            let angle = TAU*a as f64/angular as f64;
            let distance = (position[0]-radius*eta*angle.cos()).hypot(
                position[1]-radius*eta*angle.sin()).hypot(position[2]);
            let phase = omega*distance/gas.sound_speed;
            green = green + C64::new(phase.cos(), phase.sin()).scale(1.0/distance);
        }
        for (n, value) in pressure.iter_mut().enumerate() {
            let shape = bessel(ROOTS[n]*eta, 0)/bessel(ROOTS[n], 0).abs();
            *value = *value + green.scale(weight*eta*shape);
        }
    }
    pressure.map(|p| p*C64::new(0.0, -omega*gas.density/TAU
        *2.0/(3.0*radial as f64*angular as f64)))
}

#[test]
fn zero_net_flow_modes_reach_the_actual_time_receiver_with_coherent_near_field_phase() {
    let gas = air();
    let radius = 0.04;
    let rate = 48_000;
    let basis = AxisymmetricMouthBasis::new(ModeCount(3)).unwrap();
    let (surface, defect, _) = surface(radius, &basis, 32, 128).unwrap();
    assert!(defect < 1e-12);
    for n in 0..3 {
        let net: f64 = surface.iter().map(|p| p.area_m2*p.mode_shape[n]).sum();
        assert!((net-if n == 0 {1.0} else {0.0}).abs() < 2e-13);
    }
    // The source has ZERO total flow. Collapsing its modes to a piston or
    // discarding their signs would therefore fail this pressure comparison.
    let flows = [C64::ZERO, C64::new(1e-7, 0.3e-7), C64::new(-0.2e-7, 0.4e-7)];
    for position in [[0.0,0.0,0.06], [0.07,0.015,0.06]] {
        for frequency in [600.0, 4300.0] {
            let (mut receiver, report) = BaffledPressure::circular_modal_outlet_multirate_differentiating(
                radius, ModeCount(3), rate, rate, spec(position), &gas, options(), || false).unwrap();
            assert_eq!(report.radial_nodes, [32,64,128]);
            assert_eq!(report.angular_points, [64,128,256]);
            let transfers = reference(radius, position, frequency, &gas);
            let expected = transfers.iter().zip(flows)
                .fold(C64::ZERO, |sum, (&h,u)| sum+h*u);
            assert!(expected.abs() > 1e-8, "zero net flow still radiates at finite distance");
            let latency = receiver.filter_delay_samples();
            let warmup = receiver.delay_samples.1+256;
            let (mut error, mut energy) = (0.0, 0.0);
            for sample in 0..warmup+2048 {
                let phase = TAU*frequency*sample as f64/f64::from(rate);
                let frame = flows.map(|u| u.re*phase.cos()+u.im*phase.sin());
                let actual = receiver.step(&frame).unwrap();
                if sample >= warmup {
                    let phase = TAU*frequency*(sample as f64-latency)/f64::from(rate);
                    let target = expected.re*phase.cos()+expected.im*phase.sin();
                    error += (actual-target).powi(2);
                    energy += target*target;
                }
            }
            assert!((error/energy).sqrt() < 8e-4,
                "position={position:?}, frequency={frequency}, relative error={}", (error/energy).sqrt());
        }
    }
}

#[test]
fn original_signed_radial_modes_match_the_analytic_far_field_lobes() {
    // D_n(q)=2 sign(J0(gamma_n)) q J1(q)/(q²-gamma_n²), q=ka sin(theta).
    // Its removable value at q=gamma_n is |J0(gamma_n)|; at q=0 all
    // higher modes vanish only in the FAR-FIELD limit, not at finite range.
    let gas = air();
    let (radius, distance, ka) = (0.04, 40.0, 8.0);
    let frequency = ka*gas.sound_speed/(TAU*radius);
    let basis = AxisymmetricMouthBasis::new(ModeCount(3)).unwrap();
    let (samples, _, _) = surface(radius, &basis, 64, 256).unwrap();
    for q in [0.0, 1.7, ROOTS[1], 5.2, 7.5] {
        let sine = q/ka;
        let position = [distance*sine, 0.0, distance*(1.0-sine*sine).sqrt()];
        let field = super::super::harmonic::Prepared::new(&samples, 3, &[position], frequency,
            RayleighMedium::from(&gas), super::super::harmonic::Options::default(), || false).unwrap();
        let phase = -ka*distance/radius;
        let remove_propagation = C64::new(phase.cos(), phase.sin())
            *C64::new(0.0, distance/(frequency*gas.density));
        for (n, &gamma) in ROOTS.iter().enumerate() {
            let exact = if q == 0.0 {if n == 0 {1.0} else {0.0}}
                else if q == gamma {bessel(gamma, 0).abs()}
                else {2.0*bessel(gamma, 0).signum()*q*bessel(q, 1)/(q*q-gamma*gamma)};
            let actual = field.transfer(0,n).unwrap().pressure_pa*remove_propagation;
            assert!((actual-C64::from_re(exact)).abs() < 5e-3,
                "mode={n}, q={q}, actual={actual:?}, far field={exact}");
            if n == 0 && q == 5.2 { assert!(actual.re < 0.0, "phase reversal across a piston lobe"); }
        }
    }
}

#[test]
fn full_axis_work_and_unresolved_modes_refuse_without_a_running_receiver() {
    let gas = air();
    let point = spec([0.07,0.015,0.06]);
    let build = |spec, options| BaffledPressure::circular_modal_outlet_multirate_differentiating(
        0.04, ModeCount(3), 48_000, 48_000, spec, &gas, options, || false);
    let (_, report) = build(point, options()).unwrap();
    let at_cap = ModalOutletOptions { maximum_kernel_terms: report.kernel_terms, ..options() };
    assert!(build(point, at_cap).is_ok());
    assert!(build(point, ModalOutletOptions { maximum_kernel_terms: report.kernel_terms-1, ..options() }).is_err());
    assert!(build(CircularOutletReceiver { radial_rings: 4, ..point }, options()).is_err());
    assert!(build(CircularOutletReceiver { angular_points: 33, ..point }, options()).is_err());
    assert!(build(point, ModalOutletOptions { relative_kernel_tolerance: 0.0, ..options() }).is_err());
    assert!(BaffledPressure::circular_modal_outlet_multirate_differentiating(
        0.04, ModeCount(3), 48_000, 48_000, point, &gas, options(), || true).is_err());
}
