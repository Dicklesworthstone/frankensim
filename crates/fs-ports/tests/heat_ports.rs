//! G0: a heat port carries entropy `dt q / T` into its owner, with the sign
//! of the heat `q` it carries at temperature `T`; power and stream ports
//! carry none on the wire.

use core::num::NonZeroUsize;

use fs_ports::component::{ContractViolation, PortLayout, PortRole, PortSpec, Receipt, audit_step};
use fs_ports::{PortKind, StableId};

/// A small generator of the same numbers on every build (splitmix64), as
/// fs-system's load-line tests draw theirs.
struct Draw(u64);

impl Draw {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn between(&mut self, low: f64, high: f64) -> f64 {
        low + (high - low) * (self.next() >> 11) as f64 / 9_007_199_254_740_992.0
    }
}

#[test]
fn a_heat_port_carries_entropy_with_the_sign_of_its_heat() {
    let mut draw = Draw(20_261_115);
    for role in [PortRole::EffortIn, PortRole::FlowIn] {
        let port = PortSpec::heat(StableId::new("heat").expect("id"), role).expect("heat port");
        // The temperature is the effort and the heat rate the flow.
        let sides = |kelvin: f64, watts: f64| {
            if role == PortRole::EffortIn {
                ([kelvin], [watts])
            } else {
                ([watts], [kelvin])
            }
        };
        for case in 0..400 {
            let kelvin = draw.between(1.0, 3000.0);
            let watts = draw.between(-1.0e3, 1.0e3);
            let dt = draw.between(1.0e-6, 1.0);
            let (input, output) = sides(kelvin, watts);
            let entropy = port.entropy_in_j_per_k(dt, &input, &output);
            let energy = port.energy_in_j(dt, &input, &output);
            // Heat in brings entropy in; heat out takes it out.
            assert_eq!(
                (entropy > 0.0, entropy < 0.0),
                (watts > 0.0, watts < 0.0),
                "{role:?} case {case}: {entropy:e} J/K for {watts} W at {kelvin} K"
            );
            // Its size is the heat that crossed over the temperature it
            // crossed at. The heat dt q is rounded alike on both sides; the
            // port's division and the product with T round once each, at
            // most eps of the heat together, and 2 eps holds that.
            assert!(
                (entropy * kelvin - energy).abs() <= 2.0 * f64::EPSILON * energy.abs(),
                "{role:?} case {case}: {entropy:e} J/K at {kelvin} K against {energy:e} J"
            );
        }
        // No heat, no entropy, and not a rounding of it.
        let (input, output) = sides(300.0, 0.0);
        assert_eq!(
            port.entropy_in_j_per_k(0.01, &input, &output).to_bits(),
            0.0_f64.to_bits()
        );
        // At a temperature that is not above zero the entropy is NaN, which
        // every audit refuses.
        for kelvin in [0.0, -0.0, -300.0, f64::NAN] {
            let (input, output) = sides(kelvin, 5.0);
            assert!(
                port.entropy_in_j_per_k(0.01, &input, &output).is_nan(),
                "{kelvin} K"
            );
        }
    }
    // That NaN is refused where the step's audit sums it: a store taking
    // 5 W at 0 K fails, and the same step at 300 K passes with dt q / T in.
    let store = [
        PortSpec::heat(StableId::new("heat").expect("id"), PortRole::FlowIn).expect("heat port"),
    ];
    let layout = PortLayout::new(&store).expect("layout");
    let (dt, watts) = (0.01, 5.0);
    let received = Receipt::lossless(dt * watts, 0.0);
    let audit = |kelvin: f64, receipt: &Receipt| {
        audit_step(
            &store,
            &layout,
            dt,
            &[watts],
            &[kelvin],
            receipt,
            0.0,
            dt * watts,
        )
    };
    assert!(matches!(
        audit(0.0, &received),
        Err(ContractViolation::NonFinite {
            quantity: "port entropy"
        })
    ));
    let passed = audit(300.0, &received).expect("audited at 300 K");
    assert_eq!(
        passed.port_entropy_in_j_per_k.to_bits(),
        (dt * watts / 300.0).to_bits()
    );
    // A receipt's own NaN entropy is refused before the ports are summed.
    let mut unknown = received;
    unknown.entropy_generated_j_per_k = f64::NAN;
    assert!(matches!(
        audit(300.0, &unknown),
        Err(ContractViolation::NonFinite {
            quantity: "entropy_generated_j_per_k"
        })
    ));

    // Other carriers carry none, whatever crosses them.
    let power = PortSpec::power(
        StableId::new("p").expect("id"),
        PortKind::MechanicalForceVelocity,
        PortRole::EffortIn,
        NonZeroUsize::MIN,
    )
    .expect("power port");
    assert_eq!(
        power.entropy_in_j_per_k(0.01, &[3.0], &[2.0]).to_bits(),
        0.0_f64.to_bits()
    );
    let stream = PortSpec::stream(StableId::new("s").expect("id"), PortRole::FlowIn, 0)
        .expect("stream port");
    let crossing = [1.0e-3, 300.0];
    let contents = [1.0e5, 3.0e5];
    assert_eq!(
        stream
            .entropy_in_j_per_k(0.01, &crossing, &contents)
            .to_bits(),
        0.0_f64.to_bits()
    );
}
