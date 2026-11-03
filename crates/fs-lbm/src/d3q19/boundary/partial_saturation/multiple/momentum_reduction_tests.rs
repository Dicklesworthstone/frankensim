use super::*;

#[test]
fn opposite_tractions_keep_a_weak_resolved_center_couple_and_its_work() {
    let weak = 2.0_f64.powi(-62);
    let samples = [
        ([-1.0, -1.0, 0.0], [1.0 / 64.0, 0.0, 0.0]),
        ([0.0, -1.0, 0.0], [weak, 0.0, 0.0]),
        ([1.0, -1.0, 0.0], [-1.0 / 64.0, 0.0, 0.0]),
        ([-1.0, 1.0, 0.0], [1.0 / 64.0, 0.0, 0.0]),
        ([0.0, 1.0, 0.0], [-weak, 0.0, 0.0]),
        ([1.0, 1.0, 0.0], [-1.0 / 64.0, 0.0, 0.0]),
    ];
    let angular_velocity = 1.0 / 64.0;
    for reverse in [false, true] {
        let mut exchange = Exchange::default();
        for index in 0..samples.len() {
            let (arm, impulse) = samples[if reverse {
                samples.len() - 1 - index
            } else {
                index
            }];
            let offset = [arm[0] + 8.0, arm[1] - 4.0, arm[2] + 2.0];
            let velocity = [-angular_velocity * arm[1], angular_velocity * arm[0], 0.0];
            exchange.add(0.0, impulse, offset, arm, velocity, 1.0);
        }
        let result = exchange.finish();
        assert_eq!(result.fluid_population_impulse, [0.0; 3]);
        assert_eq!(
            result.body_angular_impulse_about_center,
            [0.0, 0.0, -2.0 * weak]
        );
        assert_eq!(result.body_angular_impulse, [0.0, 0.0, -2.0 * weak]);
        assert_eq!(result.body_work, -2.0 * weak * angular_velocity);
    }
}
