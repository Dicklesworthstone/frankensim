//! Diagnostic probe (E3.4b work aid): logs the generated per-capsule
//! radii for the registered 1903 rotor at both resolutions. No
//! acceptance logic lives here — it exists so cover-certificate
//! regressions carry per-item evidence in CI logs.
//! Repro: cargo test -p fs-flyer --test capsprobe_tmp -- --nocapture

use fs_flyer::aircraft::wright_rotor_v1;
use fs_flyer::sweptevents::capsules_at_resolution;

#[test]
fn log_generated_radii() {
    let rotor = wright_rotor_v1();
    for n in [8usize, 12] {
        for c in capsules_at_resolution(&rotor, n) {
            println!(
                "{{\"suite\":\"fs-flyer-capsprobe\",\"n\":{n},\"r0\":{},\"r1\":{},\"rail\":{},\"radius\":{}}}",
                c.r0_over_r, c.r1_over_r, c.rail_f, c.radius_m
            );
        }
    }
}
