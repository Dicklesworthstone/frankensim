//! G0/G1/G3 conservative rod compression stability, independent of Newton's
//! ability to converge to a straight stationary configuration.

use fs_ga::{So3, So3Tangent, Vec3};
use fs_solid::{Rod, RodSection, TipLoad};

fn column(length: f64, segments: usize, force_scale: f64, compression: f64) -> (Rod, TipLoad) {
    let section = RodSection {
        ea: 1e5 * force_scale,
        ga: 1e5 * force_scale,
        gj: force_scale * length * length,
        ei: force_scale * length * length,
    };
    let mut rod = Rod::straight(length, segments, section);
    // This is the exact straight equilibrium of the original axial law.
    // A straight equilibrium exists both below and above the buckling load.
    let stretch = 1.0 - compression / section.ea;
    for position in &mut rod.positions {
        position[0] *= stretch;
    }
    (
        rod,
        TipLoad {
            force: [-compression, 0.0, 0.0],
            moment: [0.0; 3],
        },
    )
}

#[test]
fn g1_cantilever_stability_converges_to_euler_load() {
    // A clamped-free Euler column has Pcr = pi² EI / (4 L²).
    // Continuum Mechanics of Beam and Plate Flexure, Example 3.1:
    // https://vbn.aau.dk/ws/files/197084122/Continuum_Mechanics_of_Beam_and_Plate_Flexure_reduced_.pdf
    // Equal EA and GA cancel the axial-shortening/shear correction in this
    // Cosserat law's straight-state second variation. EI/L² is one here.
    // The geodesic-midpoint discretization must converge to that continuum
    // threshold, not the four-times-larger pinned-pinned threshold.
    let euler = std::f64::consts::PI.powi(2) / 4.0;
    let mut previous_error = f64::INFINITY;
    for segments in [4, 8, 16] {
        let mut low = 0.95 * euler;
        let mut high = 1.05 * euler;
        for _ in 0..20 {
            let pressure = 0.5 * (low + high);
            let (rod, load) = column(1.0, segments, 1.0, pressure);
            if rod.is_stable_under_dead_load(&load, 1e-7).unwrap() {
                low = pressure;
            } else {
                high = pressure;
            }
        }
        let critical = 0.5 * (low + high);
        let error = critical / euler - 1.0;
        assert!(
            error > 0.0 && error < 0.027,
            "segments={segments}, error={error:e}"
        );
        if previous_error.is_finite() {
            let reduction = previous_error / error;
            assert!(
                (3.8..4.2).contains(&reduction),
                "refinement ratio {reduction}"
            );
        }
        previous_error = error;
    }
}

#[test]
fn g3_stable_and_unstable_classification_survives_rescaling_and_rotation() {
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.4, -0.7, 0.2))).unwrap();
    for (length, force_scale) in [(1e-2, 1e-2), (1.0, 1.0), (100.0, 100.0)] {
        let euler = force_scale * std::f64::consts::PI.powi(2) / 4.0;
        for (factor, expected) in [(0.7, true), (1.3, false)] {
            let (mut rod, mut load) = column(length, 8, force_scale, factor * euler);
            let tolerance = 1e-7 * force_scale * length.max(1.0);
            assert_eq!(
                rod.is_stable_under_dead_load(&load, tolerance).unwrap(),
                expected
            );
            for (position, frame) in rod.positions.iter_mut().zip(&mut rod.frames) {
                let point = rotation
                    .rotate(Vec3::new(position[0], position[1], position[2]))
                    .unwrap();
                *position = [point.x, point.y, point.z];
                *frame = rotation.compose(*frame).unwrap();
            }
            let force = rotation
                .rotate(Vec3::new(load.force[0], load.force[1], load.force[2]))
                .unwrap();
            load.force = [force.x, force.y, force.z];
            assert_eq!(
                rod.is_stable_under_dead_load(&load, tolerance).unwrap(),
                expected
            );
        }
    }
}

#[test]
fn g0_stability_refuses_unconverged_states_couples_and_invalid_sections() {
    let (rod, load) = column(1.0, 8, 1.0, 1.0);
    assert!(rod.is_stable_under_dead_load(&load, 1e-7).unwrap());
    let wrong_force = TipLoad {
        force: [-1.1, 0.0, 0.0],
        ..load
    };
    assert!(rod.is_stable_under_dead_load(&wrong_force, 1e-7).is_err());
    let couple = TipLoad {
        moment: [0.0, 0.0, 1e-12],
        ..load
    };
    assert!(rod.is_stable_under_dead_load(&couple, 1e-7).is_err());
    for tolerance in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(rod.is_stable_under_dead_load(&load, tolerance).is_err());
    }
    let mut bad_section = rod.clone();
    bad_section.section.ei = 0.0;
    assert!(bad_section.is_stable_under_dead_load(&load, 1e-7).is_err());
    let mut bad_geometry = rod;
    bad_geometry.positions[3][0] = f64::NAN;
    assert!(bad_geometry.is_stable_under_dead_load(&load, 1e-7).is_err());
    let (too_large, load) = column(1.0, 65, 1.0, 1.0);
    assert!(too_large.is_stable_under_dead_load(&load, 1e-7).is_err());
}
