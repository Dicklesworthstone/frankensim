use super::*;
use crate::render::plate::impact::squeeze::sliding::{ShearedFilm, WallVelocityPorts};
use crate::render::plate::impact::squeeze::{
    FilmCell, FilmChannel, FilmLimits, GapPort, ResistiveFilm,
};
use crate::render::plate::impact::{ImpactBody, ImpactConfig};
use fs_exec::CancelGate;
use fs_tribo::resistive_film::isothermal::AmbientGas;

const M1: f64 = 0.2;
const M2: f64 = 0.4;
const GAP: f64 = 0.001;
const WIDTH: f64 = 0.02;
const LENGTH: f64 = 0.03;
const AREA: f64 = WIDTH * LENGTH;
const K: f64 = 0.05;
const DT: f64 = 0.0001;
const LAUNCH: f64 = 0.01;

fn bare() -> ImpactSystem {
    let (a, _) = ImpactBody::free_mass(M1, 0.0, LAUNCH).unwrap();
    let (b, _) = ImpactBody::free_mass(M2, 0.0, 0.0).unwrap();
    ImpactSystem::new(
        vec![a, b],
        vec![],
        vec![],
        vec![],
        ImpactConfig {
            dt_s: DT,
            max_steps: 128,
            maximum_energy_j: 1.0,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-8,
            maximum_generalized_force: 1.0,
        },
    )
    .unwrap()
}
fn graph(minimum_gap: f64) -> ResistiveFilm {
    let gap = GapPort {
        reference_m: GAP,
        closure: vec![1.0 / M1.sqrt(), -1.0 / M2.sqrt()],
    };
    ResistiveFilm::new(
        vec![FilmCell {
            area_m2: AREA,
            gap: gap.clone(),
        }],
        vec![FilmChannel {
            from: 0,
            to: None,
            width_m: WIDTH,
            length_m: LENGTH,
            gap,
        }],
        2,
        K,
        FilmLimits {
            minimum_cell_gap_m: minimum_gap,
            maximum_gap_m: 0.01,
            maximum_pressure_pa: 1e5,
        },
    )
    .unwrap()
}
fn film(n: f64, minimum_gap: f64) -> PowerLawFilm {
    PowerLawFilm::new(
        &graph(minimum_gap),
        vec![RegularizedPowerLaw {
            consistency_pa_s_n: K,
            flow_index: n,
            minimum_shear_rate_per_s: 1.0,
        }],
        PowerLawControls {
            maximum_iterations: 64,
            maximum_line_search_steps: 32,
            absolute_flow_tolerance_m3_s: 1e-15,
            relative_flow_tolerance: 1e-10,
            absolute_power_tolerance_w: 1e-13,
            relative_power_tolerance: 1e-10,
        },
    )
    .unwrap()
}
fn pair() -> ImpactSystem {
    bare().with_power_law_film(film(0.5, 1e-6)).unwrap()
}
fn momentum(x: &[f64]) -> f64 {
    M1.sqrt() * x[1] + M2.sqrt() * x[3]
}
fn relative(x: &[f64]) -> f64 {
    x[1] / M1.sqrt() - x[3] / M2.sqrt()
}
fn closure(x: &[f64]) -> f64 {
    x[0] / M1.sqrt() - x[2] / M2.sqrt()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}

// Independent stress-profile quadrature, without the production flux formula.
// For n=1/2, gamma(tau) is linear below tau0=K and quadratic above it.
// Splitting at that physical core transition makes Simpson exact apart from
// floating-point summation: y*gamma is quadratic or cubic on each interval.
fn profile_flow(pressure: f64, gap: f64) -> f64 {
    if pressure == 0.0 {
        return 0.0;
    }
    let half = 0.5 * gap;
    let core = (K * LENGTH / pressure).min(half);
    let integrand = |y: f64| {
        let stress = pressure * y / LENGTH;
        let ratio = stress / K;
        let shear = if stress <= K { ratio } else { ratio * ratio };
        y * shear
    };
    let integrate = |lo: f64, hi: f64| {
        let step = (hi - lo) / 16.0;
        let mut sum = integrand(lo) + integrand(hi);
        for i in 1..16 {
            sum += if i % 2 == 0 { 2.0 } else { 4.0 } * integrand(lo + i as f64 * step);
        }
        step * sum / 3.0
    };
    2.0 * WIDTH * (integrate(0.0, core) + integrate(core, half))
}

// One scalar discrete-gradient equation, eliminated using conserved total
// momentum and the two independently known free masses. Pressure is recovered
// by bisection of actual profile volume flux, not a native-film evaluation.
fn reference_tick(displacement: f64, before: f64) -> (f64, f64, f64, f64) {
    let inverse_mass = 1.0 / M1 + 1.0 / M2;
    let mut lo = 0.0;
    let mut hi = before / (DT * inverse_mass * AREA);
    for _ in 0..72 {
        let p = 0.5 * (lo + hi);
        let after = before - DT * inverse_mass * AREA * p;
        let midpoint_velocity = 0.5 * (before + after);
        let midpoint_gap = GAP - displacement - 0.5 * DT * midpoint_velocity;
        if profile_flow(p, midpoint_gap) > AREA * midpoint_velocity {
            hi = p;
        } else {
            lo = p;
        }
    }
    let pressure = 0.5 * (lo + hi);
    let after = before - DT * inverse_mass * AREA * pressure;
    let midpoint_velocity = 0.5 * (before + after);
    let next_displacement = displacement + DT * midpoint_velocity;
    let heat = DT * AREA * pressure * midpoint_velocity;
    (next_displacement, after, heat, pressure)
}

#[test]
fn g1_above_core_thinning_drives_joint_motion_and_actual_profile_work() {
    let mut system = pair();
    assert_eq!(system.power_law_film_count(), 1);
    assert_eq!(system.squeeze_film_count(), 0);
    assert_eq!(system.sliding_film_count(), 0);
    assert_eq!(system.power_law_film(0).unwrap().laws()[0].flow_index, 0.5);
    assert!(system.power_law_film(1).is_none());
    let initial = system.stored_energy_j();
    let p0 = momentum(system.state());
    let mut displacement = 0.0;
    let mut rate = LAUNCH;
    let mut loss = 0.0;
    for tick in 1..=32 {
        let (next, velocity, heat, pressure) = reference_tick(displacement, rate);
        // Every step materially exceeds the physical core transition.
        let midpoint_gap = GAP - 0.5 * (displacement + next);
        assert!(pressure * midpoint_gap / (2.0 * LENGTH) > 4.0 * K);
        let frame = system
            .step(&[0.0; 2], &CancelGate::new_clock_free())
            .unwrap();
        near(closure(system.state()), next, 2e-11);
        near(relative(system.state()), velocity, 2e-10);
        near(momentum(system.state()), p0, 2e-11);
        near(frame.dissipated_energy_j, heat, 2e-12);
        assert_eq!(frame.sample, tick);
        near(frame.time_s, tick as f64 * DT, 1e-15);
        displacement = next;
        rate = velocity;
        loss += frame.dissipated_energy_j;
    }
    assert!(system.state()[3] > 0.0);
    near(system.stored_energy_j() + loss, initial, 2e-11);
}

#[test]
fn g3_complete_force_and_tangent_cover_the_original_mechanical_prefix() {
    let system = pair();
    let mut x = system.state().to_vec();
    x.extend([0.3, -0.2]);
    let mut effort = system.state().to_vec();
    effort.extend([0.7, -0.8]);
    let mut force = [99.0; 6];
    assert!(system.dissipative_flow_into(&x, &effort, &mut force));
    assert_eq!(force[0], 0.0);
    assert_eq!(force[2], 0.0);
    assert_eq!(&force[4..], &[0.0; 2]);
    near(M1.sqrt() * force[1] + M2.sqrt() * force[3], 0.0, 1e-15);

    let dx = [1e-4 * M1.sqrt(), 0.0, -5e-5 * M2.sqrt(), 0.0, 0.0, 0.0];
    let de = [0.0, 1e-3 * M1.sqrt(), 0.0, -2e-3 * M2.sqrt(), 0.0, 0.0];
    let mut tangent = [99.0; 6];
    assert!(system.dissipative_flow_tangent_into(&x, &effort, &dx, &de, &mut tangent));
    let epsilon = 1e-5;
    let xp: Vec<_> = x.iter().zip(dx).map(|(a, d)| a + epsilon * d).collect();
    let xm: Vec<_> = x.iter().zip(dx).map(|(a, d)| a - epsilon * d).collect();
    let ep: Vec<_> = effort
        .iter()
        .zip(de)
        .map(|(a, d)| a + epsilon * d)
        .collect();
    let em: Vec<_> = effort
        .iter()
        .zip(de)
        .map(|(a, d)| a - epsilon * d)
        .collect();
    let mut plus = [0.0; 6];
    let mut minus = [0.0; 6];
    assert!(system.dissipative_flow_into(&xp, &ep, &mut plus));
    assert!(system.dissipative_flow_into(&xm, &em, &mut minus));
    for i in 0..6 {
        let finite_difference = (plus[i] - minus[i]) / (2.0 * epsilon);
        near(
            tangent[i],
            finite_difference,
            2e-8 * finite_difference.abs() + 1e-10,
        );
    }
    assert_eq!(&tangent[4..], &[0.0; 2]);
}

#[test]
fn g3_prepared_images_retry_atomically_and_keep_the_exact_newtonian_limit() {
    let mut reference = pair();
    for _ in 0..32 {
        reference
            .step(&[0.0; 2], &CancelGate::new_clock_free())
            .unwrap();
    }
    for analytic in [false, true] {
        let mut system = if analytic {
            pair().prepare_analytic().unwrap()
        } else {
            pair().prepare().unwrap()
        };
        let initial = system.stored_energy_j();
        let mut loss = 0.0;
        for tick in 0..32 {
            if tick == 8 {
                let before = system.state().to_vec();
                let samples = system.samples();
                let cancelled = CancelGate::new_clock_free();
                cancelled.request();
                assert!(system.step(&[0.0; 2], &cancelled).is_err());
                assert!(
                    system
                        .step(&[2.0, 0.0], &CancelGate::new_clock_free())
                        .is_err()
                );
                assert_eq!(system.state(), before);
                assert_eq!(system.samples(), samples);
            }
            loss += system
                .step(&[0.0; 2], &CancelGate::new_clock_free())
                .unwrap()
                .dissipated_energy_j;
        }
        for (&a, &b) in system.state().iter().zip(reference.state()) {
            near(a, b, 2e-10);
        }
        near(system.stored_energy_j() + loss, initial, 2e-11);

        let mut original = bare()
            .with_squeeze_film(graph(1e-6))
            .unwrap()
            .prepare()
            .unwrap();
        let mut limit = bare()
            .with_power_law_film(film(1.0, 1e-6))
            .unwrap()
            .prepare()
            .unwrap();
        original.set_analytic_newton(analytic);
        limit.set_analytic_newton(analytic);
        for _ in 0..16 {
            let expected = original
                .step(&[0.0; 2], &CancelGate::new_clock_free())
                .unwrap();
            let actual = limit
                .step(&[0.0; 2], &CancelGate::new_clock_free())
                .unwrap();
            assert_eq!(limit.state(), original.state());
            assert_eq!(actual.dissipated_energy_j, expected.dissipated_energy_j);
        }
    }
}

fn mixed() -> ImpactSystem {
    let native = graph(1e-6);
    let sliding = ShearedFilm::new(
        native.clone(),
        vec![WallVelocityPorts {
            lower: vec![0.0; 2],
            upper: vec![0.0; 2],
        }],
    )
    .unwrap();
    pair()
        .with_squeeze_film(native.clone())
        .unwrap()
        .with_squeeze_film(native)
        .unwrap()
        .with_sliding_film(sliding)
        .unwrap()
}

#[test]
fn g4_domain_and_shared_island_refusals_leave_motion_and_time_unchanged() {
    assert_eq!(mixed().incompressible_film_count(), 4);
    assert!(mixed().with_power_law_film(film(0.5, 1e-6)).is_err());
    assert!(mixed().with_squeeze_film(graph(1e-6)).is_err());
    let sliding = ShearedFilm::new(
        graph(1e-6),
        vec![WallVelocityPorts {
            lower: vec![0.0; 2],
            upper: vec![0.0; 2],
        }],
    )
    .unwrap();
    assert!(mixed().with_sliding_film(sliding).is_err());
    let gas = AmbientGas {
        pressure_pa: 101325.0,
        temperature_k: 300.0,
        specific_gas_constant_j_kg_k: 287.05,
    };
    assert!(
        pair()
            .with_compressible_squeeze_film(graph(1e-6), gas)
            .is_err()
    );
    assert!(
        bare()
            .with_compressible_squeeze_film(graph(1e-6), gas)
            .unwrap()
            .with_power_law_film(film(0.5, 1e-6))
            .is_err()
    );

    // The midpoint remains admissible but this first endpoint closes past
    // the declared cell-gap domain. Original acceptance must reject the tick.
    let minimum_gap = GAP - 0.75 * DT * LAUNCH;
    for analytic in [false, true] {
        let mut narrow = bare()
            .with_power_law_film(film(0.5, minimum_gap))
            .unwrap()
            .prepare()
            .unwrap();
        narrow.set_analytic_newton(analytic);
        let before = narrow.state().to_vec();
        let energy = narrow.stored_energy_j();
        assert!(
            narrow
                .step(&[0.0; 2], &CancelGate::new_clock_free())
                .is_err()
        );
        assert_eq!(narrow.state(), before);
        assert_eq!(narrow.samples(), 0);
        assert_eq!(narrow.stored_energy_j(), energy);
    }
}
