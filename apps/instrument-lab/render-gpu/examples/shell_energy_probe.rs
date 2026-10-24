//! Native time-refinement diagnostic for an isolated, lossless curved annulus.
//! No contact, support, radiation or prescribed damping acts. Compare initial
//! velocity with a one-millisecond force pulse starting from exact rest; energy
//! is measured for twenty milliseconds after the pulse has completely ended.
//! Reports numerical energy drift; this is not an acoustic acceptance test.
use fs_vkplate::{Body, Drive, Edge, Instrument, Nonlinear, PlateSpec};

struct ResultRow {
    rate: f64,
    initial: f64,
    final_energy: f64,
    minimum: f64,
    maximum: f64,
    w: Vec<f64>,
    mass: Vec<f64>,
}

fn run(rate: f64, nonlinear: bool, driven: bool) -> ResultRow {
    let mut body = Body::plate_without_surface_loss(
        "isolated curved annulus",
        &PlateSpec {
            radius: 0.1,
            hole: 0.02,
            grid: 0.01,
            edge: Edge::Free,
            thickness: &|_| 0.001,
            rest_height: &|r| 0.2 * r * r,
            young: 1e11,
            poisson: 0.3,
            density: 8000.0,
            tension: 0.0,
            bending: true,
            nonlinear: if nonlinear {
                Nonlinear::VonKarman
            } else {
                Nonlinear::None
            },
            t60: [(100.0, 1.0), (1000.0, 1.0)], // Ignored by this constructor.
            support: None,
            muffle: None,
        },
    )
    .expect("valid fixed physical shell");
    let velocity: Vec<f64> = body
        .xy
        .iter()
        .map(|p| 0.05 * (-(p[0] * p[0] + p[1] * p[1]) / 0.05_f64.powi(2)).exp())
        .collect();
    let total_mass = body.mass.iter().sum::<f64>();
    let mean = body
        .mass
        .iter()
        .zip(&velocity)
        .map(|(m, v)| m * v)
        .sum::<f64>()
        / total_mass;
    // Zero net translation, the same velocity field at every time resolution.
    if !driven {
        for (previous, velocity) in body.wp.iter_mut().zip(&velocity) {
            *previous = -(velocity - mean) / rate;
        }
    }
    body.sync();
    let mut instrument = Instrument::new(rate);
    if driven {
        // A triangular pulse has area duration/2. Its fixed nodal impulse is
        // mass times the same zero-translation velocity field used above.
        for (node, (mass, velocity)) in body.mass.iter().zip(&velocity).enumerate() {
            instrument.drives.push(Drive {
                body: 0,
                node,
                knots: vec![
                    (0.0, 0.0),
                    (0.0005, 2.0 * mass * (velocity - mean) / 0.001),
                    (0.001, 0.0),
                ],
            });
        }
    }
    instrument.add(body);
    let stability = instrument.stability().0;
    assert!(
        stability.is_finite() && stability <= 0.6,
        "clock not admitted: {stability}"
    );
    if driven {
        for _ in 0..(rate * 0.001).round() as usize {
            instrument.step(&mut []);
        }
    }
    let initial = instrument.bodies[0].discrete_energy(1.0 / rate);
    assert!(initial.is_finite() && initial > 0.0);
    let (mut minimum, mut maximum) = (initial, initial);
    for step in 0..(rate * 0.02).round() as usize {
        instrument.step(&mut []);
        // Compare the same physical 1 ms boundaries at every rate.
        if (step + 1) % (rate / 1000.0).round() as usize == 0 {
            let energy = instrument.bodies[0].discrete_energy(1.0 / rate);
            assert!(energy.is_finite());
            minimum = minimum.min(energy);
            maximum = maximum.max(energy);
        }
    }
    let body = &instrument.bodies[0];
    assert!(body.w.iter().all(|x| x.is_finite()));
    ResultRow {
        rate,
        initial,
        final_energy: body.discrete_energy(1.0 / rate),
        minimum,
        maximum,
        w: body.w.clone(),
        mass: body.mass.clone(),
    }
}

fn main() {
    for driven in [false, true] {
        for nonlinear in [false, true] {
            report(nonlinear, driven);
        }
    }
}

fn report(nonlinear: bool, driven: bool) {
    let rows: Vec<_> = [48_000.0, 96_000.0, 192_000.0, 384_000.0]
        .map(|rate| run(rate, nonlinear, driven))
        .into();
    let finest = rows.last().unwrap();
    let norm = finest
        .mass
        .iter()
        .zip(&finest.w)
        .map(|(m, w)| m * w * w)
        .sum::<f64>();
    assert!(norm > 0.0);
    for row in &rows {
        let relative_state = row
            .mass
            .iter()
            .zip(&row.w)
            .zip(&finest.w)
            .map(|((m, a), b)| m * (a - b).powi(2))
            .sum::<f64>()
            .sqrt()
            / norm.sqrt();
        println!(
            "{{\"driven_from_rest\":{},\"nonlinear\":{},\"rate_hz\":{},\"nodes\":{},\"free_duration_s\":0.02,\"initial_j\":{:.17e},\"final_j\":{:.17e},\"final_relative_drift\":{:.17e},\"minimum_relative_drift\":{:.17e},\"maximum_relative_drift\":{:.17e},\"final_state_relative_to_finest\":{:.17e}}}",
            driven,
            nonlinear,
            row.rate,
            row.w.len(),
            row.initial,
            row.final_energy,
            row.final_energy / row.initial - 1.0,
            row.minimum / row.initial - 1.0,
            row.maximum / row.initial - 1.0,
            relative_state,
        );
    }
}
