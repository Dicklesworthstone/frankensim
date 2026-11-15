//! fs-system's conductance between two heat traces, `elements::ConductanceMatrix`, and the field
//! heat ports it joins: conformance and tangents over a sweep, its heat and its transpose, the
//! one-point case against `Conductance` to the bit, the entropy it generates at every drawn
//! field and at two fields whose nodal sum is negative, its refusals, and two traces of lumps
//! stepped through the kernel against the monolithic backward-Euler solve.
//! Run with `cargo test -p fs-system --test conductance_matrix -- --nocapture`.

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use fs_ports::StableId;
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt};
use fs_system::elements::{Conductance, ConductanceMatrix, InterfacePoint};
use fs_system::{AdmissionError, HeatSink, RouteTemperature, SystemBuilder};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 7,
            kernel_id: 52,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

/// Deterministic numbers in `[lo, hi)` (splitmix64), so a sweep is the same on every run.
struct Numbers(u64);

impl Numbers {
    fn next(&mut self, lo: f64, hi: f64) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        lo + (hi - lo) * ((z >> 11) as f64 / (1_u64 << 53) as f64)
    }
}

/// An interface point from each side's `(point, weight)` pairs.
fn point(a: &[(usize, f64)], b: &[(usize, f64)], conductance_w_k: f64) -> InterfacePoint {
    let side = |pairs: &[(usize, f64)]| {
        let mut slots = [(pairs[0].0, 0.0); 3];
        slots[..pairs.len()].copy_from_slice(pairs);
        slots
    };
    InterfacePoint {
        a: side(a),
        b: side(b),
        conductance_w_k,
    }
}

/// One step of `element` at the temperatures `kelvin` through the conformance kit: its loads and
/// its receipt.
fn stepped(element: &ConductanceMatrix, dt_s: f64, kelvin: &[f64]) -> (Vec<f64>, Receipt) {
    let case = Case {
        state: &(),
        stale: &(),
        dt_s,
        inputs: kelvin,
    };
    with_cx(|cx| {
        check_port_tangent(element, &case, 1.0e-8, 1.0e-9, cx)
            .unwrap_or_else(|failure| panic!("tangent at {kelvin:?}: {failure}"));
        match check_component(element, &case, cx) {
            Ok(Observed::Stepped {
                receipt, outputs, ..
            }) => (outputs, receipt),
            other => panic!("conformance at {kelvin:?}: {other:?}"),
        }
    })
}

/// What the kernel's nodal rule makes of a step: minus the entropy `dt sum Q_i / T_i` that
/// entered the element through its two ports.
fn nodal_generation(element: &ConductanceMatrix, dt_s: f64, kelvin: &[f64], loads: &[f64]) -> f64 {
    let split = element.widths()[0];
    let ports = element.ports();
    -(ports[0].entropy_in_j_per_k(dt_s, &kelvin[..split], &loads[..split])
        + ports[1].entropy_in_j_per_k(dt_s, &kelvin[split..], &loads[split..]))
}

/// A matching contact on the unit square cut into two triangles along its diagonal 0-2, points
/// 0 to 3 counter-clockwise: the degree-two quadrature of each triangle, three interface points
/// with the weights 2/3, 1/6, 1/6, each a third of the triangle's area times `conductance_w_m2_k`.
fn square_contact(conductance_w_m2_k: f64) -> ConductanceMatrix {
    let mut points = Vec::new();
    for face in [[0, 1, 2], [0, 2, 3]] {
        for q in 0..3 {
            let weights = [0, 1, 2].map(|i| (face[i], if i == q { 2.0 / 3.0 } else { 1.0 / 6.0 }));
            points.push(point(&weights, &weights, conductance_w_m2_k * 0.5 / 3.0));
        }
    }
    ConductanceMatrix::new(4, 4, points).unwrap()
}

#[test]
fn the_matrix_conforms_keeps_its_heat_and_is_its_own_transpose() {
    // Three points of a and two of b meet at four interface points.
    let element = ConductanceMatrix::new(
        3,
        2,
        vec![
            point(&[(0, 1.0)], &[(0, 1.0)], 4.0),
            point(&[(1, 1.0)], &[(0, 0.5), (1, 0.5)], 3.0),
            point(&[(2, 0.25), (1, 0.75)], &[(1, 1.0)], 2.0),
            point(&[(0, 0.2), (1, 0.3), (2, 0.5)], &[(0, 0.6), (1, 0.4)], 1.5),
        ],
    )
    .unwrap();
    assert_eq!((element.widths(), element.points().len()), ([3, 2], 4));
    let matrix = element.matrix();
    let mut numbers = Numbers(20_261_115);
    for case in 0..200 {
        // Any temperatures: neither side is the hotter one.
        let kelvin: Vec<f64> = (0..5).map(|_| numbers.next(250.0, 900.0)).collect();
        let dt_s = numbers.next(1.0e-4, 10.0);
        let (loads, receipt) = stepped(&element, dt_s, &kelvin);
        // The loads are K T, entry by entry: either side is a sum of at most a dozen products
        // of the size of K's entries times the temperatures.
        for (row, load) in loads.iter().enumerate() {
            let by_hand: f64 = (0..5).map(|j| matrix[row * 5 + j] * kelvin[j]).sum();
            let size: f64 = (0..5)
                .map(|j| (matrix[row * 5 + j] * kelvin[j]).abs())
                .sum();
            assert!(
                (load - by_hand).abs() <= 32.0 * f64::EPSILON * size,
                "case {case} row {row}: {load} against {by_hand}"
            );
        }
        assert!(receipt.entropy_generated_j_per_k >= 0.0, "case {case}");
        assert_eq!(receipt.stored_change_j.to_bits(), 0.0_f64.to_bits());
    }
    // K is symmetric, and a uniform temperature moves no heat.
    for row in 0..5 {
        for column in 0..row {
            let (upper, lower) = (matrix[column * 5 + row], matrix[row * 5 + column]);
            assert!((upper - lower).abs() <= 8.0 * f64::EPSILON * upper.abs().max(1.0));
        }
        let entries = &matrix[row * 5..row * 5 + 5];
        let size: f64 = entries.iter().map(|entry| entry.abs()).sum();
        assert!(entries.iter().sum::<f64>().abs() <= 32.0 * f64::EPSILON * size);
    }
}

#[test]
fn one_point_on_each_side_is_the_conductance_to_the_bit() {
    let matrix = ConductanceMatrix::new(1, 1, vec![point(&[(0, 1.0)], &[(0, 1.0)], 2.5)]).unwrap();
    let scalar = Conductance::new(2.5).unwrap();
    let mut numbers = Numbers(11);
    for _ in 0..100 {
        let kelvin = [numbers.next(250.0, 900.0), numbers.next(250.0, 900.0)];
        let dt_s = numbers.next(1.0e-3, 5.0);
        let (loads, receipt) = stepped(&matrix, dt_s, &kelvin);
        let case = Case {
            state: &(),
            stale: &(),
            dt_s,
            inputs: &kelvin,
        };
        let Ok(Observed::Stepped {
            receipt: scalar_receipt,
            outputs,
            ..
        }) = with_cx(|cx| check_component(&scalar, &case, cx))
        else {
            panic!("the conductance did not step at {kelvin:?}");
        };
        for (load, expected) in loads.iter().zip(&outputs) {
            assert_eq!(load.to_bits(), expected.to_bits(), "{kelvin:?}");
        }
        assert_eq!(
            receipt.entropy_generated_j_per_k.to_bits(),
            scalar_receipt.entropy_generated_j_per_k.to_bits(),
            "{kelvin:?}"
        );
    }
}

/// Heat crosses at every interface point from the hotter interpolated temperature to the
/// colder, so the element generates entropy at any field. The kernel's nodal sum over the same
/// loads does not: it is printed beside it.
#[test]
fn a_contact_generates_entropy_at_any_field_where_the_nodal_sum_does_not() {
    let contact = square_contact(1250.0);
    let mut numbers = Numbers(20_261_115);
    let (mut smallest, mut negative, mut lowest) = (f64::INFINITY, 0, 0.0_f64);
    for case in 0..2000 {
        // Nothing smooth about it: every point of either side at its own temperature.
        let kelvin: Vec<f64> = (0..8).map(|_| numbers.next(250.0, 1500.0)).collect();
        let (loads, receipt) = stepped(&contact, 0.25, &kelvin);
        let generated = receipt.entropy_generated_j_per_k;
        assert!(
            generated > 0.0,
            "case {case}: {generated} J/K at {kelvin:?}"
        );
        let nodal = nodal_generation(&contact, 0.25, &kelvin, &loads);
        negative += usize::from(nodal < 0.0);
        lowest = lowest.min(nodal / generated);
        smallest = smallest.min(generated);
    }
    eprintln!(
        "2000 drawn fields: least entropy generated {smallest:e} J/K; the nodal sum is negative in {negative} of them, at worst {lowest:e} of the generation"
    );
    // Two fields the nodal sum calls a destruction of entropy, both taken.
    // A line of two points a side under Simpson's rule (ends 1, middle 4): neighbours at 1200 K
    // and 300 K whose jumps oppose, 4 K against -1 K.
    let line = ConductanceMatrix::new(
        2,
        2,
        vec![
            point(&[(0, 1.0)], &[(0, 1.0)], 1.0),
            point(&[(0, 0.5), (1, 0.5)], &[(0, 0.5), (1, 0.5)], 4.0),
            point(&[(1, 1.0)], &[(1, 1.0)], 1.0),
        ],
    )
    .unwrap();
    // One point of a fed in equal parts by two points of b: 314.5 K from 300 K and 330 K.
    let averaged =
        ConductanceMatrix::new(1, 2, vec![point(&[(0, 1.0)], &[(0, 0.5), (1, 0.5)], 1.0)]).unwrap();
    for (name, element, kelvin) in [
        (
            "opposed jumps at 1200 K beside 300 K",
            &line,
            vec![1202.0, 299.5, 1198.0, 300.5],
        ),
        (
            "314.5 K from 300 K and 330 K",
            &averaged,
            vec![314.5, 300.0, 330.0],
        ),
    ] {
        let (loads, receipt) = stepped(element, 1.0, &kelvin);
        let nodal = nodal_generation(element, 1.0, &kelvin, &loads);
        eprintln!(
            "{name}: generated {:e} J/K, nodal sum {nodal:e} J/K",
            receipt.entropy_generated_j_per_k
        );
        assert!(
            receipt.entropy_generated_j_per_k > 0.0 && nodal < 0.0,
            "{name}"
        );
    }
}

#[test]
fn the_matrix_refuses_what_is_not_an_interface_between_two_traces() {
    let refused =
        |a: usize, b: usize, points: Vec<InterfacePoint>, what: &str| match ConductanceMatrix::new(
            a, b, points,
        ) {
            Err(error) => assert!(format!("{error:?}").contains(what), "{error:?} for {what}"),
            Ok(_) => panic!("{what}: accepted"),
        };
    let pair = |conductance_w_k: f64| point(&[(0, 1.0)], &[(0, 1.0)], conductance_w_k);
    refused(0, 1, vec![pair(1.0)], "port width");
    refused(1, 1, Vec::new(), "interface points");
    refused(
        1,
        1,
        vec![point(&[(1, 1.0)], &[(0, 1.0)], 1.0)],
        "weight's point",
    );
    refused(
        1,
        2,
        vec![point(&[(0, 1.0)], &[(0, 1.5), (1, -0.5)], 1.0)],
        "not negative",
    );
    refused(
        1,
        1,
        vec![point(&[(0, f64::NAN)], &[(0, 1.0)], 1.0)],
        "not negative",
    );
    // Weights that do not sum to one would make heat of a uniform temperature.
    refused(
        1,
        2,
        vec![point(&[(0, 1.0)], &[(0, 0.5), (1, 0.4)], 1.0)],
        "weights' sum",
    );
    for conductance_w_k in [0.0, -1.0, f64::INFINITY] {
        refused(1, 1, vec![pair(conductance_w_k)], "interface conductance");
    }
}

/// Lumps of heat capacity behind one field heat port: point `i` receives the heat load into
/// lump `i` [W] and returns that lump's END temperature, backward Euler's reading of it.
#[derive(Debug, Clone)]
struct Lumps {
    capacity_j_k: Vec<f64>,
    ports: [PortSpec; 1],
}

impl Lumps {
    fn new(capacity_j_k: &[f64]) -> Self {
        let width = NonZeroUsize::new(capacity_j_k.len()).expect("a lump");
        let id = StableId::new("heat").expect("id");
        Self {
            capacity_j_k: capacity_j_k.to_vec(),
            ports: [PortSpec::heat_field(id, PortRole::FlowIn, width).expect("a field heat port")],
        }
    }
}

impl Component for Lumps {
    type State = Vec<f64>;
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &Vec<f64>) -> Result<f64, Cancelled> {
        Ok(self
            .capacity_j_k
            .iter()
            .zip(state)
            .map(|(c, t)| c * t)
            .sum())
    }

    fn propose(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Vec<f64>,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        new.clone_from(old);
        let mut change = 0.0;
        for (i, capacity) in self.capacity_j_k.iter().enumerate() {
            new[i] = old[i] + dt_s * inputs[i] / capacity;
            outputs[i] = new[i];
            change += capacity * (new[i] - old[i]);
        }
        let scale = self.stored_energy_j(old)?.abs() + self.stored_energy_j(new)?.abs();
        Ok(Proposal::Step(Receipt::lossless(
            change,
            64.0 * f64::EPSILON * scale,
        )))
    }
}

impl PortTangent for Lumps {
    fn jvp(
        &self,
        _old: &Vec<f64>,
        dt_s: f64,
        _inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        for (i, capacity) in self.capacity_j_k.iter().enumerate() {
            d_outputs[i] = dt_s * d_inputs[i] / capacity;
        }
        Ok(())
    }

    fn vjp(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        // Diagonal, so its own transpose.
        self.jvp(old, dt_s, inputs, bar_outputs, bar_inputs, cx)
    }
}

/// `a x = b` by elimination with partial pivoting, for the monolithic step below.
fn solve(mut a: Vec<f64>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for column in 0..n {
        let pivot = (column..n)
            .max_by(|&i, &j| a[i * n + column].abs().total_cmp(&a[j * n + column].abs()))
            .expect("a row");
        for k in 0..n {
            a.swap(column * n + k, pivot * n + k);
        }
        b.swap(column, pivot);
        for row in (column + 1)..n {
            let factor = a[row * n + column] / a[column * n + column];
            for k in column..n {
                a[row * n + k] -= factor * a[column * n + k];
            }
            b[row] -= factor * b[column];
        }
    }
    for row in (0..n).rev() {
        for k in (row + 1)..n {
            b[row] -= a[row * n + k] * b[k];
        }
        b[row] /= a[row * n + row];
    }
    b
}

#[test]
fn two_traces_of_lumps_through_a_matrix_step_as_the_monolithic_backward_euler() {
    let capacity = [40.0, 55.0, 70.0, 25.0];
    // Each pair of lumps through its own point, and both pairs through one between them:
    // the block of a is [[0.8, 0.3], [0.3, 0.5]] W/K.
    let wall = ConductanceMatrix::new(
        2,
        2,
        vec![
            point(&[(0, 1.0)], &[(0, 1.0)], 0.5),
            point(&[(1, 1.0)], &[(1, 1.0)], 0.2),
            point(&[(0, 0.5), (1, 0.5)], &[(0, 0.5), (1, 0.5)], 1.2),
        ],
    )
    .unwrap();
    let conductance = wall.matrix();
    let mut kelvin = vec![420.0, 380.0, 300.0, 330.0];
    let dt_s = 2.0;
    let mut builder = SystemBuilder::new();
    let a = builder
        .add_with_tangent("a", Lumps::new(&capacity[..2]), kelvin[..2].to_vec())
        .unwrap();
    let b = builder
        .add_with_tangent("b", Lumps::new(&capacity[2..]), kelvin[2..].to_vec())
        .unwrap();
    let wall = builder.add_with_tangent("wall", wall, ()).unwrap();
    let (wall_a, wall_b) = (
        builder.port(wall, "a").unwrap(),
        builder.port(wall, "b").unwrap(),
    );
    // A heat route moves one heat rate: a trace is refused by name.
    let sink = HeatSink {
        name: StableId::new("air").expect("id"),
        temperature_k: 293.15,
    };
    match builder.heat_route(wall_a, &[], Some(sink), RouteTemperature::Fixed(293.15)) {
        Err(AdmissionError::HeatRouteWidth { port, width: 2 }) => assert_eq!(port, "wall.a"),
        other => panic!("a route from a trace: {other:?}"),
    }
    builder
        .bond(builder.port(a, "heat").unwrap(), wall_a)
        .unwrap();
    builder
        .bond(wall_b, builder.port(b, "heat").unwrap())
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(wall_a, &kelvin[..2]));
    assert!(system.guess(wall_b, &kelvin[2..]));
    with_cx(|cx| {
        for step in 0..50 {
            system.step(dt_s, cx).unwrap();
            // (C / dt + K) T' = (C / dt) T, all four lumps at once.
            let mut matrix = conductance.clone();
            let mut right = vec![0.0; 4];
            for i in 0..4 {
                matrix[i * 4 + i] += capacity[i] / dt_s;
                right[i] = capacity[i] / dt_s * kelvin[i];
            }
            kelvin = solve(matrix, right);
            let state_a = system.state::<Vec<f64>>(a).unwrap();
            let state_b = system.state::<Vec<f64>>(b).unwrap();
            for (kernel, monolithic) in state_a.iter().chain(state_b).zip(&kelvin) {
                assert!(
                    (kernel - monolithic).abs() <= 1.0e-11 * monolithic,
                    "step {step}: {kernel} against {monolithic} K"
                );
            }
        }
    });
    let ledger = system.ledger();
    let stored: f64 = capacity.iter().sum::<f64>() * 400.0;
    assert!(ledger.energy_residual_j().abs() <= 1.0e-12 * stored);
    let generated = ledger.nodes[2].entropy_generated_j_per_k.value();
    assert!(generated > 0.0 && ledger.entropy_generated_j_per_k() >= generated * (1.0 - 1.0e-9));
}
