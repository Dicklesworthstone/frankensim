//! fs-system's conductance between two heat traces, `elements::ConductanceMatrix`, and the field
//! heat ports it joins: conformance and tangents over a sweep, its heat and its transpose, the
//! one-point case against `Conductance` to the bit, the entropy it generates at every drawn
//! field and at two fields whose nodal sum is negative, its refusals, and two traces of lumps
//! stepped through the kernel against the monolithic backward-Euler solve.
//! Run with `cargo test -p fs-system --test conductance_matrix -- --nocapture`.

use core::num::NonZeroUsize;
use std::sync::Arc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use fs_ports::StableId;
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt};
use fs_system::elements::{
    Conductance, ConductanceLaw, ConductanceMatrix, ElementError, GuardRow, InterfacePoint,
};
use fs_system::{AdmissionError, HeatSink, NodeId, RouteTemperature, System, SystemBuilder};

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

/// A law authored for these tests, not a material's: the factor
/// `1 + a (T_a / 300)^2 + b T_b / 300` while both temperatures lie between `floor_k` and
/// `ceiling_k`, and no answer outside.
#[derive(Debug, Clone, Copy)]
struct Authored {
    a: f64,
    b: f64,
    floor_k: f64,
    ceiling_k: f64,
}

impl ConductanceLaw for Authored {
    fn factor(&self, kelvin_a: f64, kelvin_b: f64) -> Result<(f64, [f64; 2]), String> {
        let inside = |kelvin: f64| (self.floor_k..=self.ceiling_k).contains(&kelvin);
        if !(inside(kelvin_a) && inside(kelvin_b)) {
            return Err(format!(
                "{kelvin_a} K and {kelvin_b} K are not both between {} K and {} K",
                self.floor_k, self.ceiling_k
            ));
        }
        let (x, y) = (kelvin_a / 300.0, kelvin_b / 300.0);
        Ok((
            1.0 + self.a * x * x + self.b * y,
            [2.0 * self.a * x / 300.0, self.b / 300.0],
        ))
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend([self.a, self.b, self.floor_k, self.ceiling_k].map(f64::to_bits));
    }
}

/// The law of the sweeps below: its factor lies between 0.6 and 4.6 from 250 K to 900 K.
const SWEPT: Authored = Authored {
    a: 0.4,
    b: -0.2,
    floor_k: 200.0,
    ceiling_k: 1000.0,
};

/// The element of the first test: three points of a and two of b at four interface points.
fn four_points() -> ConductanceMatrix {
    ConductanceMatrix::new(
        3,
        2,
        vec![
            point(&[(0, 1.0)], &[(0, 1.0)], 4.0),
            point(&[(1, 1.0)], &[(0, 0.5), (1, 0.5)], 3.0),
            point(&[(2, 0.25), (1, 0.75)], &[(1, 1.0)], 2.0),
            point(&[(0, 0.2), (1, 0.3), (2, 0.5)], &[(0, 0.6), (1, 0.4)], 1.5),
        ],
    )
    .unwrap()
}

#[test]
fn a_law_driven_matrix_conforms_and_its_tangent_and_adjoint_are_transposes() {
    let constant = four_points();
    let element = four_points().following(Arc::new(SWEPT));
    assert_eq!(
        (constant.class(), element.class()),
        (Class::ExactLinear, Class::Smooth)
    );
    let mut numbers = Numbers(20_261_116);
    let mut worst = 0.0_f64;
    for case in 0..200 {
        let kelvin: Vec<f64> = (0..5).map(|_| numbers.next(250.0, 900.0)).collect();
        let dt_s = numbers.next(1.0e-4, 10.0);
        // The kit: the receipt's balance, and the tangent against differences.
        let (loads, receipt) = stepped(&element, dt_s, &kelvin);
        // The loads by hand: every point's own conductance times the law's factor at the two
        // temperatures interpolated there, spread by the same weights.
        let mut by_hand = [0.0_f64; 5];
        let mut size = [0.0_f64; 5];
        for point in element.points() {
            let at = |side: &[(usize, f64); 3], offset: usize| -> f64 {
                side.iter().map(|&(i, w)| w * kelvin[offset + i]).sum()
            };
            let (t_a, t_b) = (at(&point.a, 0), at(&point.b, 3));
            let rate = point.conductance_w_k * SWEPT.factor(t_a, t_b).unwrap().0 * (t_a - t_b);
            for &(i, w) in &point.a {
                by_hand[i] += w * rate;
                size[i] += (w * rate).abs();
            }
            for &(i, w) in &point.b {
                by_hand[3 + i] -= w * rate;
                size[3 + i] += (w * rate).abs();
            }
        }
        for row in 0..5 {
            assert!(
                (loads[row] - by_hand[row]).abs() <= 32.0 * f64::EPSILON * size[row],
                "case {case} row {row}: {} against {}",
                loads[row],
                by_hand[row]
            );
        }
        assert!(receipt.entropy_generated_j_per_k >= 0.0, "case {case}");
        // <J d, bar> = <d, J^T bar> for drawn d and bar: the adjoint is the tangent's transpose.
        let d: Vec<f64> = (0..5).map(|_| numbers.next(-1.0, 1.0)).collect();
        let bar: Vec<f64> = (0..5).map(|_| numbers.next(-1.0, 1.0)).collect();
        let (mut forward, mut backward) = (vec![0.0; 5], vec![0.0; 5]);
        with_cx(|cx| {
            element
                .jvp(&(), dt_s, &kelvin, &d, &mut forward, cx)
                .unwrap();
            element
                .vjp(&(), dt_s, &kelvin, &bar, &mut backward, cx)
                .unwrap();
        });
        let dot = |x: &[f64], y: &[f64]| -> f64 { x.iter().zip(y).map(|(x, y)| x * y).sum() };
        let scale: f64 = forward
            .iter()
            .zip(&bar)
            .map(|(x, y)| (x * y).abs())
            .sum::<f64>()
            + backward
                .iter()
                .zip(&d)
                .map(|(x, y)| (x * y).abs())
                .sum::<f64>();
        let gap = (dot(&forward, &bar) - dot(&d, &backward)).abs() / scale;
        worst = worst.max(gap);
        assert!(
            gap <= 64.0 * f64::EPSILON,
            "case {case}: transpose gap {gap:e}"
        );
    }
    eprintln!("largest relative gap between <J d, bar> and <d, J^T bar> over 200 cases: {worst:e}");
}

#[test]
fn a_law_of_one_is_the_constant_matrix_to_the_bit_and_a_law_is_in_the_identity() {
    /// The factor one with no slope.
    #[derive(Debug)]
    struct One;
    impl ConductanceLaw for One {
        fn factor(&self, _: f64, _: f64) -> Result<(f64, [f64; 2]), String> {
            Ok((1.0, [0.0, 0.0]))
        }
        fn identity_words(&self, words: &mut Vec<u64>) {
            words.push(1);
        }
    }
    let constant = four_points();
    let unit = four_points().following(Arc::new(One));
    let mut numbers = Numbers(20_261_117);
    for case in 0..100 {
        let kelvin: Vec<f64> = (0..5).map(|_| numbers.next(250.0, 900.0)).collect();
        let dt_s = numbers.next(1.0e-4, 10.0);
        let (plain, plain_receipt) = stepped(&constant, dt_s, &kelvin);
        let (scaled, scaled_receipt) = stepped(&unit, dt_s, &kelvin);
        assert_eq!(bits(&plain), bits(&scaled), "case {case}");
        assert_eq!(
            bits(&[
                plain_receipt.entropy_generated_j_per_k,
                plain_receipt.energy_tolerance_j
            ]),
            bits(&[
                scaled_receipt.entropy_generated_j_per_k,
                scaled_receipt.energy_tolerance_j
            ]),
            "case {case}"
        );
    }
    // A constant element's identity is its interface points alone; a law's words follow them.
    let words = |element: &ConductanceMatrix| {
        let mut words = Vec::new();
        element.identity_words(&mut words);
        words
    };
    let (plain, one, swept) = (
        words(&constant),
        words(&unit),
        words(&four_points().following(Arc::new(SWEPT))),
    );
    // After the law's words comes the number of guard rows, here none.
    assert_eq!(plain.len(), 1 + 4 * 13);
    assert_eq!(
        (&one[..plain.len()], &one[plain.len()..]),
        (&plain[..], &[1, 0][..])
    );
    assert_eq!(&swept[..plain.len()], &plain[..]);
    assert_eq!(swept.len(), plain.len() + 5);
    assert_ne!(one, swept);
    // A guard row's twelve words follow that number, and two rows differ.
    let guarded = |row: GuardRow| {
        let element = four_points().following(Arc::new(SWEPT));
        words(&element.guarded_at(vec![row]).unwrap())
    };
    let (first, second) = (guarded(facing(0, 1)), guarded(facing(2, 1)));
    let law_end = swept.len() - 1;
    assert_eq!(
        (&first[..law_end], first[law_end], first.len()),
        (&swept[..law_end], 1, swept.len() + 12)
    );
    assert_ne!(first, second);
}

#[test]
fn a_law_with_no_answer_refuses_the_step_and_nothing_is_clamped() {
    let narrow = Authored {
        ceiling_k: 400.0,
        ..SWEPT
    };
    let element = four_points().following(Arc::new(narrow));
    let propose = |element: &ConductanceMatrix, kelvin: &[f64]| {
        let mut loads = vec![0.0; 5];
        with_cx(|cx| element.propose(&(), 1.0, kelvin, &mut (), &mut loads, cx))
    };
    // The law is asked at the interface points, at the temperatures interpolated there. A trace
    // point above the ceiling that reaches every interface point through a mean below it is
    // answered; the first interface point reads trace point 0 alone, and there the law has none.
    assert!(propose(&element, &[300.0, 320.0, 340.0, 360.0, 380.0]).is_ok());
    assert!(propose(&element, &[300.0, 320.0, 450.0, 360.0, 380.0]).is_ok());
    match propose(&element, &[450.0, 320.0, 340.0, 360.0, 380.0]) {
        Err(ElementError::Law { reason }) => {
            assert!(
                reason.contains("are not both between 200 K and 400 K"),
                "{reason}"
            );
        }
        other => panic!("a state outside the law: {other:?}"),
    }
    // A law that answers with a factor that is not positive is refused by the element.
    let negative = Authored { a: -4.0, ..SWEPT };
    match propose(
        &four_points().following(Arc::new(negative)),
        &[600.0, 600.0, 600.0, 300.0, 300.0],
    ) {
        Err(ElementError::Law { reason }) => {
            assert!(reason.contains("is not positive and finite"), "{reason}");
        }
        other => panic!("a factor below zero: {other:?}"),
    }
    // Through the kernel: the step is refused and the lumps keep their temperatures.
    let start = [450.0, 380.0, 300.0, 330.0];
    let wall = ConductanceMatrix::new(
        2,
        2,
        vec![
            point(&[(0, 1.0)], &[(0, 1.0)], 0.5),
            point(&[(1, 1.0)], &[(1, 1.0)], 0.2),
        ],
    )
    .unwrap()
    .following(Arc::new(narrow));
    let (mut system, lumps) = lumps_through(wall, &start);
    let refusal = with_cx(|cx| system.step(2.0, cx)).expect_err("450 K is above the law");
    let said = format!("{refusal:?}");
    eprintln!("the kernel's refusal: {said}");
    assert!(
        said.contains("are not both between 200 K and 400 K"),
        "{said}"
    );
    assert_eq!(bits(&lump_kelvin(&system, lumps)), bits(&start));
}

/// Two traces of two lumps each, at `start` [K], bonded through `wall`: the sealed system and
/// its two lump nodes.
fn lumps_through(wall: ConductanceMatrix, start: &[f64; 4]) -> (System, [NodeId; 2]) {
    let mut builder = SystemBuilder::new();
    let a = builder
        .add_with_tangent("a", Lumps::new(&[40.0, 55.0]), start[..2].to_vec())
        .unwrap();
    let b = builder
        .add_with_tangent("b", Lumps::new(&[70.0, 25.0]), start[2..].to_vec())
        .unwrap();
    let wall = builder.add_with_tangent("wall", wall, ()).unwrap();
    let (wall_a, wall_b) = (
        builder.port(wall, "a").unwrap(),
        builder.port(wall, "b").unwrap(),
    );
    builder
        .bond(builder.port(a, "heat").unwrap(), wall_a)
        .unwrap();
    builder
        .bond(wall_b, builder.port(b, "heat").unwrap())
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(wall_a, &start[..2]));
    assert!(system.guess(wall_b, &start[2..]));
    (system, [a, b])
}

/// The four lumps' temperatures [K], trace a's first.
fn lump_kelvin(system: &System, lumps: [NodeId; 2]) -> Vec<f64> {
    lumps
        .iter()
        .flat_map(|&node| system.state::<Vec<f64>>(node).unwrap().clone())
        .collect()
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

/// A guard row with the weight one on point `a` of trace a and point `b` of trace b.
fn facing(a: usize, b: usize) -> GuardRow {
    GuardRow {
        a: [(a, 1.0), (a, 0.0), (a, 0.0)],
        b: [(b, 1.0), (b, 0.0), (b, 0.0)],
    }
}

#[test]
fn a_guard_row_refuses_a_trace_point_that_the_interface_points_average_away() {
    let narrow = Authored {
        ceiling_k: 400.0,
        ..SWEPT
    };
    // The one interface point reads the mean of each trace's two lumps. Lump 0 of a is above
    // the law's ceiling and the mean of a, 390 K, is below it.
    let start = [450.0, 330.0, 300.0, 330.0];
    let wall = || {
        let mean = [(0, 0.5), (1, 0.5)];
        ConductanceMatrix::new(2, 2, vec![point(&mean, &mean, 0.5)])
            .unwrap()
            .following(Arc::new(narrow))
    };
    // Without guard rows the law is asked at the mean alone, and the step is taken.
    let (mut unguarded, lumps) = lumps_through(wall(), &start);
    with_cx(|cx| unguarded.step(2.0, cx)).expect("the mean is inside the law");
    assert_ne!(bits(&lump_kelvin(&unguarded, lumps)), bits(&start));
    // Guarded at the two pairs of lumps that face each other, the same state is refused
    // through the kernel with the law's own reason, and the lumps keep their temperatures.
    let guarded = wall().guarded_at(vec![facing(0, 0), facing(1, 1)]).unwrap();
    let (mut system, lumps) = lumps_through(guarded, &start);
    let refusal = with_cx(|cx| system.step(2.0, cx)).expect_err("lump 0 of a is above the law");
    let said = format!("{refusal:?}");
    eprintln!("the kernel's refusal at a guard row: {said}");
    assert!(
        said.contains("450 K and 300 K are not both between 200 K and 400 K"),
        "{said}"
    );
    assert_eq!(bits(&lump_kelvin(&system, lumps)), bits(&start));
    // Inside the law at every lump, the guarded wall steps as the unguarded one, to the bit:
    // a guard row moves no heat.
    let inside = [390.0, 330.0, 300.0, 330.0];
    let (mut plain, plain_lumps) = lumps_through(wall(), &inside);
    let guarded = wall().guarded_at(vec![facing(0, 0), facing(1, 1)]).unwrap();
    let (mut watched, watched_lumps) = lumps_through(guarded, &inside);
    with_cx(|cx| {
        for _ in 0..10 {
            plain.step(2.0, cx).unwrap();
            watched.step(2.0, cx).unwrap();
        }
    });
    assert_eq!(
        bits(&lump_kelvin(&plain, plain_lumps)),
        bits(&lump_kelvin(&watched, watched_lumps))
    );
    // Rows guard a law: without one they are refused, and a row's weights are an interface
    // point's.
    let lawless = ConductanceMatrix::new(2, 2, vec![point(&[(0, 1.0)], &[(0, 1.0)], 0.5)]).unwrap();
    for (case, built, name) in [
        (
            "no law",
            lawless.guarded_at(vec![facing(0, 0)]),
            "guard rows",
        ),
        (
            "outside",
            wall().guarded_at(vec![facing(2, 0)]),
            "weight's point",
        ),
        (
            "sum",
            wall().guarded_at(vec![GuardRow {
                a: [(0, 0.5), (1, 0.4), (1, 0.0)],
                ..facing(0, 0)
            }]),
            "weights' sum",
        ),
    ] {
        match built {
            Err(ElementError::Parameter { name: refused, .. }) => {
                assert!(refused.contains(name), "{case}: {refused}");
            }
            other => panic!("{case}: {other:?}"),
        }
    }
}

/// The loads [W] into the wall of the kernel tests under `law`, written out point by point:
/// each pair of lumps through its own point, and the two pairs' means through one between them.
fn wall_loads(kelvin: &[f64], law: &Authored) -> [f64; 4] {
    let rate = |conductance_w_k: f64, t_a: f64, t_b: f64| {
        conductance_w_k * law.factor(t_a, t_b).expect("inside the law").0 * (t_a - t_b)
    };
    let own = [
        rate(0.5, kelvin[0], kelvin[2]),
        rate(0.2, kelvin[1], kelvin[3]),
    ];
    let shared = rate(
        1.2,
        0.5 * (kelvin[0] + kelvin[1]),
        0.5 * (kelvin[2] + kelvin[3]),
    );
    [
        own[0] + 0.5 * shared,
        own[1] + 0.5 * shared,
        -own[0] - 0.5 * shared,
        -own[1] - 0.5 * shared,
    ]
}

/// One backward-Euler step of the four lumps through that wall, all at once:
/// `C (T' - T) / dt + L(T') = 0` by Newton, its Jacobian by central differences of the loads
/// written out in `wall_loads`. Nothing of the element's own tangent is used.
fn backward_euler_newton(capacity: &[f64; 4], old: &[f64], dt_s: f64, law: &Authored) -> Vec<f64> {
    let residual = |new: &[f64]| -> Vec<f64> {
        let loads = wall_loads(new, law);
        (0..4)
            .map(|i| capacity[i] * (new[i] - old[i]) / dt_s + loads[i])
            .collect()
    };
    let mut new = old.to_vec();
    for _ in 0..50 {
        let now = residual(&new);
        let mut jacobian = vec![0.0; 16];
        for column in 0..4 {
            let mut shifted = new.clone();
            shifted[column] += 1.0e-4;
            let up = residual(&shifted);
            shifted[column] -= 2.0e-4;
            let down = residual(&shifted);
            for row in 0..4 {
                jacobian[row * 4 + column] = (up[row] - down[row]) / 2.0e-4;
            }
        }
        let update = solve(jacobian, now);
        for (value, change) in new.iter_mut().zip(&update) {
            *value -= change;
        }
        if update.iter().all(|change| change.abs() <= 1.0e-12) {
            return new;
        }
    }
    panic!("the monolithic Newton solve did not converge");
}

#[test]
fn two_traces_of_lumps_through_a_law_driven_matrix_step_as_the_backward_euler_newton_solve() {
    // Declared before the first run: the kernel's lumps within 1e-9 of the monolithic ones,
    // relative, at every step.
    const CLAIMED: f64 = 1.0e-9;
    let capacity = [40.0, 55.0, 70.0, 25.0];
    let points = vec![
        point(&[(0, 1.0)], &[(0, 1.0)], 0.5),
        point(&[(1, 1.0)], &[(1, 1.0)], 0.2),
        point(&[(0, 0.5), (1, 0.5)], &[(0, 0.5), (1, 0.5)], 1.2),
    ];
    let wall = ConductanceMatrix::new(2, 2, points)
        .unwrap()
        .following(Arc::new(SWEPT));
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
    builder
        .bond(builder.port(a, "heat").unwrap(), wall_a)
        .unwrap();
    builder
        .bond(wall_b, builder.port(b, "heat").unwrap())
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(wall_a, &kelvin[..2]));
    assert!(system.guess(wall_b, &kelvin[2..]));
    // The same four lumps with the wall's constant conductances, for what the law changes.
    let mut unscaled = kelvin.clone();
    let constant = Authored {
        a: 0.0,
        b: 0.0,
        ..SWEPT
    };
    let monolithic = |old: &[f64], law: &Authored| backward_euler_newton(&capacity, old, dt_s, law);
    let mut worst = 0.0_f64;
    with_cx(|cx| {
        for step in 0..50 {
            system.step(dt_s, cx).unwrap();
            kelvin = monolithic(&kelvin, &SWEPT);
            unscaled = monolithic(&unscaled, &constant);
            let state_a = system.state::<Vec<f64>>(a).unwrap();
            let state_b = system.state::<Vec<f64>>(b).unwrap();
            for (kernel, reference) in state_a.iter().chain(state_b).zip(&kelvin) {
                let gap = (kernel - reference).abs() / reference;
                worst = worst.max(gap);
                assert!(
                    gap <= CLAIMED,
                    "step {step}: {kernel} against {reference} K"
                );
            }
        }
    });
    // The law is felt: the same lumps through the constant wall end elsewhere.
    let moved = kelvin
        .iter()
        .zip(&unscaled)
        .map(|(with, without)| (with - without).abs())
        .fold(0.0, f64::max);
    eprintln!(
        "law-driven wall against the monolithic Newton solve over 50 steps: largest relative difference {worst:e}; the law moves a lump by up to {moved} K from the constant wall's"
    );
    assert!(moved > 0.1, "the law changed the lumps by {moved} K only");
    let ledger = system.ledger();
    let stored: f64 = capacity.iter().sum::<f64>() * 400.0;
    assert!(ledger.energy_residual_j().abs() <= 1.0e-12 * stored);
    let generated = ledger.nodes[2].entropy_generated_j_per_k.value();
    assert!(generated > 0.0 && ledger.entropy_generated_j_per_k() >= generated * (1.0 - 1.0e-9));
}
