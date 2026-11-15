//! fs-system's conductance between two heat traces, `elements::ConductanceMatrix`, and the field
//! heat ports it joins: conformance and tangents over a sweep, its heat and its transpose, the
//! one-point case against `Conductance`, the entropy it generates for smooth fields and the
//! refusal where it would destroy some, its refusals, and two traces of lumps stepped through
//! the kernel against the monolithic backward-Euler solve.
//! Run with `cargo test -p fs-system --test conductance_matrix`.

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use fs_ports::StableId;
use fs_ports::component::conformance::{
    Case, Failure, Observed, check_component, check_port_tangent,
};
use fs_ports::component::{
    Class, Component, ContractViolation, PortRole, PortSpec, PortTangent, Proposal, Receipt,
};
use fs_system::elements::{Conductance, ConductanceMatrix, ElementError};
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

/// A matching contact on the unit square cut into two triangles along its diagonal 0-2, points
/// 0 to 3 counter-clockwise: `D = [I, -I]` and `W` the consistent P1 face mass
/// `A (1 + delta) / 12` of the two faces times `conductance_w_m2_k`.
fn square_contact(conductance_w_m2_k: f64) -> ConductanceMatrix {
    let mut weights = vec![0.0; 16];
    for face in [[0, 1, 2], [0, 2, 3]] {
        for &row in &face {
            for &column in &face {
                let mass = if row == column { 0.5 / 6.0 } else { 0.5 / 12.0 };
                weights[row * 4 + column] += conductance_w_m2_k * mass;
            }
        }
    }
    let mut jump = vec![0.0; 4 * 8];
    for point in 0..4 {
        jump[point * 8 + point] = 1.0;
        jump[point * 8 + 4 + point] = -1.0;
    }
    ConductanceMatrix::new(4, 4, jump, weights).unwrap()
}

#[test]
fn the_matrix_conforms_keeps_its_heat_and_is_its_own_transpose() {
    // Three jumps over three points of a and two of b; every row sums to zero.
    let jump = vec![
        1.0, 0.0, 0.0, -1.0, 0.0, //
        0.0, 1.0, 0.0, -0.5, -0.5, //
        0.0, 0.0, 1.0, 0.0, -1.0,
    ];
    let weights = vec![4.0, 1.0, 0.0, 1.0, 3.0, 0.5, 0.0, 0.5, 2.0];
    let element = ConductanceMatrix::new(3, 2, jump, weights).unwrap();
    assert_eq!((element.widths(), element.rows()), ([3, 2], 3));
    let matrix = element.matrix();
    let mut numbers = Numbers(20_261_115);
    for case in 0..200 {
        // Side a hotter than side b at every point: each jump is positive, each jump of 1 / T
        // negative, and no entry of W is, so no term of the entropy sum can be negative.
        let mut kelvin: Vec<f64> = (0..3).map(|_| numbers.next(500.0, 900.0)).collect();
        kelvin.extend((0..2).map(|_| numbers.next(250.0, 400.0)));
        let dt_s = numbers.next(1.0e-4, 10.0);
        let (loads, receipt) = stepped(&element, dt_s, &kelvin);
        // The loads are K T, entry by entry, and their sum is no heat at all.
        let size: f64 = loads.iter().map(|load| load.abs()).sum();
        for (row, load) in loads.iter().enumerate() {
            let by_hand: f64 = (0..5).map(|j| matrix[row * 5 + j] * kelvin[j]).sum();
            assert!(
                (load - by_hand).abs() <= 64.0 * f64::EPSILON * 5.0 * 900.0 * 10.0,
                "case {case} row {row}: {load} against {by_hand}"
            );
        }
        assert!(
            loads.iter().sum::<f64>().abs() <= 32.0 * f64::EPSILON * size,
            "case {case}: net {} W of {size} W",
            loads.iter().sum::<f64>()
        );
        assert!(receipt.entropy_generated_j_per_k >= 0.0, "case {case}");
        assert_eq!(receipt.stored_change_j.to_bits(), 0.0_f64.to_bits());
    }
    // K is symmetric, and a uniform temperature moves no heat.
    for row in 0..5 {
        for column in 0..row {
            let (upper, lower) = (matrix[column * 5 + row], matrix[row * 5 + column]);
            assert!((upper - lower).abs() <= 8.0 * f64::EPSILON * upper.abs().max(1.0));
        }
        assert!(matrix[row * 5..row * 5 + 5].iter().sum::<f64>().abs() <= 32.0 * f64::EPSILON);
    }
}

#[test]
fn one_point_on_each_side_is_the_conductance() {
    let matrix = ConductanceMatrix::new(1, 1, vec![1.0, -1.0], vec![2.5]).unwrap();
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
        let heat = 2.5 * (kelvin[0] - kelvin[1]);
        for (load, expected) in loads.iter().zip(&outputs) {
            assert!((load - expected).abs() <= 4.0 * f64::EPSILON * heat.abs());
        }
        // The matrix differences 1 / T_a and 1 / T_b, each rounded at half an eps of itself, so
        // its sum errs by eps T / |T_a - T_b| of the conductance's closed form, and a few eps.
        let generated = scalar_receipt.entropy_generated_j_per_k;
        let hottest = kelvin[0].max(kelvin[1]);
        let allowed = 8.0 * f64::EPSILON * (1.0 + hottest / (kelvin[0] - kelvin[1]).abs());
        assert!(
            (receipt.entropy_generated_j_per_k - generated).abs() <= allowed * generated,
            "{} against {generated} J/K",
            receipt.entropy_generated_j_per_k
        );
    }
}

/// Smooth, here: on each side a temperature linear over the square, its corner values between
/// 300 K and 400 K, plus up to half a kelvin at every point. The two sides are drawn apart, so
/// the jump may change sign across the patch.
#[test]
fn a_contact_generates_entropy_for_smooth_fields_and_is_refused_where_it_would_destroy_it() {
    let contact = square_contact(1250.0);
    let mut numbers = Numbers(20_261_115);
    let mut smallest = f64::INFINITY;
    for case in 0..2000 {
        let mut kelvin = Vec::with_capacity(8);
        for _side in 0..2 {
            let (origin, along_x, along_y) = (
                numbers.next(300.0, 350.0),
                numbers.next(0.0, 25.0),
                numbers.next(0.0, 25.0),
            );
            for (x, y) in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] {
                kelvin.push(origin + along_x * x + along_y * y + numbers.next(-0.5, 0.5));
            }
        }
        let (loads, receipt) = stepped(&contact, 0.25, &kelvin);
        // Against the entropy the same loads would carry at 350 K.
        let scale = 0.25 * loads.iter().map(|load| load.abs()).sum::<f64>() / 350.0;
        assert!(receipt.entropy_generated_j_per_k >= 0.0, "case {case}");
        smallest = smallest.min(receipt.entropy_generated_j_per_k / scale);
    }
    eprintln!(
        "smallest entropy generated over 2000 smooth fields: {smallest:e} of dt sum |Q| / 350 K"
    );
    // Not smooth: neighbours at 1200 K and 300 K whose jumps oppose, 4 K against -1 K. The
    // consistent mass joins them, the nodal rule divides each load by its own temperature, and
    // with 1 / (T_a T_b) sixteen times apart between the two points the sum is negative.
    let line = ConductanceMatrix::new(
        2,
        2,
        vec![1.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, -1.0],
        vec![2.0, 1.0, 1.0, 2.0],
    )
    .unwrap();
    let kelvin = [1202.0, 299.5, 1198.0, 300.5];
    let case = Case {
        state: &(),
        stale: &(),
        dt_s: 1.0,
        inputs: &kelvin,
    };
    match with_cx(|cx| check_component(&line, &case, cx)) {
        Err(Failure::Contract(ContractViolation::NegativeProduction { .. })) => {}
        other => panic!("entropy destroyed at {kelvin:?} was not refused: {other:?}"),
    }
}

#[test]
fn the_matrix_refuses_what_is_not_a_conductance_between_two_traces() {
    let refused = |a: usize, b: usize, jump: &[f64], weights: &[f64], what: &str| {
        match ConductanceMatrix::new(a, b, jump.to_vec(), weights.to_vec()) {
            Err(error) => assert!(format!("{error:?}").contains(what), "{error:?} for {what}"),
            Ok(_) => panic!("{what}: accepted"),
        }
    };
    refused(0, 1, &[1.0], &[1.0], "port width");
    refused(1, 1, &[1.0, -1.0], &[1.0, 0.0], "weights entries");
    refused(1, 1, &[1.0, -1.0, 0.0], &[1.0], "jump entries");
    refused(1, 1, &[1.0, f64::NAN], &[1.0], "jump entry");
    refused(1, 1, &[1.0, -1.0], &[f64::INFINITY], "weights entry");
    // A row that does not sum to zero would make heat of a uniform temperature.
    refused(1, 1, &[1.0, -0.999], &[1.0], "jump row sum");
    let two = [1.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, -1.0];
    assert!(matches!(
        ConductanceMatrix::new(2, 2, two.to_vec(), vec![2.0, 1.0, 1.5, 2.0]),
        Err(ElementError::Asymmetric { row: 1, column: 0 })
    ));
    assert!(matches!(
        ConductanceMatrix::new(2, 2, two.to_vec(), vec![1.0, 2.0, 2.0, 1.0]),
        Err(ElementError::NotPositiveDefinite { .. })
    ));
    refused(1, 1, &[1.0, -1.0], &[-1.0], "NotPositiveDefinite");
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
    let weights = vec![0.8, 0.3, 0.3, 0.5];
    let jump = vec![1.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, -1.0];
    let wall = ConductanceMatrix::new(2, 2, jump, weights).unwrap();
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
