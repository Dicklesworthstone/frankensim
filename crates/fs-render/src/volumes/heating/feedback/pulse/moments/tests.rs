use super::*;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1.0e-300),
        "{a:e} != {b:e}"
    );
}

fn f(b: f64, x: f64) -> f64 {
    if b == 0.0 {
        1.0
    } else {
        x / (b + (1.0 - b) * x)
    }
}

fn quadrature(function: &impl Fn(f64) -> f64, left: f64, right: f64) -> f64 {
    fn split(
        function: &impl Fn(f64) -> f64,
        a: f64,
        b: f64,
        values: [f64; 3],
        whole: f64,
        tolerance: f64,
        depth: usize,
    ) -> f64 {
        let middle = 0.5 * (a + b);
        let q1 = function(0.5 * (a + middle));
        let q3 = function(0.5 * (middle + b));
        let low = (middle - a) * (values[0] + 4.0 * q1 + values[1]) / 6.0;
        let high = (b - middle) * (values[1] + 4.0 * q3 + values[2]) / 6.0;
        let difference = low + high - whole;
        if difference.abs() <= 15.0 * tolerance || depth == 0 {
            low + high + difference / 15.0
        } else {
            split(
                function,
                a,
                middle,
                [values[0], q1, values[1]],
                low,
                tolerance * 0.5,
                depth - 1,
            ) + split(
                function,
                middle,
                b,
                [values[1], q3, values[2]],
                high,
                tolerance * 0.5,
                depth - 1,
            )
        }
    }
    let values = [
        function(left),
        function(0.5 * (left + right)),
        function(right),
    ];
    let whole = (right - left) * (values[0] + 4.0 * values[1] + values[2]) / 6.0;
    split(
        function,
        left,
        right,
        values,
        whole,
        1.0e-13 * (right - left),
        40,
    ) / (right - left)
}

#[test]
fn linear_constant_and_saturation_moments() {
    for r in [0.0, 0.2, 0.75, 1.0 - f64::EPSILON, 1.0] {
        close(first_moment(1.0, r), (1.0 + r) / 2.0, 8.0e-15);
        close(
            product_moment(1.0, 1.0, r),
            (1.0 + r + r * r) / 3.0,
            8.0e-15,
        );
        assert_eq!(first_moment(0.0, r), 1.0);
        assert_eq!(product_moment(0.0, 0.0, r), 1.0);
        close(product_moment(0.0, 0.3, r), first_moment(0.3, r), 1.0e-15);
        for tiny in [1.0e-20, 1.0e-100, f64::MIN_POSITIVE, f64::from_bits(1)] {
            close(first_moment(tiny, r), 1.0, 8.0e-15);
            close(product_moment(tiny, tiny, r), 1.0, 8.0e-15);
        }
    }
}

#[test]
fn independent_time_quadrature_and_nearly_equal_shapes() {
    for a in [
        1.0e-4_f64,
        0.01,
        0.2,
        0.5,
        0.75,
        1.0 - 2.0_f64.powi(-40),
        1.0,
    ] {
        for b in [a, a * (1.0 - 1.0e-12), 0.3, 1.0] {
            for r in [0.0, 0.17, 0.9] {
                close(
                    first_moment(a, r),
                    quadrature(&|x| f(a, x), r, 1.0),
                    2.0e-11,
                );
                close(
                    product_moment(a, b, r),
                    quadrature(&|x| f(a, x) * f(b, x), r, 1.0),
                    2.0e-11,
                );
                close(product_moment(a, b, r), product_moment(b, a, r), 1.0e-15);
            }
        }
    }
}

const CELLS: [(f64, f64, f64); 4] = [
    (0.0, 0.03, 1.0),
    (0.1, 0.0, 0.4),
    (2.0, 0.4, 0.6),
    (0.02, 0.02, 2.0),
];

fn partition(power: f64, (alpha, beta, length): (f64, f64, f64)) -> (f64, f64, f64, f64) {
    let optical_depth = alpha * length;
    let opacity = -(-optical_depth).exp_m1();
    let phi = if optical_depth == 0.0 {
        1.0
    } else {
        opacity / optical_depth
    };
    let q = beta * power * length * phi;
    let linear_removed = power * opacity;
    let linear_remaining = power * (-optical_depth).exp();
    let removed = power * (opacity + q) / (1.0 + q);
    let remaining = linear_remaining / (1.0 + q);
    (removed, remaining, linear_removed, linear_remaining)
}

fn cascade(low: f64, high: f64) -> ([f64; 4], f64) {
    let mut peak = high;
    let mut shape = 1.0;
    let mut heat = [0.0; 4];
    let mut transmitted = 0.0;
    for (i, cell) in CELLS.into_iter().enumerate() {
        let (removed, remaining, linear_removed, linear_remaining) = partition(peak, cell);
        let result = mean_cell(
            low / high,
            shape,
            removed,
            remaining,
            linear_removed,
            linear_remaining,
        )
        .unwrap();
        heat[i] = result.absorbed_w;
        shape = result.next_shape;
        peak = remaining;
        transmitted = result.transmitted_w;
    }
    (heat, transmitted)
}

#[test]
fn heterogeneous_cascade_matches_independent_time_integrals() {
    let (heat, transmitted) = cascade(3.0, 9.0);
    for i in 0..=CELLS.len() {
        let expected = quadrature(
            &|x| {
                let mut power = 3.0 + 6.0 * x;
                let mut removed = 0.0;
                for cell in CELLS.into_iter().take((i + 1).min(CELLS.len())) {
                    let (loss, remaining, _, _) = partition(power, cell);
                    removed = loss;
                    power = remaining;
                }
                if i == CELLS.len() { power } else { removed }
            },
            0.0,
            1.0,
        );
        close(
            if i == CELLS.len() {
                transmitted
            } else {
                heat[i]
            },
            expected,
            2.0e-12,
        );
    }
    close(heat.iter().sum::<f64>() + transmitted, 6.0, 4.0e-15);
}

#[test]
fn interval_subdivision_and_narrow_nonzero_ramps() {
    let (whole, out) = cascade(0.0, 20.0);
    for pieces in [2, 7, 32] {
        let mut split = [0.0; 4];
        let mut split_out = 0.0;
        for part in 0..pieces {
            let a = 20.0 * f64::from(part) / f64::from(pieces);
            let b = 20.0 * f64::from(part + 1) / f64::from(pieces);
            let (heat, transmitted) = cascade(a, b);
            for i in 0..4 {
                split[i] += heat[i] / f64::from(pieces);
            }
            split_out += transmitted / f64::from(pieces);
        }
        for i in 0..4 {
            close(split[i], whole[i], 3.0e-13);
        }
        close(split_out, out, 3.0e-13);
    }
    let (near, near_out) = cascade(20.0 * (1.0 - f64::EPSILON), 20.0);
    let (constant, constant_out) = cascade(20.0, 20.0);
    for i in 0..4 {
        close(near[i], constant[i], 3.0e-15);
    }
    close(near_out, constant_out, 3.0e-15);
}

#[test]
fn weak_loss_opaque_limit_and_invalid_partitions() {
    let weak = mean_cell(0.2, 1.0, 1.0e-30, 1.0e200, 0.0, 1.0e200).unwrap();
    close(weak.absorbed_w, 1.0e-30 * 1.24 / 3.0, 5.0e-15);
    let saturated = mean_cell(0.0, 1.0, 1.0e300, 1.0e-300, 0.0, 1.0e300).unwrap();
    close(saturated.absorbed_w, 0.5e300, 5.0e-15);
    close(saturated.transmitted_w, 1.0e-300, 5.0e-15);
    let opaque = mean_cell(0.0, 1.0, 1.0e300, 0.0, 1.0e300, 0.0).unwrap();
    close(opaque.absorbed_w, 0.5e300, 5.0e-15);
    assert_eq!(opaque.transmitted_w, 0.0);
    assert!(mean_cell(-1.0, 1.0, 1.0, 1.0, 0.5, 1.5).is_err());
    assert!(mean_cell(0.0, 1.0, 0.1, 2.0, 0.2, 1.0).is_err());
    assert!(mean_cell(0.0, f64::NAN, 0.1, 0.0, 0.0, 0.0).is_err());
}
