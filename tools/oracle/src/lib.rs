//! Independent MPFR audits for FrankenSim's deterministic elementary functions
//! and outward-rounded interval arithmetic.
//!
//! The crate is a standalone development tool. It is deliberately excluded
//! from the production workspace and makes no runtime-authority claim.

use std::fmt::Write as _;

use fs_ivl::Interval;
use fs_math::{det, ulp_distance};
use rug::{Float, float::Round, ops::Pow};

/// Minimum admitted oracle precision.
pub const MIN_PRECISION_BITS: u32 = 200;
/// Default MPFR precision used by the DSR lane.
pub const DEFAULT_PRECISION_BITS: u32 = 256;
/// Default deterministic sample count per audited family.
pub const DEFAULT_SAMPLES: usize = 4_096;
/// Stable seed for the full audit input generator.
pub const AUDIT_SEED: u64 = 0xF51A_0A0C_1E20_0260;

/// Configuration for one deterministic oracle run.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AuditConfig {
    /// Number of generated samples per family, in addition to fixed edges.
    pub samples: usize,
    /// MPFR result precision in bits.
    pub precision_bits: u32,
}

impl Default for AuditConfig {
    fn default() -> Self {
        Self {
            samples: DEFAULT_SAMPLES,
            precision_bits: DEFAULT_PRECISION_BITS,
        }
    }
}

impl AuditConfig {
    /// Refuse configurations that could weaken the external-oracle contract.
    pub fn validate(self) -> Result<Self, String> {
        if self.samples == 0 {
            return Err("oracle sample count must be nonzero".to_owned());
        }
        if self.precision_bits < MIN_PRECISION_BITS {
            return Err(format!(
                "oracle precision must be at least {MIN_PRECISION_BITS} bits"
            ));
        }
        Ok(self)
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct UlpHistogram {
    bins: [u64; 9],
}

impl UlpHistogram {
    const LABELS: [&'static str; 9] = ["0", "1", "2", "3", "4", "5-8", "9-16", "17-64", "65+"];

    const fn new() -> Self {
        Self { bins: [0; 9] }
    }

    fn observe(&mut self, ulps: u64) {
        let index = match ulps {
            0..=4 => usize::try_from(ulps).expect("0..=4 fits usize"),
            5..=8 => 5,
            9..=16 => 6,
            17..=64 => 7,
            _ => 8,
        };
        self.bins[index] += 1;
    }

    fn render_json(&self) -> String {
        let mut out = String::from("{");
        for (index, (label, count)) in Self::LABELS.iter().zip(self.bins).enumerate() {
            if index != 0 {
                out.push(',');
            }
            write!(&mut out, "\"{label}\":{count}").expect("String writes cannot fail");
        }
        out.push('}');
        out
    }
}

/// One per-function accuracy verdict.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct UlpAuditRow {
    function: &'static str,
    samples: u64,
    skipped: u64,
    failures: u64,
    max_ulp: u64,
    budget_at_argmax: u64,
    max_declared_budget: u64,
    arg0_bits: u64,
    arg1_bits: Option<u64>,
    histogram: UlpHistogram,
    first_failure: Option<String>,
}

impl UlpAuditRow {
    fn new(function: &'static str) -> Self {
        Self {
            function,
            samples: 0,
            skipped: 0,
            failures: 0,
            max_ulp: 0,
            budget_at_argmax: 0,
            max_declared_budget: 0,
            arg0_bits: 0,
            arg1_bits: None,
            histogram: UlpHistogram::new(),
            first_failure: None,
        }
    }

    fn observe(
        &mut self,
        arg0: f64,
        arg1: Option<f64>,
        actual: f64,
        reference: &Float,
        budget: u64,
    ) {
        let expected = reference.to_f64_round(Round::Nearest);
        self.samples += 1;
        self.max_declared_budget = self.max_declared_budget.max(budget);

        let ulps = if actual.is_nan() || expected.is_nan() {
            if actual.is_nan() && expected.is_nan() {
                0
            } else {
                u64::MAX
            }
        } else if actual.is_infinite() || expected.is_infinite() {
            if actual.to_bits() == expected.to_bits() {
                0
            } else {
                u64::MAX
            }
        } else {
            ulp_distance(actual, expected)
        };

        if self.samples == 1 || ulps > self.max_ulp {
            self.max_ulp = ulps;
            self.budget_at_argmax = budget;
            self.arg0_bits = arg0.to_bits();
            self.arg1_bits = arg1.map(f64::to_bits);
        }
        self.histogram.observe(ulps);
        if ulps > budget {
            self.failures += 1;
            if self.first_failure.is_none() {
                let inputs = arg1.map_or_else(
                    || format!("{:016x}", arg0.to_bits()),
                    |second| format!("{:016x},{:016x}", arg0.to_bits(), second.to_bits()),
                );
                self.first_failure = Some(format!(
                    "inputs={inputs},actual=0x{:016x},expected=0x{:016x},ulp={ulps},budget={budget}",
                    actual.to_bits(),
                    expected.to_bits()
                ));
            }
        }
    }

    fn skip(&mut self) {
        self.skipped += 1;
    }

    /// Whether every observed result stayed inside its declared budget.
    #[must_use]
    pub const fn passed(&self) -> bool {
        self.failures == 0
    }

    /// Render the row as deterministic JSON Lines output.
    #[must_use]
    pub fn render_json(&self, precision_bits: u32) -> String {
        let arg1 = self
            .arg1_bits
            .map_or_else(|| "null".to_owned(), |bits| format!("\"{bits:016x}\""));
        let failure = self
            .first_failure
            .as_ref()
            .map_or_else(|| "null".to_owned(), |text| format!("\"{text}\""));
        format!(
            "{{\"check\":\"high-precision-ulp\",\"function\":\"{}\",\"status\":\"{}\",\"precision_bits\":{precision_bits},\"samples\":{},\"skipped\":{},\"max_ulp\":{},\"budget_at_argmax\":{},\"max_declared_budget\":{},\"arg0_bits\":\"{:016x}\",\"arg1_bits\":{arg1},\"histogram\":{},\"failures\":{},\"first_failure\":{failure}}}",
            self.function,
            if self.passed() { "pass" } else { "fail" },
            self.samples,
            self.skipped,
            self.max_ulp,
            self.budget_at_argmax,
            self.max_declared_budget,
            self.arg0_bits,
            self.histogram.render_json(),
            self.failures,
        )
    }
}

/// One interval-containment verdict.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct IntervalAuditRow {
    operation: &'static str,
    samples: u64,
    failures: u64,
    first_failure: Option<String>,
}

impl IntervalAuditRow {
    const fn new(operation: &'static str) -> Self {
        Self {
            operation,
            samples: 0,
            failures: 0,
            first_failure: None,
        }
    }

    fn observe(&mut self, interval: Interval, reference: &Float, inputs: &[f64]) {
        self.samples += 1;
        let lo = Float::with_val(reference.prec(), interval.lo());
        let hi = Float::with_val(reference.prec(), interval.hi());
        let contains = lo <= *reference && *reference <= hi;
        if !contains {
            self.failures += 1;
            if self.first_failure.is_none() {
                let mut rendered_inputs = String::new();
                for (index, input) in inputs.iter().enumerate() {
                    if index != 0 {
                        rendered_inputs.push(',');
                    }
                    write!(&mut rendered_inputs, "{:016x}", input.to_bits())
                        .expect("String writes cannot fail");
                }
                self.first_failure = Some(format!(
                    "inputs={rendered_inputs},lo=0x{:016x},hi=0x{:016x},reference={}",
                    interval.lo().to_bits(),
                    interval.hi().to_bits(),
                    reference.to_string_radix(16, Some(70))
                ));
            }
        }
    }

    /// Whether every high-precision point result was enclosed.
    #[must_use]
    pub const fn passed(&self) -> bool {
        self.failures == 0
    }

    /// Render the row as deterministic JSON Lines output.
    #[must_use]
    pub fn render_json(&self, precision_bits: u32) -> String {
        let failure = self
            .first_failure
            .as_ref()
            .map_or_else(|| "null".to_owned(), |text| format!("\"{text}\""));
        format!(
            "{{\"check\":\"high-precision-interval\",\"operation\":\"{}\",\"status\":\"{}\",\"precision_bits\":{precision_bits},\"samples\":{},\"failures\":{},\"first_failure\":{failure}}}",
            self.operation,
            if self.passed() { "pass" } else { "fail" },
            self.samples,
            self.failures,
        )
    }
}

/// Complete deterministic audit report.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AuditReport {
    config: AuditConfig,
    ulp_rows: Vec<UlpAuditRow>,
    interval_rows: Vec<IntervalAuditRow>,
}

impl AuditReport {
    /// True only when every ULP and interval row passes.
    #[must_use]
    pub fn passed(&self) -> bool {
        self.ulp_rows.iter().all(UlpAuditRow::passed)
            && self.interval_rows.iter().all(IntervalAuditRow::passed)
    }

    /// Render all rows plus a terminal summary as deterministic JSON Lines.
    #[must_use]
    pub fn render_json_lines(&self) -> Vec<String> {
        let mut lines = Vec::with_capacity(self.ulp_rows.len() + self.interval_rows.len() + 1);
        lines.extend(
            self.ulp_rows
                .iter()
                .map(|row| row.render_json(self.config.precision_bits)),
        );
        lines.extend(
            self.interval_rows
                .iter()
                .map(|row| row.render_json(self.config.precision_bits)),
        );
        let total_samples = self
            .ulp_rows
            .iter()
            .map(|row| row.samples)
            .chain(self.interval_rows.iter().map(|row| row.samples))
            .sum::<u64>();
        let failures = self
            .ulp_rows
            .iter()
            .map(|row| row.failures)
            .chain(self.interval_rows.iter().map(|row| row.failures))
            .sum::<u64>();
        lines.push(format!(
            "{{\"check\":\"high-precision-oracle-summary\",\"status\":\"{}\",\"precision_bits\":{},\"generated_samples_per_family\":{},\"seed\":\"{AUDIT_SEED:016x}\",\"ulp_families\":{},\"interval_families\":{},\"total_observations\":{total_samples},\"failures\":{failures},\"no_claims\":[\"finite-sample-family\",\"no-exhaustive-binary64-proof\",\"no-run-authentication\",\"no-cross-isa-proof\"]}}",
            if self.passed() { "pass" } else { "fail" },
            self.config.precision_bits,
            self.config.samples,
            self.ulp_rows.len(),
            self.interval_rows.len(),
        ));
        lines
    }
}

#[derive(Clone, Copy)]
struct Lcg(u64);

impl Lcg {
    const fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }

    fn signed_unit(&mut self) -> f64 {
        let mantissa = self.next_u64() >> 11;
        let unit = mantissa as f64 * (1.0 / ((1_u64 << 53) as f64));
        unit.mul_add(2.0, -1.0)
    }

    fn positive_finite_bits(&mut self) -> f64 {
        loop {
            let bits = self.next_u64() & 0x7fff_ffff_ffff_ffff;
            let value = f64::from_bits(bits);
            if value.is_finite() && value > 0.0 {
                return value;
            }
        }
    }

    fn finite_bits(&mut self) -> f64 {
        loop {
            let value = f64::from_bits(self.next_u64());
            if value.is_finite() {
                return value;
            }
        }
    }
}

fn mp(value: f64, precision_bits: u32) -> Float {
    Float::with_val(precision_bits, value)
}

fn unary_edges(function: &str) -> &'static [f64] {
    const GENERAL: &[f64] = &[
        -0.0,
        0.0,
        f64::from_bits(1),
        -f64::from_bits(1),
        f64::MIN_POSITIVE,
        -f64::MIN_POSITIVE,
        -1.0,
        1.0,
        -0.5,
        0.5,
        -2.0,
        2.0,
    ];
    const POSITIVE: &[f64] = &[
        f64::from_bits(1),
        f64::MIN_POSITIVE,
        0.25,
        0.5,
        1.0,
        2.0,
        4.0,
        f64::MAX,
    ];
    const EXP: &[f64] = &[
        -745.0,
        -709.0,
        -1.0,
        -f64::EPSILON,
        -0.0,
        0.0,
        f64::EPSILON,
        1.0,
        709.0,
    ];
    const ERF: &[f64] = &[-6.0, -3.5, -3.0, -1.5, -0.0, 0.0, 1.5, 3.0, 3.5, 6.0];
    match function {
        "exp" | "expm1" => EXP,
        "ln" | "sqrt" => POSITIVE,
        "erf" => ERF,
        _ => GENERAL,
    }
}

fn generated_unary_input(function: &str, rng: &mut Lcg) -> f64 {
    match function {
        "exp" | "expm1" => rng.signed_unit().mul_add(727.0, -18.0),
        "ln" | "sqrt" => rng.positive_finite_bits(),
        "sin" | "cos" => rng.finite_bits(),
        "tan" => rng.finite_bits(),
        "erf" => rng.signed_unit() * 6.0,
        "tanh" => rng.signed_unit() * 30.0,
        _ => unreachable!("registered unary function"),
    }
}

fn audit_unary<B>(
    config: AuditConfig,
    function: &'static str,
    ours: fn(f64) -> f64,
    oracle: fn(f64, u32) -> Float,
    budget: B,
    seed: u64,
) -> UlpAuditRow
where
    B: Fn(f64) -> u64,
{
    let mut row = UlpAuditRow::new(function);
    let mut rng = Lcg::new(seed);
    for &input in unary_edges(function) {
        row.observe(
            input,
            None,
            ours(input),
            &oracle(input, config.precision_bits),
            budget(input),
        );
    }
    for _ in 0..config.samples {
        let input = generated_unary_input(function, &mut rng);
        let actual = ours(input);
        let reference = oracle(input, config.precision_bits);
        row.observe(input, None, actual, &reference, budget(input));
    }
    row
}

fn audit_atan2(config: AuditConfig) -> UlpAuditRow {
    let mut row = UlpAuditRow::new("atan2");
    let edges = [
        (-0.0, -1.0),
        (0.0, -1.0),
        (-1.0, 0.0),
        (1.0, 0.0),
        (1.0, 1.0),
        (1.0, -1.0),
        (-1.0, -1.0),
        (-1.0, 1.0),
    ];
    for (y, x) in edges {
        let reference = mp(y, config.precision_bits).atan2(&mp(x, config.precision_bits));
        row.observe(
            y,
            Some(x),
            det::atan2(y, x),
            &reference,
            det::ATAN_ULP_BUDGET + 1,
        );
    }
    let mut rng = Lcg::new(AUDIT_SEED ^ 0xA7A2);
    for _ in 0..config.samples {
        let y = rng.finite_bits();
        let x = rng.finite_bits();
        if x == 0.0 && y == 0.0 {
            row.skip();
            continue;
        }
        let reference = mp(y, config.precision_bits).atan2(&mp(x, config.precision_bits));
        row.observe(
            y,
            Some(x),
            det::atan2(y, x),
            &reference,
            det::ATAN_ULP_BUDGET + 1,
        );
    }
    row
}

fn audit_pow(config: AuditConfig) -> UlpAuditRow {
    let mut row = UlpAuditRow::new("pow");
    let edges = [
        (0.5, -2.0),
        (0.5, 0.5),
        (1.0, f64::MAX),
        (2.0, -10.0),
        (2.0, 10.0),
        (10.0, 3.25),
        (-2.0, -3.0),
        (-2.0, 3.0),
        (-2.0, 4.0),
    ];
    for (x, y) in edges {
        let reference = mp(x, config.precision_bits).pow(mp(y, config.precision_bits));
        row.observe(
            x,
            Some(y),
            det::pow(x, y),
            &reference,
            det::pow_ulp_budget(x, y),
        );
    }
    let mut rng = Lcg::new(AUDIT_SEED ^ 0xB0B0);
    for sample_index in 0..config.samples {
        let (x, y) = if sample_index % 2 == 0 {
            (
                rng.signed_unit().abs().mul_add(20.0, 0.001),
                rng.signed_unit() * 30.0,
            )
        } else {
            (rng.positive_finite_bits(), rng.signed_unit() * 4.0)
        };
        let reference = mp(x, config.precision_bits).pow(mp(y, config.precision_bits));
        let expected = reference.to_f64_round(Round::Nearest);
        if expected == 0.0 || !expected.is_finite() {
            row.skip();
            continue;
        }
        row.observe(
            x,
            Some(y),
            det::pow(x, y),
            &reference,
            det::pow_ulp_budget(x, y),
        );
    }
    row
}

fn exact_add(x: f64, y: f64, precision_bits: u32) -> Float {
    let mut value = mp(x, precision_bits);
    value += mp(y, precision_bits);
    value
}

fn exact_sub(x: f64, y: f64, precision_bits: u32) -> Float {
    let mut value = mp(x, precision_bits);
    value -= mp(y, precision_bits);
    value
}

fn exact_mul(x: f64, y: f64, precision_bits: u32) -> Float {
    let mut value = mp(x, precision_bits);
    value *= mp(y, precision_bits);
    value
}

fn exact_div(x: f64, y: f64, precision_bits: u32) -> Float {
    let mut value = mp(x, precision_bits);
    value /= mp(y, precision_bits);
    value
}

#[allow(clippy::too_many_lines)] // One table-like pass keeps all interval operations auditable together.
fn audit_intervals(config: AuditConfig) -> Vec<IntervalAuditRow> {
    let mut add = IntervalAuditRow::new("add");
    let mut sub = IntervalAuditRow::new("sub");
    let mut mul = IntervalAuditRow::new("mul");
    let mut div = IntervalAuditRow::new("div");
    let mut neg = IntervalAuditRow::new("neg");
    let mut abs = IntervalAuditRow::new("abs");
    let mut hull = IntervalAuditRow::new("hull");
    let mut intersect = IntervalAuditRow::new("intersect");
    let mut width = IntervalAuditRow::new("width-upper-bound");
    let mut exp = IntervalAuditRow::new("exp");
    let mut ln = IntervalAuditRow::new("ln");
    let mut sqrt = IntervalAuditRow::new("sqrt");
    let mut sin = IntervalAuditRow::new("sin");
    let mut cos = IntervalAuditRow::new("cos");
    let mut tanh = IntervalAuditRow::new("tanh");
    let mut rng = Lcg::new(AUDIT_SEED ^ 0x01A7_E2A1);

    for _ in 0..config.samples {
        let x0 = rng.signed_unit() * 1.0e150;
        let x1 = rng.signed_unit() * 1.0e150;
        let y0 = rng.signed_unit() * 1.0e150;
        let y1 = rng.signed_unit() * 1.0e150;
        let ix = Interval::new(x0.min(x1), x0.max(x1));
        let iy = Interval::new(y0.min(y1), y0.max(y1));
        let x_points = [ix.lo(), ix.midpoint(), ix.hi()];
        let y_points = [iy.lo(), iy.midpoint(), iy.hi()];

        width.observe(
            Interval::new(0.0, ix.width()),
            &exact_sub(ix.hi(), ix.lo(), config.precision_bits),
            &[ix.lo(), ix.hi()],
        );
        let combined_hull = ix.hull(iy);
        for point in x_points.into_iter().chain(y_points) {
            hull.observe(
                combined_hull,
                &mp(point, config.precision_bits),
                &[ix.lo(), ix.hi(), iy.lo(), iy.hi(), point],
            );
        }
        if let Some(overlap) = ix.intersect(iy) {
            for point in [overlap.lo(), overlap.midpoint(), overlap.hi()] {
                intersect.observe(
                    overlap,
                    &mp(point, config.precision_bits),
                    &[ix.lo(), ix.hi(), iy.lo(), iy.hi(), point],
                );
            }
        }

        for x in x_points {
            neg.observe(-ix, &-mp(x, config.precision_bits), &[ix.lo(), ix.hi(), x]);
            abs.observe(
                ix.abs(),
                &mp(x, config.precision_bits).abs(),
                &[ix.lo(), ix.hi(), x],
            );
            for y in y_points {
                let inputs = [ix.lo(), ix.hi(), iy.lo(), iy.hi(), x, y];
                add.observe(ix + iy, &exact_add(x, y, config.precision_bits), &inputs);
                sub.observe(ix - iy, &exact_sub(x, y, config.precision_bits), &inputs);
                mul.observe(ix * iy, &exact_mul(x, y, config.precision_bits), &inputs);
                if y != 0.0 {
                    div.observe(ix / iy, &exact_div(x, y, config.precision_bits), &inputs);
                }
            }
        }

        let exp0 = rng.signed_unit().mul_add(727.0, -18.0);
        let exp1 = rng.signed_unit().mul_add(727.0, -18.0);
        let exp_input = Interval::new(exp0.min(exp1), exp0.max(exp1));
        for x in [exp_input.lo(), exp_input.midpoint(), exp_input.hi()] {
            exp.observe(
                exp_input.exp(),
                &mp(x, config.precision_bits).exp(),
                &[exp_input.lo(), exp_input.hi(), x],
            );
        }

        let positive0 = rng.positive_finite_bits();
        let positive1 = rng.positive_finite_bits();
        let positive_input = Interval::new(positive0.min(positive1), positive0.max(positive1));
        for x in [
            positive_input.lo(),
            positive_input.midpoint(),
            positive_input.hi(),
        ] {
            ln.observe(
                positive_input.ln(),
                &mp(x, config.precision_bits).ln(),
                &[positive_input.lo(), positive_input.hi(), x],
            );
            sqrt.observe(
                positive_input.sqrt(),
                &mp(x, config.precision_bits).sqrt(),
                &[positive_input.lo(), positive_input.hi(), x],
            );
        }

        let trig0 = rng.signed_unit() * det::TRIG_DOMAIN;
        let trig1 = rng.signed_unit() * det::TRIG_DOMAIN;
        let trig_input = Interval::new(trig0.min(trig1), trig0.max(trig1));
        for x in [trig_input.lo(), trig_input.midpoint(), trig_input.hi()] {
            sin.observe(
                trig_input.sin(),
                &mp(x, config.precision_bits).sin(),
                &[trig_input.lo(), trig_input.hi(), x],
            );
            cos.observe(
                trig_input.cos(),
                &mp(x, config.precision_bits).cos(),
                &[trig_input.lo(), trig_input.hi(), x],
            );
        }

        let tanh0 = rng.signed_unit() * 30.0;
        let tanh1 = rng.signed_unit() * 30.0;
        let tanh_input = Interval::new(tanh0.min(tanh1), tanh0.max(tanh1));
        for x in [tanh_input.lo(), tanh_input.midpoint(), tanh_input.hi()] {
            tanh.observe(
                tanh_input.tanh(),
                &mp(x, config.precision_bits).tanh(),
                &[tanh_input.lo(), tanh_input.hi(), x],
            );
        }
    }

    vec![
        add, sub, mul, div, neg, abs, hull, intersect, width, exp, ln, sqrt, sin, cos, tanh,
    ]
}

/// Execute the complete deterministic high-precision comparison lane.
pub fn run_audit(config: AuditConfig) -> Result<AuditReport, String> {
    let config = config.validate()?;
    let constant_budget = |budget| move |_| budget;
    let ulp_rows = vec![
        audit_unary(
            config,
            "exp",
            det::exp,
            |x, p| mp(x, p).exp(),
            constant_budget(det::EXP_ULP_BUDGET),
            AUDIT_SEED ^ 0x01,
        ),
        audit_unary(
            config,
            "expm1",
            det::expm1,
            |x, p| mp(x, p).exp_m1(),
            constant_budget(det::EXPM1_ULP_BUDGET),
            AUDIT_SEED ^ 0x02,
        ),
        audit_unary(
            config,
            "ln",
            det::ln,
            |x, p| mp(x, p).ln(),
            constant_budget(det::LN_ULP_BUDGET),
            AUDIT_SEED ^ 0x03,
        ),
        audit_unary(
            config,
            "sin",
            det::sin,
            |x, p| mp(x, p).sin(),
            |x| {
                if x.abs() <= det::TRIG_DOMAIN {
                    det::SIN_ULP_BUDGET
                } else {
                    det::SIN_LARGE_ULP_BUDGET
                }
            },
            AUDIT_SEED ^ 0x04,
        ),
        audit_unary(
            config,
            "cos",
            det::cos,
            |x, p| mp(x, p).cos(),
            |x| {
                if x.abs() <= det::TRIG_DOMAIN {
                    det::COS_ULP_BUDGET
                } else {
                    det::SIN_LARGE_ULP_BUDGET
                }
            },
            AUDIT_SEED ^ 0x05,
        ),
        audit_unary(
            config,
            "tan",
            det::tan,
            |x, p| mp(x, p).tan(),
            constant_budget(det::TAN_ULP_BUDGET),
            AUDIT_SEED ^ 0x06,
        ),
        audit_atan2(config),
        audit_unary(
            config,
            "erf",
            det::erf,
            |x, p| mp(x, p).erf(),
            constant_budget(det::ERF_ULP_BUDGET),
            AUDIT_SEED ^ 0x07,
        ),
        audit_pow(config),
        audit_unary(
            config,
            "sqrt",
            det::sqrt,
            |x, p| mp(x, p).sqrt(),
            constant_budget(0),
            AUDIT_SEED ^ 0x08,
        ),
        audit_unary(
            config,
            "tanh",
            det::tanh,
            |x, p| mp(x, p).tanh(),
            constant_budget(det::TANH_ULP_BUDGET),
            AUDIT_SEED ^ 0x09,
        ),
    ];
    let interval_rows = audit_intervals(config);
    Ok(AuditReport {
        config,
        ulp_rows,
        interval_rows,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn configuration_refuses_weakened_or_empty_runs() {
        assert!(
            AuditConfig {
                samples: 0,
                precision_bits: DEFAULT_PRECISION_BITS,
            }
            .validate()
            .is_err()
        );
        assert!(
            AuditConfig {
                samples: 1,
                precision_bits: MIN_PRECISION_BITS - 1,
            }
            .validate()
            .is_err()
        );
    }

    #[test]
    fn mpfr_known_answers_pin_comparison_plumbing() {
        let p = DEFAULT_PRECISION_BITS;
        assert_eq!(mp(0.0, p).exp().to_f64().to_bits(), 1.0f64.to_bits());
        assert_eq!(mp(1.0, p).ln().to_f64().to_bits(), 0.0f64.to_bits());
        assert_eq!(mp(0.0, p).sin().to_f64().to_bits(), 0.0f64.to_bits());
        assert_eq!(mp(4.0, p).sqrt().to_f64().to_bits(), 2.0f64.to_bits());
        assert_eq!(mp(0.0, p).erf().to_f64().to_bits(), 0.0f64.to_bits());
    }

    #[test]
    fn comparison_plumbing_detects_a_corrupted_result() {
        let mut row = UlpAuditRow::new("corruption-probe");
        let reference = mp(1.0, DEFAULT_PRECISION_BITS);
        row.observe(
            1.0,
            None,
            f64::from_bits(1.0f64.to_bits() + 1),
            &reference,
            0,
        );
        assert!(!row.passed());
        assert_eq!(row.failures, 1);
    }

    #[test]
    fn interval_comparison_uses_exact_bounds_not_rounded_reference_only() {
        let mut row = IntervalAuditRow::new("containment-probe");
        let reference = mp(1.5, DEFAULT_PRECISION_BITS);
        row.observe(Interval::new(1.0, 2.0), &reference, &[1.5]);
        assert!(row.passed());
        row.observe(Interval::new(0.0, 1.0), &reference, &[1.5]);
        assert!(!row.passed());
    }

    #[test]
    fn small_reports_are_byte_stable() {
        let config = AuditConfig {
            samples: 16,
            precision_bits: MIN_PRECISION_BITS,
        };
        let first = run_audit(config)
            .expect("valid small audit")
            .render_json_lines();
        let second = run_audit(config)
            .expect("valid small audit")
            .render_json_lines();
        assert_eq!(first, second);
        assert!(run_audit(config).expect("valid small audit").passed());
    }
}
