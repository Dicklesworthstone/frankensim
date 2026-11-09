//! G1/G3: actual dispersive-load power linewidths from the original Maxwell law.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{
    HarmonicPowerParameterJet, HarmonicPowerSample, ResonanceControls, refine_power_bandwidth,
    refine_power_bandwidth_with_parameter_derivative,
    refine_power_extremum_with_parameter_derivative,
};
use fs_modal::harmonic::{
    HarmonicFrequencyControls, HarmonicFrequencyParameterControls, HarmonicPortControls,
    HarmonicPortFrequencyParameterDirection, HarmonicPortParameterDirection,
    HarmonicResponseControls, HarmonicResponseError, PreparedHarmonicResponse,
};
use fs_modal::quadratic_modes::QuadraticParameterDirection;

const M: f64 = 1.0;
const C: f64 = 0.0;
const K: f64 = 1.0;
const B: f64 = 1.0;
const KB: f64 = 11.0 / 8.0;
const TAU: f64 = 1.0;
const FORCE: C64 = C64::new(1.2, -0.7);

#[derive(Clone, Copy, Default)]
struct Direction {
    mass: f64,
    damping: f64,
    stiffness: f64,
    attachment: f64,
    load_stiffness: f64,
    relaxation: f64,
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 127,
                kernel_id: 41,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn query_controls() -> ResonanceControls {
    ResonanceControls {
        absolute_angular_tolerance_rad_s: 1e-12,
        relative_angular_tolerance: 1e-12,
        maximum_evaluations: 300,
    }
}

fn law(w: f64, d: Direction) -> [C64; 5] {
    let s = C64::new(1.0, -w * TAU);
    [
        C64::from_re(KB * TAU) / s,
        C64::new(0.0, KB * TAU * TAU) / (s * s),
        C64::from_re(-2.0 * KB * TAU * TAU * TAU) / (s * s * s),
        C64::from_re(d.load_stiffness * TAU) / s + C64::from_re(KB * d.relaxation) / (s * s),
        C64::new(0.0, d.load_stiffness * TAU * TAU) / (s * s)
            + C64::new(0.0, 2.0 * KB * TAU * d.relaxation) / (s * s * s),
    ]
}

fn owner(w: f64, cx: &Cx<'_>) -> Result<PreparedHarmonicResponse, HarmonicResponseError> {
    PreparedHarmonicResponse::try_with_port_impedance(
        &[M],
        &[C],
        &[K],
        1,
        &[B],
        1,
        &[law(w, Direction::default())[0]],
        w,
        HarmonicResponseControls {
            maximum_coordinates: 2,
            maximum_payload_bytes: 1 << 20,
            maximum_relative_residual: 2e-11,
        },
        HarmonicPortControls {
            maximum_ports: 2,
            maximum_projection_terms: 100_000,
            maximum_payload_bytes: 1 << 20,
        },
        cx,
    )
}

fn first(w: f64, cx: &Cx<'_>) -> Result<HarmonicPowerSample, HarmonicResponseError> {
    owner(w, cx)?
        .solve_with_frequency_derivative(
            &[FORCE],
            None,
            Some(&[law(w, Direction::default())[1]]),
            HarmonicFrequencyControls {
                maximum_additional_payload_bytes: 1 << 20,
                maximum_work: 1_000_000,
            },
            cx,
        )?
        .port_power_sample()
}

fn jet(
    w: f64,
    d: Direction,
    cx: &Cx<'_>,
) -> Result<HarmonicPowerParameterJet, HarmonicResponseError> {
    let z = law(w, d);
    owner(w, cx)?
        .solve_with_frequency_parameter_and_port_derivative(
            &[FORCE],
            QuadraticParameterDirection {
                mass: &[d.mass],
                damping: &[d.damping],
                stiffness: &[d.stiffness],
            },
            None,
            None,
            None,
            None,
            &[z[1]],
            &[z[2]],
            HarmonicPortFrequencyParameterDirection {
                parameter: HarmonicPortParameterDirection {
                    force_map: &[d.attachment],
                    impedance: &[z[3]],
                },
                impedance_angular_parameter_derivative: &[z[4]],
            },
            HarmonicFrequencyParameterControls {
                maximum_additional_payload_bytes: 1 << 20,
                maximum_work: 1_000_000,
            },
            cx,
        )?
        .port_power_parameter_jet()
}

/// Clear the ORIGINAL Maxwell memory denominator s=1-i*w*tau:
/// q=F*s / [(K-M*w^2-i*C*w)*s-i*w*B^2*Kb*tau].
/// Its load power is N*x/D(x), x=w^2, D a real cubic. This reference
/// never uses a mechanical solve, impedance jet or resonance query.
fn coefficients(d: Direction) -> ([f64; 4], [f64; 4], f64, f64) {
    let a = M + C * TAU;
    let at = d.mass + d.damping * TAU + C * d.relaxation;
    let b = C + TAU * (K + B * B * KB);
    let bt = d.damping
        + d.relaxation * (K + B * B * KB)
        + TAU * (d.stiffness + 2.0 * B * d.attachment * KB + B * B * d.load_stiffness);
    let c = M * TAU;
    let ct = d.mass * TAU + M * d.relaxation;
    let polynomial = [K * K, b * b - 2.0 * K * a, a * a - 2.0 * b * c, c * c];
    let tangent = [
        2.0 * K * d.stiffness,
        2.0 * b * bt - 2.0 * (d.stiffness * a + K * at),
        2.0 * a * at - 2.0 * (bt * c + b * ct),
        2.0 * c * ct,
    ];
    let n = 0.5 * FORCE.norm_sq() * B * B * KB * TAU;
    let nt = n * (2.0 * d.attachment / B + d.load_stiffness / KB + d.relaxation / TAU);
    (polynomial, tangent, n, nt)
}

fn polynomial(a: [f64; 4], x: f64) -> f64 {
    ((a[3] * x + a[2]) * x + a[1]) * x + a[0]
}

/// Return exact-model P, P_w and P_theta using the independent real cubic.
fn reference(w: f64, d: Direction) -> [f64; 3] {
    let (a, at, n, nt) = coefficients(d);
    let x = w * w;
    let den = polynomial(a, x);
    let den_t = polynomial(at, x);
    let slope_numerator = a[0] - a[2] * x * x - 2.0 * a[3] * x * x * x;
    [
        n * x / den,
        2.0 * w * n * slope_numerator / (den * den),
        nt * x / den - n * x * den_t / (den * den),
    ]
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn maxwell_memory_moves_actual_asymmetric_power_widths_and_quality_factors() {
    with_cx(|cx| {
        // At these original physical coefficients:
        // D(x)=1+(233/64)x-(15/4)x^2+x^3, with exact stationary x=2.
        // Half peak power has the cubic x^3-(15/4)x^2+(151/64)x+1=0.
        // These two positive-root square roots were evaluated independently
        // with 90-digit decimal arithmetic, not by the numerical query owner.
        let peak_w = 2.0_f64.sqrt();
        let lower_w = 1.126_258_885_627_684_5;
        let upper_w = 1.663_281_746_433_824;
        let width = upper_w - lower_w;
        assert!(((peak_w - lower_w) - (upper_w - peak_w)).abs() > 0.03);
        for d in [
            Direction {
                relaxation: 1.0,
                ..Direction::default()
            },
            Direction {
                attachment: 1.0,
                ..Direction::default()
            },
            Direction {
                load_stiffness: 1.0,
                ..Direction::default()
            },
            Direction {
                mass: 1.0,
                ..Direction::default()
            },
            Direction {
                damping: 1.0,
                ..Direction::default()
            },
            Direction {
                stiffness: 1.0,
                ..Direction::default()
            },
        ] {
            let peak = refine_power_extremum_with_parameter_derivative(
                (1.2, 1.6),
                query_controls(),
                cx,
                |w| first(w, cx),
                |w| jet(w, d, cx),
            )
            .unwrap();
            let band = refine_power_bandwidth_with_parameter_derivative(
                &peak,
                (0.5, peak.extremum.angular_frequency_rad_s),
                (peak.extremum.angular_frequency_rad_s, 2.5),
                query_controls(),
                cx,
                |w| first(w, cx),
                |w| jet(w, d, cx),
            )
            .unwrap();
            let (a, at, _, _) = coefficients(d);
            let x = 2.0;
            let g_t = at[0] - at[2] * x * x - 2.0 * at[3] * x * x * x;
            let g_x = -2.0 * a[2] * x - 6.0 * a[3] * x * x;
            let peak_t = -g_t / (2.0 * peak_w * g_x);
            let p = reference(peak_w, d);
            let lower = reference(lower_w, d);
            let upper = reference(upper_w, d);
            let lower_t = (0.5 * p[2] - lower[2]) / lower[1];
            let upper_t = (0.5 * p[2] - upper[2]) / upper[1];
            let width_t = upper_t - lower_t;
            let quality_t = peak_t / width - peak_w * width_t / (width * width);
            near(peak.extremum.angular_frequency_rad_s, peak_w, 3e-12);
            near(peak.angular_frequency_derivative_rad_s, peak_t, 2e-9);
            near(peak.average_power_derivative, p[2], 2e-9);
            near(band.bandwidth.lower_angular_frequency_rad_s, lower_w, 4e-12);
            near(band.bandwidth.upper_angular_frequency_rad_s, upper_w, 4e-12);
            near(band.bandwidth.angular_bandwidth_rad_s, width, 6e-12);
            near(
                band.bandwidth.response_quality_factor,
                peak_w / width,
                1e-10,
            );
            near(band.lower_angular_frequency_derivative_rad_s, lower_t, 2e-9);
            near(band.upper_angular_frequency_derivative_rad_s, upper_t, 2e-9);
            near(band.angular_bandwidth_derivative_rad_s, width_t, 3e-9);
            near(band.response_quality_factor_derivative, quality_t, 1e-8);
            near(lower[0] / p[0], 0.5, 2e-14);
            near(upper[0] / p[0], 0.5, 2e-14);
        }

        // The non-parameter entry measures the same original power line.
        let peak = refine_power_extremum_with_parameter_derivative(
            (1.2, 1.6),
            query_controls(),
            cx,
            |w| first(w, cx),
            |w| jet(w, Direction::default(), cx),
        )
        .unwrap();
        let band = refine_power_bandwidth(
            &peak.extremum,
            (0.5, peak.extremum.angular_frequency_rad_s),
            (peak.extremum.angular_frequency_rad_s, 2.5),
            query_controls(),
            cx,
            |w| first(w, cx),
        )
        .unwrap();
        near(band.angular_bandwidth_rad_s, width, 6e-12);

        // Flow magnitude retains the Maxwell |1-i*w*tau|^2 factor which
        // cancels from absorbed power. Their stationary conditions differ.
        let response = owner(peak_w, cx)
            .unwrap()
            .solve_with_frequency_derivative(
                &[FORCE],
                None,
                Some(&[law(peak_w, Direction::default())[1]]),
                HarmonicFrequencyControls {
                    maximum_additional_payload_bytes: 1 << 20,
                    maximum_work: 1_000_000,
                },
                cx,
            )
            .unwrap();
        let flow = response.port_flow_transfer(&[1.0]).unwrap();
        let flow_slope = (flow.amplitude.conj() * flow.angular_derivative).re;
        assert!(flow_slope > 0.1);
        assert!(
            response
                .port_power_sample()
                .unwrap()
                .angular_derivative
                .abs()
                < 2e-13
        );
    });
}
