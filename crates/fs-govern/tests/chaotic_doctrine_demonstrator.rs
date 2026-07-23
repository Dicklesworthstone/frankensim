//! G0/G3/G5 doctrine fixture for chaotic-system evidence routing.
//!
//! One fixed-step Lorenz-63 discretization supplies two deliberately different
//! exhibits. Outward-rounded interval propagation over the discrete RK4 map
//! retains its width trajectory until it becomes decision-useless, while a
//! bounded odd long-run observable over a symmetry-randomized synthetic
//! population feeds an anytime confidence sequence. The former becomes a
//! typed `NoUsefulBound`; it neither poisons nor upgrades the latter route.
//!
//! The interval exhibit encloses the stated discrete RK4 map only. It is not a
//! validated tube for the continuous ODE. The statistical exhibit validates
//! replay plumbing against a synthetic zero-mean truth; it does not validate
//! Lorenz-63 as a physical model, prove pseudorandom independence, or establish
//! coverage beyond the declared sub-Gaussian sampling law. The long-horizon
//! statistic remains `Estimated`, while the separately scoped local stability
//! and linearized-event intervals are `Verified` only for their stated
//! numerics. Routing creates no scientific authority.

use fs_eproc::GaussianMixtureCs;
use fs_evidence::{
    BoundInterval, BoundOutcome, ClaimClass, Color, ColorRank, EvidenceRegime, NoUsefulBoundCause,
    UsefulnessCriterion, validate_color_payload,
};
use fs_govern::{
    CERTIFICATE_REGIME_SCHEMA_VERSION, CLAIM_ROUTER_SCHEMA_VERSION, ChaosBasis, ClaimExtent,
    ClaimRequest, ClaimRouteDecision, ClaimRouteRefusalCause, DecisionNeed, DynamicsProfile,
    route_claim,
};
use fs_ivl::Interval;
use fs_rand::StreamKey;
use fs_report::{LabNotebook, no_useful_bound_markdown};

const FIXTURE_VERSION: &str = "frankensim-chaotic-doctrine-demonstrator-v1";
const SEED: u64 = 0xE09_0000_0000_0004;
const STREAM_KERNEL: u32 = 0xE094;
const DT: f64 = 0.01;
const LONG_HORIZON_STEPS: usize = 2_000;
const BURN_IN_STEPS: usize = 500;
const POPULATION_SIZE: usize = 32;
const PREDICTABILITY_HORIZON: f64 = 1.0;
const USEFUL_WIDTH: f64 = 16.0;
const MAX_ENCLOSURE_STEPS: usize = 2_000;
const SIGMA: f64 = 1.0;
const RHO: f64 = 8.0;
const ALPHA: f64 = 0.05;

#[derive(Debug, Clone, Copy, PartialEq)]
struct LorenzParameters {
    sigma: f64,
    rho: f64,
    beta: f64,
}

impl LorenzParameters {
    const CLASSIC: Self = Self {
        sigma: 10.0,
        rho: 28.0,
        beta: 8.0 / 3.0,
    };
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct PointState {
    x: f64,
    y: f64,
    z: f64,
}

impl PointState {
    fn add_scaled(self, derivative: Self, scale: f64) -> Self {
        Self {
            x: self.x + scale * derivative.x,
            y: self.y + scale * derivative.y,
            z: self.z + scale * derivative.z,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct IntervalState {
    x: Interval,
    y: Interval,
    z: Interval,
}

impl IntervalState {
    fn add_scaled(self, derivative: Self, scale: Interval) -> Self {
        Self {
            x: self.x + scale * derivative.x,
            y: self.y + scale * derivative.y,
            z: self.z + scale * derivative.z,
        }
    }

    fn widest_component(self) -> Interval {
        [self.x, self.y, self.z]
            .into_iter()
            .max_by(|left, right| left.width().total_cmp(&right.width()))
            .expect("three components")
    }
}

fn point_rhs(state: PointState, parameters: LorenzParameters) -> PointState {
    PointState {
        x: parameters.sigma * (state.y - state.x),
        y: state.x * (parameters.rho - state.z) - state.y,
        z: state.x * state.y - parameters.beta * state.z,
    }
}

fn point_rk4_step(state: PointState, parameters: LorenzParameters, dt: f64) -> PointState {
    let k1 = point_rhs(state, parameters);
    let k2 = point_rhs(state.add_scaled(k1, dt * 0.5), parameters);
    let k3 = point_rhs(state.add_scaled(k2, dt * 0.5), parameters);
    let k4 = point_rhs(state.add_scaled(k3, dt), parameters);
    PointState {
        x: state.x + dt * (k1.x + 2.0 * k2.x + 2.0 * k3.x + k4.x) / 6.0,
        y: state.y + dt * (k1.y + 2.0 * k2.y + 2.0 * k3.y + k4.y) / 6.0,
        z: state.z + dt * (k1.z + 2.0 * k2.z + 2.0 * k3.z + k4.z) / 6.0,
    }
}

fn interval_rhs(state: IntervalState, parameters: LorenzParameters) -> IntervalState {
    let sigma = Interval::point(parameters.sigma);
    let rho = Interval::point(parameters.rho);
    let beta = Interval::point(parameters.beta);
    IntervalState {
        x: sigma * (state.y - state.x),
        y: state.x * (rho - state.z) - state.y,
        z: state.x * state.y - beta * state.z,
    }
}

fn interval_rk4_step(state: IntervalState, parameters: LorenzParameters, dt: f64) -> IntervalState {
    let half_dt = Interval::point(dt * 0.5);
    let full_dt = Interval::point(dt);
    let k1 = interval_rhs(state, parameters);
    let k2 = interval_rhs(state.add_scaled(k1, half_dt), parameters);
    let k3 = interval_rhs(state.add_scaled(k2, half_dt), parameters);
    let k4 = interval_rhs(state.add_scaled(k3, full_dt), parameters);
    let two = Interval::point(2.0);
    let sixth_dt = Interval::point(dt / 6.0);
    IntervalState {
        x: state.x + sixth_dt * (k1.x + two * k2.x + two * k3.x + k4.x),
        y: state.y + sixth_dt * (k1.y + two * k2.y + two * k3.y + k4.y),
        z: state.z + sixth_dt * (k1.z + two * k2.z + two * k3.z + k4.z),
    }
}

fn initial_point(sample: usize, orientation: f64) -> PointState {
    let offset = sample as f64 * 0.000_125;
    PointState {
        x: orientation * (0.1 + offset),
        y: orientation * (0.05 - 0.5 * offset),
        z: 0.02 + 0.25 * offset,
    }
}

fn initial_box() -> IntervalState {
    let radius = 1.0e-12;
    IntervalState {
        x: Interval::new(0.1 - radius, 0.1 + radius),
        y: Interval::new(0.05 - radius, 0.05 + radius),
        z: Interval::new(0.02 - radius, 0.02 + radius),
    }
}

fn bounded_odd_score(value: f64) -> f64 {
    value / (1.0 + value.abs())
}

fn long_run_score(mut state: PointState, parameters: LorenzParameters) -> f64 {
    let mut sum = 0.0;
    for step in 0..LONG_HORIZON_STEPS {
        state = point_rk4_step(state, parameters, DT);
        assert!(
            state.x.is_finite() && state.y.is_finite() && state.z.is_finite(),
            "fixture left the finite RK4 domain"
        );
        if step >= BURN_IN_STEPS {
            sum += bounded_odd_score(state.x);
        }
    }
    sum / (LONG_HORIZON_STEPS - BURN_IN_STEPS) as f64
}

#[derive(Debug, Clone, PartialEq)]
struct EnclosureExhibit {
    width_bits: Vec<u64>,
    terminal_step: usize,
    terminal_interval: BoundInterval,
    outcome: BoundOutcome,
}

fn enclosure_exhibit(parameters: LorenzParameters) -> EnclosureExhibit {
    let mut state = initial_box();
    let mut width_bits = Vec::new();
    let mut terminal_step = 0;
    for step in 1..=MAX_ENCLOSURE_STEPS {
        state = interval_rk4_step(state, parameters, DT);
        let widest = state.widest_component();
        let width = widest.width();
        width_bits.push(width.to_bits());
        terminal_step = step;
        if !width.is_finite() || width > USEFUL_WIDTH {
            break;
        }
    }
    let widest = state.widest_component();
    let terminal_interval =
        BoundInterval::try_new(widest.lo(), widest.hi()).expect("ordered interval projection");
    let criterion = UsefulnessCriterion::try_new(
        "resolve one exact Lorenz-63 state beyond the admitted predictability horizon",
        "state-coordinate",
        USEFUL_WIDTH,
    )
    .expect("valid usefulness criterion");
    let outcome = BoundOutcome::classify(
        terminal_interval,
        criterion,
        NoUsefulBoundCause::LipschitzBlowup,
        ClaimClass::LongHorizonMeanLoad,
    );
    assert!(
        outcome.no_useful_bound().is_some(),
        "fixture must reach the declared useless-width boundary"
    );
    EnclosureExhibit {
        width_bits,
        terminal_step,
        terminal_interval,
        outcome,
    }
}

fn width_curve_discovery_digest(exhibit: &EnclosureExhibit) -> u64 {
    const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;
    let mut digest = FNV_OFFSET;
    for byte in FIXTURE_VERSION
        .as_bytes()
        .iter()
        .copied()
        .chain(DT.to_bits().to_le_bytes())
        .chain(USEFUL_WIDTH.to_bits().to_le_bytes())
        .chain((exhibit.terminal_step as u64).to_le_bytes())
        .chain(exhibit.terminal_interval.lower().to_bits().to_le_bytes())
        .chain(exhibit.terminal_interval.upper().to_bits().to_le_bytes())
        .chain(
            exhibit
                .width_bits
                .iter()
                .flat_map(|bits| bits.to_le_bytes()),
        )
    {
        digest ^= u64::from(byte);
        digest = digest.wrapping_mul(FNV_PRIME);
    }
    digest
}

fn route_exact_claim(duration: f64, unit: &str, predictability: f64) -> ClaimRouteDecision {
    let request = ClaimRequest::try_new(
        "chaotic-doctrine/exact-long-trajectory",
        ClaimClass::ExactLongChaoticTrajectory,
        "one exact Lorenz-63 RK4 state trajectory",
        ClaimExtent::try_long_horizon(duration, unit).expect("valid long horizon"),
        DecisionNeed::try_new(
            "resolve one exact state-coordinate trajectory",
            "state-coordinate",
            USEFUL_WIDTH,
        )
        .expect("valid decision need"),
        DynamicsProfile::new(
            false,
            true,
            false,
            ChaosBasis::try_declared(predictability, unit).expect("valid chaos declaration"),
        ),
        vec![
            "classic Lorenz-63 parameters are caller-declared".to_string(),
            "predictability horizon is a declared routing assumption".to_string(),
            "outward-rounded RK4 intervals enclose only the discrete map".to_string(),
        ],
    )
    .expect("valid exact-trajectory request");
    route_claim(request)
}

fn route_mean_claim(duration: f64, unit: &str) -> ClaimRouteDecision {
    let request = ClaimRequest::try_new(
        "chaotic-doctrine/long-run-odd-score",
        ClaimClass::LongHorizonMeanLoad,
        "mean bounded odd score over a symmetry-randomized Lorenz-63 population",
        ClaimExtent::try_long_horizon(duration, unit).expect("valid long horizon"),
        DecisionNeed::try_new(
            "distinguish the synthetic population mean from zero",
            "dimensionless-score",
            0.25,
        )
        .expect("valid decision need"),
        DynamicsProfile::new(
            false,
            true,
            false,
            ChaosBasis::try_declared(PREDICTABILITY_HORIZON, unit)
                .expect("valid chaos declaration"),
        ),
        vec![
            "one Philox orientation draw per deterministic base state".to_string(),
            "Lorenz symmetry makes each bounded odd score conditionally mean zero".to_string(),
            "score support is contained in [-1,1], so sigma=1 is sub-Gaussian".to_string(),
            "fixed-step RK4 and the synthetic model remain Estimated".to_string(),
        ],
    )
    .expect("valid mean request");
    route_claim(request)
}

fn route_local_stability_claim() -> ClaimRouteDecision {
    let request = ClaimRequest::try_new(
        "chaotic-doctrine/origin-local-stability",
        ClaimClass::LocalStability,
        "dominant real eigenvalue of the continuous Lorenz-63 origin Jacobian",
        ClaimExtent::Local,
        DecisionNeed::try_new(
            "classify the origin as locally linearly stable or unstable",
            "inverse-simulation-second",
            1.0e-9,
        )
        .expect("valid decision need"),
        DynamicsProfile::new(
            false,
            true,
            false,
            ChaosBasis::try_declared(PREDICTABILITY_HORIZON, "simulation-second")
                .expect("valid chaos declaration"),
        ),
        vec![
            "continuous Lorenz-63 origin Jacobian at classic parameters".to_string(),
            "positive enclosed eigenvalue establishes local linear instability only".to_string(),
            "no global attraction, nonlinear robustness, or trajectory claim".to_string(),
        ],
    )
    .expect("valid local-stability request");
    route_claim(request)
}

fn route_short_event_claim() -> ClaimRouteDecision {
    let request = ClaimRequest::try_new(
        "chaotic-doctrine/linearized-mode-threshold-event",
        ClaimClass::RootOrEventTime,
        "first tenfold-amplitude crossing of the origin's unstable linearized mode",
        ClaimExtent::try_finite_horizon(0.5, "simulation-second").expect("valid short horizon"),
        DecisionNeed::try_new(
            "bracket the local linearized-mode threshold crossing",
            "simulation-second",
            1.0e-9,
        )
        .expect("valid decision need"),
        DynamicsProfile::new(
            false,
            true,
            false,
            ChaosBasis::try_declared(PREDICTABILITY_HORIZON, "simulation-second")
                .expect("valid chaos declaration"),
        ),
        vec![
            "event model is the continuous Lorenz-63 origin linearization".to_string(),
            "initial modal amplitude and threshold ratio are exact declared scalars".to_string(),
            "event bracket does not transfer to the nonlinear long trajectory".to_string(),
        ],
    )
    .expect("valid short-event request");
    route_claim(request)
}

#[derive(Debug, Clone, PartialEq)]
struct LocalCertifiedExhibit {
    dominant_eigenvalue: BoundInterval,
    crossing_time: BoundInterval,
    stability_color: Color,
    event_color: Color,
}

fn local_certified_exhibit(parameters: LorenzParameters) -> LocalCertifiedExhibit {
    let sigma = Interval::point(parameters.sigma);
    let rho = Interval::point(parameters.rho);
    let one = Interval::point(1.0);
    let two = Interval::point(2.0);
    let four = Interval::point(4.0);
    let trace_magnitude = sigma + one;
    let discriminant = (sigma - one) * (sigma - one) + four * sigma * rho;
    let dominant = (-trace_magnitude + discriminant.sqrt()) / two;
    assert!(
        dominant.lo() > 0.0,
        "positive enclosed eigenvalue establishes local linear instability"
    );
    let dominant_eigenvalue =
        BoundInterval::try_new(dominant.lo(), dominant.hi()).expect("ordered spectral enclosure");
    let crossing = Interval::point(10.0).ln() / dominant;
    let crossing_time =
        BoundInterval::try_new(crossing.lo(), crossing.hi()).expect("ordered event-time bracket");
    assert!(
        crossing_time.lower() > 0.0 && crossing_time.upper() < 0.5,
        "linearized-mode event must stay inside the routed short horizon"
    );
    let stability_color = Color::Verified {
        lo: dominant_eigenvalue.lower(),
        hi: dominant_eigenvalue.upper(),
    };
    let event_color = Color::Verified {
        lo: crossing_time.lower(),
        hi: crossing_time.upper(),
    };
    validate_color_payload(&stability_color).expect("valid spectral enclosure color");
    validate_color_payload(&event_color).expect("valid event-time enclosure color");
    LocalCertifiedExhibit {
        dominant_eigenvalue,
        crossing_time,
        stability_color,
        event_color,
    }
}

#[derive(Debug, Clone, PartialEq)]
struct StatisticalExhibit {
    orientation_words: Vec<u64>,
    score_bits: Vec<u64>,
    interval_bits: Vec<(u64, u64)>,
    center: f64,
    radius: f64,
    color: Color,
}

fn statistical_exhibit(parameters: LorenzParameters) -> StatisticalExhibit {
    let mut stream = StreamKey {
        seed: SEED,
        kernel: STREAM_KERNEL,
        tile: 0,
    }
    .stream();
    let mut confidence_sequence = GaussianMixtureCs::new(SIGMA, RHO, ALPHA);
    let mut orientation_words = Vec::with_capacity(POPULATION_SIZE);
    let mut score_bits = Vec::with_capacity(POPULATION_SIZE);
    let mut interval_bits = Vec::with_capacity(POPULATION_SIZE);
    for sample in 0..POPULATION_SIZE {
        let word = stream.next_u64();
        let orientation = if word & 1 == 0 { -1.0 } else { 1.0 };
        let score = long_run_score(initial_point(sample, orientation), parameters);
        assert!((-1.0..=1.0).contains(&score));
        confidence_sequence.observe(score);
        let (center, radius) = confidence_sequence
            .interval()
            .expect("one observation establishes an interval");
        assert!(
            center - radius <= 0.0 && 0.0 <= center + radius,
            "the finite synthetic replay must retain its known zero mean"
        );
        orientation_words.push(word);
        score_bits.push(score.to_bits());
        interval_bits.push((center.to_bits(), radius.to_bits()));
    }
    let (center, radius) = confidence_sequence.interval().expect("non-empty study");
    let color = Color::Estimated {
        estimator: "lorenz-rk4-symmetry-randomized-time-mean-score-v1".to_string(),
        dispersion: radius,
    };
    validate_color_payload(&color)
        .expect("estimated result has a valid bounded identity and spread");
    StatisticalExhibit {
        orientation_words,
        score_bits,
        interval_bits,
        center,
        radius,
        color,
    }
}

#[derive(Debug, Clone, PartialEq)]
struct DoctrineStudy {
    enclosure: EnclosureExhibit,
    statistics: StatisticalExhibit,
    local: LocalCertifiedExhibit,
    exact_route: ClaimRouteDecision,
    mean_route: ClaimRouteDecision,
    stability_route: ClaimRouteDecision,
    event_route: ClaimRouteDecision,
    report_markdown: String,
    report_hash: u64,
}

fn build_report(
    enclosure: &EnclosureExhibit,
    statistics: &StatisticalExhibit,
    local: &LocalCertifiedExhibit,
    exact_route: &ClaimRouteDecision,
    mean_route: &ClaimRouteDecision,
    stability_route: &ClaimRouteDecision,
    event_route: &ClaimRouteDecision,
) -> LabNotebook {
    let refusal = enclosure
        .outcome
        .no_useful_bound()
        .expect("paired report requires the refusal exhibit");
    let mut notebook = LabNotebook::new(
        "Chaotic-system doctrine demonstrator",
        SEED,
        FIXTURE_VERSION,
    );
    let width_curve = enclosure
        .width_bits
        .iter()
        .enumerate()
        .map(|(index, bits)| format!("{}:0x{bits:016x}", index + 1))
        .collect::<Vec<_>>()
        .join(",");
    let final_claims_table = format!(
        "| Claim | Doctrine row | Evidence/result | Color | Boundary |\n\
         | --- | --- | --- | --- | --- |\n\
         | Exact long trajectory | CR-08 | NoUsefulBound | none | discrete RK4 interval width became decision-useless |\n\
         | Long-run bounded odd score | CR-05 | GaussianMixtureCs | Estimated | synthetic sub-Gaussian sampling law; model unvalidated |\n\
         | Origin local stability | CR-04 | dominant eigenvalue in [{lambda_lo}, {lambda_hi}] 1/s | Verified | local continuous-model linearization only |\n\
         | Linearized-mode event | CR-01 | crossing time in [{event_lo}, {event_hi}] s | Verified | linearized mode only; no nonlinear transfer |\n\
         | Conserved quantity | not requested | not applicable | none | Lorenz-63 is dissipative; no invariant asserted |",
        lambda_lo = local.dominant_eigenvalue.lower(),
        lambda_hi = local.dominant_eigenvalue.upper(),
        event_lo = local.crossing_time.lower(),
        event_hi = local.crossing_time.upper(),
    );
    notebook
        .prose(
            "Paired result: the exact long chaotic trajectory is refused as \
             NoUsefulBound, while the distinct long-horizon mean observable \
             routes to statistical/model evidence. Neither route mints authority.",
        )
        .prose(exact_route.render_record())
        .prose(no_useful_bound_markdown(refusal))
        .prose(mean_route.render_record())
        .prose(stability_route.render_record())
        .prose(event_route.render_record())
        .prose(format!("width-curve-step-to-f64-bits={width_curve}"))
        .prose(format!(
            "width-stop=criterion-crossed step={} threshold={} state-coordinate cause={}",
            enclosure.terminal_step,
            USEFUL_WIDTH,
            refusal.cause().code(),
        ))
        .prose(final_claims_table)
        .prose(
            "Long-horizon statistical color: Estimated. The two Verified rows \
             certify only their retained interval numerics for the stated \
             continuous-model linearization. The physical model, pseudorandom \
             independence, nonlinear event transfer, and cross-ISA \
             floating-point replay are not validated here.",
        )
        .metric(
            "interval blowup step",
            enclosure.terminal_step as f64,
            "RK4-step",
        )
        .metric(
            "interval blowup time",
            enclosure.terminal_step as f64 * DT,
            "simulation-second",
        )
        .metric(
            "terminal interval width",
            refusal.width_achieved(),
            "state-coordinate",
        )
        .metric(
            "confidence-sequence center",
            statistics.center,
            "dimensionless-score",
        )
        .metric(
            "confidence-sequence radius",
            statistics.radius,
            "dimensionless-score",
        )
        .metric(
            "dominant origin eigenvalue lower bound",
            local.dominant_eigenvalue.lower(),
            "inverse-simulation-second",
        )
        .metric(
            "dominant origin eigenvalue upper bound",
            local.dominant_eigenvalue.upper(),
            "inverse-simulation-second",
        )
        .metric(
            "linearized-mode event lower bound",
            local.crossing_time.lower(),
            "simulation-second",
        )
        .metric(
            "linearized-mode event upper bound",
            local.crossing_time.upper(),
            "simulation-second",
        )
        .step(
            "five-explicits",
            vec![
                "units=simulation-second,state-coordinate,dimensionless-score".to_string(),
                format!("seeds=philox:{SEED}:kernel:{STREAM_KERNEL}"),
                format!(
                    "budgets=useful-width:{USEFUL_WIDTH},alpha:{ALPHA},max-enclosure-steps:{MAX_ENCLOSURE_STEPS}"
                ),
                format!(
                    "versions=fixture:{FIXTURE_VERSION},claim-router-schema:{CLAIM_ROUTER_SCHEMA_VERSION},certificate-regime-schema:{CERTIFICATE_REGIME_SCHEMA_VERSION}"
                ),
                "capabilities=fs-ivl-interval-rk4,fs-eproc-gaussian-mixture-cs,fs-rand-philox,fs-report-lab-notebook".to_string(),
            ],
        )
        .step(
            "lorenz-rk4-interval-width-trajectory",
            vec![
                format!("dt={DT}"),
                format!("max-steps={MAX_ENCLOSURE_STEPS}"),
                format!("useful-width={USEFUL_WIDTH}"),
                "arithmetic=fs-ivl-outward-rounded".to_string(),
                "scope=discrete-map-only".to_string(),
            ],
        )
        .step(
            "lorenz-rk4-symmetry-population",
            vec![
                format!("seed={SEED}"),
                format!("kernel={STREAM_KERNEL}"),
                format!("population={POPULATION_SIZE}"),
                format!("steps={LONG_HORIZON_STEPS}"),
                format!("burn-in={BURN_IN_STEPS}"),
                format!("sigma={SIGMA}"),
                format!("rho={RHO}"),
                format!("alpha={ALPHA}"),
            ],
        );
    notebook
}

fn run_study(parameters: LorenzParameters) -> DoctrineStudy {
    let enclosure = enclosure_exhibit(parameters);
    let statistics = statistical_exhibit(parameters);
    let local = local_certified_exhibit(parameters);
    let duration = LONG_HORIZON_STEPS as f64 * DT;
    let exact_route = route_exact_claim(duration, "simulation-second", PREDICTABILITY_HORIZON);
    let mean_route = route_mean_claim(duration, "simulation-second");
    let stability_route = route_local_stability_claim();
    let event_route = route_short_event_claim();
    let notebook = build_report(
        &enclosure,
        &statistics,
        &local,
        &exact_route,
        &mean_route,
        &stability_route,
        &event_route,
    );
    DoctrineStudy {
        enclosure,
        statistics,
        local,
        exact_route,
        mean_route,
        stability_route,
        event_route,
        report_markdown: notebook.render_markdown(),
        report_hash: notebook.content_hash(),
    }
}

#[test]
fn g0_paired_routes_keep_refusal_and_statistical_evidence_distinct() {
    let study = run_study(LorenzParameters::CLASSIC);
    let exact_refusal = study
        .exact_route
        .refusal()
        .expect("exact long chaotic path must refuse");
    assert_eq!(exact_refusal.row_id(), "CR-08");
    assert_eq!(
        exact_refusal.required_evidence(),
        EvidenceRegime::NoUsefulBound
    );
    assert_eq!(
        exact_refusal.cause(),
        &ClaimRouteRefusalCause::ExactLongChaoticTrajectoryHasNoUsefulRoute {
            requested: 20.0,
            predictability_horizon: PREDICTABILITY_HORIZON,
            unit: "simulation-second".to_string(),
        }
    );
    assert_eq!(
        study
            .enclosure
            .outcome
            .no_useful_bound()
            .expect("width blowup is typed")
            .cause(),
        NoUsefulBoundCause::LipschitzBlowup
    );

    let mean_route = study
        .mean_route
        .routed()
        .expect("the distinct mean claim must route");
    assert_eq!(mean_route.row_id(), "CR-05");
    assert_eq!(
        mean_route.evidence(),
        EvidenceRegime::StatisticalObservableWithModelEvidence
    );
    assert_eq!(study.statistics.color.rank(), ColorRank::Estimated);
    let stability_route = study
        .stability_route
        .routed()
        .expect("local-stability claim must route");
    assert_eq!(stability_route.row_id(), "CR-04");
    assert_eq!(
        stability_route.evidence(),
        EvidenceRegime::SpectralOrLyapunovCertificate
    );
    assert_eq!(study.local.stability_color.rank(), ColorRank::Verified);
    let event_route = study.event_route.routed().expect("short event must route");
    assert_eq!(event_route.row_id(), "CR-01");
    assert_eq!(
        event_route.evidence(),
        EvidenceRegime::IntervalRootOrTaylorEnclosure
    );
    assert_eq!(study.local.event_color.rank(), ColorRank::Verified);
    assert!(study.report_markdown.contains("### NoUsefulBound"));
    assert!(study.report_markdown.contains("doctrine-row=CR-05"));
    assert!(
        study
            .report_markdown
            .contains("Long-horizon statistical color: Estimated")
    );
    assert!(
        study
            .report_markdown
            .contains("Neither route mints authority")
    );
    assert!(
        study
            .report_markdown
            .contains("| Origin local stability | CR-04 |")
    );
    assert!(
        study
            .report_markdown
            .contains("| Linearized-mode event | CR-01 |")
    );
    assert!(
        study
            .report_markdown
            .contains("| Conserved quantity | not requested | not applicable |")
    );
    assert!(
        study
            .report_markdown
            .contains("width-stop=criterion-crossed")
    );
}

#[test]
fn g3_time_unit_rescaling_and_parameter_perturbation_preserve_doctrine() {
    let seconds = route_exact_claim(20.0, "simulation-second", 1.0);
    let milliseconds = route_exact_claim(20_000.0, "simulation-millisecond", 1_000.0);
    let seconds_refusal = seconds.refusal().expect("seconds refusal");
    let milliseconds_refusal = milliseconds.refusal().expect("milliseconds refusal");
    assert_eq!(seconds_refusal.row_id(), milliseconds_refusal.row_id());
    assert_eq!(
        seconds_refusal.required_evidence(),
        milliseconds_refusal.required_evidence()
    );
    assert_eq!(
        seconds_refusal.suggested_reformulation(),
        milliseconds_refusal.suggested_reformulation()
    );
    assert_eq!(
        seconds_refusal.cause().code(),
        milliseconds_refusal.cause().code()
    );

    let classic = enclosure_exhibit(LorenzParameters::CLASSIC);
    let perturbed = enclosure_exhibit(LorenzParameters {
        rho: 28.5,
        ..LorenzParameters::CLASSIC
    });
    assert_ne!(
        classic.width_bits, perturbed.width_bits,
        "parameter perturbation must move the numerical exhibit"
    );
    for exhibit in [&classic, &perturbed] {
        assert_eq!(
            exhibit
                .outcome
                .no_useful_bound()
                .expect("both trajectories cross the usefulness boundary")
                .suggested_reformulation(),
            ClaimClass::LongHorizonMeanLoad
        );
    }
}

#[test]
fn synthetic_symmetry_truth_is_covered_at_every_retained_checkpoint() {
    for sample in 0..8 {
        let positive = long_run_score(initial_point(sample, 1.0), LorenzParameters::CLASSIC);
        let negative = long_run_score(initial_point(sample, -1.0), LorenzParameters::CLASSIC);
        assert_eq!(
            positive.to_bits() ^ (1_u64 << 63),
            negative.to_bits(),
            "the represented RK4 map must preserve the Lorenz sign symmetry"
        );
    }

    let statistics = statistical_exhibit(LorenzParameters::CLASSIC);
    assert_eq!(statistics.interval_bits.len(), POPULATION_SIZE);
    for (center_bits, radius_bits) in statistics.interval_bits {
        let center = f64::from_bits(center_bits);
        let radius = f64::from_bits(radius_bits);
        assert!(center - radius <= 0.0 && 0.0 <= center + radius);
    }
}

#[test]
fn g5_same_isa_replay_is_bit_stable() {
    let first = run_study(LorenzParameters::CLASSIC);
    let second = run_study(LorenzParameters::CLASSIC);
    assert_eq!(first, second);
    assert_ne!(first.report_hash, 0);
    assert_eq!(
        first.enclosure.width_bits.len(),
        first.enclosure.terminal_step
    );
    println!(
        "{{\"suite\":\"fs-govern/chaotic-doctrine-demonstrator-v1\",\"curve_discovery_digest\":\"0x{:016x}\",\"terminal_step\":{},\"curve_points\":{},\"report_hash\":\"0x{:016x}\",\"status\":\"discovery-not-frozen\"}}",
        width_curve_discovery_digest(&first.enclosure),
        first.enclosure.terminal_step,
        first.enclosure.width_bits.len(),
        first.report_hash,
    );
}
