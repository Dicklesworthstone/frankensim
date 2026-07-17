//! Retained exact conformal-shift and surrogate-escalation evidence.
//!
//! The battery exhausts every leave-one-out rank in two deliberately separated
//! residual regimes.  Mondrian calibration has exact nominal coverage in both
//! buckets, while pooling still has nominal marginal coverage and silently
//! undercovers the shifted bucket.  The same pooled width is cross-checked
//! through `fs-surrogate`, then the production certify-or-escalate policy is
//! exercised for calibrated, unseen, drifted, and recalibrated queries.
//!
//! This is a finite deterministic counterexample and replay receipt, not a
//! theorem for arbitrary distributions, shifts, ties, buckets, or callers.
#![cfg(feature = "conformal-hardening")]

use fs_eproc::hardening::{BucketBand, DriftMonitor, ExchangeabilityCard, MondrianConformal};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, EventKind, Severity};
use fs_rand::StreamKey;
use fs_surrogate::{ConformalBand, Decision, certify_or_escalate, conformal_band};

const SUITE: &str = "fs-eproc-conformal-shift";
const CASE: &str = "exact-mondrian-surrogate-escalation-v1";
const INPUT_SEED: u64 = 0xC0F0_4A11_5EED_0001;
const ROTATION_KERNEL: u32 = 0xE719;
const ROTATION_TILE: u32 = 0;

const REGIME_COUNT: usize = 2;
const ORBIT_SIZE: usize = 20;
const MIN_CALIBRATION: usize = ORBIT_SIZE - 1;
const ALPHA_NUMERATOR: u64 = 1;
const ALPHA_DENOMINATOR: u64 = 10;
const ALPHA: f64 = ALPHA_NUMERATOR as f64 / ALPHA_DENOMINATOR as f64;
const EXPECTED_MONDRIAN_HITS: usize = 18;
const EXPECTED_POOLED_SMOOTH_HITS: usize = 20;
const EXPECTED_POOLED_SHOCK_HITS: usize = 16;
const EXPECTED_POOLED_TOTAL_HITS: usize = 36;

const DRIFT_TRAIN_SIZE: usize = 64;
const DRIFT_ALPHA: f64 = 0.05;
const DRIFT_BUDGET: usize = 64;
const SHIFTED_QUERY: f64 = 2.0;
const FS_SURROGATE_API_SCHEMA: &str = "conformal-band+certify-or-escalate-v1";

const _: () = assert!(REGIME_COUNT == 2);
const _: () = assert!(ORBIT_SIZE == 20 && MIN_CALIBRATION == 19);
const _: () = assert!(ALPHA_NUMERATOR == 1 && ALPHA_DENOMINATOR == 10);
const _: () = assert!(EXPECTED_MONDRIAN_HITS * 10 == ORBIT_SIZE * 9);
const _: () = assert!(EXPECTED_POOLED_TOTAL_HITS == 2 * EXPECTED_MONDRIAN_HITS);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Regime {
    Smooth,
    Shock,
}

impl Regime {
    const fn name(self) -> &'static str {
        match self {
            Self::Smooth => "smooth",
            Self::Shock => "shock",
        }
    }

    const fn index(self) -> usize {
        match self {
            Self::Smooth => 0,
            Self::Shock => 1,
        }
    }

    fn residual(self, rank: usize) -> f64 {
        let one_based = u32::try_from(rank + 1).expect("fixture rank fits u32");
        match self {
            Self::Smooth => f64::from(one_based),
            Self::Shock => 100.0 + f64::from(one_based),
        }
    }
}

const REGIMES: [Regime; REGIME_COUNT] = [Regime::Smooth, Regime::Shock];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RetainedBand {
    Calibrated { half_width_bits: u64, n: usize },
    Refused { have: usize, need: usize },
}

impl RetainedBand {
    fn from_bucket(band: BucketBand) -> Self {
        match band {
            BucketBand::Calibrated { half_width, n } => Self::Calibrated {
                half_width_bits: half_width.to_bits(),
                n,
            },
            BucketBand::Refused { have, need } => Self::Refused { have, need },
        }
    }

    fn half_width(&self) -> Option<f64> {
        match *self {
            Self::Calibrated {
                half_width_bits, ..
            } => Some(f64::from_bits(half_width_bits)),
            Self::Refused { .. } => None,
        }
    }

    fn half_width_bits(&self) -> Option<u64> {
        match *self {
            Self::Calibrated {
                half_width_bits, ..
            } => Some(half_width_bits),
            Self::Refused { .. } => None,
        }
    }

    fn covers(&self, residual: f64) -> bool {
        self.half_width().is_some_and(|width| residual <= width)
    }

    fn policy_band(&self) -> ConformalBand {
        ConformalBand {
            half_width: self.half_width().unwrap_or(f64::INFINITY),
            alpha: ALPHA,
        }
    }

    fn decision_tolerance(&self) -> f64 {
        self.half_width()
            .filter(|width| width.is_finite())
            .unwrap_or(0.0)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum PolicyRecord {
    UseSurrogate { band_half_width_bits: u64 },
    Escalate { reason: String },
}

impl PolicyRecord {
    fn from_decision(decision: Decision) -> Self {
        match decision {
            Decision::UseSurrogate { band_half_width } => Self::UseSurrogate {
                band_half_width_bits: band_half_width.to_bits(),
            },
            Decision::Escalate { reason } => Self::Escalate { reason },
        }
    }

    fn uses_exact_width(&self, expected_bits: u64) -> bool {
        matches!(
            self,
            Self::UseSurrogate {
                band_half_width_bits
            } if *band_half_width_bits == expected_bits
        )
    }

    fn escalates_with(&self, expected: &str) -> bool {
        matches!(self, Self::Escalate { reason } if reason == expected)
    }
}

fn policy_for(band: &RetainedBand, in_validity_domain: bool) -> PolicyRecord {
    PolicyRecord::from_decision(certify_or_escalate(
        &band.policy_band(),
        in_validity_domain,
        band.decision_tolerance(),
    ))
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct HoldoutRecord {
    regime: Regime,
    rank: usize,
    residual_bits: u64,
    bucket_names_canonical: bool,
    bucket_band: RetainedBand,
    marginal_band: RetainedBand,
    surrogate_half_width_bits: u64,
    surrogate_alpha_bits: u64,
    bucket_covered: bool,
    marginal_covered: bool,
    bucket_policy: PolicyRecord,
    pooled_policy: PolicyRecord,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Campaign {
    regime_rotation: usize,
    selector_draws: u64,
    records: Vec<HoldoutRecord>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CoverageCounts {
    mondrian_smooth: usize,
    mondrian_shock: usize,
    pooled_smooth: usize,
    pooled_shock: usize,
}

impl CoverageCounts {
    const fn pooled_total(self) -> usize {
        self.pooled_smooth + self.pooled_shock
    }
}

impl Campaign {
    fn coverage_counts(&self) -> CoverageCounts {
        let count = |regime: Regime, mondrian: bool| {
            self.records
                .iter()
                .filter(|record| {
                    record.regime == regime
                        && if mondrian {
                            record.bucket_covered
                        } else {
                            record.marginal_covered
                        }
                })
                .count()
        };
        CoverageCounts {
            mondrian_smooth: count(Regime::Smooth, true),
            mondrian_shock: count(Regime::Shock, true),
            pooled_smooth: count(Regime::Smooth, false),
            pooled_shock: count(Regime::Shock, false),
        }
    }
}

fn run_campaign() -> Campaign {
    let mut selector = StreamKey {
        seed: INPUT_SEED,
        kernel: ROTATION_KERNEL,
        tile: ROTATION_TILE,
    }
    .stream();
    let regime_rotation =
        usize::try_from(selector.next_below(2)).expect("two-regime selector result fits usize");
    let selector_draws = selector.index();
    let rank_rotation = regime_rotation * (ORBIT_SIZE / REGIME_COUNT);
    let regime_order = if regime_rotation == 0 {
        REGIMES
    } else {
        [Regime::Shock, Regime::Smooth]
    };

    let mut records = Vec::with_capacity(REGIME_COUNT * ORBIT_SIZE);
    for regime in regime_order {
        for rank_offset in 0..ORBIT_SIZE {
            let rank = (rank_offset + rank_rotation) % ORBIT_SIZE;
            let residual = regime.residual(rank);
            let mut calibrator = MondrianConformal::new(MIN_CALIBRATION);
            let mut pooled_residuals = Vec::with_capacity(REGIME_COUNT * ORBIT_SIZE - 1);
            for candidate_regime in REGIMES {
                for candidate_rank in 0..ORBIT_SIZE {
                    if candidate_regime == regime && candidate_rank == rank {
                        continue;
                    }
                    let candidate_residual = candidate_regime.residual(candidate_rank);
                    calibrator.add(candidate_regime.name(), candidate_residual);
                    pooled_residuals.push(candidate_residual);
                }
            }

            let bucket_names_canonical = calibrator
                .bucket_names()
                .iter()
                .map(String::as_str)
                .eq(["shock", "smooth"]);
            let bucket_band = RetainedBand::from_bucket(calibrator.band(regime.name(), ALPHA));
            let marginal_band = RetainedBand::from_bucket(calibrator.marginal_band(ALPHA));
            let surrogate_band = conformal_band(&pooled_residuals, ALPHA);
            let bucket_covered = bucket_band.covers(residual);
            let marginal_covered = marginal_band.covers(residual);
            let bucket_policy = policy_for(&bucket_band, true);
            let pooled_policy = PolicyRecord::from_decision(certify_or_escalate(
                &surrogate_band,
                true,
                surrogate_band.half_width,
            ));

            records.push(HoldoutRecord {
                regime,
                rank,
                residual_bits: residual.to_bits(),
                bucket_names_canonical,
                bucket_band,
                marginal_band,
                surrogate_half_width_bits: surrogate_band.half_width.to_bits(),
                surrogate_alpha_bits: surrogate_band.alpha.to_bits(),
                bucket_covered,
                marginal_covered,
                bucket_policy,
                pooled_policy,
            });
        }
    }

    Campaign {
        regime_rotation,
        selector_draws,
        records,
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PolicyTrace {
    known_band: RetainedBand,
    known_policy: PolicyRecord,
    unseen_band: RetainedBand,
    unseen_policy: PolicyRecord,
    control_fired: bool,
    control_observations: usize,
    control_log_e_bits: u64,
    shifted_detection: u64,
    shifted_validity_scale_bits: u64,
    shifted_log_e_bits: u64,
    shifted_policy: PolicyRecord,
    ignore_drift_policy: PolicyRecord,
    recalibrated_band: RetainedBand,
    recalibrated_policy: PolicyRecord,
    exchangeability_card_json: String,
}

impl PolicyTrace {
    fn pass(&self) -> bool {
        let control_log_e = f64::from_bits(self.control_log_e_bits);
        let shifted_log_e = f64::from_bits(self.shifted_log_e_bits);
        let shifted_scale = f64::from_bits(self.shifted_validity_scale_bits);
        let known_width = self.known_band.half_width_bits();
        let recalibrated_width = self.recalibrated_band.half_width_bits();
        matches!(
            self.known_band,
            RetainedBand::Calibrated {
                n: MIN_CALIBRATION,
                ..
            }
        ) && known_width.is_some_and(|bits| self.known_policy.uses_exact_width(bits))
            && self.unseen_band
                == RetainedBand::Refused {
                    have: 0,
                    need: MIN_CALIBRATION,
                }
            && self
                .unseen_policy
                .escalates_with("conformal band is unbounded")
            && !self.control_fired
            && self.control_observations == DRIFT_TRAIN_SIZE
            && control_log_e.is_finite()
            && (1..=usize_u64(DRIFT_BUDGET)).contains(&self.shifted_detection)
            && shifted_log_e.is_finite()
            && shifted_scale.is_finite()
            && shifted_scale > 0.0
            && shifted_scale < 1.0
            && self
                .shifted_policy
                .escalates_with("query outside the surrogate's validity domain")
            && known_width.is_some_and(|bits| self.ignore_drift_policy.uses_exact_width(bits))
            && matches!(
                self.recalibrated_band,
                RetainedBand::Calibrated {
                    n: MIN_CALIBRATION,
                    ..
                }
            )
            && recalibrated_width
                .is_some_and(|bits| self.recalibrated_policy.uses_exact_width(bits))
            && self.exchangeability_card_json
                == "{\"bucketing\":\"regime-class\",\"drift_alpha\":0.05,\"fcr_budget\":0.1,\"refresh_policy\":\"detect-refuse-then-recalibrate\"}"
    }
}

fn drift_training_value(index: usize) -> f64 {
    f64::from(u32::try_from(index).expect("drift fixture index fits u32"))
        / f64::from(u32::try_from(DRIFT_TRAIN_SIZE).expect("training size fits u32"))
}

fn run_policy_trace() -> PolicyTrace {
    let mut calibrator = MondrianConformal::new(MIN_CALIBRATION);
    for rank in 0..MIN_CALIBRATION {
        calibrator.add(Regime::Smooth.name(), Regime::Smooth.residual(rank));
    }
    let known_band = RetainedBand::from_bucket(calibrator.band(Regime::Smooth.name(), ALPHA));
    let known_policy = policy_for(&known_band, true);
    let unseen_band = RetainedBand::from_bucket(calibrator.band(Regime::Shock.name(), ALPHA));
    let unseen_policy = policy_for(&unseen_band, false);

    let training: Vec<_> = (0..DRIFT_TRAIN_SIZE).map(drift_training_value).collect();
    let mut control_monitor = DriftMonitor::new(training.clone(), DRIFT_ALPHA);
    let mut control_fired = false;
    let mut control_observations = 0usize;
    for offset in 0..(DRIFT_TRAIN_SIZE / 2) {
        for index in [
            DRIFT_TRAIN_SIZE / 2 - 1 - offset,
            DRIFT_TRAIN_SIZE / 2 + offset,
        ] {
            let verdict = control_monitor.observe(drift_training_value(index));
            control_fired |= verdict.drifted;
            control_observations += 1;
        }
    }
    let control_log_e_bits = control_monitor.log_e_value().to_bits();

    let mut shifted_monitor = DriftMonitor::new(training, DRIFT_ALPHA);
    let mut shifted_detection = 0u64;
    let mut shifted_validity_scale_bits = 1.0f64.to_bits();
    for _ in 0..DRIFT_BUDGET {
        let verdict = shifted_monitor.observe(SHIFTED_QUERY);
        if verdict.drifted {
            shifted_detection = verdict.samples_at_detection;
            shifted_validity_scale_bits = verdict.validity_scale.to_bits();
            break;
        }
    }
    let shifted_log_e_bits = shifted_monitor.log_e_value().to_bits();
    let shifted_policy = policy_for(&known_band, shifted_detection == 0);
    let ignore_drift_policy = policy_for(&known_band, true);

    let mut recalibrator = MondrianConformal::new(MIN_CALIBRATION);
    for rank in 0..MIN_CALIBRATION {
        recalibrator.add(Regime::Shock.name(), Regime::Shock.residual(rank));
    }
    let recalibrated_band =
        RetainedBand::from_bucket(recalibrator.band(Regime::Shock.name(), ALPHA));
    let recalibrated_policy = policy_for(&recalibrated_band, true);
    let exchangeability_card_json = ExchangeabilityCard {
        bucketing: "regime-class".to_string(),
        drift_alpha: DRIFT_ALPHA,
        fcr_budget: ALPHA,
        refresh_policy: "detect-refuse-then-recalibrate".to_string(),
    }
    .to_json();

    PolicyTrace {
        known_band,
        known_policy,
        unseen_band,
        unseen_policy,
        control_fired,
        control_observations,
        control_log_e_bits,
        shifted_detection,
        shifted_validity_scale_bits,
        shifted_log_e_bits,
        shifted_policy,
        ignore_drift_policy,
        recalibrated_band,
        recalibrated_policy,
        exchangeability_card_json,
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Evidence {
    campaign: Campaign,
    policy: PolicyTrace,
}

fn run_evidence() -> Evidence {
    Evidence {
        campaign: run_campaign(),
        policy: run_policy_trace(),
    }
}

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("fixture cardinality fits u64")
}

fn config_identity() -> ReplayIdentity {
    IdentityBuilder::new("fs-eproc-conformal-shift-config-v1")
        .str("units", "dimensionless-absolute-residual")
        .u64("input-seed", INPUT_SEED)
        .u64("rotation-kernel", u64::from(ROTATION_KERNEL))
        .u64("rotation-tile", u64::from(ROTATION_TILE))
        .str("rotation-draw", "next-below-two-once")
        .str("regime-order", "smooth,shock-cyclic")
        .str("rank-order", "rank+rotation*10-modulo-20")
        .u64(
            "stream-semantics-version",
            u64::from(fs_rand::STREAM_SEMANTICS_VERSION),
        )
        .u64("regime-count", usize_u64(REGIME_COUNT))
        .u64("orbit-size-per-regime", usize_u64(ORBIT_SIZE))
        .u64("minimum-calibration", usize_u64(MIN_CALIBRATION))
        .u64("alpha-numerator", ALPHA_NUMERATOR)
        .u64("alpha-denominator", ALPHA_DENOMINATOR)
        .f64_bits("alpha", ALPHA)
        .str("smooth-residual-lattice", "rank-one-through-rank-twenty")
        .str("shock-residual-lattice", "one-hundred-plus-smooth-rank")
        .str(
            "holdout-rule",
            "exhaustive-leave-one-out-with-other-bucket-full",
        )
        .str("quantile-rule", "ceil((n+1)*(1-alpha))-th-order-statistic")
        .u64(
            "expected-mondrian-hits-per-bucket",
            usize_u64(EXPECTED_MONDRIAN_HITS),
        )
        .u64(
            "expected-pooled-smooth-hits",
            usize_u64(EXPECTED_POOLED_SMOOTH_HITS),
        )
        .u64(
            "expected-pooled-shock-hits",
            usize_u64(EXPECTED_POOLED_SHOCK_HITS),
        )
        .u64(
            "expected-pooled-total-hits",
            usize_u64(EXPECTED_POOLED_TOTAL_HITS),
        )
        .u64("drift-training-size", usize_u64(DRIFT_TRAIN_SIZE))
        .str("drift-training-grid", "i/64-for-i-zero-through-63")
        .str("drift-control-order", "31,32,30,33,...,0,63")
        .f64_bits("drift-alpha", DRIFT_ALPHA)
        .u64("drift-observation-budget", usize_u64(DRIFT_BUDGET))
        .f64_bits("shifted-query", SHIFTED_QUERY)
        .str("refusal-policy", "unbounded-band-forces-escalation")
        .str("drift-policy", "drifted-means-outside-validity-domain")
        .str("refresh-policy", "detect-refuse-then-recalibrate")
        .str("pooled-mutant", "pool-buckets-and-ignore-drift")
        .str("fs-surrogate-api-schema", FS_SURROGATE_API_SCHEMA)
        .str(
            "capabilities",
            "safe-rust;exact-finite-orbit;production-policy-seam",
        )
        .str("execution-context", "synchronous-direct-test-no-Cx")
        .flag("conformal-hardening-feature", true)
        .str("fs-eproc-version", fs_eproc::VERSION)
        .str("fs-math-version", fs_math::VERSION)
        .str("fs-rand-version", fs_rand::VERSION)
        .str("fs-obs-version", fs_obs::VERSION)
        .finish()
}

fn bind_band(mut builder: IdentityBuilder, role: &str, band: &RetainedBand) -> IdentityBuilder {
    builder = builder.str("band-role", role);
    match *band {
        RetainedBand::Calibrated { half_width_bits, n } => builder
            .str("band-kind", "calibrated")
            .u64("band-half-width-bits", half_width_bits)
            .u64("band-calibration-count", usize_u64(n)),
        RetainedBand::Refused { have, need } => builder
            .str("band-kind", "refused")
            .u64("band-have", usize_u64(have))
            .u64("band-need", usize_u64(need)),
    }
}

fn bind_policy(mut builder: IdentityBuilder, role: &str, policy: &PolicyRecord) -> IdentityBuilder {
    builder = builder.str("policy-role", role);
    match policy {
        PolicyRecord::UseSurrogate {
            band_half_width_bits,
        } => builder
            .str("policy-decision", "use-surrogate")
            .u64("policy-half-width-bits", *band_half_width_bits),
        PolicyRecord::Escalate { reason } => builder
            .str("policy-decision", "escalate")
            .str("policy-reason", reason),
    }
}

fn campaign_identity(config: &ReplayIdentity, campaign: &Campaign) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-eproc-conformal-shift-campaign-v1")
        .child("config", config)
        .u64("regime-rotation", usize_u64(campaign.regime_rotation))
        .u64("selector-draws", campaign.selector_draws)
        .u64("record-count", usize_u64(campaign.records.len()));
    for record in &campaign.records {
        builder = builder
            .str("regime", record.regime.name())
            .u64("rank", usize_u64(record.rank))
            .u64("residual-bits", record.residual_bits)
            .flag("bucket-names-canonical", record.bucket_names_canonical);
        builder = bind_band(builder, "mondrian", &record.bucket_band);
        builder = bind_band(builder, "pooled", &record.marginal_band);
        builder = builder
            .u64(
                "surrogate-half-width-bits",
                record.surrogate_half_width_bits,
            )
            .u64("surrogate-alpha-bits", record.surrogate_alpha_bits)
            .flag("bucket-covered", record.bucket_covered)
            .flag("marginal-covered", record.marginal_covered);
        builder = bind_policy(builder, "mondrian", &record.bucket_policy);
        builder = bind_policy(builder, "pooled", &record.pooled_policy);
    }
    builder.finish()
}

fn policy_identity(config: &ReplayIdentity, policy: &PolicyTrace) -> ReplayIdentity {
    let mut builder =
        IdentityBuilder::new("fs-eproc-conformal-shift-policy-v1").child("config", config);
    builder = bind_band(builder, "known", &policy.known_band);
    builder = bind_policy(builder, "known", &policy.known_policy);
    builder = bind_band(builder, "unseen", &policy.unseen_band);
    builder = bind_policy(builder, "unseen", &policy.unseen_policy);
    builder = builder
        .flag("control-fired", policy.control_fired)
        .u64(
            "control-observations",
            usize_u64(policy.control_observations),
        )
        .u64("control-log-e-bits", policy.control_log_e_bits)
        .u64("shifted-detection", policy.shifted_detection)
        .u64(
            "shifted-validity-scale-bits",
            policy.shifted_validity_scale_bits,
        )
        .u64("shifted-log-e-bits", policy.shifted_log_e_bits);
    builder = bind_policy(builder, "shifted", &policy.shifted_policy);
    builder = bind_policy(builder, "ignore-drift-mutant", &policy.ignore_drift_policy);
    builder = bind_band(builder, "recalibrated", &policy.recalibrated_band);
    builder = bind_policy(builder, "recalibrated", &policy.recalibrated_policy);
    builder
        .str("exchangeability-card", &policy.exchangeability_card_json)
        .finish()
}

fn first_campaign_mismatch(left: &Campaign, right: &Campaign) -> Option<String> {
    if left.regime_rotation != right.regime_rotation {
        return Some(format!(
            "regime-rotation:{}!={}",
            left.regime_rotation, right.regime_rotation
        ));
    }
    if left.selector_draws != right.selector_draws {
        return Some(format!(
            "selector-draws:{}!={}",
            left.selector_draws, right.selector_draws
        ));
    }
    if left.records.len() != right.records.len() {
        return Some(format!(
            "record-count:{}!={}",
            left.records.len(),
            right.records.len()
        ));
    }
    left.records
        .iter()
        .zip(&right.records)
        .enumerate()
        .find_map(|(index, (a, b))| {
            (a != b).then(|| format!("record[{index}]:left={a:?};right={b:?}"))
        })
}

fn campaign_accounting_mismatch(campaign: &Campaign) -> Option<String> {
    if campaign.regime_rotation >= REGIME_COUNT {
        return Some(format!("regime-rotation={}", campaign.regime_rotation));
    }
    if campaign.selector_draws != 1 {
        return Some(format!("selector-draws={}!=1", campaign.selector_draws));
    }
    if campaign.records.len() != REGIME_COUNT * ORBIT_SIZE {
        return Some(format!(
            "record-count={}!=expected-{}",
            campaign.records.len(),
            REGIME_COUNT * ORBIT_SIZE
        ));
    }

    let mut seen = [[0u8; ORBIT_SIZE]; REGIME_COUNT];
    for (index, record) in campaign.records.iter().enumerate() {
        let regime_position = index / ORBIT_SIZE;
        let rank_position = index % ORBIT_SIZE;
        let expected_regime = if campaign.regime_rotation == 0 {
            REGIMES[regime_position]
        } else {
            REGIMES[REGIME_COUNT - 1 - regime_position]
        };
        let expected_rank =
            (rank_position + campaign.regime_rotation * (ORBIT_SIZE / REGIME_COUNT)) % ORBIT_SIZE;
        if record.regime != expected_regime || record.rank != expected_rank {
            return Some(format!(
                "record[{index}]-schedule={}/{}!=expected-{}/{}",
                record.regime.name(),
                record.rank,
                expected_regime.name(),
                expected_rank
            ));
        }
        if record.rank >= ORBIT_SIZE {
            return Some(format!("record[{index}]-rank={}", record.rank));
        }
        let slot = &mut seen[record.regime.index()][record.rank];
        *slot = slot.saturating_add(1);
        if *slot != 1 {
            return Some(format!(
                "record[{index}]-duplicate={}/{}",
                record.regime.name(),
                record.rank
            ));
        }
        if record.residual_bits != record.regime.residual(record.rank).to_bits() {
            return Some(format!("record[{index}]-residual-mismatch"));
        }
        if !record.bucket_names_canonical {
            return Some(format!("record[{index}]-bucket-order"));
        }
        if !matches!(
            record.bucket_band,
            RetainedBand::Calibrated {
                n: MIN_CALIBRATION,
                half_width_bits
            } if f64::from_bits(half_width_bits).is_finite()
        ) {
            return Some(format!(
                "record[{index}]-bucket-band={:?}",
                record.bucket_band
            ));
        }
        if !matches!(
            record.marginal_band,
            RetainedBand::Calibrated {
                n,
                half_width_bits
            } if n == REGIME_COUNT * ORBIT_SIZE - 1
                && f64::from_bits(half_width_bits).is_finite()
        ) {
            return Some(format!(
                "record[{index}]-marginal-band={:?}",
                record.marginal_band
            ));
        }
        let Some(marginal_bits) = record.marginal_band.half_width().map(f64::to_bits) else {
            return Some(format!("record[{index}]-marginal-refused"));
        };
        if record.surrogate_half_width_bits != marginal_bits {
            return Some(format!(
                "record[{index}]-surrogate-bits=0x{:016x}!=marginal-0x{marginal_bits:016x}",
                record.surrogate_half_width_bits
            ));
        }
        if record.surrogate_alpha_bits != ALPHA.to_bits() {
            return Some(format!(
                "record[{index}]-surrogate-alpha=0x{:016x}!=0x{:016x}",
                record.surrogate_alpha_bits,
                ALPHA.to_bits()
            ));
        }
        let expected_residual = f64::from_bits(record.residual_bits);
        if record.bucket_covered != record.bucket_band.covers(expected_residual)
            || record.marginal_covered != record.marginal_band.covers(expected_residual)
        {
            return Some(format!("record[{index}]-stored-coverage-flag"));
        }
        let Some(bucket_bits) = record.bucket_band.half_width_bits() else {
            return Some(format!("record[{index}]-bucket-refused"));
        };
        if !record.bucket_policy.uses_exact_width(bucket_bits)
            || !record
                .pooled_policy
                .uses_exact_width(record.surrogate_half_width_bits)
        {
            return Some(format!("record[{index}]-unexpected-policy-width"));
        }
    }
    seen.iter().enumerate().find_map(|(regime, ranks)| {
        ranks
            .iter()
            .position(|&count| count != 1)
            .map(|rank| format!("seen[{regime}][{rank}]={}", ranks[rank]))
    })
}

fn first_accounting_mismatch(original: &Evidence, replay: &Evidence) -> Option<String> {
    campaign_accounting_mismatch(&original.campaign)
        .map(|mismatch| format!("original:{mismatch}"))
        .or_else(|| {
            campaign_accounting_mismatch(&replay.campaign)
                .map(|mismatch| format!("replay:{mismatch}"))
        })
        .or_else(|| (!original.policy.pass()).then(|| "original-policy-trace".to_string()))
        .or_else(|| (!replay.policy.pass()).then(|| "replay-policy-trace".to_string()))
}

#[derive(Debug, Clone, Copy)]
#[allow(clippy::struct_excessive_bools)]
struct Verdicts {
    conditional_coverage: bool,
    pooled_counterexample: bool,
    surrogate_seam: bool,
    drift_and_refresh: bool,
    accounting: bool,
    replay: bool,
}

impl Verdicts {
    const fn pass(self) -> bool {
        self.conditional_coverage
            && self.pooled_counterexample
            && self.surrogate_seam
            && self.drift_and_refresh
            && self.accounting
            && self.replay
    }
}

#[allow(clippy::too_many_arguments)]
fn result_identity(
    config: &ReplayIdentity,
    campaign: &ReplayIdentity,
    campaign_replay: &ReplayIdentity,
    policy: &ReplayIdentity,
    policy_replay: &ReplayIdentity,
    counts: CoverageCounts,
    evidence_mismatch: Option<&str>,
    accounting_mismatch: Option<&str>,
    verdicts: Verdicts,
) -> ReplayIdentity {
    IdentityBuilder::new("fs-eproc-conformal-shift-result-v1")
        .child("config", config)
        .child("campaign", campaign)
        .child("campaign-replay", campaign_replay)
        .child("policy", policy)
        .child("policy-replay", policy_replay)
        .u64("mondrian-smooth-hits", usize_u64(counts.mondrian_smooth))
        .u64("mondrian-shock-hits", usize_u64(counts.mondrian_shock))
        .u64("pooled-smooth-hits", usize_u64(counts.pooled_smooth))
        .u64("pooled-shock-hits", usize_u64(counts.pooled_shock))
        .u64("pooled-total-hits", usize_u64(counts.pooled_total()))
        .str(
            "first-evidence-mismatch",
            evidence_mismatch.unwrap_or("none"),
        )
        .str(
            "first-accounting-mismatch",
            accounting_mismatch.unwrap_or("none"),
        )
        .flag("conditional-coverage-pass", verdicts.conditional_coverage)
        .flag("pooled-counterexample-pass", verdicts.pooled_counterexample)
        .flag("surrogate-seam-pass", verdicts.surrogate_seam)
        .flag("drift-and-refresh-pass", verdicts.drift_and_refresh)
        .flag("accounting-pass", verdicts.accounting)
        .flag("replay-pass", verdicts.replay)
        .finish()
}

fn optional_json_string(value: Option<&str>) -> String {
    value.map_or_else(
        || "null".to_string(),
        |value| {
            let mut escaped = String::with_capacity(value.len() + 2);
            escaped.push('"');
            for character in value.chars() {
                match character {
                    '"' => escaped.push_str("\\\""),
                    '\\' => escaped.push_str("\\\\"),
                    '\n' => escaped.push_str("\\n"),
                    '\r' => escaped.push_str("\\r"),
                    '\t' => escaped.push_str("\\t"),
                    other => escaped.push(other),
                }
            }
            escaped.push('"');
            escaped
        },
    )
}

fn retained_band_json(band: &RetainedBand) -> String {
    match *band {
        RetainedBand::Calibrated { half_width_bits, n } => format!(
            "{{\"kind\":\"calibrated\",\"half_width_bits\":\"0x{half_width_bits:016x}\",\"n\":{n}}}"
        ),
        RetainedBand::Refused { have, need } => {
            format!("{{\"kind\":\"refused\",\"have\":{have},\"need\":{need}}}")
        }
    }
}

fn policy_json(policy: &PolicyRecord) -> String {
    match policy {
        PolicyRecord::UseSurrogate {
            band_half_width_bits,
        } => format!(
            "{{\"decision\":\"use-surrogate\",\"band_half_width_bits\":\"0x{band_half_width_bits:016x}\"}}"
        ),
        PolicyRecord::Escalate { reason } => {
            let reason_json = optional_json_string(Some(reason));
            format!("{{\"decision\":\"escalate\",\"reason\":{reason_json}}}")
        }
    }
}

fn coverage(count: usize) -> f64 {
    f64::from(u32::try_from(count).expect("coverage count fits u32"))
        / f64::from(u32::try_from(ORBIT_SIZE).expect("orbit size fits u32"))
}

fn emit_case(emitter: &mut Emitter, case: &str, pass: bool, detail: String) {
    let event = emitter.emit(
        if pass {
            Severity::Info
        } else {
            Severity::Error
        },
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass,
            detail,
            seed: INPUT_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("conformal-shift verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("conformal-shift verdict must use the fs-obs schema");
    println!("{line}");
}

fn emit_benchmark(emitter: &mut Emitter, metric: &str, value: f64) {
    let event = emitter.emit(
        Severity::Info,
        EventKind::BenchmarkResult {
            kernel: CASE.to_string(),
            metric: metric.to_string(),
            value,
            machine: 0,
        },
        None,
    );
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("conformal-shift row must use the fs-obs schema");
    println!("{line}");
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn emit_receipt(
    emitter: &mut Emitter,
    config: &ReplayIdentity,
    result: &ReplayIdentity,
    campaign: &ReplayIdentity,
    campaign_replay: &ReplayIdentity,
    policy: &ReplayIdentity,
    policy_replay: &ReplayIdentity,
    evidence: &Evidence,
    counts: CoverageCounts,
    evidence_mismatch: Option<&str>,
    accounting_mismatch: Option<&str>,
    verdicts: Verdicts,
) {
    let evidence_mismatch_json = optional_json_string(evidence_mismatch);
    let accounting_mismatch_json = optional_json_string(accounting_mismatch);
    let unseen_band_json = retained_band_json(&evidence.policy.unseen_band);
    let unseen_policy_json = policy_json(&evidence.policy.unseen_policy);
    let alpha_bits = ALPHA.to_bits();
    let conditional_coverage_pass = verdicts.conditional_coverage;
    let pooled_counterexample_pass = verdicts.pooled_counterexample;
    let surrogate_seam_pass = verdicts.surrogate_seam;
    let drift_and_refresh_pass = verdicts.drift_and_refresh;
    let accounting_pass = verdicts.accounting;
    let replay_pass = verdicts.replay;
    let event = emitter.emit(
        if verdicts.pass() {
            Severity::Info
        } else {
            Severity::Error
        },
        EventKind::Custom {
            name: "exact-conformal-shift-and-escalation".to_string(),
            json: format!(
                "{{\"config_identity\":\"{}\",\"result_identity\":\"{}\",\
                 \"campaign_identity\":\"{}\",\"campaign_replay_identity\":\"{}\",\
                 \"policy_identity\":\"{}\",\"policy_replay_identity\":\"{}\",\
                 \"units\":\"dimensionless-absolute-residual\",\
                 \"input_seed\":{INPUT_SEED},\"rotation_kernel\":{ROTATION_KERNEL},\
                 \"rotation_tile\":{ROTATION_TILE},\"regime_rotation\":{},\
                 \"selector_draws\":{},\"stream_semantics_version\":{},\
                 \"alpha\":{ALPHA},\"alpha_bits\":\"0x{alpha_bits:016x}\",\
                 \"orbit_size_per_regime\":{ORBIT_SIZE},\
                 \"minimum_calibration\":{MIN_CALIBRATION},\
                 \"coverage\":{{\"mondrian_smooth\":{},\"mondrian_shock\":{},\
                 \"pooled_smooth\":{},\"pooled_shock\":{},\"pooled_total\":{}}},\
                 \"counts\":{{\"mondrian_smooth\":{},\"mondrian_shock\":{},\
                 \"pooled_smooth\":{},\"pooled_shock\":{},\"pooled_total\":{}}},\
                 \"drift\":{{\"control_fired\":{},\"control_observations\":{},\
                 \"control_log_e_bits\":\"0x{:016x}\",\
                 \"shifted_detection\":{},\"observation_budget\":{DRIFT_BUDGET},\
                 \"validity_scale_bits\":\"0x{:016x}\",\
                 \"shifted_log_e_bits\":\"0x{:016x}\"}},\
                 \"unseen_bucket\":{{\"band\":{unseen_band_json},\
                 \"policy\":{unseen_policy_json}}},\
                 \"mutant\":\"pool-buckets-and-ignore-drift\",\
                 \"mutant_caught\":{},\"exchangeability_card\":{},\
                 \"verdicts\":{{\"conditional_coverage\":{conditional_coverage_pass},\
                 \"pooled_counterexample\":{pooled_counterexample_pass},\
                 \"surrogate_seam\":{surrogate_seam_pass},\
                 \"drift_and_refresh\":{drift_and_refresh_pass},\
                 \"accounting\":{accounting_pass},\"replay\":{replay_pass}}},\
                 \"first_evidence_mismatch\":{evidence_mismatch_json},\
                 \"first_accounting_mismatch\":{accounting_mismatch_json},\
                 \"versions\":{{\"fs_eproc\":\"{}\",\"fs_math\":\"{}\",\
                 \"fs_rand\":\"{}\",\"fs_obs\":\"{}\",\
                 \"fs_surrogate_api\":\"{FS_SURROGATE_API_SCHEMA}\"}},\
                 \"no_claims\":[\"arbitrary-distributions-or-ties\",\
                 \"all-alphas-or-sample-sizes\",\"all-bucket-schemes\",\
                 \"all-shift-magnitudes-or-projections\",\"weighted-conformal\",\
                 \"automatic-caller-wiring\",\"cross-ISA-execution\",\
                 \"Cx-or-cancellation\",\"performance\"],\"pass\":{}}}",
                config.hex(),
                result.hex(),
                campaign.hex(),
                campaign_replay.hex(),
                policy.hex(),
                policy_replay.hex(),
                evidence.campaign.regime_rotation,
                evidence.campaign.selector_draws,
                fs_rand::STREAM_SEMANTICS_VERSION,
                coverage(counts.mondrian_smooth),
                coverage(counts.mondrian_shock),
                coverage(counts.pooled_smooth),
                coverage(counts.pooled_shock),
                f64::from(u32::try_from(counts.pooled_total()).expect("count fits u32"))
                    / f64::from(u32::try_from(REGIME_COUNT * ORBIT_SIZE).expect("total fits u32")),
                counts.mondrian_smooth,
                counts.mondrian_shock,
                counts.pooled_smooth,
                counts.pooled_shock,
                counts.pooled_total(),
                evidence.policy.control_fired,
                evidence.policy.control_observations,
                evidence.policy.control_log_e_bits,
                evidence.policy.shifted_detection,
                evidence.policy.shifted_validity_scale_bits,
                evidence.policy.shifted_log_e_bits,
                verdicts.pooled_counterexample && verdicts.drift_and_refresh,
                evidence.policy.exchangeability_card_json,
                fs_eproc::VERSION,
                fs_math::VERSION,
                fs_rand::VERSION,
                fs_obs::VERSION,
                verdicts.pass(),
            ),
        },
        None,
    );
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("conformal-shift receipt must use the fs-obs schema");
    println!("{line}");
}

#[test]
#[allow(clippy::too_many_lines)]
fn exact_mondrian_shift_and_surrogate_escalation_replay() {
    let config = config_identity();
    let original = run_evidence();
    let replay = run_evidence();
    let campaign = campaign_identity(&config, &original.campaign);
    let campaign_replay = campaign_identity(&config, &replay.campaign);
    let policy = policy_identity(&config, &original.policy);
    let policy_replay = policy_identity(&config, &replay.policy);
    let counts = original.campaign.coverage_counts();
    let evidence_mismatch =
        first_campaign_mismatch(&original.campaign, &replay.campaign).or_else(|| {
            (original.policy != replay.policy).then(|| {
                format!(
                    "policy:left={:?};right={:?}",
                    original.policy, replay.policy
                )
            })
        });
    let accounting_mismatch = first_accounting_mismatch(&original, &replay);

    let conditional_coverage = counts.mondrian_smooth == EXPECTED_MONDRIAN_HITS
        && counts.mondrian_shock == EXPECTED_MONDRIAN_HITS;
    let pooled_counterexample = counts.pooled_smooth == EXPECTED_POOLED_SMOOTH_HITS
        && counts.pooled_shock == EXPECTED_POOLED_SHOCK_HITS
        && counts.pooled_total() == EXPECTED_POOLED_TOTAL_HITS
        && counts.pooled_total() == 2 * EXPECTED_MONDRIAN_HITS
        && counts.pooled_shock < EXPECTED_MONDRIAN_HITS;
    let surrogate_seam = original.campaign.records.iter().all(|record| {
        record.marginal_band.half_width_bits() == Some(record.surrogate_half_width_bits)
            && record.surrogate_alpha_bits == ALPHA.to_bits()
            && record
                .bucket_band
                .half_width_bits()
                .is_some_and(|bits| record.bucket_policy.uses_exact_width(bits))
            && record
                .pooled_policy
                .uses_exact_width(record.surrogate_half_width_bits)
    }) && original
        .policy
        .unseen_policy
        .escalates_with("conformal band is unbounded")
        && original
            .policy
            .known_band
            .half_width_bits()
            .is_some_and(|bits| original.policy.known_policy.uses_exact_width(bits));
    let verdicts = Verdicts {
        conditional_coverage,
        pooled_counterexample,
        surrogate_seam,
        drift_and_refresh: original.policy.pass(),
        accounting: accounting_mismatch.is_none(),
        replay: evidence_mismatch.is_none()
            && original == replay
            && campaign.root() == campaign_replay.root()
            && policy.root() == policy_replay.root(),
    };
    let result = result_identity(
        &config,
        &campaign,
        &campaign_replay,
        &policy,
        &policy_replay,
        counts,
        evidence_mismatch.as_deref(),
        accounting_mismatch.as_deref(),
        verdicts,
    );

    let mut emitter = Emitter::new(SUITE, CASE);
    emit_receipt(
        &mut emitter,
        &config,
        &result,
        &campaign,
        &campaign_replay,
        &policy,
        &policy_replay,
        &original,
        counts,
        evidence_mismatch.as_deref(),
        accounting_mismatch.as_deref(),
        verdicts,
    );
    emit_benchmark(
        &mut emitter,
        "mondrian_smooth_conditional_coverage",
        coverage(counts.mondrian_smooth),
    );
    emit_benchmark(
        &mut emitter,
        "mondrian_shock_conditional_coverage",
        coverage(counts.mondrian_shock),
    );
    emit_benchmark(
        &mut emitter,
        "pooled_shock_conditional_coverage",
        coverage(counts.pooled_shock),
    );
    emit_benchmark(
        &mut emitter,
        "pooled_marginal_coverage",
        f64::from(u32::try_from(counts.pooled_total()).expect("count fits u32"))
            / f64::from(u32::try_from(REGIME_COUNT * ORBIT_SIZE).expect("total fits u32")),
    );
    emit_benchmark(
        &mut emitter,
        "shifted_drift_detection_samples",
        f64::from(
            u32::try_from(original.policy.shifted_detection)
                .expect("drift detection budget fits u32"),
        ),
    );

    emit_case(
        &mut emitter,
        "exact-mondrian-conditional-coverage",
        verdicts.conditional_coverage,
        format!(
            "config={}; result={}; exhaustive leave-one-out hits smooth/shock={}/{} and {}/{}; nominal=18/20",
            config.hex(),
            result.hex(),
            counts.mondrian_smooth,
            ORBIT_SIZE,
            counts.mondrian_shock,
            ORBIT_SIZE,
        ),
    );
    emit_case(
        &mut emitter,
        "pooled-marginal-hides-shock-undercoverage",
        verdicts.pooled_counterexample,
        format!(
            "config={}; result={}; pooled hits smooth/shock/total={}/{}, {}/{}, {}/{}; overall 36/40 is nominal while shock 16/20 is not",
            config.hex(),
            result.hex(),
            counts.pooled_smooth,
            ORBIT_SIZE,
            counts.pooled_shock,
            ORBIT_SIZE,
            counts.pooled_total(),
            REGIME_COUNT * ORBIT_SIZE,
        ),
    );
    emit_case(
        &mut emitter,
        "fs-surrogate-width-and-policy-seam",
        verdicts.surrogate_seam,
        format!(
            "config={}; result={}; all 40 pooled widths are bit-equal across fs-eproc and fs-surrogate; calibrated finite bands authorize at their exact tolerance; unseen unbounded band escalates",
            config.hex(),
            result.hex(),
        ),
    );
    emit_case(
        &mut emitter,
        "drift-refusal-and-recalibration-policy",
        verdicts.drift_and_refresh,
        format!(
            "config={}; result={}; control fired={}; shifted detection={}/{}; correct drift wiring escalates, ignore-drift mutant authorizes, unseen bucket refuses, and 19 shock residuals recalibrate; card={}",
            config.hex(),
            result.hex(),
            original.policy.control_fired,
            original.policy.shifted_detection,
            DRIFT_BUDGET,
            original.policy.exchangeability_card_json,
        ),
    );
    emit_case(
        &mut emitter,
        "complete-accounting-and-bitwise-replay",
        verdicts.accounting && verdicts.replay,
        format!(
            "config={}; result={}; campaign roots={}/{}; policy roots={}/{}; first evidence mismatch={evidence_mismatch:?}; first accounting mismatch={accounting_mismatch:?}",
            config.hex(),
            result.hex(),
            campaign.hex(),
            campaign_replay.hex(),
            policy.hex(),
            policy_replay.hex(),
        ),
    );

    assert!(
        verdicts.pass(),
        "conformal-shift battery failed: verdicts={verdicts:?}; counts={counts:?}; evidence mismatch={evidence_mismatch:?}; accounting mismatch={accounting_mismatch:?}; policy={:?}",
        original.policy,
    );
}
