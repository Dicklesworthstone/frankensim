//! Measured budget lane (music bead `frankensim-music-v8-root-3ez8g.2.2`).
//!
//! D25: a live-default image needs a MEASURED budget row — samples/sec,
//! state count, headroom at 48 kHz, machine fingerprint, build profile,
//! repeatability band — never timing prose. This lane renders pinned
//! fixtures through the block API (so numbers are product-shaped, not
//! microbenchmark-shaped) and emits content-addressed budget rows the
//! claims registry references as `budget_row`.
//!
//! Why a test binary and not production code: wall-clock time must never
//! enter a production path or a receipt (host-stable bytes doctrine), and
//! the roofline precedent keeps measurement in dev lanes. The row bytes
//! carry the timing RESULTS as data with the build profile stamped, and
//! the ADMISSIBILITY rule is structural: `admissible_for_registry()` is
//! true only for `release` rows with dispersion inside the authored band —
//! debug numbers are logged but can never become registry evidence.
//!
//! The lane runs by default in this suite at a small rep count (it is a
//! schema + machinery gate); the `--ignored` heavy variant does the real
//! 32-rep measurement for row minting. Repeatability is measured as the
//! roofline convention's relative IQR dispersion ((p75-p25)/median) —
//! "a benchmark without variance bars is folklore." Host-quiescence
//! caveat (recorded per the perf-certification memory): tight bands on a
//! shared machine are uncertifiable; the band travels IN the row so a
//! reviewer sees the noise, and rows minted on a noisy host are honestly
//! wide rather than quietly lucky.

use fs_couple::reed_bore::ReedSolverMode;
use fs_couple::render::{ReedBoreVoice, RenderContext, RenderVoice};
use fs_couple::thin_plate::PlateBank;
use fs_duct::{Duct, Segment, Termination};
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::BeatingReed;
use fs_substrate::CapabilityProbe;

const RATE: u32 = 48_000;
const BLOCK: usize = 512;

/// One measured budget row: everything the registry's `budget_row`
/// reference must resolve to. Canonical line encoding mirrors the
/// bake-off receipt discipline (fixed order, `{:e}` floats, tabs).
struct BudgetRow {
    fixture: &'static str,
    image: &'static str,
    states: usize,
    block_len: usize,
    sample_rate_hz: u32,
    build_profile: &'static str,
    machine_fingerprint: u64,
    reps: usize,
    median_samples_per_sec: f64,
    dispersion: f64,
    headroom_at_48k: f64,
}

impl BudgetRow {
    fn canonical(&self) -> String {
        format!(
            "frankensim-budget-row-v1\nfixture\t{}\nimage\t{}\nstates\t{}\nblock\t{}\nrate\t{}\n\
             profile\t{}\nmachine\t{:016x}\nreps\t{}\nmedian-samples-per-sec\t{:e}\n\
             dispersion\t{:e}\nheadroom-48k\t{:e}\nadmissible\t{}\n",
            self.fixture,
            self.image,
            self.states,
            self.block_len,
            self.sample_rate_hz,
            self.build_profile,
            self.machine_fingerprint,
            self.reps,
            self.median_samples_per_sec,
            self.dispersion,
            self.headroom_at_48k,
            self.admissible_for_registry()
        )
    }

    /// D25 admissibility is structural: release build AND a dispersion
    /// inside the authored band. Debug rows are diagnostics forever.
    fn admissible_for_registry(&self) -> bool {
        self.build_profile == "release" && self.dispersion < 0.20
    }
}

fn build_profile() -> &'static str {
    if cfg!(debug_assertions) {
        "debug"
    } else {
        "release"
    }
}

fn reed_fixture(mode: ReedSolverMode) -> RenderContext {
    let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air");
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.0022,
            length: 0.50,
        }],
    };
    let reed = BeatingReed {
        rest_opening_m: 4.0e-4,
        width_m: 0.013,
        closing_pressure_pa: 6_000.0,
        blowing_pressure_pa: 2_800.0,
        attack_s: 0.008,
        mass_kg: 0.0,
        stiffness_n_m: 0.0,
    };
    let mut voice = ReedBoreVoice::new(
        &duct,
        &air,
        reed,
        Termination::UnflangedOpen,
        PlateBank::default(),
        1.0,
        RATE,
        RATE as usize,
        None,
    )
    .expect("voice admits");
    voice.set_solver_mode(mode);
    RenderContext::new(vec![RenderVoice::ReedBore(voice)], BLOCK)
}

/// Measure one fixture: warm-up, timed reps of `blocks_per_rep` blocks,
/// median + relative-IQR dispersion (the roofline convention).
fn measure(context: &mut RenderContext, reps: usize, blocks_per_rep: usize) -> (f64, f64) {
    let mut block = vec![0.0; BLOCK];
    // Warm-up: one rep untimed.
    for _ in 0..blocks_per_rep {
        context.block(&mut block).expect("warm-up");
    }
    let mut times: Vec<f64> = Vec::with_capacity(reps);
    for _ in 0..reps {
        let start = std::time::Instant::now();
        for _ in 0..blocks_per_rep {
            context.block(&mut block).expect("timed block");
        }
        times.push(start.elapsed().as_secs_f64());
    }
    times.sort_by(f64::total_cmp);
    let median = times[times.len() / 2];
    let p25 = times[times.len() / 4];
    let p75 = times[(3 * times.len()) / 4];
    let dispersion = if median > 0.0 {
        (p75 - p25) / median
    } else {
        f64::INFINITY
    };
    let samples = (blocks_per_rep * BLOCK) as f64;
    (samples / median, dispersion)
}

fn mint_row(reps: usize) -> BudgetRow {
    let mut context = reed_fixture(ReedSolverMode::Strict);
    let (samples_per_sec, dispersion) = measure(&mut context, reps, 32);
    BudgetRow {
        fixture: "massless-reed 2.2mm radiation-damped pipe, empty plate bank",
        image: "wind-reed/char-line",
        states: 1, // scalar reed state; the FIR line history is capacity, not dynamical state count
        block_len: BLOCK,
        sample_rate_hz: RATE,
        build_profile: build_profile(),
        machine_fingerprint: CapabilityProbe::run().fingerprint(),
        reps,
        median_samples_per_sec: samples_per_sec,
        dispersion,
        headroom_at_48k: samples_per_sec / f64::from(RATE),
    }
}

#[test]
fn budget_lane_measures_and_encodes_a_row() {
    // Schema/machinery gate at a small rep count: the row must encode,
    // carry a finite positive rate, a finite dispersion, the machine
    // fingerprint, and the correct admissibility semantics for THIS
    // build profile. The real minting run is the ignored heavy variant.
    let row = mint_row(5);
    assert!(row.median_samples_per_sec.is_finite() && row.median_samples_per_sec > 0.0);
    assert!(row.dispersion.is_finite() && row.dispersion >= 0.0);
    assert!(row.machine_fingerprint != 0, "fingerprint must be real");
    let bytes = row.canonical();
    println!("{bytes}");
    // Structural admissibility: a debug row can NEVER be admissible,
    // whatever its numbers — the rule is profile-first.
    if row.build_profile == "debug" {
        assert!(
            !row.admissible_for_registry(),
            "debug rows are diagnostics forever (D25)"
        );
    }
    // The canonical bytes carry every load-bearing field.
    for field in [
        "fixture\t",
        "image\t",
        "machine\t",
        "median-samples-per-sec\t",
        "dispersion\t",
        "headroom-48k\t",
        "admissible\t",
    ] {
        assert!(bytes.contains(field), "canonical row missing {field:?}");
    }
    // Repeatability of the MACHINERY (not a perf claim): a second
    // measurement on the same fixture class must produce the same
    // schema and a rate within a very loose sanity band (10x), so a
    // wildly broken timer or fixture is caught even in debug.
    let again = mint_row(5);
    let ratio = row.median_samples_per_sec / again.median_samples_per_sec;
    assert!(
        (0.1..=10.0).contains(&ratio),
        "same-fixture rates differ by {ratio}x; the lane itself is broken"
    );
}

#[test]
#[ignore = "heavy minting run: 32 reps in release via the recorded heavy-run recipe"]
fn mint_registry_budget_row() {
    // The real row: 32 reps. Run in release (RCH_MIN_LOCAL_TIME_MS
    // recipe, --release); the printed canonical bytes are the artifact to
    // commit under data/budget-rows/ and reference from the registry's
    // budget_row field once a live-default promotion is reviewed.
    let row = mint_row(32);
    println!("{}", row.canonical());
    assert!(
        row.build_profile == "release",
        "minting runs must be release builds; debug numbers are never admissible"
    );
    assert!(
        row.headroom_at_48k > 1.0,
        "fixture renders below real time on this machine ({}x); a live-default \
         claim would be false here",
        row.headroom_at_48k
    );
}

// --- Fusion #1 before/after rows (bead frankensim-kigr6) ----------------
//
// The aperture-solve fusion (2s4i5) already flipped its pricing; this
// lane's job for kigr6 is the BEFORE row pair for the FIR convolution
// fusion candidate: the SAME massless-reed fixture rendered through
// both aperture solver modes so the char-image budget delta from the
// junction solve is isolated from everything else. Strict is the
// certification default; FastNewton rows are the declared-fast-mode
// comparator.

fn mint_fusion_pair(reps: usize) -> (BudgetRow, BudgetRow, f64) {
    let mut strict_ctx = reed_fixture(ReedSolverMode::Strict);
    let (strict_rate, strict_disp) = measure(&mut strict_ctx, reps, 32);
    let mut fast_ctx = reed_fixture(ReedSolverMode::FastNewton);
    let (fast_rate, fast_disp) = measure(&mut fast_ctx, reps, 32);
    let base = BudgetRow {
        fixture: "massless-reed 2.2mm radiation-damped pipe, empty plate bank",
        image: "wind-reed/char-line",
        states: 1,
        block_len: BLOCK,
        sample_rate_hz: RATE,
        build_profile: build_profile(),
        machine_fingerprint: CapabilityProbe::run().fingerprint(),
        reps,
        median_samples_per_sec: strict_rate,
        dispersion: strict_disp,
        headroom_at_48k: strict_rate / f64::from(RATE),
    };
    let fused = BudgetRow {
        fixture: base.fixture,
        image: "wind-reed/char-line+fast-newton-junction",
        states: base.states,
        block_len: base.block_len,
        sample_rate_hz: base.sample_rate_hz,
        build_profile: base.build_profile,
        machine_fingerprint: base.machine_fingerprint,
        reps: base.reps,
        median_samples_per_sec: fast_rate,
        dispersion: fast_disp,
        headroom_at_48k: fast_rate / f64::from(RATE),
    };
    let ratio = fast_rate / strict_rate.max(f64::MIN_POSITIVE);
    (base, fused, ratio)
}

fn fusion_row_canonical(label: &str, row: &BudgetRow, ratio: Option<f64>) -> String {
    let mut s = format!(
        "frankensim-fusion-budget-row-v2\nlane\t{klabel}\n{rcanonical}",
        klabel = label,
        rcanonical = {
            let _ = label;
            row.canonical()
        }
    )
    .replace("frankensim-budget-row-v1", "");
    if let Some(r) = ratio {
        s.push_str(&format!("ratio-vs-strict\t{r:e}\n"));
    }
    s
}

#[test]
fn fusion_before_after_rows_encode_and_stay_sane() {
    // Machinery gate at small reps (mirrors budget_lane_measures...):
    // both modes must encode finite rows and the FAST mode must not be
    // catastrophically slower than STRICT on the same host (10x sanity
    // band; a real regression shows up as a budget-row finding, and a
    // genuine speedup shows up in the release mint).
    let (strict_row, fast_row, ratio) = mint_fusion_pair(5);
    for row in [&strict_row, &fast_row] {
        assert!(row.median_samples_per_sec.is_finite() && row.median_samples_per_sec > 0.0);
        assert!(row.dispersion.is_finite());
        let bytes = row.canonical();
        assert!(bytes.contains("image\twind-reed/char-line"));
    }
    println!(
        "{}",
        fusion_row_canonical("strict-before", &strict_row, None)
    );
    println!(
        "{}",
        fusion_row_canonical("fast-newton-after", &fast_row, Some(ratio))
    );
    assert!(
        (0.1..=10.0).contains(&ratio),
        "fast-newton junction rate is {ratio}x strict on the same host; \\
         outside the 10x sanity band either the lane or the solver is broken"
    );
}

#[test]
#[ignore = "heavy minting run: 32 reps in release via the recorded heavy-run recipe"]
fn mint_fusion_registry_rows() {
    let (strict_row, fast_row, ratio) = mint_fusion_pair(32);
    assert!(strict_row.build_profile == "release");
    println!(
        "{}",
        fusion_row_canonical("strict-before", &strict_row, None)
    );
    println!(
        "{}",
        fusion_row_canonical("fast-newton-after", &fast_row, Some(ratio))
    );
}
