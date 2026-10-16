//! First fs-modalid consumer: identify the existing coupled solver's FRF.
//! Authored modal fixture, not a measured bar/body specimen or a 3-D claim.
//! The two cavity modes have nonzero overlap with distinct structural modes;
//! the remaining four structural modes are orthogonal to that interface.

use fs_couple::vibroacoustic::{
    AcousticMedium, CavityModes, StructuralModes, VibroacousticModel, assemble_coupling,
};
use fs_math::{c64::C64, det};
use fs_modalid::{
    FrfData, Identification, IdentifyOptions, ModalIdError, ModeCompareOptions, StabilizationRule,
    compare_mode_tables, estimate_snr, identify,
};
use fs_rand::StreamKey;
use std::sync::OnceLock;
use std::time::Instant;

const TAU: f64 = core::f64::consts::TAU;
const STRUCTURAL_HZ: [f64; 6] = [120.0, 260.0, 410.0, 570.0, 760.0, 980.0];
const CAVITY_HZ: [f64; 2] = [180.0, 340.0];
const COUPLING_RATIO: [f64; 2] = [0.35, 0.30];
const ETA_S: f64 = 0.01;
const ETA_A: f64 = 0.005;
const GRID_POINTS: usize = 4096;
const NOISE_SEED: u64 = 0x4d37_3457_3151_5001;
const NOISY_TEST: &str = "modalid_qp::noisy_frf_inside_intervals_and_snr_refusal";

struct Fixture {
    structure: StructuralModes,
    cavity: CavityModes,
    coupling: Vec<f64>,
    model: VibroacousticModel,
}

fn fixture(shift_cents: f64) -> Fixture {
    let mut omegas: Vec<f64> = STRUCTURAL_HZ.iter().map(|f| TAU * f).collect();
    // The sixth mode is interface-orthogonal. Its deliberate shift must not
    // move the other seven modes through a hidden global parameter change.
    omegas[5] *= det::exp(shift_cents * core::f64::consts::LN_2 / 1200.0);
    let structure = StructuralModes {
        omegas,
        shapes: (0..6)
            .map(|r| (0..6).map(|i| if r == i { 1.0 } else { 0.0 }).collect())
            .collect(),
        loss_factor: ETA_S,
    };
    let air = AcousticMedium::air();
    let omega_a: Vec<f64> = CAVITY_HZ.iter().map(|f| TAU * f).collect();
    // Pressure basis amplitudes are arbitrary: scale psi and its norm together
    // so s=rho*c²/Lambda=1/omega_a². This balances the existing unbalanced dense
    // pencil without changing physical poles. Lambda is a mode norm, not volume.
    let cavity = CavityModes {
        lambdas: omega_a
            .iter()
            .map(|w| air.rho0 * air.c0 * air.c0 * w * w)
            .collect(),
        interface: (0..2)
            .map(|q| {
                (0..6)
                    .map(|i| {
                        if i == q {
                            COUPLING_RATIO[q] * omega_a[q] * omega_a[q] / 0.01
                        } else {
                            0.0
                        }
                    })
                    .collect()
            })
            .collect(),
        omegas: omega_a,
        loss_factor: ETA_A,
        rho0: air.rho0,
        c0: air.c0,
    };
    let coupling = assemble_coupling(&structure, &cavity, &[0.01; 6]).unwrap();
    let model = VibroacousticModel::try_new(&structure, &cavity, coupling.clone(), None).unwrap();
    Fixture {
        structure,
        cavity,
        coupling,
        model,
    }
}

fn options(snr_floor: f64) -> IdentifyOptions {
    // Eight true conjugate pairs require order16, plus two higher orders for
    // the declared three-consecutive-order stabilization gate. No defaults.
    IdentifyOptions {
        min_order: 16,
        max_order: 20,
        order_step: 2,
        snr_floor,
        rule: StabilizationRule {
            freq_tol: 0.001,
            damp_tol: 0.05,
            required_runs: 3,
        },
        window_tau: None,
    }
}

fn frf_data(f: &Fixture, noise: f64) -> FrfData {
    let omega: Vec<f64> = (0..GRID_POINTS)
        .map(|i| TAU * (80.0 + 1120.0 * i as f64 / (GRID_POINTS - 1) as f64))
        .collect();
    let mut channels = vec![Vec::with_capacity(GRID_POINTS); 9];
    let force = [C64::ONE; 6];
    for &w in &omega {
        let response = f.model.frf(w, &force).unwrap();
        // The reference observes every mode; individual structural/pressure
        // modal coordinates supply distinct shapes on these same nine channels.
        channels[0].push(response.b.iter().copied().fold(C64::ZERO, |a, b| a + b));
        for (r, &value) in response.b.iter().enumerate() {
            channels[1 + r].push(value);
        }
        for (q, &value) in response.a.iter().enumerate() {
            channels[7 + q].push(value);
        }
    }
    let kernel = NOISY_TEST.bytes().fold(0x811c_9dc5_u32, |h, b| {
        (h ^ u32::from(b)).wrapping_mul(0x0100_0193)
    });
    for (channel, h) in channels.iter_mut().enumerate() {
        let mut stream = StreamKey {
            seed: NOISE_SEED,
            kernel,
            tile: u32::try_from(channel).unwrap(),
        }
        .stream();
        for value in h {
            let phase = TAU * stream.next_f64();
            let magnitude = noise * value.abs();
            *value = *value + C64::new(det::cos(phase), det::sin(phase)).scale(magnitude);
        }
    }
    // The solver is engineering e^-iwt; ingest must conjugate to fs-vfit's axis.
    FrfData::new(omega, channels.into_iter().map(|h| (h, None)).collect()).unwrap()
}

fn clean_identification() -> &'static Identification {
    static CLEAN: OnceLock<Identification> = OnceLock::new();
    CLEAN.get_or_init(|| identify(&frf_data(&fixture(0.0), 0.0), &options(10.0)).unwrap())
}

/// Independent complex eigenvalue oracle for this block-coupled fixture.
/// For x=omega², each two-mode block has
/// x²-(ds+da+s C²)x+ds*da=0, ds=omega_s²(1-i eta_s), likewise da.
/// Four remaining roots are ds alone. This is not eta/2 or an FRF peak pick.
fn coupled_poles(f: &Fixture) -> Vec<(f64, f64)> {
    let mut squares = Vec::new();
    for q in 0..2 {
        let ds = C64::new(1.0, -ETA_S).scale(f.structure.omegas[q] * f.structure.omegas[q]);
        let da = C64::new(1.0, -ETA_A).scale(f.cavity.omegas[q] * f.cavity.omegas[q]);
        let s = f.cavity.rho0 * f.cavity.c0 * f.cavity.c0 / f.cavity.lambdas[q];
        let c = f.coupling[q * 2 + q];
        let trace = ds + da + C64::from_re(s * c * c);
        let discriminant = (trace * trace - (ds * da).scale(4.0)).sqrt();
        squares.extend([
            (trace - discriminant).scale(0.5),
            (trace + discriminant).scale(0.5),
        ]);
    }
    for &w in &f.structure.omegas[2..] {
        squares.push(C64::new(w * w, -ETA_S * w * w));
    }
    let mut poles: Vec<_> = squares
        .into_iter()
        .map(|x| {
            let w = x.sqrt();
            assert!(w.re > 0.0 && w.im < 0.0);
            (w.abs() / TAU, -w.im / w.abs())
        })
        .collect();
    poles.sort_by(|a, b| a.0.total_cmp(&b.0));
    poles
}

fn cents(candidate: f64, reference: f64) -> f64 {
    1200.0 * det::ln(candidate / reference) / core::f64::consts::LN_2
}

/// G1: recover every actual coupled mode, its damping, and no extra mode.
#[test]
fn modal_id_recovers_simulated_coupled_modes() {
    let started = Instant::now();
    let f = fixture(0.0);
    let undamped = f.model.undamped_natural_frequencies().unwrap();
    let poles = coupled_poles(&f);
    let id = clean_identification();
    assert_eq!(undamped.len(), 8);
    assert_eq!(
        id.modes.len(),
        undamped.len(),
        "no missing or spurious accepted mode: {:?}",
        id.modes
    );
    assert!(
        (undamped[0] / f.structure.omegas[0] - 1.0).abs() > 0.01,
        "coupling must move a pole"
    );
    for (rung, ((mode, &w), &(pole_hz, zeta))) in
        id.modes.iter().zip(&undamped).zip(&poles).enumerate()
    {
        let df = (mode.frequency_hz / (w / TAU) - 1.0).abs();
        let dz = (mode.damping_ratio / zeta - 1.0).abs();
        assert!(df < 0.0005, "mode{rung}: {df:e}, {mode:?}");
        assert!(
            dz < 0.05,
            "mode{rung}: damping error{dz:e}, coupled-pole zeta{zeta:e}"
        );
        // Keep the complex-pole frequency visible alongside the conservative one.
        assert!((mode.frequency_hz / pole_hz - 1.0).abs() < 0.0005);
        println!(
            "{{\"case\":\"modal_id_recovers_simulated_coupled_modes\",\"inputs\":\"analytic six-structure/two-cavity model; eta_s=.01 eta_a=.005; 4096 FRF samples\",\"axis\":\"coupled mode\",\"rung\":{rung},\"observed\":[{df:e},{dz:e}],\"oracle\":\"undamped_natural_frequencies and complex quadratic eigenvalues\",\"tol\":[0.0005,0.05],\"verdict\":\"pass\",\"elapsed_ms\":{}}}",
            started.elapsed().as_millis()
        );
    }
}

/// G3: compare identified QP mode tables and retain the deliberately bad pair.
#[test]
fn qp_mode_table_comparison_flags_a_perturbed_model() {
    let started = Instant::now();
    let reference = clean_identification();
    let opts = ModeCompareOptions {
        mac_floor: 0.95,
        max_cents: 1.0,
    };
    let same = compare_mode_tables(&reference.modes, &reference.modes, &opts).unwrap();
    assert_eq!(same.pairs.len(), 8);
    assert!(same.unmatched_reference.is_empty() && same.unmatched_candidate.is_empty());
    assert!(
        same.pairs
            .iter()
            .all(|p| p.within_tolerance && p.cents.abs() <= 1.0 && p.mac.unwrap() >= 0.95)
    );
    let changed = identify(&frf_data(&fixture(20.0), 0.0), &options(10.0)).unwrap();
    let diff = compare_mode_tables(&reference.modes, &changed.modes, &opts).unwrap();
    assert_eq!(diff.pairs.len(), 8);
    assert!(diff.unmatched_reference.is_empty() && diff.unmatched_candidate.is_empty());
    let rejected: Vec<_> = diff.pairs.iter().filter(|p| !p.within_tolerance).collect();
    assert_eq!(
        rejected.len(),
        1,
        "single structural perturbation: {diff:?}"
    );
    assert!((rejected[0].cents - 20.0).abs() <= 1.0);
    for (rung, pair) in diff.pairs.iter().enumerate() {
        assert!(pair.mac.unwrap() >= 0.95);
        if pair.within_tolerance {
            assert!(pair.cents.abs() <= 1.0);
        }
        println!(
            "{{\"case\":\"qp_mode_table_comparison_flags_a_perturbed_model\",\"inputs\":\"analytic coupled FRF; sixth structural mode +20c\",\"axis\":\"mode pair\",\"rung\":{rung},\"observed\":{},\"oracle\":\"zero shift except one +20c pole\",\"mac\":{},\"within_tolerance\":{},\"tol\":1.0,\"verdict\":\"pass\",\"elapsed_ms\":{}}}",
            pair.cents,
            pair.mac.unwrap(),
            pair.within_tolerance,
            started.elapsed().as_millis()
        );
    }
}

/// G2: seeded 1% relative complex noise and the exact estimated-SNR boundary.
#[test]
fn noisy_frf_inside_intervals_and_snr_refusal() {
    let started = Instant::now();
    let f = fixture(0.0);
    let data = frf_data(&f, 0.01);
    let snr = estimate_snr(&data);
    assert!(snr > 0.0 && snr.is_finite());
    let above = identify(&data, &options(snr.next_up()));
    assert!(
        matches!(above, Err(ModalIdError::SnrTooLow { snr: got, floor }) if got == snr && floor > snr)
    );
    let id = identify(&data, &options(snr.next_down())).unwrap();
    let truth = f.model.undamped_natural_frequencies().unwrap();
    assert_eq!(
        id.modes.len(),
        truth.len(),
        "no noise-induced missing/spurious mode"
    );
    for (rung, (mode, &w)) in id.modes.iter().zip(&truth).enumerate() {
        let target = w / TAU;
        let delta = (mode.frequency_hz - target).abs();
        let shift = cents(mode.frequency_hz, target).abs();
        let inside = mode.frequency_ci_hz.is_finite() && delta <= mode.frequency_ci_hz;
        assert!(
            inside || shift <= 3.0,
            "mode{rung}: {mode:?}, truth{target}, cents{shift}"
        );
        println!(
            "{{\"case\":\"noisy_frf_inside_intervals_and_snr_refusal\",\"inputs\":\"analytic coupled FRF; relative complex noise .01; Philox seed {NOISE_SEED}; test-name stream\",\"axis\":\"coupled mode\",\"rung\":{rung},\"observed\":{shift:e},\"frequency_ci_hz\":{},\"oracle\":\"undamped coupled frequencies\",\"snr\":{snr:e},\"tol\":\"inside split-sample interval or <=3c; SNR next_up refuses, next_down identifies\",\"verdict\":\"pass\",\"elapsed_ms\":{}}}",
            if mode.frequency_ci_hz.is_finite() {
                mode.frequency_ci_hz.to_string()
            } else {
                "null".to_owned()
            },
            started.elapsed().as_millis()
        );
    }
}
