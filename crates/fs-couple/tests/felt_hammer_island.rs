//! Wool-felt hammer–string strike island (music bead
//! `frankensim-music-t-piano-felt-87zbd`): the new `fs_material::WoolFelt`
//! law driving a modal string through an explicit contact island — the
//! piano-attack MECHANISM demonstrated at fixture scale.
//!
//! Composition: a hammer mass flies at the string; while the felt overlap
//! is positive, the contact force is `A_c · σ_felt(overlap / t_felt)`
//! through the committed hysteretic state (loading rides the stiffening
//! envelope; unloading rides the steeper crush-anchored path, so the
//! strike dissipates by construction). The string is three exact-ZOH
//! modal coordinates (`fs_couple::modal_acoustic_time`); the island
//! substeps the hammer–felt ODE inside each audio sample (the plan's
//! "substepped hysteretic island").
//!
//! Physics gates, each the felt's OWN signature rather than generic
//! oscillation:
//! 1. CONTACT DURATION falls with impact velocity (the stiffening p > 1
//!    envelope: harder blows meet a stiffer felt) — the classic measured
//!    hammer behavior, and the OPPOSITE of a linear spring whose contact
//!    time is amplitude-independent.
//! 2. SPECTRAL TILT: the harder strike puts proportionally MORE energy
//!    into the upper string modes than the soft strike (mode-3 to mode-1
//!    energy ratio rises with velocity), while the LINEAR-SPRING control
//!    on the same fixture shows a much smaller shift — the hysteretic
//!    stiffening law is thereby audible in the spectrum, not asserted.
//! 3. ENERGY HONESTY: the felt's loading/unloading loop dissipates a
//!    strictly positive fraction of the hammer's kinetic energy, and the
//!    total ledger (string energy + returned hammer KE + felt loss)
//!    closes against the initial KE within a stated band.
//!
//! Provenance honesty: the felt parameters here are AUTHORED,
//! Estimate-labeled fixture values (the licensed coupon hunt found the
//! Stulov/Chabassier sources but PDF retention was blocked — recorded on
//! the corpus absence row `acoustic-absent-felt-coupon`); nothing in this
//! test claims measured-felt authority. The mechanism gates are exactly
//! the ones a future coupon-backed card must also pass.

use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_material::{Uniaxial, WoolFelt};

const RATE: u32 = 48_000;
/// Island substeps per audio sample (contact is stiff and brief).
const SUBSTEPS: usize = 64;
/// Felt pad thickness [m] — overlap/thickness is the felt strain.
const T_FELT: f64 = 8.0e-3;
/// Contact patch area [m^2].
const A_CONTACT: f64 = 1.0e-4;
/// Hammer mass [kg] (piano-hammer scale).
const M_HAMMER: f64 = 8.0e-3;

fn string_model() -> ModalAcousticTimeModel {
    let wave_speed = (60.0f64 / 6.0e-4).sqrt();
    let modes = (1..=3)
        .map(|k| ModalAcousticMode {
            angular_frequency_rad_s: k as f64 * core::f64::consts::PI * wave_speed / 0.65,
            damping_ratio: 5.0e-4 * k as f64,
            pressure_per_modal_velocity: fs_math::c64::C64::new(1.0, 0.0),
        })
        .collect();
    ModalAcousticTimeModel::try_new(RATE, modes, ModalAcousticTimeBudget::audible_reference())
        .expect("string modes admit")
}

/// Contact-force law for the strike: felt (hysteretic) or the
/// linear-spring control matched to the felt's SECANT stiffness at a
/// reference overlap, so the two fixtures are comparable.
enum Hammer {
    Felt {
        law: WoolFelt,
        state: <WoolFelt as Uniaxial>::State,
    },
    LinearSpring {
        k: f64,
    },
}

impl Hammer {
    fn force(&self, overlap: f64) -> f64 {
        if overlap <= 0.0 {
            return 0.0;
        }
        match self {
            Self::Felt { law, state } => {
                let strain = (overlap / T_FELT).min(0.65);
                A_CONTACT * law.stress(strain, state)
            }
            Self::LinearSpring { k } => k * overlap,
        }
    }

    fn commit(&mut self, overlap: f64) {
        if let Self::Felt { law, state } = self {
            let strain = (overlap.max(0.0) / T_FELT).min(0.65);
            *state = law.update_state(strain, state);
        }
    }
}

struct StrikeOutcome {
    /// Samples with positive contact force.
    contact_samples: usize,
    /// Modal energies after the hammer separates [J].
    modal_energy_j: Vec<f64>,
    /// Hammer kinetic energy after separation [J].
    hammer_ke_out_j: f64,
    /// Energy dissipated in the felt loop [J] (ledger residual).
    felt_loss_j: f64,
}

/// Strike the string at `v0` and integrate until separation (+ margin).
fn strike(mut hammer: Hammer, v0: f64) -> StrikeOutcome {
    let mut string = string_model();
    // Mode shapes at the strike point (normalized sin(k pi x0/L), x0/L = 0.12
    // — near-bridge-ish striking point).
    let phi: Vec<f64> = (1..=3)
        .map(|k| fs_math::det::sin(k as f64 * core::f64::consts::PI * 0.12))
        .collect();
    let dt = 1.0 / f64::from(RATE);
    let h = dt / SUBSTEPS as f64;
    // Hammer starts just touching, flying at the string.
    let mut y_hammer = 0.0f64; // hammer tip minus string rest line
    let mut v_hammer = v0;
    let ke_in = 0.5 * M_HAMMER * v0 * v0;
    let mut work_into_string = 0.0f64;
    let mut contact_samples = 0usize;
    let mut frames = 0usize;
    loop {
        // String displacement at the strike point from modal states.
        let string_disp: f64 = string
            .states()
            .iter()
            .zip(&phi)
            .map(|(s, p)| s.displacement_m_sqrt_kg * p)
            .sum();
        // Substep the hammer against a frozen string displacement (the
        // string moves at audio rate; the island resolves the stiff felt).
        let mut force_avg = 0.0f64;
        let mut peak_overlap = 0.0f64;
        for _ in 0..SUBSTEPS {
            let overlap = y_hammer - string_disp;
            let f = hammer.force(overlap);
            force_avg += f;
            peak_overlap = peak_overlap.max(overlap);
            // Semi-implicit Euler on the hammer.
            v_hammer -= h * f / M_HAMMER;
            y_hammer += h * v_hammer;
        }
        force_avg /= SUBSTEPS as f64;
        hammer.commit(peak_overlap.max(y_hammer - string_disp));
        if force_avg > 0.0 {
            contact_samples += 1;
        }
        // Drive the string with the sample-held modal generalized forces.
        let generalized: Vec<f64> = phi.iter().map(|p| force_avg * p).collect();
        let frame = string.step(&generalized).expect("string step");
        work_into_string += frame.input_work_j;
        frames += 1;
        // Separation: hammer moving away and no contact for a margin.
        let separated = v_hammer < 0.0 && (y_hammer - string_disp) < -1.0e-4;
        if separated || frames > 4 * RATE as usize / 100 {
            let modal_energy_j = frame.modal_energy_j.clone();
            let hammer_ke_out_j = 0.5 * M_HAMMER * v_hammer * v_hammer;
            // Felt loss = KE_in - KE_out - work delivered to the string
            // (the hammer's whole energy budget; string damping loss during
            // contact is inside the string's own ledger, not the felt's).
            let felt_loss_j = ke_in - hammer_ke_out_j - work_into_string;
            return StrikeOutcome {
                contact_samples,
                modal_energy_j,
                hammer_ke_out_j,
                felt_loss_j,
            };
        }
    }
}

fn felt() -> WoolFelt {
    // AUTHORED fixture parameters (Estimate-labeled; see module doc):
    // exponent p = 2.5 in the published hammer range, q = 3.2, 25% crush.
    WoolFelt::new(4.0e5, 0.2, 2.5, 3.2, 0.25, 0.8).expect("felt admits")
}

/// Linear control matched to the felt's secant stiffness at 30% strain.
fn matched_spring() -> Hammer {
    let law = felt();
    let s = law.envelope(0.3).0;
    Hammer::LinearSpring {
        k: A_CONTACT * s / (0.3 * T_FELT),
    }
}

#[test]
fn contact_duration_falls_with_velocity_unlike_the_linear_spring() {
    let soft = strike(
        Hammer::Felt {
            law: felt(),
            state: felt().initial_state(),
        },
        0.8,
    );
    let hard = strike(
        Hammer::Felt {
            law: felt(),
            state: felt().initial_state(),
        },
        4.0,
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"felt-contact-duration\",\"soft_samples\":{},\
         \"hard_samples\":{}}}",
        soft.contact_samples, hard.contact_samples
    );
    assert!(
        hard.contact_samples < soft.contact_samples,
        "stiffening felt must shorten contact with velocity ({} vs {})",
        hard.contact_samples,
        soft.contact_samples
    );
    // The linear-spring control: contact time is amplitude-independent
    // (half-period of the spring-mass pair) — within one sample.
    let spring_soft = strike(matched_spring(), 0.8);
    let spring_hard = strike(matched_spring(), 4.0);
    assert!(
        (spring_soft.contact_samples as i64 - spring_hard.contact_samples as i64).abs() <= 1,
        "linear control must hold contact time ({} vs {})",
        spring_soft.contact_samples,
        spring_hard.contact_samples
    );
}

#[test]
fn spectral_tilt_rises_with_velocity_beyond_the_linear_control() {
    let tilt = |o: &StrikeOutcome| o.modal_energy_j[2] / o.modal_energy_j[0].max(1.0e-30);
    let felt_soft = strike(
        Hammer::Felt {
            law: felt(),
            state: felt().initial_state(),
        },
        0.8,
    );
    let felt_hard = strike(
        Hammer::Felt {
            law: felt(),
            state: felt().initial_state(),
        },
        4.0,
    );
    let spring_soft = strike(matched_spring(), 0.8);
    let spring_hard = strike(matched_spring(), 4.0);
    let felt_shift = tilt(&felt_hard) / tilt(&felt_soft).max(1.0e-30);
    let spring_shift = tilt(&spring_hard) / tilt(&spring_soft).max(1.0e-30);
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"felt-spectral-tilt\",\"felt_shift\":{felt_shift:.4},\
         \"spring_shift\":{spring_shift:.4}}}"
    );
    assert!(
        felt_shift > 1.2,
        "harder felt strikes must brighten the spectrum (shift {felt_shift})"
    );
    assert!(
        felt_shift > 1.5 * spring_shift,
        "the tilt must come from the FELT nonlinearity, not the fixture: felt shift \
         {felt_shift} vs linear control {spring_shift}"
    );
}

#[test]
fn the_felt_loop_dissipates_and_the_ledger_closes() {
    let v0 = 2.5;
    let out = strike(
        Hammer::Felt {
            law: felt(),
            state: felt().initial_state(),
        },
        v0,
    );
    let ke_in = 0.5 * M_HAMMER * v0 * v0;
    let string_total: f64 = out.modal_energy_j.iter().sum();
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"felt-energy-ledger\",\"ke_in\":{ke_in:.6},\
         \"ke_out\":{:.6},\"string\":{string_total:.6},\"felt_loss\":{:.6}}}",
        out.hammer_ke_out_j, out.felt_loss_j
    );
    assert!(
        out.felt_loss_j > 0.02 * ke_in,
        "the hysteretic loop must dissipate a real fraction of the strike \
         (loss {} of {ke_in})",
        out.felt_loss_j
    );
    // Ledger closure: everything is accounted for. The band covers the
    // frozen-string-during-substep splitting error of the island.
    let accounted = out.hammer_ke_out_j + string_total + out.felt_loss_j;
    assert!(
        (accounted - ke_in).abs() / ke_in < 0.15,
        "energy ledger must close within the island's splitting band \
         ({accounted} vs {ke_in})"
    );
    // The linear control returns essentially everything it borrowed
    // (its only loss channel is string damping during contact).
    let spring = strike(matched_spring(), v0);
    assert!(
        spring.felt_loss_j < 0.5 * out.felt_loss_j,
        "the elastic control must dissipate far less than the felt \
         ({} vs {})",
        spring.felt_loss_j,
        out.felt_loss_j
    );
}

#[test]
fn prony_chain_covers_the_hammer_band_and_refuses_outside() {
    // DONE-WHEN 3 of the felt bead: the rate-dependent half of felt runs
    // through the EXISTING visco chain — an authored (Estimate-labeled)
    // FractionalZener lowered via lower_to_prony to a GeneralizedMaxwell
    // whose certified band covers the hammer contact spectrum
    // (~50 Hz .. 5 kHz), with the out-of-band refusal firing by name.
    use fs_material::visco::{FractionalZener, ViscoError, lower_to_prony};
    let fz = FractionalZener::new(4.0e5, 1.2e6, 0.35, 1.0e-3).expect("felt-ish FZ admits");
    let lowered = lower_to_prony(&fz, 50.0, 5_000.0, 12, 0.05).expect("lowering certifies");
    let (lo, hi) = lowered.band;
    let tau = core::f64::consts::TAU;
    assert!(
        lo <= tau * 50.0 && hi >= tau * 5_000.0,
        "certified band [{lo}, {hi}] rad/s must cover the hammer spectrum"
    );
    assert!(
        lowered.sup_rel_err <= 0.05,
        "certificate honest: sup error {}",
        lowered.sup_rel_err
    );
    // In-band evaluation works; out-of-band refuses BY NAME.
    lowered
        .modulus_checked(tau * 500.0)
        .expect("mid-band evaluates");
    match lowered.modulus_checked(tau * 50_000.0) {
        Err(ViscoError::OutOfBand { .. }) => {}
        other => panic!("out-of-band must refuse by name, got {other:?}"),
    }
    // The runtime island steps at audio dt with a non-negative
    // dissipation ledger over a strain cycle.
    let mut state = lowered.model.state();
    let dt = 1.0 / 48_000.0;
    for k in 0..960 {
        let t = f64::from(k) * dt;
        let strain = 0.05 * fs_math::det::sin(tau * 220.0 * t);
        lowered.model.step(&mut state, strain, dt);
    }
    assert!(
        state.dissipated >= 0.0,
        "Prony dissipation ledger stays non-negative ({})",
        state.dissipated
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"felt-prony-chain\",\"band_rad_s\":[{lo:.1},{hi:.1}],\
         \"sup_rel_err\":{:.4},\"dissipated_j\":{:.3e}}}",
        lowered.sup_rel_err, state.dissipated
    );
}

// ---------------------------------------------------------------------------
// DONE-WHEN 4 seam: the island geometry is SAMPLED (fs-contact FiniteGap
// chart-to-gap, certified brackets), and the string is stepped under
// `fs_phs::step` (discrete-gradient, energy-consistent) rather than a
// bespoke integrator — the felt island rides the existing machinery.
// ---------------------------------------------------------------------------

/// Hammer crown radius [m] (felt outer surface, AUTHORED fixture value).
const R_CROWN: f64 = 8.0e-3;

/// Sample the felt-crowned hammer surface with the certified chart-to-gap
/// sampler and return the receipt. The crown is a `SphereChart` tangent to
/// the string line at the apex — authored, Estimate-labeled geometry, but
/// the SAMPLING is the real fs-contact pipeline with rigorous brackets.
fn sampled_crown() -> fs_contact::normal_patch::FiniteGapChartSamplingReceipt {
    use fs_contact::normal_patch::{
        FiniteGapChartEvidenceRequirement, FiniteGapChartSamplingRequest, FiniteGapContactFrame,
        FiniteGapGrid, sample_finite_gap_from_chart,
    };
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    use fs_geom::fixtures::SphereChart;
    use fs_geom::{Point3, Vec3};

    let chart = SphereChart {
        center: Point3::new(0.0, 0.0, R_CROWN),
        radius: R_CROWN,
    };
    let request = FiniteGapChartSamplingRequest {
        chart: &chart,
        source_geometry_id: "music/felt-hammer-crown/r8mm/v1",
        grid: FiniteGapGrid {
            cells_x: 5,
            cells_y: 5,
            cell_width_m: 5.0e-4,
            cell_depth_m: 5.0e-4,
        },
        frame: FiniteGapContactFrame {
            surface_point_m: Point3::new(0.0, 0.0, 0.0),
            outward_normal: Vec3::new(0.0, 0.0, -1.0),
            tangent_x: Vec3::new(1.0, 0.0, 0.0),
            tangent_y: Vec3::new(0.0, 1.0, 0.0),
        },
        evidence_requirement: FiniteGapChartEvidenceRequirement::RigorousEnclosure,
        outside_probe_m: 1.0e-3,
        maximum_inward_search_m: 5.0e-3,
        root_tolerance_m: 1.0e-9,
        maximum_bisection_steps: 64,
    };
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 7,
                kernel_id: 11,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        sample_finite_gap_from_chart(&request, &cx).expect("crown sampling certifies")
    })
}

#[test]
fn the_crown_gap_field_is_certified_and_matches_the_sphere_sagitta() {
    use fs_contact::normal_patch::FiniteGapChartSamplingAuthority;
    let receipt = sampled_crown();
    assert_eq!(
        receipt.authority,
        FiniteGapChartSamplingAuthority::Enclosure,
        "crown sampling must be rigorous, not nominal"
    );
    // Every sampled midpoint sits within its own bracket of the analytic
    // sphere sagitta r - sqrt(r^2 - x^2 - y^2) (cell-center convention
    // mirrors the sampler's: (i + 0.5 - cells/2) * width).
    let g = receipt.grid;
    for (index, &mid) in receipt.undeformed_gap_m.iter().enumerate() {
        let ix = index % g.cells_x;
        let iy = index / g.cells_x;
        let x = (ix as f64 + 0.5 - 0.5 * g.cells_x as f64) * g.cell_width_m;
        let y = (iy as f64 + 0.5 - 0.5 * g.cells_y as f64) * g.cell_depth_m;
        let sagitta = R_CROWN - fs_math::det::sqrt(R_CROWN * R_CROWN - x * x - y * y);
        assert!(
            (mid - sagitta).abs() <= receipt.maximum_surface_bracket_m + 1.0e-12,
            "cell ({ix},{iy}): sampled {mid} vs analytic {sagitta}"
        );
    }
    // The apex cells carry the smallest standoff (crown geometry).
    let apex = apex_standoff(&receipt);
    assert!(
        receipt
            .undeformed_gap_m
            .iter()
            .all(|&h| h >= apex - 1.0e-12),
        "apex must be the minimal sampled height"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"felt-crown-gap-field\",\"cells\":{},\
         \"apex_standoff_m\":{apex:.3e},\"max_bracket_m\":{:.3e}}}",
        receipt.undeformed_gap_m.len(),
        receipt.maximum_surface_bracket_m
    );
}

/// Minimal sampled surface height — the apex standoff the island starts from.
fn apex_standoff(receipt: &fs_contact::normal_patch::FiniteGapChartSamplingReceipt) -> f64 {
    receipt
        .undeformed_gap_m
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min)
}

#[test]
fn the_phs_string_island_ledger_closes_with_sampled_geometry() {
    // One string mode as a genuine PortHamiltonian (x = [q, p], H =
    // (omega^2 q^2 + p^2)/2, force port through the mode shape), stepped by
    // the EXISTING discrete-gradient `fs_phs::step`; the initial hammer
    // standoff comes from the certified crown sampling above. The
    // discrete-gradient ledger (delta_h = supplied - dissipated) must close
    // at machine precision EVERY step, and the strike-level energy budget
    // must close: KE_in = KE_out + work-into-string + felt loss.
    use fs_phs::{PortHamiltonian, QuadraticStorage, step};

    let wave_speed = fs_math::det::sqrt(60.0f64 / 6.0e-4);
    let omega = core::f64::consts::PI * wave_speed / 0.65;
    let zeta = 5.0e-4;
    let phi = fs_math::det::sin(core::f64::consts::PI * 0.12);
    let storage = QuadraticStorage::new(vec![omega * omega, 0.0, 0.0, 1.0], 2).expect("Q admits");
    let sys = PortHamiltonian::new(
        2,
        1,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, 2.0 * zeta * omega],
        vec![0.0, phi],
        Box::new(storage),
    )
    .expect("string mode admits as a PHS");

    let standoff = apex_standoff(&sampled_crown());
    let mut hammer = Hammer::Felt {
        law: felt(),
        state: <WoolFelt as Uniaxial>::initial_state(&felt()),
    };
    let v0 = 2.0f64;
    let dt = 1.0 / f64::from(RATE);
    let h = dt / SUBSTEPS as f64;
    let mut x = vec![0.0f64, 0.0];
    let mut y_hammer = -standoff;
    let mut v_hammer = v0;
    let ke_in = 0.5 * M_HAMMER * v0 * v0;
    let mut supplied_total = 0.0f64;
    let mut string_h = 0.0f64;
    let mut contacted = false;
    for _ in 0..(RATE as usize / 25) {
        let string_disp = phi * x[0];
        let mut force_avg = 0.0f64;
        let mut peak_overlap = 0.0f64;
        for _ in 0..SUBSTEPS {
            let overlap = y_hammer - string_disp;
            let f = hammer.force(overlap);
            force_avg += f;
            peak_overlap = peak_overlap.max(overlap);
            v_hammer -= h * f / M_HAMMER;
            y_hammer += h * v_hammer;
        }
        force_avg /= SUBSTEPS as f64;
        hammer.commit(peak_overlap.max(y_hammer - string_disp));
        if force_avg > 0.0 {
            contacted = true;
        }
        let record = step(&sys, &x, &[force_avg], dt).expect("PHS step");
        // Discrete-gradient honesty, every step, machine precision.
        let residual = (record.delta_h - (record.supplied - record.dissipated)).abs();
        assert!(
            residual <= 1.0e-10 * record.supplied.abs().max(1.0),
            "PHS ledger must close per step (residual {residual})"
        );
        supplied_total += record.supplied;
        string_h += record.delta_h;
        x = record.x;
        if v_hammer < 0.0 && (y_hammer - phi * x[0]) < -1.0e-4 {
            break;
        }
    }
    assert!(contacted, "the hammer must actually strike");
    let ke_out = 0.5 * M_HAMMER * v_hammer * v_hammer;
    let felt_loss = ke_in - ke_out - supplied_total;
    assert!(
        felt_loss > 0.02 * ke_in,
        "felt loop must dissipate a real fraction (loss {felt_loss:.3e} of {ke_in:.3e})"
    );
    assert!(
        felt_loss < ke_in && ke_out < ke_in,
        "ledger sanity: no energy minted"
    );
    assert!(
        string_h > 0.0 && string_h <= supplied_total + 1.0e-12,
        "string retains at most what the port supplied (H {string_h:.3e} vs {supplied_total:.3e})"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"felt-phs-island\",\"standoff_m\":{standoff:.3e},\
         \"ke_in_j\":{ke_in:.4e},\"ke_out_j\":{ke_out:.4e},\"work_into_string_j\":{supplied_total:.4e},\
         \"felt_loss_j\":{felt_loss:.4e}}}"
    );
}

/// Composed check for the INGESTED felt geometry (music bead
/// `frankensim-music-v8-root-3ez8g.3.4`): a `fs_query::felt` thickness
/// field drives the whole certified chain — chart -> `AllowEstimate`
/// gap sampling -> convergence-checked dense response curve ->
/// `FiniteGapPoint` — with the authority held at Estimate end to end
/// and the constant-thickness control matching the closed-form
/// cylinder sagitta.
#[test]
fn the_ingested_field_drives_the_certified_response_chain() {
    use fs_contact::normal_patch::{
        FiniteGapChartEvidenceRequirement, FiniteGapChartSamplingAuthority,
        FiniteGapChartSamplingRequest, FiniteGapContactFrame, FiniteGapGrid,
        FiniteGapResponseCurveRequest, NormalPatchLaw, build_finite_gap_response_curve,
        sample_finite_gap_from_chart,
    };
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    use fs_geom::{Point3, Vec3};
    use fs_query::{FeltCoordinateUnits, FeltStation, FeltThicknessChart, FeltThicknessField};

    const CORE_R: f64 = 8.0e-3;
    const T_APEX: f64 = 9.0e-3;
    let field = |crown_c: f64, tag: &str| -> FeltThicknessField {
        let stations: Vec<FeltStation> = (0..41)
            .map(|i| {
                let theta = -1.0 + 0.05 * f64::from(i);
                FeltStation {
                    coordinate: theta,
                    thickness_m: T_APEX - crown_c * theta * theta,
                    half_width_m: 0.0,
                }
            })
            .collect();
        FeltThicknessField::try_new(
            stations,
            FeltCoordinateUnits::RadiansAroundCore,
            tag,
            "authored fixture (this test)",
            "authored analytic profile; Estimate by authorship",
        )
        .expect("field admits")
    };
    let grid = FiniteGapGrid {
        cells_x: 7,
        cells_y: 7,
        cell_width_m: 6.0e-4,
        cell_depth_m: 6.0e-4,
    };
    let sample = |chart: &FeltThicknessChart, grid: FiniteGapGrid, id: &str| {
        let request = FiniteGapChartSamplingRequest {
            chart,
            source_geometry_id: id,
            grid,
            frame: FiniteGapContactFrame {
                surface_point_m: Point3::new(0.0, 0.0, 0.0),
                outward_normal: Vec3::new(0.0, 0.0, -1.0),
                tangent_x: Vec3::new(1.0, 0.0, 0.0),
                tangent_y: Vec3::new(0.0, 1.0, 0.0),
            },
            evidence_requirement: FiniteGapChartEvidenceRequirement::AllowEstimate,
            outside_probe_m: 1.0e-3,
            maximum_inward_search_m: 8.0e-3,
            root_tolerance_m: 1.0e-9,
            maximum_bisection_steps: 64,
        };
        let gate = CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 7,
                    kernel_id: 13,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            sample_finite_gap_from_chart(&request, &cx).expect("ingested chart samples")
        })
    };

    // CONTROL: constant thickness -> the outer surface is a cylinder of
    // radius R0 = core + t and the gap field must match the closed-form
    // sagitta R0 - sqrt(R0^2 - x^2), y-invariant (an oracle wholly
    // independent of the chart code).
    const WIDTH_R: f64 = 9.0e-3;
    let cylinder = FeltThicknessChart::try_new(
        field(0.0, "music/felt-crown/constant-control/v1"),
        CORE_R,
        6.0e-3,
        WIDTH_R,
    )
    .expect("control chart");
    let control = sample(&cylinder, grid, "music/felt-ingest/control/v1");
    assert_eq!(control.authority, FiniteGapChartSamplingAuthority::Estimate);
    let r0 = CORE_R + T_APEX;
    for iy in 0..7usize {
        for ix in 0..7usize {
            let x = (ix as f64 + 0.5 - 3.5) * 6.0e-4;
            let y = (iy as f64 + 0.5 - 3.5) * 6.0e-4;
            // Exact closed form: the width sag shrinks the RADIUS, so
            // the surface along the vertical ray sits at
            // R0 - sqrt((R0 - sag)^2 - x^2).
            let shrunk = r0 - y * y / (2.0 * WIDTH_R);
            let expected = r0 - (shrunk * shrunk - x * x).sqrt();
            let got = control.undeformed_gap_m[iy * 7 + ix];
            assert!(
                (got - expected).abs() < 1.0e-8,
                "cylinder sagitta at cell ({ix},{iy}): got {got:.3e} expected {expected:.3e}"
            );
        }
    }

    // CROWNED: the same core with a real crown curves TIGHTER, so every
    // off-apex column shows a LARGER gap than the constant control.
    let crowned = FeltThicknessChart::try_new(
        field(4.0e-3, "music/felt-crown/authored-parabolic/v1"),
        CORE_R,
        6.0e-3,
        WIDTH_R,
    )
    .expect("crowned chart");
    let receipt = crowned.field_receipt();
    let sampling = sample(&crowned, grid, "music/felt-ingest/crowned/v1");
    assert_eq!(
        sampling.authority,
        FiniteGapChartSamplingAuthority::Estimate
    );
    let apex = sampling.undeformed_gap_m[3 * 7 + 3];
    assert!(apex.abs() < 1.0e-6, "the apex touches the tangent plane");
    let edge_crowned = sampling.undeformed_gap_m[3 * 7];
    let edge_control = control.undeformed_gap_m[3 * 7];
    assert!(
        edge_crowned > edge_control,
        "crowning must open the off-apex gap ({edge_crowned:.3e} vs {edge_control:.3e})"
    );

    // THE CHAIN: receipt -> convergence-checked response curve ->
    // FiniteGapPoint, authority never promoted.
    let build = |sampling: &fs_contact::normal_patch::FiniteGapChartSamplingReceipt,
                 grid: FiniteGapGrid| {
        let request = FiniteGapResponseCurveRequest {
            gap_source_identity: sampling.identity,
            gap_source_authority: sampling.authority,
            grid,
            undeformed_gap_m: sampling.undeformed_gap_m.clone(),
            reduced_modulus_pa: 5.0e8,
            // Dense near zero: the F ~ delta^(3/2) ramp makes the first
            // trapezoid the worst integral segment, and the curve's own
            // energy gate refuses a too-coarse ladder (it did).
            approach_nodes_m: vec![
                0.0, 2.0e-6, 5.0e-6, 1.0e-5, 1.5e-5, 2.2e-5, 3.0e-5, 4.0e-5, 5.0e-5, 6.0e-5,
            ],
            maximum_active_set_iterations: 200,
            complementarity_tolerance_m: 1.0e-12,
            boundary_clearance_cells: 1,
            absolute_energy_tolerance_j: 1.0e-12,
            relative_energy_tolerance: 0.20,
        };
        let gate = CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 7,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            build_finite_gap_response_curve(&request, &cx).expect("response curve certifies")
        })
    };
    let curve = build(&sampling, grid);
    assert_eq!(
        curve.gap_source_authority,
        FiniteGapChartSamplingAuthority::Estimate
    );
    assert_eq!(curve.gap_source_identity, sampling.identity);
    let response = curve.evaluate(4.0e-5).expect("interpolates");
    assert!(response.normal_force_n > 0.0 && response.normal_force_n.is_finite());
    let mut semiaxes = response.equivalent_pressure_semiaxes_m;
    if semiaxes[1] > semiaxes[0] {
        semiaxes.swap(0, 1);
    }
    let law = NormalPatchLaw::FiniteGapPoint {
        response_identity: curve.identity,
        reference_radius_m: 1.16e-2,
        elastic_force_n: response.normal_force_n,
        elastic_tangent_n_per_m: response.normal_tangent_n_per_m,
        reversible_energy_j: response.reversible_energy_j,
        peak_pressure_pa: response.peak_pressure_pa,
        equivalent_pressure_semiaxes_m: semiaxes,
        dissipation_s_per_m: 0.0,
    };
    // The binding: every value in the law came from ONE evaluate call
    // on the identity-carrying curve.
    match law {
        NormalPatchLaw::FiniteGapPoint {
            response_identity,
            elastic_force_n,
            ..
        } => {
            assert_eq!(response_identity, curve.identity);
            assert!((elastic_force_n - response.normal_force_n).abs() == 0.0);
        }
        _ => unreachable!(),
    }

    // REFINEMENT: a finer grid over the same ingested field agrees on
    // the force within the family's own convergence class.
    let fine_grid = FiniteGapGrid {
        cells_x: 11,
        cells_y: 11,
        cell_width_m: 3.8e-4,
        cell_depth_m: 3.8e-4,
    };
    let fine_sampling = sample(&crowned, fine_grid, "music/felt-ingest/crowned-fine/v1");
    let fine_curve = build(&fine_sampling, fine_grid);
    let fine = fine_curve.evaluate(4.0e-5).expect("fine interpolates");
    let rel = ((fine.normal_force_n - response.normal_force_n) / fine.normal_force_n).abs();
    // Measured 2026-08-17: rel = 0.016 (16.56 N vs 16.30 N); band =
    // 3x headroom.
    assert!(
        rel < 0.05,
        "grid refinement must converge (coarse {:.4} N vs fine {:.4} N, rel {rel:.3})",
        response.normal_force_n,
        fine.normal_force_n
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"felt-ingest-chain\",\"verdict\":\"pass\",\
         \"field_digest\":\"{:#018x}\",\"force_n\":{:.4},\"fine_force_n\":{:.4},\
         \"refinement_rel\":{rel:.4},\"authority\":\"Estimate\"}}",
        receipt.digest, response.normal_force_n, fine.normal_force_n
    );
}
