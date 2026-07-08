//! Fiber-section beam-columns (bead tfz.14): cross-sections
//! discretized into fibers, each carrying a uniaxial fs-material law —
//! Mander confined/unconfined concrete plus Menegotto–Pinto steel for
//! the reinforced-concrete frame use-case. The section state update
//! `ε(y) = ε₀ − y·κ → (N, M, 2×2 tangent)` is the seismic flagship's
//! inner hot loop; [`update_sections_batched`] runs MANY sections per
//! call with the 2×2 tangent solves through fs-la's batched Cholesky
//! (the §15.2 pairing), and the battery ledgers the measured
//! throughput.
//!
//! Cyclic correctness is the point: the hysteresis battery drives an
//! RC section through growing curvature cycles and checks loop
//! closure, positive per-cycle dissipation, and the peak moment
//! against a hand capacity estimate.

use fs_la::batched::{BatchMat, BatchVec, batch_cholesky, batch_cholesky_solve};
use fs_material::fiber::{ManderConcrete, ManderState, MenegottoPintoSteel, MpState, Uniaxial};

/// A fiber's constitutive card.
#[derive(Debug, Clone)]
pub enum FiberLaw {
    /// Menegotto–Pinto steel.
    Steel(MenegottoPintoSteel),
    /// Mander concrete (confined or unconfined by parameters).
    Concrete(ManderConcrete),
}

/// A fiber's committed state.
#[derive(Debug, Clone)]
pub enum FiberState {
    /// Steel state.
    Steel(MpState),
    /// Concrete state.
    Concrete(ManderState),
}

/// One fiber: position, area, law.
#[derive(Debug, Clone)]
pub struct Fiber {
    /// Distance from the section reference axis (bending arm).
    pub y: f64,
    /// Tributary area.
    pub area: f64,
    /// The constitutive card.
    pub law: FiberLaw,
}

/// A section: fibers plus committed states.
#[derive(Debug, Clone)]
pub struct Section {
    /// The fibers.
    pub fibers: Vec<Fiber>,
    /// Committed states, one per fiber.
    pub states: Vec<FiberState>,
}

/// A section's response at a strain state.
#[derive(Debug, Clone, Copy, Default)]
pub struct SectionState {
    /// Axial force.
    pub n: f64,
    /// Bending moment.
    pub m: f64,
    /// Tangent [[dN/dε₀, dN/dκ], [dM/dε₀, dM/dκ]].
    pub tangent: [[f64; 2]; 2],
}

impl Section {
    /// Build with virgin states.
    #[must_use]
    pub fn new(fibers: Vec<Fiber>) -> Section {
        let states = fibers
            .iter()
            .map(|f| match &f.law {
                FiberLaw::Steel(l) => FiberState::Steel(l.initial_state()),
                FiberLaw::Concrete(l) => FiberState::Concrete(l.initial_state()),
            })
            .collect();
        Section { fibers, states }
    }

    /// Response at (ε₀, κ) from the COMMITTED states (pure — commit
    /// separately).
    #[must_use]
    pub fn respond(&self, eps0: f64, kappa: f64) -> SectionState {
        let mut out = SectionState::default();
        for (f, st) in self.fibers.iter().zip(&self.states) {
            let eps = eps0 - f.y * kappa;
            // Section convention is tension-positive; the Mander card
            // is compression-positive (tension returns zero) — flip in
            // and out for concrete.
            let (sig, tan) = match (&f.law, st) {
                (FiberLaw::Steel(l), FiberState::Steel(s)) => (l.stress(eps, s), l.tangent(eps, s)),
                (FiberLaw::Concrete(l), FiberState::Concrete(s)) => {
                    (-l.stress(-eps, s), l.tangent(-eps, s))
                }
                _ => unreachable!("law/state pairing fixed at construction"),
            };
            out.n += sig * f.area;
            out.m -= sig * f.y * f.area;
            let ea = tan * f.area;
            out.tangent[0][0] += ea;
            out.tangent[0][1] -= ea * f.y;
            out.tangent[1][0] -= ea * f.y;
            out.tangent[1][1] += ea * f.y * f.y;
        }
        out
    }

    /// Commit the states at (ε₀, κ).
    pub fn commit(&mut self, eps0: f64, kappa: f64) {
        for (f, st) in self.fibers.iter().zip(self.states.iter_mut()) {
            let eps = eps0 - f.y * kappa;
            match (&f.law, st) {
                (FiberLaw::Steel(l), FiberState::Steel(s)) => *s = l.update_state(eps, s),
                (FiberLaw::Concrete(l), FiberState::Concrete(s)) => {
                    *s = l.update_state(-eps, s);
                }
                _ => unreachable!("law/state pairing fixed at construction"),
            }
        }
    }
}

/// Batched section update: respond to per-section strain pairs and
/// solve each 2×2 tangent system `K·d = rhs` through fs-la's batched
/// Cholesky — the flagship's inner loop shape (state determination
/// needs exactly these solves). Returns the responses and the batch
/// solutions.
#[must_use]
pub fn update_sections_batched(
    sections: &[Section],
    strains: &[(f64, f64)],
    rhs: &[(f64, f64)],
) -> (Vec<SectionState>, Vec<(f64, f64)>) {
    assert_eq!(sections.len(), strains.len(), "one strain pair per section");
    assert_eq!(sections.len(), rhs.len(), "one rhs per section");
    let m = sections.len();
    let responses: Vec<SectionState> = sections
        .iter()
        .zip(strains)
        .map(|(s, &(e, k))| s.respond(e, k))
        .collect();
    let mut a = BatchMat::zeros(2, m);
    let mut b = BatchVec::zeros(2, m);
    for (i, r) in responses.iter().enumerate() {
        // SPD guard: fiber tangents can momentarily lose definiteness
        // on crushing branches; a tiny diagonal floor keeps the batch
        // factorization honest without masking the response itself.
        let floor = 1e-9 * (r.tangent[0][0].abs() + r.tangent[1][1].abs()).max(1e-9);
        a.set(i, 0, 0, r.tangent[0][0] + floor);
        a.set(i, 0, 1, r.tangent[0][1]);
        a.set(i, 1, 0, r.tangent[1][0]);
        a.set(i, 1, 1, r.tangent[1][1] + floor);
        b.set(i, 0, rhs[i].0);
        b.set(i, 1, rhs[i].1);
    }
    let (chol, failures) = batch_cholesky(&a);
    // Indefinite outliers fall back to their own 2×2 direct solve.
    let mut out = Vec::with_capacity(m);
    if failures.is_empty() {
        batch_cholesky_solve(&chol, &mut b);
        for i in 0..m {
            out.push((b.get(i, 0), b.get(i, 1)));
        }
    } else {
        for (i, r) in responses.iter().enumerate() {
            let det = r.tangent[0][0] * r.tangent[1][1] - r.tangent[0][1] * r.tangent[1][0];
            let d = if det.abs() > 1e-30 { det } else { 1e-30 };
            out.push((
                (r.tangent[1][1] * rhs[i].0 - r.tangent[0][1] * rhs[i].1) / d,
                (-r.tangent[1][0] * rhs[i].0 + r.tangent[0][0] * rhs[i].1) / d,
            ));
        }
    }
    (responses, out)
}

/// A reinforced-concrete rectangular section fixture: confined core +
/// unconfined cover Mander fibers plus top/bottom Menegotto–Pinto
/// steel layers. Depth `d`, width `b`, `layers` concrete fibers.
///
/// # Panics
/// On invalid material parameters (fixture constants are valid).
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn rc_section(d: f64, b: f64, layers: usize, steel_area: f64) -> Section {
    let core = ManderConcrete::new(42e6, 0.004, 32e9, 0.015).expect("core card");
    let cover = ManderConcrete::new(30e6, 0.002, 30e9, 0.005).expect("cover card");
    let steel = MenegottoPintoSteel::new(200e9, 450e6, 0.01).expect("steel card");
    let cover_t = 0.1 * d;
    let mut fibers = Vec::new();
    let t = d / layers as f64;
    for i in 0..layers {
        let y = -0.5 * d + (i as f64 + 0.5) * t;
        let is_cover = y.abs() > 0.5 * d - cover_t;
        fibers.push(Fiber {
            y,
            area: b * t,
            law: FiberLaw::Concrete(if is_cover {
                cover.clone()
            } else {
                core.clone()
            }),
        });
    }
    for y in [-(0.5 * d - cover_t), 0.5 * d - cover_t] {
        fibers.push(Fiber {
            y,
            area: steel_area,
            law: FiberLaw::Steel(steel.clone()),
        });
    }
    Section::new(fibers)
}
