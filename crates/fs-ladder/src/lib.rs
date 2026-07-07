//! fs-ladder — the fidelity-ladder registry (plan addendum, Proposal 3).
//! Layer: L3 (FLUX-adjacent; fidelity rungs are a physics concept, but this
//! crate is a pure abstraction with no numerical dependencies — concrete
//! transfer operators are injected).
//!
//! Each physics kernel declares its RUNGS — an ordered ladder of fidelity
//! levels (e.g. correlation → RANS → LES; linear-elastic → hyperelastic →
//! plasticity) — with typed PROLONGATION (coarse→fine) and RESTRICTION
//! (fine→coarse) operators between ADJACENT rungs. This is the ONE shared
//! substrate that many addendum proposals walk, so no single one should own
//! it:
//! - Proposal 9's coarse-rung proposer + verifier (prolongate rung k−1 → k);
//! - Proposal 8's planner ladder-walk (climb / descend rungs);
//! - Proposal D's Goodhart guard step (i) "re-solve at rung k+1";
//! - Proposal 11's tolerance-band-extreme re-solves;
//! - Proposal 3's discrepancy probes (adjacent-rung evaluation).
//!
//! Rungs are TOTALLY ORDERED per kernel (index 0 = coarsest/cheapest). The
//! registry owns the rung DECLARATIONS and the adjacency/ordering; the actual
//! numerical [`Transfer`] between two rungs is a pluggable trait object a
//! consumer supplies (fs-feec transfer operators, a correlation model, …).
//! One concrete demonstrator ships — [`Refine1d`], a 1D coarsen/refine-by-2
//! whose `restrict ∘ prolongate = identity` makes the ladder immediately real
//! and its G0 approximation property testable.
//!
//! Determinism: rung resolution and transfer application are pure functions
//! (no RNG, no I/O), so a replayed ladder walk reproduces bit-identical
//! results — load-bearing because transfer outputs feed verified-color
//! certificates.

use std::collections::BTreeMap;
use std::error::Error;
use std::fmt;

/// A numerical transfer between two adjacent rungs. `prolongate` maps a
/// coarse-rung state to the next finer rung; `restrict` maps a fine-rung
/// state to the next coarser rung. Implementations bring the actual operator;
/// the ladder only sequences them.
pub trait Transfer {
    /// Coarse (rung k) → fine (rung k+1).
    fn prolongate(&self, coarse: &[f64]) -> Vec<f64>;
    /// Fine (rung k+1) → coarse (rung k).
    fn restrict(&self, fine: &[f64]) -> Vec<f64>;
}

/// One fidelity rung: its position, a name, a relative-cost hint, and a note.
#[derive(Debug, Clone, PartialEq)]
pub struct Rung {
    /// Position in the ladder (0 = coarsest/cheapest).
    pub index: u32,
    /// Human/agent-readable rung name (e.g. `"RANS"`).
    pub name: String,
    /// Relative cost hint vs the coarsest rung (advisory; the planner learns
    /// real costs from telemetry).
    pub relative_cost: f64,
    /// What this rung models / its validity note.
    pub note: String,
}

/// The finer/coarser neighbours of a rung (either may be absent at the ends).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdjacentRungs {
    /// The next coarser rung's index, or `None` at the bottom.
    pub coarser: Option<u32>,
    /// The next finer rung's index, or `None` at the top.
    pub finer: Option<u32>,
}

/// A structured ladder error (a refusal that teaches).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LadderError {
    /// No rung with this index exists in the kernel's ladder.
    NoSuchRung {
        /// The kernel.
        kernel: String,
        /// The requested index.
        index: u32,
        /// The number of rungs that exist.
        len: u32,
    },
    /// `prolongate` was called at the finest rung (there is no finer rung).
    AtTop {
        /// The kernel.
        kernel: String,
        /// The (top) index.
        index: u32,
    },
    /// `restrict` was called at the coarsest rung (there is no coarser rung).
    AtBottom {
        /// The kernel.
        kernel: String,
        /// The (bottom) index.
        index: u32,
    },
    /// No kernel with this name is registered.
    NoKernel {
        /// The requested kernel.
        kernel: String,
    },
}

impl fmt::Display for LadderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LadderError::NoSuchRung { kernel, index, len } => write!(
                f,
                "kernel '{kernel}' has no rung {index} (it has {len} rung(s), indices 0..{len}); \
                 fix: request a rung in range"
            ),
            LadderError::AtTop { kernel, index } => write!(
                f,
                "kernel '{kernel}' rung {index} is the finest — cannot prolongate above the top; \
                 fix: stop climbing, or add a finer rung"
            ),
            LadderError::AtBottom { kernel, index } => write!(
                f,
                "kernel '{kernel}' rung {index} is the coarsest — cannot restrict below the bottom; \
                 fix: stop descending, or add a coarser rung"
            ),
            LadderError::NoKernel { kernel } => {
                write!(f, "no kernel '{kernel}' is registered; fix: register its ladder first")
            }
        }
    }
}

impl Error for LadderError {}

/// A per-kernel fidelity ladder: an ordered list of rungs plus the transfer
/// between each adjacent pair (`transfers[k]` bridges rung `k` and `k+1`).
pub struct Ladder {
    kernel: String,
    rungs: Vec<Rung>,
    transfers: Vec<Box<dyn Transfer>>,
}

impl Ladder {
    /// Start a ladder for `kernel` with its coarsest (base) rung.
    #[must_use]
    pub fn new(kernel: &str, base_name: &str, base_cost: f64, base_note: &str) -> Ladder {
        Ladder {
            kernel: kernel.to_string(),
            rungs: vec![Rung {
                index: 0,
                name: base_name.to_string(),
                relative_cost: base_cost,
                note: base_note.to_string(),
            }],
            transfers: Vec::new(),
        }
    }

    /// Append the next finer rung, reached from the current top by `transfer`.
    /// Keeps `transfers.len() == rungs.len() - 1` by construction.
    #[must_use]
    pub fn then(
        mut self,
        transfer: Box<dyn Transfer>,
        name: &str,
        cost: f64,
        note: &str,
    ) -> Ladder {
        let index = self.rungs.len() as u32;
        self.rungs.push(Rung {
            index,
            name: name.to_string(),
            relative_cost: cost,
            note: note.to_string(),
        });
        self.transfers.push(transfer);
        self
    }

    /// The kernel name.
    #[must_use]
    pub fn kernel(&self) -> &str {
        &self.kernel
    }

    /// The number of rungs.
    #[must_use]
    pub fn len(&self) -> u32 {
        self.rungs.len() as u32
    }

    /// Is the ladder empty? (Never true via the public builder — a ladder
    /// always has at least its base rung.)
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.rungs.is_empty()
    }

    /// The coarsest rung (index 0).
    #[must_use]
    pub fn bottom(&self) -> &Rung {
        &self.rungs[0]
    }

    /// The finest rung.
    #[must_use]
    pub fn top(&self) -> &Rung {
        &self.rungs[self.rungs.len() - 1]
    }

    /// The rung at `index`.
    ///
    /// # Errors
    /// [`LadderError::NoSuchRung`] if out of range.
    pub fn rung(&self, index: u32) -> Result<&Rung, LadderError> {
        self.rungs.get(index as usize).ok_or_else(|| LadderError::NoSuchRung {
            kernel: self.kernel.clone(),
            index,
            len: self.len(),
        })
    }

    /// The finer/coarser neighbours of a rung (empty in the correct direction
    /// at the ends).
    ///
    /// # Errors
    /// [`LadderError::NoSuchRung`] if `index` is out of range.
    pub fn adjacent_rungs(&self, index: u32) -> Result<AdjacentRungs, LadderError> {
        if index >= self.len() {
            return Err(LadderError::NoSuchRung {
                kernel: self.kernel.clone(),
                index,
                len: self.len(),
            });
        }
        Ok(AdjacentRungs {
            coarser: index.checked_sub(1),
            finer: if index + 1 < self.len() { Some(index + 1) } else { None },
        })
    }

    /// Prolongate a state from rung `from` up to rung `from + 1`.
    ///
    /// # Errors
    /// [`LadderError::NoSuchRung`] if `from` is out of range;
    /// [`LadderError::AtTop`] if `from` is already the finest rung.
    pub fn prolongate(&self, from: u32, coarse: &[f64]) -> Result<Vec<f64>, LadderError> {
        if from >= self.len() {
            return Err(LadderError::NoSuchRung {
                kernel: self.kernel.clone(),
                index: from,
                len: self.len(),
            });
        }
        let t = self.transfers.get(from as usize).ok_or_else(|| LadderError::AtTop {
            kernel: self.kernel.clone(),
            index: from,
        })?;
        Ok(t.prolongate(coarse))
    }

    /// Restrict a state from rung `from` down to rung `from - 1`.
    ///
    /// # Errors
    /// [`LadderError::NoSuchRung`] if `from` is out of range;
    /// [`LadderError::AtBottom`] if `from` is already the coarsest rung.
    pub fn restrict(&self, from: u32, fine: &[f64]) -> Result<Vec<f64>, LadderError> {
        if from >= self.len() {
            return Err(LadderError::NoSuchRung {
                kernel: self.kernel.clone(),
                index: from,
                len: self.len(),
            });
        }
        let Some(below) = from.checked_sub(1) else {
            return Err(LadderError::AtBottom {
                kernel: self.kernel.clone(),
                index: from,
            });
        };
        Ok(self.transfers[below as usize].restrict(fine))
    }
}

/// A registry of per-kernel fidelity ladders — the ONE service every
/// consumer queries.
#[derive(Default)]
pub struct LadderRegistry {
    ladders: BTreeMap<String, Ladder>,
}

impl LadderRegistry {
    /// An empty registry.
    #[must_use]
    pub fn new() -> LadderRegistry {
        LadderRegistry {
            ladders: BTreeMap::new(),
        }
    }

    /// Register a kernel's ladder (replacing any prior one for that kernel).
    pub fn register(&mut self, ladder: Ladder) {
        self.ladders.insert(ladder.kernel.clone(), ladder);
    }

    /// A kernel's ladder.
    ///
    /// # Errors
    /// [`LadderError::NoKernel`] if the kernel is not registered.
    pub fn ladder(&self, kernel: &str) -> Result<&Ladder, LadderError> {
        self.ladders.get(kernel).ok_or_else(|| LadderError::NoKernel {
            kernel: kernel.to_string(),
        })
    }

    /// Registered kernel names, sorted.
    #[must_use]
    pub fn kernels(&self) -> Vec<&str> {
        self.ladders.keys().map(String::as_str).collect()
    }

    /// A registry seeded with the conjugate-heat-transfer (electronics
    /// cooling) ladder from Proposal 7 — the correlation-based bottom rung
    /// makes the fidelity ladder immediately real. Rungs:
    /// `correlation-Nu` (Nusselt correlation) → `RANS` → `LES`.
    #[must_use]
    pub fn cht() -> LadderRegistry {
        let mut r = LadderRegistry::new();
        let ladder = Ladder::new(
            "cht",
            "correlation-Nu",
            1.0,
            "cheap bottom rung: forced-convection Nusselt correlation",
        )
        .then(Box::new(Refine1d), "RANS", 40.0, "steady RANS CFD")
        .then(Box::new(Refine1d), "LES", 2000.0, "large-eddy simulation");
        r.register(ladder);
        r
    }
}

/// A 1D coarsen/refine-by-2 transfer demonstrator. Prolongation is linear
/// interpolation (`n → 2n-1` points: fine `2i` copies coarse `i`, fine `2i+1`
/// averages neighbours); restriction is injection (`2n-1 → n`, coarse `i` =
/// fine `2i`). Hence `restrict ∘ prolongate = identity` on the coarse space
/// (the declared approximation property) and `prolongate ∘ restrict` is an
/// idempotent projection.
pub struct Refine1d;

impl Transfer for Refine1d {
    fn prolongate(&self, coarse: &[f64]) -> Vec<f64> {
        if coarse.is_empty() {
            return Vec::new();
        }
        if coarse.len() == 1 {
            return vec![coarse[0]];
        }
        let mut fine = Vec::with_capacity(coarse.len() * 2 - 1);
        for i in 0..coarse.len() - 1 {
            fine.push(coarse[i]);
            fine.push(0.5 * (coarse[i] + coarse[i + 1]));
        }
        fine.push(coarse[coarse.len() - 1]);
        fine
    }

    fn restrict(&self, fine: &[f64]) -> Vec<f64> {
        if fine.is_empty() {
            return Vec::new();
        }
        // injection: take every other sample (indices 0, 2, 4, …).
        fine.iter().step_by(2).copied().collect()
    }
}
