//! Bounded dispatch of the exact native multi-backstress axial restrictions.
//!
//! Each active variant keeps its original branch count and ordered histories.
//! Responses and commits neither allocate nor pad inactive constitutive branches.
//! Enum storage is sized for the largest variant; this is not compact per-N storage.

use super::FiberLaw;
use fs_material::MaterialError;
use fs_material::fiber::Uniaxial;
use fs_material::plastic::chaboche::axial::{UniaxialChaboche, UniaxialChabocheState};
use fs_material::plastic::chaboche::{ChabocheError, ChabocheJ2};

/// Two through eight independent native backstress branches. The existing
/// single-branch FiberLaw::Chaboche variant remains the N=1 representation.
#[derive(Clone, Debug, PartialEq)]
pub enum ChabocheFiberLaw {
    /// Two ordered recovery branches.
    Two(UniaxialChaboche<2>),
    /// Three ordered recovery branches.
    Three(UniaxialChaboche<3>),
    /// Four ordered recovery branches.
    Four(UniaxialChaboche<4>),
    /// Five ordered recovery branches.
    Five(UniaxialChaboche<5>),
    /// Six ordered recovery branches.
    Six(UniaxialChaboche<6>),
    /// Seven ordered recovery branches.
    Seven(UniaxialChaboche<7>),
    /// Eight ordered recovery branches.
    Eight(UniaxialChaboche<8>),
}

/// Exact native axial histories, with independent state in every branch.
#[derive(Clone, Debug, PartialEq)]
pub enum ChabocheFiberState {
    /// Two ordered histories.
    Two(UniaxialChabocheState<2>),
    /// Three ordered histories.
    Three(UniaxialChabocheState<3>),
    /// Four ordered histories.
    Four(UniaxialChabocheState<4>),
    /// Five ordered histories.
    Five(UniaxialChabocheState<5>),
    /// Six ordered histories.
    Six(UniaxialChabocheState<6>),
    /// Seven ordered histories.
    Seven(UniaxialChabocheState<7>),
    /// Eight ordered histories.
    Eight(UniaxialChabocheState<8>),
}

macro_rules! paired {
    ($law:expr, $state:expr, |$l:ident, $s:ident| $body:expr, $otherwise:expr) => {
        match ($law, $state) {
            (ChabocheFiberLaw::Two($l), ChabocheFiberState::Two($s)) => $body,
            (ChabocheFiberLaw::Three($l), ChabocheFiberState::Three($s)) => $body,
            (ChabocheFiberLaw::Four($l), ChabocheFiberState::Four($s)) => $body,
            (ChabocheFiberLaw::Five($l), ChabocheFiberState::Five($s)) => $body,
            (ChabocheFiberLaw::Six($l), ChabocheFiberState::Six($s)) => $body,
            (ChabocheFiberLaw::Seven($l), ChabocheFiberState::Seven($s)) => $body,
            (ChabocheFiberLaw::Eight($l), ChabocheFiberState::Eight($s)) => $body,
            _ => $otherwise,
        }
    };
}

impl ChabocheFiberLaw {
    /// Number of physical recovery branches evaluated by this variant.
    #[must_use]
    pub const fn branch_count(&self) -> usize {
        match self {
            Self::Two(_) => 2,
            Self::Three(_) => 3,
            Self::Four(_) => 4,
            Self::Five(_) => 5,
            Self::Six(_) => 6,
            Self::Seven(_) => 7,
            Self::Eight(_) => 8,
        }
    }

    pub(super) fn initial_state(&self) -> ChabocheFiberState {
        match self {
            Self::Two(l) => ChabocheFiberState::Two(l.initial_state()),
            Self::Three(l) => ChabocheFiberState::Three(l.initial_state()),
            Self::Four(l) => ChabocheFiberState::Four(l.initial_state()),
            Self::Five(l) => ChabocheFiberState::Five(l.initial_state()),
            Self::Six(l) => ChabocheFiberState::Six(l.initial_state()),
            Self::Seven(l) => ChabocheFiberState::Seven(l.initial_state()),
            Self::Eight(l) => ChabocheFiberState::Eight(l.initial_state()),
        }
    }

    pub(super) fn response(
        &self,
        strain: f64,
        state: &ChabocheFiberState,
    ) -> Result<(f64, f64), ChabocheError> {
        paired!(
            self,
            state,
            |law, old| law
                .response(strain, old)
                .map(|r| (r.stress_pa, r.tangent_pa)),
            Err(ChabocheError::Invalid(
                "fiber law and backstress-history branch counts differ"
            ))
        )
    }

    pub(super) fn commit(&self, strain: f64, state: &mut ChabocheFiberState) {
        paired!(
            self,
            state,
            |law, old| *old = law.update_state(strain, old),
            unreachable!("law/state pairing fixed at construction")
        )
    }
}

impl FiberLaw {
    /// Admit an authored native axial law with one through eight branches.
    /// Construction preserves exact coefficients and ordering. Subsequent
    /// responses dispatch directly to the original const-generic owner.
    pub fn chaboche_branches<const N: usize>(
        law: UniaxialChaboche<N>,
    ) -> Result<Self, MaterialError> {
        Ok(match N {
            1 => Self::Chaboche(select::<N, 1>(&law)?),
            2 => Self::ChabocheBranches(ChabocheFiberLaw::Two(select::<N, 2>(&law)?)),
            3 => Self::ChabocheBranches(ChabocheFiberLaw::Three(select::<N, 3>(&law)?)),
            4 => Self::ChabocheBranches(ChabocheFiberLaw::Four(select::<N, 4>(&law)?)),
            5 => Self::ChabocheBranches(ChabocheFiberLaw::Five(select::<N, 5>(&law)?)),
            6 => Self::ChabocheBranches(ChabocheFiberLaw::Six(select::<N, 6>(&law)?)),
            7 => Self::ChabocheBranches(ChabocheFiberLaw::Seven(select::<N, 7>(&law)?)),
            8 => Self::ChabocheBranches(ChabocheFiberLaw::Eight(select::<N, 8>(&law)?)),
            _ => {
                return Err(MaterialError::Parameters {
                    what: "fiber Chaboche laws require one through eight recovery branches".into(),
                });
            }
        })
    }
}

// Rust cannot refine a const-generic N from a runtime match. Rebuild the same
// admitted owner once at construction, copying each ordered parameter exactly.
fn select<const N: usize, const M: usize>(
    law: &UniaxialChaboche<N>,
) -> Result<UniaxialChaboche<M>, MaterialError> {
    if N != M {
        return Err(MaterialError::Parameters {
            what: "fiber Chaboche branch count must match the native law".into(),
        });
    }
    let original = law.backbone();
    UniaxialChaboche::new(ChabocheJ2::new(
        original.backbone().clone(),
        std::array::from_fn(|i| original.branches()[i]),
    )?)
}
