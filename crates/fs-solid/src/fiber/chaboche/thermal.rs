//! Exact one-through-eight-branch dispatch for sourced thermal Chaboche fibers.
//!
//! Each variant retains the original prepared interval and the native const-N
//! history. Cloning shares material preparation; mechanical trials neither
//! allocate branch arrays nor combine, pad or reset independent backstresses.

use super::super::FiberLaw;
use fs_blake3::ContentHash;
use fs_material::IsotropicElastic;
use fs_material::fiber::Uniaxial;
use fs_material::plastic::chaboche::ChabocheError;
use fs_material::plastic::energy::StoredEnergy;
use fs_material::plastic::thermal::chaboche::ThermalChabocheError;
use fs_material::plastic::thermal::chaboche::axial::{
    UniaxialThermalChaboche, UniaxialThermalChabocheState,
};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::Voigt;

/// Original prepared source interval, with its actual ordered branch count.
#[derive(Clone, Debug)]
pub enum ThermalChabocheFiberLaw {
    One(UniaxialThermalChaboche<1>),
    Two(UniaxialThermalChaboche<2>),
    Three(UniaxialThermalChaboche<3>),
    Four(UniaxialThermalChaboche<4>),
    Five(UniaxialThermalChaboche<5>),
    Six(UniaxialThermalChaboche<6>),
    Seven(UniaxialThermalChaboche<7>),
    Eight(UniaxialThermalChaboche<8>),
}

/// Accepted temperature and all original native axial branch histories.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermalChabocheFiberState {
    One(UniaxialThermalChabocheState<1>),
    Two(UniaxialThermalChabocheState<2>),
    Three(UniaxialThermalChabocheState<3>),
    Four(UniaxialThermalChabocheState<4>),
    Five(UniaxialThermalChabocheState<5>),
    Six(UniaxialThermalChabocheState<6>),
    Seven(UniaxialThermalChabocheState<7>),
    Eight(UniaxialThermalChabocheState<8>),
}

macro_rules! each_law {
    ($value:expr, |$law:ident| $body:expr) => {
        match $value {
            ThermalChabocheFiberLaw::One($law) => $body,
            ThermalChabocheFiberLaw::Two($law) => $body,
            ThermalChabocheFiberLaw::Three($law) => $body,
            ThermalChabocheFiberLaw::Four($law) => $body,
            ThermalChabocheFiberLaw::Five($law) => $body,
            ThermalChabocheFiberLaw::Six($law) => $body,
            ThermalChabocheFiberLaw::Seven($law) => $body,
            ThermalChabocheFiberLaw::Eight($law) => $body,
        }
    };
}
macro_rules! each_state {
    ($value:expr, |$state:ident| $body:expr) => {
        match $value {
            ThermalChabocheFiberState::One($state) => $body,
            ThermalChabocheFiberState::Two($state) => $body,
            ThermalChabocheFiberState::Three($state) => $body,
            ThermalChabocheFiberState::Four($state) => $body,
            ThermalChabocheFiberState::Five($state) => $body,
            ThermalChabocheFiberState::Six($state) => $body,
            ThermalChabocheFiberState::Seven($state) => $body,
            ThermalChabocheFiberState::Eight($state) => $body,
        }
    };
}
macro_rules! paired {
    ($law:expr, $state:expr, |$l:ident, $s:ident| $body:expr) => {
        match ($law, $state) {
            (ThermalChabocheFiberLaw::One($l), ThermalChabocheFiberState::One($s)) => $body,
            (ThermalChabocheFiberLaw::Two($l), ThermalChabocheFiberState::Two($s)) => $body,
            (ThermalChabocheFiberLaw::Three($l), ThermalChabocheFiberState::Three($s)) => $body,
            (ThermalChabocheFiberLaw::Four($l), ThermalChabocheFiberState::Four($s)) => $body,
            (ThermalChabocheFiberLaw::Five($l), ThermalChabocheFiberState::Five($s)) => $body,
            (ThermalChabocheFiberLaw::Six($l), ThermalChabocheFiberState::Six($s)) => $body,
            (ThermalChabocheFiberLaw::Seven($l), ThermalChabocheFiberState::Seven($s)) => $body,
            (ThermalChabocheFiberLaw::Eight($l), ThermalChabocheFiberState::Eight($s)) => $body,
            _ => Err(branch_mismatch()),
        }
    };
}
fn branch_mismatch() -> ThermalChabocheError {
    ChabocheError::Invalid("thermal fiber law and original branch-history counts differ").into()
}

impl ThermalChabocheFiberLaw {
    pub fn branch_count(&self) -> usize {
        match self {
            Self::One(_) => 1,
            Self::Two(_) => 2,
            Self::Three(_) => 3,
            Self::Four(_) => 4,
            Self::Five(_) => 5,
            Self::Six(_) => 6,
            Self::Seven(_) => 7,
            Self::Eight(_) => 8,
        }
    }
    pub fn old_temperature_k(&self) -> f64 {
        each_law!(self, |law| law.old_temperature_k())
    }
    pub fn temperature_k(&self) -> f64 {
        each_law!(self, |law| law.temperature_k())
    }
    pub fn duration_s(&self) -> f64 {
        each_law!(self, |law| law.duration_s())
    }
    pub fn density_kg_m3(&self) -> f64 {
        each_law!(self, |law| law.density_kg_m3())
    }
    pub fn caloric_curve_identity(&self) -> ContentHash {
        each_law!(self, |law| law.caloric_curve_identity())
    }
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        each_law!(self, |law| law.material())
    }
    pub fn elastic(&self) -> &IsotropicElastic {
        each_law!(self, |law| &law.tensor_law().law().backbone().backbone().elastic)
    }
    pub fn yield_stress_pa(&self) -> f64 {
        each_law!(self, |law| law.tensor_law().law().backbone().backbone().yield_stress)
    }
    pub fn initial_state(&self) -> ThermalChabocheFiberState {
        match self {
            Self::One(law) => ThermalChabocheFiberState::One(law.initial_state()),
            Self::Two(law) => ThermalChabocheFiberState::Two(law.initial_state()),
            Self::Three(law) => ThermalChabocheFiberState::Three(law.initial_state()),
            Self::Four(law) => ThermalChabocheFiberState::Four(law.initial_state()),
            Self::Five(law) => ThermalChabocheFiberState::Five(law.initial_state()),
            Self::Six(law) => ThermalChabocheFiberState::Six(law.initial_state()),
            Self::Seven(law) => ThermalChabocheFiberState::Seven(law.initial_state()),
            Self::Eight(law) => ThermalChabocheFiberState::Eight(law.initial_state()),
        }
    }
    pub fn tensor_strain(
        &self,
        strain: f64,
        state: &ThermalChabocheFiberState,
    ) -> Result<Voigt, ThermalChabocheError> {
        paired!(self, state, |law, old| law.tensor_strain(strain, old))
    }
    pub fn stored_energy(
        &self,
        strain: f64,
        state: &ThermalChabocheFiberState,
    ) -> Result<StoredEnergy, ThermalChabocheError> {
        paired!(self, state, |law, old| law.stored_energy(strain, old))
    }
    pub fn response(
        &self,
        strain: f64,
        state: &ThermalChabocheFiberState,
    ) -> Result<(f64, f64), ThermalChabocheError> {
        paired!(self, state, |law, old| law
            .response(strain, old)
            .map(|response| (response.stress_pa, response.tangent_pa)))
    }
    pub(crate) fn commit(
        &self,
        strain: f64,
        state: &mut ThermalChabocheFiberState,
    ) -> Result<(), ThermalChabocheError> {
        paired!(self, state, |law, old| {
            *old = law.response(strain, old)?.state;
            Ok(())
        })
    }
    /// Keep the original branch ordering, source choices and caloric reference.
    pub fn admit_successor(&self, next: &Self) -> Result<(), ThermalChabocheError> {
        let result = match (self, next) {
            (Self::One(old), Self::One(next)) => old.admit_successor(next),
            (Self::Two(old), Self::Two(next)) => old.admit_successor(next),
            (Self::Three(old), Self::Three(next)) => old.admit_successor(next),
            (Self::Four(old), Self::Four(next)) => old.admit_successor(next),
            (Self::Five(old), Self::Five(next)) => old.admit_successor(next),
            (Self::Six(old), Self::Six(next)) => old.admit_successor(next),
            (Self::Seven(old), Self::Seven(next)) => old.admit_successor(next),
            (Self::Eight(old), Self::Eight(next)) => old.admit_successor(next),
            _ => return Err(branch_mismatch()),
        };
        result.map_err(ThermalChabocheError::from)
    }
}

impl ThermalChabocheFiberState {
    pub fn branch_count(&self) -> usize {
        match self {
            Self::One(_) => 1,
            Self::Two(_) => 2,
            Self::Three(_) => 3,
            Self::Four(_) => 4,
            Self::Five(_) => 5,
            Self::Six(_) => 6,
            Self::Seven(_) => 7,
            Self::Eight(_) => 8,
        }
    }
    pub fn temperature_k(&self) -> f64 {
        each_state!(self, |state| state.temperature_k())
    }
    pub fn plastic_strain(&self) -> f64 {
        each_state!(self, |state| state.mechanical.plastic_strain())
    }
    pub fn equivalent_plastic_strain(&self) -> f64 {
        each_state!(self, |state| state.mechanical.equivalent_plastic_strain())
    }
    /// Read a physical branch translation without allocating or advancing time.
    pub fn backstress_pa(&self, branch: usize) -> Option<f64> {
        each_state!(self, |state| state.mechanical.tensor_history().backstress_pa
            .get(branch)
            .map(|value| 1.5 * value[0]))
    }
}

macro_rules! from_native {
    ($n:literal, $variant:ident) => {
        impl From<UniaxialThermalChaboche<$n>> for ThermalChabocheFiberLaw {
            fn from(value: UniaxialThermalChaboche<$n>) -> Self {
                Self::$variant(value)
            }
        }
        impl From<UniaxialThermalChaboche<$n>> for FiberLaw {
            fn from(value: UniaxialThermalChaboche<$n>) -> Self {
                Self::ThermalChaboche(ThermalChabocheFiberLaw::$variant(value))
            }
        }
    };
}
from_native!(1, One);
from_native!(2, Two);
from_native!(3, Three);
from_native!(4, Four);
from_native!(5, Five);
from_native!(6, Six);
from_native!(7, Seven);
from_native!(8, Eight);

impl From<ThermalChabocheFiberLaw> for FiberLaw {
    fn from(law: ThermalChabocheFiberLaw) -> Self {
        Self::ThermalChaboche(law)
    }
}
