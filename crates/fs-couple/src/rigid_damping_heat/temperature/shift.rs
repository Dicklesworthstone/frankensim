//! Calibrated Kelvin viscosity shifts on the original passive spatial graph.
//!
//! The material owner's convention is a_T = tau(T)/tau(T_ref). A common shift
//! of a thermorheologically simple Kelvin law therefore gives C(T) = a_T C_ref.
//! Full anisotropy, cross-coupling and exact null directions stay in the original
//! native matrix. No interface temperature, universal calibration, material
//! source validity, Maxwell memory or elastic parameter work is inferred.

use super::TemperatureDampingLaw;
use fs_material::visco::{ViscoError, temperature::TemperatureShift};
use fs_mbd::spatial_loads::internal_impulse::viscous::{
    BushingMatrixError, PassiveBushingMatrix, network::ViscousNetworkLink,
};

/// Failure leaves the reference calibration and every physical state untouched.
#[derive(Debug)]
pub enum KelvinTemperatureShiftError {
    Invalid(&'static str),
    /// The supplied connector no longer has its original reference form.
    Reference {
        connector: usize,
    },
    Shift(ViscoError),
    Matrix(BushingMatrixError),
}

impl core::fmt::Display for KelvinTemperatureShiftError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(message) => write!(f, "Kelvin damping calibration: {message}"),
            Self::Reference { connector } => {
                write!(f, "connector {connector} changed its reference damping law")
            }
            Self::Shift(error) => error.fmt(f),
            Self::Matrix(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for KelvinTemperatureShiftError {}

/// Borrowed original reference laws and one explicit calibration per connector.
///
/// References are declared at each calibration's own reference temperature.
/// Retain this binding and the original links across continuation calls. The
/// accepted endpoint matrices are physical outputs, never new reference laws.
/// An accidental attempt to rescale an already shifted output refuses.
///
/// The caller qualifies temperature superposition for these same materials and
/// the fixed attachment geometry. Coefficients and calibrated support come only
/// from the existing WLF/Arrhenius owner. Preparation performs no allocation or
/// factorization; evaluation scales the original passive factors in fixed-size
/// storage. The thermal consumer caches unchanged-temperature evaluations.
#[derive(Clone, Copy, Debug)]
pub struct CalibratedKelvinDamping<'a> {
    reference: &'a [ViscousNetworkLink],
    shifts: &'a [TemperatureShift],
}

impl<'a> CalibratedKelvinDamping<'a> {
    pub fn try_new(
        reference: &'a [ViscousNetworkLink],
        shifts: &'a [TemperatureShift],
    ) -> Result<Self, KelvinTemperatureShiftError> {
        if reference.len() != shifts.len() {
            return Err(KelvinTemperatureShiftError::Invalid(
                "one temperature calibration per original connector required",
            ));
        }
        Ok(Self { reference, shifts })
    }

    pub fn reference_links(&self) -> &[ViscousNetworkLink] {
        self.reference
    }

    pub fn shifts(&self) -> &[TemperatureShift] {
        self.shifts
    }
}

impl TemperatureDampingLaw for CalibratedKelvinDamping<'_> {
    type Error = KelvinTemperatureShiftError;

    fn at_temperature(
        &self,
        connector: usize,
        reference: &PassiveBushingMatrix,
        temperature_k: f64,
    ) -> Result<PassiveBushingMatrix, Self::Error> {
        let original = self
            .reference
            .get(connector)
            .ok_or(KelvinTemperatureShiftError::Reference { connector })?;
        if original.damping != *reference {
            return Err(KelvinTemperatureShiftError::Reference { connector });
        }
        let state = self.shifts[connector]
            .evaluate(temperature_k)
            .map_err(KelvinTemperatureShiftError::Shift)?;
        original
            .damping
            .scaled(state.factor())
            .map_err(KelvinTemperatureShiftError::Matrix)
    }
}
