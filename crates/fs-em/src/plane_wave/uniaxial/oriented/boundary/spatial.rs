//! Prepared original-carrier fields and their coherent spatial frequency jets.

use super::super::super::dispersion::kernel::Z;
use super::super::{
    OrientedCrystalField, OrientedCrystalFieldSide, ReflectionError, Result, kernel,
};
use super::{OrientedCrystalBoundaryFrequencyDispersion, OrientedCrystalVectorFrequencyJet};
use fs_math::{c64::C64, dd::Dd};

/// One original physical-point carrier field and its complete Cartesian jets.
///
/// The carrier's coherent Poynting vector is unchanged. Component logarithmic
/// phases belong to their chosen Cartesian projection, not to a common pulse
/// delay for the interfering incident/reflected or ordinary/extraordinary waves.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrientedCrystalPointFrequencyDispersion {
    carrier: OrientedCrystalField,
    vector: OrientedCrystalVectorFrequencyJet,
}

impl OrientedCrystalPointFrequencyDispersion {
    /// Exactly the original prepared field, including coherent Poynting [W/m²].
    pub const fn carrier_field(self) -> OrientedCrystalField {
        self.carrier
    }

    /// Complete E/H frequency jets at this fixed physical point.
    pub const fn vector_frequency_jet(self) -> OrientedCrystalVectorFrequencyJet {
        self.vector
    }
}

impl OrientedCrystalBoundaryFrequencyDispersion<'_, '_> {
    /// Evaluate local frequency derivatives at a fixed Cartesian point [m].
    /// The native half-space selection and exp(-i k.r) convention are retained.
    ///
    /// Both phase branches and complete boundary polarization jets contribute.
    /// FixedAngle follows the original incident bulk law's signed kx derivative;
    /// FixedTransverseWavenumber has exactly zero kx derivatives. Reflected kz
    /// derivatives have the opposite sign to the original incident branch.
    ///
    /// Preparation already retained all original material and wavevector jets.
    /// This query performs no material evaluation, root solve, boundary solve,
    /// finite-frequency difference, frequency fit or allocation. Common complex
    /// phase jets are removed before coherent DD interference; the original
    /// phase reduction and weighted evanescent decay retain weak phase/decay and
    /// derivative tails even when a carrier or one excitation vanishes.
    ///
    /// An unrepresentable weighted frequency polynomial before attenuation,
    /// or relative-decay underflow of a field needed by a later common-phase
    /// product, refuses explicitly. A raw zero field component retains its derivatives,
    /// but has no logarithmic phase/delay/GDD; no shared vector pulse is inferred.
    pub fn field_frequency_dispersion(
        &self,
        side: OrientedCrystalFieldSide,
        point_m: [f64; 3],
    ) -> Result<OrientedCrystalPointFrequencyDispersion> {
        use OrientedCrystalFieldSide::{
            Extraordinary, Incident, Ordinary, Reflected, TotalIncidentSide, Transmitted,
        };
        let boundary = self.boundary_response();
        let carrier = kernel::field(boundary, side, point_m).map_err(ReflectionError::Invalid)?;
        let incident_normal = if matches!(side, Incident | Reflected | TotalIncidentSide) {
            self.incident_normal_derivatives
                .map_err(ReflectionError::Invalid)?
        } else {
            [Z::real(Dd::ZERO); 2]
        };
        let normal_jets = self.normal.normal_derivative_jets;
        let derivatives = [
            self.incident,
            self.reflected,
            self.transmitted[0],
            self.transmitted[1],
        ]
        .map(|field| kernel::boundary::FieldDerivatives {
            electric: field.scaled_electric_derivatives_v_m(),
            magnetic: field.scaled_magnetic_derivatives_a_m(),
        });
        let incident = kernel::FrequencyWave {
            normal: C64::from_re(boundary.incident_normal_wavenumber_rad_m),
            normal_derivatives: incident_normal,
            electric: &boundary.incident_electric_v_m,
            magnetic: &boundary.incident_magnetic_a_m,
            derivatives: &derivatives[0],
        };
        let reflected = kernel::FrequencyWave {
            normal: -incident.normal,
            normal_derivatives: incident_normal.map(|value| value.scale(-Dd::ONE)),
            electric: &boundary.reflected_electric_v_m,
            magnetic: &boundary.reflected_magnetic_a_m,
            derivatives: &derivatives[1],
        };
        let transmitted = |mode: usize| kernel::FrequencyWave {
            normal: boundary.transmitted[mode].normal_wavenumber_rad_m,
            normal_derivatives: normal_jets[mode],
            electric: &boundary.transmitted[mode].electric_v_m,
            magnetic: &boundary.transmitted[mode].magnetic_a_m,
            derivatives: &derivatives[mode + 2],
        };
        let waves = match side {
            Incident => [Some(incident), None],
            Reflected => [Some(reflected), None],
            TotalIncidentSide => [Some(incident), Some(reflected)],
            Ordinary => [Some(transmitted(0)), None],
            Extraordinary => [Some(transmitted(1)), None],
            Transmitted => [Some(transmitted(0)), Some(transmitted(1))],
        };
        let derivatives = kernel::field_frequency_derivatives(
            self.omega_rad_s(),
            boundary.tangential_wavenumber_rad_m,
            self.tangential_derivatives,
            side,
            point_m,
            waves,
        )
        .map_err(ReflectionError::Invalid)?;
        Ok(OrientedCrystalPointFrequencyDispersion {
            vector: OrientedCrystalVectorFrequencyJet::from_scaled(
                self.omega_rad_s(),
                carrier.electric_v_m,
                carrier.magnetic_a_m,
                derivatives,
            )?,
            carrier,
        })
    }
}
