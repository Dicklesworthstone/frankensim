//! Local spectral propagation of the original refracted crystal eigenmode.
//!
//! A normal-transfer delay depends on what is held fixed while frequency
//! changes. The same explicit constraint as isotropic Fresnel dispersion is
//! retained here. The normal response keeps exp(-i*kz*z) propagation separate
//! from the source-bound Fresnel jets in [`UniaxialBoundaryFrequencyDispersion`].
//! Neither supplies finite-crystal residence time or a shared delay for a
//! coherent superposition of ordinary and extraordinary fields.

use super::kernel::Incidence;
use super::{
    DielectricResponse, Polarization, PrincipalPlaneUniaxialInterface, ReflectionError,
    UniaxialInterfaceResponse, UniaxialOpticalResponse,
};
use crate::plane_wave::dispersion::{
    BulkPropagation, IncidenceFrequencyConstraint, propagation_at_distance,
};
use fs_math::c64::C64;

pub(super) mod kernel;

mod boundary;
pub use boundary::UniaxialBoundaryFrequencyDispersion;

type Result<T> = core::result::Result<T, ReflectionError>;

/// Original principal-plane carrier plus its normal spatial-frequency jet.
/// The source and optic axis remain borrowed unchanged. Local derivatives
/// preserve the caller's fixed-angle or fixed-transverse-wavenumber path.
#[derive(Clone, Copy)]
pub struct UniaxialFrequencyDispersion<'s, 'm> {
    source: &'s PrincipalPlaneUniaxialInterface<'m>,
    omega_rad_s: f64,
    constraint: IncidenceFrequencyConstraint,
    incidence: Incidence,
    indices: [f64; 3],
    carrier: UniaxialInterfaceResponse,
    normal_derivatives: [C64; 2],
}

impl<'s, 'm> UniaxialFrequencyDispersion<'s, 'm> {
    pub const fn source(self) -> &'s PrincipalPlaneUniaxialInterface<'m> {
        self.source
    }
    pub const fn angular_frequency_rad_s(self) -> f64 {
        self.omega_rad_s
    }
    pub const fn frequency_constraint(self) -> IncidenceFrequencyConstraint {
        self.constraint
    }
    /// Unchanged native Snell/Fresnel fields, phase roots and normal powers.
    pub const fn carrier(self) -> UniaxialInterfaceResponse {
        self.carrier
    }
    /// Complex d(kz)/domega [s/m] at fixed transverse observation coordinate.
    /// This is not slowness along the crystal's generally different energy ray.
    pub const fn transmitted_normal_wavenumber_derivative_s_m(self) -> C64 {
        self.normal_derivatives[0]
    }
    /// Complex d²(kz)/domega² [s²/m]. Even nondispersive crystal laws have
    /// geometric normal-transfer dispersion at fixed transverse wavenumber.
    pub const fn transmitted_normal_wavenumber_second_derivative_s2_m(self) -> C64 {
        self.normal_derivatives[1]
    }
    /// Native carrier propagation through a nonnegative normal distance [m],
    /// with local phase/group delay, GDD and evanescent attenuation slopes.
    /// The original Fresnel boundary amplitude is not multiplied into this
    /// transfer; its separate phase jet is available through
    /// [`Self::boundary_frequency_dispersion`].
    ///
    /// Signed slopes are retained, including tilted negative phase roots.
    /// They do not establish signal-front velocity or broadband pulse motion.
    pub fn transmitted_propagation(self, normal_distance_m: f64) -> Result<BulkPropagation> {
        propagation_at_distance(
            self.carrier.transmitted_normal_wavenumber_rad_m,
            self.omega_rad_s,
            self.normal_derivatives,
            normal_distance_m,
        )
    }
}

impl<'a> PrincipalPlaneUniaxialInterface<'a> {
    /// Differentiate the original transmitted normal branch analytically.
    /// The crystal's original two-sided spectral jets retain both principal
    /// slopes/curvatures, support, interpolation gaps and extinction refusals.
    /// Fixed-angle queries also require the incident law's original lossless
    /// frequency jet; fixed transverse wavenumber needs no incident derivative.
    ///
    /// Exact critical roots have no admitted square-root jet. Their carrier
    /// response remains available. Preparation allocates nothing and performs
    /// no frequency differencing or second Maxwell boundary solve. DD
    /// contractions retain cancellation residuals, but near mode coalescence
    /// the accuracy is still limited by the resolved native carrier.
    pub fn frequency_dispersion(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        polarization: Polarization,
        constraint: IncidenceFrequencyConstraint,
    ) -> Result<UniaxialFrequencyDispersion<'_, 'a>> {
        self.prepare_frequency_dispersion(
            omega_rad_s,
            Incidence::Cosine(incidence_cosine),
            polarization,
            constraint,
        )
    }

    /// Explicit signed carrier kx/k0 without an angle round-trip. FixedAngle
    /// retains the incident direction defined by this original carrier;
    /// FixedTransverseWavenumber retains the dimensional kx instead.
    pub fn frequency_dispersion_at_tangential_index(
        &self,
        omega_rad_s: f64,
        tangential_index: f64,
        polarization: Polarization,
        constraint: IncidenceFrequencyConstraint,
    ) -> Result<UniaxialFrequencyDispersion<'_, 'a>> {
        self.prepare_frequency_dispersion(
            omega_rad_s,
            Incidence::TangentialIndex(tangential_index),
            polarization,
            constraint,
        )
    }

    fn prepare_frequency_dispersion(
        &self,
        omega_rad_s: f64,
        incidence: Incidence,
        polarization: Polarization,
        constraint: IncidenceFrequencyConstraint,
    ) -> Result<UniaxialFrequencyDispersion<'_, 'a>> {
        let (carrier, indices) =
            self.evaluate_with_indices(omega_rad_s, incidence, polarization)?;
        let derivatives =
            material_derivatives(self.incident, self.crystal, omega_rad_s, constraint)?;
        prepare_normal(
            self,
            omega_rad_s,
            incidence,
            constraint,
            carrier,
            indices,
            &derivatives,
        )
    }
}

fn prepare_normal<'s, 'm>(
    source: &'s PrincipalPlaneUniaxialInterface<'m>,
    omega_rad_s: f64,
    incidence: Incidence,
    constraint: IncidenceFrequencyConstraint,
    carrier: UniaxialInterfaceResponse,
    indices: [f64; 3],
    derivatives: &MaterialDerivatives,
) -> Result<UniaxialFrequencyDispersion<'s, 'm>> {
    let normal_derivatives = kernel::normal_derivatives(kernel::Input {
        omega: omega_rad_s,
        indices,
        relative_permittivity_derivatives: derivatives.relative_permittivity,
        principal_wavenumber_derivatives: derivatives.principal_wavenumber,
        incidence,
        axis: [source.optic_axis[0], 0.0, source.optic_axis[1]],
        mode: carrier.mode,
        kind: carrier.wave_kind,
        normal_wavenumber: carrier.transmitted_normal_wavenumber_rad_m,
        constraint,
    })
    .map_err(ReflectionError::Invalid)?;
    Ok(UniaxialFrequencyDispersion {
        source,
        omega_rad_s,
        constraint,
        incidence,
        indices,
        carrier,
        normal_derivatives,
    })
}

pub(super) struct MaterialDerivatives {
    pub(super) relative_permittivity: [[f64; 2]; 3],
    pub(super) principal_wavenumber: [[f64; 2]; 2],
}

// One source evaluation supplies both eigenbranches. This is also used by the
// arbitrary-axis owner; it does not reinterpret spectral data through a fit.
pub(super) fn material_derivatives(
    incident: DielectricResponse<'_>,
    crystal: UniaxialOpticalResponse<'_>,
    omega: f64,
    constraint: IncidenceFrequencyConstraint,
) -> Result<MaterialDerivatives> {
    let jets = crystal.transparent_frequency_jets(omega)?;
    let incident_relative = if constraint == IncidenceFrequencyConstraint::FixedAngle {
        let bulk = incident.bulk_dispersion(omega)?;
        relative_permittivity(
            bulk.permittivity_f_m(),
            bulk.scaled_permittivity_derivatives_f_m(),
        )?
    } else {
        [0.0; 2]
    };
    let ordinary = relative_permittivity(jets[0].epsilon, jets[0].epsilon_derivatives)?;
    let extraordinary = relative_permittivity(jets[1].epsilon, jets[1].epsilon_derivatives)?;
    Ok(MaterialDerivatives {
        relative_permittivity: [incident_relative, ordinary, extraordinary],
        principal_wavenumber: jets.map(|jet| [jet.first_s_m.re, jet.second_scaled_s_m.re]),
    })
}

fn relative_permittivity(epsilon: C64, derivatives: [C64; 2]) -> Result<[f64; 2]> {
    if !epsilon.re.is_finite()
        || epsilon.re <= 0.0
        || epsilon.im != 0.0
        || derivatives.iter().any(|value| value.im != 0.0)
    {
        return Err(ReflectionError::Invalid(
            "crystal normal dispersion requires lossless original frequency jets",
        ));
    }
    let result = derivatives.map(|value| value.re / epsilon.re);
    if result.iter().any(|value| !value.is_finite())
        || result
            .iter()
            .zip(derivatives)
            .any(|(&ratio, value)| ratio == 0.0 && value.re != 0.0)
    {
        return Err(ReflectionError::Invalid(
            "crystal relative permittivity derivatives are unrepresentable",
        ));
    }
    Ok(result)
}
