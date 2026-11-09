//! Original-source Fresnel amplitudes of the principal-plane crystal port.
//!
//! Normal propagation and boundary phase are distinct contributions. The
//! original source, authored incidence coordinate and unchanged native carrier
//! remain attached to both; no neighboring frequency or fitted law is used.

use super::{
    Incidence, MaterialDerivatives, Polarization, PrincipalPlaneUniaxialInterface, ReflectionError,
    Result, UniaxialFrequencyDispersion, UniaxialInterfaceResponse, kernel, material_derivatives,
    prepare_normal,
};
use crate::plane_wave::dispersion::{
    FresnelAmplitudeDispersion, IncidenceFrequencyConstraint, PowerFractionDispersion,
};
use crate::plane_wave::uniaxial::UniaxialWaveKind;

/// Unchanged crystal carrier, its normal propagation jet and the local
/// reflection/transmission field jets at the physical interface z=0.
/// The fixed tangential electric axes and normalized power-conjugate magnetic
/// field are exactly those of the original principal-plane Maxwell owner.
#[derive(Clone, Copy)]
pub struct UniaxialBoundaryFrequencyDispersion<'s, 'm> {
    normal: UniaxialFrequencyDispersion<'s, 'm>,
    amplitudes: [FresnelAmplitudeDispersion; 3],
}

impl<'s, 'm> UniaxialBoundaryFrequencyDispersion<'s, 'm> {
    pub const fn source(self) -> &'s PrincipalPlaneUniaxialInterface<'m> {
        self.normal.source()
    }
    pub const fn angular_frequency_rad_s(self) -> f64 {
        self.normal.angular_frequency_rad_s()
    }
    pub const fn frequency_constraint(self) -> IncidenceFrequencyConstraint {
        self.normal.frequency_constraint()
    }
    pub const fn carrier(self) -> UniaxialInterfaceResponse {
        self.normal.carrier()
    }
    /// Same original mode and coordinate path, including complex evanescent
    /// propagation. Boundary phase is not included in this separate transfer.
    pub const fn normal_dispersion(self) -> UniaxialFrequencyDispersion<'s, 'm> {
        self.normal
    }
    /// Raw complex reflection derivatives remain defined at a Brewster zero.
    /// The amplitude owner's logarithmic query returns None there.
    pub const fn reflected_tangential_e(self) -> FresnelAmplitudeDispersion {
        self.amplitudes[0]
    }
    pub const fn transmitted_tangential_e(self) -> FresnelAmplitudeDispersion {
        self.amplitudes[1]
    }
    /// H_transmitted/(Y_incident E_incident), including the incident
    /// admittance's frequency dependence under either coordinate constraint.
    pub const fn transmitted_conjugate_h(self) -> FresnelAmplitudeDispersion {
        self.amplitudes[2]
    }

    /// Native reflected normal power and its complete first frequency slope.
    /// Total internal reflection has exactly zero power slope while its
    /// complex field can have a nonzero phase delay and GDD.
    pub fn reflected_power_fraction(self) -> Result<PowerFractionDispersion> {
        let r = self.amplitudes[0];
        let carrier = self.normal.carrier;
        let result = if carrier.wave_kind == UniaxialWaveKind::Evanescent {
            PowerFractionDispersion::from_scaled(
                carrier.reflected_power_fraction,
                self.normal.omega_rad_s,
                0.0,
            )
        } else {
            PowerFractionDispersion::from_fields(
                carrier.reflected_power_fraction,
                self.normal.omega_rad_s,
                r.amplitude(),
                r.amplitude(),
                r.scaled_derivatives()[0],
                r.scaled_derivatives()[0],
            )
        };
        result.ok_or(ReflectionError::Invalid(
            "crystal reflected Poynting-flux derivative is unrepresentable",
        ))
    }

    /// Native transmitted normal power and its complete first frequency
    /// slope, retaining both electric and magnetic admittance derivatives.
    /// An evanescent half-space remains an exact zero-normal-power port.
    pub fn transmitted_power_fraction(self) -> Result<PowerFractionDispersion> {
        let [_, e, h] = self.amplitudes;
        let carrier = self.normal.carrier;
        let result = if carrier.wave_kind == UniaxialWaveKind::Evanescent {
            PowerFractionDispersion::from_scaled(
                carrier.transmitted_power_fraction,
                self.normal.omega_rad_s,
                0.0,
            )
        } else {
            PowerFractionDispersion::from_fields(
                carrier.transmitted_power_fraction,
                self.normal.omega_rad_s,
                e.amplitude(),
                h.amplitude(),
                e.scaled_derivatives()[0],
                h.scaled_derivatives()[0],
            )
        };
        result.ok_or(ReflectionError::Invalid(
            "crystal transmitted Poynting-flux derivative is unrepresentable",
        ))
    }
}

impl<'a> PrincipalPlaneUniaxialInterface<'a> {
    /// Differentiate the original crystal Fresnel boundary and normal
    /// propagation analytically under the explicit local frequency path.
    /// Both constraints require the incident frequency jet because its port
    /// admittance varies even when dimensional kx is held fixed.
    ///
    /// The crystal and incident jets must stay lossless through both derivative
    /// orders. Original endpoints, interpolation cusps and gaps refuse. An
    /// exact critical crystal root refuses; carrier-only response remains
    /// available. Preparation allocates nothing and solves one native boundary.
    /// A local phase delay is not a signal speed or finite-crystal residence time.
    pub fn boundary_frequency_dispersion(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        polarization: Polarization,
        constraint: IncidenceFrequencyConstraint,
    ) -> Result<UniaxialBoundaryFrequencyDispersion<'_, 'a>> {
        self.prepare_boundary_frequency_dispersion(
            omega_rad_s,
            Incidence::Cosine(incidence_cosine),
            polarization,
            constraint,
        )
    }

    /// Signed carrier kx/k0 counterpart, retaining the authored coordinate
    /// without converting through a rounded angle or changing the optic axis.
    pub fn boundary_frequency_dispersion_at_tangential_index(
        &self,
        omega_rad_s: f64,
        tangential_index: f64,
        polarization: Polarization,
        constraint: IncidenceFrequencyConstraint,
    ) -> Result<UniaxialBoundaryFrequencyDispersion<'_, 'a>> {
        self.prepare_boundary_frequency_dispersion(
            omega_rad_s,
            Incidence::TangentialIndex(tangential_index),
            polarization,
            constraint,
        )
    }

    fn prepare_boundary_frequency_dispersion(
        &self,
        omega_rad_s: f64,
        incidence: Incidence,
        polarization: Polarization,
        constraint: IncidenceFrequencyConstraint,
    ) -> Result<UniaxialBoundaryFrequencyDispersion<'_, 'a>> {
        let (carrier, indices) =
            self.evaluate_with_indices(omega_rad_s, incidence, polarization)?;
        // FixedAngle requests the incident jet from the shared material
        // helper; the actual caller constraint still governs both kernels.
        let derivatives = material_derivatives(
            self.incident,
            self.crystal,
            omega_rad_s,
            IncidenceFrequencyConstraint::FixedAngle,
        )?;
        let normal = prepare_normal(
            self,
            omega_rad_s,
            incidence,
            constraint,
            carrier,
            indices,
            &derivatives,
        )?;
        prepare(normal, &derivatives)
    }
}

impl<'s, 'm> UniaxialFrequencyDispersion<'s, 'm> {
    /// Add Fresnel boundary jets to this exact existing normal-propagation
    /// response. Its original incidence and native carrier remain unchanged;
    /// no Maxwell boundary or outgoing-root solve is repeated.
    ///
    /// This additional query requires the incident law's original frequency
    /// jet under both constraints. Normal-only fixed-kx queries retain their
    /// broader admission and remain independently usable.
    pub fn boundary_frequency_dispersion(
        self,
    ) -> Result<UniaxialBoundaryFrequencyDispersion<'s, 'm>> {
        let derivatives = material_derivatives(
            self.source.incident,
            self.source.crystal,
            self.omega_rad_s,
            IncidenceFrequencyConstraint::FixedAngle,
        )?;
        prepare(self, &derivatives)
    }
}

fn prepare<'s, 'm>(
    normal: UniaxialFrequencyDispersion<'s, 'm>,
    derivatives: &MaterialDerivatives,
) -> Result<UniaxialBoundaryFrequencyDispersion<'s, 'm>> {
    let carrier = normal.carrier;
    let scaled = kernel::boundary::amplitude_derivatives(kernel::boundary::Input {
        indices: normal.indices,
        relative_permittivity_derivatives: derivatives.relative_permittivity,
        incidence: normal.incidence,
        axis: normal.source.optic_axis,
        constraint: normal.constraint,
        carrier,
    })
    .map_err(ReflectionError::Invalid)?;
    let omega = normal.omega_rad_s;
    let amplitudes = [
        FresnelAmplitudeDispersion::try_from_scaled(carrier.reflected_tangential_e, omega, scaled)?,
        FresnelAmplitudeDispersion::try_from_scaled(
            carrier.transmitted_tangential_e,
            omega,
            scaled,
        )?,
        FresnelAmplitudeDispersion::try_from_scaled(
            carrier.transmitted_conjugate_h,
            omega,
            scaled.map(core::ops::Neg::neg),
        )?,
    ];
    Ok(UniaxialBoundaryFrequencyDispersion { normal, amplitudes })
}
