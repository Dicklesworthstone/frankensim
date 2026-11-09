//! Original-source vector Fresnel jets for arbitrary optic-axis orientation.
//!
//! The incident tangential electric phasors are held fixed. Both outgoing
//! polarizations, their coupling, and the complete reflected Cartesian E/H
//! fields are differentiated analytically at the physical interface.

use super::super::UniaxialMode;
use super::super::dispersion::kernel::Z;
use super::super::dispersion::material_derivatives;
use super::{
    OrientedCrystalFrequencyDispersion, OrientedCrystalResponse, OrientedUniaxialInterface,
    PreparedOrientedCrystalField, ReflectionError, Result, dispersion, kernel,
};
use crate::plane_wave::dispersion::{
    FresnelAmplitudeDispersion, FresnelLogarithmicDispersion, IncidenceFrequencyConstraint,
};
use fs_math::{c64::C64, dd::Dd};

mod spatial;
pub use spatial::OrientedCrystalPointFrequencyDispersion;

/// Full Cartesian peak-field jet at one carrier frequency.
///
/// The zeroth-order vectors are the unchanged native fields. Derivatives are
/// with respect to angular frequency under the owner's explicit constraint.
/// The component logarithmic queries describe the chosen Cartesian projection;
/// they do not assign one phase or pulse delay to an interfering vector field.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrientedCrystalVectorFrequencyJet {
    electric: [FresnelAmplitudeDispersion; 3],
    magnetic: [FresnelAmplitudeDispersion; 3],
}

impl OrientedCrystalVectorFrequencyJet {
    /// Original complex electric field [V/m].
    pub fn electric_v_m(self) -> [C64; 3] {
        self.electric.map(FresnelAmplitudeDispersion::amplitude)
    }

    /// Original complex magnetic field [A/m].
    pub fn magnetic_a_m(self) -> [C64; 3] {
        self.magnetic.map(FresnelAmplitudeDispersion::amplitude)
    }

    /// dE/domega [V s/m], including polarization rotation and mixing.
    pub fn electric_derivative_v_s_m(self) -> [C64; 3] {
        self.electric.map(FresnelAmplitudeDispersion::derivative_s)
    }

    /// d²E/domega² [V s²/m].
    pub fn electric_second_derivative_v_s2_m(self) -> [C64; 3] {
        self.electric
            .map(FresnelAmplitudeDispersion::second_derivative_s2)
    }

    /// dH/domega [A s/m], including the changing wave admittance.
    pub fn magnetic_derivative_a_s_m(self) -> [C64; 3] {
        self.magnetic.map(FresnelAmplitudeDispersion::derivative_s)
    }

    /// d²H/domega² [A s²/m].
    pub fn magnetic_second_derivative_a_s2_m(self) -> [C64; 3] {
        self.magnetic
            .map(FresnelAmplitudeDispersion::second_derivative_s2)
    }

    /// [omega E', omega² E''] in original electric-field units [V/m].
    pub fn scaled_electric_derivatives_v_m(self) -> [[C64; 3]; 2] {
        core::array::from_fn(|order| {
            self.electric
                .map(|component| component.scaled_derivatives()[order])
        })
    }

    /// [omega H', omega² H''] in original magnetic-field units [A/m].
    pub fn scaled_magnetic_derivatives_a_m(self) -> [[C64; 3]; 2] {
        core::array::from_fn(|order| {
            self.magnetic
                .map(|component| component.scaled_derivatives()[order])
        })
    }

    /// Separate logarithmic jets of Ex, Ey and Ez. Each entry is None at an
    /// exact native component zero or an unrepresentable quotient. Its raw
    /// first/second field derivatives remain available, including at a zero.
    pub fn electric_component_logarithmic_dispersion(
        self,
    ) -> [Option<FresnelLogarithmicDispersion>; 3] {
        self.electric
            .map(FresnelAmplitudeDispersion::logarithmic_dispersion)
    }

    /// Separate logarithmic jets of Hx, Hy and Hz, with the same zero handling.
    pub fn magnetic_component_logarithmic_dispersion(
        self,
    ) -> [Option<FresnelLogarithmicDispersion>; 3] {
        self.magnetic
            .map(FresnelAmplitudeDispersion::logarithmic_dispersion)
    }

    fn from_scaled(
        omega: f64,
        electric: [C64; 3],
        magnetic: [C64; 3],
        derivatives: kernel::boundary::FieldDerivatives,
    ) -> Result<Self> {
        let e = |j| {
            FresnelAmplitudeDispersion::try_from_scaled(
                electric[j],
                omega,
                [derivatives.electric[0][j], derivatives.electric[1][j]],
            )
        };
        let h = |j| {
            FresnelAmplitudeDispersion::try_from_scaled(
                magnetic[j],
                omega,
                [derivatives.magnetic[0][j], derivatives.magnetic[1][j]],
            )
        };
        Ok(Self {
            electric: [e(0)?, e(1)?, e(2)?],
            magnetic: [h(0)?, h(1)?, h(2)?],
        })
    }
}

/// Source-bound native boundary and full incident/reflected/modal E/H jets.
/// A zero carrier modal excitation does not force its spectral derivative to
/// vanish: frequency-dependent polarization coupling can excite that mode.
#[derive(Clone, Copy)]
pub struct OrientedCrystalBoundaryFrequencyDispersion<'s, 'm> {
    normal: OrientedCrystalFrequencyDispersion<'s, 'm>,
    incident: OrientedCrystalVectorFrequencyJet,
    reflected: OrientedCrystalVectorFrequencyJet,
    transmitted: [OrientedCrystalVectorFrequencyJet; 2],
    tangential_derivatives: [Dd; 2],
    // An unrepresentable incident phase jet need not invalidate the boundary
    // fields or transmitted spatial queries. The incident-side point query
    // reports this cached refusal without reevaluating the source.
    incident_normal_derivatives: core::result::Result<[Z; 2], &'static str>,
}

impl<'s, 'm> OrientedCrystalBoundaryFrequencyDispersion<'s, 'm> {
    pub const fn source(&self) -> &'s OrientedUniaxialInterface<'m> {
        self.normal.source()
    }

    pub const fn omega_rad_s(&self) -> f64 {
        self.normal.omega_rad_s()
    }

    pub const fn frequency_constraint(&self) -> IncidenceFrequencyConstraint {
        self.normal.frequency_constraint()
    }

    pub const fn boundary_response(&self) -> &OrientedCrystalResponse {
        self.normal.boundary_response()
    }

    /// Separate native normal propagation. These modal propagation jets do not
    /// include the vector boundary amplitudes returned by this owner.
    pub const fn normal_dispersion(&self) -> OrientedCrystalFrequencyDispersion<'s, 'm> {
        self.normal
    }

    /// Fixed authored tangential E, with the incident longitudinal E and H
    /// following the original incident law and chosen frequency constraint.
    pub const fn incident_field(&self) -> OrientedCrystalVectorFrequencyJet {
        self.incident
    }

    /// Complete reflected vector, retaining conversion between the two
    /// tangential electric polarizations and complex total-reflection phase.
    pub const fn reflected_field(&self) -> OrientedCrystalVectorFrequencyJet {
        self.reflected
    }

    /// Original ordinary or extraordinary eigenfield and its complete jets.
    /// At admitted exact local isotropy, the native slots remain TE then TM.
    pub const fn transmitted_field(&self, mode: UniaxialMode) -> OrientedCrystalVectorFrequencyJet {
        self.transmitted[match mode {
            UniaxialMode::Ordinary => 0,
            UniaxialMode::Extraordinary => 1,
        }]
    }

    /// Coherent sum of the two transmitted vectors at z=0, including their
    /// frequency-dependent interference. No modal phase averaging is used.
    pub fn transmitted_boundary_field(&self) -> Result<OrientedCrystalVectorFrequencyJet> {
        let [a, b] = self.transmitted;
        let sum = |a: [C64; 3], b: [C64; 3]| core::array::from_fn(|j| a[j] + b[j]);
        OrientedCrystalVectorFrequencyJet::from_scaled(
            self.omega_rad_s(),
            sum(a.electric_v_m(), b.electric_v_m()),
            sum(a.magnetic_a_m(), b.magnetic_a_m()),
            kernel::boundary::FieldDerivatives {
                electric: core::array::from_fn(|order| {
                    sum(
                        a.scaled_electric_derivatives_v_m()[order],
                        b.scaled_electric_derivatives_v_m()[order],
                    )
                }),
                magnetic: core::array::from_fn(|order| {
                    sum(
                        a.scaled_magnetic_derivatives_a_m()[order],
                        b.scaled_magnetic_derivatives_a_m()[order],
                    )
                }),
            },
        )
    }
}

impl<'a> OrientedUniaxialInterface<'a> {
    /// Analytic full vector Fresnel jets of one original arbitrary-axis carrier.
    /// The complex tangential incident [Ex,Ey] is frequency independent.
    ///
    /// Both constraints require the incident and principal materials' original
    /// two-sided lossless frequency jets. Original domains, endpoints and cusps
    /// remain authoritative. Either exactly critical branch refuses, even when
    /// unexcited. At equal principal phase indices with an out-of-plane optic
    /// axis, both local principal jets must coincide so that the native TE/TM
    /// convention is an admitted spectral basis.
    ///
    /// One native boundary solve records a fixed-size local basis witness.
    /// Eigenvector and boundary derivatives are analytic; no neighboring
    /// frequency, fitted material, frequency-dependent gauge normalization or
    /// per-query allocation is used.
    pub fn boundary_frequency_dispersion(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        incident_tangential_e_v_m: [C64; 2],
        constraint: IncidenceFrequencyConstraint,
    ) -> Result<OrientedCrystalBoundaryFrequencyDispersion<'_, 'a>> {
        self.prepare_field(omega_rad_s, incidence_cosine, incident_tangential_e_v_m)?
            .boundary_frequency_dispersion(constraint)
    }

    /// Signed carrier kx/k0 without an angle round trip. FixedAngle follows
    /// that original incident direction; FixedTransverseWavenumber holds its
    /// dimensional transverse wavenumber as frequency changes.
    pub fn boundary_frequency_dispersion_at_tangential_index(
        &self,
        omega_rad_s: f64,
        tangential_index: f64,
        incident_tangential_e_v_m: [C64; 2],
        constraint: IncidenceFrequencyConstraint,
    ) -> Result<OrientedCrystalBoundaryFrequencyDispersion<'_, 'a>> {
        self.prepare_field_at_tangential_index(
            omega_rad_s,
            tangential_index,
            incident_tangential_e_v_m,
        )?
        .boundary_frequency_dispersion(constraint)
    }
}

impl<'s, 'm> PreparedOrientedCrystalField<'s, 'm> {
    /// Differentiate this exact prepared boundary and its native local basis.
    /// The original incidence, eigenvectors, coefficients and illumination were
    /// retained by the one preparation solve. No boundary or carrier-root solve
    /// is repeated, and the prepared coherent carrier field remains unchanged.
    ///
    /// Fixed-kx normal propagation alone retains its broader original material
    /// admission. This additional boundary query requires the incident jet too.
    pub fn boundary_frequency_dispersion(
        &self,
        constraint: IncidenceFrequencyConstraint,
    ) -> Result<OrientedCrystalBoundaryFrequencyDispersion<'s, 'm>> {
        let source = self.source;
        let original = source.resolve_frequency(self.omega_rad_s)?;
        let derivatives = material_derivatives(
            source.incident,
            source.crystal,
            original.omega,
            IncidenceFrequencyConstraint::FixedAngle,
        )?;
        let normal = dispersion::prepare_with_derivatives(
            source,
            original,
            self.incidence,
            self.boundary,
            constraint,
            &derivatives,
        )?;
        let fields = kernel::boundary::derivatives(kernel::boundary::Input {
            original: original.input(
                source.optic_axis,
                self.incidence,
                self.boundary_witness.incident_electric,
            ),
            boundary: &self.boundary,
            witness: &self.boundary_witness,
            relative_permittivity_derivatives: derivatives.relative_permittivity,
            principal_wavenumber_derivatives: derivatives.principal_wavenumber,
            constraint,
        })
        .map_err(ReflectionError::Invalid)?;
        let field = |electric, magnetic, derivatives| {
            OrientedCrystalVectorFrequencyJet::from_scaled(
                original.omega,
                electric,
                magnetic,
                derivatives,
            )
        };
        let incident_bulk = derivatives
            .incident_wavenumber
            .ok_or(ReflectionError::Invalid(
                "original incident bulk frequency derivatives were not retained",
            ))?;
        let tangential_derivatives = if constraint == IncidenceFrequencyConstraint::FixedAngle {
            let sine = Dd::from_f64(self.boundary_witness.basis.p)
                / Dd::from_f64(self.boundary_witness.basis.i);
            incident_bulk.map(|value| Dd::from_f64(value) * sine)
        } else {
            [Dd::ZERO; 2]
        };
        let incident_normal_derivatives = super::super::dispersion::kernel::normal_frequency_jet(
            super::super::dispersion::kernel::Input {
                omega: original.omega,
                indices: [original.ni; 3],
                relative_permittivity_derivatives: [derivatives.relative_permittivity[0]; 3],
                principal_wavenumber_derivatives: [[
                    incident_bulk[0],
                    incident_bulk[1] * original.omega,
                ]; 2],
                incidence: self.incidence,
                axis: [0.0, 0.0, 1.0],
                mode: UniaxialMode::Ordinary,
                kind: super::super::UniaxialWaveKind::Propagating,
                normal_wavenumber: C64::from_re(self.boundary.incident_normal_wavenumber_rad_m),
                constraint,
            },
        );
        Ok(OrientedCrystalBoundaryFrequencyDispersion {
            normal,
            tangential_derivatives,
            incident_normal_derivatives,
            incident: field(
                self.boundary.incident_electric_v_m,
                self.boundary.incident_magnetic_a_m,
                fields.incident,
            )?,
            reflected: field(
                self.boundary.reflected_electric_v_m,
                self.boundary.reflected_magnetic_a_m,
                fields.reflected,
            )?,
            transmitted: [
                field(
                    self.boundary.transmitted[0].electric_v_m,
                    self.boundary.transmitted[0].magnetic_a_m,
                    fields.transmitted[0],
                )?,
                field(
                    self.boundary.transmitted[1].electric_v_m,
                    self.boundary.transmitted[1].magnetic_a_m,
                    fields.transmitted[1],
                )?,
            ],
        })
    }
}
