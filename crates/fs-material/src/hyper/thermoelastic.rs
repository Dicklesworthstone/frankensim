//! Finite thermoelastic derivatives of an explicitly prescribed natural map.
//!
//! With H=G(T)^-1, Fe=F H and Jg=det(G), the original-volume Helmholtz
//! potential is Psi(F,T)=Jg We(Fe). Native energy, Piola stress and Hessian
//! supply its exact chain rule. The ordinary entry holds elastic coefficients
//! fixed; the source child additionally carries the original Neo-Hookean
//! coefficient derivatives without replacing this natural-map chain rule.
//!
//! Mechanical entropy is -Psi_T, internal energy is Psi-T Psi_T, and the
//! mechanical contribution to fixed-F heat capacity is -T Psi_TT. These
//! are additions to a separately admitted caloric reference, not a complete
//! caloric law, inferred heat, or a temperature/history evolution algorithm.
//!
//! For beta=dP/dT and positive total reference-volume capacity c_F, eliminating
//! temperature from ds=(c_F/T)dT-beta:dF at ds=0 gives
//! A_S=A_T+(T/c_F) beta (x) beta. This is a full compressible 3-D local
//! tangent. The viscoelastic child uses the same chain rule with accepted
//! viscous metrics held fixed. Constrained pressure and history evolution
//! require their own thermodynamic owners. Source smoothness and the actual
//! meaning of the caloric capacity remain explicit.

use super::acoustic::HyperelasticAcousticState;
use super::reference::ReferenceDeformation;
use super::{Hyperelastic, det3};
use crate::MaterialError;
use fs_math::dd::Dd;

/// Exact elastic-volume plane stress with moving natural volume and thermal work.
pub mod incompressible_plane_stress;
/// Finite-compliance transverse relaxation of the same thermal potential.
pub mod plane_stress;
/// Original elastic-source derivatives for the admitted Neo-Hookean model.
pub mod source;

const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

fn invalid(what: &str) -> MaterialError {
    MaterialError::Parameters {
        what: format!("finite thermoelasticity: {what}"),
    }
}
fn unrepresentable() -> MaterialError {
    MaterialError::State {
        what: "finite thermoelastic response is unrepresentable".into(),
    }
}
fn product(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3).fold(0.0, |sum, k| {
            a[3 * (ij / 3) + k].mul_add(b[3 * k + ij % 3], sum)
        })
    })
}
fn trace(a: &[f64; 9]) -> f64 {
    a[0] + a[4] + a[8]
}
fn dot(a: &[f64; 9], b: &[f64; 9]) -> f64 {
    a.iter().zip(b).fold(0.0, |sum, (x, y)| x.mul_add(*y, sum))
}

fn elastic_deformation_rate(
    elastic_deformation: &[f64; 9],
    reference_jet: &ThermalReferenceJet,
    deformation_rate: &[f64; 9],
    temperature_rate_k_s: f64,
) -> [f64; 9] {
    let h = reference_jet.reference.elastic_deformation(&IDENTITY);
    let thermal_direction = product(elastic_deformation, &reference_jet.first_derivative_per_k);
    let numerator = core::array::from_fn(|i| {
        (-temperature_rate_k_s).mul_add(thermal_direction[i], deformation_rate[i])
    });
    product(&numerator, &h)
}

// Partial derivatives at fixed Fe, supplied only by the admitted source owner.
// The natural-reference derivatives below remain with this shared chain rule.
#[derive(Default)]
struct ElasticTemperatureDerivatives {
    energy_j_m3_k: f64,
    energy_j_m3_k2: f64,
    piola_pa_k: [f64; 9],
}

// An original-volume physical potential at fixed natural reference/history.
// The originating constitutive owner, rather than this chain-rule kernel,
// admits its parameters and supplies its exact gradient and Hessian.
pub(crate) struct ReferencePotentialDerivatives {
    pub(crate) free_energy_j_m3: f64,
    pub(crate) first_piola_pa: [f64; 9],
    pub(crate) isothermal_tangent_pa: [[f64; 9]; 9],
}

/// A local natural deformation and its first two temperature derivatives.
/// All entries are row-major and use the same fixed original/intermediate
/// axes. The owner supplying this jet must establish a smooth local source
/// neighborhood; finite entries alone do not establish that physical claim.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalReferenceJet {
    reference: ReferenceDeformation,
    first_derivative_per_k: [f64; 9],
    second_derivative_per_k2: [f64; 9],
}

impl ThermalReferenceJet {
    pub fn try_new(
        reference: ReferenceDeformation,
        first_derivative_per_k: [f64; 9],
        second_derivative_per_k2: [f64; 9],
    ) -> Result<Self, MaterialError> {
        if first_derivative_per_k
            .iter()
            .chain(&second_derivative_per_k2)
            .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "finite natural-map temperature derivatives required",
            ));
        }
        Ok(Self {
            reference,
            first_derivative_per_k,
            second_derivative_per_k2,
        })
    }

    pub const fn reference(&self) -> &ReferenceDeformation {
        &self.reference
    }
    pub const fn first_derivative_per_k(&self) -> &[f64; 9] {
        &self.first_derivative_per_k
    }
    pub const fn second_derivative_per_k2(&self) -> &[f64; 9] {
        &self.second_derivative_per_k2
    }
}

/// Immutable original-volume mechanical thermodynamics at one total F and T.
/// Internal variables, when supplied by the native viscoelastic owner, remain
/// frozen throughout these reversible derivatives; no relaxation heat or
/// history update is represented by this snapshot.
/// Preparation admits the elastic domain, including the original native
/// constitutive switching conventions. It does not certify a two-sided
/// tangent at a fibre activation surface or an arbitrary source kink.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteThermoelasticState {
    // Only elastic endpoints retain the direct stable increment owner. A
    // fixed-Cv potential must not fall back to its equilibrium branch alone.
    law: Option<Hyperelastic>,
    deformation_gradient: [f64; 9],
    elastic_deformation_gradient: [f64; 9],
    temperature_k: f64,
    reference_jet: ThermalReferenceJet,
    free_energy_j_m3: f64,
    entropy_j_m3_k: f64,
    internal_energy_j_m3: f64,
    mechanical_heat_capacity_j_m3_k: f64,
    first_piola_pa: [f64; 9],
    isothermal_tangent_pa: [[f64; 9]; 9],
    stress_temperature_derivative_pa_k: [f64; 9],
}

/// Instantaneous reversible rates of one admitted thermal snapshot.
/// Deformation gradients are dimensionless, Piola stress is per ORIGINAL
/// area, and energy/entropy densities are per ORIGINAL volume. These are
/// directional derivatives, not a time integrator or dissipative heat.
/// Accepted internal variables remain fixed; temperature-dependent elastic
/// coefficients and the natural map retain the snapshot's source derivatives.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteThermoelasticRate {
    deformation_gradient_rate_per_s: [f64; 9],
    elastic_deformation_gradient_rate_per_s: [f64; 9],
    temperature_rate_k_s: f64,
    first_piola_rate_pa_s: [f64; 9],
    free_energy_rate_w_m3: f64,
    entropy_rate_w_m3_k: f64,
    internal_energy_rate_w_m3: f64,
}

impl FiniteThermoelasticRate {
    /// Total Fdot in row-major order [1/s].
    pub const fn deformation_gradient_rate_per_s(&self) -> &[f64; 9] {
        &self.deformation_gradient_rate_per_s
    }
    /// Elastic Fedot=(Fdot-Fe*Gdot)*G^-1 in the original natural axes [1/s].
    pub const fn elastic_deformation_gradient_rate_per_s(&self) -> &[f64; 9] {
        &self.elastic_deformation_gradient_rate_per_s
    }
    pub const fn temperature_rate_k_s(&self) -> f64 {
        self.temperature_rate_k_s
    }
    /// Material time derivative of first Piola stress, not an objective
    /// spatial stress rate. A changing observer also contributes Qdot*P.
    pub const fn first_piola_rate_pa_s(&self) -> &[f64; 9] {
        &self.first_piola_rate_pa_s
    }
    /// Mechanical Helmholtz rate P:Fdot-S*Tdot [W/m³].
    pub const fn free_energy_rate_w_m3(&self) -> f64 {
        self.free_energy_rate_w_m3
    }
    /// Reversible mechanical entropy rate [W/(m³ K)]; no caloric or
    /// irreversible entropy production is included.
    pub const fn entropy_rate_w_m3_k(&self) -> f64 {
        self.entropy_rate_w_m3_k
    }
    /// Mechanical internal-energy rate P:Fdot+T*Sdot [W/m³].
    pub const fn internal_energy_rate_w_m3(&self) -> f64 {
        self.internal_energy_rate_w_m3
    }
}

impl FiniteThermoelasticState {
    /// Evaluate the original full elastic energy and its exact thermal chain
    /// rule. The elastic stretch envelope applies to Fe=F G^-1; total F has
    /// its own orientation/representability gate. Caloric capacity is not
    /// needed to observe mechanical entropy/internal energy.
    pub fn try_new(
        law: &Hyperelastic,
        total_deformation_gradient: &[f64; 9],
        temperature_k: f64,
        reference_jet: ThermalReferenceJet,
    ) -> Result<Self, MaterialError> {
        Self::with_elastic_temperature_derivatives(
            law,
            total_deformation_gradient,
            temperature_k,
            reference_jet,
            |_| Ok(ElasticTemperatureDerivatives::default()),
        )
    }

    fn with_elastic_temperature_derivatives(
        law: &Hyperelastic,
        total_deformation_gradient: &[f64; 9],
        temperature_k: f64,
        reference_jet: ThermalReferenceJet,
        derivatives: impl FnOnce(&[f64; 9]) -> Result<ElasticTemperatureDerivatives, MaterialError>,
    ) -> Result<Self, MaterialError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(invalid("positive finite absolute temperature required"));
        }
        let jacobian = det3(total_deformation_gradient);
        if total_deformation_gradient
            .iter()
            .any(|value| !value.is_finite())
            || !jacobian.is_finite()
            || jacobian <= 0.0
        {
            return Err(invalid(
                "finite orientation-preserving total deformation required",
            ));
        }
        // Public native parameters are mutable, so admit them once before
        // constructing this immutable thermal/acoustic snapshot.
        let law = Hyperelastic::new(law.model.clone(), law.stretch_limit)?;
        let reference = reference_jet.reference;
        let elastic_deformation_gradient =
            reference.elastic_deformation(total_deformation_gradient);
        let elastic_piola = law.piola(&elastic_deformation_gradient)?;
        let elastic_tangent = law.tangent(&elastic_deformation_gradient)?;
        let elastic_energy = law.energy(&elastic_deformation_gradient);
        let partial = derivatives(&elastic_deformation_gradient)?;
        if [partial.energy_j_m3_k, partial.energy_j_m3_k2]
            .iter()
            .chain(&partial.piola_pa_k)
            .any(|value| !value.is_finite())
        {
            return Err(unrepresentable());
        }
        let mut state = Self::with_reference_potential(
            total_deformation_gradient,
            temperature_k,
            reference_jet,
            ReferencePotentialDerivatives {
                free_energy_j_m3: reference.volume_ratio() * elastic_energy,
                first_piola_pa: reference.pullback_piola(&elastic_piola),
                isothermal_tangent_pa: reference.pullback_tangent(&elastic_tangent),
            },
            partial,
        )?;
        state.law = Some(law);
        Ok(state)
    }

    /// Shared chain rule for a native physical potential whose coefficients
    /// and internal variables are fixed. Not an admission seam for a finite-
    /// interval incremental potential, which is not physical storage.
    pub(crate) fn from_reference_potential(
        total_deformation_gradient: &[f64; 9],
        temperature_k: f64,
        reference_jet: ThermalReferenceJet,
        potential: ReferencePotentialDerivatives,
    ) -> Result<Self, MaterialError> {
        Self::with_reference_potential(
            total_deformation_gradient,
            temperature_k,
            reference_jet,
            potential,
            ElasticTemperatureDerivatives::default(),
        )
    }

    #[allow(clippy::too_many_lines)]
    fn with_reference_potential(
        total_deformation_gradient: &[f64; 9],
        temperature_k: f64,
        reference_jet: ThermalReferenceJet,
        potential: ReferencePotentialDerivatives,
        partial: ElasticTemperatureDerivatives,
    ) -> Result<Self, MaterialError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(invalid("positive finite absolute temperature required"));
        }
        let jacobian = det3(total_deformation_gradient);
        if total_deformation_gradient
            .iter()
            .any(|value| !value.is_finite())
            || !jacobian.is_finite()
            || jacobian <= 0.0
        {
            return Err(invalid(
                "finite orientation-preserving total deformation required",
            ));
        }
        let ReferencePotentialDerivatives {
            free_energy_j_m3,
            first_piola_pa,
            isothermal_tangent_pa,
        } = potential;
        let reference = reference_jet.reference;
        let elastic_deformation_gradient =
            reference.elastic_deformation(total_deformation_gradient);
        let jg = reference.volume_ratio();

        // Express the natural-map path as an equivalent deformation path at
        // the CURRENT fixed G, so original Psi, P and A can be used directly.
        // D=H G_T, C=H G_TT; v=-F D, w=F(2D²-C), Jg_T/Jg=tr(D).
        // This never changes the native history binding or pushes it to a
        // different frame. In particular, H and G_T need not commute.
        let h = reference.elastic_deformation(&IDENTITY);
        let d = product(&h, &reference_jet.first_derivative_per_k);
        let c = product(&h, &reference_jet.second_derivative_per_k2);
        let d_squared = product(&d, &d);
        let trace_d = trace(&d);
        let relative_volume_tt = trace_d * trace_d - trace(&d_squared) + trace(&c);
        let equivalent_t = product(total_deformation_gradient, &d).map(|value| -value);
        let curvature = core::array::from_fn(|ij| 2.0 * d_squared[ij] - c[ij]);
        let equivalent_tt = product(total_deformation_gradient, &curvature);
        let piola_t: [f64; 9] =
            core::array::from_fn(|row| dot(&isothermal_tangent_pa[row], &equivalent_t));
        let work_t = dot(&first_piola_pa, &equivalent_t);
        let mut psi_t = trace_d * free_energy_j_m3 + work_t;
        let mut psi_tt = relative_volume_tt * free_energy_j_m3
            + 2.0 * trace_d * work_t
            + dot(&equivalent_t, &piola_t)
            + dot(&first_piola_pa, &equivalent_tt);
        let variable_coefficients = partial.energy_j_m3_k != 0.0
            || partial.energy_j_m3_k2 != 0.0
            || partial.piola_pa_k.iter().any(|value| *value != 0.0);
        if variable_coefficients {
            // Psi=Jg W(Fe,T): include the partial thermal response and both
            // mixed reference/coefficient terms, without re-differencing G.
            psi_t += jg * partial.energy_j_m3_k;
            psi_tt += 2.0 * trace_d * jg * partial.energy_j_m3_k
                + 2.0
                    * dot(
                        &reference.pullback_piola(&partial.piola_pa_k),
                        &equivalent_t,
                    )
                + jg * partial.energy_j_m3_k2;
        }
        let entropy_j_m3_k = -psi_t;
        let internal_energy_j_m3 = temperature_k.mul_add(entropy_j_m3_k, free_energy_j_m3);
        let mechanical_heat_capacity_j_m3_k = -temperature_k * psi_tt;

        // d(Jg Pe H^T)/dT = tr(D) P + A:v - P D^T.
        // The complete original tangent carries every shear/cross component.
        let mut stress_temperature_derivative_pa_k = core::array::from_fn(|ij| {
            let (i, j) = (ij / 3, ij % 3);
            let inverse_work = (0..3).fold(0.0, |sum, k| {
                first_piola_pa[3 * i + k].mul_add(d[3 * j + k], sum)
            });
            trace_d.mul_add(first_piola_pa[ij], piola_t[ij] - inverse_work)
        });
        if variable_coefficients {
            let coefficient_rate = reference.pullback_piola(&partial.piola_pa_k);
            for (derivative, increment) in stress_temperature_derivative_pa_k
                .iter_mut()
                .zip(coefficient_rate)
            {
                *derivative += increment;
            }
        }
        if [
            trace_d,
            relative_volume_tt,
            free_energy_j_m3,
            entropy_j_m3_k,
            internal_energy_j_m3,
            mechanical_heat_capacity_j_m3_k,
        ]
        .iter()
        .chain(first_piola_pa.iter())
        .chain(isothermal_tangent_pa.iter().flatten())
        .chain(stress_temperature_derivative_pa_k.iter())
        .any(|value| !value.is_finite())
        {
            return Err(unrepresentable());
        }
        Ok(Self {
            law: None,
            deformation_gradient: *total_deformation_gradient,
            elastic_deformation_gradient,
            temperature_k,
            reference_jet,
            free_energy_j_m3,
            entropy_j_m3_k,
            internal_energy_j_m3,
            mechanical_heat_capacity_j_m3_k,
            first_piola_pa,
            isothermal_tangent_pa,
            stress_temperature_derivative_pa_k,
        })
    }

    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub const fn reference_jet(&self) -> &ThermalReferenceJet {
        &self.reference_jet
    }
    pub const fn deformation_gradient(&self) -> &[f64; 9] {
        &self.deformation_gradient
    }
    pub const fn elastic_deformation_gradient(&self) -> &[f64; 9] {
        &self.elastic_deformation_gradient
    }
    /// Mechanical Helmholtz potential per ORIGINAL volume [J/m³].
    pub const fn free_energy_j_m3(&self) -> f64 {
        self.free_energy_j_m3
    }
    /// Mechanical entropy contribution per ORIGINAL volume [J/(m³ K)].
    pub const fn entropy_j_m3_k(&self) -> f64 {
        self.entropy_j_m3_k
    }
    /// Mechanical internal energy Psi+T*S per ORIGINAL volume [J/m³].
    pub const fn internal_energy_j_m3(&self) -> f64 {
        self.internal_energy_j_m3
    }
    /// Direct change of the original-volume mechanical Helmholtz potential.
    /// Native source-amplitude and deformation increments, including changing
    /// natural volume, are retained independently of large absolute storage.
    /// The caller owns the physical source/trajectory relationship of the two
    /// states; this operation does not label their difference as work or heat.
    ///
    /// Both endpoints must retain the native elastic increment owner. Fixed-
    /// history viscoelastic snapshots explicitly refuse: their separate stable
    /// Cv/reference increment owner is not supplied by this thermal snapshot.
    pub fn free_energy_increment_to(&self, after: &Self) -> Result<f64, MaterialError> {
        let (Some(before_law), Some(after_law)) = (&self.law, &after.law) else {
            return Err(invalid(
                "stable thermal energy increments require native elastic endpoints; fixed-history viscoelastic increments are unsupported",
            ));
        };
        let before_f = &self.elastic_deformation_gradient;
        let after_f = &after.elastic_deformation_gradient;
        let increment = self.reference_jet.reference.original_density_increment_to(
            &after.reference_jet.reference,
            before_law.energy(before_f),
            after_law.energy(after_f),
            before_law.energy_increment_to(before_f, after_law, after_f),
        );
        if increment.is_finite() {
            Ok(increment)
        } else {
            Err(unrepresentable())
        }
    }
    /// Mechanical ΔU=ΔPsi+Δ(T*S) on the original endpoint potentials.
    /// The symmetric product identity avoids subtracting two large rounded
    /// internal energies. Entropy values retain their original source-derived
    /// numerical precision; this is not a new constitutive derivative model.
    /// As for `free_energy_increment_to`, fixed-history viscoelastic endpoints
    /// refuse until their own direct stable energy increment is supplied.
    pub fn internal_energy_increment_to(&self, after: &Self) -> Result<f64, MaterialError> {
        let before_t = Dd::from_f64(self.temperature_k);
        let after_t = Dd::from_f64(after.temperature_k);
        let before_s = Dd::from_f64(self.entropy_j_m3_k);
        let after_s = Dd::from_f64(after.entropy_j_m3_k);
        let half = Dd::from_f64(0.5);
        let change = (Dd::from_f64(self.free_energy_increment_to(after)?)
            + (half * before_t + half * after_t) * (after_s - before_s)
            + (half * before_s + half * after_s) * (after_t - before_t))
            .to_f64();
        if change.is_finite() {
            Ok(change)
        } else {
            Err(unrepresentable())
        }
    }
    /// Mechanical addition to fixed-total-F capacity [J/(m³ K)], which may
    /// be negative. Positivity belongs to the complete caloric+mechanical sum.
    pub const fn mechanical_heat_capacity_j_m3_k(&self) -> f64 {
        self.mechanical_heat_capacity_j_m3_k
    }
    pub const fn first_piola_pa(&self) -> &[f64; 9] {
        &self.first_piola_pa
    }
    pub const fn isothermal_tangent_pa(&self) -> &[[f64; 9]; 9] {
        &self.isothermal_tangent_pa
    }
    /// Complete mixed derivative dP/dT at fixed TOTAL F [Pa/K].
    pub const fn stress_temperature_derivative_pa_k(&self) -> &[f64; 9] {
        &self.stress_temperature_derivative_pa_k
    }

    /// Differentiate this physical potential along prescribed total Fdot
    /// and Tdot, using its original A=dP/dF, beta=dP/dT and heat capacity.
    ///
    /// Pdot=A:Fdot+beta*Tdot and Sdot=-beta:Fdot+(c_F/T)*Tdot retain both
    /// temperature-dependent coefficients and the full noncommuting G(T)
    /// chain rule. Fedot removes Gdot=G_T*Tdot before applying G^-1, so
    /// thermal expansion is not reported as elastic strain rate.
    ///
    /// This is a local derivative, not a finite step or a source-band proof.
    /// Caloric storage, irreversible heating and internal-variable rates are
    /// absent. In particular a fixed-history viscoelastic snapshot remains a
    /// fixed-history derivative even when Tdot is nonzero.
    pub fn rate(
        &self,
        deformation_gradient_rate_per_s: [f64; 9],
        temperature_rate_k_s: f64,
    ) -> Result<FiniteThermoelasticRate, MaterialError> {
        if !temperature_rate_k_s.is_finite()
            || deformation_gradient_rate_per_s
                .iter()
                .any(|value| !value.is_finite())
        {
            return Err(invalid("finite deformation and temperature rates required"));
        }
        let elastic_deformation_gradient_rate_per_s = elastic_deformation_rate(
            &self.elastic_deformation_gradient,
            &self.reference_jet,
            &deformation_gradient_rate_per_s,
            temperature_rate_k_s,
        );
        let first_piola_rate_pa_s = core::array::from_fn(|i| {
            self.stress_temperature_derivative_pa_k[i].mul_add(
                temperature_rate_k_s,
                dot(
                    &self.isothermal_tangent_pa[i],
                    &deformation_gradient_rate_per_s,
                ),
            )
        });
        let stress_power = dot(&self.first_piola_pa, &deformation_gradient_rate_per_s);
        let entropy_rate_w_m3_k = (self.mechanical_heat_capacity_j_m3_k / self.temperature_k)
            .mul_add(
                temperature_rate_k_s,
                -dot(
                    &self.stress_temperature_derivative_pa_k,
                    &deformation_gradient_rate_per_s,
                ),
            );
        let free_energy_rate_w_m3 =
            (-self.entropy_j_m3_k).mul_add(temperature_rate_k_s, stress_power);
        let internal_energy_rate_w_m3 = self
            .temperature_k
            .mul_add(entropy_rate_w_m3_k, stress_power);
        if [
            free_energy_rate_w_m3,
            entropy_rate_w_m3_k,
            internal_energy_rate_w_m3,
        ]
        .iter()
        .chain(&elastic_deformation_gradient_rate_per_s)
        .chain(&first_piola_rate_pa_s)
        .any(|value| !value.is_finite())
        {
            return Err(unrepresentable());
        }
        Ok(FiniteThermoelasticRate {
            deformation_gradient_rate_per_s,
            elastic_deformation_gradient_rate_per_s,
            temperature_rate_k_s,
            first_piola_rate_pa_s,
            free_energy_rate_w_m3,
            entropy_rate_w_m3_k,
            internal_energy_rate_w_m3,
        })
    }

    /// Freeze the local reversible tangent at constant total entropy.
    ///
    /// The supplied capacity is the caloric dU/dT per ORIGINAL volume and
    /// excludes this state's mechanical contribution. A source owner must
    /// establish its thermodynamic reference and a smooth local neighborhood;
    /// an arbitrary tabulated Cp or Cv is not interchangeable with this input.
    /// The original inertial density supplies rho_current=rho0/det(F).
    pub fn adiabatic_acoustics(
        &self,
        reference_density_kg_m3: f64,
        caloric_heat_capacity_j_m3_k: f64,
    ) -> Result<AdiabaticHyperelasticState, MaterialError> {
        if !reference_density_kg_m3.is_finite()
            || reference_density_kg_m3 <= 0.0
            || !caloric_heat_capacity_j_m3_k.is_finite()
            || caloric_heat_capacity_j_m3_k <= 0.0
        {
            return Err(invalid(
                "positive finite reference density and caloric capacity required",
            ));
        }
        let fixed_deformation_heat_capacity_j_m3_k =
            caloric_heat_capacity_j_m3_k + self.mechanical_heat_capacity_j_m3_k;
        if !fixed_deformation_heat_capacity_j_m3_k.is_finite()
            || fixed_deformation_heat_capacity_j_m3_k <= 0.0
        {
            return Err(invalid(
                "total heat capacity at fixed total deformation must be positive",
            ));
        }
        let beta = self.stress_temperature_derivative_pa_k;
        let mut tangent = self.isothermal_tangent_pa;
        let mut temperature_gradient_k = [0.0; 9];
        if beta.iter().any(|value| *value != 0.0) {
            let scale = self.temperature_k / fixed_deformation_heat_capacity_j_m3_k;
            if !scale.is_finite() || scale <= 0.0 {
                return Err(unrepresentable());
            }
            for i in 0..9 {
                temperature_gradient_k[i] = scale * beta[i];
                if !temperature_gradient_k[i].is_finite()
                    || (beta[i] != 0.0 && temperature_gradient_k[i] == 0.0)
                {
                    return Err(unrepresentable());
                }
                for j in 0..=i {
                    let correction = temperature_gradient_k[i] * beta[j];
                    if !correction.is_finite()
                        || (temperature_gradient_k[i] != 0.0 && beta[j] != 0.0 && correction == 0.0)
                    {
                        return Err(unrepresentable());
                    }
                    // A single symmetric outer correction leaves every
                    // original tangent entry unchanged apart from that term.
                    tangent[i][j] += correction;
                    if i != j {
                        tangent[j][i] += correction;
                    }
                }
            }
        }
        let acoustic = HyperelasticAcousticState::from_tangent(
            &self.deformation_gradient,
            reference_density_kg_m3,
            tangent,
        )?;
        Ok(AdiabaticHyperelasticState {
            source: self.clone(),
            caloric_heat_capacity_j_m3_k,
            fixed_deformation_heat_capacity_j_m3_k,
            temperature_gradient_k,
            acoustic,
        })
    }
}

/// Native finite-deformation acoustic snapshot with explicit thermodynamic
/// inputs. Directional instabilities remain signed in the shared wave owner;
/// positive local heat capacity is not a global ellipticity or source-band proof.
#[derive(Clone, Debug, PartialEq)]
pub struct AdiabaticHyperelasticState {
    source: FiniteThermoelasticState,
    caloric_heat_capacity_j_m3_k: f64,
    fixed_deformation_heat_capacity_j_m3_k: f64,
    temperature_gradient_k: [f64; 9],
    acoustic: HyperelasticAcousticState,
}

impl AdiabaticHyperelasticState {
    pub const fn source(&self) -> &FiniteThermoelasticState {
        &self.source
    }
    pub const fn caloric_heat_capacity_j_m3_k(&self) -> f64 {
        self.caloric_heat_capacity_j_m3_k
    }
    pub const fn fixed_deformation_heat_capacity_j_m3_k(&self) -> f64 {
        self.fixed_deformation_heat_capacity_j_m3_k
    }
    /// Infinitesimal isentropic dT/dF. Contract with a TOTAL deformation
    /// increment; this gradient does not advance a finite temperature path.
    pub const fn temperature_gradient_k(&self) -> &[f64; 9] {
        &self.temperature_gradient_k
    }
    pub const fn isentropic_tangent_pa(&self) -> &[[f64; 9]; 9] {
        self.acoustic.tangent_pa()
    }
    pub const fn acoustic(&self) -> &HyperelasticAcousticState {
        &self.acoustic
    }
}
