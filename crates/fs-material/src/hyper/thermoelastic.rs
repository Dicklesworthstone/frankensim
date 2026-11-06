//! Finite thermoelastic derivatives of an explicitly prescribed natural map.
//!
//! With H=G(T)^-1, Fe=F H and Jg=det(G), the original-volume Helmholtz
//! potential is Psi(F,T)=Jg We(Fe). Native energy, Piola stress and Hessian
//! supply its exact chain rule. All elastic coefficients are held fixed;
//! the declared G, dG/dT and d2G/dT2 own the temperature dependence.
//!
//! Mechanical entropy is -Psi_T, internal energy is Psi-T Psi_T, and the
//! mechanical contribution to fixed-F heat capacity is -T Psi_TT. These
//! are additions to a separately admitted caloric reference, not a complete
//! caloric law, inferred heat, or a temperature/history evolution algorithm.
//!
//! For beta=dP/dT and positive total reference-volume capacity c_F, eliminating
//! temperature from ds=(c_F/T)dT-beta:dF at ds=0 gives
//! A_S=A_T+(T/c_F) beta (x) beta. This is a full compressible 3-D local
//! tangent. Incompressible/mixed pressure, plane stress and irreversible
//! histories require their own constrained thermodynamic owners. Source
//! smoothness and the actual meaning of the caloric capacity remain explicit.

use super::acoustic::HyperelasticAcousticState;
use super::reference::ReferenceDeformation;
use super::{Hyperelastic, det3};
use crate::MaterialError;

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
/// Preparation admits the elastic domain, including the original native
/// constitutive switching conventions. It does not certify a two-sided
/// tangent at a fibre activation surface or an arbitrary source kink.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteThermoelasticState {
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

impl FiniteThermoelasticState {
    /// Evaluate the original full elastic energy and its exact thermal chain
    /// rule. The elastic stretch envelope applies to Fe=F G^-1; total F has
    /// its own orientation/representability gate. Caloric capacity is not
    /// needed to observe mechanical entropy/internal energy.
    #[allow(clippy::too_many_lines)]
    pub fn try_new(
        law: &Hyperelastic,
        total_deformation_gradient: &[f64; 9],
        temperature_k: f64,
        reference_jet: ThermalReferenceJet,
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
        let jg = reference.volume_ratio();
        let free_energy_j_m3 = jg * elastic_energy;
        let first_piola_pa = reference.pullback_piola(&elastic_piola);
        let isothermal_tangent_pa = reference.pullback_tangent(&elastic_tangent);

        // Reuse the prepared inverse through its existing kinematic owner.
        // B=G_T H, C=G_TT H; Fe_T=-Fe B and Fe_TT=Fe(2B²-C).
        // Cyclic trace invariance gives Jg_T/Jg=tr(B).
        let h = reference.elastic_deformation(&IDENTITY);
        let b = product(&reference_jet.first_derivative_per_k, &h);
        let c = product(&reference_jet.second_derivative_per_k2, &h);
        let b_squared = product(&b, &b);
        let trace_b = trace(&b);
        let jg_t = jg * trace_b;
        let jg_tt = jg * (trace_b * trace_b - trace(&b_squared) + trace(&c));
        let elastic_t = product(&elastic_deformation_gradient, &b).map(|value| -value);
        let curvature = core::array::from_fn(|ij| 2.0 * b_squared[ij] - c[ij]);
        let elastic_tt = product(&elastic_deformation_gradient, &curvature);
        let elastic_piola_t: [f64; 9] =
            core::array::from_fn(|row| dot(&elastic_tangent[row], &elastic_t));
        let elastic_work_t = dot(&elastic_piola, &elastic_t);
        let psi_t = jg_t * elastic_energy + jg * elastic_work_t;
        let psi_tt = jg_tt * elastic_energy
            + 2.0 * jg_t * elastic_work_t
            + jg * (dot(&elastic_t, &elastic_piola_t) + dot(&elastic_piola, &elastic_tt));
        let entropy_j_m3_k = -psi_t;
        let internal_energy_j_m3 = temperature_k.mul_add(entropy_j_m3_k, free_energy_j_m3);
        let mechanical_heat_capacity_j_m3_k = -temperature_k * psi_tt;

        // d(Jg Pe H^T)/dT = tr(B) P + Jg (Pe_T - Pe B^T) H^T.
        // The same native pullback carries every shear/cross component.
        let pulled_rate = reference.pullback_piola(&core::array::from_fn(|ij| {
            let (i, j) = (ij / 3, ij % 3);
            let inverse_work = (0..3).fold(0.0, |sum, k| {
                elastic_piola[3 * i + k].mul_add(b[3 * j + k], sum)
            });
            elastic_piola_t[ij] - inverse_work
        }));
        let stress_temperature_derivative_pa_k =
            core::array::from_fn(|ij| trace_b.mul_add(first_piola_pa[ij], pulled_rate[ij]));
        if [
            elastic_energy,
            jg_t,
            jg_tt,
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
