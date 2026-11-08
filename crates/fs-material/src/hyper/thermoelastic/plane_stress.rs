//! Reversible thermodynamics on the original finite-compliance plane branch.
//!
//! The native transverse solve eliminates z=[F02,F12,F22] at P_z=0 while
//! a=[F00,F01,F10,F11] and T are held. Its same positive Hessian block gives
//! z_T=-A_zz^-1 beta_z, beta_bar=beta_a-A_az A_zz^-1 beta_z, and
//! Psi_bar,TT=Psi_TT-beta_z^T A_zz^-1 beta_z. Entropy and internal energy
//! keep their full-state values; fixed-a heat capacity acquires the positive
//! contribution T beta_z^T A_zz^-1 beta_z.
//!
//! Every density is per ORIGINAL reference volume and every gradient is total
//! F. A membrane owner multiplies by its ORIGINAL reference thickness and uses
//! the existing Green-strain chain rule; relaxed thickness is not inventory.
//! This is a local reversible snapshot, not a heat source or an evolution law.
//! Exact-volume pressure constraints need a different constrained thermal
//! derivation and are explicitly refused by the snapshot-reduction entry.

use super::{
    FiniteThermoelasticRate, FiniteThermoelasticState, ThermalReferenceJet, invalid,
    unrepresentable,
};
use crate::MaterialError;
use crate::hyper::Hyperelastic;
use crate::hyper::mixed::VolumetricTreatment;
use crate::hyper::plane_stress::{
    PlaneStressKinematicResponse, PlaneStressResponse, PlaneStressSettings, PlaneStressState,
    PositiveBlock, plane_stress_kinematic_response_from, plane_stress_state,
    plane_stress_state_with_seed,
};
use crate::hyper::reference::ReferencedHyperelastic;

const IN_PLANE: [usize; 4] = [0, 1, 3, 4];
const TRANSVERSE: [usize; 3] = [2, 5, 8];

/// Mechanical thermal storage and mixed derivatives on one stable native
/// plane-stress branch. The retained full snapshot supplies original F, G(T),
/// temperature and physical energy; no local material history is advanced.
#[derive(Clone, Debug)]
pub struct FiniteThermoelasticPlaneStressState {
    full_state: FiniteThermoelasticState,
    kinematics: PlaneStressKinematicResponse,
    mechanical_heat_capacity_j_m3_k: f64,
    stress_temperature_derivative_pa_k: [f64; 4],
    transverse_temperature_derivative_per_k: [f64; 3],
}

impl FiniteThermoelasticPlaneStressState {
    /// Relax all three transverse tractions with the ORIGINAL native solver,
    /// then evaluate the full thermal snapshot and its exact Schur reduction.
    /// Elastic coefficients stay fixed; the supplied natural-map jet retains
    /// its original material frame and all first/second thermal derivatives.
    ///
    /// det(G)/det(a) is only an initial thickness guess, never a volume
    /// constraint. Native domain, stable-branch and iteration gates apply.
    pub fn try_new(
        law: &Hyperelastic,
        in_plane: &[f64; 4],
        temperature_k: f64,
        reference_jet: ThermalReferenceJet,
        settings: PlaneStressSettings,
    ) -> Result<Self, MaterialError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(invalid("positive finite absolute temperature required"));
        }
        // Public native parameters may have changed since their constructor.
        let law = Hyperelastic::new(law.model.clone(), law.stretch_limit)?;
        let material = ReferencedHyperelastic::new(&law, *reference_jet.reference());
        let determinant = in_plane[0] * in_plane[3] - in_plane[1] * in_plane[2];
        let thickness = reference_jet.reference().volume_ratio() / determinant;
        let plane = if thickness.is_finite() && thickness > 0.0 {
            plane_stress_state_with_seed(&material, in_plane, [0.0, 0.0, thickness], settings)?
        } else {
            // An unavailable guess does not replace native geometry/domain
            // admission or manufacture a clipped transverse deformation.
            plane_stress_state(&material, in_plane, settings)?
        };
        let full_state = FiniteThermoelasticState::try_new(
            &law,
            &plane.deformation_gradient,
            temperature_k,
            reference_jet,
        )?;
        Self::condense(full_state, plane)
    }

    /// Reduce an already relaxed physical full-state snapshot without another
    /// transverse solve. This also accepts source-coefficient derivatives and
    /// native fixed-history viscoelastic snapshots of the same Helmholtz energy.
    /// An incremental dissipative potential is not a physical thermal snapshot.
    ///
    /// The native chart requires F20=F21=0, positive in-plane orientation and
    /// positive total thickness. Transverse traction is checked using the
    /// current native compliance scale; A_zz must be positive definite. The
    /// response reports zero iterations because this operation performs none.
    /// No field is projected, no pressure constraint is inferred, and the input
    /// state remains the full physical thermodynamic owner.
    pub fn try_from_state(
        full_state: FiniteThermoelasticState,
        volumetric: VolumetricTreatment,
        settings: PlaneStressSettings,
    ) -> Result<Self, MaterialError> {
        if volumetric != VolumetricTreatment::OriginalMaterial {
            return Err(invalid(
                "exact-volume thermoelastic plane stress requires a constrained thermal owner",
            ));
        }
        if !settings.relative_tolerance.is_finite()
            || settings.relative_tolerance <= 0.0
            || settings.relative_tolerance >= 1.0
            || settings.max_iterations == 0
            || settings.max_backtracks == 0
        {
            return Err(invalid("valid native plane-stress settings required"));
        }
        let f = *full_state.deformation_gradient();
        let determinant = f[0] * f[4] - f[1] * f[3];
        if f[6] != 0.0
            || f[7] != 0.0
            || f[8] <= 0.0
            || !determinant.is_finite()
            || determinant <= 0.0
        {
            return Err(invalid(
                "native total-deformation plane chart required without projection",
            ));
        }
        let a = full_state.isothermal_tangent_pa();
        let block = PositiveBlock::new(TRANSVERSE.map(|i| TRANSVERSE.map(|j| a[i][j])))?;
        let gradient_norm = f[0].hypot(f[1]).hypot(f[3]).hypot(f[4]);
        let scale = block.traction_scale()? * gradient_norm.max(1.0);
        if !scale.is_finite() || scale <= 0.0 {
            return Err(unrepresentable());
        }
        let p = full_state.first_piola_pa();
        let residual = (p[2] / scale).hypot(p[5] / scale).hypot(p[8] / scale);
        if !residual.is_finite() || residual > settings.relative_tolerance {
            return Err(invalid(
                "full thermal snapshot is not in transverse equilibrium",
            ));
        }
        let plane = PlaneStressState {
            deformation_gradient: f,
            energy_density_pa: full_state.free_energy_j_m3(),
            piola: IN_PLANE.map(|i| p[i]),
            relative_transverse_residual: residual,
            iterations: 0,
        };
        Self::condense(full_state, plane)
    }

    fn condense(
        full_state: FiniteThermoelasticState,
        plane: PlaneStressState,
    ) -> Result<Self, MaterialError> {
        let a = full_state.isothermal_tangent_pa();
        let beta = full_state.stress_temperature_derivative_pa_k();
        let block = PositiveBlock::new(TRANSVERSE.map(|i| TRANSVERSE.map(|j| a[i][j])))?;
        let extension = block.solve(TRANSVERSE.map(|i| beta[i]))?;
        let curvature_reduction = TRANSVERSE
            .iter()
            .zip(extension)
            .fold(0.0, |sum, (&i, value)| beta[i].mul_add(value, sum));
        let mechanical_heat_capacity_j_m3_k = full_state.temperature_k().mul_add(
            curvature_reduction,
            full_state.mechanical_heat_capacity_j_m3_k(),
        );
        let stress_temperature_derivative_pa_k = IN_PLANE.map(|i| {
            beta[i]
                - TRANSVERSE
                    .iter()
                    .zip(extension)
                    .fold(0.0, |sum, (&j, value)| a[i][j].mul_add(value, sum))
        });
        if !curvature_reduction.is_finite()
            || curvature_reduction < 0.0
            || !mechanical_heat_capacity_j_m3_k.is_finite()
            || stress_temperature_derivative_pa_k
                .iter()
                .any(|x| !x.is_finite())
        {
            return Err(unrepresentable());
        }
        let kinematics = plane_stress_kinematic_response_from(plane, *a)?;
        Ok(Self {
            full_state,
            kinematics,
            mechanical_heat_capacity_j_m3_k,
            stress_temperature_derivative_pa_k,
            transverse_temperature_derivative_per_k: extension.map(|x| -x),
        })
    }

    pub const fn full_state(&self) -> &FiniteThermoelasticState {
        &self.full_state
    }
    /// Native total F, in-plane Piola stress and isothermal Schur tangent.
    pub const fn response(&self) -> &PlaneStressResponse {
        &self.kinematics.response
    }
    /// Prepared transverse sensitivity to in-plane F at fixed temperature
    /// and history. Use `rate` when the temperature is changing.
    pub const fn kinematics(&self) -> &PlaneStressKinematicResponse {
        &self.kinematics
    }
    /// Reversible full material rates on the current plane-stress branch.
    ///
    /// The caller supplies [Fdot00,Fdot01,Fdot10,Fdot11] and Tdot. The same
    /// native Hessian determines zdot=z_a*adot+z_T*Tdot, including transverse
    /// shear and thickness. It therefore enforces
    /// A_zz*zdot+A_za*adot+beta_z*Tdot=0, rather than holding temperature
    /// fixed when lifting an in-plane rate.
    ///
    /// The original thermal owner then supplies full Piola, elastic-F,
    /// entropy and energy rates. In-plane stress rates occupy [0,1,3,4] of
    /// the returned full row-major array. Transverse stress-rate roundoff is
    /// retained, not overwritten with zeros. At exact transverse equilibrium
    /// these give Sdot=-beta_bar:adot+(c_a/T)*Tdot and
    /// Psidot=P_a:adot-S*Tdot; no pressure, caloric law or history evolution
    /// is inferred. Only matrix contractions are performed after preparation.
    pub fn rate(
        &self,
        in_plane_rate_per_s: [f64; 4],
        temperature_rate_k_s: f64,
    ) -> Result<FiniteThermoelasticRate, MaterialError> {
        if !temperature_rate_k_s.is_finite() {
            return Err(invalid("finite temperature rate required"));
        }
        let mut rate = self
            .kinematics
            .deformation_gradient_rate(in_plane_rate_per_s)?;
        for (i, component) in TRANSVERSE.into_iter().enumerate() {
            rate[component] = self.transverse_temperature_derivative_per_k[i]
                .mul_add(temperature_rate_k_s, rate[component]);
        }
        self.full_state.rate(rate, temperature_rate_k_s)
    }
    pub const fn temperature_k(&self) -> f64 {
        self.full_state.temperature_k()
    }
    pub const fn free_energy_j_m3(&self) -> f64 {
        self.full_state.free_energy_j_m3()
    }
    /// Envelope derivative -Psi_bar,T at fixed in-plane total deformation.
    pub const fn entropy_j_m3_k(&self) -> f64 {
        self.full_state.entropy_j_m3_k()
    }
    pub const fn internal_energy_j_m3(&self) -> f64 {
        self.full_state.internal_energy_j_m3()
    }
    /// -T Psi_bar,TT at fixed in-plane total deformation [J/(m³ K)]. This
    /// mechanical contribution may be negative; a caloric owner must supply
    /// its own reference capacity and admit the positive complete sum.
    pub const fn mechanical_heat_capacity_j_m3_k(&self) -> f64 {
        self.mechanical_heat_capacity_j_m3_k
    }
    /// dP_in-plane/dT after transverse relaxation, in native [00,01,10,11]
    /// order [Pa/K]. It is also minus the in-plane gradient of entropy.
    pub const fn stress_temperature_derivative_pa_k(&self) -> &[f64; 4] {
        &self.stress_temperature_derivative_pa_k
    }
    /// d[F02,F12,F22]/dT at fixed [F00,F01,F10,F11] [1/K].
    pub const fn transverse_temperature_derivative_per_k(&self) -> &[f64; 3] {
        &self.transverse_temperature_derivative_per_k
    }
}
