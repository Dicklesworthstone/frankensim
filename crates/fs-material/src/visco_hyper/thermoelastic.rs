//! Physical finite-strain thermodynamics at fixed accepted viscous memory.
//!
//! For the original law, Psi(F,T,Cv)=det(G(T))*W(F G(T)^-1,Cv).
//! The native frozen potential at ZERO elapsed time is this physical stored
//! energy: all branch metrics stay fixed and no relaxation takes place. Its
//! original stress and exact tangent supply the shared thermal chain rule.
//! A positive-duration incremental potential is never interpreted as Psi.
//!
//! Equilibrium and branch energy coefficients are held fixed. Temperature may
//! change the natural map and, independently, positive relaxation times, but
//! the latter do not enter reversible derivatives at fixed Cv. This snapshot
//! does not evolve temperature/history, supply a complete caloric law, or
//! convert reference work into heat. Exact incompressibility and remainder
//! histories retain their separate pressure/constraint thermodynamic owners.

use super::{Hyperelastic, MAX_FINITE_STRAIN_BRANCHES, ViscoHyperelastic, ViscoState, parameters};
use crate::MaterialError;
use crate::hyper::cards::FiniteStrainPotential;
use crate::hyper::mixed::VolumetricTreatment;
use crate::hyper::plane_stress::{PlaneStressSettings, plane_stress_state_with_seed};
use crate::hyper::thermoelastic::plane_stress::FiniteThermoelasticPlaneStressState;
use crate::hyper::thermoelastic::{
    FiniteThermoelasticState, ReferencePotentialDerivatives, ThermalReferenceJet,
};

impl ViscoHyperelastic {
    /// Reversible full-compressible thermodynamics at total F and absolute T,
    /// holding every accepted inverse viscous metric fixed in its original
    /// natural-coordinate frame. The deformation may be a trial different
    /// from the accepted endpoint; the complete accepted history is borrowed
    /// unchanged and no physical interval is advanced.
    ///
    /// The jet's current G must equal this law's reference exactly, and the
    /// accepted history must retain that native binding. A prescribed finite
    /// change of G must first use the existing reference-change transaction;
    /// this derivative query never silently rebinds or resets its memory.
    ///
    /// The result contains mechanical Helmholtz energy, entropy, internal
    /// energy, fixed-F heat-capacity contribution and dP/dT per ORIGINAL
    /// volume. It also retains the unchanged instantaneous native P and dP/dF.
    /// The heat-capacity contribution may be negative and must be combined
    /// with the separately admitted caloric capacity. No dissipative heat or
    /// finite-frequency viscoelastic acoustic response is claimed.
    ///
    /// Native parameter, history, elastic-domain and representability gates
    /// apply. As for ordinary native stepping, accepted public fields must be
    /// unmodified outputs of the original owner; this query cannot certify an
    /// externally edited or imported physical ledger.
    pub fn frozen_thermoelastic_state(
        &self,
        accepted: &ViscoState,
        total_deformation_gradient: &[f64; 9],
        temperature_k: f64,
        reference_jet: ThermalReferenceJet,
    ) -> Result<FiniteThermoelasticState, MaterialError> {
        self.validate_frozen_thermoelastic(temperature_k, &reference_jet)?;
        let frozen = self.freeze(accepted, 0.0)?;
        let (first_piola_pa, free_energy_j_m3) =
            frozen.stress_and_potential(total_deformation_gradient)?;
        let isothermal_tangent_pa = frozen.consistent_tangent(total_deformation_gradient)?;
        FiniteThermoelasticState::from_reference_potential(
            total_deformation_gradient,
            temperature_k,
            reference_jet,
            ReferencePotentialDerivatives {
                free_energy_j_m3,
                first_piola_pa,
                isothermal_tangent_pa,
            },
        )
    }

    /// Relax the original finite-compliance transverse tractions at fixed
    /// in-plane F, T and accepted Cv, then condense the same thermal potential.
    /// This is a reversible instantaneous material query, not a proposed time
    /// step: all trial states borrow the same metrics and zero elapsed time.
    ///
    /// The reduced heat-capacity contribution and dP_in_plane/dT include the
    /// native transverse response to temperature. They differ from holding
    /// the full F fixed. No exact-volume pressure or remainder state is
    /// admitted by this full-energy entry.
    ///
    /// A supplied transverse column seeds the original native solve. `None`
    /// retains the accepted shears and uses det(G)/det(F_in_plane) for the
    /// initial thickness; the native solver owns fallback and convergence.
    pub fn frozen_thermoelastic_plane_stress_state(
        &self,
        accepted: &ViscoState,
        in_plane: &[f64; 4],
        temperature_k: f64,
        reference_jet: ThermalReferenceJet,
        initial_transverse: Option<[f64; 3]>,
        settings: PlaneStressSettings,
    ) -> Result<FiniteThermoelasticPlaneStressState, MaterialError> {
        self.validate_frozen_thermoelastic(temperature_k, &reference_jet)?;
        let frozen = self.freeze(accepted, 0.0)?;
        let seed = initial_transverse.unwrap_or_else(|| {
            let determinant = in_plane[0] * in_plane[3] - in_plane[1] * in_plane[2];
            [
                accepted.deformation_gradient[2],
                accepted.deformation_gradient[5],
                self.reference_deformation().volume_ratio() / determinant,
            ]
        });
        let plane = plane_stress_state_with_seed(&frozen, in_plane, seed, settings)?;
        let full = self.frozen_thermoelastic_state(
            accepted,
            &plane.deformation_gradient,
            temperature_k,
            reference_jet,
        )?;
        FiniteThermoelasticPlaneStressState::try_from_state(
            full,
            VolumetricTreatment::OriginalMaterial,
            settings,
        )
    }

    fn validate_frozen_thermoelastic(
        &self,
        temperature_k: f64,
        reference_jet: &ThermalReferenceJet,
    ) -> Result<(), MaterialError> {
        if reference_jet.reference() != self.reference_deformation() {
            return Err(parameters(
                "frozen thermodynamics requires the native current reference deformation".into(),
            ));
        }
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(parameters(
                "frozen thermodynamics requires positive finite absolute temperature".into(),
            ));
        }
        // Public native parameters may have been edited after construction.
        // Reuse their original validators without cloning or rebinding Cv.
        Hyperelastic::new(
            self.equilibrium.model.clone(),
            self.equilibrium.stretch_limit,
        )?;
        if self.branches.len() > MAX_FINITE_STRAIN_BRANCHES {
            return Err(parameters(
                "frozen thermodynamic branch count exceeds its original cap".into(),
            ));
        }
        for branch in &self.branches {
            branch.energy.validate()?;
            if !branch.tau_s.is_finite() || branch.tau_s <= 0.0 {
                return Err(parameters(
                    "frozen thermodynamics requires finite positive relaxation times".into(),
                ));
            }
        }
        Ok(())
    }
}
