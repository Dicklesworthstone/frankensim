//! Native finite-deformation ductility using Kirchhoff stress, tau = P F^T.
//!
//! A first Piola tensor is not a symmetric spatial stress. Its raw components
//! must never drive triaxiality. The positive Jacobian factor between Kirchhoff
//! and Cauchy stress cancels from mean stress / equivalent stress.

/// Pointwise initiation on native finite plastic and thermoplastic sheets.
pub mod sheet;

use fs_material::plastic::finite::FiniteJ2Material;
use fs_material::plastic::finite::thermal::FiniteThermalJ2Material;

use super::*;
use crate::hyper3::HyperTetSettings;
use crate::hyper3::plastic::{
    FinitePlasticStepControls, FinitePlasticTetError, FinitePlasticTetSolid, FinitePlasticTetStep,
};
use crate::hyper3::plastic::dynamics::{
    FinitePlasticMotionControls, FinitePlasticMotionStep, FinitePlasticTetMotion,
};
use crate::hyper3::plastic::hex::{
    FiniteThermoplasticHexError, FiniteThermoplasticHexSolid, FiniteThermoplasticHexStep,
};
use crate::hyper3::plastic::thermal::{
    FiniteThermoplasticControls, FiniteThermoplasticTetError, FiniteThermoplasticTetSolid,
    FiniteThermoplasticTetStep,
};
use crate::IncrementSettings;

fn check_thermal_families<M: FiniteThermalJ2Material>(
    families: &[M], cx: &Cx<'_>,
) -> Result<(), DuctileSolidError> {
    checkpoint(cx)?;
    for (element, family) in families.iter().enumerate() {
        if element % 256 == 0 { checkpoint(cx)?; }
        if family.thermal_finite_j2().is_norton_creep() {
            return Err(DuctileSolidError::Input(
                "pure Norton creep history is not plastic strain; a separate creep-rupture calibration is required",
            ));
        }
    }
    Ok(())
}

fn kirchhoff(p: &[f64; 9], f: &[f64; 9]) -> Result<Voigt, DuctileSolidError> {
    let mut tau = [0.0; 9];
    for i in 0..3 {
        for j in 0..3 {
            tau[3 * i + j] = (0..3).map(|k| p[3 * i + k] * f[3 * j + k]).sum();
        }
    }
    if tau.iter().any(|x| !x.is_finite()) {
        return Err(DuctileSolidError::Input("Kirchhoff stress is not representable"));
    }
    let scale = tau.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    if scale > 0.0 {
        for (i, j) in [(1, 3), (5, 7), (6, 2)] {
            if (tau[i] / scale - tau[j] / scale).abs() > 4096.0 * f64::EPSILON {
                return Err(DuctileSolidError::Input("native Piola/deformation pair gives nonsymmetric spatial stress"));
            }
        }
    }
    Ok([
        tau[0], tau[4], tau[8],
        0.5 * tau[1] + 0.5 * tau[3],
        0.5 * tau[5] + 0.5 * tau[7],
        0.5 * tau[6] + 0.5 * tau[2],
    ])
}

impl<'a, L: DuctileInitiationLaw> DuctileHistory<'a, L> {
    fn trial_finite(
        &self, states: &[FiniteJ2State], piola_pa: &[[f64; 9]], deformation: &[[f64; 9]],
        conditions: DuctileStepConditions<'_>, cx: &Cx<'_>,
    ) -> Result<Self, DuctileSolidError> {
        if states.len() != self.points.len()
            || piola_pa.len() != states.len() || deformation.len() != states.len()
        {
            return Err(DuctileSolidError::Input("native finite-strain point arrays do not match ductile history"));
        }
        self.trial_stress(states, |point| kirchhoff(&piola_pa[point], &deformation[point]),
            conditions, cx)
    }

    /// One original finite-J2 tetrahedral equilibrium increment. Reconstruct
    /// the spatial stress from the SAME returned Piola/deformation pair.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_finite_tet<'m, M: FiniteJ2Material>(
        &self, solid: &FinitePlasticTetSolid<'m, M>,
        prescribed_dofs_m: &[(usize, f64)], nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3], conditions: DuctileStepConditions<'_>,
        settings: HyperTetSettings, controls: FinitePlasticStepControls, cx: &Cx<'_>,
    ) -> Result<DuctileAccepted<'a, FinitePlasticTetStep<'m, M>, L>, DuctileAdvanceError<FinitePlasticTetError>> {
        self.prepare(solid.element_states(), conditions, cx)?;
        let mechanics = solid.advance(prescribed_dofs_m, nodal_forces_n, body_force_n_m3,
            settings, controls, cx).map_err(DuctileAdvanceError::Mechanics)?;
        let next = &mechanics.solid;
        let ductility = self.trial_finite(next.element_states(), next.element_first_piola_pa(),
            next.element_deformation(), conditions, cx)?;
        Ok(DuctileAccepted { mechanics, ductility })
    }

    /// Native finite-J2 backward-Euler motion. The same physical duration
    /// drives inertia and the interval-average equivalent plastic strain rate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_finite_tet_motion<'m, M: FiniteJ2Material>(
        &self, motion: &FinitePlasticTetMotion<'m, M>,
        prescribed_dofs_m: &[(usize, f64)], nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3], conditions: DuctileStepConditions<'_>,
        settings: HyperTetSettings, controls: FinitePlasticMotionControls, cx: &Cx<'_>,
    ) -> Result<DuctileAccepted<'a, FinitePlasticMotionStep<'m, M>, L>, DuctileAdvanceError<FinitePlasticTetError>> {
        self.prepare(motion.solid().element_states(), conditions, cx)?;
        let mechanics = motion.advance(prescribed_dofs_m, nodal_forces_n, body_force_n_m3,
            conditions.duration_s, settings, controls, cx).map_err(DuctileAdvanceError::Mechanics)?;
        let next = mechanics.motion.solid();
        let ductility = self.trial_finite(next.element_states(), next.element_first_piola_pa(),
            next.element_deformation(), conditions, cx)?;
        Ok(DuctileAccepted { mechanics, ductility })
    }

    /// Native finite thermoplastic tetrahedra, including separately accumulated
    /// coupled creep. Native alpha is plastic; creep_alpha never consumes this
    /// initiation allowance. Pure-creep families refuse because their alpha is
    /// creep accumulation. Actual endpoint material temperatures are used.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_thermoplastic_tet<'m, M: FiniteThermalJ2Material>(
        &self, solid: &FiniteThermoplasticTetSolid<'m, M>, temperatures_k: &[f64],
        prescribed_dofs_m: &[(usize, f64)], nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3], duration_s: f64, settings: HyperTetSettings,
        controls: FiniteThermoplasticControls, cx: &Cx<'_>,
    ) -> Result<DuctileAccepted<'a, FiniteThermoplasticTetStep<'m, M>, L>, DuctileAdvanceError<FiniteThermoplasticTetError>> {
        check_thermal_families(solid.families(), cx)?;
        let conditions = DuctileStepConditions {
            duration_s, temperatures_k: DuctileTemperatures::PerPoint(temperatures_k),
        };
        self.prepare(solid.element_states(), conditions, cx)?;
        let mechanics = solid.advance(temperatures_k, prescribed_dofs_m, nodal_forces_n,
            body_force_n_m3, duration_s, settings, controls, cx)
            .map_err(DuctileAdvanceError::Mechanics)?;
        let next = &mechanics.solid;
        let ductility = self.trial_finite(next.element_states(), next.element_first_piola_pa(),
            next.element_deformation(), DuctileStepConditions {
                duration_s, temperatures_k: DuctileTemperatures::PerPoint(next.temperatures_k()),
            }, cx)?;
        Ok(DuctileAccepted { mechanics, ductility })
    }

    /// Native fully integrated/B-bar thermoplastic hexes. Temperature is supplied
    /// per element to the original solver and evaluated at ALL eight Gauss points.
    /// Use the assumed deformation matching stored Piola stress, not a mismatched
    /// compatible deformation when the native B-bar formulation is active.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_thermoplastic_hex<'m, M: FiniteThermalJ2Material>(
        &self, solid: &FiniteThermoplasticHexSolid<'m, M>, temperatures_k: &[f64],
        prescribed_dofs_m: &[(usize, f64)], nodal_forces_n: &[(usize, f64)],
        duration_s: f64, settings: IncrementSettings,
        controls: FinitePlasticStepControls, cx: &Cx<'_>,
    ) -> Result<DuctileAccepted<'a, FiniteThermoplasticHexStep<'m, M>, L>, DuctileAdvanceError<FiniteThermoplasticHexError>> {
        check_thermal_families(solid.families(), cx)?;
        let conditions = DuctileStepConditions {
            duration_s, temperatures_k: DuctileTemperatures::PerHexElement(temperatures_k),
        };
        self.prepare(solid.quadrature_states(), conditions, cx)?;
        let mechanics = solid.advance(temperatures_k, prescribed_dofs_m, nodal_forces_n,
            duration_s, settings, controls, cx).map_err(DuctileAdvanceError::Mechanics)?;
        let next = &mechanics.solid;
        let ductility = self.trial_finite(next.quadrature_states(), next.quadrature_first_piola_pa(),
            next.quadrature_assumed_deformation(), conditions, cx)?;
        Ok(DuctileAccepted { mechanics, ductility })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kirchhoff_uses_rotated_spatial_stress_and_rejects_skew() {
        let f = [0.0, -1.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0];
        let p = [0.0, -9.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 3.0];
        assert_eq!(kirchhoff(&p, &f).unwrap(), [9.0, 3.0, 3.0, 0.0, 0.0, 0.0]);
        let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        assert!(kirchhoff(&[1.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], &identity).is_err());
    }
}
