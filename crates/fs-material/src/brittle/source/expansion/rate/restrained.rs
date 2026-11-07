//! Timed thermal loading with physically restrained TOTAL strain.
//!
//! Coefficients transfer first at the original mechanical history. The total
//! strain is then held while the original integrated free strain follows its
//! declared affine path over positive physical time. The native timed law
//! supplies the complete mean stress, intrinsic fracture and bond-viscous work.
//! This is a restrained material experiment, not the free-expansion solve of
//! a deformable body. An affine free-strain path does not imply affine T when
//! alpha(T) varies, and no caloric heat input or entropy is inferred.

use super::super::{BrittleSourceError, ResolvedThermalBrittle, check_work};
use crate::brittle::rate::BrittleRatePathStep;
use crate::brittle::temperature::BrittleTemperatureStep;
use crate::brittle::{BrittleState, finite, positive};

/// Original coefficient preparation and timed mechanical path, with their
/// combined restrained thermal work. The final accepted state is in
/// `path.endpoint.equilibrium.state`; coefficient preparation is intermediate.
#[derive(Clone, Debug, PartialEq)]
pub struct BrittleThermalRateStep<'a> {
    pub coefficient: BrittleTemperatureStep<'a>,
    pub path: BrittleRatePathStep<'a>,
    /// Minus the TOTAL path-mean stress trace times the free-strain change.
    pub eigenstrain_work_j_m3: f64,
    /// Original intrinsic coefficient- and strain-induced fracture combined.
    /// Physical viscous loss remains in `path.endpoint.viscous_dissipated_j_m3`.
    pub fracture_dissipated_j_m3: f64,
    /// Final storage minus storage before either stage, J/m3.
    pub stored_change_j_m3: f64,
    /// Eigenstrain port minus exact native mechanical-path work, J/m3.
    pub eigenstrain_work_defect_j_m3: f64,
    /// Parameter/eigenstrain work minus storage, fracture and viscous loss.
    pub balance_defect_j_m3: f64,
}

impl ResolvedThermalBrittle {
    /// Advance an explicitly restrained thermal experiment over physical time.
    ///
    /// The original total tensor strain is held. Both coefficient and strain
    /// stages preserve the original source, fixed expansion reference and
    /// irreversible fracture owner. Complete mean stress includes viscosity in
    /// the eigenstrain work port; its positive loss is never charged again as
    /// intrinsic fracture. Every refusal leaves the supplied history unchanged.
    ///
    /// This declares an affine FREE-strain trajectory under fixed new
    /// coefficients. It does not infer an endpoint temperature derivative or
    /// replace a solid's simultaneous displacement/eigenstrain solve.
    pub fn transfer_temperature_over_time<'next>(
        &self,
        state: &BrittleState<'_>,
        next: &'next Self,
        dt_s: f64,
        energy_tolerance_j_m3: f64,
    ) -> Result<BrittleThermalRateStep<'next>, BrittleSourceError> {
        positive(dt_s)?;
        positive(energy_tolerance_j_m3)?;
        let total = self.total_strain(state)?;
        let coefficient = self.transfer_coefficients(state, next, energy_tolerance_j_m3)?;
        let path = next.trial_total_over_time_with_mean_stress(&total, &coefficient.state, dt_s)?;
        let delta_free = finite(next.free_linear_strain() - self.free_linear_strain())?;
        let trace =
            finite(path.mean_stress_pa[0] + path.mean_stress_pa[1] + path.mean_stress_pa[2])?;
        let eigenstrain_work_j_m3 = finite(-delta_free * trace)?;
        let eigenstrain_work_defect_j_m3 = finite(eigenstrain_work_j_m3 - path.path_work_j_m3)?;
        check_work(eigenstrain_work_defect_j_m3, energy_tolerance_j_m3)?;
        let endpoint = &path.endpoint.equilibrium;
        let fracture_dissipated_j_m3 =
            finite(coefficient.damage_dissipated_j_m3 + endpoint.energy.dissipated)?;
        let stored_change_j_m3 =
            finite(endpoint.state.stored_energy_j_m3() - state.stored_energy_j_m3())?;
        let balance_defect_j_m3 = finite(
            (((coefficient.parameter_work_j_m3 + eigenstrain_work_j_m3) - stored_change_j_m3)
                - fracture_dissipated_j_m3)
                - path.endpoint.viscous_dissipated_j_m3,
        )?;
        check_work(balance_defect_j_m3, energy_tolerance_j_m3)?;
        if endpoint.state.damage() < state.damage()
            || endpoint.state.dissipated_j_m3() < state.dissipated_j_m3()
            || fracture_dissipated_j_m3 < 0.0
            || !path.endpoint.viscous_dissipated_j_m3.is_finite()
            || path.endpoint.viscous_dissipated_j_m3 < 0.0
        {
            return Err(BrittleSourceError::Invalid(
                "timed thermal loading cannot heal fracture or create negative viscous loss",
            ));
        }
        Ok(BrittleThermalRateStep {
            coefficient,
            path,
            eigenstrain_work_j_m3,
            fracture_dissipated_j_m3,
            stored_change_j_m3,
            eigenstrain_work_defect_j_m3,
            balance_defect_j_m3,
        })
    }
}
