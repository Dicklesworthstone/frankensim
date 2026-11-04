//! Native deformable debonding, finite-body cooling and explicit loss heating.
//!
//! Mechanics and thermal transport use the accepted coefficients and contact
//! set over one straight mechanical segment. Endpoint fracture/reclosure then
//! changes the NEXT thermal interval, and explicitly converted physical losses
//! enter the original caloric cells once. This is a staggered split, not an
//! event-localized or simultaneous temperature-dependent constitutive solve.

use fs_conduction::transient::enthalpy::deposition::fracture::{
    CohesiveHeatPartition, Error as FractureCoolingError, FractureCooling, FractureHeatReport,
};
use fs_exec::CancelGate;
use fs_solid::cohesive::dynamics::{CohesiveMotion, DynamicsError, DynamicsStep};
use fs_thermal::enthalpy::{DepositionError, EnergyAudit, StepConfig};

/// Native thermal controls and absolute budgets for conversion and total energy.
/// Mechanical bounds remain with the original `CohesiveMotion`.
#[derive(Clone, Copy, Debug)]
pub struct CohesiveCoolingControls {
    pub thermal: StepConfig,
    pub heat_energy_tolerance_j: f64,
    pub combined_energy_tolerance_j: f64,
}

#[derive(Debug)]
pub enum CohesiveCoolingError {
    Invalid(&'static str),
    Mechanical(DynamicsError),
    Cooling(FractureCoolingError),
    Energy { residual_j: f64, tolerance_j: f64 },
}
impl From<DynamicsError> for CohesiveCoolingError {
    fn from(error: DynamicsError) -> Self {
        Self::Mechanical(error)
    }
}
impl From<FractureCoolingError> for CohesiveCoolingError {
    fn from(error: FractureCoolingError) -> Self {
        Self::Cooling(error)
    }
}
impl From<DepositionError> for CohesiveCoolingError {
    fn from(error: DepositionError) -> Self {
        Self::Cooling(error.into())
    }
}
impl core::fmt::Display for CohesiveCoolingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "deformable fracture cooling refused: {self:?}")
    }
}
impl std::error::Error for CohesiveCoolingError {}

/// Both original owners advance together only after every native/combined gate.
pub struct CohesiveCoolingCandidate<'s, 'a, 't> {
    pub mechanical: DynamicsStep<'s, 'a>,
    pub cooling: FractureCooling<'t, 'a>,
    /// Heat transported with the PREVIOUS accepted closure held fixed.
    pub transport: EnergyAudit,
    /// Endpoint closure and heat, including both unconverted physical losses.
    pub heat: FractureHeatReport,
    /// Independently reduced delta(K + U_bulk + U_bond + H) plus nonthermal
    /// losses, minus applied/prescribed work, birth energy and external heat.
    pub combined_energy_residual_j: f64,
    /// Transport residual + heat-transfer residual - mechanical input residual.
    pub summed_energy_residual_j: f64,
}

/// Advance actual native inertia/deformation through debonding and reclosure.
///
/// Every mechanical patch must retain the same ID, exact original area, law and
/// complete accepted point as its cooling interface. The caller declares the
/// fixed-reference coincidence and ordered A/B thermal-cell/body/material
/// correspondence: a reduced jump map cannot establish that geometry itself.
/// All three mechanical, thermal and synchronized-history clocks must agree.
/// `target_time_s` is exact; the represented elapsed time drives BOTH viscosity
/// and thermal transport, including at nonbinary absolute starting times.
///
/// Original mechanics advances once; transport then advances with the OLD
/// closure. The existing cooling owner verifies the same straight material path
/// and deposits only the explicit fractions of NEW fracture and viscosity loss.
/// This replay validates history and energy; it is not another physical loss.
/// Nonthermal remainders stay in `heat`. Nonzero damaged-area friction refuses:
/// that channel requires its own declared receiving/heat owner.
///
/// Opening/reclosure affects the next interval. Refine/split the step to resolve
/// events and feedback; coefficients remain frozen and no temperature-dependent
/// source refresh is inferred. The cooling policy retains its original slip,
/// frame and area limits, independent bonded/reclosed conductances, caloric
/// domains, events and native thermal work bounds. All candidates stay private
/// until the final shared cancellation and energy gates succeed.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_deformable_cohesive_cooling<'s, 'a, 't>(
    motion: &CohesiveMotion<'s, 'a>,
    cooling: &FractureCooling<'t, 'a>,
    mean_force_n: &[f64],
    endpoint_prescribed_jumps_m: &[[f64; 3]],
    target_time_s: f64,
    partitions: &[CohesiveHeatPartition],
    controls: CohesiveCoolingControls,
    gate: &CancelGate,
) -> Result<CohesiveCoolingCandidate<'s, 'a, 't>, CohesiveCoolingError> {
    let check = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    check()?;
    let dt_s = target_time_s - motion.time_s();
    let patches = motion.state().system().patches();
    if motion.time_s() != cooling.thermal().time_s()
        || motion.time_s() != cooling.last_synchronized_time_s()
        || !target_time_s.is_finite()
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || patches.is_empty()
        || patches.len() != cooling.interfaces().len()
        || !controls.thermal.max_dt_s.is_finite()
        || controls.thermal.max_dt_s <= 0.0
        || [
            controls.thermal.energy_tolerance_j,
            controls.heat_energy_tolerance_j,
            controls.combined_energy_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(CohesiveCoolingError::Invalid(
            "complete synchronized original owners and finite positive interval required",
        ));
    }
    for (patch, point) in patches.iter().zip(motion.state().points()) {
        check()?;
        let interface =
            cooling
                .interfaces()
                .get(&patch.id())
                .ok_or(CohesiveCoolingError::Invalid(
                    "unknown original fracture-cooling patch",
                ))?;
        if interface.reference.area_m2 != patch.area_m2()
            || interface.point != *point
            || !core::ptr::eq(interface.point.law(), point.law())
        {
            return Err(CohesiveCoolingError::Invalid(
                "fracture cooling must retain the original patch area, law and complete history",
            ));
        }
    }
    let mechanical = motion.advance_to_time_controlled(
        mean_force_n,
        endpoint_prescribed_jumps_m,
        target_time_s,
        &mut || gate.is_requested(),
    )?;
    if mechanical.friction_dissipated_j != 0.0 {
        return Err(CohesiveCoolingError::Invalid(
            "nonzero residual-friction loss needs its own explicit thermal transfer owner",
        ));
    }
    let (transported, transport) = cooling.advance_to(target_time_s, controls.thermal, check)?;
    let proposed: Vec<_> = patches
        .iter()
        .zip(mechanical.motion.state().points())
        .map(|(patch, point)| (patch.id(), *point))
        .collect();
    let (next, heat) = transported.synchronize_over_time_with_heat(
        target_time_s,
        dt_s,
        &proposed,
        partitions,
        motion.dynamics_config().solve.local_energy_tolerance_j_m2,
        controls.heat_energy_tolerance_j,
        check,
    )?;

    // The authoritative native mechanical and cooling replays must agree on
    // BOTH losses. Equal endpoints alone cannot establish a heat/work account.
    for difference in [
        heat.contact.fracture_dissipation_j - mechanical.cohesive_energy_j.dissipated,
        heat.contact.viscous_dissipation_j - mechanical.viscous_dissipated_j,
    ] {
        require_energy(difference, controls.combined_energy_tolerance_j)?;
    }
    let mut balance = EnergySum::default();
    for term in [
        mechanical.kinetic_after_j - mechanical.kinetic_before_j,
        mechanical.bulk_stored_after_j - mechanical.bulk_stored_before_j,
        mechanical.cohesive_stored_after_j - mechanical.cohesive_stored_before_j,
        heat.nonthermal_fracture_energy_j,
        heat.nonthermal_viscous_energy_j,
        -mechanical.external_work_j,
        -mechanical.prescribed_jump_work_j,
        -transport.incoming_energy_j,
        -transport.external_heat_j,
    ] {
        balance.add(term)?;
    }
    // Independently inspect represented caloric endpoints. Subtract per-cell
    // enthalpy coordinates BEFORE summation to retain small increments beside
    // large formation-energy references. Births contribute their actual energy.
    for (&id, cell) in next.thermal().cells() {
        check()?;
        let delta_h = cooling
            .thermal()
            .cells()
            .get(&id)
            .map_or(cell.specific_h_j_kg, |old| {
                cell.specific_h_j_kg - old.specific_h_j_kg
            });
        balance.add(cell.mass_kg * delta_h)?;
    }
    let combined = balance.value()?;
    let summed = transport.residual_j + heat.transfer_residual_j - mechanical.energy_residual_j;
    for residual in [combined, summed, combined - summed] {
        require_energy(residual, controls.combined_energy_tolerance_j)?;
    }
    check()?;
    Ok(CohesiveCoolingCandidate {
        mechanical,
        cooling: next,
        transport,
        heat,
        combined_energy_residual_j: combined,
        summed_energy_residual_j: summed,
    })
}

fn require_energy(residual_j: f64, tolerance_j: f64) -> Result<(), CohesiveCoolingError> {
    if !residual_j.is_finite() || residual_j.abs() > tolerance_j {
        Err(CohesiveCoolingError::Energy {
            residual_j,
            tolerance_j,
        })
    } else {
        Ok(())
    }
}

#[derive(Default)]
struct EnergySum {
    high: f64,
    low: f64,
}
impl EnergySum {
    fn add(&mut self, value: f64) -> Result<(), CohesiveCoolingError> {
        let (high, error) = fs_math::eft::two_sum(self.high, value);
        let low = self.low + error;
        if !high.is_finite() || !low.is_finite() {
            return Err(DepositionError::Arithmetic.into());
        }
        self.high = high;
        self.low = low;
        Ok(())
    }
    fn value(&self) -> Result<f64, CohesiveCoolingError> {
        let value = self.high + self.low;
        if !value.is_finite() {
            Err(DepositionError::Arithmetic.into())
        } else {
            Ok(value)
        }
    }
}
