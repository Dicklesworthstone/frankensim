//! One irreversible magnetic history coupled to multiple physical windings.
//!
//! The existing uniform core is an explicitly chosen reference winding. With
//! signed turns ratios r and positive leakage D, total linkages obey
//! lambda_i = D_i i_i + r_i lambda_ref and I_ref = sum(r_i i_i).
//! Eliminating winding currents gives a positive linear load line on the SAME
//! original core law. Its pinning heat and internal polarizations occur once.
//!
//! The work-conjugate mean current integrates the exact piecewise-linear
//! magnetic path along a straight chord in TOTAL winding linkage. Endpoint
//! current averages alone fail when the common core crosses a pinning or
//! saturation knee. Reversals of the common drive must be separate increments.
//! Geometry and material coefficients are fixed; copper/eddy heat and motion
//! belong to their original circuit/field/mechanical owners.

use super::{
    FieldStorageError, HystereticCoreStep, MagneticHysteresisState, UniformHystereticCore, invalid,
    positive,
};
use fs_matdb::QueryPoint;
use std::sync::Arc;

/// One physical winding coupled to the existing reference-core coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedHystereticWinding {
    /// Signed physical turns / reference-core turns. Both must be explicit.
    pub turns_ratio: f64,
    /// Positive fixed physical leakage inductance [H], never regularization.
    pub leakage_h: f64,
}

#[derive(Clone, Debug, PartialEq)]
struct Geometry {
    windings: Vec<SharedHystereticWinding>,
    drive_weights: Vec<f64>,
    loading_per_h: f64,
    load_line_h_m: f64,
    linkage_per_t: f64,
    current_per_field: f64,
    volume_m3: f64,
}

/// Source-backed reference core and explicit leakage/turns for 1..=16 windings.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedHystereticCore {
    core: UniformHystereticCore,
    geometry: Arc<Geometry>,
}

/// One accepted common polarization history and all physical total linkages.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedHystereticState {
    geometry: Arc<Geometry>,
    core: MagneticHysteresisState,
    flux_linkages_wb: Vec<f64>,
    currents_a: Vec<f64>,
    stored_energy_j: f64,
}
impl SharedHystereticState {
    pub fn flux_linkages_wb(&self) -> &[f64] {
        &self.flux_linkages_wb
    }
    pub fn currents_a(&self) -> &[f64] {
        &self.currents_a
    }
    pub const fn core_state(&self) -> &MagneticHysteresisState {
        &self.core
    }
    /// One common core plus all explicit leakage energy [J].
    pub const fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }
}

/// Transactional candidate; native circuit work uses `mean_currents_a`.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedHystereticStep {
    pub state: SharedHystereticState,
    /// Original core result, retaining exactly one pinning heat/source account.
    pub core: HystereticCoreStep,
    pub mean_currents_a: Vec<f64>,
    /// Stable change of complete core and leakage energy [J].
    pub stored_change_j: f64,
    pub stored_energy_j: f64,
    /// Exact monotone-path winding work = stored change + core pinning heat [J].
    pub magnetic_work_j: f64,
}

fn finite(value: f64) -> Result<f64, FieldStorageError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("shared magnetic arithmetic is unrepresentable"))
    }
}

impl SharedHystereticCore {
    pub fn new(
        core: UniformHystereticCore,
        windings: Vec<SharedHystereticWinding>,
    ) -> Result<Self, FieldStorageError> {
        if windings.is_empty() || windings.len() > 16 {
            return Err(invalid(
                "one to sixteen physical shared-core windings required",
            ));
        }
        let mut drive_weights = Vec::with_capacity(windings.len());
        let mut loading_per_h = 0.0;
        for winding in &windings {
            if !winding.turns_ratio.is_finite()
                || winding.turns_ratio == 0.0
                || !positive(winding.leakage_h)
                || !positive(winding.leakage_h.recip())
            {
                return Err(invalid(
                    "shared windings require signed nonzero ratios and positive physical leakage",
                ));
            }
            let weight = finite(winding.turns_ratio / winding.leakage_h)?;
            let loading = finite(weight * winding.turns_ratio)?;
            if !positive(loading) {
                return Err(invalid("shared magnetic load line underflows"));
            }
            loading_per_h = finite(loading_per_h + loading)?;
            drive_weights.push(weight);
        }
        let load_line_h_m = (core.current_per_field / core.linkage_per_t) / loading_per_h;
        if !positive(load_line_h_m) {
            return Err(invalid("shared magnetic load line is unrepresentable"));
        }
        let geometry = Arc::new(Geometry {
            windings,
            drive_weights,
            loading_per_h,
            load_line_h_m,
            linkage_per_t: core.linkage_per_t,
            current_per_field: core.current_per_field,
            volume_m3: core.volume_m3,
        });
        Ok(Self { core, geometry })
    }

    pub const fn core(&self) -> &UniformHystereticCore {
        &self.core
    }
    pub fn windings(&self) -> &[SharedHystereticWinding] {
        &self.geometry.windings
    }

    /// Virgin initialization. Prepared/remanent specimens use state_from_core.
    pub fn initial_state(&self) -> SharedHystereticState {
        SharedHystereticState {
            geometry: Arc::clone(&self.geometry),
            core: self.core.initial_state(),
            flux_linkages_wb: vec![0.0; self.windings().len()],
            currents_a: vec![0.0; self.windings().len()],
            stored_energy_j: 0.0,
        }
    }

    /// Bind explicitly prepared core memory and TOTAL winding linkage. The
    /// magnetic balance is checked without changing any polarization/current.
    pub fn state_from_core(
        &self,
        core_state: MagneticHysteresisState,
        flux_linkages_wb: Vec<f64>,
    ) -> Result<SharedHystereticState, FieldStorageError> {
        let core = self.core.hold_candidate(&core_state)?;
        self.state_at(&core, flux_linkages_wb)
    }

    /// Recheck the same source coefficients at a new actual query, retaining
    /// fixed geometry and all history. Changing fit parameters still refuses.
    pub fn at_source_query(&self, point: &QueryPoint) -> Result<Self, FieldStorageError> {
        Ok(Self {
            core: self.core.at_source_query(point)?,
            geometry: Arc::clone(&self.geometry),
        })
    }

    pub fn hold_candidate(
        &self,
        previous: &SharedHystereticState,
    ) -> Result<SharedHystereticStep, FieldStorageError> {
        self.candidate_linkages(previous, previous.flux_linkages_wb())
    }

    /// Frozen-memory active-set tangent d(lambda)/di [H], row-major in winding
    /// order. At pinning/saturation thresholds the original material's selected
    /// one-sided tangent is retained; this is not a complex AC hysteresis law.
    pub fn incremental_inductance_h(
        &self,
        state: &SharedHystereticState,
    ) -> Result<Vec<f64>, FieldStorageError> {
        if state.geometry != self.geometry {
            return Err(invalid(
                "shared magnetic tangent requires its original winding geometry",
            ));
        }
        let core = self.core.hold_candidate(state.core_state())?;
        let n = self.windings().len();
        let mut matrix = vec![0.0; n * n];
        for (i, a) in self.windings().iter().enumerate() {
            for (j, b) in self.windings()[..=i].iter().enumerate() {
                let value = finite(
                    finite(a.turns_ratio * core.incremental_inductance_h)?
                        .mul_add(b.turns_ratio, if i == j { a.leakage_h } else { 0.0 }),
                )?;
                if value == 0.0 {
                    return Err(invalid("nonzero shared magnetic tangent underflows"));
                }
                matrix[i * n + j] = value;
                matrix[j * n + i] = value;
            }
        }
        Ok(matrix)
    }

    pub fn candidate_linkages(
        &self,
        previous: &SharedHystereticState,
        linkages: &[f64],
    ) -> Result<SharedHystereticStep, FieldStorageError> {
        if previous.geometry != self.geometry
            || linkages.len() != self.windings().len()
            || linkages.iter().any(|x| !x.is_finite())
        {
            return Err(invalid(
                "shared magnetic state must retain its geometry and complete finite linkages",
            ));
        }
        let old_core = self.core.hold_candidate(previous.core_state())?;
        let mut drive = 0.0;
        let mut delta_drive = 0.0;
        for ((&weight, &after), &before) in self
            .geometry
            .drive_weights
            .iter()
            .zip(linkages)
            .zip(previous.flux_linkages_wb())
        {
            drive = finite(weight.mul_add(after, drive))?;
            delta_drive = finite(weight.mul_add(finite(after - before)?, delta_drive))?;
        }
        let core = if delta_drive == 0.0 {
            old_core.clone()
        } else {
            // I_ref + k*lambda_ref = drive, with lambda_ref=b*B and
            // I_ref=a*H. Divide by k*b and reuse the original material inverse.
            let target = finite((drive / self.geometry.loading_per_h) / self.core.linkage_per_t)?;
            self.core.finish(self.core.law.candidate_linear_flux(
                previous.core_state(),
                target,
                self.geometry.load_line_h_m,
            )?)?
        };
        let state = self.state_at(&core, linkages.to_vec())?;
        let delta_core = finite(core.flux_linkage_wb - old_core.flux_linkage_wb)?;
        if delta_drive != 0.0 && delta_core != 0.0 && delta_drive.signum() != delta_core.signum() {
            return Err(invalid(
                "shared magnetic inverse changed the monotone drive direction",
            ));
        }
        let mut mean_currents_a: Vec<_> = previous
            .currents_a()
            .iter()
            .zip(state.currents_a())
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect();
        let mut mean_mmf = 0.0;
        let mut leakage_change = 0.0;
        for ((winding, &middle), (&before, &after)) in self
            .windings()
            .iter()
            .zip(&mean_currents_a)
            .zip(previous.flux_linkages_wb().iter().zip(linkages))
        {
            mean_mmf = finite(winding.turns_ratio.mul_add(middle, mean_mmf))?;
            // D*delta(i) = delta(lambda) - r*delta(lambda_ref). Evaluate
            // directly in the conserved coordinates: subtracting reconstructed
            // near-equal currents loses resolvable small-step leakage work.
            let leakage_increment =
                finite((-winding.turns_ratio).mul_add(delta_core, finite(after - before)?))?;
            leakage_change = finite(middle.mul_add(leakage_increment, leakage_change))?;
        }
        // Integrate Phi ds by parts: mean(Phi)=Phi_mid+
        // (Iref_mid*dPhi-Wcore)/ds. Applying the corresponding correction
        // directly to currents avoids cancellation in large common fluxes.
        // The endpoint winding MMF includes only constitutive solve roundoff;
        // using it here closes the exact retained leakage/core energy account.
        let defect = finite(mean_mmf.mul_add(delta_core, -core.magnetic_work_j))?;
        if delta_drive == 0.0 {
            if delta_core != 0.0 || defect != 0.0 {
                return Err(invalid(
                    "zero shared magnetic drive cannot advance core memory",
                ));
            }
        } else {
            let correction = finite(defect / delta_drive)?;
            for (&weight, current) in self.geometry.drive_weights.iter().zip(&mut mean_currents_a) {
                *current = finite((-weight).mul_add(correction, *current))?;
            }
        }
        let core_change = finite(self.core.volume_m3 * core.material.energy_change_j_m3)?;
        let stored_change_j = finite(leakage_change + core_change)?;
        let magnetic_work_j = finite(stored_change_j + core.core_heat_j)?;
        let mut paired_work = 0.0;
        let mut scale = magnetic_work_j.abs() + stored_change_j.abs() + core.core_heat_j.abs();
        for ((&current, &after), &before) in mean_currents_a
            .iter()
            .zip(linkages)
            .zip(previous.flux_linkages_wb())
        {
            let contribution = finite(current * (after - before))?;
            paired_work = finite(paired_work + contribution)?;
            scale = finite(scale + contribution.abs())?;
        }
        if (paired_work - magnetic_work_j).abs() > 2048.0 * f64::EPSILON * scale {
            return Err(invalid(
                "shared magnetic winding work does not close core and leakage storage plus pinning heat",
            ));
        }
        Ok(SharedHystereticStep {
            stored_energy_j: state.stored_energy_j(),
            state,
            core,
            mean_currents_a,
            stored_change_j,
            magnetic_work_j,
        })
    }

    fn state_at(
        &self,
        core: &HystereticCoreStep,
        linkages: Vec<f64>,
    ) -> Result<SharedHystereticState, FieldStorageError> {
        if linkages.len() != self.windings().len() || linkages.iter().any(|x| !x.is_finite()) {
            return Err(invalid("complete finite shared winding linkages required"));
        }
        let mut currents = Vec::with_capacity(linkages.len());
        let mut mmf = 0.0;
        let mut scale = core.current_a.abs();
        let mut stored_energy_j = core.stored_energy_j;
        for (winding, &linkage) in self.windings().iter().zip(&linkages) {
            let current = finite(
                (-winding.turns_ratio).mul_add(core.flux_linkage_wb, linkage) / winding.leakage_h,
            )?;
            mmf = finite(winding.turns_ratio.mul_add(current, mmf))?;
            scale = finite(scale + finite(winding.turns_ratio * current)?.abs())?;
            stored_energy_j =
                finite((0.5 * winding.leakage_h * current).mul_add(current, stored_energy_j))?;
            currents.push(current);
        }
        // Include the common flux subtraction scale: large coupling may make
        // small leakage currents sensitive to the admitted core inverse error.
        scale = finite(scale + finite(self.geometry.loading_per_h * core.flux_linkage_wb)?.abs())?;
        if !stored_energy_j.is_finite()
            || stored_energy_j < 0.0
            || (mmf - core.current_a).abs() > 256.0 * f64::EPSILON * scale
        {
            return Err(invalid(
                "prepared common core and winding linkages violate magnetic balance",
            ));
        }
        Ok(SharedHystereticState {
            geometry: Arc::clone(&self.geometry),
            core: core.material.state.clone(),
            flux_linkages_wb: linkages,
            currents_a: currents,
            stored_energy_j,
        })
    }
}
