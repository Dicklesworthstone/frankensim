//! Time-resolved molecular impulse on the ORIGINAL adaptive chamber controller.
//!
//! Each trial compares the coarse interval's load with the SUM of its two fine
//! intervals. Only accepted fine intervals contribute to the returned integral;
//! coarse probes and rejected trials never add impulse, heat or phase mass.
//! The same original snapshots own all streams and loads on every evaluation.
//!
//! Frames, exposed areas and the moment origin are fixed throughout the call.
//! Supports remain stationary: this computes their resolved impulse history,
//! not moving-body dynamics or a bulk-gas momentum state. The error estimate is
//! local in time, not a rigorous bound on the complete integral. Original mass,
//! species and thermal conservation budgets still apply to the entire duration.

use super::{
    ClosedGasParcel, CondensedPhaseStep, GasChamberSurface, GasSurfaceError, MaxwellReflection,
    MolecularChamberImpulse, MolecularSurfaceFrame, MolecularSurfaceImpulse, SharpMeltingCaloric,
    Sum, ThermochemicalGasModel, phase_reports, poll,
};
use crate::sublimation::gas::chamber::advance::{
    GasChamberAdvanceCandidate, GasChamberAdvanceError, GasChamberAdvanceSolve,
    advance_history_observed,
};
use std::collections::BTreeMap;

/// Independent local temporal-error scales for momentum and original states.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MomentumChamberAdvanceSolve {
    pub chamber: GasChamberAdvanceSolve,
    /// Applies separately to each gross beam and total normal impulse on EACH
    /// surface, and to the norm of the complete resultant [N s]. A carrier's
    /// large reflection load cannot hide an inaccurate emission recoil.
    pub absolute_impulse_tolerance_n_s: f64,
    /// Applies to each surface's and the complete resultant's moment norm about
    /// the caller's fixed origin [N m s].
    pub absolute_angular_impulse_tolerance_n_m_s: f64,
    /// Multiplies the larger coarse/fine interval load before addition to the
    /// absolute scale.
    /// Accumulated history never enlarges a subsequent interval's error scale.
    pub relative_impulse_tolerance: f64,
}

/// Original endpoint and conservative stream history, with accepted-only loads.
#[derive(Clone, Debug, PartialEq)]
pub struct MomentumChamberAdvanceCandidate<'catalog> {
    history: GasChamberAdvanceCandidate<'catalog>,
    phases: BTreeMap<u64, CondensedPhaseStep>,
    momentum: MolecularChamberImpulse,
    maximum_accepted_impulse_error_ratio: f64,
}
impl<'catalog> MomentumChamberAdvanceCandidate<'catalog> {
    /// Continue from these original surface/gas states. The original endpoint's
    /// step receipt is its LAST half step; audit() covers the complete window.
    pub const fn history(&self) -> &GasChamberAdvanceCandidate<'catalog> {
        &self.history
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.history.gas()
    }
    /// Whole-window constituent changes, present only with sharp fusion charts.
    pub fn phases(&self) -> &BTreeMap<u64, CondensedPhaseStep> {
        &self.phases
    }
    /// Complete-duration integral, including full carrier reflection pressure.
    /// Every vector uses the same fixed frame and explicit moment origin.
    pub const fn momentum(&self) -> &MolecularChamberImpulse {
        &self.momentum
    }
    pub const fn maximum_accepted_impulse_error_ratio(&self) -> f64 {
        self.maximum_accepted_impulse_error_ratio
    }
    #[allow(clippy::type_complexity)]
    pub fn into_parts(
        self,
    ) -> (
        GasChamberAdvanceCandidate<'catalog>,
        BTreeMap<u64, CondensedPhaseStep>,
        MolecularChamberImpulse,
        f64,
    ) {
        (
            self.history,
            self.phases,
            self.momentum,
            self.maximum_accepted_impulse_error_ratio,
        )
    }
}

/// Advance fixed stationary surfaces with independently resolved force/torque.
/// Optional Maxwell fractions retain their original reflection and heat model.
/// Three original kernel calls per trial and O(surface count) retained sums use
/// the existing surface/trial ceilings. No time-series grows with accepted steps.
/// A failed, exhausted or cancelled window publishes no physical state or load.
#[allow(clippy::too_many_arguments)]
pub fn advance_chamber_adaptive_with_momentum<'catalog>(
    surfaces: &[GasChamberSurface<'_>],
    frames: &[MolecularSurfaceFrame],
    origin_m: [f64; 3],
    scattering: Option<&BTreeMap<u64, MaxwellReflection>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: MomentumChamberAdvanceSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<MomentumChamberAdvanceCandidate<'catalog>, GasChamberAdvanceError> {
    advance(
        surfaces, frames, origin_m, scattering, None, gas, catalog, duration_s, controls, cancelled,
    )
}

/// The same controller with original solid/liquid inventories and latent heat.
/// Exposed area is the explicitly supplied SOLID area, not a bulk phase fraction.
/// Both constituent-mass and molecular-load errors participate in acceptance.
#[allow(clippy::too_many_arguments)]
pub fn advance_partitioned_chamber_adaptive_with_momentum<'catalog>(
    surfaces: &[GasChamberSurface<'_>],
    calorics: &BTreeMap<u64, SharpMeltingCaloric<'_>>,
    frames: &[MolecularSurfaceFrame],
    origin_m: [f64; 3],
    scattering: Option<&BTreeMap<u64, MaxwellReflection>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: MomentumChamberAdvanceSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<MomentumChamberAdvanceCandidate<'catalog>, GasChamberAdvanceError> {
    advance(
        surfaces,
        frames,
        origin_m,
        scattering,
        Some(calorics),
        gas,
        catalog,
        duration_s,
        controls,
        cancelled,
    )
}

#[derive(Default)]
struct SurfaceSums {
    normal: [Sum; 4],
    impulse: [Sum; 3],
    angular: [Sum; 3],
}
fn normal(row: &MolecularSurfaceImpulse) -> [f64; 4] {
    [
        row.outgoing_normal_impulse_n_s,
        row.sticking_normal_impulse_n_s,
        row.reflected_normal_impulse_n_s,
        row.total_normal_impulse_n_s,
    ]
}
impl SurfaceSums {
    fn add(&mut self, row: &MolecularSurfaceImpulse) -> Result<(), GasSurfaceError> {
        for (sum, value) in self.normal.iter_mut().zip(normal(row)) {
            sum.add(value)?;
        }
        for axis in 0..3 {
            self.impulse[axis].add(row.impulse_on_surface_n_s[axis])?;
            self.angular[axis].add(row.angular_impulse_on_surface_n_m_s[axis])?;
        }
        Ok(())
    }
    fn finish(&self) -> Result<MolecularSurfaceImpulse, GasSurfaceError> {
        Ok(MolecularSurfaceImpulse {
            outgoing_normal_impulse_n_s: self.normal[0].value()?,
            sticking_normal_impulse_n_s: self.normal[1].value()?,
            reflected_normal_impulse_n_s: self.normal[2].value()?,
            total_normal_impulse_n_s: self.normal[3].value()?,
            impulse_on_surface_n_s: [
                self.impulse[0].value()?,
                self.impulse[1].value()?,
                self.impulse[2].value()?,
            ],
            angular_impulse_on_surface_n_m_s: [
                self.angular[0].value()?,
                self.angular[1].value()?,
                self.angular[2].value()?,
            ],
        })
    }
}

fn ratio(
    error: f64,
    scale: f64,
    absolute: f64,
    relative: f64,
) -> Result<f64, GasChamberAdvanceError> {
    let scale = absolute + relative * scale;
    let ratio = error / scale;
    if !scale.is_finite() || scale <= 0.0 || !ratio.is_finite() {
        return Err(GasChamberAdvanceError::Invalid(
            "molecular impulse error scale is unrepresentable",
        ));
    }
    Ok(ratio)
}
fn norm(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}
fn vector_ratio(
    a: [f64; 3],
    first: [f64; 3],
    second: [f64; 3],
    absolute: f64,
    relative: f64,
) -> Result<f64, GasChamberAdvanceError> {
    let fine = core::array::from_fn(|axis| first[axis] + second[axis]);
    let difference = core::array::from_fn(|axis| a[axis] - fine[axis]);
    ratio(
        norm(difference),
        norm(a).max(norm(fine)),
        absolute,
        relative,
    )
}
fn compare(
    coarse: &MolecularChamberImpulse,
    first: &MolecularChamberImpulse,
    second: &MolecularChamberImpulse,
    controls: MomentumChamberAdvanceSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, GasChamberAdvanceError> {
    if coarse.surfaces.len() != first.surfaces.len()
        || coarse.surfaces.len() != second.surfaces.len()
    {
        return Err(GasChamberAdvanceError::Invalid(
            "molecular impulse trial roster changed",
        ));
    }
    let absolute = controls.absolute_impulse_tolerance_n_s;
    let angular = controls.absolute_angular_impulse_tolerance_n_m_s;
    let relative = controls.relative_impulse_tolerance;
    let mut maximum = 0.0_f64;
    for (id, a) in &coarse.surfaces {
        poll(cancelled)?;
        let b = first
            .surfaces
            .get(id)
            .ok_or(GasChamberAdvanceError::Invalid(
                "missing first molecular impulse",
            ))?;
        let c = second
            .surfaces
            .get(id)
            .ok_or(GasChamberAdvanceError::Invalid(
                "missing second molecular impulse",
            ))?;
        for ((a, b), c) in normal(a).into_iter().zip(normal(b)).zip(normal(c)) {
            let fine = b + c;
            maximum = maximum.max(ratio(
                (a - fine).abs(),
                a.abs().max(fine.abs()),
                absolute,
                relative,
            )?);
        }
        maximum = maximum.max(vector_ratio(
            a.angular_impulse_on_surface_n_m_s,
            b.angular_impulse_on_surface_n_m_s,
            c.angular_impulse_on_surface_n_m_s,
            angular,
            relative,
        )?);
    }
    maximum = maximum.max(vector_ratio(
        coarse.impulse_on_surfaces_n_s,
        first.impulse_on_surfaces_n_s,
        second.impulse_on_surfaces_n_s,
        absolute,
        relative,
    )?);
    maximum = maximum.max(vector_ratio(
        coarse.angular_impulse_on_surfaces_n_m_s,
        first.angular_impulse_on_surfaces_n_m_s,
        second.angular_impulse_on_surfaces_n_m_s,
        angular,
        relative,
    )?);
    Ok(maximum)
}

#[allow(clippy::too_many_arguments)]
fn advance<'catalog>(
    surfaces: &[GasChamberSurface<'_>],
    frames: &[MolecularSurfaceFrame],
    origin: [f64; 3],
    scattering: Option<&BTreeMap<u64, MaxwellReflection>>,
    calorics: Option<&BTreeMap<u64, SharpMeltingCaloric<'_>>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration: f64,
    controls: MomentumChamberAdvanceSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MomentumChamberAdvanceCandidate<'catalog>, GasChamberAdvanceError> {
    poll(&mut cancelled)?;
    if [
        controls.absolute_impulse_tolerance_n_s,
        controls.absolute_angular_impulse_tolerance_n_m_s,
    ]
    .iter()
    .any(|x| !x.is_finite() || *x <= 0.0)
        || !controls.relative_impulse_tolerance.is_finite()
        || controls.relative_impulse_tolerance < 0.0
        || surfaces.len() > controls.chamber.step.maximum_surfaces
        || frames.len() != surfaces.len()
    {
        return Err(GasChamberAdvanceError::Invalid(
            "bounded frames and positive physical impulse error controls required",
        ));
    }
    let mut sums = BTreeMap::new();
    for frame in frames {
        poll(&mut cancelled)?;
        if sums
            .insert(frame.surface_id, SurfaceSums::default())
            .is_some()
        {
            return Err(GasChamberAdvanceError::Invalid(
                "duplicate molecular history frame",
            ));
        }
    }
    let mut maximum_accepted_impulse_error_ratio = 0.0_f64;
    let history = advance_history_observed(
        surfaces,
        None,
        None,
        calorics,
        gas,
        catalog,
        duration,
        controls.chamber,
        |rows, gas, dt, cancelled| {
            let (candidate, reflection, _, impulse) = super::advance(
                rows,
                frames,
                origin,
                scattering,
                calorics,
                gas,
                catalog,
                dt,
                controls.chamber.step,
                cancelled,
            )?
            .into_parts();
            Ok((candidate, reflection, impulse))
        },
        |coarse, first, second, cancelled| compare(coarse, first, second, controls, cancelled),
        |first, second, error, cancelled| {
            for part in [first, second] {
                for (id, row) in &part.surfaces {
                    poll(cancelled)?;
                    sums.get_mut(id)
                        .ok_or(GasChamberAdvanceError::Invalid(
                            "accepted molecular history roster changed",
                        ))?
                        .add(row)?;
                }
            }
            maximum_accepted_impulse_error_ratio = maximum_accepted_impulse_error_ratio.max(error);
            Ok(())
        },
        &mut cancelled,
    )?;
    let phases = match calorics {
        Some(laws) => phase_reports(
            surfaces,
            laws,
            history.endpoint(),
            controls.chamber.step.surface.gas.mass_tolerance_kg,
            &mut cancelled,
        )?,
        None => BTreeMap::new(),
    };
    let mut loads = BTreeMap::new();
    let mut impulse: [Sum; 3] = core::array::from_fn(|_| Sum::default());
    let mut angular: [Sum; 3] = core::array::from_fn(|_| Sum::default());
    for (id, sum) in sums {
        poll(&mut cancelled)?;
        let row = sum.finish()?;
        for axis in 0..3 {
            impulse[axis].add(row.impulse_on_surface_n_s[axis])?;
            angular[axis].add(row.angular_impulse_on_surface_n_m_s[axis])?;
        }
        loads.insert(id, row);
    }
    poll(&mut cancelled)?;
    Ok(MomentumChamberAdvanceCandidate {
        momentum: MolecularChamberImpulse {
            duration_s: history.audit().elapsed_s,
            origin_m: origin,
            surfaces: loads,
            impulse_on_surfaces_n_s: [
                impulse[0].value()?,
                impulse[1].value()?,
                impulse[2].value()?,
            ],
            angular_impulse_on_surfaces_n_m_s: [
                angular[0].value()?,
                angular[1].value()?,
                angular[2].value()?,
            ],
        },
        history,
        phases,
        maximum_accepted_impulse_error_ratio,
    })
}

#[cfg(test)]
mod tests;
