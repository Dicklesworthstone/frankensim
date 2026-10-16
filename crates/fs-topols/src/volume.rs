//! Hard material-area projection of a bilinear level set, without remeshing.
//!
//! Free nodal values receive one common offset; explicitly fixed nodal values
//! do not move. Increasing the offset can only remove material from the exact
//! bilinear field. Area is measured by the SAME cut-quadrature functional used
//! by the elasticity optimizer. The tolerance therefore concerns that numerical
//! area, not a certified continuum volume. Flat fields, pinned regions and a
//! bounded shift can make the requested area unattainable; they never receive a
//! fabricated success. Neither a refusal nor interruption mutates the input.

use crate::optimize::material_volume;
use crate::GridSdf;
use fs_cutfem::{CutFemError, Quadtree};
use std::convert::Infallible;
use std::ops::ControlFlow;

/// Explicit normalized-area, displacement and work limits for projection.
#[derive(Debug, Clone, Copy)]
pub struct VolumeProjectionSettings {
    /// Target material area in the normalized unit square, strictly in (0, 1).
    pub target: f64,
    /// Absolute tolerance on the cut-quadrature area; strictly in (0, 1).
    pub tolerance: f64,
    /// Maximum magnitude of the common offset to FREE nodal values.
    pub max_shift: f64,
    /// Maximum complete area evaluations, including the origin and bracket.
    pub max_evaluations: usize,
}

/// Evidence for the exact field published by a successful projection.
#[derive(Debug, Clone, Copy)]
pub struct VolumeProjectionReport {
    /// Area after imposing fixed nodal values, before shifting free nodes.
    pub initial_volume: f64,
    /// Measured area of the field actually returned to the caller.
    pub volume: f64,
    /// Common offset applied to free nodes of the original input.
    pub shift: f64,
    /// Number of complete area evaluations, including bracket evaluations.
    pub evaluations: usize,
}

/// Cooperative checkpoints; one area evaluation is not internally preemptible.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VolumeProjectionStage {
    /// Before staging geometry or allocating the quadtree.
    Prepare,
    /// Before an area evaluation; the count is one-based.
    Measure(usize),
    /// After the numerical area gate, before changing caller geometry.
    Publish,
}

fn refused(what: impl Into<String>) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: what.into() }
}

fn validate(
    phi: &GridSdf,
    level: u32,
    fixed: &[(usize, f64)],
    settings: VolumeProjectionSettings,
) -> Result<(), CutFemError> {
    if !(1..=8).contains(&level) || phi.n() != (1usize << level) {
        return Err(refused("volume projection requires a matching grid level in 1..=8"));
    }
    if phi.nodes().iter().any(|value| !value.is_finite()) {
        return Err(refused("volume projection requires finite nodal values"));
    }
    if !(settings.target.is_finite() && settings.target > 0.0 && settings.target < 1.0
        && settings.tolerance.is_finite() && settings.tolerance > 0.0
        && settings.tolerance < settings.target.min(1.0 - settings.target))
    {
        return Err(refused("volume projection needs target in (0, 1) and a positive tolerance smaller than both target and 1-target"));
    }
    if !(settings.max_shift.is_finite() && settings.max_shift > 0.0)
        || !(3..=128).contains(&settings.max_evaluations)
    {
        return Err(refused("volume projection requires finite positive max_shift and 3..=128 area evaluations"));
    }
    if fixed.len() > phi.nodes().len() {
        return Err(refused("fixed-node count exceeds the level-set lattice"));
    }
    let mut previous = None;
    for &(index, value) in fixed {
        if index >= phi.nodes().len() || !value.is_finite()
            || previous.is_some_and(|old| index <= old)
        {
            return Err(refused("fixed nodes must have finite values and strictly increasing, unique in-range indices"));
        }
        previous = Some(index);
    }
    Ok(())
}

/// Project onto a fixed numerical material-area equality transactionally.
///
/// `fixed` contains row-major nodal indices and their prescribed values, in
/// strictly increasing index order. Preserving nodes preserves their bilinear
/// edges exactly; it does not implicitly prescribe an entire adjacent cell.
/// All other nodes receive a uniform offset in `[-max_shift, max_shift]`.
///
/// # Errors
/// Refuses malformed inputs, an unattained bracket, non-finite arithmetic,
/// exhausted evaluation budget or a non-representable bisection midpoint.
/// Caller geometry is unchanged on every error.
pub fn project_material_volume(
    phi: &mut GridSdf,
    level: u32,
    fixed: &[(usize, f64)],
    settings: VolumeProjectionSettings,
) -> Result<VolumeProjectionReport, CutFemError> {
    match project_material_volume_controlled(
        phi, level, fixed, settings,
        |_| ControlFlow::<Infallible>::Continue(()),
    )? {
        ControlFlow::Continue(report) => Ok(report),
        ControlFlow::Break(never) => match never {},
    }
}

/// Controlled counterpart of [`project_material_volume`].
///
/// `Break` at ANY callback (including `Publish`) leaves the exact input bits
/// untouched. No wall-clock bound is claimed for an individual area evaluation.
pub fn project_material_volume_controlled<B>(
    phi: &mut GridSdf,
    level: u32,
    fixed: &[(usize, f64)],
    settings: VolumeProjectionSettings,
    mut control: impl FnMut(VolumeProjectionStage) -> ControlFlow<B>,
) -> Result<ControlFlow<B, VolumeProjectionReport>, CutFemError> {
    validate(phi, level, fixed, settings)?;
    if let ControlFlow::Break(reason) = control(VolumeProjectionStage::Prepare) {
        return Ok(ControlFlow::Break(reason));
    }
    let grid = Quadtree::uniform(level);
    let mut staged = phi.clone();
    let mut lo = -settings.max_shift;
    let mut hi = settings.max_shift;
    let mut low_volume = None;
    let mut high_volume = None;
    let mut initial_volume = 0.0;
    let mut shift = 0.0;
    let mut best_error = f64::INFINITY;
    let mut best_volume = 0.0;

    for evaluation in 1..=settings.max_evaluations {
        if let ControlFlow::Break(reason) = control(VolumeProjectionStage::Measure(evaluation)) {
            return Ok(ControlFlow::Break(reason));
        }
        let mut fixed_cursor = 0;
        for (index, (value, original)) in staged.nodes_mut().iter_mut()
            .zip(phi.nodes()).enumerate()
        {
            if fixed.get(fixed_cursor).is_some_and(|&(node, _)| node == index) {
                *value = fixed[fixed_cursor].1;
                fixed_cursor += 1;
            } else {
                // Preserve signed zero when the required offset is exactly zero.
                *value = if shift == 0.0 { *original } else { *original + shift };
            }
            if !value.is_finite() {
                return Err(refused("volume projection offset produced a non-finite node"));
            }
        }
        let volume = material_volume(&grid, &staged);
        if !(volume.is_finite() && (0.0..=1.0 + 1e-10).contains(&volume)) {
            return Err(refused("volume projection received an invalid numerical material area"));
        }
        if evaluation == 1 {
            initial_volume = volume;
        }
        let error = (volume - settings.target).abs();
        if error < best_error {
            best_error = error;
            best_volume = volume;
        }
        if error <= settings.tolerance {
            if let ControlFlow::Break(reason) = control(VolumeProjectionStage::Publish) {
                return Ok(ControlFlow::Break(reason));
            }
            *phi = staged;
            return Ok(ControlFlow::Continue(VolumeProjectionReport {
                initial_volume, volume, shift, evaluations: evaluation,
            }));
        }
        match evaluation {
            1 => shift = lo,
            2 => {
                low_volume = Some(volume);
                shift = hi;
            }
            _ => {
                if evaluation == 3 {
                    high_volume = Some(volume);
                    if low_volume.is_none_or(|area| area < settings.target)
                        || volume > settings.target
                    {
                        return Err(refused(format!(
                            "volume target {} is not bracketed within the declared shift: areas {:?}..{}; fixed regions or max_shift limit attainment",
                            settings.target, low_volume, volume,
                        )));
                    }
                } else {
                    // Inclusion is exact for the represented fields, but the
                    // area measurement is numerical. Never silently accept a
                    // significantly non-monotone quadrature response.
                    if low_volume.is_some_and(|area| volume > area + settings.tolerance)
                        || high_volume.is_some_and(|area| volume < area - settings.tolerance)
                    {
                        return Err(refused("volume projection area evaluations violate the monotone bracket"));
                    }
                    if volume > settings.target {
                        lo = shift;
                        low_volume = Some(volume);
                    } else {
                        hi = shift;
                        high_volume = Some(volume);
                    }
                }
                let midpoint = 0.5 * lo + 0.5 * hi;
                if midpoint == lo || midpoint == hi {
                    return Err(refused("volume projection exhausted representable offsets before meeting the area tolerance"));
                }
                shift = midpoint;
            }
        }
    }
    Err(refused(format!(
        "volume projection exhausted {} area evaluations: best area {best_volume}, target {}, residual {best_error}",
        settings.max_evaluations, settings.target,
    )))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn controls(target: f64) -> VolumeProjectionSettings {
        VolumeProjectionSettings { target, tolerance: 1e-6, max_shift: 2.0, max_evaluations: 64 }
    }

    fn bits(phi: &GridSdf) -> Vec<u64> {
        phi.nodes().iter().map(|value| value.to_bits()).collect()
    }

    #[test]
    fn g0_plane_projection_grows_and_shrinks_the_actual_material() {
        for target in [0.25, 0.75] {
            let mut phi = GridSdf::from_fn(8, &|_, y| y - 0.5);
            let report = project_material_volume(&mut phi, 3, &[], controls(target)).unwrap();
            assert!((report.volume - target).abs() <= 1e-6);
            assert!((report.shift - (0.5 - target)).abs() <= 1e-6);
            let measured = material_volume(&Quadtree::uniform(3), &phi);
            assert_eq!(measured.to_bits(), report.volume.to_bits());
        }
    }

    #[test]
    fn g3_fixed_boundary_nodes_and_deterministic_replay() {
        let origin = GridSdf::from_fn(8, &|_, y| y - 0.5);
        let fixed: Vec<_> = (0..=8).map(|j| (j * 9, origin.node(0, j))).collect();
        let mut first = origin.clone();
        let mut second = origin.clone();
        let a = project_material_volume(&mut first, 3, &fixed, controls(0.4)).unwrap();
        let b = project_material_volume(&mut second, 3, &fixed, controls(0.4)).unwrap();
        assert_eq!(bits(&first), bits(&second));
        assert_eq!(a.shift.to_bits(), b.shift.to_bits());
        for &(index, value) in &fixed {
            assert_eq!(first.nodes()[index].to_bits(), value.to_bits());
        }
        assert!((a.volume - 0.4).abs() <= 1e-6);
    }

    #[test]
    fn g4_unattainable_or_exhausted_projection_never_changes_geometry() {
        let origin = GridSdf::from_fn(8, &|_, y| y - 0.5);
        let fixed: Vec<_> = origin.nodes().iter().copied().enumerate().collect();
        let mut phi = origin.clone();
        assert!(project_material_volume(&mut phi, 3, &fixed, controls(0.25)).is_err());
        assert_eq!(bits(&phi), bits(&origin));
        let short = VolumeProjectionSettings { max_evaluations: 3, ..controls(0.3) };
        assert!(project_material_volume(&mut phi, 3, &[], short).is_err());
        assert_eq!(bits(&phi), bits(&origin));
        for fixed in [vec![(0, 0.0), (0, 1.0)], vec![(100, 0.0)], vec![(0, f64::NAN)]] {
            assert!(project_material_volume(&mut phi, 3, &fixed, controls(0.3)).is_err());
            assert_eq!(bits(&phi), bits(&origin));
        }
    }

    #[test]
    fn g4_cancellation_before_measurement_and_publication_is_transactional() {
        for stop in [VolumeProjectionStage::Prepare, VolumeProjectionStage::Measure(2), VolumeProjectionStage::Publish] {
            let origin = GridSdf::from_fn(8, &|_, y| y - 0.5);
            let mut phi = origin.clone();
            let result = project_material_volume_controlled(&mut phi, 3, &[], controls(0.25), |stage| {
                if stage == stop { ControlFlow::Break(stage) } else { ControlFlow::Continue(()) }
            }).unwrap();
            assert!(matches!(result, ControlFlow::Break(stage) if stage == stop));
            assert_eq!(bits(&phi), bits(&origin));
        }
    }
}
