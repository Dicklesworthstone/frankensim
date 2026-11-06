//! Boundary Fourier/thermal-resistance baths on the original volume timeline.

use std::collections::BTreeMap;

use fs_material::thermal::conduction::phase::PhaseConductivitySample;

use super::*;
use crate::d3q19::Face3;

/// A prescribed thermal bath covering one complete exposed thermal-cell face.
///
/// The original Cartesian geometry owns the area and centre-to-surface distance
/// (half the local width). The external resistance is additional to that cell's
/// own source-conducting half cell. It may represent a declared layer, contact
/// or heat-transfer film, without adding storage or borrowing another material.
/// Fluid-only entries admit fluid receivers adjoining the exterior or an
/// unbound fixed-solid voxel. Implicit conjugate entries also admit explicitly
/// bound finite solids. A bath never replaces a face between two active cells.
/// Each `(cell, face)` may occur once; unlisted faces remain insulated.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalVolumeBathFace {
    /// Original fluid-cell index, or explicitly bound finite-solid cell in an
    /// implicit conjugate entry. An unbound fixed-solid receiver is refused.
    pub cell: usize,
    /// Exposed side in the global Cartesian frame.
    pub face: Face3,
    /// Held external absolute temperature [K], finite and strictly positive.
    /// This is not an incoming-fluid enthalpy or a finite bath's state.
    pub temperature_k: f64,
    /// Finite nonnegative resistance per area [m² K/W]. Zero imposes the
    /// temperature at the actual surface; positive values act in series with
    /// the receiving cell's half cell. Omit the face to prescribe insulation.
    pub resistance_m2_k_w: f64,
}

pub(super) fn prepare<'a>(
    domain: &VoxelDomain,
    baths: &'a [ThermalVolumeBathFace],
    maximum_faces: usize,
    gate: &CancelGate,
) -> Result<BTreeMap<(usize, usize), &'a ThermalVolumeBathFace>, ChtError> {
    prepare_for_active_cells(
        domain,
        baths,
        maximum_faces,
        |cell| domain.is_fluid(cell),
        gate,
    )
}

/// Shared exclusive face admission. The caller owns which thermal inventories
/// exist; geometry alone never makes an unbound solid a bath receiver.
pub(super) fn prepare_for_active_cells<'a>(
    domain: &VoxelDomain,
    baths: &'a [ThermalVolumeBathFace],
    maximum_faces: usize,
    active: impl Fn(usize) -> bool,
    gate: &CancelGate,
) -> Result<BTreeMap<(usize, usize), &'a ThermalVolumeBathFace>, ChtError> {
    if baths.len() > maximum_faces {
        return Err(refuse(
            "conductive face-work budget exhausted by thermal baths",
        ));
    }
    let mut ordered = BTreeMap::new();
    for (index, bath) in baths.iter().enumerate() {
        poll_every(gate, index)?;
        let key = admit_face(domain, bath, &active)?;
        if ordered.insert(key, bath).is_some() {
            return Err(refuse("duplicate thermal bath on one boundary face"));
        }
    }
    Ok(ordered)
}

/// Admit the same physical face for scalar and source-owned boundary laws.
pub(super) fn admit_face(
    domain: &VoxelDomain,
    bath: &ThermalVolumeBathFace,
    active: impl Fn(usize) -> bool,
) -> Result<(usize, usize), ChtError> {
    let side = bath.face as usize;
    if bath.cell >= domain.cell_count()
        || !active(bath.cell)
        || domain.neighbor(bath.cell, side).is_some_and(&active)
    {
        return Err(refuse(
            "thermal bath requires a boundary fluid-cell face or an exposed bound finite-solid face",
        ));
    }
    finite_positive("flow.conductive_bath_temperature", bath.temperature_k)?;
    if !bath.resistance_m2_k_w.is_finite() || bath.resistance_m2_k_w < 0.0 {
        return Err(refuse(
            "thermal bath resistance must be finite and nonnegative",
        ));
    }
    Ok((bath.cell, side))
}

#[derive(Default)]
pub(super) struct Transfers {
    pub(super) heat_j: [f64; 6],
    pub(super) absolute_heat_j: f64,
}

// The SAME piecewise-linear original caloric law supplies its global slope
// over the source-admitted h interval. Looking only at the cell's current
// slope would miss sensible heating after a latent plateau is exhausted.
fn maximum_slope(
    conductivity: &PreparedPhaseConductivity<'_>,
    gate: &CancelGate,
) -> Result<f64, ChtError> {
    let [lo, hi] = conductivity.enthalpy_bounds_j_kg();
    let knots = conductivity.caloric_curve().knots();
    let first = knots
        .partition_point(|k| k.specific_enthalpy_j_kg <= lo)
        .saturating_sub(1)
        .min(knots.len() - 2);
    let end = knots
        .partition_point(|k| k.specific_enthalpy_j_kg < hi)
        .min(knots.len() - 1);
    let mut maximum = 0.0_f64;
    for (index, pair) in knots[first..=end].windows(2).enumerate() {
        poll_every(gate, index)?;
        let change = pair[1].temperature_k - pair[0].temperature_k;
        let slope = change / (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg);
        if !slope.is_finite() || slope < 0.0 || (change > 0.0 && slope == 0.0) {
            return Err(refuse(
                "original caloric bath stability bound is unrepresentable",
            ));
        }
        maximum = maximum.max(slope);
    }
    Ok(maximum)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn apply(
    domain: &VoxelDomain,
    conductivity: &PreparedPhaseConductivity<'_>,
    baths: &BTreeMap<(usize, usize), &ThermalVolumeBathFace>,
    samples: &[Option<PhaseConductivitySample>],
    mass: &[f64],
    dt: f64,
    internal_faces: usize,
    maximum_faces: usize,
    increments: &mut [f64],
    mixing_weights: &mut [Sum],
    assembly_defect: &mut Sum,
    gate: &CancelGate,
) -> Result<Transfers, ChtError> {
    if baths.is_empty() {
        return Ok(Transfers::default());
    }
    if internal_faces
        .checked_add(baths.len())
        .is_none_or(|n| n > maximum_faces)
    {
        return Err(refuse(
            "conductive face-work budget exhausted by thermal baths",
        ));
    }
    let slope = maximum_slope(conductivity, gate)?;
    let mut heat = [Sum::default(); 6];
    let mut absolute = Sum::default();
    for (index, (&(cell, side), bath)) in baths.iter().enumerate() {
        poll_every(gate, index)?;
        let sample = samples[cell].expect("admitted fluid bath and original phase state");
        let k = sample.conductivity_w_m_k();
        let delta_t = bath.temperature_k - sample.temperature_k();
        if k == 0.0 {
            continue;
        }
        let axis = side / 2;
        let bulk = (0.5 * domain.widths(cell)[axis]) / k;
        finite_positive("flow.conductive_bath_half_resistance", bulk)?;
        let resistance = bath.resistance_m2_k_w;
        let largest = bulk.max(resistance);
        let g = (domain.face_area(cell, axis) / largest) / (1.0 + bulk.min(resistance) / largest);
        finite_positive("flow.conductive_bath_conductance", g)?;
        if slope > 0.0 {
            let weight = (dt * g / mass[cell]) * slope;
            finite_positive("flow.conductive_bath_enthalpy_weight", weight)?;
            mixing_weights[cell].add(weight);
            if !mixing_weights[cell].value.is_finite() || mixing_weights[cell].value > 1.0 {
                return Err(refuse(
                    "explicit conductive enthalpy mixing weight exceeds one",
                ));
            }
        }
        // An initially equilibrated bath still contributes stiffness to the
        // explicit interval, including when other authored heat changes h.
        if delta_t == 0.0 {
            continue;
        }
        let q = dt * (g * delta_t);
        finite("flow.conductive_bath_heat", q)?;
        if q == 0.0 {
            return Err(refuse("nonzero conductive bath heat underflowed"));
        }
        let (sum, residual) = accounted_add(increments[cell], q)?;
        increments[cell] = sum;
        assembly_defect.add(residual.abs());
        heat[side].add(q);
        absolute.add(q.abs());
        finite("flow.conductive_bath_side_heat", heat[side].value)?;
        finite("flow.conductive_bath_absolute_heat", absolute.value)?;
    }
    Ok(Transfers {
        heat_j: heat.map(|sum| sum.value),
        absolute_heat_j: absolute.value,
    })
}
