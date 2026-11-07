//! Conserved original-reference density through the existing liquid handoff.
//!
//! A current affine tetrahedron carries rho_reference/J. Native VOF keeps its
//! own uniform physical density: the original source is never rewritten as a
//! numerical pressure or a different geometric fill.

use super::{
    LiquidTransferMoments, TetLiquidTransferError, TetrahedralLiquidState, add, cross, dot,
    field_gradient, invalid, scale, sub,
};
use fs_exec::Cx;
use fs_mesh::overlap::{DensityMoments, tetrahedron_box_density_moments};

/// Original material inventory on the SAME tetrahedral connectivity.
///
/// Reference positions may differ from the current liquid geometry. Each
/// tetrahedron retains its four positive reference densities in its original
/// corner order; adjacent elements may have different original inventories.
/// The supplied current geometry and velocity must already describe liquid
/// without stored elastic/plastic history. This field does not authorize a
/// phase change or discard such history.
#[derive(Clone, Copy, Debug)]
pub struct TetrahedralLiquidReferenceDensity<'a> {
    pub positions_m: &'a [[f64; 3]],
    /// Original mass per reference volume [kg/m3], affine within each tet.
    pub corner_density_kg_m3: &'a [[f64; 4]],
}

/// Explicit physical accuracy for the reference-to-native density transition.
#[derive(Clone, Copy, Debug)]
pub struct TetReferenceDensityTransferControls {
    /// Maximum absolute difference between ANY current rho_reference/J corner
    /// and the destination's original physical density [kg/m3]. Zero requires
    /// exact represented agreement. A nonzero allowance is reported, never
    /// converted into pressure or silently repaired by changing fill.
    pub maximum_density_defect_kg_m3: f64,
    /// Sum of absolute per-cell native occupied volume minus geometric
    /// intersection volume [m3], including actual population readback.
    pub maximum_volume_defect_m3: f64,
}

/// Original and deformed density inventories, about the report's grid origin.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TetReferenceDensityTransferReport {
    /// Exact affine-density reference integrals before any deformation.
    pub reference: DensityMoments,
    /// Affine push-forward of the SAME mass and first/second moments.
    /// Its volume is the actual current geometry, not mass/native_density.
    pub deformed: DensityMoments,
    pub maximum_current_density_defect_kg_m3: f64,
    /// Sum over tets of |clipped current mass - original reference mass| [kg].
    pub overlap_mass_partition_defect_kg: f64,
    /// Sum over tets of absolute geometric overlap partition defects [m3].
    pub overlap_volume_partition_defect_m3: f64,
    /// Sum over receiving cells of |actual native mass - overlap mass| [kg].
    /// Together with the partition defect this must fit the existing mass bound.
    pub sum_absolute_cell_mass_defect_kg: f64,
    /// Actual native occupied geometry: sum(fill_fraction * cell_size^3).
    /// Numerical pressure-density departures are audited by mass separately.
    pub native_occupied_volume_m3: f64,
    pub native_minus_deformed_volume_m3: f64,
    pub sum_absolute_cell_volume_defect_m3: f64,
}

#[derive(Clone, Copy)]
pub(super) struct Binding<'a> {
    pub field: TetrahedralLiquidReferenceDensity<'a>,
    pub controls: TetReferenceDensityTransferControls,
}

impl Binding<'_> {
    pub(super) fn admit(
        self,
        source: TetrahedralLiquidState<'_>,
    ) -> Result<(), TetLiquidTransferError> {
        if self.field.positions_m.len() != source.positions_m.len()
            || self.field.corner_density_kg_m3.len() != source.tetrahedra.len()
            || [
                self.controls.maximum_density_defect_kg_m3,
                self.controls.maximum_volume_defect_m3,
            ]
            .into_iter()
            .any(|value| !value.is_finite() || value < 0.0)
        {
            return Err(invalid(
                "invalid original reference-density field or physical accuracy bounds",
            ));
        }
        Ok(())
    }

    pub(super) fn tetrahedron(
        self,
        index: usize,
        indices: [u32; 4],
        current: [[f64; 3]; 4],
        current_volume: f64,
        native_density: f64,
        cx: &Cx<'_>,
    ) -> Result<Tetrahedron, TetLiquidTransferError> {
        let reference = indices.map(|i| self.field.positions_m[i as usize]);
        let density = self.field.corner_density_kg_m3[index];
        if reference.iter().flatten().any(|x| !x.is_finite()) {
            return Err(invalid("nonfinite original reference geometry"));
        }
        let lower = core::array::from_fn(|axis| {
            reference
                .iter()
                .map(|p| p[axis])
                .fold(f64::INFINITY, f64::min)
        });
        let upper = core::array::from_fn(|axis| {
            reference
                .iter()
                .map(|p| p[axis])
                .fold(f64::NEG_INFINITY, f64::max)
        });
        let original =
            tetrahedron_box_density_moments(reference, density, lower, upper, reference[0], cx)?;
        if !original.volume.is_finite() || original.volume <= 0.0 {
            return Err(invalid(
                "original reference tetrahedron has no positive volume",
            ));
        }
        let deformation = field_gradient(reference, current)?;
        let magnitude = deformation
            .iter()
            .flatten()
            .fold(0.0_f64, |a, b| a.max(b.abs()));
        if !magnitude.is_finite() || magnitude <= 0.0 {
            return Err(invalid("unrepresentable reference-to-current deformation"));
        }
        let normalized = deformation.map(|row| row.map(|value| value / magnitude));
        let signed = dot(normalized[0], cross(normalized[1], normalized[2]));
        if !signed.is_finite() || signed <= 0.0 {
            return Err(invalid(
                "reference-to-current tetrahedron is inverted or singular",
            ));
        }
        let jacobian = current_volume / original.volume;
        if !jacobian.is_finite() || jacobian <= 0.0 {
            return Err(invalid("unrepresentable current/reference volume ratio"));
        }
        let current_density = density.map(|rho| rho / jacobian);
        let maximum_defect = current_density
            .iter()
            .map(|rho| (rho - native_density).abs())
            .fold(0.0_f64, f64::max);
        if current_density
            .iter()
            .any(|rho| !rho.is_finite() || *rho <= 0.0)
            || !maximum_defect.is_finite()
            || maximum_defect > self.controls.maximum_density_defect_kg_m3
        {
            return Err(invalid(
                "deformed reference density is incompatible with the native liquid density; pressure and fill cannot replace physical density",
            ));
        }
        let deformed = transform(original, deformation, [0.0; 3], current_volume)?;
        Ok(Tetrahedron {
            original,
            deformed,
            reference_origin: reference[0],
            current_origin: current[0],
            reference_density: density,
            current_density,
            maximum_defect,
        })
    }
}

pub(super) struct Tetrahedron {
    pub original: DensityMoments,
    /// About current corner zero; the mass is the original reference integral.
    pub deformed: DensityMoments,
    pub reference_origin: [f64; 3],
    pub current_origin: [f64; 3],
    pub reference_density: [f64; 4],
    pub current_density: [f64; 4],
    pub maximum_defect: f64,
}

impl Tetrahedron {
    pub(super) fn add_report(
        &self,
        report: &mut TetReferenceDensityTransferReport,
        grid_origin: [f64; 3],
    ) -> Result<(), TetLiquidTransferError> {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        accumulate(
            &mut report.reference,
            transform(
                self.original,
                identity,
                sub(self.reference_origin, grid_origin),
                self.original.volume,
            )?,
        )?;
        accumulate(
            &mut report.deformed,
            transform(
                self.deformed,
                identity,
                sub(self.current_origin, grid_origin),
                self.deformed.volume,
            )?,
        )?;
        report.maximum_current_density_defect_kg_m3 = report
            .maximum_current_density_defect_kg_m3
            .max(self.maximum_defect);
        Ok(())
    }

    pub(super) fn add_source(
        &self,
        total: &mut LiquidTransferMoments,
        velocity: [f64; 3],
        gradient: [[f64; 3]; 3],
        grid_origin: [f64; 3],
    ) {
        let m = self.deformed;
        let offset = sub(self.current_origin, grid_origin);
        let gradient_first = gradient.map(|row| dot(row, m.first));
        let momentum = add(scale(velocity, m.mass), gradient_first);
        let mut angular = add(cross(offset, momentum), cross(m.first, velocity));
        for axis in 0..3 {
            let b = (axis + 1) % 3;
            let c = (axis + 2) % 3;
            for k in 0..3 {
                angular[axis] += m.second[b][k] * gradient[c][k] - m.second[c][k] * gradient[b][k];
            }
        }
        let quadratic: f64 = gradient
            .into_iter()
            .map(|row| {
                (0..3)
                    .map(|i| {
                        (0..3)
                            .map(|j| row[i] * m.second[i][j] * row[j])
                            .sum::<f64>()
                    })
                    .sum::<f64>()
            })
            .sum();
        total.mass_kg += m.mass;
        total.first_mass_moment_kg_m = add(
            total.first_mass_moment_kg_m,
            add(m.first, scale(offset, m.mass)),
        );
        total.momentum_kg_m_s = add(total.momentum_kg_m_s, momentum);
        total.angular_momentum_kg_m2_s = add(total.angular_momentum_kg_m2_s, angular);
        total.kinetic_energy_j += 0.5 * m.mass * dot(velocity, velocity)
            + dot(velocity, gradient_first)
            + 0.5 * quadratic;
    }

    /// Original P1 enthalpy against original P1 density: a quadratic integral.
    pub(super) fn specific_field_integral(&self, values: [f64; 4]) -> f64 {
        let scale = self.reference_density.into_iter().fold(0.0_f64, f64::max);
        let weights = self.reference_density.map(|rho| rho / scale);
        let sum = weights.into_iter().sum::<f64>();
        let total = values.into_iter().sum::<f64>();
        self.original.mass / 5.0
            * (total + (0..4).map(|i| weights[i] / sum * values[i]).sum::<f64>())
    }
}

/// Push forward mass integrals; density is divided by J, not multiplied by it.
fn transform(
    original: DensityMoments,
    map: [[f64; 3]; 3],
    shift: [f64; 3],
    volume: f64,
) -> Result<DensityMoments, TetLiquidTransferError> {
    let first = map.map(|row| dot(row, original.first));
    let mut second = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in i..3 {
            let mut value =
                original.mass * shift[i] * shift[j] + shift[i] * first[j] + shift[j] * first[i];
            for a in 0..3 {
                for b in 0..3 {
                    value += map[i][a] * original.second[a][b] * map[j][b];
                }
            }
            second[i][j] = value;
            second[j][i] = value;
        }
    }
    let result = DensityMoments {
        volume,
        mass: original.mass,
        first: add(first, scale(shift, original.mass)),
        second,
    };
    finite(result)
        .then_some(result)
        .ok_or_else(|| invalid("unrepresentable transported original density moments"))
}

fn finite(moments: DensityMoments) -> bool {
    moments.volume.is_finite()
        && moments.volume > 0.0
        && moments.mass.is_finite()
        && moments.mass > 0.0
        && moments
            .first
            .iter()
            .chain(moments.second.iter().flatten())
            .all(|v| v.is_finite())
}

fn accumulate(
    total: &mut DensityMoments,
    value: DensityMoments,
) -> Result<(), TetLiquidTransferError> {
    total.volume += value.volume;
    total.mass += value.mass;
    total.first = add(total.first, value.first);
    for i in 0..3 {
        for j in 0..3 {
            total.second[i][j] += value.second[i][j];
        }
    }
    if !finite(*total) {
        return Err(invalid(
            "unrepresentable combined original density inventory",
        ));
    }
    Ok(())
}

/// Repeat the exact same density/volume audit after thermal stress attachment.
/// Native initialization can round populations again; report the final cells.
pub(super) fn audit_native(
    flow: &super::MaterialFreeSurface3<'_>,
    report: &mut TetReferenceDensityTransferReport,
    cells: &std::collections::BTreeMap<[usize; 3], (f64, f64)>,
    controls: TetReferenceDensityTransferControls,
    mass_tolerance: f64,
    cx: &Cx<'_>,
) -> Result<(), TetLiquidTransferError> {
    let dx = flow.scaling.cell_size_m();
    let cell_volume = dx * dx * dx;
    report.native_occupied_volume_m3 = 0.0;
    report.sum_absolute_cell_mass_defect_kg = 0.0;
    report.sum_absolute_cell_volume_defect_m3 = 0.0;
    for (&coordinates, &(mass, volume)) in cells {
        super::checkpoint(cx)?;
        let native_mass = super::native_cell_mass_kg(flow, coordinates)?;
        let [x, y, z] = coordinates;
        let sample = flow
            .sample_at_si(x, y, z)?
            .ok_or_else(|| invalid("final reference-density readback lost an occupied cell"))?;
        let native_volume = sample.fill_fraction * cell_volume;
        if !native_volume.is_finite() || native_volume <= 0.0 {
            return Err(invalid("final native occupied geometry is unrepresentable"));
        }
        report.native_occupied_volume_m3 += native_volume;
        report.sum_absolute_cell_mass_defect_kg += (native_mass - mass).abs();
        report.sum_absolute_cell_volume_defect_m3 += (native_volume - volume).abs();
    }
    report.native_minus_deformed_volume_m3 =
        report.native_occupied_volume_m3 - report.deformed.volume;
    if ![
        report.native_occupied_volume_m3,
        report.native_minus_deformed_volume_m3,
        report.sum_absolute_cell_mass_defect_kg,
        report.sum_absolute_cell_volume_defect_m3,
    ]
    .into_iter()
    .all(f64::is_finite)
        || report.sum_absolute_cell_mass_defect_kg + report.overlap_mass_partition_defect_kg
            > mass_tolerance
        || report.sum_absolute_cell_volume_defect_m3 + report.overlap_volume_partition_defect_m3
            > controls.maximum_volume_defect_m3
    {
        return Err(invalid(
            "final native density cannot preserve the original mass and geometric volume within their explicit bounds",
        ));
    }
    Ok(())
}
