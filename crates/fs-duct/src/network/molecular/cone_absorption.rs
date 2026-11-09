//! Positive loss forms of the ORIGINAL frozen-radius spherical wave slices.
//!
//! With r=r0+s*x, q=r*p/r_ref and Z_ref=z_specific/A_ref, the native
//! equations give U*r_ref/r = [(1-s/(i*k*r))*F*exp(i*k*x)
//! -(1+s/(i*k*r))*B*exp(i*k*(L-x))]/Z_ref. This spherical induction
//! term is essential to viscous work. Pressure-weighted losses reduce exactly
//! to the shared positive q integral; no signed power residual is called heat.

use super::absorption::physical_channels;
use super::*;
use crate::compose::pressure_wave_integral_weights;
use crate::molecular::local::absorption::MolecularCone;
use crate::network::field::absorption::{
    CylinderLossOperator, admit_cells, divide, multiply, norm_squared, observe_loss_forms,
};
use crate::network::field::{MAX_FIELD_PROBES, NetworkCylinderAbsorptionCell, NetworkFieldSection};
use fs_math::{c64::C64, dd::Dd};

#[derive(Clone, Copy, Debug)]
pub(super) struct ConeSlice<'source> {
    pub(super) start_m: f64,
    pub(super) end_m: f64,
    pub(super) native: MolecularCone<'source>,
}

/// Work and refinement limits for positive spherical-induction integration.
#[derive(Clone, Copy, Debug)]
pub struct ConeAbsorptionControls {
    pub maximum_cells: usize,
    /// Per-cell composite four-node Gauss panels, at most 8192.
    pub maximum_panels_per_cell: usize,
    /// Aggregate node evaluations, including discarded coarse rules, at most
    /// 4,194,304. The same bound covers all cells in this preparation.
    pub maximum_quadrature_evaluations: usize,
    /// Successive viscous Gram-entry change / largest entry. This is a
    /// refinement estimate at matrix scale, not a relative guarantee for a
    /// nearly cancelling excitation or a continuum/material error bound.
    pub relative_tolerance: f64,
}

impl Default for ConeAbsorptionControls {
    fn default() -> Self {
        Self {
            maximum_cells: MAX_FIELD_PROBES,
            maximum_panels_per_cell: 4096,
            maximum_quadrature_evaluations: 1_048_576,
            relative_tolerance: 1e-8,
        }
    }
}

/// One physical frustum contained in a single ORIGINAL native cone slice.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularNetworkConeAbsorptionCell {
    pub probe_indices: [usize; 2],
    pub segment: usize,
    pub slice: usize,
    /// Coordinates relative to the authored segment inlet.
    pub start_m: f64,
    pub end_m: f64,
    pub inlet_radius_m: f64,
    pub outlet_radius_m: f64,
    pub volume_m3: f64,
}

/// A posteriori change of the viscous form, independently of the source drive.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ConeAbsorptionQuadrature {
    pub panels: usize,
    pub evaluations: usize,
    /// Twice the largest absolute entry change [m]. Multiplying this by
    /// |F/Z_ref|²+|B/Z_ref|² bounds the CHANGE between the last two integrated
    /// velocity forms. It is an estimate of quadrature error, not an enclosure.
    pub operator_change_m: f64,
    pub relative_change: f64,
}

/// Prepared coherent gas-loss and separate mechanical-wall-work forms.
/// Warm observations use constant-size positive factors and original source
/// contractions only; there are no new graph, material or radial solves.
#[derive(Debug)]
pub struct PreparedMolecularNetworkConeAbsorption<'fields, 'network> {
    fields: &'fields PreparedMolecularNetworkFields<'network>,
    cells: Vec<MolecularNetworkConeAbsorptionCell>,
    operators: Vec<ConeOperator>,
    quadrature: Vec<ConeAbsorptionQuadrature>,
}

#[derive(Clone, Copy, Debug)]
enum ConeOperator {
    Zero,
    Uniform(CylinderLossOperator<4>),
    Spherical {
        impedance: C64,
        end_impedance: [C64; 2],
        induction: [C64; 2],
        radius_ratio: [f64; 2],
        pressure_weights: [f64; 2],
        losses: [f64; 4],
        viscous: PositiveGram,
    },
}

impl<'network> PreparedMolecularNetworkFields<'network> {
    /// Integrate physical losses between ordered probes inside one original
    /// spherical cone slice. Whole cones are covered by cells at their native
    /// slice cuts; a cell crossing an unprobed cut refuses. Equal endpoints
    /// give exact zero. Overlapping cells are never automatically summed.
    ///
    /// The native frozen midpoint Cp, diffusion spectra and mechanical wall
    /// admittance are retained. Area normalization varies exactly as in its
    /// spherical propagator; the radial/material approximation is not refined
    /// here. Fourier/molecular/wall pressure forms are analytic. Viscosity
    /// retains spherical induction through a positive, bounded Gauss form.
    /// Only viscosity is declared direct caloric heat; availability losses and
    /// mechanical wall work retain the cylinder channel carrier's meanings.
    ///
    /// # Errors
    /// Original identity/order/slice, invalid controls, work exhaustion,
    /// unresolved quadrature, unrepresentable positive form or cancellation.
    /// Profiles require a native within-Magnus-panel field law and refuse.
    pub fn prepare_cone_absorption(
        &self,
        probe_pairs: &[[usize; 2]],
        controls: ConeAbsorptionControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularNetworkConeAbsorption<'_, 'network>, DuctError> {
        poll(&mut cancelled, 0)?;
        admit_cells(probe_pairs.len(), controls.maximum_cells)?;
        if controls.maximum_panels_per_cell < 2
            || controls.maximum_panels_per_cell > 8192
            || controls.maximum_quadrature_evaluations > 4_194_304
            || !controls.relative_tolerance.is_finite()
            || !(128.0 * f64::EPSILON..=0.1).contains(&controls.relative_tolerance)
        {
            return Err(bad(
                "invalid spherical absorption work or refinement controls",
            ));
        }
        let mut result = PreparedMolecularNetworkConeAbsorption {
            fields: self,
            cells: Vec::with_capacity(probe_pairs.len()),
            operators: Vec::with_capacity(probe_pairs.len()),
            quadrature: Vec::with_capacity(probe_pairs.len()),
        };
        let rule = crate::dg::gauss_legendre_rule(4);
        let mut evaluations = 0;
        let mut coefficients = BTreeMap::new();
        for (index, &pair) in probe_pairs.iter().enumerate() {
            poll(&mut cancelled, index)?;
            let (cell, source, length) = self.cone_cell(pair)?;
            let mut report = ConeAbsorptionQuadrature::default();
            let operator = if length == 0.0 {
                ConeOperator::Zero
            } else {
                let edge = self.fields().locations()[pair[0]].probe.edge_id.as_str();
                let key = (edge, cell.segment, cell.slice);
                let losses = if let Some(&value) = coefficients.get(&key) {
                    value
                } else {
                    let gas = source.original.loss_weights()?;
                    let value = [
                        gas[0],
                        gas[1],
                        gas[2],
                        source.original.mechanical_admittance_per_metre.re,
                    ];
                    if value.iter().any(|v| !v.is_finite() || *v < 0.0) {
                        return Err(bad(
                            "spherical absorption needs passive original loss forms",
                        ));
                    }
                    coefficients.insert(key, value);
                    value
                };
                let wave = source.original.wave;
                let r_ref = source.original.radius_m;
                let area_ref = core::f64::consts::PI * r_ref * r_ref;
                if (source.outlet_radius_m - source.inlet_radius_m).abs()
                    < 1e-12 * source.inlet_radius_m
                {
                    // Match the original TwoPort's exact cylindrical branch.
                    let ratio = (source.inlet_radius_m / r_ref).powi(2);
                    ConeOperator::Uniform(CylinderLossOperator::new(
                        wave.wavenumber,
                        wave.specific_impedance.scale(1.0 / (area_ref * ratio)),
                        length,
                        [
                            losses[0] / ratio,
                            losses[1] * ratio,
                            losses[2] * ratio,
                            losses[3] * ratio,
                        ],
                    )?)
                } else {
                    let radii = [cell.inlet_radius_m, cell.outlet_radius_m];
                    let slope = (source.outlet_radius_m - source.inlet_radius_m) / source.length_m;
                    let inverse_ik = (C64::new(0.0, 1.0) * wave.wavenumber).recip();
                    let induction = radii.map(|r| inverse_ik.scale(slope / r));
                    let (plus, minus) = pressure_wave_integral_weights(wave.wavenumber, length)?;
                    let viscous = if losses[0] == 0.0 {
                        PositiveGram::ZERO
                    } else {
                        let (form, observed) = prepare_viscous(
                            wave.wavenumber,
                            length,
                            radii,
                            slope,
                            controls,
                            &rule,
                            &mut evaluations,
                            &mut cancelled,
                        )?;
                        report = observed;
                        form
                    };
                    ConeOperator::Spherical {
                        impedance: wave.specific_impedance.scale(1.0 / area_ref),
                        end_impedance: radii.map(|r| {
                            wave.specific_impedance
                                .scale(1.0 / (core::f64::consts::PI * r * r))
                        }),
                        induction,
                        radius_ratio: radii.map(|r| r / r_ref),
                        pressure_weights: [plus, minus],
                        losses,
                        viscous,
                    }
                }
            };
            result.cells.push(cell);
            result.operators.push(operator);
            result.quadrature.push(report);
        }
        poll(&mut cancelled, result.cells.len())?;
        Ok(result)
    }

    fn cone_cell(
        &self,
        pair: [usize; 2],
    ) -> Result<(MolecularNetworkConeAbsorptionCell, MolecularCone<'_>, f64), DuctError> {
        let coordinate = |index: usize| -> Result<(&str, usize, f64), DuctError> {
            let location = self
                .fields()
                .locations()
                .get(index)
                .ok_or_else(|| bad("unknown spherical absorption probe"))?;
            let segment = match location.probe.section {
                NetworkFieldSection::Inlet { segment }
                | NetworkFieldSection::Outlet { segment }
                | NetworkFieldSection::Cone { segment, .. } => segment,
                _ => return Err(bad("spherical absorption needs an original cone probe")),
            };
            let edge = self
                .network()
                .network()
                .edges
                .iter()
                .find(|e| e.id == location.probe.edge_id)
                .expect("original field edge exists");
            let Some(Segment::Cone { length, .. }) = edge.duct.segments.get(segment) else {
                return Err(bad("spherical absorption needs an authored cone"));
            };
            let x = match location.probe.section {
                NetworkFieldSection::Inlet { .. } => 0.0,
                NetworkFieldSection::Outlet { .. } => *length,
                NetworkFieldSection::Cone { offset_m, .. } => offset_m,
                _ => unreachable!(),
            };
            Ok((&location.probe.edge_id, segment, x))
        };
        let (edge, segment, start_m) = coordinate(pair[0])?;
        let (other, last, end_m) = coordinate(pair[1])?;
        if edge != other || segment != last || start_m > end_m {
            return Err(bad(
                "spherical cell endpoints must be ordered in one original cone",
            ));
        }
        let retained = self
            .cones
            .get(edge)
            .ok_or_else(|| bad("spherical absorption needs retained original slices"))?;
        let (&(_, slice), source) = retained
            .range((segment, 0)..=(segment, usize::MAX))
            .find(|(_, source)| start_m >= source.start_m && end_m <= source.end_m)
            .ok_or_else(|| {
                bad("spherical cell crosses a native slice; include its original cut probes")
            })?;
        let local = |x: f64| (x - source.start_m) / (source.end_m - source.start_m);
        let a = local(start_m);
        let b = local(end_m);
        let native = source.native;
        let radius =
            |t: f64| native.inlet_radius_m + (native.outlet_radius_m - native.inlet_radius_m) * t;
        let (r0, r1) = (radius(a), radius(b));
        let length = native.length_m * (b - a);
        let volume_m3 = (core::f64::consts::PI / 3.0) * length * (r0 * r0 + r0 * r1 + r1 * r1);
        if !volume_m3.is_finite() || (end_m > start_m && (volume_m3 <= 0.0 || length <= 0.0)) {
            return Err(bad("spherical cell volume is unrepresentable"));
        }
        Ok((
            MolecularNetworkConeAbsorptionCell {
                probe_indices: pair,
                segment,
                slice,
                start_m,
                end_m,
                inlet_radius_m: r0,
                outlet_radius_m: r1,
                volume_m3,
            },
            native,
            length,
        ))
    }
}

impl<'fields, 'network> PreparedMolecularNetworkConeAbsorption<'fields, 'network> {
    pub fn fields(&self) -> &'fields PreparedMolecularNetworkFields<'network> {
        self.fields
    }
    pub fn cells(&self) -> &[MolecularNetworkConeAbsorptionCell] {
        &self.cells
    }
    pub fn quadrature(&self) -> &[ConeAbsorptionQuadrature] {
        &self.quadrature
    }
    pub fn sample(
        &self,
        excitation: &[C64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<MolecularNetworkCylinderAbsorption>, DuctError> {
        self.observe_coherence(Some(excitation), &[], cancelled)
    }

    /// Mean and covariance factors use the same canonical source coordinates
    /// and phase-preserving contractions as original cylindrical losses.
    pub fn observe_coherence(
        &self,
        mean: Option<&[C64]>,
        factors: &[&[C64]],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<MolecularNetworkCylinderAbsorption>, DuctError> {
        observe_loss_forms(
            self.fields.fields(),
            self.cells.len(),
            mean,
            factors,
            cancelled,
            |index, excitation| {
                self.operators[index].channels(self.fields.fields(), &self.cells[index], excitation)
            },
            physical_channels,
        )
    }
}

impl ConeOperator {
    fn channels(
        &self,
        fields: &crate::network::field::PreparedNetworkFields,
        cell: &MolecularNetworkConeAbsorptionCell,
        excitation: &[C64],
    ) -> Result<[Dd; 4], DuctError> {
        let Self::Spherical {
            impedance,
            end_impedance,
            induction,
            radius_ratio,
            pressure_weights: [plus, minus],
            losses,
            viscous,
        } = self
        else {
            return match self {
                Self::Zero => Ok([Dd::ZERO; 4]),
                Self::Uniform(operator) => operator.channels(
                    fields,
                    &NetworkCylinderAbsorptionCell {
                        probe_indices: cell.probe_indices,
                        segment: cell.segment,
                        start_m: cell.start_m,
                        end_m: cell.end_m,
                        volume_m3: cell.volume_m3,
                    },
                    excitation,
                ),
                Self::Spherical { .. } => unreachable!(),
            };
        };
        let [p0, u0] = fields.parts(cell.probe_indices[0], excitation)?;
        let [p1, u1] = fields.parts(cell.probe_indices[1], excitation)?;
        let u0 = multiply(u0, end_impedance[0]);
        let u1 = multiply(u1, end_impedance[1]);
        let a0 = multiply(p0, induction[0]);
        let a1 = multiply(p1, induction[1]);
        // Forward at the inlet, backward at the outlet: no growing travel.
        let f =
            std::array::from_fn(|i| (p0[i] + u0[i] + a0[i]) * Dd::from_f64(0.5 * radius_ratio[0]));
        let b =
            std::array::from_fn(|i| (p1[i] - u1[i] - a1[i]) * Dd::from_f64(0.5 * radius_ratio[1]));
        let sum = [f[0] + b[0], f[1] + b[1]];
        let difference = [f[0] - b[0], f[1] - b[1]];
        let pressure = Dd::from_f64(*plus) * norm_squared(sum)
            + Dd::from_f64(*minus) * norm_squared(difference);
        let mut watts = [Dd::ZERO; 4];
        if losses[0] != 0.0 {
            watts[0] = viscous.evaluate(divide(f, *impedance), divide(b, *impedance))
                * Dd::from_f64(0.5 * losses[0]);
        }
        for channel in 1..4 {
            watts[channel] = pressure * Dd::from_f64(0.5 * losses[channel]);
        }
        Ok(watts)
    }
}

// Positive Gauss Gram, then DD Cholesky: constant warm cost and no negative
// cancellation of large signed work. Zero/rank-one forms remain exact limits.
#[derive(Clone, Copy, Debug)]
struct PositiveGram {
    swapped: bool,
    root: Dd,
    cross: [Dd; 2],
    residual: Dd,
}
impl PositiveGram {
    const ZERO: Self = Self {
        swapped: false,
        root: Dd::ZERO,
        cross: [Dd::ZERO; 2],
        residual: Dd::ZERO,
    };
    fn new(matrix: [Dd; 4]) -> Result<Self, DuctError> {
        let [mut a, mut d, re, mut im] = matrix;
        let swapped = a.lt(d);
        if swapped {
            std::mem::swap(&mut a, &mut d);
            im = -im;
        }
        if a == Dd::ZERO && d == Dd::ZERO && re == Dd::ZERO && im == Dd::ZERO {
            return Ok(Self::ZERO);
        }
        let root = a.sqrt();
        let cross = [re / root, im / root];
        let residual = d - norm_squared(cross);
        if matrix
            .iter()
            .chain([root, cross[0], cross[1], residual].iter())
            .any(|v| !v.hi.is_finite() || !v.lo.is_finite())
            || residual.lt(Dd::ZERO)
        {
            return Err(bad("positive spherical viscous factor is unrepresentable"));
        }
        Ok(Self {
            swapped,
            root,
            cross,
            residual,
        })
    }
    fn evaluate(&self, mut f: [Dd; 2], mut b: [Dd; 2]) -> Dd {
        if self.swapped {
            std::mem::swap(&mut f, &mut b);
        }
        let value = [
            self.root * f[0] + self.cross[0] * b[0] - self.cross[1] * b[1],
            self.root * f[1] + self.cross[0] * b[1] + self.cross[1] * b[0],
        ];
        norm_squared(value) + self.residual * norm_squared(b)
    }
}

fn prepare_viscous(
    k: C64,
    length: f64,
    radii: [f64; 2],
    slope: f64,
    controls: ConeAbsorptionControls,
    rule: &(Vec<f64>, Vec<f64>),
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(PositiveGram, ConeAbsorptionQuadrature), DuctError> {
    // Do not accept an aliased oscillation or a missed opaque endpoint layer
    // merely because two coarse rules both happen to sample zero.
    let needed = (k.abs() * length).max(2.0 * (radii[1] - radii[0]).abs() / radii[0].min(radii[1]));
    if !needed.is_finite() || needed > controls.maximum_panels_per_cell as f64 {
        return Err(bad(
            "spherical induction resolution exceeds its panel budget",
        ));
    }
    let mut panels = (needed.ceil() as usize).max(1);
    let initial_evaluations = *evaluations;
    let inverse_ik = (C64::new(0.0, 1.0) * k).recip();
    let mut integrate = |panels: usize| -> Result<[Dd; 4], DuctError> {
        let work = panels * rule.0.len();
        if work
            > controls
                .maximum_quadrature_evaluations
                .saturating_sub(*evaluations)
        {
            return Err(bad(
                "spherical absorption quadrature evaluation budget exceeded",
            ));
        }
        *evaluations += work;
        let h = length / panels as f64;
        let mut matrix = [Dd::ZERO; 4];
        for panel in 0..panels {
            poll(cancelled, panel)?;
            for (&node, &weight) in rule.0.iter().zip(&rule.1) {
                let x = h * (panel as f64 + 0.5 * (node + 1.0));
                let radius = radii[0] + (radii[1] - radii[0]) * (x / length);
                let alpha = inverse_ik.scale(slope / radius);
                let f = crate::exp_i(k.scale(x)) * (C64::ONE - alpha);
                let b = crate::exp_i(k.scale(length - x)) * (-C64::ONE - alpha);
                let f = [Dd::from_f64(f.re), Dd::from_f64(f.im)];
                let b = [Dd::from_f64(b.re), Dd::from_f64(b.im)];
                let w = Dd::from_f64(0.5 * h * weight);
                let entries = [
                    norm_squared(f),
                    norm_squared(b),
                    f[0] * b[0] + f[1] * b[1],
                    f[0] * b[1] - f[1] * b[0],
                ];
                for i in 0..4 {
                    matrix[i] = matrix[i] + w * entries[i];
                }
            }
        }
        if matrix
            .iter()
            .any(|v| !v.hi.is_finite() || !v.lo.is_finite())
        {
            return Err(bad("spherical viscous integral is unrepresentable"));
        }
        Ok(matrix)
    };
    let mut previous = integrate(panels)?;
    loop {
        panels *= 2;
        if panels > controls.maximum_panels_per_cell {
            return Err(bad(
                "spherical viscous quadrature did not converge within its panel budget",
            ));
        }
        let current = integrate(panels)?;
        let entry_maximum = |matrix: [Dd; 4]| {
            matrix[0]
                .abs()
                .to_f64()
                .max(matrix[1].abs().to_f64())
                .max(C64::new(matrix[2].to_f64(), matrix[3].to_f64()).abs())
        };
        // The complex off-diagonal counts as ONE entry: bounding its real
        // and imaginary components separately would miss a sqrt(2) factor.
        let change = entry_maximum(std::array::from_fn(|i| current[i] - previous[i]));
        let scale = entry_maximum(current);
        if scale == 0.0 {
            return Err(bad("positive spherical viscous integral underflowed"));
        }
        if change <= controls.relative_tolerance * scale {
            return Ok((
                PositiveGram::new(current)?,
                ConeAbsorptionQuadrature {
                    panels,
                    evaluations: *evaluations - initial_evaluations,
                    operator_change_m: 2.0 * change,
                    relative_change: change / scale,
                },
            ));
        }
        previous = current;
    }
}
