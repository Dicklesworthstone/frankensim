//! Positive cycle-mean viscous and thermal work in original cylinders.
//!
//! In exp(-i omega t), p'=i k Zc U and U'=i k p/Zc. Therefore
//! -d Re(p U*)/dx / 2 = [Im(k Zc)|U|^2 + Im(k/Zc)|p|^2]/2.
//! The two weights are read directly from the original Bessel series/shunt
//! law before its square roots, preserving small losses and their mechanisms.
//! Exact positive wave integrals retain coherent pressure and velocity nodes;
//! finite-window boundary work is never used as a heat source.

use super::*;
use crate::WallDiffusionInputs;
use crate::compose::pressure_wave_integral_weights;

/// One bounded physical interval in an original authored cylinder.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkCylinderAbsorptionCell {
    /// Indices into the borrowed field owner's original probe list, in order.
    pub probe_indices: [usize; 2],
    pub segment: usize,
    /// Offsets from THIS segment's inlet, not a sum of preceding lengths.
    pub start_m: f64,
    pub end_m: f64,
    pub volume_m3: f64,
}

/// Irreversible cycle-mean acoustic power in one complete cross-section cell.
/// The channel split does not specify a wall/fluid temperature or heat partition.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NetworkCylinderAbsorption {
    pub viscous_power_w: f64,
    pub thermal_power_w: f64,
    pub absorbed_power_w: f64,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct CylinderSection {
    segment: usize,
    offset_m: f64,
    radius_m: f64,
    diffusion: Option<WallDiffusionInputs>,
}

pub(super) fn cylinder_section(
    probe: &NetworkFieldProbe,
    edge: &NetworkEdge,
    diffusion: Option<WallDiffusionInputs>,
) -> Option<CylinderSection> {
    let segment = match probe.section {
        NetworkFieldSection::Inlet { segment }
        | NetworkFieldSection::Outlet { segment }
        | NetworkFieldSection::Cylinder { segment, .. } => segment,
        _ => return None,
    };
    let Some(&Segment::Cylinder { radius, length }) = edge.duct.segments.get(segment) else {
        return None;
    };
    let offset_m = match probe.section {
        NetworkFieldSection::Inlet { .. } => 0.0,
        NetworkFieldSection::Outlet { .. } => length,
        NetworkFieldSection::Cylinder { offset_m, .. } => offset_m,
        _ => unreachable!("only cylinder-compatible sections admitted"),
    };
    Some(CylinderSection {
        segment,
        offset_m,
        radius_m: radius,
        diffusion,
    })
}

#[derive(Clone, Copy, Debug)]
pub(in crate::network) struct CylinderLossOperator<const N: usize> {
    impedance: C64,
    plus_m: f64,
    minus_m: f64,
    // Channel zero weights |U|^2; all remaining channels weight |p|^2.
    weights: [f64; N],
}

type CellOperator = CylinderLossOperator<2>;

impl<const N: usize> CylinderLossOperator<N> {
    pub(in crate::network) const ZERO: Self = Self {
        impedance: C64::ONE,
        plus_m: 0.0,
        minus_m: 0.0,
        weights: [0.0; N],
    };

    pub(in crate::network) fn new(
        k: C64,
        impedance: C64,
        length: f64,
        weights: [f64; N],
    ) -> Result<Self, DuctError> {
        if !finite(impedance)
            || impedance == C64::ZERO
            || weights.into_iter().any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(bad(
                "cylinder absorption needs finite passive original loss forms",
            ));
        }
        let (plus_m, minus_m) = pressure_wave_integral_weights(k, length)?;
        Ok(Self {
            impedance,
            plus_m,
            minus_m,
            weights,
        })
    }
}

pub(in crate::network) fn admit_cells(count: usize, maximum_cells: usize) -> Result<(), DuctError> {
    if count == 0 || count > maximum_cells || maximum_cells > MAX_FIELD_PROBES {
        Err(bad(
            "cylinder absorption needs a bounded nonempty cell list",
        ))
    } else {
        Ok(())
    }
}

/// Positive cell loss forms bound to the SAME solved graph, source coordinates,
/// clock, original material snapshot and rigid-wall propagation model.
/// Preparation holds at most 1024 cell records and constant-size operators;
/// source transfers and bounded per-probe material records remain borrowed.
/// Repeated coherent or factor queries perform contractions only.
#[derive(Clone, Debug)]
pub struct PreparedNetworkCylinderAbsorption<'a> {
    fields: &'a PreparedNetworkFields,
    cells: Vec<NetworkCylinderAbsorptionCell>,
    operators: Vec<CellOperator>,
}

impl PreparedNetworkFields {
    /// Prepare exact positive cycle-mean heating between pairs of existing
    /// probes in the SAME original cylinder. Inlet/outlet probes are accepted.
    /// Cells may overlap; their output is never automatically summed or deposited.
    /// Equal endpoints are exact zero. No quadrature or remeshing is introduced.
    ///
    /// Bessel uses original viscosity, thermal conductivity and actual Cp/Cv.
    /// Lossless produces exact zero. WideTube's first-order k/Z approximation
    /// does not define this independently sourced mechanism split and refuses.
    /// Cones, profile panels, holes, junctions and exterior loads need their own
    /// volume-loss law and are not assigned artificial cylindrical heat.
    ///
    /// # Errors
    /// Empty/excessive cell lists, invalid or reversed probe pairs, different
    /// original cylinders, unsupported losses, unrepresentable weights, or
    /// cancellation. Both cell arrays fit the caller's maximum_cells and the
    /// hard probe cap, including the temporary per-cylinder coefficient cache;
    /// no partial preparation is published.
    pub fn prepare_cylinder_absorption(
        &self,
        probe_pairs: &[[usize; 2]],
        maximum_cells: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedNetworkCylinderAbsorption<'_>, DuctError> {
        poll(&mut cancelled, 0)?;
        let loss = self.loss.ok_or_else(|| {
            bad("cylinder heat needs the solved original constitutive loss forms")
        })?;
        admit_cells(probe_pairs.len(), maximum_cells)?;
        if loss == LossModel::WideTube {
            return Err(bad(
                "cylinder heat requires Bessel or exact lossless propagation",
            ));
        }
        let mut cells = Vec::with_capacity(probe_pairs.len());
        let mut operators = Vec::with_capacity(probe_pairs.len());
        // At most one original constitutive evaluation per selected cylinder;
        // this temporary cache is also bounded by the number of cells.
        let mut constitutive = std::collections::BTreeMap::new();
        for (index, &probe_indices) in probe_pairs.iter().enumerate() {
            poll(&mut cancelled, index)?;
            let cell = self.cylinder_absorption_cell(probe_indices)?;
            let [left, right] = probe_indices;
            let a = self
                .cylinders
                .get(left)
                .and_then(Option::as_ref)
                .ok_or_else(|| bad("cell start is not an original cylinder probe"))?;
            let b = self
                .cylinders
                .get(right)
                .and_then(Option::as_ref)
                .ok_or_else(|| bad("cell end is not an original cylinder probe"))?;
            let length = b.offset_m - a.offset_m;
            let operator = if loss == LossModel::Lossless || length == 0.0 {
                CellOperator::ZERO
            } else {
                let key = (self.locations[left].probe.edge_id.as_str(), a.segment);
                let (k, impedance, viscous_weight, thermal_weight) =
                    if let Some(&cached) = constitutive.get(&key) {
                        cached
                    } else {
                        let inputs = a.diffusion.as_ref().ok_or_else(|| {
                            bad("cylinder heat needs its original diffusion snapshot")
                        })?;
                        let (series, shunt, _) = crate::bessel_telegraph_pair_from_inputs(
                            inputs,
                            a.radius_m,
                            self.ports.omega(),
                        )?;
                        let (k, impedance) = crate::helmholtz_wave_from_telegraph(series, shunt);
                        if !finite(impedance) || impedance == C64::ZERO {
                            return Err(DuctError::Singular);
                        }
                        let cached = (k, impedance, series.re, shunt.re);
                        constitutive.insert(key, cached);
                        cached
                    };
                CellOperator::new(k, impedance, length, [viscous_weight, thermal_weight])?
            };
            cells.push(cell);
            operators.push(operator);
        }
        poll(&mut cancelled, cells.len())?;
        Ok(PreparedNetworkCylinderAbsorption {
            fields: self,
            cells,
            operators,
        })
    }

    /// Geometry admission shared by original equilibrium and molecular loss
    /// owners. This does not select, infer or evaluate a constitutive law.
    pub(in crate::network) fn cylinder_absorption_cell(
        &self,
        probe_indices: [usize; 2],
    ) -> Result<NetworkCylinderAbsorptionCell, DuctError> {
        let [left, right] = probe_indices;
        let a = self
            .cylinders
            .get(left)
            .and_then(Option::as_ref)
            .ok_or_else(|| bad("cell start is not an original cylinder probe"))?;
        let b = self
            .cylinders
            .get(right)
            .and_then(Option::as_ref)
            .ok_or_else(|| bad("cell end is not an original cylinder probe"))?;
        if self.locations[left].probe.edge_id != self.locations[right].probe.edge_id
            || a.segment != b.segment
            || a.offset_m > b.offset_m
        {
            return Err(bad(
                "cell endpoints must be ordered in the same original cylinder",
            ));
        }
        let length = b.offset_m - a.offset_m;
        let volume_m3 = self.locations[left].area_m2 * length;
        if !volume_m3.is_finite() || (length > 0.0 && volume_m3 == 0.0) {
            return Err(bad("cylinder cell volume is unrepresentable"));
        }
        Ok(NetworkCylinderAbsorptionCell {
            probe_indices,
            segment: a.segment,
            start_m: a.offset_m,
            end_m: b.offset_m,
            volume_m3,
        })
    }
}

impl PreparedNetworkCylinderAbsorption<'_> {
    pub fn fields(&self) -> &PreparedNetworkFields {
        self.fields
    }

    pub fn cells(&self) -> &[NetworkCylinderAbsorptionCell] {
        &self.cells
    }

    /// Combine the complete coherent source vector before either loss square.
    pub fn sample(
        &self,
        excitation: &[C64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<NetworkCylinderAbsorption>, DuctError> {
        self.observe_coherence(Some(excitation), &[], cancelled)
    }

    /// Cycle-mean heating for a mean source vector and covariance factors.
    /// Source = mean + sum(factor_r z_r), E[z_r]=0 and
    /// E[z_r conj(z_s)]=delta_rs. Each factor is coherent internally; independent
    /// factors contribute positive powers. No Gaussian or proper-process law
    /// is required for this cycle mean. All entries use the borrowed owner's
    /// canonical source units and clock. No cross-frequency average is inferred.
    pub fn observe_coherence(
        &self,
        mean: Option<&[C64]>,
        factors: &[&[C64]],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<NetworkCylinderAbsorption>, DuctError> {
        observe_channels(
            self.fields,
            &self.cells,
            &self.operators,
            mean,
            factors,
            cancelled,
            |watts, total| NetworkCylinderAbsorption {
                viscous_power_w: watts[0],
                thermal_power_w: watts[1],
                absorbed_power_w: total,
            },
        )
    }
}

// One source contraction and coherent axial integration implementation serves
// both equilibrium and molecular owners. The latter supplies its own original
// positive constitutive coefficients, not a substituted equilibrium loss law.
pub(in crate::network) fn observe_channels<const N: usize, T>(
    fields: &PreparedNetworkFields,
    cells: &[NetworkCylinderAbsorptionCell],
    operators: &[CylinderLossOperator<N>],
    mean: Option<&[C64]>,
    factors: &[&[C64]],
    cancelled: impl FnMut() -> bool,
    output: impl FnMut([f64; N], f64) -> T,
) -> Result<Vec<T>, DuctError> {
    observe_loss_forms(
        fields,
        cells.len(),
        mean,
        factors,
        cancelled,
        |index, values| operators[index].channels(fields, &cells[index], values),
        output,
    )
}

/// Same canonical-source contractions and positive factor sum for native
/// cylinder and spherical-cone forms. The supplied form owns its geometry.
pub(in crate::network) fn observe_loss_forms<const N: usize, T>(
    fields: &PreparedNetworkFields,
    count: usize,
    mean: Option<&[C64]>,
    factors: &[&[C64]],
    mut cancelled: impl FnMut() -> bool,
    mut channels: impl FnMut(usize, &[C64]) -> Result<[Dd; N], DuctError>,
    mut output: impl FnMut([f64; N], f64) -> T,
) -> Result<Vec<T>, DuctError> {
    poll(&mut cancelled, 0)?;
    if factors.len() > MAX_FIELD_FACTORS {
        return Err(bad("cylinder absorption factor cap exceeded"));
    }
    if let Some(values) = mean {
        fields.validate_excitation(values)?;
    }
    for values in factors {
        fields.validate_excitation(values)?;
    }
    let mut result = Vec::with_capacity(count);
    for index in 0..count {
        poll(&mut cancelled, index)?;
        let mut watts = [Dd::ZERO; N];
        if let Some(values) = mean {
            watts = channels(index, values)?;
        }
        for (factor, values) in factors.iter().enumerate() {
            poll(&mut cancelled, factor)?;
            let extra = channels(index, values)?;
            for channel in 0..N {
                watts[channel] = watts[channel] + extra[channel];
            }
        }
        let total = watts.iter().fold(Dd::ZERO, |sum, &v| sum + v).to_f64();
        let watts = watts.map(Dd::to_f64);
        if watts
            .into_iter()
            .chain([total])
            .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(bad("cylinder cell heating is unrepresentable"));
        }
        result.push(output(watts, total));
    }
    poll(&mut cancelled, result.len())?;
    Ok(result)
}

impl<const N: usize> CylinderLossOperator<N> {
    pub(in crate::network) fn channels(
        &self,
        fields: &PreparedNetworkFields,
        cell: &NetworkCylinderAbsorptionCell,
        excitation: &[C64],
    ) -> Result<[Dd; N], DuctError> {
        let op = self;
        if op.weights.into_iter().all(|v| v == 0.0) {
            return Ok([Dd::ZERO; N]);
        }
        let [left, right] = cell.probe_indices;
        let [p_left, u_left] = fields.parts(left, excitation)?;
        let [p_right, u_right] = fields.parts(right, excitation)?;
        let zu_left = multiply(u_left, op.impedance);
        let zu_right = multiply(u_right, op.impedance);
        let half = Dd::from_f64(0.5);
        // Anchor the forward wave at the cell inlet and the backward wave at
        // its outlet. Neither amplitude is reconstructed with growing travel.
        let f = [
            (p_left[0] + zu_left[0]) * half,
            (p_left[1] + zu_left[1]) * half,
        ];
        let b = [
            (p_right[0] - zu_right[0]) * half,
            (p_right[1] - zu_right[1]) * half,
        ];
        let sum = [f[0] + b[0], f[1] + b[1]];
        let difference = [f[0] - b[0], f[1] - b[1]];
        let plus = Dd::from_f64(op.plus_m);
        let minus = Dd::from_f64(op.minus_m);
        let mut result = [Dd::ZERO; N];
        result[0] = if op.weights[0] == 0.0 {
            Dd::ZERO
        } else {
            (minus * norm_squared(divide(sum, op.impedance))
                + plus * norm_squared(divide(difference, op.impedance)))
                * Dd::from_f64(op.weights[0])
                * half
        };
        if op.weights[1..].iter().any(|&v| v != 0.0) {
            let pressure = plus * norm_squared(sum) + minus * norm_squared(difference);
            for channel in 1..N {
                if op.weights[channel] != 0.0 {
                    result[channel] = pressure * Dd::from_f64(op.weights[channel]) * half;
                }
            }
        }
        Ok(result)
    }
}

pub(in crate::network) fn multiply(z: [Dd; 2], w: C64) -> [Dd; 2] {
    let re = Dd::from_f64(w.re);
    let im = Dd::from_f64(w.im);
    [z[0] * re - z[1] * im, z[0] * im + z[1] * re]
}

pub(in crate::network) fn divide(z: [Dd; 2], w: C64) -> [Dd; 2] {
    let scale = w.re.abs().max(w.im.abs());
    let re = Dd::from_f64(w.re / scale);
    let im = Dd::from_f64(w.im / scale);
    let denominator = re * re + im * im;
    let scale = Dd::from_f64(scale);
    [
        (z[0] * re + z[1] * im) / denominator / scale,
        (z[1] * re - z[0] * im) / denominator / scale,
    ]
}

pub(in crate::network) fn norm_squared(z: [Dd; 2]) -> Dd {
    z[0] * z[0] + z[1] * z[1]
}
