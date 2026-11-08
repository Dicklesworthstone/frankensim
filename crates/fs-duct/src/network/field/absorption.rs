//! Positive cycle-mean viscous and thermal work in original cylinders.
//!
//! In exp(-i omega t), p'=i k Zc U and U'=i k p/Zc. Therefore
//! -d Re(p U*)/dx / 2 = [Im(k Zc)|U|^2 + Im(k/Zc)|p|^2]/2.
//! The two weights are read directly from the original Bessel series/shunt
//! law before its square roots, preserving small losses and their mechanisms.
//! Exact positive wave integrals retain coherent pressure and velocity nodes;
//! finite-window boundary work is never used as a heat source.

use super::*;
use crate::compose::pressure_wave_integral_weights;
use crate::WallDiffusionInputs;

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
struct CellOperator {
    impedance: C64,
    plus_m: f64,
    minus_m: f64,
    viscous_weight: f64,
    thermal_weight: f64,
}

impl CellOperator {
    const ZERO: Self = Self {
        impedance: C64::ONE,
        plus_m: 0.0,
        minus_m: 0.0,
        viscous_weight: 0.0,
        thermal_weight: 0.0,
    };
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
        if probe_pairs.is_empty()
            || probe_pairs.len() > maximum_cells
            || maximum_cells > MAX_FIELD_PROBES
        {
            return Err(bad("cylinder absorption needs a bounded nonempty cell list"));
        }
        if self.loss == LossModel::WideTube {
            return Err(bad("cylinder heat requires Bessel or exact lossless propagation"));
        }
        let mut cells = Vec::with_capacity(probe_pairs.len());
        let mut operators = Vec::with_capacity(probe_pairs.len());
        // At most one original constitutive evaluation per selected cylinder;
        // this temporary cache is also bounded by the number of cells.
        let mut constitutive = std::collections::BTreeMap::new();
        for (index, &probe_indices) in probe_pairs.iter().enumerate() {
            poll(&mut cancelled, index)?;
            let [left, right] = probe_indices;
            let a = self.cylinders.get(left).and_then(Option::as_ref)
                .ok_or_else(|| bad("cell start is not an original cylinder probe"))?;
            let b = self.cylinders.get(right).and_then(Option::as_ref)
                .ok_or_else(|| bad("cell end is not an original cylinder probe"))?;
            if self.locations[left].probe.edge_id != self.locations[right].probe.edge_id
                || a.segment != b.segment
                || a.offset_m > b.offset_m
            {
                return Err(bad("cell endpoints must be ordered in the same original cylinder"));
            }
            let length = b.offset_m - a.offset_m;
            let volume_m3 = self.locations[left].area_m2 * length;
            if !volume_m3.is_finite() || (length > 0.0 && volume_m3 == 0.0) {
                return Err(bad("cylinder cell volume is unrepresentable"));
            }
            let operator = if self.loss == LossModel::Lossless || length == 0.0 {
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
                            inputs, a.radius_m, self.ports.omega(),
                        )?;
                        let (k, impedance) = crate::helmholtz_wave_from_telegraph(series, shunt);
                        if !finite(impedance) || impedance == C64::ZERO {
                            return Err(DuctError::Singular);
                        }
                        let cached = (k, impedance, series.re, shunt.re);
                        constitutive.insert(key, cached);
                        cached
                    };
                let (plus_m, minus_m) = pressure_wave_integral_weights(k, length)?;
                CellOperator { impedance, plus_m, minus_m, viscous_weight, thermal_weight }
            };
            cells.push(NetworkCylinderAbsorptionCell {
                probe_indices,
                segment: a.segment,
                start_m: a.offset_m,
                end_m: b.offset_m,
                volume_m3,
            });
            operators.push(operator);
        }
        poll(&mut cancelled, cells.len())?;
        Ok(PreparedNetworkCylinderAbsorption { fields: self, cells, operators })
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
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<NetworkCylinderAbsorption>, DuctError> {
        poll(&mut cancelled, 0)?;
        if factors.len() > MAX_FIELD_FACTORS {
            return Err(bad("cylinder absorption factor cap exceeded"));
        }
        if let Some(values) = mean {
            self.fields.validate_excitation(values)?;
        }
        for values in factors {
            self.fields.validate_excitation(values)?;
        }
        let mut result = Vec::with_capacity(self.cells.len());
        for (index, (cell, operator)) in self.cells.iter().zip(&self.operators).enumerate() {
            poll(&mut cancelled, index)?;
            let mut watts = [Dd::ZERO; 2];
            if let Some(values) = mean {
                watts = self.channels(cell, operator, values)?;
            }
            for (factor, values) in factors.iter().enumerate() {
                poll(&mut cancelled, factor)?;
                let extra = self.channels(cell, operator, values)?;
                watts[0] = watts[0] + extra[0];
                watts[1] = watts[1] + extra[1];
            }
            let output = NetworkCylinderAbsorption {
                viscous_power_w: watts[0].to_f64(),
                thermal_power_w: watts[1].to_f64(),
                absorbed_power_w: (watts[0] + watts[1]).to_f64(),
            };
            if [output.viscous_power_w, output.thermal_power_w, output.absorbed_power_w]
                .into_iter().any(|v| !v.is_finite() || v < 0.0)
            {
                return Err(bad("cylinder cell heating is unrepresentable"));
            }
            result.push(output);
        }
        poll(&mut cancelled, result.len())?;
        Ok(result)
    }

    fn channels(
        &self,
        cell: &NetworkCylinderAbsorptionCell,
        operator: &CellOperator,
        excitation: &[C64],
    ) -> Result<[Dd; 2], DuctError> {
        let op = operator;
        if op.viscous_weight == 0.0 && op.thermal_weight == 0.0 {
            return Ok([Dd::ZERO; 2]);
        }
        let [left, right] = cell.probe_indices;
        let [p_left, u_left] = self.fields.parts(left, excitation)?;
        let [p_right, u_right] = self.fields.parts(right, excitation)?;
        let zu_left = multiply(u_left, op.impedance);
        let zu_right = multiply(u_right, op.impedance);
        let half = Dd::from_f64(0.5);
        // Anchor the forward wave at the cell inlet and the backward wave at
        // its outlet. Neither amplitude is reconstructed with growing travel.
        let f = [(p_left[0] + zu_left[0]) * half, (p_left[1] + zu_left[1]) * half];
        let b = [(p_right[0] - zu_right[0]) * half, (p_right[1] - zu_right[1]) * half];
        let sum = [f[0] + b[0], f[1] + b[1]];
        let difference = [f[0] - b[0], f[1] - b[1]];
        let plus = Dd::from_f64(op.plus_m);
        let minus = Dd::from_f64(op.minus_m);
        let viscous = if op.viscous_weight == 0.0 {
            Dd::ZERO
        } else {
            (minus * norm_squared(divide(sum, op.impedance))
                + plus * norm_squared(divide(difference, op.impedance)))
                * Dd::from_f64(op.viscous_weight) * half
        };
        let thermal = if op.thermal_weight == 0.0 {
            Dd::ZERO
        } else {
            (plus * norm_squared(sum) + minus * norm_squared(difference))
                * Dd::from_f64(op.thermal_weight) * half
        };
        Ok([viscous, thermal])
    }
}

fn multiply(z: [Dd; 2], w: C64) -> [Dd; 2] {
    let re = Dd::from_f64(w.re);
    let im = Dd::from_f64(w.im);
    [z[0] * re - z[1] * im, z[0] * im + z[1] * re]
}

fn divide(z: [Dd; 2], w: C64) -> [Dd; 2] {
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

fn norm_squared(z: [Dd; 2]) -> Dd {
    z[0] * z[0] + z[1] * z[1]
}
