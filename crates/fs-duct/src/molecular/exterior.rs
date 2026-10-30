//! Reciprocal exterior feedback at every molecular chimney and the terminal.
//!
//! The original decaying wave blocks are reduced to an exposed-port scattering
//! operator. Only the aperture count enters the dense closure; axial resolution
//! and wall diffusion poles do not enlarge it. Bore spans are composed before
//! adjoining the next aperture. An exterior owner supplies the complete acoustic
//! matrix, including all self/mutual terms, in physical pressure/volume-flow
//! coordinates. No radiation law, exterior gas or end correction is inferred.
//!
//! The compact junction and vent mixture retain their original model scope.
//! A vent's physical open-mouth flow is sigma times the full-radius branch flow.
//! Consequently its chimney impedance scales by 1/sigma. To reproduce an old
//! independent scalar mouth load Z, put Z/sigma on this matrix's diagonal. This
//! convention preserves reciprocal work with mutual entries: p = Z_rad U,
//! where U is the actual weighted outward flow, without weighting power twice.

use super::*;
use crate::network::{CoupledTermination, MAX_NETWORK_PORTS};
use fs_la::eigen_complex::lu_complex;

/// Frequency-prepared molecular duct with every open/vented mouth and terminal
/// exposed. The inlet is scattering port 0; `exterior_labels` follows it.
/// Scattering uses one positive real reference impedance for all ports.
#[derive(Debug, Clone)]
pub struct MolecularPortScattering {
    omega: f64,
    reference: f64,
    labels: Vec<String>,
    scattering: Vec<C64>,
    min_shear: f64,
    composition_condition: f64,
    maximum_condition: f64,
}

/// The closed radiation/duct solution for unit volume flow at the inlet.
/// Exterior arrays follow `MolecularPortScattering::exterior_labels`, even if
/// the supplied matrix uses a different label order. Individual mouth work may
/// be negative under mutual loading; only the summed exterior work is passive.
#[derive(Debug, Clone)]
pub struct MolecularExteriorResponse {
    /// Inlet pressure per inlet volume flow [Pa s/m3].
    pub z_in: C64,
    /// Mean mouth pressures per inlet volume flow [Pa s/m3].
    pub exterior_pressures_over_u_in: Vec<C64>,
    /// Physical flows directed OUT of the duct, per inlet flow.
    pub exterior_flows_over_u_in: Vec<C64>,
    /// Total mean exterior work / squared inlet flow [Pa s/m3].
    pub exterior_power_per_squared_inlet_flow: f64,
    /// Input work minus exterior work, without clipping [Pa s/m3].
    /// Includes molecular and wall heat in bore, open chimneys and sealed stubs.
    pub dissipated_power_per_squared_inlet_flow: f64,
    /// Infinity-norm condition of the dimensionless row-equilibrated closure.
    pub condition_estimate: f64,
}

impl MolecularDuct {
    /// Prepare all open/vented mouth ports in authored order, followed by the
    /// terminal/bell port. `labels` must name every one uniquely. Closed tone
    /// holes retain their physical sealed stub but have no exterior port.
    /// The hard exterior cap is `network::MAX_NETWORK_PORTS` (including bell).
    pub fn prepare_exterior_ports(
        &self,
        omega: f64,
        labels: &[String],
        options: ComposeOptions,
    ) -> Result<MolecularPortScattering, MolecularDuctError> {
        self.prepare_exterior_ports_controlled(omega, labels, options, &mut || false)
    }

    pub fn prepare_exterior_ports_controlled(
        &self,
        omega: f64,
        labels: &[String],
        options: ComposeOptions,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularPortScattering, MolecularDuctError> {
        poll(cancelled, 0)?;
        options.validate()?;
        let mouths = self
            .pieces
            .iter()
            .filter(|piece| matches!(piece, Piece::ToneHole(hole) if hole.site.open_fraction > 0.0))
            .count();
        if labels.len() != mouths + 1
            || labels.len() > MAX_NETWORK_PORTS
            || labels
                .iter()
                .enumerate()
                .any(|(i, label)| label.trim().is_empty() || labels[..i].contains(label))
        {
            return Err(bad("molecular exterior needs unique labels for each open mouth and terminal within the port cap").into());
        }
        if !positive(omega) || self.blocks() > options.maximum_blocks {
            return Err(bad(
                "molecular exterior needs finite positive frequency within the block cap",
            )
            .into());
        }
        let (waves, min_shear) = self.waves_controlled(omega, cancelled)?;
        let reference = self.reference_impedance;
        let mut ports = Ports {
            n: 2,
            s: TwoPort::through(reference)?.scattering.to_vec(),
            condition: 1.0,
        };
        let mut pending = TwoPort::through(reference)?;
        for (i, piece) in self.pieces.iter().enumerate() {
            poll(cancelled, i)?;
            if let Piece::ToneHole(hole) = piece {
                let (series, junction, branch) = hole.exposed_parts(
                    &waves[hole.point],
                    self.state_at_point(hole.point).gas().density,
                    omega,
                    reference,
                    options,
                )?;
                let (next, condition) = pending.star(series, options.maximum_condition)?;
                pending = next;
                ports.condition = ports.condition.max(condition);
                if let Some(branch) = branch {
                    ports.connect(ports.n - 1, &pending.scattering, 2, options, cancelled)?;
                    ports.connect(ports.n - 1, &junction, 3, options, cancelled)?;
                    ports.connect(ports.n - 2, &branch.scattering, 2, options, cancelled)?;
                    pending = series;
                } else {
                    let shunt = TwoPort::new(
                        junction.try_into().expect("closed junction has two ports"),
                        reference,
                    )?;
                    let (next, condition) = pending.star(shunt, options.maximum_condition)?;
                    pending = next;
                    ports.condition = ports.condition.max(condition);
                    let (next, condition) = pending.star(series, options.maximum_condition)?;
                    pending = next;
                    ports.condition = ports.condition.max(condition);
                }
            } else {
                let block = self.piece_wave_block(piece, &waves)?;
                let (next, condition) = pending.star(block, options.maximum_condition)?;
                pending = next;
                ports.condition = ports.condition.max(condition);
            }
        }
        ports.connect(ports.n - 1, &pending.scattering, 2, options, cancelled)?;
        Ok(MolecularPortScattering {
            omega,
            reference,
            labels: labels.to_vec(),
            scattering: ports.s,
            min_shear,
            composition_condition: ports.condition,
            maximum_condition: options.maximum_condition,
        })
    }
}

impl MolecularPortScattering {
    pub fn omega(&self) -> f64 {
        self.omega
    }
    pub fn reference_impedance(&self) -> f64 {
        self.reference
    }
    /// Inlet first, then exterior labels; row-major square scattering matrix.
    pub fn scattering(&self) -> &[C64] {
        &self.scattering
    }
    pub fn exterior_labels(&self) -> &[String] {
        &self.labels
    }
    pub fn min_shear_number(&self) -> f64 {
        self.min_shear
    }
    pub fn composition_condition(&self) -> f64 {
        self.composition_condition
    }

    /// Close every exterior port by one passive acoustic matrix, p = Z U.
    /// A singular/rank-deficient Z, including pressure release, is admissible:
    /// constitutive equations are retained and Z is never inverted. The
    /// existing network admission checks its actual Hermitian part. Matrices
    /// from BEM can be passed after the owner's acoustic port conversion.
    pub fn response(
        &self,
        load: &CoupledTermination,
    ) -> Result<MolecularExteriorResponse, MolecularDuctError> {
        self.response_controlled(load, &mut || false)
    }

    pub fn response_controlled(
        &self,
        load: &CoupledTermination,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularExteriorResponse, MolecularDuctError> {
        poll(cancelled, 0)?;
        crate::network::admit_coupled(load, self.omega)?;
        if load.labels.len() != self.labels.len() {
            return Err(bad(
                "molecular exterior matrix must cover every exposed mouth and terminal",
            )
            .into());
        }
        let mut order = Vec::with_capacity(self.labels.len());
        for label in &self.labels {
            let index = load
                .labels
                .iter()
                .position(|entry| entry == label)
                .ok_or_else(|| {
                    bad("molecular exterior matrix label does not match an exposed port")
                })?;
            if order.contains(&index) {
                return Err(bad("molecular exterior matrix labels must be unique").into());
            }
            order.push(index);
        }
        let n = self.labels.len() + 1;
        let m = n - 1;
        let s = &self.scattering;
        let mut a = vec![C64::ZERO; n * n];
        // Incoming/outgoing waves are divided by sqrt(R): p=R(a+b), U=a-b.
        // Unit inlet flow, then R(a+b)=Z(b-a) at each exterior port.
        for j in 0..n {
            a[j] = (if j == 0 { C64::ONE } else { C64::ZERO }) - s[j];
        }
        for i in 0..m {
            poll(cancelled, i)?;
            let scale = load.z[order[i] * m..(order[i] + 1) * m]
                .iter()
                .map(|z| z.abs())
                .fold(self.reference, f64::max);
            for j in 0..n {
                let eye = if j == i + 1 { C64::ONE } else { C64::ZERO };
                let mut value = (eye + s[(i + 1) * n + j]).scale(self.reference / scale);
                for k in 0..m {
                    let delta = if j == k + 1 { C64::ONE } else { C64::ZERO };
                    value = value
                        + load.z[order[i] * m + order[k]].scale(1.0 / scale)
                            * (delta - s[(k + 1) * n + j]);
                }
                a[(i + 1) * n + j] = value;
            }
        }
        let (incoming, condition) = solve_source(a, n, self.maximum_condition, cancelled)?;
        let mut outgoing = vec![C64::ZERO; n];
        for i in 0..n {
            poll(cancelled, i)?;
            for j in 0..n {
                outgoing[i] = outgoing[i] + s[i * n + j] * incoming[j];
            }
        }
        let z_in = (incoming[0] + outgoing[0]).scale(self.reference);
        let mut pressures = Vec::with_capacity(m);
        let mut flows = Vec::with_capacity(m);
        let mut exterior_work = 0.0;
        for i in 1..n {
            let pressure = (incoming[i] + outgoing[i]).scale(self.reference);
            let flow = outgoing[i] - incoming[i];
            exterior_work += 0.5 * (pressure * flow.conj()).re;
            pressures.push(pressure);
            flows.push(flow);
        }
        if !finite(z_in)
            || !exterior_work.is_finite()
            || pressures.iter().chain(&flows).any(|&z| !finite(z))
        {
            return Err(DuctError::Singular.into());
        }
        Ok(MolecularExteriorResponse {
            z_in,
            exterior_pressures_over_u_in: pressures,
            exterior_flows_over_u_in: flows,
            exterior_power_per_squared_inlet_flow: exterior_work,
            dissipated_power_per_squared_inlet_flow: 0.5 * z_in.re - exterior_work,
            condition_estimate: condition,
        })
    }
}

// Join a resolved 2-port to any exposed port, or a 3-port junction to the
// downstream port. Only a single scalar feedback loop is eliminated. Keeping
// every mouth exposed preserves cross-loading, reciprocity and dark channels.
struct Ports {
    n: usize,
    s: Vec<C64>,
    condition: f64,
}
impl Ports {
    fn connect(
        &mut self,
        port: usize,
        right: &[C64],
        rn: usize,
        options: ComposeOptions,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DuctError> {
        let (denominator, condition) = compose::feedback(self.s[port * self.n + port], right[0]);
        if !condition.is_finite() || condition > options.maximum_condition {
            return Err(DuctError::Singular);
        }
        self.condition = self.condition.max(condition);
        let map: Vec<(bool, usize)> = if rn == 2 {
            (0..self.n)
                .map(|i| if i == port { (false, 1) } else { (true, i) })
                .collect()
        } else {
            debug_assert_eq!(port, self.n - 1);
            (0..port)
                .map(|i| (true, i))
                .chain([(false, 1), (false, 2)])
                .collect()
        };
        let n = map.len();
        let mut result = vec![C64::ZERO; n * n];
        for (i, &(left_i, row)) in map.iter().enumerate() {
            poll(cancelled, i)?;
            for (j, &(left_j, col)) in map.iter().enumerate() {
                result[i * n + j] = match (left_i, left_j) {
                    (true, true) => {
                        self.s[row * self.n + col]
                            + self.s[row * self.n + port]
                                * right[0]
                                * (self.s[port * self.n + col] / denominator)
                    }
                    (true, false) => self.s[row * self.n + port] * (right[col] / denominator),
                    (false, true) => right[row * rn] * (self.s[port * self.n + col] / denominator),
                    (false, false) => {
                        right[row * rn + col]
                            + right[row * rn]
                                * self.s[port * self.n + port]
                                * (right[col] / denominator)
                    }
                };
            }
        }
        if result.iter().any(|&z| !finite(z)) {
            return Err(DuctError::Singular);
        }
        self.n = n;
        self.s = result;
        Ok(())
    }
}

// One LU, one source solve; inverse-column row sums give the actual normalized
// infinity-norm condition without storing another dense matrix. N <= 33.
fn solve_source(
    mut a: Vec<C64>,
    n: usize,
    limit: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<C64>, f64), DuctError> {
    let mut source = 0.0;
    for (i, row) in a.chunks_exact_mut(n).enumerate() {
        let scale = row.iter().map(|z| z.abs()).fold(0.0, f64::max);
        if !positive(scale) || row.iter().any(|&z| !finite(z)) {
            return Err(DuctError::Singular);
        }
        for z in row {
            *z = z.scale(1.0 / scale);
        }
        if i == 0 {
            source = 1.0 / scale;
        }
    }
    poll(cancelled, 0)?;
    let lu = lu_complex(&a, n).map_err(|_| DuctError::Singular)?;
    let norm = a
        .chunks_exact(n)
        .map(|row| row.iter().map(|z| z.abs()).sum::<f64>())
        .fold(0.0, f64::max);
    let mut row_sums = vec![0.0; n];
    let mut solution = Vec::new();
    let mut column = vec![C64::ZERO; n];
    for j in 0..n {
        poll(cancelled, j)?;
        column.fill(C64::ZERO);
        column[j] = C64::ONE;
        lu.solve(&mut column);
        for i in 0..n {
            row_sums[i] += column[i].abs();
        }
        if j == 0 {
            solution = column.iter().map(|z| z.scale(source)).collect();
        }
    }
    let condition = norm * row_sums.into_iter().fold(0.0, f64::max);
    if !condition.is_finite() || condition > limit || solution.iter().any(|&z| !finite(z)) {
        return Err(DuctError::Singular);
    }
    Ok((solution, condition))
}

#[cfg(test)]
#[path = "exterior_tests.rs"]
mod tests;
