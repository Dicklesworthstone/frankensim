//! Reversible electromagnetic and original polarization energy at one film point.
//!
//! Complete coherent fields are summed before squaring. Endpoint changes use
//! analytic field increments, delta U = delta q * (2 q + delta q) / 2, rather
//! than subtracting two stored energies. A proper source removes one common
//! carrier; a real source restores that carrier and its increment explicitly.

use super::*;
use fs_fft::harmonic_window::phase_rad_s_dd_from_clock_dd;
use integrated_heat::ScaledDd;

type Complex = [ScaledDd; 2];
type Vector = [Complex; 3];

/// Preparation bounds for an original-field energy observer. Reconstruction
/// retains the original layer's point-line and FFT costs. Coordinate terms
/// count three axes for E, H, each restoring pole and each Lorentz kinetic pole.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumExposureLayerEnergyLimits {
    pub maximum_point_line_pairs: usize,
    pub maximum_fft_work: u64,
    pub maximum_coordinate_frequency_terms: usize,
    /// Original layer workspace plus retained fields/banks and one temporary
    /// frequency-point bank. Borrowed sources and allocator overhead excluded.
    pub maximum_workspace_bytes: usize,
}

/// Reversible energy density [J/m^3], or its signed endpoint change when returned
/// by `change_from_clock_dd`. The total includes every named subaccount once.
/// Thermal energy, parameter work and switch-on transients are not represented.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AngularSpectrumExposureLayerEnergyDensity {
    pub electromagnetic_j_m3: Dd,
    pub polarization_restoring_j_m3: Dd,
    pub polarization_kinetic_j_m3: Dd,
    pub total_j_m3: Dd,
}

/// Original settled fields at one original grid point and one selected depth.
/// Queries allocate no memory, perform no new field/material solve, and have
/// linear work in the retained coordinate-frequency terms. Each independent
/// source factor keeps its original real/proper statistics. Source uncertainty
/// and spatial resolution remain those of the prepared layer.
pub struct PreparedAngularSpectrumExposureLayerEnergy<'p, 'e, 'b, 's, 'm> {
    layer: &'p PreparedAngularSpectrumExposureLayer<'e, 'b, 's, 'm>,
    depth_m: f64,
    point_index: usize,
    banks: Vec<EnergyBank>,
    lines: Vec<EnergyLine>,
    coordinate_frequency_terms: usize,
    workspace_bytes: usize,
}

impl<'p, 'e, 'b, 's, 'm> PreparedAngularSpectrumExposureLayerEnergy<'p, 'e, 'b, 's, 'm> {
    pub fn layer(&self) -> &'p PreparedAngularSpectrumExposureLayer<'e, 'b, 's, 'm> {
        self.layer
    }
    pub fn depth_m(&self) -> f64 {
        self.depth_m
    }
    pub fn point_index(&self) -> usize {
        self.point_index
    }
    pub fn grid(&self) -> AngularSpectrumGrid {
        self.layer.anchor.grid()
    }
    pub fn coordinate_frequency_terms(&self) -> usize {
        self.coordinate_frequency_terms
    }
    pub fn preparation_workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Expected instantaneous original-field stored energy. The DD source epoch
    /// is not rounded to f64; proper factors retain their physical half-factor.
    pub fn evaluate_from_clock_dd(
        &self,
        time_s: Dd,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureLayerEnergyDensity, AngularSpectrumError> {
        energy(&self.banks, &self.lines, time_s, None, &mut cancelled)
    }

    /// U(start+duration)-U(start), retaining the start and positive duration
    /// independently. Analytic increments precede products, so neither a short
    /// accepted step nor a small storage change is lost in large endpoint values.
    /// This signed reversible change is not deposited heat or source work.
    pub fn change_from_clock_dd(
        &self,
        start_s: Dd,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureLayerEnergyDensity, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(bad("film energy interval needs positive finite duration"));
        }
        energy(
            &self.banks,
            &self.lines,
            start_s,
            Some(duration_s),
            &mut cancelled,
        )
    }
}

struct EnergyLine {
    omega: f64,
    field: Vector,
}

struct EnergyBank {
    start: usize,
    end: usize,
    proper: bool,
    account: usize,
    scale: ScaledDd,
}

struct PointField {
    omega: f64,
    cache: usize,
    field: [C64; 6],
}

#[derive(Clone, Copy)]
enum Coordinate {
    Electric,
    Magnetic,
    Restoring(usize),
    Kinetic(usize),
}

impl<'e, 'b, 's, 'm> PreparedAngularSpectrumExposureLayer<'e, 'b, 's, 'm> {
    /// Retain the original coherent E/H and normalized restoring/kinetic fields
    /// at one point. The material's existing causal reservoir law authorizes
    /// U=epsilon_inf E^2/2+mu H^2/2+sum(P^2/Delta)/2
    ///   +sum(P_dot^2/(omega0^2 Delta))/2.
    /// No complex-permittivity energy shortcut or new material fit is used.
    pub fn point_energy<'p>(
        &'p self,
        depth_m: f64,
        point_index: usize,
        limits: AngularSpectrumExposureLayerEnergyLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAngularSpectrumExposureLayerEnergy<'p, 'e, 'b, 's, 'm>, AngularSpectrumError>
    {
        poll(&mut cancelled)?;
        if !depth_m.is_finite() || depth_m < 0.0 || depth_m > self.thickness_m {
            return Err(bad("film energy depth is outside the selected layer"));
        }
        if point_index >= self.anchor.modes.len() {
            return Err(bad("film energy point is outside the original grid"));
        }
        let coordinates = 2usize
            .checked_add(self.law.pole_count())
            .and_then(|n| n.checked_add(self.law.lorentz_poles.len()))
            .ok_or_else(overflow)?;
        let mut lines_count = 0usize;
        let mut banks_count = 0usize;
        let mut maximum_frequencies = 0usize;
        for (bank, _) in ensemble::banks(self.mean, self.factors, self.statistics) {
            poll(&mut cancelled)?;
            let frequencies = bank.frequency_count();
            maximum_frequencies = maximum_frequencies.max(frequencies);
            lines_count = frequencies
                .checked_mul(coordinates)
                .and_then(|n| n.checked_add(lines_count))
                .ok_or_else(overflow)?;
            if frequencies != 0 {
                banks_count = banks_count.checked_add(coordinates).ok_or_else(overflow)?;
            }
        }
        let coordinate_frequency_terms = lines_count.checked_mul(3).ok_or_else(overflow)?;
        let workspace_bytes = lines_count
            .checked_mul(core::mem::size_of::<EnergyLine>())
            .and_then(|n| {
                n.checked_add(banks_count.checked_mul(core::mem::size_of::<EnergyBank>())?)
            })
            .and_then(|n| {
                n.checked_add(maximum_frequencies.checked_mul(core::mem::size_of::<PointField>())?)
            })
            .and_then(|n| n.checked_add(self.workspace_bytes))
            .ok_or_else(overflow)?;
        if self.point_line_pairs > limits.maximum_point_line_pairs
            || self.field_fft_work > limits.maximum_fft_work
            || coordinate_frequency_terms > limits.maximum_coordinate_frequency_terms
            || workspace_bytes > limits.maximum_workspace_bytes
        {
            return Err(bad("film energy exceeds its work or workspace budget"));
        }
        let mut points = allocated(maximum_frequencies)?;
        let mut lines = allocated(lines_count)?;
        let mut banks = allocated(banks_count)?;
        for (bank, proper) in ensemble::banks(self.mean, self.factors, self.statistics) {
            poll(&mut cancelled)?;
            points.clear();
            for group in 0..bank.frequency_count() {
                poll(&mut cancelled)?;
                let omega = bank.frequencies[group];
                let cache = self
                    .frequencies
                    .binary_search_by(|value| value.omega.total_cmp(&omega))
                    .map_err(|_| bad("original film frequency cache is missing"))?;
                let layer = &self.frequencies[cache].layer;
                let fields =
                    bank.frequency_field_with(group, &mut cancelled, |line, cancelled| {
                        line.modal_layer_field(layer.as_ref(), self.layer_index, depth_m, cancelled)
                    })?;
                points.push(PointField {
                    omega,
                    cache,
                    field: fields[point_index],
                });
            }
            if points.is_empty() {
                continue;
            }
            let coordinates = core::iter::once(Coordinate::Electric)
                .chain(core::iter::once(Coordinate::Magnetic))
                .chain((0..self.law.pole_count()).map(Coordinate::Restoring))
                .chain(
                    (self.law.debye_poles.len()..self.law.pole_count()).map(Coordinate::Kinetic),
                );
            for coordinate in coordinates {
                poll(&mut cancelled)?;
                let start = lines.len();
                for point in &points {
                    poll(&mut cancelled)?;
                    let (offset, transfer) = match coordinate {
                        Coordinate::Electric => (0, C64::ONE),
                        Coordinate::Magnetic => (3, C64::ONE),
                        Coordinate::Restoring(pole) => {
                            (0, self.frequencies[point.cache].poles[pole].storage)
                        }
                        Coordinate::Kinetic(pole) => {
                            (0, self.frequencies[point.cache].poles[pole].kinetic_storage)
                        }
                    };
                    lines.push(EnergyLine {
                        omega: point.omega,
                        field: core::array::from_fn(|c| {
                            filtered(point.field[c + offset], transfer)
                        }),
                    });
                }
                let (account, coefficient) = match coordinate {
                    Coordinate::Electric => (0, self.law.epsilon_infinity),
                    Coordinate::Magnetic => (0, self.law.permeability),
                    Coordinate::Restoring(_) => (1, 1.0),
                    Coordinate::Kinetic(_) => (2, 1.0),
                };
                banks.push(EnergyBank {
                    start,
                    end: lines.len(),
                    proper,
                    account,
                    scale: ScaledDd::from_f64(coefficient) * ScaledDd::from_f64(0.5),
                });
            }
        }
        poll(&mut cancelled)?;
        Ok(PreparedAngularSpectrumExposureLayerEnergy {
            layer: self,
            depth_m,
            point_index,
            banks,
            lines,
            coordinate_frequency_terms,
            workspace_bytes,
        })
    }
}

fn overflow() -> AngularSpectrumError {
    bad("film energy work or storage count overflows")
}

fn allocated<T>(count: usize) -> Result<Vec<T>, AngularSpectrumError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| bad("film energy allocation failed"))?;
    Ok(values)
}

// Conjugate native +i*omega*t amplitudes into the shared negative-time law.
// Exponents stay separate until field, material and physical half-scales meet.
fn filtered(value: C64, transfer: C64) -> Complex {
    let [a, b, c, d] = [value.re, value.im, transfer.re, transfer.im].map(ScaledDd::from_f64);
    [a * c - b * d, -(a * d + b * c)]
}

fn rotate(value: Complex, phase: [Dd; 2]) -> Complex {
    multiply(value, phase.map(ScaledDd::from_dd))
}

fn multiply(value: Complex, [c, s]: Complex) -> Complex {
    [value[0] * c - value[1] * s, value[0] * s + value[1] * c]
}

fn phase(omega: Dd, time: Dd) -> Result<[Dd; 2], AngularSpectrumError> {
    phase_rad_s_dd_from_clock_dd(omega, time)
        .map_err(|_| bad("film energy source phase is unrepresentable"))
}

// exp(-i omega duration)-1, avoiding cos(theta)-1 near a zero increment.
fn increment(omega: Dd, duration: f64) -> Result<Complex, AngularSpectrumError> {
    let half_omega = omega * Dd::from_f64(0.5);
    if omega != Dd::ZERO && half_omega * Dd::from_f64(duration) == Dd::ZERO {
        return Err(bad("film energy phase increment is unrepresentable"));
    }
    let half = phase(half_omega, Dd::from_f64(duration))?.map(ScaledDd::from_dd);
    let two = ScaledDd::from_f64(2.0);
    Ok([-two * half[1] * half[1], two * half[0] * half[1]])
}

fn energy(
    banks: &[EnergyBank],
    lines: &[EnergyLine],
    time: Dd,
    duration: Option<f64>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AngularSpectrumExposureLayerEnergyDensity, AngularSpectrumError> {
    poll(cancelled)?;
    if !time.hi.is_finite() || !time.lo.is_finite() {
        return Err(bad("film energy source clock is nonfinite"));
    }
    let time = Dd::from_pair(time.hi, time.lo);
    if !time.hi.is_finite() || !time.lo.is_finite() {
        return Err(bad("film energy source clock is unrepresentable"));
    }
    let mut accounts = [ScaledDd::ZERO; 3];
    let two = ScaledDd::from_f64(2.0);
    for bank in banks {
        poll(cancelled)?;
        let fields = &lines[bank.start..bank.end];
        let Some(first) = fields.first() else {
            continue;
        };
        let reference = Dd::from_f64(first.omega);
        let mut value = [[ScaledDd::ZERO; 2]; 3];
        let mut delta = value;
        for line in fields {
            poll(cancelled)?;
            let omega = Dd::from_f64(line.omega) - reference;
            let rotation = phase(omega, time)?;
            let difference = duration.map(|dt| increment(omega, dt)).transpose()?;
            for c in 0..3 {
                let q = rotate(line.field[c], rotation);
                for i in 0..2 {
                    value[c][i] = value[c][i] + q[i];
                }
                if let Some(difference) = difference {
                    let dq = multiply(q, difference);
                    for i in 0..2 {
                        delta[c][i] = delta[c][i] + dq[i];
                    }
                }
            }
        }
        let scale = bank.scale * ScaledDd::from_f64(if bank.proper { 0.5 } else { 1.0 });
        if bank.proper {
            for c in 0..3 {
                for i in 0..2 {
                    let product = if duration.is_some() {
                        delta[c][i] * (two * value[c][i] + delta[c][i])
                    } else {
                        value[c][i] * value[c][i]
                    };
                    accounts[bank.account] = accounts[bank.account] + scale * product;
                }
            }
        } else {
            let carrier = phase(reference, time)?;
            let carrier_delta = duration.map(|dt| increment(reference, dt)).transpose()?;
            for c in 0..3 {
                let q = rotate(value[c], carrier)[0];
                let product = if let Some(change) = carrier_delta {
                    // delta[(F+delta F)*carrier] = carrier *
                    // [delta F + (F+delta F)*(delta carrier / carrier)].
                    let end = [value[c][0] + delta[c][0], value[c][1] + delta[c][1]];
                    let shifted = multiply(end, change);
                    let dq = rotate(
                        [delta[c][0] + shifted[0], delta[c][1] + shifted[1]],
                        carrier,
                    )[0];
                    dq * (two * q + dq)
                } else {
                    q * q
                };
                accounts[bank.account] = accounts[bank.account] + scale * product;
            }
        }
    }
    let physical = |value: ScaledDd| {
        value
            .physical()
            .map_err(|_| bad("film stored energy or its change is unrepresentable"))
    };
    let electromagnetic_j_m3 = physical(accounts[0])?;
    let polarization_restoring_j_m3 = physical(accounts[1])?;
    let polarization_kinetic_j_m3 = physical(accounts[2])?;
    let total_j_m3 = physical(accounts[0] + accounts[1] + accounts[2])?;
    poll(cancelled)?;
    Ok(AngularSpectrumExposureLayerEnergyDensity {
        electromagnetic_j_m3,
        polarization_restoring_j_m3,
        polarization_kinetic_j_m3,
        total_j_m3,
    })
}

#[cfg(test)]
mod tests;
