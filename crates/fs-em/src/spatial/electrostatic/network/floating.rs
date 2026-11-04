//! Charge-constrained floating electrodes on the original tetrahedral field.
//!
//! A floating electrode keeps its TOTAL free charge while its equipotential
//! changes. It can carry nonzero charges on individual faces even when neutral.
//! Partitioning Q=C V gives V_f=C_ff^-1(Q_f-C_fr V_r). The Schur matrix describes
//! differential retained charge at fixed Q_f; it does not discard trapped charge
//! or the energy Q_f.C_ff^-1.Q_f/2. The complete state satisfies
//! Q_r=C_reduced V_r+C_rf C_ff^-1 Q_f and
//! U=V_r.C_reduced.V_r/2+Q_f.C_ff^-1.Q_f/2.
//! Native field recovery checks every state.

use super::{
    ElectrostaticControls, Result, SpatialCapacitanceMatrix, SpatialCapacitanceResponse,
    SpatialElectricField, add, check, finite, invalid, positive,
};
use fs_la::factor::{Cholesky, cholesky};

/// Explicit bounds for the dense reduction, beyond native field extraction.
#[derive(Clone, Copy, Debug)]
pub struct SpatialFloatingCapacitanceControls {
    pub maximum_floating_terminals: usize,
    /// Conservative scalar-product work: f^3 + 2*r*f^2 + r^2*f + n^2.
    /// The native factorization and triangular solves are polled immediately
    /// before and after; their sizes are admitted before allocation.
    pub maximum_dense_work: usize,
    /// Absolute [F] allowance for charge reduction, reciprocity and gauge checks.
    pub coefficient_tolerance_f: f64,
}

/// Reusable floating-terminal factor and Schur reduction of the frozen store.
/// Floating order fixes the charge-input order; retained terminals remain in
/// original order. Subsequent states need no new finite-element field solve.
pub struct SpatialFloatingCapacitance<'s, 'a> {
    source: &'s SpatialCapacitanceMatrix<'a>,
    retained: Vec<usize>,
    floating: Vec<usize>,
    references: Vec<usize>,
    factor: Cholesky,
    scales: Vec<f64>,
    transfer: Vec<f64>,
    reduced_f: Vec<f64>,
    maximum_coefficient_residual_f: f64,
}

/// One mixed voltage/total-charge state, with the original native field.
/// Physical arithmetic keeps the component voltage origin separate so a small
/// charge-induced voltage survives a large, independently chosen common bias.
pub struct SpatialFloatingCapacitanceResponse<'a> {
    /// Complete original terminal order, relative to each component's first
    /// retained terminal. These are the voltages supplied to native field recovery.
    pub component_relative_voltage_v: Vec<f64>,
    /// Physical voltage origin for each original terminal. Adding this to the
    /// relative observation can round away a small signal; fields never do so.
    pub component_reference_voltage_v: Vec<f64>,
    /// Independently recovered weak free-charge reactions, original terminal order.
    pub terminal_charge_c: Vec<f64>,
    /// Affine charge offset C_rf C_ff^-1 Q_f, in retained-terminal order.
    /// Total retained charge is C_reduced V_r plus this offset.
    pub retained_charge_offset_c: Vec<f64>,
    /// Q_f.C_ff^-1.Q_f/2 [J], a part of the native field's stored energy.
    /// This is not an additional store to add to `field().energy_j()`.
    pub trapped_charge_energy_j: f64,
    pub maximum_floating_charge_residual_c: f64,
    pub maximum_reduced_charge_residual_c: f64,
    /// Native field energy minus one half of all original terminal V.Q work.
    pub terminal_energy_residual_j: f64,
    /// Native field energy minus reduced driven energy and trapped-charge energy.
    pub reduced_energy_residual_j: f64,
    field: SpatialElectricField<'a>,
}

impl<'a> SpatialFloatingCapacitanceResponse<'a> {
    /// Original material/geometry field. Nodal voltage observations use the
    /// component-relative origins above, preserving charge-induced small fields.
    pub fn field(&self) -> &SpatialElectricField<'a> {
        &self.field
    }
    pub fn into_field(self) -> SpatialElectricField<'a> {
        self.field
    }
}

impl<'a> SpatialCapacitanceMatrix<'a> {
    /// Prepare a nonempty proper subset of floating electrodes. Every electrical
    /// component carrying terminals must retain a voltage reference, including
    /// components joined through a shared equipotential electrode. Physical
    /// topology admits that reference BEFORE numerical factorization; a gauge
    /// nullspace cannot become a spurious charge store through roundoff.
    ///
    /// The original scalar/tensor field owns the matrix. Positive mutual entries
    /// do not require a pair-capacitor approximation or prevent this reduction.
    pub fn prepare_floating(
        &self,
        floating_terminals: &[usize],
        controls: SpatialFloatingCapacitanceControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialFloatingCapacitance<'_, 'a>> {
        check(cancelled)?;
        let n = self.faces.len();
        let f = floating_terminals.len();
        if f == 0
            || f >= n
            || f > controls.maximum_floating_terminals
            || !positive(controls.coefficient_tolerance_f)
            || floating_terminals
                .iter()
                .enumerate()
                .any(|(i, &t)| t >= n || floating_terminals[..i].contains(&t))
        {
            return Err(invalid(
                "bounded unique floating terminals and retained voltage references required",
            ));
        }
        let r = n - f;
        let work = f.checked_mul(f).and_then(|ff| {
            f.checked_mul(ff)?
                .checked_add(2usize.checked_mul(r)?.checked_mul(ff)?)?
                .checked_add(r.checked_mul(r)?.checked_mul(f)?)?
                .checked_add(n.checked_mul(n)?)
        });
        if work.is_none_or(|work| work > controls.maximum_dense_work) {
            return Err(invalid(
                "floating capacitance factorization exceeds dense work allowance",
            ));
        }
        let retained: Vec<_> = (0..n).filter(|t| !floating_terminals.contains(t)).collect();
        let mut component_reference = vec![None; n];
        for &terminal in &retained {
            check(cancelled)?;
            component_reference[self.references[terminal]].get_or_insert(terminal);
        }
        let references =
            self.references
                .iter()
                .map(|&component| {
                    check(cancelled)?;
                    component_reference[component].ok_or_else(|| invalid(
                "every floating electrical component needs a retained voltage reference",
            ))
                })
                .collect::<Result<Vec<_>>>()?;
        let scales = floating_terminals
            .iter()
            .map(|&terminal| {
                check(cancelled)?;
                let value = self.matrix_f[terminal * n + terminal];
                if positive(value) {
                    Ok(value.sqrt())
                } else {
                    Err(invalid(
                        "floating electrode has no representably positive electric storage",
                    ))
                }
            })
            .collect::<Result<Vec<_>>>()?;
        // Dimensionless diagonal equilibration admits picofarad devices using
        // the same native factorization as large capacitors, without a farad cutoff.
        let mut block = vec![0.0; f * f];
        for (i, &a) in floating_terminals.iter().enumerate() {
            check(cancelled)?;
            for (j, &b) in floating_terminals[..=i].iter().enumerate() {
                let value = finite((self.matrix_f[a * n + b] / scales[i]) / scales[j])?;
                block[i * f + j] = value;
                block[j * f + i] = value;
            }
        }
        check(cancelled)?;
        let factor = cholesky(&block, f)
            .map_err(|_| invalid("floating capacitance block is not positive definite"))?;
        check(cancelled)?;
        let mut transfer = vec![0.0; f * r];
        for (column, &terminal) in retained.iter().enumerate() {
            check(cancelled)?;
            let mut rhs = floating_terminals
                .iter()
                .enumerate()
                .map(|(i, &t)| finite(-self.matrix_f[t * n + terminal] / scales[i]))
                .collect::<Result<Vec<_>>>()?;
            check(cancelled)?;
            factor.solve(&mut rhs);
            check(cancelled)?;
            for i in 0..f {
                transfer[i * r + column] = finite(rhs[i] / scales[i])?;
            }
        }
        let mut reduced_f = vec![0.0; r * r];
        let mut maximum_coefficient_residual_f = 0.0_f64;
        for (column, &terminal) in retained.iter().enumerate() {
            for &t in floating_terminals {
                check(cancelled)?;
                let (mut residual, mut correction) = (self.matrix_f[t * n + terminal], 0.0);
                for (j, &other) in floating_terminals.iter().enumerate() {
                    add(
                        &mut residual,
                        &mut correction,
                        finite(self.matrix_f[t * n + other] * transfer[j * r + column])?,
                    )?;
                }
                maximum_coefficient_residual_f = maximum_coefficient_residual_f.max(residual.abs());
            }
            for (row, &t) in retained.iter().enumerate() {
                check(cancelled)?;
                let (mut value, mut correction) = (self.matrix_f[t * n + terminal], 0.0);
                for (j, &other) in floating_terminals.iter().enumerate() {
                    add(
                        &mut value,
                        &mut correction,
                        finite(self.matrix_f[t * n + other] * transfer[j * r + column])?,
                    )?;
                }
                reduced_f[row * r + column] = value;
            }
        }
        for (i, &terminal) in retained.iter().enumerate() {
            let (mut gauge, mut correction) = (0.0, 0.0);
            for (j, &other) in retained.iter().enumerate() {
                check(cancelled)?;
                add(&mut gauge, &mut correction, reduced_f[i * r + j])?;
                maximum_coefficient_residual_f = maximum_coefficient_residual_f
                    .max(finite(reduced_f[i * r + j] - reduced_f[j * r + i])?.abs());
                if self.references[terminal] != self.references[other] {
                    maximum_coefficient_residual_f =
                        maximum_coefficient_residual_f.max(reduced_f[i * r + j].abs());
                }
            }
            maximum_coefficient_residual_f = maximum_coefficient_residual_f.max(gauge.abs());
        }
        if maximum_coefficient_residual_f > controls.coefficient_tolerance_f {
            return Err(invalid(
                "floating capacitance reduction, reciprocity or component gauges failed",
            ));
        }
        check(cancelled)?;
        Ok(SpatialFloatingCapacitance {
            source: self,
            retained,
            floating: floating_terminals.to_vec(),
            references,
            factor,
            scales,
            transfer,
            reduced_f,
            maximum_coefficient_residual_f,
        })
    }
}

impl<'s, 'a> SpatialFloatingCapacitance<'s, 'a> {
    pub fn source(&self) -> &'s SpatialCapacitanceMatrix<'a> {
        self.source
    }
    pub fn retained_terminals(&self) -> &[usize] {
        &self.retained
    }
    pub fn floating_terminals(&self) -> &[usize] {
        &self.floating
    }
    /// First retained voltage-reference terminal for each original terminal.
    pub fn reference_terminals(&self) -> &[usize] {
        &self.references
    }
    /// Row-major Maxwell matrix [F], in retained order. This is differential
    /// charge at fixed floating charge. Nonzero Q_f also supplies affine retained
    /// charge and positive trapped-charge storage, available through `solve`.
    pub fn reduced_matrix_f(&self) -> &[f64] {
        &self.reduced_f
    }
    pub fn maximum_coefficient_residual_f(&self) -> f64 {
        self.maximum_coefficient_residual_f
    }

    /// Solve mixed prescribed retained voltages and TOTAL floating free charges.
    /// The native response independently verifies the original field reactions
    /// and storage. Its exact geometry/material owner, charge/energy tolerances
    /// and force support remain in force; no new FE solve is performed.
    pub fn solve(
        &self,
        retained_voltage_v: &[f64],
        floating_charge_c: &[f64],
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialFloatingCapacitanceResponse<'a>> {
        check(cancelled)?;
        let r = self.retained.len();
        let n = self.source.faces.len();
        if retained_voltage_v.len() != r
            || floating_charge_c.len() != self.floating.len()
            || retained_voltage_v
                .iter()
                .chain(floating_charge_c)
                .any(|v| !v.is_finite())
            || !positive(controls.charge_tolerance_c)
            || !positive(controls.energy_tolerance_j)
        {
            return Err(invalid(
                "complete finite retained voltages, floating charges and positive physical tolerances required",
            ));
        }
        let mut voltage = vec![0.0; n];
        for (&terminal, &value) in self.retained.iter().zip(retained_voltage_v) {
            voltage[terminal] = value;
        }
        let component_reference_voltage_v: Vec<_> = self
            .references
            .iter()
            .map(|&reference| voltage[reference])
            .collect();
        for (&terminal, &value) in self.retained.iter().zip(retained_voltage_v) {
            check(cancelled)?;
            voltage[terminal] = finite(value - component_reference_voltage_v[terminal])?;
        }
        let mut charged_voltage = floating_charge_c
            .iter()
            .zip(&self.scales)
            .map(|(&charge, &scale)| finite(charge / scale))
            .collect::<Result<Vec<_>>>()?;
        check(cancelled)?;
        self.factor.solve(&mut charged_voltage);
        check(cancelled)?;
        let (mut trapped_charge_energy_j, mut trapped_correction) = (0.0, 0.0);
        for (i, &terminal) in self.floating.iter().enumerate() {
            check(cancelled)?;
            charged_voltage[i] = finite(charged_voltage[i] / self.scales[i])?;
            add(
                &mut trapped_charge_energy_j,
                &mut trapped_correction,
                finite(0.5 * floating_charge_c[i] * charged_voltage[i])?,
            )?;
            let (mut value, mut correction) = (charged_voltage[i], 0.0);
            for (j, &retained) in self.retained.iter().enumerate() {
                add(
                    &mut value,
                    &mut correction,
                    finite(self.transfer[i * r + j] * voltage[retained])?,
                )?;
            }
            voltage[terminal] = value;
        }
        let SpatialCapacitanceResponse {
            field,
            terminal_charge_c,
            terminal_energy_residual_j,
            ..
        } = self.source.response(&voltage, controls, cancelled)?;
        let mut maximum_floating_charge_residual_c = 0.0_f64;
        for (&terminal, &charge) in self.floating.iter().zip(floating_charge_c) {
            check(cancelled)?;
            maximum_floating_charge_residual_c = maximum_floating_charge_residual_c
                .max(finite(terminal_charge_c[terminal] - charge)?.abs());
        }
        let mut retained_charge_offset_c = Vec::with_capacity(r);
        let mut maximum_reduced_charge_residual_c = 0.0_f64;
        let (mut driven_energy, mut energy_correction) = (0.0, 0.0);
        for (row, &terminal) in self.retained.iter().enumerate() {
            check(cancelled)?;
            let (mut offset, mut correction) = (0.0, 0.0);
            for (i, &other) in self.floating.iter().enumerate() {
                add(
                    &mut offset,
                    &mut correction,
                    finite(self.source.matrix_f[terminal * n + other] * charged_voltage[i])?,
                )?;
            }
            let (mut driven, mut correction) = (0.0, 0.0);
            for (column, &other) in self.retained.iter().enumerate() {
                add(
                    &mut driven,
                    &mut correction,
                    finite(self.reduced_f[row * r + column] * voltage[other])?,
                )?;
            }
            maximum_reduced_charge_residual_c = maximum_reduced_charge_residual_c
                .max(finite(terminal_charge_c[terminal] - finite(driven + offset)?)?.abs());
            add(
                &mut driven_energy,
                &mut energy_correction,
                finite(0.5 * voltage[terminal] * driven)?,
            )?;
            retained_charge_offset_c.push(offset);
        }
        let reduced_energy_residual_j =
            finite(field.energy_j() - finite(driven_energy + trapped_charge_energy_j)?)?;
        if trapped_charge_energy_j < 0.0
            || maximum_floating_charge_residual_c > controls.charge_tolerance_c
            || maximum_reduced_charge_residual_c > controls.charge_tolerance_c
            || reduced_energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(invalid(
                "native floating charge, retained charge offset or reduced/trapped storage failed",
            ));
        }
        check(cancelled)?;
        Ok(SpatialFloatingCapacitanceResponse {
            field,
            component_relative_voltage_v: voltage,
            component_reference_voltage_v,
            terminal_charge_c,
            retained_charge_offset_c,
            trapped_charge_energy_j,
            maximum_floating_charge_residual_c,
            maximum_reduced_charge_residual_c,
            terminal_energy_residual_j,
            reduced_energy_residual_j,
        })
    }
}
