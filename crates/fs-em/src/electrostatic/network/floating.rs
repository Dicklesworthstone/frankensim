//! Floating conductors retain their prescribed TOTAL free charge while their
//! equipotential adjusts. They are not grounded or assigned zero local D.n.
//!
//! Partitioning Q=C V into retained r and floating f gives
//! V_f=C_ff^-1(Q_f-C_fr V_r) and C_reduced=C_rr-C_rf C_ff^-1 C_fr.
//! A nonzero Q_f also induces an affine retained-charge offset; the reduced
//! matrix alone must never replace that charged physical state. References:
//! [Ansys Q3D: grounding and floating conductors](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/Q3DExtractor/Content/Q3D/TNMatrixReductionGroundFloat.htm)
//! and [COMSOL: floating potential](https://doc.comsol.com/6.4/doc/com.comsol.help.semicond/semicond_ug_acdc.7.20.html).

use super::{PlanarCapacitanceMatrix, add, invalid};
use crate::electrostatic::{ElectricField, ElectrostaticError, finite, poll, positive};
use fs_la::factor::{Cholesky, cholesky};

/// Bounds for the dense terminal reduction, additional to field extraction.
#[derive(Clone, Copy, Debug)]
pub struct FloatingCapacitanceControls {
    pub maximum_floating_terminals: usize,
    /// Conservative scalar-product work: f^3 + 2*r*f^2 + r^2*f + n^2.
    /// The native dense Cholesky operation is bounded by this allowance;
    /// cancellation is checked immediately before and after it.
    pub maximum_dense_work: usize,
    /// Absolute [F] allowance for reduction, reciprocity and gauge residuals.
    pub coefficient_tolerance_f: f64,
}

/// Reusable charge-constrained reduction retaining the original material field.
/// Preparation factors only the floating terminal block; subsequent states
/// require triangular solves and original field recovery, not new FE solves.
pub struct FloatingCapacitance<'s, 'm> {
    source: &'s PlanarCapacitanceMatrix<'m>,
    retained: Vec<usize>,
    floating: Vec<usize>,
    factor: Cholesky,
    scales: Vec<f64>,
    transfer: Vec<f64>,
    reduced_f: Vec<f64>,
    maximum_coefficient_residual_f: f64,
}

/// One actual-amplitude mixed voltage/charge solution on all original terminals.
pub struct FloatingCapacitanceResponse<'m> {
    /// Complete original terminal order, relative to `reference_terminal`.
    /// Keeping the reference separate avoids losing small floating voltages
    /// when a large, physically irrelevant common potential is supplied.
    pub terminal_voltage_v: Vec<f64>,
    pub reference_terminal: usize,
    pub reference_voltage_v: f64,
    /// Independently recovered original weak electrode reactions, in coulombs.
    pub terminal_charge_c: Vec<f64>,
    pub maximum_floating_charge_residual_c: f64,
    /// Original field energy minus one half of all terminal voltage-charge work.
    pub terminal_energy_residual_j: f64,
    field: ElectricField<'m>,
}

impl<'m> FloatingCapacitanceResponse<'m> {
    /// Native spatial field, including material identity, energy and force API.
    /// Its nodal voltage reference remains the original FIRST terminal, as in
    /// `PlanarCapacitanceMatrix::field_at_voltages`.
    pub fn field(&self) -> &ElectricField<'m> {
        &self.field
    }
}

impl<'m> PlanarCapacitanceMatrix<'m> {
    /// Prepare the physical Schur reduction for a nonempty proper subset of
    /// floating electrodes. Their order fixes subsequent charge-input order;
    /// retained terminals keep their original relative order.
    ///
    /// Every dielectric component, including components joined by a common
    /// equipotential terminal, needs a retained voltage reference. An entirely
    /// floating component refuses BEFORE factorization: roundoff must not turn
    /// its gauge nullspace into an apparently invertible capacitance.
    /// The frozen field can contain anisotropic dielectrics and positive mutual
    /// entries; a passive pair-capacitor realization is not required here.
    pub fn prepare_floating(
        &self,
        floating_terminals: &[usize],
        controls: FloatingCapacitanceControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FloatingCapacitance<'_, 'm>, ElectrostaticError> {
        poll(cancelled)?;
        let n = self.terminals.len();
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
        admit_references(self, &retained, cancelled)?;
        let scales: Vec<_> = floating_terminals
            .iter()
            .map(|&t| {
                let value = self.matrix_f[t * n + t];
                if positive(value) {
                    Ok(value.sqrt())
                } else {
                    Err(invalid(
                        "floating electrode has no positive electric storage",
                    ))
                }
            })
            .collect::<Result<_, _>>()?;
        // Diagonal equilibration avoids an absolute pivot threshold in farads.
        // Picofarad devices and large capacitors follow the same equations.
        let mut block = vec![0.0; f * f];
        for (i, &a) in floating_terminals.iter().enumerate() {
            poll(cancelled)?;
            for (j, &b) in floating_terminals.iter().enumerate() {
                block[i * f + j] = finite((self.matrix_f[a * n + b] / scales[i]) / scales[j])?;
            }
        }
        poll(cancelled)?;
        let factor = cholesky(&block, f)
            .map_err(|_| invalid("floating capacitance block is not positive definite"))?;
        poll(cancelled)?;
        let mut transfer = vec![0.0; f * r];
        for (column, &terminal) in retained.iter().enumerate() {
            poll(cancelled)?;
            let mut rhs = floating_terminals
                .iter()
                .enumerate()
                .map(|(i, &t)| finite(-self.matrix_f[t * n + terminal] / scales[i]))
                .collect::<Result<Vec<_>, _>>()?;
            factor.solve(&mut rhs);
            for i in 0..f {
                transfer[i * r + column] = finite(rhs[i] / scales[i])?;
            }
        }
        let mut reduced_f = vec![0.0; r * r];
        let mut maximum_coefficient_residual_f: f64 = 0.0;
        for (column, &terminal) in retained.iter().enumerate() {
            for &t in floating_terminals {
                poll(cancelled)?;
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
                poll(cancelled)?;
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
        for i in 0..r {
            let (mut gauge, mut correction) = (0.0, 0.0);
            for j in 0..r {
                poll(cancelled)?;
                add(&mut gauge, &mut correction, reduced_f[i * r + j])?;
                maximum_coefficient_residual_f = maximum_coefficient_residual_f
                    .max(finite(reduced_f[i * r + j] - reduced_f[j * r + i])?.abs());
            }
            maximum_coefficient_residual_f = maximum_coefficient_residual_f.max(gauge.abs());
        }
        if maximum_coefficient_residual_f > controls.coefficient_tolerance_f {
            return Err(invalid(
                "floating capacitance reduction, reciprocity or gauge failed",
            ));
        }
        poll(cancelled)?;
        Ok(FloatingCapacitance {
            source: self,
            retained,
            floating: floating_terminals.to_vec(),
            factor,
            scales,
            transfer,
            reduced_f,
            maximum_coefficient_residual_f,
        })
    }
}

impl<'s, 'm> FloatingCapacitance<'s, 'm> {
    pub fn source(&self) -> &'s PlanarCapacitanceMatrix<'m> {
        self.source
    }
    pub fn retained_terminals(&self) -> &[usize] {
        &self.retained
    }
    pub fn floating_terminals(&self) -> &[usize] {
        &self.floating
    }
    /// Row-major Maxwell matrix on retained electrodes, in farads. It describes
    /// DIFFERENTIAL charge at fixed floating charge, or total charge when every
    /// floating charge is zero. Use `solve` to retain nonzero-charge offsets and
    /// their energy. A zero-voltage charged device can still store energy.
    pub fn reduced_matrix_f(&self) -> &[f64] {
        &self.reduced_f
    }
    pub fn maximum_coefficient_residual_f(&self) -> f64 {
        self.maximum_coefficient_residual_f
    }

    /// Solve mixed prescribed voltages and total floating charges. The original
    /// native field independently checks actual charge and energy, including
    /// errors in the extracted basis and scaled dense solve. Failed admission
    /// or cancellation returns no partial state and leaves preparation reusable.
    ///
    /// The result supplies all original electrode voltages for initializing the
    /// existing native capacitor network; omitting the floating nodes from that
    /// circuit without their charge constraints would change the physics.
    pub fn solve(
        &self,
        retained_voltage_v: &[f64],
        floating_charge_c: &[f64],
        charge_tolerance_c: f64,
        energy_tolerance_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FloatingCapacitanceResponse<'m>, ElectrostaticError> {
        poll(cancelled)?;
        let r = self.retained.len();
        if retained_voltage_v.len() != r
            || floating_charge_c.len() != self.floating.len()
            || retained_voltage_v
                .iter()
                .chain(floating_charge_c)
                .any(|v| !v.is_finite())
            || !positive(charge_tolerance_c)
            || !positive(energy_tolerance_j)
        {
            return Err(invalid(
                "complete finite retained voltages, floating charges and positive physical tolerances required",
            ));
        }
        let reference_voltage_v = retained_voltage_v[0];
        let differences = retained_voltage_v
            .iter()
            .map(|&v| finite(v - reference_voltage_v))
            .collect::<Result<Vec<_>, _>>()?;
        let mut charged_voltage = floating_charge_c
            .iter()
            .zip(&self.scales)
            .map(|(&q, &s)| finite(q / s))
            .collect::<Result<Vec<_>, _>>()?;
        self.factor.solve(&mut charged_voltage);
        poll(cancelled)?;
        let mut voltage = vec![0.0; self.source.terminals.len()];
        for (&t, &v) in self.retained.iter().zip(&differences) {
            voltage[t] = v;
        }
        for (i, &t) in self.floating.iter().enumerate() {
            poll(cancelled)?;
            let (mut value, mut correction) = (finite(charged_voltage[i] / self.scales[i])?, 0.0);
            for (j, &v) in differences.iter().enumerate() {
                add(
                    &mut value,
                    &mut correction,
                    finite(self.transfer[i * r + j] * v)?,
                )?;
            }
            voltage[t] = value;
        }
        let field = self.source.field_at_voltages(
            &voltage,
            charge_tolerance_c,
            energy_tolerance_j,
            cancelled,
        )?;
        let mut terminal_charge_c = Vec::with_capacity(voltage.len());
        for terminal in &self.source.terminals {
            let (mut charge, mut correction) = (0.0, 0.0);
            for &v in terminal {
                poll(cancelled)?;
                add(
                    &mut charge,
                    &mut correction,
                    field.nodal_boundary_charge_c()[v].ok_or_else(|| {
                        invalid("floating terminal lost its native charge reaction")
                    })?,
                )?;
            }
            terminal_charge_c.push(charge);
        }
        let mut maximum_floating_charge_residual_c: f64 = 0.0;
        for (&t, &charge) in self.floating.iter().zip(floating_charge_c) {
            maximum_floating_charge_residual_c = maximum_floating_charge_residual_c
                .max(finite(terminal_charge_c[t] - charge)?.abs());
        }
        let (mut work, mut correction) = (0.0, 0.0);
        for (&v, &q) in voltage.iter().zip(&terminal_charge_c) {
            poll(cancelled)?;
            add(&mut work, &mut correction, finite(0.5 * v * q)?)?;
        }
        let terminal_energy_residual_j = finite(field.energy_j() - work)?;
        if maximum_floating_charge_residual_c > charge_tolerance_c
            || terminal_energy_residual_j.abs() > energy_tolerance_j
        {
            return Err(invalid(
                "native floating charge or terminal-field energy failed",
            ));
        }
        poll(cancelled)?;
        Ok(FloatingCapacitanceResponse {
            terminal_voltage_v: voltage,
            reference_terminal: self.retained[0],
            reference_voltage_v,
            terminal_charge_c,
            maximum_floating_charge_residual_c,
            terminal_energy_residual_j,
            field,
        })
    }
}

// Physical connectivity, rather than a numerical pivot threshold, owns gauge
// admission. A terminal can itself join disjoint dielectric mesh components.
fn admit_references(
    source: &PlanarCapacitanceMatrix<'_>,
    retained: &[usize],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ElectrostaticError> {
    let mesh = source.problem.mesh();
    let terminals: Vec<_> = source.terminals.iter().map(Vec::as_slice).collect();
    let parent = crate::scalar_potential::terminal_components(mesh, &terminals, cancelled)?;
    let mut anchored = vec![false; parent.len()];
    for &terminal in retained {
        anchored[parent[source.terminals[terminal][0]]] = true;
    }
    for &component in &parent {
        poll(cancelled)?;
        if !anchored[component] {
            return Err(invalid(
                "every floating dielectric component needs a retained voltage reference",
            ));
        }
    }
    Ok(())
}
