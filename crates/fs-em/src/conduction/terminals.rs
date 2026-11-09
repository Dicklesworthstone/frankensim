//! Direct simultaneous ideal-electrode controls on the native sparse operator.
//!
//! Each electrode contributes one equipotential coordinate. The original
//! stiffness is identified under that constraint and solved once at the actual
//! drive, without a unit-voltage probe or a dense terminal-response matrix.
//! Native cell/contact recovery and independent terminal work remain binding.

use super::*;
use fs_math::dd::Dd;

/// Held controls in the original terminal order; currents enter the conductor.
#[derive(Clone, Copy, Debug)]
pub enum ConductionTerminalControl<'a> {
    /// Complete finite electrode potentials [V].
    Voltages(&'a [f64]),
    /// Complete conserved total currents [A]. Terminal zero fixes the gauge.
    /// Zero current floats an equipotential electrode, permitting face transfer.
    Currents(&'a [f64]),
    /// Current values follow the supplied unique proper nonempty subset.
    /// Voltage values follow its complement in original terminal order.
    Mixed {
        current_terminals: &'a [usize],
        currents_a: &'a [f64],
        voltages_v: &'a [f64],
    },
}

/// Admission for the extra sparse identification/solve workspace.
/// Existing mesh, material and checked-CG bounds remain independently active.
#[derive(Clone, Copy, Debug)]
pub struct ConductionTerminalLimits {
    pub maximum_terminals: usize,
    /// Conservative scalar/index count: 32*vertices + 8*original_nonzeros
    /// + 8*terminals. Prescribed-motion methods also charge their documented
    /// preparation buffers. This is not an allocator byte-peak guarantee.
    pub maximum_workspace_values: usize,
}

/// Actual ordered electrode observations and independently authored source work.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductionTerminalReport {
    /// Physical voltage differences [V] from original terminal zero. Small
    /// differences survive even when absolute field observations share a huge bias.
    pub voltage_v: Vec<f64>,
    /// Native total inward electrode reactions [A], in original order.
    pub current_a: Vec<f64>,
    /// Authored current controls times solved differences, plus voltage-source
    /// reactions times held differences, with each electrode counted once [W].
    pub supplied_power_w: f64,
    pub maximum_current_residual_a: f64,
    /// Native bulk/contact heat minus independent supplied power [W]. With
    /// prescribed motion, include signed mechanical output alongside heat.
    pub power_residual_w: f64,
}

fn invalid(reason: &'static str) -> ConductionError {
    ConductionError::Invalid(reason)
}

pub(crate) struct Controls {
    pub(crate) voltages: Vec<Option<f64>>,
    pub(crate) currents: Vec<Option<f64>>,
}

impl ConductionTerminalControl<'_> {
    pub(crate) fn prepare(
        self,
        count: usize,
        tolerance_a: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Controls, ConductionError> {
        let mut result = Controls {
            voltages: vec![None; count],
            currents: vec![None; count],
        };
        match self {
            Self::Voltages(values) => {
                if values.len() != count || values.iter().any(|v| !v.is_finite()) {
                    return Err(invalid("complete finite electrode voltages required"));
                }
                for (out, &value) in result.voltages.iter_mut().zip(values) {
                    poll(cancelled)?;
                    *out = Some(value);
                }
            }
            Self::Currents(values) => {
                if values.len() != count || values.iter().any(|v| !v.is_finite()) {
                    return Err(invalid("complete finite electrode currents required"));
                }
                let mut sum = Dd::from_f64(0.0);
                for (out, &value) in result.currents.iter_mut().zip(values) {
                    poll(cancelled)?;
                    *out = Some(value);
                    sum = sum + Dd::from_f64(value);
                }
                if finite(sum.to_f64())?.abs() > tolerance_a {
                    return Err(invalid("terminal currents must sum to zero"));
                }
                result.voltages[0] = Some(0.0);
            }
            Self::Mixed {
                current_terminals,
                currents_a,
                voltages_v,
            } => {
                if current_terminals.is_empty()
                    || current_terminals.len() >= count
                    || currents_a.len() != current_terminals.len()
                    || voltages_v.len() != count - current_terminals.len()
                    || currents_a.iter().chain(voltages_v).any(|v| !v.is_finite())
                {
                    return Err(invalid(
                        "complete mixed controls and a proper current subset required",
                    ));
                }
                for (&terminal, &current) in current_terminals.iter().zip(currents_a) {
                    poll(cancelled)?;
                    if terminal >= count || result.currents[terminal].is_some() {
                        return Err(invalid("unique in-range current terminals required"));
                    }
                    result.currents[terminal] = Some(current);
                }
                let mut voltages = voltages_v.iter();
                for terminal in 0..count {
                    poll(cancelled)?;
                    if result.currents[terminal].is_none() {
                        result.voltages[terminal] = voltages.next().copied();
                    }
                }
            }
        }
        Ok(result)
    }
}

impl<C: PlanarOhmicResponse> PlanarConduction<C> {
    /// Solve simultaneous voltage/current/floating electrodes once at their
    /// actual amplitude, using sparse identification of equipotential vertices.
    ///
    /// Terminals are disjoint unions of complete exterior traces. A terminal
    /// may join disconnected physical pieces. All terminals must share one
    /// conducting connection through bulk, contacts or ideal electrodes; every
    /// body needs an admitted connection to those controls. Unlisted surfaces
    /// insulate. Each axisymmetric electrode requires positive physical area.
    ///
    /// A current control specifies total reaction, never a uniform boundary
    /// density. Original scalar/tensor material and contact laws recover local
    /// current and heat. Both authored terminal currents and source work are
    /// checked independently after the original algebraic/current/work gates.
    pub fn solve_terminals(
        &self,
        terminals: &[&[usize]],
        control: ConductionTerminalControl<'_>,
        limits: ConductionTerminalLimits,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(ConductionField<'_, C>, ConductionTerminalReport), ConductionError> {
        poll(cancelled)?;
        let n = self.mesh.vertices().len();
        let count = terminals.len();
        if count < 2 || count > n || count > limits.maximum_terminals {
            return Err(invalid("bounded physical terminal count required"));
        }
        let mut nonzeros = 0usize;
        for vertex in 0..n {
            poll(cancelled)?;
            nonzeros = nonzeros
                .checked_add(self.stiffness.row(vertex).0.len())
                .ok_or_else(|| invalid("terminal workspace size overflow"))?;
        }
        let workspace = n
            .checked_mul(32)
            .and_then(|v| v.checked_add(nonzeros.checked_mul(8)?))
            .and_then(|v| v.checked_add(count.checked_mul(8)?));
        if workspace.is_none_or(|v| v > limits.maximum_workspace_values) {
            return Err(invalid("direct terminal workspace exceeds its bound"));
        }
        self.validate_controls(0, controls)?;
        let authored = control.prepare(count, controls.current_tolerance_a, cancelled)?;
        if self.mesh.metric().is_axisymmetric() {
            axisymmetric::admit_terminal_traces(&self.mesh, terminals, cancelled)?;
        }
        // Complete trace admission uses temporary zero potentials only as
        // geometry labels; no material is queried at a fictitious field.
        let boundary = scalar_potential::terminal_boundary(
            &self.mesh,
            terminals,
            &vec![0.0; count],
            cancelled,
        )?;
        let mut prescribed = self.prescribed_boundary(&boundary, cancelled)?;
        let contact_links: Vec<_> = self
            .contacts
            .iter()
            .map(|contact| [contact.first_vertices[0], contact.second_vertices[0]])
            .collect();
        let mut links = terminals.to_vec();
        links.extend(contact_links.iter().map(|pair| pair.as_slice()));
        let components = scalar_potential::terminal_components(&self.mesh, &links, cancelled)?;
        if terminals
            .iter()
            .any(|terminal| components[terminal[0]] != components[terminals[0][0]])
        {
            return Err(invalid(
                "terminals lack a common bulk/contact/electrode connection",
            ));
        }

        let mut groups = vec![usize::MAX; n];
        for (terminal, vertices) in terminals.iter().enumerate() {
            for &vertex in *vertices {
                poll(cancelled)?;
                groups[vertex] = terminal;
            }
        }
        let mut group_count = count;
        for group in &mut groups {
            poll(cancelled)?;
            if *group == usize::MAX {
                *group = group_count;
                group_count += 1;
            }
        }
        let mut rhs = vec![0.0; group_count];
        let mut held = vec![None; group_count];
        for terminal in 0..count {
            poll(cancelled)?;
            rhs[terminal] = authored.currents[terminal].unwrap_or(0.0);
            held[terminal] = authored.voltages[terminal];
        }
        let solved = scalar_potential::solve_grouped_relative(
            &self.stiffness,
            &groups,
            &rhs,
            &held,
            controls.linear,
            cancelled,
        )?;
        for (terminal, vertices) in terminals.iter().enumerate() {
            for &vertex in *vertices {
                poll(cancelled)?;
                prescribed[vertex] = Some(match authored.voltages[terminal] {
                    Some(value) => value,
                    None => finite(solved.references[vertex] + solved.relative[vertex])?,
                });
            }
        }
        let field = self.recover_field(
            solved.relative,
            solved.references,
            &prescribed,
            vec![None; n],
            controls.current_tolerance_a,
            controls.power_tolerance_w,
            solved.report,
            cancelled,
        )?;
        let reference = terminals[0][0];
        let relative = field.relative_voltage_v();
        let references = field.reference_voltage_v();
        let mut report = ConductionTerminalReport {
            voltage_v: Vec::with_capacity(count),
            current_a: Vec::with_capacity(count),
            supplied_power_w: 0.0,
            maximum_current_residual_a: 0.0,
            power_residual_w: 0.0,
        };
        let mut work = Dd::from_f64(0.0);
        for (terminal, vertices) in terminals.iter().enumerate() {
            poll(cancelled)?;
            let vertex = vertices[0];
            let voltage = (Dd::from_f64(relative[vertex]) - Dd::from_f64(relative[reference]))
                + (Dd::from_f64(references[vertex]) - Dd::from_f64(references[reference]));
            let mut current = Dd::from_f64(0.0);
            for &vertex in *vertices {
                poll(cancelled)?;
                current = current
                    + Dd::from_f64(
                        field.nodal_boundary_current_a()[vertex]
                            .ok_or_else(|| invalid("terminal lost its native current reaction"))?,
                    );
            }
            let supplied = match authored.currents[terminal] {
                Some(target) => {
                    let defect = finite((current - Dd::from_f64(target)).to_f64())?;
                    report.maximum_current_residual_a =
                        report.maximum_current_residual_a.max(defect.abs());
                    Dd::from_f64(target)
                }
                None => current,
            };
            work = work + voltage * supplied;
            report.voltage_v.push(finite(voltage.to_f64())?);
            report.current_a.push(finite(current.to_f64())?);
        }
        report.supplied_power_w = finite(work.to_f64())?;
        report.power_residual_w = finite((Dd::from_f64(field.power_w()) - work).to_f64())?;
        if report.maximum_current_residual_a > controls.current_tolerance_a
            || report.power_residual_w.abs() > controls.power_tolerance_w
        {
            return Err(ConductionError::Balance {
                current_residual_a: report.maximum_current_residual_a,
                maximum_free_current_residual_a: field.maximum_free_current_residual_a(),
                power_residual_w: report.power_residual_w,
            });
        }
        if (authored.currents.iter().flatten().any(|&i| i != 0.0)
            || report.voltage_v.iter().any(|&v| v != 0.0))
            && !positive(field.power_w())
        {
            return Err(invalid(
                "nonzero terminal drive lost representable Joule heat",
            ));
        }
        poll(cancelled)?;
        Ok((field, report))
    }
}
