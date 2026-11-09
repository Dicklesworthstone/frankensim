//! Actual-amplitude sparse electrode controls on the original 3D Hodge.
//!
//! The geometry-neutral grouped scalar solve identifies each ideal electrode
//! before solving. No terminal basis fields or dense response factor are built.

use super::*;
use crate::{
    conduction::{
        ConductionError,
        terminals::{
            ConductionTerminalControl, ConductionTerminalLimits, ConductionTerminalReport,
        },
    },
    scalar_potential::{self, PotentialError},
};
use fs_math::dd::Dd;

fn invalid(reason: &'static str) -> SpatialConductionError {
    SpatialConductionError::Invalid(reason)
}

fn control_error(error: ConductionError) -> SpatialConductionError {
    match error {
        ConductionError::Invalid(reason) => invalid(reason),
        ConductionError::Cancelled => SpatialConductionError::Cancelled,
        ConductionError::Assembly(error) => SpatialConductionError::Assembly(error),
        ConductionError::Nonconvergence(report) => {
            SpatialMagneticError::Nonconvergence(report).into()
        }
        ConductionError::Balance {
            current_residual_a,
            maximum_free_current_residual_a,
            power_residual_w,
        } => SpatialConductionError::Balance {
            current_residual_a,
            maximum_free_current_residual_a,
            power_residual_w,
            hodge_power_residual_w: 0.0,
        },
    }
}

impl<C: SpatialConductorMaterial> SpatialConduction<C> {
    /// Solve simultaneous ideal voltage/current/floating electrodes at their
    /// actual amplitude with one grouped sparse native field solve.
    ///
    /// Each terminal is a disjoint union of complete original exterior faces.
    /// A shared electrode may connect separate bodies; every body and terminal
    /// must share a conducting connection through bulk, contacts or electrodes.
    /// Unlisted faces insulate. Total-current constraints determine local
    /// injection without imposing a uniform face-current density.
    ///
    /// Field voltages retain authored absolute potentials, as in `solve`.
    /// The terminal report retains differences from terminal zero separately,
    /// including drops too small to survive addition to a large common bias.
    /// All-current control uses terminal zero only as a numerical voltage gauge.
    /// Scalar/tensor bulk heat, finite-contact heat, native current/work gates
    /// and independently authored terminal current/work gates remain binding.
    pub fn solve_terminals(
        &self,
        terminals: &[&[usize]],
        control: ConductionTerminalControl<'_>,
        limits: ConductionTerminalLimits,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(SpatialConductionField, ConductionTerminalReport)> {
        self.solve_terminals_core(terminals, control, limits, None, controls, cancelled)
            .map(|(field, _, report)| (field, report))
    }

    // Both stationary and prescribed-motion controls identify the same native
    // operator. Motion contributes only its original weak nodal load; recovery
    // owns the conductor-frame heat and the separate mechanical transfer.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub(super) fn solve_terminals_core(
        &self,
        terminals: &[&[usize]],
        control: ConductionTerminalControl<'_>,
        limits: ConductionTerminalLimits,
        motion: Option<motion::SpatialConductionMotion<'_>>,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(
        SpatialConductionField,
        Option<motion::MotionPower>,
        ConductionTerminalReport,
    )> {
        check(cancelled)?;
        let n = self.mesh().vertex_count;
        let count = terminals.len();
        if count < 2 || count > n || count > limits.maximum_terminals {
            return Err(invalid("bounded physical terminal count required"));
        }
        let mut nonzeros = 0usize;
        for vertex in 0..n {
            check(cancelled)?;
            nonzeros = nonzeros
                .checked_add(self.stiffness.row(vertex).0.len())
                .ok_or_else(|| invalid("terminal workspace size overflow"))?;
        }
        let workspace = n
            .checked_mul(32)
            .and_then(|v| v.checked_add(nonzeros.checked_mul(8)?))
            .and_then(|v| v.checked_add(count.checked_mul(8)?));
        // Prepared motion retains n nodal loads and three values per cell.
        // Its DD grouped-load reduction adds two values per group, at most
        // two per original vertex. Admit all of them before preparing motion.
        let workspace = if motion.is_some() {
            workspace
                .and_then(|v| v.checked_add(n.checked_mul(3)?))
                .and_then(|v| v.checked_add(self.conductors.len().checked_mul(3)?))
        } else {
            workspace
        };
        if workspace.is_none_or(|v| v > limits.maximum_workspace_values) {
            return Err(invalid("direct terminal workspace exceeds its bound"));
        }
        let (faces, vertices) = self.electrode_faces(terminals, cancelled)?;
        // Zero values label admitted exterior geometry only. No fictitious
        // field is solved and no constitutive law is queried at a probe drive.
        let boundary = terminals::voltage_boundary(&faces, &vec![0.0; count]);
        let mut traces = self.prepare_boundary(&boundary, controls, cancelled)?;
        let authored = control
            .prepare(count, controls.current_tolerance_a, cancelled)
            .map_err(control_error)?;
        let prepared = motion
            .map(|motion| self.prepare_motion(motion, cancelled))
            .transpose()?;
        let mut groups = vec![usize::MAX; n];
        for (terminal, side) in vertices.iter().enumerate() {
            for &vertex in side {
                check(cancelled)?;
                groups[vertex] = terminal;
            }
        }
        let mut group_count = count;
        for group in &mut groups {
            check(cancelled)?;
            if *group == usize::MAX {
                *group = group_count;
                group_count += 1;
            }
        }
        let mut rhs = vec![0.0; group_count];
        let mut held = vec![None; group_count];
        for terminal in 0..count {
            check(cancelled)?;
            rhs[terminal] = authored.currents[terminal].unwrap_or(0.0);
            held[terminal] = authored.voltages[terminal];
        }
        if let Some(prepared) = &prepared {
            // P^T b_motion preserves nonuniform currents within a shared
            // electrode. The total-current equation is P^T K P V =
            // I_authored + P^T b_motion; no separate source solve is needed.
            let mut grouped = vec![Dd::from_f64(0.0); group_count];
            for (vertex, &load) in prepared.load_a.iter().enumerate() {
                check(cancelled)?;
                let group = groups[vertex];
                grouped[group] = grouped[group] + Dd::from_f64(load);
            }
            for (load, motional) in rhs.iter_mut().zip(grouped) {
                check(cancelled)?;
                *load = finite((Dd::from_f64(*load) + motional).to_f64())?;
            }
        }
        let solved = scalar_potential::solve_grouped_relative(
            &self.stiffness,
            &groups,
            &rhs,
            &held,
            controls.linear,
            cancelled,
        )
        .map_err(|error| match error {
            PotentialError::Invalid(reason) => invalid(reason),
            PotentialError::Cancelled => SpatialConductionError::Cancelled,
            PotentialError::Nonconvergence(report) => {
                SpatialMagneticError::Nonconvergence(report).into()
            }
        })?;
        let mut report = ConductionTerminalReport {
            voltage_v: Vec::with_capacity(count),
            current_a: Vec::with_capacity(count),
            supplied_power_w: 0.0,
            maximum_current_residual_a: 0.0,
            power_residual_w: 0.0,
        };
        let reference = vertices[0][0];
        let mut voltage_differences = Vec::with_capacity(count);
        for side in &vertices {
            check(cancelled)?;
            let vertex = side[0];
            let voltage = (Dd::from_f64(solved.relative[vertex])
                - Dd::from_f64(solved.relative[reference]))
                + (Dd::from_f64(solved.references[vertex])
                    - Dd::from_f64(solved.references[reference]));
            report.voltage_v.push(finite(voltage.to_f64())?);
            voltage_differences.push(voltage);
        }
        // Native recovery consumes component-relative potentials. The grouped
        // solve carries one shared reference across bodies joined by electrodes.
        // Never restore that reference before E/J, contact jumps or source work.
        for vertex in 0..n {
            check(cancelled)?;
            traces.references[self.components[vertex]] = solved.references[vertex];
        }
        if matches!(control, ConductionTerminalControl::Currents(_)) {
            traces.gauges.push(reference);
        }
        let (field, power) = self.recover_field_with_motion(
            solved.relative,
            traces,
            prepared.as_ref(),
            controls.current_tolerance_a,
            controls.power_tolerance_w,
            solved.report,
            cancelled,
        )?;
        let mut work = Dd::from_f64(0.0);
        for (terminal, side) in vertices.iter().enumerate() {
            check(cancelled)?;
            let mut current = Dd::from_f64(0.0);
            for &vertex in side {
                check(cancelled)?;
                current = current
                    + Dd::from_f64(
                        field.nodal_boundary_current_a[vertex]
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
            work = work + voltage_differences[terminal] * supplied;
            report.current_a.push(finite(current.to_f64())?);
        }
        report.supplied_power_w = finite(work.to_f64())?;
        let mechanical = power.as_ref().map_or(0.0, |power| power.mechanical_power_w);
        report.power_residual_w =
            finite(((Dd::from_f64(field.power_w) + Dd::from_f64(mechanical)) - work).to_f64())?;
        let terminal_defect = finite((Dd::from_f64(field.supplied_power_w) - work).to_f64())?;
        if report.maximum_current_residual_a > controls.current_tolerance_a
            || report.power_residual_w.abs() > controls.power_tolerance_w
            || terminal_defect.abs() > controls.power_tolerance_w
        {
            return Err(SpatialConductionError::Balance {
                current_residual_a: report.maximum_current_residual_a,
                maximum_free_current_residual_a: field.maximum_free_current_residual_a,
                power_residual_w: if terminal_defect.abs() > report.power_residual_w.abs() {
                    terminal_defect
                } else {
                    report.power_residual_w
                },
                hodge_power_residual_w: field.hodge_power_residual_w,
            });
        }
        // Motion can sustain an open-circuit voltage with exactly zero heat.
        // Nonzero prescribed current must still retain representable dissipation.
        if (authored.currents.iter().flatten().any(|&i| i != 0.0)
            || (prepared.is_none() && report.voltage_v.iter().any(|&v| v != 0.0)))
            && !terminals::positive(field.power_w)
        {
            return Err(invalid(
                "nonzero terminal drive lost representable Joule heat",
            ));
        }
        check(cancelled)?;
        Ok((field, power, report))
    }
}
