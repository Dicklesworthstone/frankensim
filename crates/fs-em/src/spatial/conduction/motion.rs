//! Prescribed three-dimensional motion in the original steady-current field.
//!
//! `J = Sigma(E + u cross B)`, with `E = -grad(V)`. The motional weak load
//! changes charge continuity; the original conductivity operator stays SPD.
//! Conductor-frame Joule heat and Lorentz output power remain distinct:
//! `integral(E.J) = integral((E + u cross B).J) + integral(u.(J cross B))`.
//!
//! Velocity and magnetic flux density are constant per original tetrahedron.
//! This fixed-domain, nonrelativistic DC approximation prescribes B and motion;
//! it does not evolve induction, finite moving interfaces or mechanical state.

use super::terminals::{terminal_current, voltage_boundary};
use super::*;
use crate::conduction::terminals::{
    ConductionTerminalControl, ConductionTerminalLimits, ConductionTerminalReport,
};

mod mechanical;
pub use crate::spatial::vector_potential::induction::mechanical::{
    SpatialLorentzLoadError, SpatialLorentzLoadLimits, SpatialLorentzLoads,
};

/// Laboratory-frame prescribed fields on the original tetrahedra. A fixed
/// Eulerian mesh requires motion to preserve the represented material/source
/// occupancy, or a separately justified instantaneous quasistatic snapshot.
#[derive(Clone, Copy, Debug)]
pub struct SpatialConductionMotion<'a> {
    pub velocity_m_s: &'a [[f64; 3]],
    pub flux_density_t: &'a [[f64; 3]],
}

pub(super) struct PreparedMotion<'a> {
    pub(super) prescribed: SpatialConductionMotion<'a>,
    pub(super) field_v_m: Vec<[f64; 3]>,
    pub(super) load_a: Vec<f64>,
}

/// One zero-electrode-voltage solve supplies the affine motional source. Its
/// relative potentials share the original operator and component gauges with
/// every retained homogeneous terminal basis.
pub(super) struct TerminalMotion<'a> {
    prepared: PreparedMotion<'a>,
    relative: Vec<f64>,
    pub(super) currents_a: Vec<f64>,
}

fn cross(a: [f64; 3], b: [f64; 3]) -> Result<[f64; 3]> {
    Ok([
        finite(a[1].mul_add(b[2], -a[2] * b[1]))?,
        finite(a[2].mul_add(b[0], -a[0] * b[2]))?,
        finite(a[0].mul_add(b[1], -a[1] * b[0]))?,
    ])
}

impl<C: SpatialConductorMaterial> SpatialConduction<C> {
    /// Solve the original terminal/current-trace problem with prescribed
    /// `u cross B`. Open components retain their native zero-current gauge:
    /// a conservative motional field is screened by voltage, not an invented
    /// short circuit. Scalar and full rotated conductivity tensors share this
    /// same solve and source ownership.
    ///
    /// Independent charge, local work, Hodge and boundary-power gates precede
    /// publication. Finite contacts retain their stationary surface law only
    /// when both adjacent original cells have zero prescribed velocity. Moving
    /// bulk elsewhere is admitted; moving/sliding interfaces still refuse.
    /// Surface electrical work is heat, with zero surface mechanical work.
    /// The caller owns magnetic consistency, motion and thermal evolution
    /// outside this quasistatic snapshot.
    pub fn solve_with_motion(
        &self,
        boundary: &SpatialConductionBoundary,
        motion: SpatialConductionMotion<'_>,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialMovingConductionField<'_, C>> {
        check(cancelled)?;
        let prepared = self.prepare_motion(motion, cancelled)?;
        let traces = self.prepare_boundary(boundary, controls, cancelled)?;
        let (relative, report) = self.solve_motion(&traces, &prepared, controls, cancelled)?;
        self.recover_moving_field(relative, traces, &prepared, controls, report, cancelled)
    }

    /// Solve simultaneous voltage/total-current/floating electrodes with the
    /// original prescribed motional weak load in one grouped sparse solve.
    ///
    /// All electrode and sparse-workspace admission is shared with
    /// `solve_terminals`. In addition to its stationary workspace, motion
    /// requires at most `3*vertices + 3*tetrahedra` scalar values for its native
    /// load, per-cell field and compensated group reduction. This bound is
    /// checked before preparing those arrays.
    ///
    /// The report retains voltage differences from terminal zero independently
    /// of absolute displayed voltages. Its signed supplied work uses authored
    /// current controls and native voltage-source reactions, and must equal
    /// conductor-frame Joule heat plus signed Lorentz output work. Generation
    /// may have negative terminal work; an open circuit may have nonzero
    /// voltage with zero current and heat. Stationary-contact admission,
    /// scalar/tensor source ownership and all native local work gates remain
    /// those of `solve_with_motion`; no moving-interface law is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_terminals_with_motion(
        &self,
        terminals: &[&[usize]],
        control: ConductionTerminalControl<'_>,
        limits: ConductionTerminalLimits,
        motion: SpatialConductionMotion<'_>,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(
        SpatialMovingConductionField<'_, C>,
        ConductionTerminalReport,
    )> {
        let (field, power, terminals) = self.solve_terminals_core(
            terminals,
            control,
            limits,
            Some(motion),
            controls,
            cancelled,
        )?;
        let power = power.ok_or(SpatialConductionError::Invalid(
            "missing original motional power recovery",
        ))?;
        check(cancelled)?;
        Ok((
            SpatialMovingConductionField {
                problem: self,
                velocity_m_s: motion.velocity_m_s.to_vec(),
                flux_density_t: motion.flux_density_t.to_vec(),
                field,
                power,
            },
            terminals,
        ))
    }

    pub(super) fn prepare_motion<'v>(
        &self,
        motion: SpatialConductionMotion<'v>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMotion<'v>> {
        check(cancelled)?;
        let cells = self.conductors.len();
        if motion.velocity_m_s.len() != cells || motion.flux_density_t.len() != cells {
            return Err(SpatialConductionError::Invalid(
                "one prescribed velocity and magnetic field per original tetrahedron required",
            ));
        }
        for cells in self.contact_cells() {
            check(cancelled)?;
            if cells
                .iter()
                .any(|&cell| motion.velocity_m_s[cell] != [0.0; 3])
            {
                return Err(SpatialConductionError::Invalid(
                    "moving electrical contacts require an owned surface work and heat law",
                ));
            }
        }
        let g = &self.geometry;
        let mut prepared = PreparedMotion {
            prescribed: motion,
            field_v_m: Vec::with_capacity(cells),
            load_a: vec![0.0; g.mesh.vertex_count],
        };
        for (cell, tet) in g.mesh.tets.iter().enumerate() {
            check(cancelled)?;
            for value in motion.velocity_m_s[cell]
                .into_iter()
                .chain(motion.flux_density_t[cell])
            {
                finite(value)?;
            }
            let motional = cross(motion.velocity_m_s[cell], motion.flux_density_t[cell])?;
            prepared.field_v_m.push(motional);
            if motional == [0.0; 3] {
                continue;
            }
            let current = self.conductivity[cell].current_density_a_m2(motional);
            for value in current {
                finite(value)?;
            }
            let volume = g.geometry.vol_signed[cell].abs();
            let mut local_load = [0.0; 4];
            for (local, load) in local_load.iter_mut().enumerate().skip(1) {
                *load = finite(volume * dot(g.geometry.grads[cell][local], current))?;
            }
            // Exact partition of unity in the same reference-node convention
            // as electric-field recovery: motion cannot inject net charge.
            local_load[0] = finite(-finite(
                finite(local_load[1] + local_load[2])? + local_load[3],
            )?)?;
            for (&vertex, load) in tet.iter().zip(local_load) {
                let vertex = vertex as usize;
                prepared.load_a[vertex] = finite(prepared.load_a[vertex] + load)?;
            }
        }
        Ok(prepared)
    }

    fn solve_motion(
        &self,
        traces: &ScalarBoundary,
        prepared: &PreparedMotion<'_>,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, Option<SolveReport>)> {
        let mut rhs = traces.flux_load_wb.clone();
        for (value, &motional) in rhs.iter_mut().zip(&prepared.load_a) {
            check(cancelled)?;
            if motional != 0.0 {
                *value = finite(*value + motional)?;
            }
        }
        Ok(solve_scalar_system(
            &self.stiffness,
            &rhs,
            &traces.prescribed,
            controls.linear,
            cancelled,
        )?)
    }

    fn recover_moving_field(
        &self,
        relative: Vec<f64>,
        traces: ScalarBoundary,
        prepared: &PreparedMotion<'_>,
        controls: ConductionControls,
        report: Option<SolveReport>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialMovingConductionField<'_, C>> {
        let (field, power) = self.recover_field_with_motion(
            relative,
            traces,
            Some(prepared),
            controls.current_tolerance_a,
            controls.power_tolerance_w,
            report,
            cancelled,
        )?;
        let power = power.ok_or(SpatialConductionError::Invalid(
            "missing original motional power recovery",
        ))?;
        check(cancelled)?;
        Ok(SpatialMovingConductionField {
            problem: self,
            velocity_m_s: prepared.prescribed.velocity_m_s.to_vec(),
            flux_density_t: prepared.prescribed.flux_density_t.to_vec(),
            field,
            power,
        })
    }

    /// The same native solve as the boundary-driven path, with every actual
    /// electrode clamped to zero. This response is a computational source;
    /// its heat and Lorentz work are never added to the accepted actual field.
    pub(super) fn terminal_motion<'v>(
        &self,
        faces: &[Vec<usize>],
        vertices: &[Vec<usize>],
        motion: SpatialConductionMotion<'v>,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TerminalMotion<'v>> {
        let prepared = self.prepare_motion(motion, cancelled)?;
        let boundary = voltage_boundary(faces, &vec![0.0; faces.len()]);
        let traces = self.prepare_boundary(&boundary, controls, cancelled)?;
        let (relative, report) = self.solve_motion(&traces, &prepared, controls, cancelled)?;
        let field = self.recover_moving_field(
            relative.clone(),
            traces,
            &prepared,
            controls,
            report,
            cancelled,
        )?;
        let currents_a = vertices
            .iter()
            .map(|vertices| terminal_current(&field.field, vertices, cancelled))
            .collect::<Result<Vec<_>>>()?;
        check(cancelled)?;
        Ok(TerminalMotion {
            prepared,
            relative,
            currents_a,
        })
    }

    /// Add ONLY potentials, then independently recover original local currents,
    /// bulk/contact heat and Lorentz work. All cross terms remain present.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn recover_terminal_motion(
        &self,
        mut relative: Vec<f64>,
        source: &TerminalMotion<'_>,
        faces: &[Vec<usize>],
        vertices: &[Vec<usize>],
        voltages: &[f64],
        expected_currents: &[f64],
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialMovingConductionField<'_, C>> {
        let boundary = voltage_boundary(faces, voltages);
        let traces = self.prepare_boundary(&boundary, controls, cancelled)?;
        for (v, value) in relative.iter_mut().enumerate() {
            check(cancelled)?;
            *value = match traces.prescribed[v] {
                Some(prescribed) => prescribed,
                None => finite(*value + source.relative[v])?,
            };
        }
        // The source solve's report is not a report for this superposed field.
        let field = self.recover_moving_field(
            relative,
            traces,
            &source.prepared,
            controls,
            None,
            cancelled,
        )?;
        let mut current_residual_a = 0.0_f64;
        let (mut work, mut correction) = (0.0, 0.0);
        for ((vertices, &current), &voltage) in vertices.iter().zip(expected_currents).zip(voltages)
        {
            check(cancelled)?;
            let defect = finite(terminal_current(&field.field, vertices, cancelled)? - current)?;
            if defect.abs() > current_residual_a.abs() {
                current_residual_a = defect;
            }
            let adjusted = finite(finite(voltage * current)? - correction)?;
            let next = finite(work + adjusted)?;
            correction = finite(finite(next - work)? - adjusted)?;
            work = next;
        }
        // Signed authored terminal work may vanish at short circuit or be
        // negative during generation even while conductor-frame heat is positive.
        let terminal_defect = finite(field.terminal_power_w() - work)?;
        let physical_defect =
            finite(finite(field.joule_power_w() + field.mechanical_power_w())? - work)?;
        let power_residual_w = if terminal_defect.abs() > physical_defect.abs() {
            terminal_defect
        } else {
            physical_defect
        };
        if current_residual_a.abs() > controls.current_tolerance_a
            || power_residual_w.abs() > controls.power_tolerance_w
        {
            return Err(SpatialConductionError::Balance {
                current_residual_a,
                maximum_free_current_residual_a: field.maximum_free_current_residual_a(),
                power_residual_w,
                hodge_power_residual_w: field.hodge_power_residual_w(),
            });
        }
        check(cancelled)?;
        Ok(field)
    }
}

#[derive(Clone, Debug)]
pub(super) struct MotionPower {
    conductor_field_v_m: Vec<[f64; 3]>,
    lorentz_force_density_n_m3: Vec<[f64; 3]>,
    cell_mechanical_power_w: Vec<f64>,
    pub(super) mechanical_power_w: f64,
    pub(super) electrical_power_w: f64,
    pub(super) contact_power_w: f64,
    pub(super) work_pairing_residual_w: f64,
}

impl MotionPower {
    pub(super) fn new(cells: usize) -> Self {
        Self {
            conductor_field_v_m: Vec::with_capacity(cells),
            lorentz_force_density_n_m3: Vec::with_capacity(cells),
            cell_mechanical_power_w: Vec::with_capacity(cells),
            mechanical_power_w: 0.0,
            electrical_power_w: 0.0,
            contact_power_w: 0.0,
            work_pairing_residual_w: 0.0,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn record(
        &mut self,
        motion: SpatialConductionMotion<'_>,
        cell: usize,
        electric: [f64; 3],
        conductor: [f64; 3],
        current: [f64; 3],
        heat: f64,
        volume: f64,
    ) -> Result<()> {
        let force = cross(current, motion.flux_density_t[cell])?;
        let mechanical = finite(volume * dot(motion.velocity_m_s[cell], force))?;
        let electrical = finite(volume * dot(electric, current))?;
        let residual = finite(finite(heat + mechanical)? - electrical)?;
        self.mechanical_power_w = finite(self.mechanical_power_w + mechanical)?;
        self.electrical_power_w = finite(self.electrical_power_w + electrical)?;
        // Opposite cellwise defects cannot hide in one cancelling global sum.
        self.work_pairing_residual_w = finite(self.work_pairing_residual_w + residual.abs())?;
        self.conductor_field_v_m.push(conductor);
        self.lorentz_force_density_n_m3.push(force);
        self.cell_mechanical_power_w.push(mechanical);
        Ok(())
    }
}

/// Accepted original-cell currents and physical power transfers. Positive
/// mechanical power is delivered TO prescribed motion. Terminal power is
/// positive INTO the conductor and reverses sign during loaded generation.
/// Only Joule power belongs in a thermal receiver's irreversible heat source.
pub struct SpatialMovingConductionField<'a, C = ResolvedConductorState> {
    problem: &'a SpatialConduction<C>,
    velocity_m_s: Vec<[f64; 3]>,
    flux_density_t: Vec<[f64; 3]>,
    field: SpatialConductionField,
    power: MotionPower,
}

impl<C: SpatialConductorMaterial> SpatialMovingConductionField<'_, C> {
    pub fn problem(&self) -> &SpatialConduction<C> {
        self.problem
    }
    /// Consume this accepted result without copying its native field arrays.
    ///
    /// The caller must retain this exact `problem` when coupling the snapshot
    /// to original cells or refreshed material sources. The snapshot has no
    /// topology-dependent load or state-update operations, and cannot be
    /// reattached to an unrelated field owner.
    pub fn into_snapshot(self) -> SpatialMovingConductionSnapshot {
        SpatialMovingConductionSnapshot {
            velocity_m_s: self.velocity_m_s,
            flux_density_t: self.flux_density_t,
            field: self.field,
            power: self.power,
        }
    }
    pub fn voltage_v(&self) -> &[f64] {
        &self.field.voltage_v
    }
    pub fn velocity_m_s(&self) -> &[[f64; 3]] {
        &self.velocity_m_s
    }
    pub fn flux_density_t(&self) -> &[[f64; 3]] {
        &self.flux_density_t
    }
    /// Laboratory `E=-grad(V)` [V/m].
    pub fn electric_field_v_m(&self) -> &[[f64; 3]] {
        &self.field.electric_field_v_m
    }
    /// Conductor-frame `E'=E+u cross B` [V/m].
    pub fn conductor_field_v_m(&self) -> &[[f64; 3]] {
        &self.power.conductor_field_v_m
    }
    pub fn current_density_a_m2(&self) -> &[[f64; 3]] {
        &self.field.current_density_a_m2
    }
    /// Actual tetrahedron-volume integral of `E'.J` [W].
    pub fn cell_joule_power_w(&self) -> &[f64] {
        &self.field.cell_power_w
    }
    /// Total irreversible heat, including stationary contact loss once [W].
    pub fn joule_power_w(&self) -> f64 {
        self.field.power_w
    }
    /// Original P1 contact jumps, currents and surface heat in contact order.
    pub fn contact_response(&self) -> &[SpatialConductionContactResponse] {
        &self.field.contacts
    }
    /// Stationary surface electrical work, entirely irreversible heat [W].
    pub fn contact_power_w(&self) -> f64 {
        self.power.contact_power_w
    }
    pub fn lorentz_force_density_n_m3(&self) -> &[[f64; 3]] {
        &self.power.lorentz_force_density_n_m3
    }
    pub fn cell_mechanical_power_w(&self) -> &[f64] {
        &self.power.cell_mechanical_power_w
    }
    pub fn mechanical_power_w(&self) -> f64 {
        self.power.mechanical_power_w
    }
    /// Signed volume integral of laboratory `E.J` [W].
    /// Stationary contact electrical work is separate in `contact_power_w`.
    pub fn electrical_power_w(&self) -> f64 {
        self.power.electrical_power_w
    }
    /// Independently recovered boundary voltage-current work [W].
    pub fn terminal_power_w(&self) -> f64 {
        self.field.supplied_power_w
    }
    pub fn nodal_boundary_current_a(&self) -> &[Option<f64>] {
        &self.field.nodal_boundary_current_a
    }
    pub fn nodal_applied_current_a(&self) -> &[f64] {
        &self.field.nodal_applied_current_a
    }
    pub fn boundary_current_a(&self) -> f64 {
        self.field.current_residual_a
    }
    pub fn maximum_free_current_residual_a(&self) -> f64 {
        self.field.maximum_free_current_residual_a
    }
    pub fn work_pairing_residual_w(&self) -> f64 {
        self.power.work_pairing_residual_w
    }
    /// Total `Joule + mechanical - terminal` defect [W].
    pub fn power_balance_residual_w(&self) -> f64 {
        self.field.power_residual_w
    }
    /// Bulk plus stationary-contact electrical work minus weak nodal work [W].
    pub fn hodge_power_residual_w(&self) -> f64 {
        self.field.hodge_power_residual_w
    }
    pub fn gauges(&self) -> &[usize] {
        &self.field.gauges
    }
    pub fn solve_report(&self) -> Option<&SolveReport> {
        self.field.linear.as_ref()
    }
}

/// Owned accepted diagnostics for a coupled candidate retaining its original
/// `SpatialConduction` owner. Constructed only by consuming a solved field;
/// storage transfers without allocating or cloning the native arrays.
///
/// Positive mechanical power goes TO the prescribed motion and terminal power
/// goes INTO the conductor. Only conductor-frame Joule power is thermal heat.
/// This diagnostic snapshot does not move geometry or update mechanical state.
#[derive(Clone, Debug)]
pub struct SpatialMovingConductionSnapshot {
    velocity_m_s: Vec<[f64; 3]>,
    flux_density_t: Vec<[f64; 3]>,
    field: SpatialConductionField,
    power: MotionPower,
}

impl SpatialMovingConductionSnapshot {
    pub fn voltage_v(&self) -> &[f64] {
        &self.field.voltage_v
    }
    pub fn velocity_m_s(&self) -> &[[f64; 3]] {
        &self.velocity_m_s
    }
    pub fn flux_density_t(&self) -> &[[f64; 3]] {
        &self.flux_density_t
    }
    /// Laboratory `E=-grad(V)` [V/m].
    pub fn electric_field_v_m(&self) -> &[[f64; 3]] {
        &self.field.electric_field_v_m
    }
    /// Conductor-frame `E'=E+u cross B` [V/m].
    pub fn conductor_field_v_m(&self) -> &[[f64; 3]] {
        &self.power.conductor_field_v_m
    }
    pub fn current_density_a_m2(&self) -> &[[f64; 3]] {
        &self.field.current_density_a_m2
    }
    pub fn cell_joule_power_w(&self) -> &[f64] {
        &self.field.cell_power_w
    }
    /// Total irreversible heat, including stationary contact loss once [W].
    pub fn joule_power_w(&self) -> f64 {
        self.field.power_w
    }
    /// Original P1 contact jumps, currents and surface heat in contact order.
    pub fn contact_response(&self) -> &[SpatialConductionContactResponse] {
        &self.field.contacts
    }
    /// Stationary surface electrical work, entirely irreversible heat [W].
    pub fn contact_power_w(&self) -> f64 {
        self.power.contact_power_w
    }
    pub fn lorentz_force_density_n_m3(&self) -> &[[f64; 3]] {
        &self.power.lorentz_force_density_n_m3
    }
    pub fn cell_mechanical_power_w(&self) -> &[f64] {
        &self.power.cell_mechanical_power_w
    }
    pub fn mechanical_power_w(&self) -> f64 {
        self.power.mechanical_power_w
    }
    pub fn electrical_power_w(&self) -> f64 {
        self.power.electrical_power_w
    }
    pub fn terminal_power_w(&self) -> f64 {
        self.field.supplied_power_w
    }
    pub fn nodal_boundary_current_a(&self) -> &[Option<f64>] {
        &self.field.nodal_boundary_current_a
    }
    pub fn nodal_applied_current_a(&self) -> &[f64] {
        &self.field.nodal_applied_current_a
    }
    pub fn boundary_current_a(&self) -> f64 {
        self.field.current_residual_a
    }
    pub fn maximum_free_current_residual_a(&self) -> f64 {
        self.field.maximum_free_current_residual_a
    }
    pub fn work_pairing_residual_w(&self) -> f64 {
        self.power.work_pairing_residual_w
    }
    pub fn power_balance_residual_w(&self) -> f64 {
        self.field.power_residual_w
    }
    pub fn hodge_power_residual_w(&self) -> f64 {
        self.field.hodge_power_residual_w
    }
    pub fn gauges(&self) -> &[usize] {
        &self.field.gauges
    }
    pub fn solve_report(&self) -> Option<&SolveReport> {
        self.field.linear.as_ref()
    }
}
