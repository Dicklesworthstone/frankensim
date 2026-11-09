//! Charge-conserving conduction in a prescribed moving material snapshot.
//!
//! For in-plane velocity `u` and imposed normal flux density `Bz`, the original
//! scalar or symmetric tensor conductivity acts on `E' = -grad(V) + u x B`.
//! The same prepared conductivity Hodge enforces current continuity and the
//! original voltage/current traces. Unloaded exterior boundaries insulate the
//! TOTAL current, including its motional part.
//!
//! Heat is `E'.J`, force is `J x B`, and signed electrical work is `E.J`:
//! `electrical power = irreversible heat + mechanical power`. These quantities
//! share the same cell current, geometry and material source; a motional EMF
//! cannot become an extra Joule source or an unbalanced terminal current.
//!
//! This is a quasistatic fixed-geometry snapshot with authored cell velocities
//! and imposed normal fields. It does not evolve induction, move the mesh,
//! infer material parameters, add Hall response, or integrate fluid/solid
//! motion. Remesh/re-resolve changed geometry or material through its owner.

use super::{
    ConductionBoundaryResponse, ConductionCellResponse, ConductionContactResponse,
    ConductionControls, ConductionCurrentTrace, ConductionError, ConductionField, PlanarConduction,
    finite, gradients, poll, scalar_potential,
};
use fs_material::conductor::ResolvedConductorState;
use fs_material::conductor::planar::PlanarConductorMaterial;
use fs_solver::SolveReport;

/// Prescribed values on each ORIGINAL triangle, in the mesh-coordinate frame.
/// Positive normal flux density points out of the oriented x-y plane.
#[derive(Clone, Copy, Debug)]
pub struct ConductionMotion<'a> {
    pub velocity_m_s: &'a [[f64; 2]],
    pub normal_flux_density_t: &'a [f64],
}

// Prepared once per frozen conductivity operator and motion snapshot. The
// terminal owner groups this same weak load; trace and electrode recovery use
// exactly the same native current, heat, force and work calculations.
pub(super) struct PreparedMotion<'a> {
    motion: ConductionMotion<'a>,
    field_v_m: Vec<[f64; 2]>,
    pub(super) load_a: Vec<f64>,
    has_field: bool,
}

impl<C: PlanarConductorMaterial> PlanarConduction<C> {
    /// Solve charge continuity with J=Sigma(-grad(V)+u cross B) on the original
    /// SPD conductivity operator. Unloaded boundaries insulate TOTAL current.
    ///
    /// A voltage reference alone lets the free potential screen conservative
    /// motional fields without inventing a closed-circuit braking current.
    /// Both original cells adjacent to a finite electrical contact must remain
    /// stationary; remote moving bulk is admitted. Stationary contact heat has
    /// no mechanical work. Geometry and occupancy remain fixed snapshots.
    pub fn solve_with_motion(
        &self,
        boundary_voltage_v: &[(usize, f64)],
        current_traces: &[ConductionCurrentTrace],
        motion: ConductionMotion<'_>,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingConductionField<'_, C>, ConductionError> {
        poll(cancelled)?;
        self.validate_controls(boundary_voltage_v.len(), controls)?;
        let prescribed = self.prescribed_boundary(boundary_voltage_v, cancelled)?;
        let applied = self.current_loads(current_traces, &prescribed, cancelled)?;
        let prepared = self.prepare_motion(motion, cancelled)?;
        let mut rhs: Vec<_> = applied.iter().map(|load| load.unwrap_or(0.0)).collect();
        if prepared.has_field {
            for (right, &motional) in rhs.iter_mut().zip(&prepared.load_a) {
                *right = finite(*right + motional)?;
            }
        }
        let solved = scalar_potential::solve_relative(
            &self.stiffness,
            &rhs,
            &prescribed,
            controls.linear,
            cancelled,
        )?;
        self.recover_moving_field(
            solved.relative,
            solved.references,
            &prescribed,
            applied,
            &prepared,
            controls,
            solved.report,
            cancelled,
        )
    }

    /// Solve direct ideal voltage/current/floating/mixed electrodes under
    /// prescribed motion with one actual-amplitude grouped sparse solve.
    ///
    /// Total current constrains an electrode's reaction; its local injection
    /// remains a result of the original bulk/contact field. Zero-current
    /// electrodes float and may transfer current between their separate traces.
    /// The native motional weak load is identified together with the stiffness,
    /// never represented by a voltage probe or scaled accepted field.
    ///
    /// The terminal report retains physical voltage drops beneath any common
    /// bias and independently pairs authored current with those solved drops.
    /// Native conductor-frame heat plus signed Lorentz output equals that
    /// signed supply, including generation and zero-current motional voltage.
    /// Extra preparation workspace is 3*vertices + 2*triangles scalar/index
    /// values beyond the stationary terminal bound, admitted before allocation.
    /// The stationary-contact and fixed-planar-geometry restrictions remain.
    pub fn solve_terminals_with_motion(
        &self,
        terminals: &[&[usize]],
        control: super::terminals::ConductionTerminalControl<'_>,
        limits: super::terminals::ConductionTerminalLimits,
        motion: ConductionMotion<'_>,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<
        (
            MovingConductionField<'_, C>,
            super::terminals::ConductionTerminalReport,
        ),
        ConductionError,
    > {
        self.solve_terminals_core(
            terminals,
            control,
            limits,
            Some(motion),
            controls,
            cancelled,
            |relative, references, prescribed, report, prepared, cancelled| {
                self.recover_moving_field(
                    relative,
                    references,
                    prescribed,
                    vec![None; self.mesh.vertices().len()],
                    prepared.ok_or(ConductionError::Invalid(
                        "native terminal motion preparation is absent",
                    ))?,
                    controls,
                    report,
                    &mut || cancelled(),
                )
            },
        )
    }
}

impl<C: super::PlanarOhmicResponse> PlanarConduction<C> {
    pub(super) fn prepare_motion<'m>(
        &self,
        motion: ConductionMotion<'m>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMotion<'m>, ConductionError> {
        poll(cancelled)?;
        if !matches!(self.mesh.metric(), super::Metric2::Planar { .. }) {
            return Err(ConductionError::Invalid(
                "planar motion requires the original planar metric",
            ));
        }
        let cells = self.mesh.faces().len();
        if motion.velocity_m_s.len() != cells || motion.normal_flux_density_t.len() != cells {
            return Err(ConductionError::Invalid(
                "one prescribed velocity and normal magnetic field per original triangle required",
            ));
        }
        if !self.contacts().is_empty() {
            let trace = self.mesh.boundary_trace().map_err(|_| {
                ConductionError::Invalid("electrical contacts require original exterior traces")
            })?;
            for edge in &trace.edges {
                poll(cancelled)?;
                if self.is_contact_edge(edge.oriented_vertices.map(|v| v as usize))
                    && motion.velocity_m_s[edge.source_face] != [0.0; 2]
                {
                    return Err(ConductionError::Invalid(
                        "moving electrical contacts require an owned surface work and heat law",
                    ));
                }
            }
        }
        let mut field_v_m = Vec::with_capacity(cells);
        let mut load_a = vec![0.0; self.mesh.vertices().len()];
        let mut has_field = false;
        for (cell, face) in self.mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            let [ux, uy] = motion.velocity_m_s[cell];
            let b = motion.normal_flux_density_t[cell];
            finite(ux)?;
            finite(uy)?;
            finite(b)?;
            let m = [finite(uy * b)?, finite(-ux * b)?];
            field_v_m.push(m);
            if m == [0.0; 2] {
                continue;
            }
            has_field = true;
            let gradients = gradients(face.map(|v| self.mesh.vertices()[v as usize]))
                .map_err(|_| ConductionError::Invalid("invalid conductor cell gradient"))?;
            let current = self.conductivity[cell].current_density_a_m2(m);
            let j = [finite(current[0])?, finite(current[1])?];
            let volume = self
                .mesh
                .face_measure(cell)
                .ok_or(ConductionError::Invalid("missing conductor cell volume"))?;
            let mut load = [0.0; 3];
            for local in 1..3 {
                load[local] =
                    finite(volume * gradients[local][0].mul_add(j[0], gradients[local][1] * j[1]))?;
            }
            // Partition of unity makes the cell source sum zero. Use the same
            // reference basis as the gauge-invariant electric-field recovery.
            load[0] = finite(-finite(load[1] + load[2])?)?;
            for local in 0..3 {
                let vertex = face[local] as usize;
                load_a[vertex] = finite(load_a[vertex] + load[local])?;
            }
        }
        Ok(PreparedMotion {
            motion,
            field_v_m,
            load_a,
            has_field,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn recover_moving_field(
        &self,
        relative: Vec<f64>,
        references: Vec<f64>,
        prescribed: &[Option<f64>],
        applied: Vec<Option<f64>>,
        prepared: &PreparedMotion<'_>,
        controls: ConductionControls,
        report: Option<SolveReport>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingConductionField<'_, C>, ConductionError> {
        let motion = prepared.motion;
        let cells = self.mesh.faces().len();
        let motion_field = prepared.has_field.then_some(prepared.field_v_m.as_slice());
        let motion_load = prepared.has_field.then_some(prepared.load_a.as_slice());
        let response = self.cell_response_with_motion(&relative, motion_field, cancelled)?;
        let contact_response = self.contact_response(&relative, cancelled)?;
        let mut contact_power_w = 0.0;
        for contact in &contact_response {
            poll(cancelled)?;
            contact_power_w = finite(contact_power_w + contact.power_w)?;
        }
        let joule_power_w = finite(response.power_w + contact_power_w)?;
        let mut conductor_field_v_m = Vec::with_capacity(cells);
        let mut lorentz_force_density_n_m3 = Vec::with_capacity(cells);
        let mut cell_mechanical_power_w = Vec::with_capacity(cells);
        let mut mechanical_power_w = 0.0;
        let mut electrical_power_w = 0.0;
        let mut work_pairing_residual_w = 0.0;
        for (cell, &m) in prepared.field_v_m.iter().enumerate() {
            poll(cancelled)?;
            let e = response.electric_field_v_m[cell];
            let j = response.current_density_a_m2[cell];
            let [ux, uy] = motion.velocity_m_s[cell];
            let b = motion.normal_flux_density_t[cell];
            let force = [finite(j[1] * b)?, finite(-j[0] * b)?];
            let volume = self
                .mesh
                .face_measure(cell)
                .ok_or(ConductionError::Invalid("missing conductor cell volume"))?;
            let mechanical = finite(volume * ux.mul_add(force[0], uy * force[1]))?;
            let electrical = finite(volume * e[0].mul_add(j[0], e[1] * j[1]))?;
            let residual = finite(finite(response.cell_power_w[cell] + mechanical)? - electrical)?;
            // An L1 check prevents opposite local work errors from cancelling.
            work_pairing_residual_w = finite(work_pairing_residual_w + residual.abs())?;
            mechanical_power_w = finite(mechanical_power_w + mechanical)?;
            electrical_power_w = finite(electrical_power_w + electrical)?;
            conductor_field_v_m.push([finite(e[0] + m[0])?, finite(e[1] + m[1])?]);
            lorentz_force_density_n_m3.push(force);
            cell_mechanical_power_w.push(mechanical);
        }
        let boundary = self.boundary_response(
            &relative,
            prescribed,
            &applied,
            motion_load,
            finite(electrical_power_w + contact_power_w)?,
            controls.current_tolerance_a,
            controls.power_tolerance_w,
            cancelled,
        )?;
        let power_balance_residual_w =
            finite(finite(joule_power_w + mechanical_power_w)? - boundary.terminal_power_w)?;
        if work_pairing_residual_w > controls.power_tolerance_w
            || power_balance_residual_w.abs() > controls.power_tolerance_w
        {
            return Err(ConductionError::Balance {
                current_residual_a: boundary.boundary_current_a,
                maximum_free_current_residual_a: boundary.maximum_free_current_residual_a,
                power_residual_w: if work_pairing_residual_w > controls.power_tolerance_w {
                    work_pairing_residual_w
                } else {
                    power_balance_residual_w
                },
            });
        }
        let mut voltage = Vec::with_capacity(relative.len());
        for v in 0..relative.len() {
            poll(cancelled)?;
            voltage.push(match prescribed[v] {
                Some(value) => value,
                None => finite(references[v] + relative[v])?,
            });
        }
        poll(cancelled)?;
        Ok(MovingConductionField {
            problem: self,
            voltage_v: voltage,
            relative_voltage_v: relative,
            reference_voltage_v: references,
            velocity_m_s: motion.velocity_m_s.to_vec(),
            normal_flux_density_t: motion.normal_flux_density_t.to_vec(),
            response,
            contact_response,
            contact_power_w,
            joule_power_w,
            conductor_field_v_m,
            lorentz_force_density_n_m3,
            cell_mechanical_power_w,
            mechanical_power_w,
            electrical_power_w,
            boundary,
            nodal_applied_current_a: applied,
            work_pairing_residual_w,
            power_balance_residual_w,
            report,
        })
    }
}

/// Original prescribed-motion observations accompanying a nonlinear field.
/// The conductor-frame field drives material response and irreversible heat.
/// Mechanical power is positive OUT to motion; terminal power is positive IN.
#[derive(Clone, Debug)]
pub struct ConductionMotionReport {
    pub velocity_m_s: Vec<[f64; 2]>,
    pub normal_flux_density_t: Vec<f64>,
    pub conductor_field_v_m: Vec<[f64; 2]>,
    pub lorentz_force_density_n_m3: Vec<[f64; 2]>,
    pub cell_mechanical_power_w: Vec<f64>,
    pub mechanical_power_w: f64,
    /// Signed bulk laboratory E.J, excluding stationary contact work.
    pub electrical_power_w: f64,
    pub contact_power_w: f64,
    pub terminal_power_w: f64,
    /// Sum of absolute local heat + mechanical - electrical defects.
    pub work_pairing_residual_w: f64,
}

/// Accepted current, heat and force with the exact original source/geometry
/// and prescribed velocity/field snapshots. Mechanical power is positive when
/// the electromagnetic force does work on the moving conductor; terminal power
/// is positive into the conductor and can be negative during generation.
pub struct MovingConductionField<'a, C = ResolvedConductorState> {
    problem: &'a PlanarConduction<C>,
    voltage_v: Vec<f64>,
    relative_voltage_v: Vec<f64>,
    reference_voltage_v: Vec<f64>,
    velocity_m_s: Vec<[f64; 2]>,
    normal_flux_density_t: Vec<f64>,
    response: ConductionCellResponse,
    contact_response: Vec<ConductionContactResponse>,
    contact_power_w: f64,
    joule_power_w: f64,
    conductor_field_v_m: Vec<[f64; 2]>,
    lorentz_force_density_n_m3: Vec<[f64; 2]>,
    cell_mechanical_power_w: Vec<f64>,
    mechanical_power_w: f64,
    electrical_power_w: f64,
    boundary: ConductionBoundaryResponse,
    nodal_applied_current_a: Vec<Option<f64>>,
    work_pairing_residual_w: f64,
    power_balance_residual_w: f64,
    report: Option<SolveReport>,
}

impl<'a, C: PlanarConductorMaterial> MovingConductionField<'a, C> {
    // Consume native arrays into the shared nonlinear source closure. This
    // projection never escapes as a stationary public solve: its caller pairs
    // the field with the mandatory physical motion account below.
    pub(super) fn into_response_parts(self) -> (ConductionField<'a, C>, ConductionMotionReport) {
        let motion = ConductionMotionReport {
            velocity_m_s: self.velocity_m_s,
            normal_flux_density_t: self.normal_flux_density_t,
            conductor_field_v_m: self.conductor_field_v_m,
            lorentz_force_density_n_m3: self.lorentz_force_density_n_m3,
            cell_mechanical_power_w: self.cell_mechanical_power_w,
            mechanical_power_w: self.mechanical_power_w,
            electrical_power_w: self.electrical_power_w,
            contact_power_w: self.contact_power_w,
            terminal_power_w: self.boundary.terminal_power_w,
            work_pairing_residual_w: self.work_pairing_residual_w,
        };
        let field = ConductionField {
            problem: self.problem,
            voltage_v: self.voltage_v,
            relative_voltage_v: self.relative_voltage_v,
            reference_voltage_v: self.reference_voltage_v,
            electric_field_v_m: self.response.electric_field_v_m,
            current_density_a_m2: self.response.current_density_a_m2,
            cell_power_w: self.response.cell_power_w,
            contact_response: self.contact_response,
            power_w: self.joule_power_w,
            nodal_boundary_current_a: self.boundary.nodal_boundary_current_a,
            nodal_applied_current_a: self.nodal_applied_current_a,
            boundary_current_a: self.boundary.boundary_current_a,
            maximum_free_current_residual_a: self.boundary.maximum_free_current_residual_a,
            power_balance_residual_w: self.power_balance_residual_w,
            report: self.report,
        };
        (field, motion)
    }

    pub fn problem(&self) -> &PlanarConduction<C> {
        self.problem
    }
    /// Absolute reporting view, with authored voltage traces preserved exactly.
    /// Current, heat and force use the retained relative solution, since a
    /// large common voltage can round away a small motional or resistive drop.
    pub fn voltage_v(&self) -> &[f64] {
        &self.voltage_v
    }
    pub fn relative_voltage_v(&self) -> &[f64] {
        &self.relative_voltage_v
    }
    pub fn reference_voltage_v(&self) -> &[f64] {
        &self.reference_voltage_v
    }
    pub fn velocity_m_s(&self) -> &[[f64; 2]] {
        &self.velocity_m_s
    }
    pub fn normal_flux_density_t(&self) -> &[f64] {
        &self.normal_flux_density_t
    }
    /// Laboratory electric field `E = -grad(V)` on each original cell [V/m].
    pub fn electric_field_v_m(&self) -> &[[f64; 2]] {
        &self.response.electric_field_v_m
    }
    /// Conductor-frame field `E' = E + u x B` used by the conductivity [V/m].
    pub fn conductor_field_v_m(&self) -> &[[f64; 2]] {
        &self.conductor_field_v_m
    }
    pub fn current_density_a_m2(&self) -> &[[f64; 2]] {
        &self.response.current_density_a_m2
    }
    /// Irreversible `E'.J` integrated over each actual extruded cell [W].
    pub fn cell_joule_power_w(&self) -> &[f64] {
        &self.response.cell_power_w
    }
    /// Total irreversible heat, including stationary contact loss once [W].
    pub fn joule_power_w(&self) -> f64 {
        self.joule_power_w
    }
    /// Original P1 contact jumps, currents and surface heat in contact order.
    pub fn contact_response(&self) -> &[ConductionContactResponse] {
        &self.contact_response
    }
    /// Stationary surface electrical work, entirely irreversible heat [W].
    pub fn contact_power_w(&self) -> f64 {
        self.contact_power_w
    }
    /// Original cell `J x B` [N/m3], without interfacial smoothing.
    pub fn lorentz_force_density_n_m3(&self) -> &[[f64; 2]] {
        &self.lorentz_force_density_n_m3
    }
    pub fn cell_mechanical_power_w(&self) -> &[f64] {
        &self.cell_mechanical_power_w
    }
    pub fn mechanical_power_w(&self) -> f64 {
        self.mechanical_power_w
    }
    /// Signed volume integral of the LABORATORY `E.J` [W].
    /// Stationary contact electrical work is separate in `contact_power_w`.
    pub fn electrical_power_w(&self) -> f64 {
        self.electrical_power_w
    }
    /// Independently recovered boundary voltage-current work, positive in [W].
    pub fn terminal_power_w(&self) -> f64 {
        self.boundary.terminal_power_w
    }
    /// Combined voltage reactions and authored inward current loads [A].
    pub fn nodal_boundary_current_a(&self) -> &[Option<f64>] {
        &self.boundary.nodal_boundary_current_a
    }
    pub fn nodal_applied_current_a(&self) -> &[Option<f64>] {
        &self.nodal_applied_current_a
    }
    pub fn boundary_current_a(&self) -> f64 {
        self.boundary.boundary_current_a
    }
    pub fn maximum_free_current_residual_a(&self) -> f64 {
        self.boundary.maximum_free_current_residual_a
    }
    /// Sum of absolute local `heat + mechanical - electrical` defects [W].
    pub fn work_pairing_residual_w(&self) -> f64 {
        self.work_pairing_residual_w
    }
    /// Total `heat + mechanical - terminal` defect [W].
    pub fn power_balance_residual_w(&self) -> f64 {
        self.power_balance_residual_w
    }
    pub fn solve_report(&self) -> Option<&SolveReport> {
        self.report.as_ref()
    }
}

/// Exact affine two-electrode response of ONE frozen conductivity and motion
/// snapshot. Preparation uses one homogeneous resistance solve and one native
/// zero-voltage moving solve. Actual-drive recovery superposes their potentials
/// and re-evaluates original cell/contact current, heat and Lorentz work.
///
/// This response is invalid after changing material, geometry, motion or B.
/// A nonlinear owner must prepare it anew on every original constitutive trial;
/// it is never an accepted nonlinear resistance or a scalable heat field.
pub struct PlanarMovingResistance<'p, 'm, C = ResolvedConductorState> {
    resistance: super::PlanarResistance<'p, C>,
    prepared: PreparedMotion<'m>,
    zero_relative_voltage_v: Vec<f64>,
    motion_current_a: f64,
    controls: ConductionControls,
}

impl<C: PlanarConductorMaterial> PlanarConduction<C> {
    /// Prepare I=V/R+I_motion on the original frozen material operator.
    /// Positive current enters the positive electrode. Both electrodes are
    /// complete disjoint exterior traces, with the original native resistance
    /// and stationary-contact admission. Every recovered actual field retains
    /// source validity, charge/current and local/global work gates.
    ///
    /// Two sparse field solves are required; later actual-voltage recovery
    /// performs no sparse solve. Their preparation buffers remain bounded by
    /// this already admitted original mesh and the caller's field/work limits.
    pub fn moving_resistance<'p, 'm>(
        &'p self,
        positive: &[usize],
        negative: &[usize],
        motion: ConductionMotion<'m>,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarMovingResistance<'p, 'm, C>, ConductionError> {
        poll(cancelled)?;
        let prepared = self.prepare_motion(motion, cancelled)?;
        let resistance = self.resistance(positive, negative, controls, cancelled)?;
        let mut prescribed = vec![None; self.mesh.vertices().len()];
        for &vertex in positive.iter().chain(negative) {
            poll(cancelled)?;
            prescribed[vertex] = Some(0.0);
        }
        let solved = scalar_potential::solve_relative(
            &self.stiffness,
            &prepared.load_a,
            &prescribed,
            controls.linear,
            cancelled,
        )?;
        let zero = self.recover_moving_field(
            solved.relative,
            solved.references,
            &prescribed,
            vec![None; self.mesh.vertices().len()],
            &prepared,
            controls,
            solved.report,
            cancelled,
        )?;
        let mut current = fs_math::dd::Dd::from_f64(0.0);
        for &vertex in positive {
            poll(cancelled)?;
            current = current
                + fs_math::dd::Dd::from_f64(zero.nodal_boundary_current_a()[vertex].ok_or(
                    ConductionError::Invalid("moving affine electrode lost its reaction"),
                )?);
        }
        let motion_current_a = finite(current.to_f64())?;
        poll(cancelled)?;
        Ok(PlanarMovingResistance {
            resistance,
            prepared,
            zero_relative_voltage_v: zero.relative_voltage_v,
            motion_current_a,
            controls,
        })
    }
}

impl<'p, C: PlanarConductorMaterial> PlanarMovingResistance<'p, '_, C> {
    /// Homogeneous resistance of the SAME frozen bulk/contact operator [ohm].
    pub fn resistance_ohm(&self) -> f64 {
        self.resistance.resistance_ohm()
    }
    /// Signed inward positive-electrode current at zero voltage [A].
    pub fn motion_current_a(&self) -> f64 {
        self.motion_current_a
    }

    /// Recover the original moving field at an actual voltage [V], without
    /// another sparse solve. Only native potentials are superposed; current,
    /// source-field admission, contact loss, heat and mechanical work are all
    /// freshly recovered. No accepted nonlinear field or heat is scaled.
    pub fn field_at_voltage(
        &self,
        voltage_v: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingConductionField<'p, C>, ConductionError> {
        use fs_math::dd::Dd;
        poll(cancelled)?;
        finite(voltage_v)?;
        let basis = self.resistance.reference_field();
        let reference = self.resistance.reference_voltage_v();
        let mut relative = Vec::with_capacity(self.zero_relative_voltage_v.len());
        let mut references = Vec::with_capacity(self.zero_relative_voltage_v.len());
        for vertex in 0..self.zero_relative_voltage_v.len() {
            poll(cancelled)?;
            relative.push(finite(
                (Dd::from_f64(self.zero_relative_voltage_v[vertex])
                    + Dd::from_f64(basis.relative_voltage_v()[vertex] / reference)
                        * Dd::from_f64(voltage_v))
                .to_f64(),
            )?);
            references.push(finite(
                (basis.reference_voltage_v()[vertex] / reference) * voltage_v,
            )?);
        }
        let problem = self.resistance.field.problem;
        let mut prescribed = vec![None; relative.len()];
        for (vertices, value) in self.resistance.terminals.iter().zip([voltage_v, 0.0]) {
            for &vertex in vertices {
                poll(cancelled)?;
                prescribed[vertex] = Some(value);
            }
        }
        let field = problem.recover_moving_field(
            relative,
            references,
            &prescribed,
            vec![None; self.zero_relative_voltage_v.len()],
            &self.prepared,
            self.controls,
            None,
            cancelled,
        )?;
        let mut actual = Dd::from_f64(0.0);
        for &vertex in &self.resistance.terminals[0] {
            poll(cancelled)?;
            actual = actual
                + Dd::from_f64(field.nodal_boundary_current_a()[vertex].ok_or(
                    ConductionError::Invalid("moving affine electrode lost its actual reaction"),
                )?);
        }
        let expected =
            Dd::from_f64(voltage_v / self.resistance_ohm()) + Dd::from_f64(self.motion_current_a);
        let current_defect = finite((actual - expected).to_f64())?;
        let power_defect = finite(
            (Dd::from_f64(field.terminal_power_w()) - Dd::from_f64(voltage_v) * expected).to_f64(),
        )?;
        if current_defect.abs() > self.controls.current_tolerance_a
            || power_defect.abs() > self.controls.power_tolerance_w
        {
            return Err(ConductionError::Balance {
                current_residual_a: current_defect,
                maximum_free_current_residual_a: field.maximum_free_current_residual_a(),
                power_residual_w: power_defect,
            });
        }
        poll(cancelled)?;
        Ok(field)
    }
}
