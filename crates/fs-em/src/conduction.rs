//! Spatial steady Ohmic conduction from complete material source snapshots.
//!
//! The conductivity Hodge solves -div(sigma grad(V))=0 on the original planar
//! mesh. E=-grad(V), J=sigma E and the actual volume integral E.J determine
//! current and irreversible heating. The dielectric field shares only scalar
//! algebra and electrode geometry; its charge/storage law is not reused here.
//!
//! This is a stationary symmetric Ohmic model with fixed geometry and material
//! state. Induction, displacement current, contacts, Hall response and nonlinear
//! field-dependent conduction require their own models. A changed temperature
//! or material must be resolved and prepared again, not substituted into an
//! already accepted field. No electrical storage or thermal evolution is added.
//! The separate [`crate::electrostatic::relaxation`] owner composes these same
//! cell currents/losses with retained dielectric charge and displacement current.

use fs_feec::{
    assembly::{incidence_to_csr, stiffness},
    weighted::{WeightedAssemblyLimits, WeightedError},
    weighted2::{CellWeight2, weighted_mass_matrix_2d},
};
use fs_material::conductor::planar::{PlanarConductivity, PlanarConductorMaterial};
use fs_material::conductor::{ELECTRICAL_CONDUCTIVITY_DIMS, ResolvedConductorState};
use fs_phs::circuit::Branch;
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::{CheckedCgConfig, SolveReport};
use fs_sparse::Csr;
use std::collections::BTreeMap;

use crate::planar::gradients;
use crate::scalar_potential::{self, PotentialError};

/// Independent algebraic, current-conservation and power acceptance limits.
#[derive(Clone, Copy, Debug)]
pub struct ConductionControls {
    pub linear: CheckedCgConfig,
    /// Absolute net terminal current and free-node current residual [A].
    pub current_tolerance_a: f64,
    /// Absolute difference between integrated Joule and terminal power [W].
    pub power_tolerance_w: f64,
}

/// Uniform inward current density on one original exterior edge. The edge's
/// extruded area and P1 trace basis supply its nodal current loads; this is not
/// a point electrode or an equipotential terminal. Vertex order is immaterial.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductionCurrentTrace {
    pub vertices: [usize; 2],
    /// Positive INTO the conductor [A/m2], negative for extraction.
    pub inward_current_density_a_m2: f64,
}

/// Material-model, geometry, numerical or physical-balance refusal.
#[derive(Debug)]
pub enum ConductionError {
    Assembly(WeightedError),
    Invalid(&'static str),
    Cancelled,
    Nonconvergence(SolveReport),
    Balance {
        current_residual_a: f64,
        maximum_free_current_residual_a: f64,
        power_residual_w: f64,
    },
}
impl core::fmt::Display for ConductionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar Ohmic conduction: {self:?}")
    }
}
impl std::error::Error for ConductionError {}
impl From<PotentialError> for ConductionError {
    fn from(error: PotentialError) -> Self {
        match error {
            PotentialError::Invalid(reason) => Self::Invalid(reason),
            PotentialError::Cancelled => Self::Cancelled,
            PotentialError::Nonconvergence(report) => Self::Nonconvergence(report),
        }
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ConductionError> {
    if cancelled() {
        Err(ConductionError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ConductionError::Invalid(
            "conductor field arithmetic is not representable",
        ))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

/// Immutable original mesh, source-bound conductors and conductivity
/// Hodge. One complete material snapshot is retained per original triangle.
/// The default scalar type preserves existing callers; oriented/mixed source
/// states use this SAME field operator, boundary solve, current and heat owner.
pub struct PlanarConduction<C = ResolvedConductorState> {
    mesh: TriComplex2,
    conductors: Vec<C>,
    conductivity: Vec<PlanarConductivity>,
    conductivity_s_m: Vec<f64>,
    stiffness: Csr,
}

// One original Ohmic cell law for steady balance and retained charge dynamics.
pub(crate) struct ConductionCellResponse {
    pub electric_field_v_m: Vec<[f64; 2]>,
    pub current_density_a_m2: Vec<[f64; 2]>,
    pub cell_power_w: Vec<f64>,
    pub power_w: f64,
}
impl<C: PlanarConductorMaterial> PlanarConduction<C> {
    /// Assemble the same native weighted FEEC operator used by other scalar
    /// fields, with the conductivity's actual SI dimensions. The canonical
    /// conductor owner selects resistivity or conductivity explicitly; material
    /// names and an assumed unit wire play no role in the spatial law.
    pub fn new(
        mesh: TriComplex2,
        conductors: Vec<C>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConductionError> {
        poll(cancelled)?;
        if !matches!(mesh.metric(), Metric2::Planar { .. }) {
            return Err(ConductionError::Invalid(
                "this conductor owner requires a planar metric",
            ));
        }
        if conductors.len() != mesh.faces().len()
            || mesh.vertices().len() > limits.max_dofs
            || mesh.edges().len() > limits.max_dofs
            || mesh.faces().len() > limits.max_cells
            || mesh
                .faces()
                .len()
                .checked_mul(9)
                .is_none_or(|n| n > limits.max_triplets)
        {
            return Err(ConductionError::Invalid(
                "one conductor per triangle and bounded assembly extents required",
            ));
        }
        let mut conductivity_s_m = Vec::with_capacity(conductors.len());
        let mut conductivity = Vec::with_capacity(conductors.len());
        let mut weights = Vec::with_capacity(conductors.len());
        let mut geometry_frame = None;
        for conductor in &conductors {
            poll(cancelled)?;
            if !conductor.has_admitted_static_sources() {
                return Err(ConductionError::Invalid(
                    "spatial DC conduction requires static Ohmic coefficients at fixed external or explicitly owned phase conditions",
                ));
            }
            if let Some(frame) = conductor.frame() {
                if geometry_frame.is_some_and(|id| id != frame.geometry_frame()) {
                    return Err(ConductionError::Invalid(
                        "all oriented conductors must use one mesh-coordinate frame",
                    ));
                }
                geometry_frame = Some(frame.geometry_frame());
            }
            let sigma = conductor.conductivity();
            let weight = match sigma {
                PlanarConductivity::Scalar(value) => {
                    if !positive(value) {
                        return Err(ConductionError::Invalid(
                            "source conductivity must be positive and representable",
                        ));
                    }
                    conductivity_s_m.push(value);
                    CellWeight2::Scalar(value)
                }
                // Existing outward-interval SPD admission belongs to the
                // original weighted FEEC owner, not a second electrical gate.
                PlanarConductivity::Tensor(matrix) => CellWeight2::Tensor(matrix),
            };
            conductivity.push(sigma);
            weights.push(weight);
        }
        let mass = weighted_mass_matrix_2d(
            &mesh,
            1,
            &weights,
            ELECTRICAL_CONDUCTIVITY_DIMS,
            limits,
            cancelled,
        )
        .map_err(|error| match error {
            WeightedError::Cancelled => ConductionError::Cancelled,
            error => ConductionError::Assembly(error),
        })?;
        poll(cancelled)?;
        let stiffness = stiffness(&incidence_to_csr(&mesh.d0()), mass.matrix());
        for row in 0..stiffness.nrows() {
            poll(cancelled)?;
            if stiffness.row(row).1.iter().any(|v| !v.is_finite()) {
                return Err(ConductionError::Invalid(
                    "unrepresentable conductor stiffness",
                ));
            }
        }
        Ok(Self {
            mesh,
            conductors,
            conductivity,
            conductivity_s_m,
            stiffness,
        })
    }

    pub fn mesh(&self) -> &TriComplex2 {
        &self.mesh
    }
    /// Complete original card, query, requirements and electrical source basis.
    pub fn conductors(&self) -> &[C] {
        &self.conductors
    }
    /// Actual scalar/tensor conductivity in the original mesh coordinates.
    pub fn conductivity_coefficients(&self) -> &[PlanarConductivity] {
        &self.conductivity
    }

    pub(crate) fn stiffness_operator(&self) -> &Csr {
        &self.stiffness
    }

    /// Solve at prescribed exterior voltage traces. Remaining exterior edges
    /// are electrically insulated; each connected component needs a complete
    /// anchored trace. No volume-current source or point electrode is inferred.
    pub fn solve(
        &self,
        boundary_voltage_v: &[(usize, f64)],
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionField<'_, C>, ConductionError> {
        self.solve_with_current_traces(boundary_voltage_v, &[], controls, cancelled)
    }

    /// Solve mixed voltage and prescribed inward-current-density boundaries.
    /// Each current record covers an actual complete exterior edge; its weak
    /// load is J_in * edge_area / 2 at each endpoint. Remaining edges insulate.
    /// Each connected component still needs a voltage trace to fix its gauge
    /// and carry the return current. A current edge cannot also have both
    /// endpoint voltages prescribed. Adjacent voltage/current edges may meet
    /// at a vertex; the voltage-source reaction excludes the applied load there.
    ///
    /// Current density is prescribed, not voltage or equipotentiality on that
    /// trace. This supports distributed current injection without reducing a
    /// heterogeneous specimen to one resistance. The original FEEC solver and
    /// cell E.J recovery own the actual spatial response and power acceptance.
    pub fn solve_with_current_traces(
        &self,
        boundary_voltage_v: &[(usize, f64)],
        current_traces: &[ConductionCurrentTrace],
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionField<'_, C>, ConductionError> {
        poll(cancelled)?;
        if boundary_voltage_v.len() > self.mesh.vertices().len()
            || !positive(controls.linear.tolerance)
            || controls.linear.tolerance >= 1.0
            || !positive(controls.current_tolerance_a)
            || !positive(controls.power_tolerance_w)
        {
            return Err(ConductionError::Invalid(
                "bounded boundary traces and positive solve/current/power tolerances required",
            ));
        }
        let prescribed =
            scalar_potential::prescribed_boundary(&self.mesh, boundary_voltage_v, cancelled)?;
        let applied = self.current_loads(current_traces, &prescribed, cancelled)?;
        let rhs: Vec<_> = applied.iter().map(|load| load.unwrap_or(0.0)).collect();
        let (voltage, report) = scalar_potential::solve(
            &self.stiffness,
            &rhs,
            &prescribed,
            controls.linear,
            cancelled,
        )?;
        self.recover_field(
            voltage,
            &prescribed,
            applied,
            controls.current_tolerance_a,
            controls.power_tolerance_w,
            report,
            cancelled,
        )
    }

    fn current_loads(
        &self,
        currents: &[ConductionCurrentTrace],
        prescribed: &[Option<f64>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<Option<f64>>, ConductionError> {
        let mut loads = vec![None; self.mesh.vertices().len()];
        if currents.is_empty() {
            return Ok(loads);
        }
        let trace = self
            .mesh
            .boundary_trace()
            .map_err(|_| ConductionError::Invalid("invalid current boundary trace"))?;
        if currents.len() > trace.edges.len() {
            return Err(ConductionError::Invalid(
                "current traces exceed exterior edge count",
            ));
        }
        let mut exterior = BTreeMap::new();
        for edge in &trace.edges {
            poll(cancelled)?;
            let mut vertices = edge.oriented_vertices.map(|v| v as usize);
            vertices.sort_unstable();
            exterior.insert(vertices, edge.global_edge);
        }
        let mut ordered = BTreeMap::new();
        for current in currents {
            poll(cancelled)?;
            let mut vertices = current.vertices;
            vertices.sort_unstable();
            if !current.inward_current_density_a_m2.is_finite()
                || !exterior.contains_key(&vertices)
                || ordered
                    .insert(vertices, current.inward_current_density_a_m2)
                    .is_some()
                || (prescribed[vertices[0]].is_some() && prescribed[vertices[1]].is_some())
            {
                return Err(ConductionError::Invalid(
                    "current traces need distinct exterior edges without a prescribed voltage trace",
                ));
            }
        }
        // Canonical edge order makes input permutation and orientation immaterial.
        for (vertices, density) in ordered {
            poll(cancelled)?;
            let area =
                self.mesh
                    .edge_measure(exterior[&vertices])
                    .ok_or(ConductionError::Invalid(
                        "current trace lacks physical area",
                    ))?;
            let half = finite(finite(area * density)? * 0.5)?;
            if density != 0.0 && half == 0.0 {
                return Err(ConductionError::Invalid(
                    "applied edge current is unrepresentable",
                ));
            }
            for vertex in vertices {
                let previous = loads[vertex].unwrap_or(0.0);
                let sum = finite(previous + half)?;
                if half != 0.0 && sum == previous {
                    return Err(ConductionError::Invalid(
                        "applied edge current disappeared in nodal sum",
                    ));
                }
                loads[vertex] = Some(sum);
            }
        }
        Ok(loads)
    }

    pub(crate) fn cell_response(
        &self,
        voltage: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionCellResponse, ConductionError> {
        if voltage.len() != self.mesh.vertices().len() || voltage.iter().any(|v| !v.is_finite()) {
            return Err(ConductionError::Invalid("complete finite native potentials required"));
        }
        let mut electric_field_v_m = Vec::with_capacity(self.conductors.len());
        let mut current_density_a_m2 = Vec::with_capacity(self.conductors.len());
        let mut cell_power_w = Vec::with_capacity(self.conductors.len());
        let mut power_w = 0.0;
        for (cell, face) in self.mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            let gradients = gradients(face.map(|v| self.mesh.vertices()[v as usize]))
                .map_err(|_| ConductionError::Invalid("invalid conductor cell gradient"))?;
            let mut e = [0.0; 2];
            for local in 1..3 {
                let delta = finite(voltage[face[local] as usize] - voltage[face[0] as usize])?;
                for axis in 0..2 {
                    e[axis] = finite((-delta).mul_add(gradients[local][axis], e[axis]))?;
                }
            }
            let current = self.conductivity[cell].current_density_a_m2(e);
            let j = [finite(current[0])?, finite(current[1])?];
            let volume = self
                .mesh
                .face_measure(cell)
                .ok_or(ConductionError::Invalid("missing conductor cell volume"))?;
            // This is irreversible POWER, with no dielectric-energy factor 1/2.
            let power = finite(volume * e[0].mul_add(j[0], e[1] * j[1]))?;
            if power < 0.0 || (e != [0.0; 2] && power == 0.0) {
                return Err(ConductionError::Invalid(
                    "conductor Joule power is not representably positive",
                ));
            }
            power_w = finite(power_w + power)?;
            electric_field_v_m.push(e);
            current_density_a_m2.push(j);
            cell_power_w.push(power);
        }
        Ok(ConductionCellResponse {
            electric_field_v_m, current_density_a_m2, cell_power_w, power_w,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn recover_field(
        &self,
        voltage: Vec<f64>,
        prescribed: &[Option<f64>],
        applied: Vec<Option<f64>>,
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        report: Option<SolveReport>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionField<'_, C>, ConductionError> {
        let n = self.mesh.vertices().len();
        let ConductionCellResponse {
            electric_field_v_m, current_density_a_m2, cell_power_w, power_w,
        } = self.cell_response(&voltage, cancelled)?;
        let mut nodal_boundary_current_a = vec![None; n];
        let mut boundary_current_a = 0.0;
        let mut maximum_free_current_residual_a: f64 = 0.0;
        let mut terminal_power_w = 0.0;
        // A common voltage cannot change terminal work. Form local differences
        // in the weak reaction as well as the field before multiplication.
        let reference = prescribed.iter().flatten().copied().next().unwrap_or(0.0);
        for v in 0..n {
            poll(cancelled)?;
            let (cols, vals) = self.stiffness.row(v);
            let mut reaction = 0.0;
            for (&w, &value) in cols.iter().zip(vals) {
                let difference = finite(voltage[w] - voltage[v])?;
                reaction = finite(value.mul_add(difference, reaction))?;
            }
            let injection = applied[v].unwrap_or(0.0);
            if prescribed[v].is_none() {
                maximum_free_current_residual_a =
                    maximum_free_current_residual_a.max(finite(reaction - injection)?.abs());
            }
            // On a voltage vertex the source supplies reaction - injection;
            // adding the authored injection gives the actual combined current.
            // On a current-only vertex use the AUTHORED load, independently of
            // the recovered field reaction, for conservation and source work.
            let terminal = if prescribed[v].is_some() {
                Some(reaction)
            } else {
                applied[v]
            };
            if let Some(current) = terminal {
                nodal_boundary_current_a[v] = Some(current);
                boundary_current_a = finite(boundary_current_a + current)?;
                let difference = finite(voltage[v] - reference)?;
                terminal_power_w = finite(difference.mul_add(current, terminal_power_w))?;
            }
        }
        let power_balance_residual_w = finite(power_w - terminal_power_w)?;
        if boundary_current_a.abs() > current_tolerance_a
            || maximum_free_current_residual_a > current_tolerance_a
            || power_balance_residual_w.abs() > power_tolerance_w
        {
            return Err(ConductionError::Balance {
                current_residual_a: boundary_current_a,
                maximum_free_current_residual_a,
                power_residual_w: power_balance_residual_w,
            });
        }
        poll(cancelled)?;
        Ok(ConductionField {
            problem: self,
            voltage_v: voltage,
            electric_field_v_m,
            current_density_a_m2,
            cell_power_w,
            power_w,
            nodal_boundary_current_a,
            nodal_applied_current_a: applied,
            boundary_current_a,
            maximum_free_current_residual_a,
            power_balance_residual_w,
            report,
        })
    }

    /// Extract the original geometry's two-terminal resistance. One volt drives
    /// disjoint complete electrode traces; other exterior edges are insulated.
    /// Integrated E.J must agree with independently summed terminal current.
    pub fn resistance(
        &self,
        positive_terminal: &[usize],
        negative_terminal: &[usize],
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarResistance<'_, C>, ConductionError> {
        let excitation_v = 1.0;
        let boundary = scalar_potential::terminal_boundary(
            &self.mesh,
            &[positive_terminal, negative_terminal],
            &[excitation_v, 0.0],
            cancelled,
        )?;
        let field = self.solve(&boundary, controls, cancelled)?;
        let mut terminal_current_a = 0.0;
        for &v in positive_terminal {
            poll(cancelled)?;
            terminal_current_a = finite(
                terminal_current_a
                    + field.nodal_boundary_current_a[v].ok_or(ConductionError::Invalid(
                        "lost original conductor terminal reaction",
                    ))?,
            )?;
        }
        // G=P/V^2=I/V, with V=1 V. Both derive from the spatial field, not a
        // fitted bulk length/area approximation or a named material preset.
        let conductance_s = finite(field.power_w / (excitation_v * excitation_v))?;
        let current_conductance_s = finite(terminal_current_a / excitation_v)?;
        let allowance_s = controls.current_tolerance_a / excitation_v
            + controls.power_tolerance_w / (excitation_v * excitation_v);
        let resistance_ohm = conductance_s.recip();
        if !positive(conductance_s)
            || !positive(resistance_ohm)
            || terminal_current_a <= 0.0
            || (current_conductance_s - conductance_s).abs() > allowance_s
        {
            return Err(ConductionError::Invalid(
                "nonpositive or unbalanced two-terminal conductance",
            ));
        }
        poll(cancelled)?;
        Ok(PlanarResistance {
            field,
            resistance_ohm,
        })
    }
}

impl PlanarConduction<ResolvedConductorState> {
    /// Original scalar conductivity view; anisotropic fields expose the full
    /// matrix through `conductivity_coefficients`, never a scalar surrogate.
    pub fn conductivity_s_m(&self) -> &[f64] {
        &self.conductivity_s_m
    }
}

/// Accepted voltage/current/heat tied to its original source and geometry.
pub struct ConductionField<'a, C = ResolvedConductorState> {
    problem: &'a PlanarConduction<C>,
    voltage_v: Vec<f64>,
    electric_field_v_m: Vec<[f64; 2]>,
    current_density_a_m2: Vec<[f64; 2]>,
    cell_power_w: Vec<f64>,
    power_w: f64,
    nodal_boundary_current_a: Vec<Option<f64>>,
    nodal_applied_current_a: Vec<Option<f64>>,
    boundary_current_a: f64,
    maximum_free_current_residual_a: f64,
    power_balance_residual_w: f64,
    report: Option<SolveReport>,
}
impl<C: PlanarConductorMaterial> ConductionField<'_, C> {
    pub fn problem(&self) -> &PlanarConduction<C> {
        self.problem
    }
    pub fn voltage_v(&self) -> &[f64] {
        &self.voltage_v
    }
    /// Original cell E=-grad(V) [V/m]; no interfacial smoothing.
    pub fn electric_field_v_m(&self) -> &[[f64; 2]] {
        &self.electric_field_v_m
    }
    /// Original cell J=sigma E [A/m2].
    pub fn current_density_a_m2(&self) -> &[[f64; 2]] {
        &self.current_density_a_m2
    }
    /// Physical irreversible heat rate, integrated over each extruded cell [W].
    pub fn cell_power_w(&self) -> &[f64] {
        &self.cell_power_w
    }
    pub fn power_w(&self) -> f64 {
        self.power_w
    }
    /// Combined voltage-source reaction and prescribed trace current [A],
    /// positive INTO the conductor. Unloaded, unconstrained vertices have None.
    /// Current-only trace vertices retain the authored weak load; voltage
    /// vertices retain the field reaction. Both contribute to terminal work.
    pub fn nodal_boundary_current_a(&self) -> &[Option<f64>] {
        &self.nodal_boundary_current_a
    }
    /// Integrated prescribed trace loads [A], including any shared voltage
    /// corner. Subtract this from a voltage vertex's combined boundary current
    /// to obtain that voltage source's own reaction. None means no current trace.
    pub fn nodal_applied_current_a(&self) -> &[Option<f64>] {
        &self.nodal_applied_current_a
    }
    pub fn boundary_current_a(&self) -> f64 {
        self.boundary_current_a
    }
    pub fn maximum_free_current_residual_a(&self) -> f64 {
        self.maximum_free_current_residual_a
    }
    /// Integrated Joule power minus sum(terminal voltage * injected current).
    pub fn power_balance_residual_w(&self) -> f64 {
        self.power_balance_residual_w
    }
    pub fn report(&self) -> Option<&SolveReport> {
        self.report.as_ref()
    }
}

/// Retained unit-voltage solution and its independently power/current-checked
/// resistance for the existing Kirchhoff/pHS circuit owner.
pub struct PlanarResistance<'a, C = ResolvedConductorState> {
    field: ConductionField<'a, C>,
    resistance_ohm: f64,
}
impl<'a, C: PlanarConductorMaterial> PlanarResistance<'a, C> {
    pub fn unit_voltage_field(&self) -> &ConductionField<'a, C> {
        &self.field
    }
    pub fn resistance_ohm(&self) -> f64 {
        self.resistance_ohm
    }
    /// This is the conductor's own dissipative branch, not a second parallel
    /// copy. Retain this reduction to reconstruct its spatial heat afterwards.
    pub fn circuit_branch(&self) -> Branch {
        Branch::Resistor {
            ohms: self.resistance_ohm,
        }
    }

    /// Recover spatial current and heat at the native circuit's actual branch
    /// voltage without another field solve. For a midpoint circuit step use
    /// its midpoint voltage; dt*cell_power then partitions that SAME resistor
    /// heat. Do not also debit it as an additional circuit loss.
    ///
    /// Negative voltages reverse current but preserve nonnegative heating.
    /// Zero represents the negative terminal's reference. Actual-drive current
    /// and power acceptance are independent of unit-field extraction budgets.
    pub fn field_at_voltage(
        &self,
        voltage_v: f64,
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionField<'a, C>, ConductionError> {
        poll(cancelled)?;
        if !voltage_v.is_finite() || !positive(current_tolerance_a) || !positive(power_tolerance_w)
        {
            return Err(ConductionError::Invalid(
                "finite branch voltage and positive actual-field tolerances required",
            ));
        }
        let mut voltage = Vec::with_capacity(self.field.voltage_v.len());
        let mut prescribed = Vec::with_capacity(self.field.voltage_v.len());
        for (&unit, current) in self
            .field
            .voltage_v
            .iter()
            .zip(&self.field.nodal_boundary_current_a)
        {
            poll(cancelled)?;
            let value = finite(unit * voltage_v)?;
            voltage.push(value);
            prescribed.push(current.map(|_| value));
        }
        self.field.problem.recover_field(
            voltage,
            &prescribed,
            vec![None; self.field.voltage_v.len()],
            current_tolerance_a,
            power_tolerance_w,
            None,
            cancelled,
        )
    }
}
