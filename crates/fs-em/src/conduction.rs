//! Spatial steady Ohmic conduction from complete material source snapshots.
//!
//! The conductivity Hodge solves -div(sigma grad(V))=0 on the original planar
//! mesh. E=-grad(V), J=sigma E and the actual volume integral E.J determine
//! current and irreversible heating. The dielectric field shares only scalar
//! algebra and electrode geometry; its charge/storage law is not reused here.
//!
//! This is a stationary, isotropic DC model with fixed geometry and material
//! state. Induction, displacement current, contacts, Hall response and nonlinear
//! field-dependent conduction require their own models. A changed temperature
//! or material must be resolved and prepared again, not substituted into an
//! already accepted field. No electrical storage or thermal evolution is added.

use fs_feec::{
    assembly::{incidence_to_csr, stiffness},
    weighted::{WeightedAssemblyLimits, WeightedError},
    weighted2::{CellWeight2, weighted_mass_matrix_2d},
};
use fs_material::conductor::{ELECTRICAL_CONDUCTIVITY_DIMS, ResolvedConductorState};
use fs_phs::circuit::Branch;
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::{CheckedCgConfig, SolveReport};
use fs_sparse::Csr;

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

/// Immutable original mesh, source-bound scalar conductors and conductivity
/// Hodge. One complete material snapshot is retained per original triangle.
pub struct PlanarConduction {
    mesh: TriComplex2,
    conductors: Vec<ResolvedConductorState>,
    conductivity_s_m: Vec<f64>,
    stiffness: Csr,
}
impl PlanarConduction {
    /// Assemble the same native weighted FEEC operator used by other scalar
    /// fields, with the conductivity's actual SI dimensions. The canonical
    /// conductor owner selects resistivity or conductivity explicitly; material
    /// names and an assumed unit wire play no role in the spatial law.
    pub fn new(
        mesh: TriComplex2,
        conductors: Vec<ResolvedConductorState>,
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
        let mut weights = Vec::with_capacity(conductors.len());
        for conductor in &conductors {
            poll(cancelled)?;
            if !conductor
                .source_property()
                .is_static_at_fixed_external_conditions()
            {
                return Err(ConductionError::Invalid(
                    "spatial DC conduction requires a static Ohmic scalar at fixed external conditions",
                ));
            }
            let sigma = conductor.resistivity_ohm_m().recip();
            if !positive(sigma) {
                return Err(ConductionError::Invalid(
                    "source conductivity must be positive and representable",
                ));
            }
            conductivity_s_m.push(sigma);
            weights.push(CellWeight2::Scalar(sigma));
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
            conductivity_s_m,
            stiffness,
        })
    }

    pub fn mesh(&self) -> &TriComplex2 {
        &self.mesh
    }
    /// Complete original card, query, requirements and electrical source basis.
    pub fn conductors(&self) -> &[ResolvedConductorState] {
        &self.conductors
    }
    pub fn conductivity_s_m(&self) -> &[f64] {
        &self.conductivity_s_m
    }

    /// Solve at prescribed exterior voltage traces. Remaining exterior edges
    /// are electrically insulated; each connected component needs a complete
    /// anchored trace. No volume-current source or point electrode is inferred.
    pub fn solve(
        &self,
        boundary_voltage_v: &[(usize, f64)],
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionField<'_>, ConductionError> {
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
        let (voltage, report) = scalar_potential::solve(
            &self.stiffness,
            &vec![0.0; self.mesh.vertices().len()],
            &prescribed,
            controls.linear,
            cancelled,
        )?;
        self.recover_field(
            voltage,
            &prescribed,
            controls.current_tolerance_a,
            controls.power_tolerance_w,
            report,
            cancelled,
        )
    }

    fn recover_field(
        &self,
        voltage: Vec<f64>,
        prescribed: &[Option<f64>],
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        report: Option<SolveReport>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionField<'_>, ConductionError> {
        let n = self.mesh.vertices().len();
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
            let j = [
                finite(self.conductivity_s_m[cell] * e[0])?,
                finite(self.conductivity_s_m[cell] * e[1])?,
            ];
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
            if prescribed[v].is_some() {
                nodal_boundary_current_a[v] = Some(reaction);
                boundary_current_a = finite(boundary_current_a + reaction)?;
                let difference = finite(voltage[v] - reference)?;
                terminal_power_w = finite(difference.mul_add(reaction, terminal_power_w))?;
            } else {
                maximum_free_current_residual_a =
                    maximum_free_current_residual_a.max(reaction.abs());
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
    ) -> Result<PlanarResistance<'_>, ConductionError> {
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

/// Accepted voltage/current/heat tied to its original source and geometry.
pub struct ConductionField<'a> {
    problem: &'a PlanarConduction,
    voltage_v: Vec<f64>,
    electric_field_v_m: Vec<[f64; 2]>,
    current_density_a_m2: Vec<[f64; 2]>,
    cell_power_w: Vec<f64>,
    power_w: f64,
    nodal_boundary_current_a: Vec<Option<f64>>,
    boundary_current_a: f64,
    maximum_free_current_residual_a: f64,
    power_balance_residual_w: f64,
    report: Option<SolveReport>,
}
impl ConductionField<'_> {
    pub fn problem(&self) -> &PlanarConduction {
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
    /// Weak terminal reactions [A], positive INTO the conductor. Free vertices
    /// have None; they are not fictitious terminals or point current sources.
    pub fn nodal_boundary_current_a(&self) -> &[Option<f64>] {
        &self.nodal_boundary_current_a
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
pub struct PlanarResistance<'a> {
    field: ConductionField<'a>,
    resistance_ohm: f64,
}
impl<'a> PlanarResistance<'a> {
    pub fn unit_voltage_field(&self) -> &ConductionField<'a> {
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
    ) -> Result<ConductionField<'a>, ConductionError> {
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
            current_tolerance_a,
            power_tolerance_w,
            None,
            cancelled,
        )
    }
}
