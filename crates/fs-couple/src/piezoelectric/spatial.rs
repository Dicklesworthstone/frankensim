//! Reciprocal small-strain piezoelectricity on conforming P1 tetrahedra.
//!
//! The native `fs-solid` stiffness/mass and `fs-em` dielectric weak form share
//! the same mesh and source-bound tensor. With `E = -grad(phi)`, the coupled
//! static equations are `[K G; G^T -K_e] [u,phi] = [f,-q]`. Both off-diagonal
//! blocks come from ONE `G = integral B^T e^T grad(N) dV`; neither a sensor gain
//! nor an independently calibrated actuator force is introduced.
//!
//! Electrical electrodes are equipotential faces. A voltage electrode fixes
//! its potential; a charge electrode solves one shared potential at a declared
//! TOTAL free charge (zero is an isolated, initially uncharged electrode).
//! Remaining exterior faces are electrically insulating. A voltage reference
//! is required on each connected body. Coefficients and reference geometry are
//! fixed. This is a static linear model, not a transient, ferroelectric,
//! breakdown, loss, acoustic-radiation or exterior-fringing-field model.

use fs_blake3::ContentHash;
use fs_em::electrostatic::tetrahedral::{
    TetDielectricBudget, TetDielectricError, TetDielectricOperator,
};
use fs_exec::Cx;
use fs_material::piezoelectric::PiezoelectricError;
use fs_material::piezoelectric::tensor::{TensorPiezoelectric, TensorPiezoelectricResponse};
use fs_rep_mesh::TetComplex;
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticAssembly, TetElasticError, TetElasticMaterial,
    TetLinearElasticProblem, TetMaterialField,
};
use fs_solver::{CsrOp, LinearOp, MinresState};
use fs_sparse::{Coo, Csr};

mod boundary;
use boundary::{BoundaryMap, prepare_boundary};

/// Electrical condition on a complete conducting surface electrode.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ElectrodeCondition {
    /// Prescribed potential [V], including the chosen voltage reference.
    Voltage(f64),
    /// Prescribed total electrode free charge [C]; potential is solved.
    Charge(f64),
}

/// One equipotential electrode, as a disjoint union of complete exterior faces.
#[derive(Clone, Copy, Debug)]
pub struct PiezoelectricElectrode<'a> {
    /// Original mesh vertices; duplicates and overlap with other electrodes refuse.
    pub nodes: &'a [usize],
    /// Voltage control or isolated fixed-charge control.
    pub condition: ElectrodeCondition,
}

/// Explicit accuracy and small-signal envelope for a static coupled solve.
#[derive(Clone, Copy, Debug)]
pub struct TetPiezoelectricControls {
    /// Relative residual of the diagonally equilibrated coupled system, in (0,1).
    pub relative_tolerance: f64,
    /// Maximum native MINRES iterations, including residual correction solves.
    pub maximum_iterations: usize,
    /// Maximum absolute residual on an unconstrained mechanical DOF [N].
    pub force_tolerance_n: f64,
    /// Maximum absolute free-node or total floating-electrode charge residual [C].
    pub charge_tolerance_c: f64,
    /// Maximum absolute error in the independent ramp-work/storage identity [J].
    pub energy_tolerance_j: f64,
    /// Maximum Frobenius norm of each infinitesimal strain tensor, in (0,1).
    pub maximum_strain_norm: f64,
    /// Maximum Euclidean magnitude of each electric field [V/m].
    pub maximum_electric_field_v_m: f64,
}

/// Source, mesh, boundary, amplitude or numerical refusal; no partial field escapes.
#[derive(Debug)]
pub enum TetPiezoelectricError {
    /// The source tensor cannot evaluate the requested state.
    Material(PiezoelectricError),
    /// Native solid assembly or support admission failed.
    Solid(TetElasticError),
    /// Native dielectric assembly or field reconstruction failed.
    Electric(TetDielectricError),
    /// Invalid complete input, unsupported condition or nonfinite arithmetic.
    Invalid(&'static str),
    /// Native iteration failed the independently recomputed scaled residual.
    Nonconvergence {
        /// Native iteration count.
        iterations: usize,
        /// Recomputed relative residual.
        relative_residual: f64,
    },
    /// Independent mechanical/electrical residual or work check failed.
    Balance(&'static str),
    /// Cancellation observed before publishing any candidate field.
    Cancelled,
}
impl From<PiezoelectricError> for TetPiezoelectricError {
    fn from(value: PiezoelectricError) -> Self {
        Self::Material(value)
    }
}
impl From<TetElasticError> for TetPiezoelectricError {
    fn from(value: TetElasticError) -> Self {
        Self::Solid(value)
    }
}
impl From<TetDielectricError> for TetPiezoelectricError {
    fn from(value: TetDielectricError) -> Self {
        Self::Electric(value)
    }
}
impl core::fmt::Display for TetPiezoelectricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tetrahedral piezoelectricity: {self:?}")
    }
}
impl std::error::Error for TetPiezoelectricError {}
type Result<T> = core::result::Result<T, TetPiezoelectricError>;
fn invalid(reason: &'static str) -> TetPiezoelectricError {
    TetPiezoelectricError::Invalid(reason)
}
fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint()
        .map_err(|_| TetPiezoelectricError::Cancelled)
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

/// Prepared, immutable native operators. Assembly is reused across voltage,
/// charge and mechanical load cases; no tensor inversion occurs in the solve.
pub struct TetPiezoelectricModel<'a> {
    nodes_m: &'a [[f64; 3]],
    tetrahedra: &'a [[usize; 4]],
    materials: &'a [TensorPiezoelectric],
    elastic_materials: Vec<TetElasticMaterial>,
    mesh_frame: ContentHash,
    budget: TetAssemblyBudget,
    solid: TetElasticAssembly,
    dielectric: TetDielectricOperator,
    coupling: Csr,
}

/// One complete solved cell; Mandel ordering is xx,yy,zz,sqrt(2)xy,yz,zx.
#[derive(Clone, Copy, Debug)]
pub struct TetPiezoelectricCell {
    /// Infinitesimal strain in the declared mesh frame.
    pub strain_mandel: [f64; 6],
    /// Stress, field, displacement and storage from the original source energy.
    pub response: TensorPiezoelectricResponse,
    /// Stored energy integrated over this reference tetrahedron [J].
    pub stored_energy_j: f64,
}

/// Electrode charge and voltage recovered from the coupled field.
#[derive(Clone, Copy, Debug)]
pub struct PiezoelectricElectrodeResponse {
    /// Shared equipotential value [V].
    pub voltage_v: f64,
    /// Total conductor charge, excluding separately supplied volume charge [C].
    pub charge_c: f64,
}

/// Accepted displacement, potential and work-conjugate reactions.
#[derive(Clone, Debug)]
pub struct TetPiezoelectricSolution {
    /// Solved displacement at every original node [m].
    pub displacement_m: Vec<[f64; 3]>,
    /// Solved potential at every original node [V].
    pub voltage_v: Vec<f64>,
    /// Cell state in original connectivity order.
    pub cells: Vec<TetPiezoelectricCell>,
    /// Internal nodal force minus applied force [N]; free entries are residuals.
    pub support_force_n: Vec<[f64; 3]>,
    /// Electric reaction minus applied volume charge [C]. Floating electrode
    /// entries sum to its specified charge; their distribution is solved.
    pub boundary_charge_c: Vec<f64>,
    /// Conducting electrode responses, in caller order.
    pub electrodes: Vec<PiezoelectricElectrodeResponse>,
    /// Positive internal elastic plus electric storage [J], not saddle enthalpy.
    pub stored_energy_j: f64,
    /// One half displacement dot total mechanical load [J].
    pub mechanical_ramp_work_j: f64,
    /// One half potential dot total free charge [J].
    pub electrical_ramp_work_j: f64,
    /// Storage minus the two proportional-load ramp works [J].
    pub energy_residual_j: f64,
    /// Maximum unconstrained-force residual [N].
    pub maximum_force_residual_n: f64,
    /// Maximum free-node/total floating-electrode charge residual [C].
    pub maximum_charge_residual_c: f64,
    /// Native iterations actually used.
    pub iterations: usize,
    /// Independently recomputed equilibrated relative residual.
    pub relative_residual: f64,
}

impl<'a> TetPiezoelectricModel<'a> {
    /// Prepare the native P1 solid/dielectric/coupling operators. `materials`
    /// contains either one uniform source or exactly one source per tetrahedron.
    /// Every tensor must already be rotated into `mesh_frame` explicitly.
    pub fn assemble(
        nodes_m: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        materials: &'a [TensorPiezoelectric],
        mesh_frame: ContentHash,
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        poll(cx)?;
        if mesh_frame == ContentHash([0; 32])
            || materials.is_empty()
            || !(materials.len() == 1 || materials.len() == tetrahedra.len())
            || materials
                .iter()
                .any(|material| material.frame_identity() != mesh_frame)
        {
            return Err(invalid(
                "complete piezoelectric materials in the declared mesh frame required",
            ));
        }
        let elastic_materials = materials
            .iter()
            .map(|material| {
                TetElasticMaterial::try_new_mandel(
                    material.density_kg_m3(),
                    *material.stiffness_mandel_pa(),
                    material.identity().0,
                )
            })
            .collect::<core::result::Result<Vec<_>, _>>()?;
        let field = if elastic_materials.len() == 1 {
            TetMaterialField::Uniform(&elastic_materials[0])
        } else {
            TetMaterialField::PerElement(&elastic_materials)
        };
        // Keep full native matrices: nonzero mechanical displacements and
        // floating electrodes are eliminated together for each load case.
        let solid = TetLinearElasticProblem {
            nodes_m,
            tetrahedra,
            materials: field,
            fixed_dofs: &[],
            budget,
        }
        .assemble(cx)?;
        let connectivity = tetrahedra
            .iter()
            .map(|tet| {
                let mut out = [0u32; 4];
                for i in 0..4 {
                    out[i] = u32::try_from(tet[i])
                        .map_err(|_| invalid("node index exceeds native complex capacity"))?;
                }
                Ok(out)
            })
            .collect::<Result<Vec<_>>>()?;
        let mesh = TetComplex::from_tets(nodes_m.len(), connectivity);
        let epsilon = (0..tetrahedra.len())
            .map(|cell| *materials[if materials.len() == 1 { 0 } else { cell }].permittivity_f_m())
            .collect::<Vec<_>>();
        let dielectric = TetDielectricOperator::assemble(
            &mesh,
            nodes_m,
            &epsilon,
            TetDielectricBudget {
                maximum_nodes: budget.maximum_nodes,
                maximum_tetrahedra: budget.maximum_tetrahedra,
                minimum_scaled_jacobian: budget.minimum_scaled_jacobian,
            },
            &mut || cx.checkpoint().is_err(),
        )?;
        let mut coupling = Coo::new(3 * nodes_m.len(), nodes_m.len());
        for (cell, tet) in tetrahedra.iter().enumerate() {
            poll(cx)?;
            let b = strain_matrix(&dielectric.gradients()[cell]);
            let material = &materials[if materials.len() == 1 { 0 } else { cell }];
            let e = material.coupling_c_m2();
            let volume = dielectric.volumes_m3()[cell];
            for local in 0..12 {
                let mut polarization = [0.0; 3];
                for i in 0..3 {
                    for a in 0..6 {
                        polarization[i] += e[i][a] * b[a][local];
                    }
                }
                for node in 0..4 {
                    let value = volume * dot3(polarization, dielectric.gradients()[cell][node]);
                    if !value.is_finite() {
                        return Err(invalid("coupled element is not representable"));
                    }
                    if value != 0.0 {
                        coupling.push(3 * tet[local / 3] + local % 3, tet[node], value);
                    }
                }
            }
        }
        poll(cx)?;
        Ok(Self {
            nodes_m,
            tetrahedra,
            materials,
            elastic_materials,
            mesh_frame,
            budget,
            solid,
            dielectric,
            coupling: coupling.assemble(),
        })
    }

    /// Original mesh coordinate frame shared by all three constitutive blocks.
    pub const fn mesh_frame(&self) -> ContentHash {
        self.mesh_frame
    }
    /// Original source tensor, uniform or in connectivity order.
    pub const fn materials(&self) -> &[TensorPiezoelectric] {
        self.materials
    }
    /// Native full mechanical stiffness and consistent mass; no electrode clamp
    /// has yet condensed the coupling into this short-circuit stiffness.
    pub const fn solid(&self) -> &TetElasticAssembly {
        &self.solid
    }
    /// Native constant-strain dielectric operator on the same original mesh.
    pub const fn dielectric(&self) -> &TetDielectricOperator {
        &self.dielectric
    }
    /// Work-conjugate coupling [C/m], rows mechanical, columns electrical.
    pub const fn coupling(&self) -> &Csr {
        &self.coupling
    }

    /// Solve displacement and electrostatics together. Mechanical constraints
    /// are `(3*node+axis, displacement_m)`. Nodal charge is a complete integrated
    /// volume-load vector; electrode total charge is supplied separately.
    /// All inputs and physical residuals are checked before returning any field.
    #[allow(clippy::too_many_arguments)]
    pub fn solve(
        &self,
        prescribed_displacements_m: &[(usize, f64)],
        electrodes: &[PiezoelectricElectrode<'_>],
        nodal_forces_n: &[[f64; 3]],
        nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<TetPiezoelectricSolution> {
        poll(cx)?;
        let n = self.nodes_m.len();
        if !positive(controls.relative_tolerance)
            || controls.relative_tolerance >= 1.0
            || controls.maximum_iterations == 0
            || !positive(controls.force_tolerance_n)
            || !positive(controls.charge_tolerance_c)
            || !positive(controls.energy_tolerance_j)
            || !positive(controls.maximum_strain_norm)
            || controls.maximum_strain_norm >= 1.0
            || !positive(controls.maximum_electric_field_v_m)
            || nodal_forces_n.len() != n
            || nodal_free_charge_c.len() != n
            || nodal_forces_n
                .iter()
                .flatten()
                .chain(nodal_free_charge_c)
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "complete finite loads and positive finite solve/amplitude controls required",
            ));
        }
        let boundary = prepare_boundary(
            self.nodes_m,
            self.tetrahedra,
            prescribed_displacements_m,
            electrodes,
            nodal_forces_n,
            nodal_free_charge_c,
            cx,
        )?;
        let fixed_dofs = prescribed_displacements_m
            .iter()
            .map(|&(dof, _)| dof)
            .collect::<Vec<_>>();
        let field = if self.elastic_materials.len() == 1 {
            TetMaterialField::Uniform(&self.elastic_materials[0])
        } else {
            TetMaterialField::PerElement(&self.elastic_materials)
        };
        TetLinearElasticProblem {
            nodes_m: self.nodes_m,
            tetrahedra: self.tetrahedra,
            materials: field,
            fixed_dofs: &fixed_dofs,
            budget: self.budget,
        }
        .validate_static_constraints(cx)?;
        let (matrix, rhs) = self.reduced_system(&boundary, cx)?;
        let (solution, iterations, relative_residual) =
            equilibrated_solve(matrix, rhs, controls, cx)?;
        let mut full = boundary.prescribed.clone();
        for (i, index) in boundary.index.iter().enumerate() {
            if let Some(index) = index {
                full[i] = solution[*index];
            }
        }
        // Restore exact user voltage values after solving only voltage differences.
        let voltage_v = full[3 * n..]
            .iter()
            .zip(&boundary.voltage_reference)
            .map(|(v, reference)| v + reference)
            .collect::<Vec<_>>();
        let displacement_m = full[..3 * n]
            .chunks_exact(3)
            .map(|u| [u[0], u[1], u[2]])
            .collect::<Vec<_>>();
        self.recover(
            &full,
            displacement_m,
            voltage_v,
            &boundary,
            electrodes,
            nodal_forces_n,
            nodal_free_charge_c,
            controls,
            iterations,
            relative_residual,
            cx,
        )
    }

    fn reduced_system(&self, boundary: &BoundaryMap, cx: &Cx<'_>) -> Result<(Csr, Vec<f64>)> {
        let mut matrix = Coo::new(boundary.rhs.len(), boundary.rhs.len());
        let mut rhs = boundary.rhs.clone();
        let mut add = |row: usize, col: usize, value: f64| -> Result<()> {
            if !value.is_finite() {
                return Err(invalid("nonfinite native matrix coefficient"));
            }
            if let Some(r) = boundary.index[row] {
                if let Some(c) = boundary.index[col] {
                    matrix.push(r, c, value);
                } else {
                    rhs[r] -= value * boundary.prescribed[col];
                }
            }
            Ok(())
        };
        let nu = 3 * self.nodes_m.len();
        for row in 0..nu {
            poll(cx)?;
            let (cols, vals) = self.solid.stiffness.row(row);
            for (&col, &value) in cols.iter().zip(vals) {
                add(row, col, value)?;
            }
            let (cols, vals) = self.coupling.row(row);
            for (&col, &value) in cols.iter().zip(vals) {
                add(row, nu + col, value)?;
                add(nu + col, row, value)?;
            }
        }
        for row in 0..self.nodes_m.len() {
            poll(cx)?;
            let (cols, vals) = self.dielectric.stiffness().row(row);
            for (&col, &value) in cols.iter().zip(vals) {
                add(nu + row, nu + col, -value)?;
            }
        }
        if rhs.iter().any(|v| !v.is_finite()) {
            return Err(invalid("reduced load is not representable"));
        }
        Ok((matrix.assemble(), rhs))
    }

    #[allow(clippy::too_many_arguments)]
    fn recover(
        &self,
        full: &[f64],
        displacement_m: Vec<[f64; 3]>,
        voltage_v: Vec<f64>,
        boundary: &BoundaryMap,
        electrodes: &[PiezoelectricElectrode<'_>],
        forces: &[[f64; 3]],
        charges: &[f64],
        controls: TetPiezoelectricControls,
        iterations: usize,
        relative_residual: f64,
        cx: &Cx<'_>,
    ) -> Result<TetPiezoelectricSolution> {
        let n = self.nodes_m.len();
        let electric = self
            .dielectric
            .electric_field_v_m(&full[3 * n..], &mut || cx.checkpoint().is_err())?;
        let mut internal = vec![[0.0; 3]; n];
        let mut nodal_charge = vec![0.0; n];
        let mut cells = Vec::with_capacity(self.tetrahedra.len());
        let mut stored_energy_j = 0.0;
        for (cell, tet) in self.tetrahedra.iter().enumerate() {
            poll(cx)?;
            let gradients = &self.dielectric.gradients()[cell];
            let b = strain_matrix(gradients);
            let mut strain = [0.0; 6];
            for a in 0..6 {
                for local in 0..12 {
                    // Relative nodal displacement removes rigid translation from
                    // small-strain reconstruction before multiplication.
                    strain[a] += b[a][local]
                        * (displacement_m[tet[local / 3]][local % 3]
                            - displacement_m[tet[0]][local % 3]);
                }
            }
            if norm(&strain) > controls.maximum_strain_norm
                || norm(&electric[cell]) > controls.maximum_electric_field_v_m
            {
                return Err(invalid(
                    "solved field leaves the declared small-signal envelope",
                ));
            }
            let material = &self.materials[if self.materials.len() == 1 { 0 } else { cell }];
            let response = material.field_response(strain, electric[cell])?;
            let volume = self.dielectric.volumes_m3()[cell];
            let energy = response.stored_energy_j_m3 * volume;
            stored_energy_j += energy;
            for local in 0..12 {
                let value = (0..6)
                    .map(|a| b[a][local] * response.stress_mandel_pa[a])
                    .sum::<f64>();
                internal[tet[local / 3]][local % 3] += volume * value;
            }
            for local in 0..4 {
                nodal_charge[tet[local]] -=
                    volume * dot3(gradients[local], response.displacement_c_m2);
            }
            cells.push(TetPiezoelectricCell {
                strain_mandel: strain,
                response,
                stored_energy_j: energy,
            });
        }
        let mut maximum_force_residual_n: f64 = 0.0;
        let mut support_force_n = internal.clone();
        for node in 0..n {
            for axis in 0..3 {
                support_force_n[node][axis] -= forces[node][axis];
                if boundary.index[3 * node + axis].is_some() {
                    maximum_force_residual_n =
                        maximum_force_residual_n.max(support_force_n[node][axis].abs());
                }
            }
        }
        let boundary_charge_c = nodal_charge
            .iter()
            .zip(charges)
            .map(|(q, b)| q - b)
            .collect::<Vec<_>>();
        let mut maximum_charge_residual_c: f64 = 0.0;
        let mut electric_residual = vec![0.0; boundary.rhs.len()];
        for node in 0..n {
            if let Some(index) = boundary.index[3 * n + node] {
                electric_residual[index] += boundary_charge_c[node];
            }
        }
        let mut electrode_responses = Vec::with_capacity(electrodes.len());
        for electrode in electrodes {
            let charge_c = electrode
                .nodes
                .iter()
                .map(|&node| boundary_charge_c[node])
                .sum::<f64>();
            if let ElectrodeCondition::Charge(charge) = electrode.condition {
                let index =
                    boundary.index[3 * n + electrode.nodes[0]].expect("admitted charge electrode");
                electric_residual[index] -= charge;
            }
            electrode_responses.push(PiezoelectricElectrodeResponse {
                voltage_v: voltage_v[electrode.nodes[0]],
                charge_c,
            });
        }
        for value in electric_residual {
            maximum_charge_residual_c = maximum_charge_residual_c.max(value.abs());
        }
        // Use the actual imposed loads on free coordinates. Reactions supply
        // work only where a displacement or voltage source is prescribed;
        // otherwise this check would merely repeat the constitutive Euler
        // identity even for a field that did not solve the boundary-value problem.
        let mut mechanical_ramp_work_j = 0.0;
        for node in 0..n {
            for axis in 0..3 {
                let load = if boundary.index[3 * node + axis].is_some() {
                    forces[node][axis]
                } else {
                    internal[node][axis]
                };
                mechanical_ramp_work_j += 0.5 * displacement_m[node][axis] * load;
            }
        }
        // Relative potentials avoid contaminating work with arbitrary gauges.
        let mut electrical_ramp_work_j = 0.5
            * full[3 * n..]
                .iter()
                .zip(charges)
                .map(|(v, q)| v * q)
                .sum::<f64>();
        for (electrode, response) in electrodes.iter().zip(&electrode_responses) {
            let charge = match electrode.condition {
                ElectrodeCondition::Voltage(_) => response.charge_c,
                ElectrodeCondition::Charge(charge) => charge,
            };
            electrical_ramp_work_j += 0.5 * full[3 * n + electrode.nodes[0]] * charge;
        }
        let energy_residual_j = stored_energy_j - mechanical_ramp_work_j - electrical_ramp_work_j;
        if [
            stored_energy_j,
            mechanical_ramp_work_j,
            electrical_ramp_work_j,
            energy_residual_j,
            maximum_force_residual_n,
            maximum_charge_residual_c,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || internal
                .iter()
                .flatten()
                .chain(&boundary_charge_c)
                .any(|v| !v.is_finite())
            || voltage_v.iter().any(|v| !v.is_finite())
        {
            return Err(invalid("solved response or work is not representable"));
        }
        if maximum_force_residual_n > controls.force_tolerance_n {
            return Err(TetPiezoelectricError::Balance(
                "unconstrained mechanical force residual exceeds tolerance",
            ));
        }
        if maximum_charge_residual_c > controls.charge_tolerance_c {
            return Err(TetPiezoelectricError::Balance(
                "free-node or electrode-total charge residual exceeds tolerance",
            ));
        }
        if energy_residual_j.abs() > controls.energy_tolerance_j {
            return Err(TetPiezoelectricError::Balance(
                "independent internal energy and ramp work do not close",
            ));
        }
        poll(cx)?;
        Ok(TetPiezoelectricSolution {
            displacement_m,
            voltage_v,
            cells,
            support_force_n,
            boundary_charge_c,
            electrodes: electrode_responses,
            stored_energy_j,
            mechanical_ramp_work_j,
            electrical_ramp_work_j,
            energy_residual_j,
            maximum_force_residual_n,
            maximum_charge_residual_c,
            iterations,
            relative_residual,
        })
    }
}

fn strain_matrix(g: &[[f64; 3]; 4]) -> [[f64; 12]; 6] {
    let mut b = [[0.0; 12]; 6];
    for (node, &[x, y, z]) in g.iter().enumerate() {
        let i = 3 * node;
        b[0][i] = x;
        b[1][i + 1] = y;
        b[2][i + 2] = z;
        b[3][i] = y * core::f64::consts::FRAC_1_SQRT_2;
        b[3][i + 1] = x * core::f64::consts::FRAC_1_SQRT_2;
        b[4][i + 1] = z * core::f64::consts::FRAC_1_SQRT_2;
        b[4][i + 2] = y * core::f64::consts::FRAC_1_SQRT_2;
        b[5][i] = z * core::f64::consts::FRAC_1_SQRT_2;
        b[5][i + 2] = x * core::f64::consts::FRAC_1_SQRT_2;
    }
    b
}
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |sum, &x| sum.hypot(x))
}

fn equilibrated_solve(
    matrix: Csr,
    rhs: Vec<f64>,
    controls: TetPiezoelectricControls,
    cx: &Cx<'_>,
) -> Result<(Vec<f64>, usize, f64)> {
    let n = rhs.len();
    if n == 0 {
        return Ok((Vec::new(), 0, 0.0));
    }
    let mut scale = vec![0.0; n];
    for (row, s) in scale.iter_mut().enumerate() {
        poll(cx)?;
        let (cols, vals) = matrix.row(row);
        let diag = cols
            .iter()
            .position(|&col| col == row)
            .map(|i| vals[i].abs())
            .unwrap_or(0.0);
        if !positive(diag) {
            return Err(invalid(
                "free coordinate has zero or nonfinite stiffness/capacitance",
            ));
        }
        *s = 1.0 / fs_math::det::sqrt(diag);
        if !positive(*s) {
            return Err(invalid("coupled diagonal scaling is not representable"));
        }
    }
    let mut balanced = Coo::new(n, n);
    for row in 0..n {
        poll(cx)?;
        let (cols, vals) = matrix.row(row);
        for (&col, &value) in cols.iter().zip(vals) {
            let entry = value * scale[row] * scale[col];
            if !entry.is_finite() {
                return Err(invalid("equilibrated matrix is not representable"));
            }
            balanced.push(row, col, entry);
        }
    }
    let rhs = rhs
        .iter()
        .zip(&scale)
        .map(|(b, s)| b * s)
        .collect::<Vec<_>>();
    if rhs.iter().any(|v| !v.is_finite()) {
        return Err(invalid("equilibrated load is not representable"));
    }
    let rhs_norm = norm(&rhs);
    if rhs_norm == 0.0 {
        return Ok((vec![0.0; n], 0, 0.0));
    }
    if !rhs_norm.is_finite() {
        return Err(invalid("equilibrated load norm is not representable"));
    }
    let operator = CsrOp::symmetric(balanced.assemble());
    let mut solution = vec![0.0; n];
    let mut residual = rhs.clone();
    let mut applied = vec![0.0; n];
    let mut iterations = 0;
    let mut relative_residual = 1.0;
    // Native recurrence estimates are not acceptance evidence. Recompute from
    // the original equilibrated operator; bounded correction solves remove
    // roundoff drift without altering the physical matrix or its boundary data.
    for _ in 0..3 {
        if relative_residual <= controls.relative_tolerance
            || iterations >= controls.maximum_iterations
        {
            break;
        }
        // Native MINRES forms dot products internally. Normalize before
        // entering it so finite, very small/large physical loads cannot make
        // its initial squared norm underflow or overflow.
        let correction_norm = norm(&residual);
        let normalized = residual
            .iter()
            .map(|value| value / correction_norm)
            .collect::<Vec<_>>();
        let mut state = MinresState::new(&operator, &normalized);
        let available = controls.maximum_iterations - iterations;
        while state.iters < available && state.rel_residual() >= controls.relative_tolerance * 0.25
        {
            poll(cx)?;
            state.run(
                &operator,
                controls.relative_tolerance * 0.25,
                (available - state.iters).min(32),
            );
            if state.x.iter().any(|v| !v.is_finite()) {
                break;
            }
        }
        iterations += state.iters;
        for (x, delta) in solution.iter_mut().zip(&state.x) {
            *x += delta * correction_norm;
        }
        operator.apply(&solution, &mut applied);
        for i in 0..n {
            residual[i] = rhs[i] - applied[i];
        }
        relative_residual = norm(&residual) / rhs_norm;
        if !relative_residual.is_finite() {
            break;
        }
    }
    if !relative_residual.is_finite() || relative_residual > controls.relative_tolerance {
        return Err(TetPiezoelectricError::Nonconvergence {
            iterations,
            relative_residual,
        });
    }
    for (x, s) in solution.iter_mut().zip(scale) {
        *x *= s;
    }
    if solution.iter().any(|v| !v.is_finite()) {
        return Err(invalid("physical solution is not representable"));
    }
    poll(cx)?;
    Ok((solution, iterations, relative_residual))
}
