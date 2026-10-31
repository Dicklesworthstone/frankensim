//! Eight-node hexahedral small-strain inelasticity, using the existing
//! isoparametric geometry and incremental Newton owners.
//!
//! Eight Gauss points retain independent native constitutive histories.
//! [`Formulation::Standard`] uses the compatible strain. [`Formulation::BBar`]
//! replaces only its volumetric part by the reference-volume-weighted element
//! mean: `epsilon_bar = dev(epsilon) + mean(trace(epsilon)) I / 3`. The same
//! projected columns enter strain, internal force and the complete algorithmic
//! tangent; the latter is never symmetrized. This is the small-strain B-bar
//! brick formulation, not a general stable displacement/pressure element.
//! See the [OpenSees formulation description](https://opensees.github.io/OpenSeesDocumentation/user/manual/model/elements/bbarBrick.html).
//!
//! Raw compatible strain and rotation remain independently bounded, so the
//! projection cannot hide an inadmissible displacement field. Returned strains
//! are the assumed strains actually supplied to the law. Material parameters
//! remain fixed during a load path; rate-dependent laws must own one physical
//! interval per increment. Histories commit only after global equilibrium.
//!
//! Geometry admission samples the corners, centre and integration points;
//! positive sampled Jacobians do not prove global injectivity of a warped hex.
//! No finite-strain, general shear-locking, limit-point, contact or dynamics
//! claim is made. The full deviatoric quadrature is retained; no one-point
//! stiffness or artificial hourglass penalty is introduced.

use super::increment::{self, IncrementControls, IncrementOperator};
use super::{
    IncrementRecord, IncrementSettings, SmallStrainLawField, SmallStrainTetError, WEIGHT, b_column,
    checkpoint,
};
use crate::linear::Formulation;
use crate::linear3::hex::{HEX27_REFERENCE_NODES, hex8_point};
use fs_exec::Cx;
use fs_feec::highorder::quad1d::gauss_legendre;
use fs_material::{SmallStrainLaw, Voigt};
use fs_sparse::{Coo, Csr};
use std::collections::BTreeSet;

/// Allocation and geometry bounds checked before quadrature preparation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HexAssemblyBudget {
    /// Maximum mesh vertices.
    pub maximum_nodes: usize,
    /// Maximum eight-node cells.
    pub maximum_elements: usize,
    /// Maximum retained integration-point histories (eight per cell).
    pub maximum_quadrature_points: usize,
    /// Maximum unconstrained scalar displacement unknowns.
    pub maximum_free_dofs: usize,
    /// Maximum local matrix entries before sparse duplicate assembly.
    pub maximum_tangent_entries: usize,
    /// Maximum physical increments and retained records.
    pub maximum_load_steps: usize,
    /// Minimum positive `det(J) / longest_corner_distance^3` at each probe.
    pub minimum_scaled_jacobian: f64,
}

impl Default for HexAssemblyBudget {
    fn default() -> Self {
        Self {
            maximum_nodes: 100_000,
            maximum_elements: 100_000,
            maximum_quadrature_points: 800_000,
            maximum_free_dofs: 300_000,
            maximum_tangent_entries: 6_400_000,
            maximum_load_steps: 10_000,
            minimum_scaled_jacobian: 1e-12,
        }
    }
}

/// A bounded three-dimensional load path on conforming trilinear hexes.
#[derive(Debug, Clone, Copy)]
pub struct SmallStrainHexProblem<'a, L: SmallStrainLaw> {
    /// Reference vertex coordinates, m.
    pub nodes_m: &'a [[f64; 3]],
    /// Canonical corners: four around z-, then the corresponding four at z+.
    /// Ordering matches the first eight [`HEX27_REFERENCE_NODES`].
    pub hexahedra: &'a [[usize; 8]],
    /// Native laws in connectivity order; coefficients are not homogenized.
    pub laws: SmallStrainLawField<'a, L>,
    /// Full compatible integration or volumetric B-bar projection.
    pub formulation: Formulation,
    /// Maximum Frobenius norm of both raw and assumed infinitesimal strain.
    pub strain_limit: f64,
    /// Maximum infinitesimal rotation-vector norm, radians.
    pub rotation_limit_rad: f64,
    /// Full-load prescribed displacements: `(3 * node + component, m)`.
    pub prescribed_m: &'a [(usize, f64)],
    /// Full-load dead nodal forces: `(3 * node + component, N)`.
    pub nodal_forces_n: &'a [(usize, f64)],
    /// Absolute load factors in physical order; reversal and unloading allowed.
    pub load_path: &'a [f64],
    /// Caller-owned allocation and geometry limits.
    pub budget: HexAssemblyBudget,
    /// Existing incremental Newton and Krylov controls.
    pub settings: IncrementSettings,
}

/// Accepted endpoint of a hexahedral material-history load path.
#[derive(Debug, Clone, PartialEq)]
pub struct SmallStrainHexSolution<S> {
    /// Nodal displacement, m.
    pub displacement_m: Vec<[f64; 3]>,
    /// Compatible or B-bar assumed strains, in element/Gauss-point order.
    pub quadrature_strain: Vec<Voigt>,
    /// Stress from the same old-history update as the accepted endpoint, Pa.
    pub quadrature_stress_pa: Vec<Voigt>,
    /// Separate committed native history at every integration point.
    pub quadrature_states: Vec<S>,
    /// Reference integration weights in the same order, m^3.
    pub quadrature_volume_m3: Vec<f64>,
    /// Prefix offsets into all quadrature arrays; each element owns eight points.
    pub element_offsets: Vec<usize>,
    /// Accepted increments, support reactions and physical force residuals.
    pub increments: Vec<IncrementRecord>,
}

#[derive(Debug, Clone, Copy)]
struct Point {
    gradients: [[f64; 3]; 8],
    volume: f64,
}

struct Element {
    nodes: [usize; 8],
    points: [Point; 8],
    mean_gradients: [[f64; 3]; 8],
}

fn invalid(what: impl Into<String>) -> SmallStrainTetError {
    SmallStrainTetError::InvalidInput { what: what.into() }
}

fn tensor_norm(tensor: &Voigt) -> f64 {
    tensor.iter().enumerate().fold(0.0_f64, |value, (i, x)| {
        let value = value.hypot(*x);
        if i >= 3 { value.hypot(*x) } else { value }
    })
}

impl Element {
    fn columns(&self, point: &Point, formulation: Formulation) -> [Voigt; 24] {
        std::array::from_fn(|dof| {
            let (node, component) = (dof / 3, dof % 3);
            let mut column = b_column(point.gradients[node], component);
            if formulation == Formulation::BBar {
                let correction =
                    (self.mean_gradients[node][component] - point.gradients[node][component]) / 3.0;
                for value in &mut column[..3] {
                    *value += correction;
                }
            }
            column
        })
    }
}

impl<L: SmallStrainLaw> SmallStrainHexProblem<'_, L> {
    fn prepare(&self, cx: &Cx<'_>) -> Result<Vec<Element>, SmallStrainTetError> {
        checkpoint(cx)?;
        let (nodes, cells) = (self.nodes_m.len(), self.hexahedra.len());
        let dofs = nodes
            .checked_mul(3)
            .ok_or_else(|| invalid("hex DOF count overflow"))?;
        let points = cells
            .checked_mul(8)
            .ok_or_else(|| invalid("hex quadrature count overflow"))?;
        let entries = cells
            .checked_mul(24 * 24)
            .ok_or_else(|| invalid("hex tangent count overflow"))?;
        if nodes == 0 || cells == 0 || self.load_path.is_empty() {
            return Err(invalid("empty hexahedral mesh or load path"));
        }
        for (name, requested, maximum) in [
            ("nodes", nodes, self.budget.maximum_nodes),
            ("elements", cells, self.budget.maximum_elements),
            (
                "quadrature points",
                points,
                self.budget.maximum_quadrature_points,
            ),
            (
                "tangent entries",
                entries,
                self.budget.maximum_tangent_entries,
            ),
            (
                "load steps",
                self.load_path.len(),
                self.budget.maximum_load_steps,
            ),
        ] {
            if requested > maximum {
                return Err(invalid(format!(
                    "hex {name}: requested {requested}, budget {maximum}"
                )));
            }
        }
        if !(self.strain_limit.is_finite() && self.strain_limit > 0.0 && self.strain_limit <= 0.05)
            || !(self.rotation_limit_rad.is_finite()
                && self.rotation_limit_rad > 0.0
                && self.rotation_limit_rad <= 0.05)
            || !(self.budget.minimum_scaled_jacobian.is_finite()
                && self.budget.minimum_scaled_jacobian > 0.0)
            || self.settings.max_newton_iterations == 0
            || !(self.settings.relative_tolerance.is_finite()
                && self.settings.relative_tolerance > 0.0)
            || !(self.settings.linear_tolerance.is_finite() && self.settings.linear_tolerance > 0.0)
        {
            return Err(invalid(
                "inadmissible hex kinematic, geometry or solver limits",
            ));
        }
        self.laws.validate_count(cells)?;
        if self
            .nodes_m
            .iter()
            .flatten()
            .chain(self.load_path)
            .any(|value| !value.is_finite())
        {
            return Err(invalid("nonfinite hex coordinate or load factor"));
        }
        let mut fixed = vec![false; dofs];
        for &(dof, value) in self.prescribed_m {
            if dof >= dofs || !value.is_finite() || fixed[dof] {
                return Err(invalid("invalid or repeated prescribed hex DOF"));
            }
            fixed[dof] = true;
        }
        if dofs - self.prescribed_m.len() > self.budget.maximum_free_dofs {
            return Err(invalid("hex free DOFs exceed budget"));
        }
        if self
            .nodal_forces_n
            .iter()
            .any(|&(dof, value)| dof >= dofs || !value.is_finite())
        {
            return Err(invalid("invalid hex nodal force"));
        }
        let (coordinate, weight) = gauss_legendre(2);
        let mut elements = Vec::with_capacity(cells);
        let mut unique = BTreeSet::new();
        let mut incident = vec![false; nodes];
        for (index, connectivity) in self.hexahedra.iter().enumerate() {
            checkpoint(cx)?;
            let mut key = *connectivity;
            key.sort_unstable();
            if key[7] >= nodes
                || key.windows(2).any(|pair| pair[0] == pair[1])
                || !unique.insert(key)
            {
                return Err(invalid(format!(
                    "hex {index} has missing/repeated nodes or duplicate connectivity"
                )));
            }
            for &node in connectivity {
                incident[node] = true;
            }
            let positions = connectivity.map(|node| self.nodes_m[node]);
            let mut samples = Vec::with_capacity(8);
            let mut volume = 0.0;
            let mut mean_gradients = [[0.0; 3]; 8];
            for i in 0..2 {
                for j in 0..2 {
                    for k in 0..2 {
                        let point = hex8_point(
                            positions,
                            [coordinate[i], coordinate[j], coordinate[k]],
                            self.budget.minimum_scaled_jacobian,
                            index,
                            samples.len(),
                        )
                        .map_err(|error| invalid(error.to_string()))?;
                        let dv = point.determinant_m3 * weight[i] * weight[j] * weight[k];
                        volume += dv;
                        for (mean, gradient) in mean_gradients.iter_mut().zip(point.gradients) {
                            for (value, g) in mean.iter_mut().zip(gradient) {
                                *value = dv.mul_add(g, *value);
                            }
                        }
                        samples.push(Point {
                            gradients: point.gradients,
                            volume: dv,
                        });
                    }
                }
            }
            for (probe, reference) in HEX27_REFERENCE_NODES[..8].iter().enumerate() {
                hex8_point(
                    positions,
                    reference.map(f64::from),
                    self.budget.minimum_scaled_jacobian,
                    index,
                    8 + probe,
                )
                .map_err(|error| invalid(error.to_string()))?;
            }
            hex8_point(
                positions,
                [0.0; 3],
                self.budget.minimum_scaled_jacobian,
                index,
                16,
            )
            .map_err(|error| invalid(error.to_string()))?;
            if !(volume.is_finite() && volume > 0.0) {
                return Err(invalid(format!(
                    "hex {index} has nonfinite integrated volume"
                )));
            }
            for mean in &mut mean_gradients {
                for value in mean {
                    *value /= volume;
                }
            }
            if mean_gradients
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
            {
                return Err(invalid(format!("hex {index} has nonfinite mean gradients")));
            }
            elements.push(Element {
                nodes: *connectivity,
                points: samples.try_into().expect("fixed eight-point Gauss rule"),
                mean_gradients,
            });
        }
        if incident
            .iter()
            .enumerate()
            .any(|(node, used)| !used && (0..3).any(|i| !fixed[3 * node + i]))
        {
            return Err(invalid("unused hex node has free displacement DOFs"));
        }
        Ok(elements)
    }

    fn check(&self, displacement: &[f64], states: &[L::State]) -> Result<(), SmallStrainTetError> {
        if displacement.len() != 3 * self.nodes_m.len()
            || states.len() != 8 * self.hexahedra.len()
            || displacement.iter().any(|value| !value.is_finite())
        {
            return Err(invalid(
                "hex displacement/history length mismatch or nonfinite displacement",
            ));
        }
        Ok(())
    }

    fn strain(
        &self,
        element: &Element,
        point: &Point,
        u: &[f64],
        e: usize,
        q: usize,
    ) -> Result<Voigt, SmallStrainTetError> {
        let mut gradient = [[0.0; 3]; 3];
        for a in 1..8 {
            for (i, row) in gradient.iter_mut().enumerate() {
                let du = u[3 * element.nodes[a] + i] - u[3 * element.nodes[0] + i];
                for (j, value) in row.iter_mut().enumerate() {
                    *value = point.gradients[a][j].mul_add(du, *value);
                }
            }
        }
        let rotation = (0.5 * (gradient[2][1] - gradient[1][2]))
            .hypot(0.5 * (gradient[0][2] - gradient[2][0]))
            .hypot(0.5 * (gradient[1][0] - gradient[0][1]));
        let mut strain = [
            gradient[0][0],
            gradient[1][1],
            gradient[2][2],
            0.5 * (gradient[0][1] + gradient[1][0]),
            0.5 * (gradient[1][2] + gradient[2][1]),
            0.5 * (gradient[2][0] + gradient[0][2]),
        ];
        let raw_norm = tensor_norm(&strain);
        if !rotation.is_finite()
            || rotation > self.rotation_limit_rad
            || !raw_norm.is_finite()
            || raw_norm > self.strain_limit
        {
            return Err(invalid(format!(
                "hex {e} point {q} exceeds raw small-strain/rotation domain"
            )));
        }
        if self.formulation == Formulation::BBar {
            let mut correction = 0.0;
            for a in 1..8 {
                for i in 0..3 {
                    let du = u[3 * element.nodes[a] + i] - u[3 * element.nodes[0] + i];
                    correction = ((element.mean_gradients[a][i] - point.gradients[a][i]) / 3.0)
                        .mul_add(du, correction);
                }
            }
            for value in &mut strain[..3] {
                *value += correction;
            }
            let projected_norm = tensor_norm(&strain);
            if !projected_norm.is_finite() || projected_norm > self.strain_limit {
                return Err(invalid(format!(
                    "hex {e} point {q} exceeds assumed small-strain domain"
                )));
            }
        }
        Ok(strain)
    }

    fn internal(
        &self,
        elements: &[Element],
        states: &[L::State],
        u: &[f64],
        prescribed: Option<&[bool]>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64), SmallStrainTetError> {
        let mut force = vec![0.0; u.len()];
        let mut bounds = prescribed.map(|_| vec![0.0; u.len()]);
        for (e, element) in elements.iter().enumerate() {
            checkpoint(cx)?;
            let law = self.laws.at(e);
            for (q, point) in element.points.iter().enumerate() {
                let strain = self.strain(element, point, u, e, q)?;
                let stress = law.stress(&strain, &states[8 * e + q]);
                if stress.iter().any(|value| !value.is_finite()) {
                    return Err(invalid(format!(
                        "hex {e} point {q} has nonfinite constitutive stress"
                    )));
                }
                let columns = element.columns(point, self.formulation);
                let stress_norm = prescribed.map_or(0.0, |_| tensor_norm(&stress));
                for (local, column) in columns.iter().enumerate() {
                    let dof = 3 * element.nodes[local / 3] + local % 3;
                    let value = (0..6).fold(0.0, |sum, i| {
                        (WEIGHT[i] * column[i]).mul_add(stress[i], sum)
                    });
                    force[dof] = point.volume.mul_add(value, force[dof]);
                    if let (Some(fixed), Some(bounds)) = (prescribed, &mut bounds)
                        && !fixed[dof]
                    {
                        bounds[dof] += point.volume * tensor_norm(column) * stress_norm;
                    }
                }
            }
        }
        let scale = bounds.as_ref().map_or(0.0, |values| {
            values.iter().fold(0.0_f64, |n, value| n.hypot(*value))
        });
        if !scale.is_finite() || force.iter().any(|value| !value.is_finite()) {
            return Err(invalid("nonfinite assembled hex force/scale"));
        }
        Ok((force, scale))
    }

    fn tangent(
        &self,
        elements: &[Element],
        states: &[L::State],
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Csr, SmallStrainTetError> {
        let mut coo = Coo::new(u.len(), u.len());
        for (e, element) in elements.iter().enumerate() {
            checkpoint(cx)?;
            let law = self.laws.at(e);
            let mut local = [[0.0; 24]; 24];
            for (q, point) in element.points.iter().enumerate() {
                let strain = self.strain(element, point, u, e, q)?;
                let tangent = law.tangent(&strain, &states[8 * e + q]);
                if tangent.iter().flatten().any(|value| !value.is_finite()) {
                    return Err(invalid(format!(
                        "hex {e} point {q} has nonfinite constitutive tangent"
                    )));
                }
                let columns = element.columns(point, self.formulation);
                let cb: [Voigt; 24] = columns.map(|column| {
                    std::array::from_fn(|i| {
                        (0..6).fold(0.0, |sum, j| tangent[i][j].mul_add(column[j], sum))
                    })
                });
                for (i, bi) in columns.iter().enumerate() {
                    for (j, cbj) in cb.iter().enumerate() {
                        let value =
                            (0..6).fold(0.0, |sum, v| (WEIGHT[v] * bi[v]).mul_add(cbj[v], sum));
                        local[i][j] = point.volume.mul_add(value, local[i][j]);
                    }
                }
            }
            for (i, row) in local.iter().enumerate() {
                for (j, value) in row.iter().enumerate() {
                    if !value.is_finite() {
                        return Err(invalid("nonfinite local hex tangent"));
                    }
                    if *value != 0.0 {
                        coo.push(
                            3 * element.nodes[i / 3] + i % 3,
                            3 * element.nodes[j / 3] + j % 3,
                            *value,
                        );
                    }
                }
            }
        }
        let matrix = coo.assemble();
        if (0..matrix.nrows()).any(|row| matrix.row(row).1.iter().any(|value| !value.is_finite())) {
            return Err(invalid("nonfinite assembled hex tangent"));
        }
        checkpoint(cx)?;
        Ok(matrix)
    }

    /// Full internal force evaluated from the supplied committed history.
    pub fn internal_force(
        &self,
        displacement_m: &[f64],
        states: &[L::State],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, SmallStrainTetError> {
        let elements = self.prepare(cx)?;
        self.check(displacement_m, states)?;
        self.internal(&elements, states, displacement_m, None, cx)
            .map(|(force, _)| force)
    }

    /// Full consistent tangent, retaining every nonsymmetric material entry.
    pub fn tangent_matrix(
        &self,
        displacement_m: &[f64],
        states: &[L::State],
        cx: &Cx<'_>,
    ) -> Result<Csr, SmallStrainTetError> {
        let elements = self.prepare(cx)?;
        self.check(displacement_m, states)?;
        self.tangent(&elements, states, displacement_m, cx)
    }

    /// Follow the absolute load path from virgin constitutive histories.
    pub fn solve(
        &self,
        cx: &Cx<'_>,
    ) -> Result<SmallStrainHexSolution<L::State>, SmallStrainTetError> {
        self.solve_candidate(None, cx)
    }

    /// Continue from accepted displacement and point histories. Input history
    /// is borrowed; any failed increment leaves the complete caller state
    /// unchanged. Each load factor remains an absolute target.
    pub fn solve_from(
        &self,
        displacement_m: &[f64],
        quadrature_states: &[L::State],
        cx: &Cx<'_>,
    ) -> Result<SmallStrainHexSolution<L::State>, SmallStrainTetError> {
        self.solve_candidate(Some((displacement_m, quadrature_states)), cx)
    }

    fn solve_candidate(
        &self,
        initial: Option<(&[f64], &[L::State])>,
        cx: &Cx<'_>,
    ) -> Result<SmallStrainHexSolution<L::State>, SmallStrainTetError> {
        let elements = self.prepare(cx)?;
        if let Some((u, states)) = initial {
            self.check(u, states)?;
        }
        let operator = HexOperator {
            problem: self,
            elements: &elements,
        };
        let accepted = increment::solve(
            &operator,
            IncrementControls {
                dofs: 3 * self.nodes_m.len(),
                prescribed_m: self.prescribed_m,
                nodal_forces_n: self.nodal_forces_n,
                load_path: self.load_path,
                settings: self.settings,
            },
            initial,
            None,
            cx,
        )?;
        let mut quadrature_strain = Vec::with_capacity(8 * elements.len());
        let mut quadrature_volume_m3 = Vec::with_capacity(8 * elements.len());
        for (e, element) in elements.iter().enumerate() {
            checkpoint(cx)?;
            for (q, point) in element.points.iter().enumerate() {
                quadrature_strain.push(self.strain(
                    element,
                    point,
                    &accepted.displacement_m,
                    e,
                    q,
                )?);
                quadrature_volume_m3.push(point.volume);
            }
        }
        Ok(SmallStrainHexSolution {
            displacement_m: accepted.displacement_m.as_chunks::<3>().0.to_vec(),
            quadrature_strain,
            quadrature_stress_pa: accepted.point_stress,
            quadrature_states: accepted.point_states,
            quadrature_volume_m3,
            element_offsets: (0..=elements.len()).map(|e| 8 * e).collect(),
            increments: accepted.increments,
        })
    }
}

struct HexOperator<'p, 'm, L: SmallStrainLaw> {
    problem: &'p SmallStrainHexProblem<'m, L>,
    elements: &'p [Element],
}

impl<L: SmallStrainLaw> IncrementOperator for HexOperator<'_, '_, L> {
    type State = L::State;

    fn initial_states(&self) -> Vec<Self::State> {
        (0..self.elements.len())
            .flat_map(|e| (0..8).map(move |_| self.problem.laws.at(e).initial_state()))
            .collect()
    }

    fn internal(
        &self,
        states: &[Self::State],
        u: &[f64],
        prescribed: Option<&[bool]>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64), SmallStrainTetError> {
        self.problem
            .internal(self.elements, states, u, prescribed, cx)
    }

    fn tangent(
        &self,
        states: &[Self::State],
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Csr, SmallStrainTetError> {
        self.problem.tangent(self.elements, states, u, cx)
    }

    fn tangent_symmetric(&self) -> bool {
        (0..self.elements.len()).all(|e| self.problem.laws.at(e).admissibility().tangent_symmetric)
    }

    fn commit(
        &self,
        states: &mut [Self::State],
        u: &[f64],
        mut final_stress: Option<&mut Vec<Voigt>>,
        cx: &Cx<'_>,
    ) -> Result<usize, SmallStrainTetError> {
        let mut evolving = 0;
        for (e, element) in self.elements.iter().enumerate() {
            checkpoint(cx)?;
            let law = self.problem.laws.at(e);
            let mut changed = false;
            for (q, point) in element.points.iter().enumerate() {
                let strain = self.problem.strain(element, point, u, e, q)?;
                let old = &states[8 * e + q];
                if let Some(stresses) = &mut final_stress {
                    let stress = law.stress(&strain, old);
                    if stress.iter().any(|value| !value.is_finite()) {
                        return Err(invalid(format!(
                            "hex {e} point {q} refused endpoint stress"
                        )));
                    }
                    stresses.push(stress);
                }
                let next = law.update_state(&strain, old);
                changed |= next != *old;
                states[8 * e + q] = next;
            }
            evolving += usize::from(changed);
        }
        Ok(evolving)
    }

    fn recoverable_trial_error(&self, error: &SmallStrainTetError) -> bool {
        matches!(error, SmallStrainTetError::InvalidInput { .. })
    }
}
