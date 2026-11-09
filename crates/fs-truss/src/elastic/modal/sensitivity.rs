//! Isolated prestressed eigenvalue derivatives in the original truss design.
//!
//! For mass-normalized phi, d(lambda) = phi^T (dK_t - lambda dM) phi.
//! The geometric term depends on the *compatible* member forces: keeping their
//! values frozen while sizing a redundant, heated or settled truss is incorrect.
//! Its displacement covector is B diag(EA/L) c, c_i = (T_i^T phi)^2/L_i.
//! One native elastic adjoint accounts for this preload redistribution.
//!
//! Reference node coordinates may vary with fixed member connectivity, global
//! support flags, force/motion components and material states. These partials
//! do not differentiate eigeniterations, source resolution, topology selection,
//! catalog snapping, density-dependent loads or a self-weight fixed point.

use super::{ModalAssemblyControls, ModalInertia, TrussModalError, TrussModalPencil, assemble};
use crate::elastic::{ElasticControls, ElasticImposedState, MemberTangent, checkpoint, invalid};
use crate::{GroundStructure, LayoutCase, TrussConstructionError};
use fs_exec::Cx;
use fs_modal::{MultiSliceReport, SliceOptions, WindowPlan};

/// One extra compatible solve and its numerical payload allowance.
///
/// Native spectral work, factors and vectors retain the separate WindowPlan
/// admission. This allowance bounds the additional gradient/adjoint arrays,
// not the already prepared pencil, source ground or OS allocator overhead.
#[derive(Debug, Clone, Copy)]
pub struct ModalSensitivityControls {
    /// Tolerance and total iteration allowance of the single elastic adjoint.
    pub elastic: ElasticControls,
    /// Conservative additional numerical payload in bytes.
    pub maximum_bytes: usize,
}

impl Default for ModalSensitivityControls {
    fn default() -> Self {
        Self {
            elastic: ElasticControls::default(),
            maximum_bytes: 512 * 1024 * 1024,
        }
    }
}

/// Fixed-geometry derivatives of signed lambda = omega^2 [s^-2].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MemberModalSensitivity {
    /// d(lambda)/dA [s^-2/m²], including compatible preload and reference mass.
    pub area_derivative: f64,
    /// d(lambda)/dE [s^-2/Pa], with all other source properties frozen.
    pub young_modulus_derivative: f64,
    /// d(lambda)/d(epsilon_0) [s^-2], positive epsilon_0 is free expansion.
    pub stress_free_strain_derivative: f64,
    /// d(lambda)/d(rho_reference) [s^-2/(kg/m³)], with nodal loads fixed.
    pub density_derivative: f64,
}

/// One inertia-isolated native mode and its original-order design covectors.
#[derive(Debug, Clone)]
pub struct ModalSensitivity {
    /// Exactly one window with exactly one mode. Its original-pencil residual,
    /// interval and factor/work diagnostics remain available without promotion.
    pub spectrum: MultiSliceReport,
    /// Signed lambda = omega^2 [s^-2]; negative values remain negative.
    pub eigenvalue_s2: f64,
    /// None means a removed member, with no topology-activation derivative.
    pub members: Vec<Option<MemberModalSensitivity>>,
    /// d(lambda)/df [s^-2/N] in original node/x-y order; zero at supports.
    /// Only variations retaining compatible equilibrium are admissible when
    /// the structure contains unloaded axial mechanisms.
    pub load_gradient: Vec<[f64; 2]>,
    /// d(lambda)/du_prescribed [s^-2/m]; zero on free coordinates.
    pub prescribed_displacement_gradient: Vec<[f64; 2]>,
    /// d(lambda)/d(reference node coordinate) [s^-2/m], in original x-y order.
    /// Both member direction and length vary; A, E, free strain, reference
    /// density, joint masses and GLOBAL force/prescribed-motion components stay
    /// fixed. Supported nodes are geometric variables too; entries are not
    /// masked by mechanical support flags. Only topology/rank-preserving
    /// variations retaining compatible axial equilibrium are admissible.
    pub node_coordinate_gradient: Vec<[f64; 2]>,
    /// d(lambda)/d(isotropic joint mass) [s^-2/kg] in original node order.
    /// Zero-mass joints admit only nonnegative mass perturbations.
    pub joint_mass_derivative: Vec<f64>,
    /// Actual physical free-force residual of the compatible elastic adjoint.
    pub adjoint_relative_equilibrium_residual: f64,
    /// Native CG iterations used by the adjoint, within its separate allowance.
    pub adjoint_iterations: usize,
}

/// Immutable design inputs and one prepared physical vibration pencil.
///
/// Borrowing prevents source arrays or the load case from changing underneath
/// the derivatives. Repeated isolated windows reuse the primal equilibrium and
/// sparse assembly; each derivative call performs exactly one elastic adjoint.
pub struct PreparedModalDesign<'a> {
    ground: &'a GroundStructure,
    case: &'a LayoutCase,
    areas_m2: &'a [f64],
    young_pa: &'a [f64],
    inertia: ModalInertia<'a>,
    pencil: TrussModalPencil,
}

impl<'a> PreparedModalDesign<'a> {
    /// Prepare the original compatible equilibrium and reference-mass pencil.
    ///
    /// The imposed arrays are consumed during preparation. The accepted
    /// response retains their complete elastic-strain/prestress contribution.
    /// # Errors
    /// Inherits the original physical and resource admission of [assemble].
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        ground: &'a GroundStructure,
        case: &'a LayoutCase,
        areas_m2: &'a [f64],
        young_pa: &'a [f64],
        imposed: ElasticImposedState<'_>,
        inertia: ModalInertia<'a>,
        controls: ModalAssemblyControls,
        cx: &Cx<'_>,
    ) -> Result<Self, TrussModalError> {
        let pencil = assemble(
            ground, case, areas_m2, young_pa, imposed, inertia, controls, cx,
        )?;
        Ok(Self {
            ground,
            case,
            areas_m2,
            young_pa,
            inertia,
            pencil,
        })
    }

    /// The same accepted equilibrium, original reference mass and modal DOFs.
    #[must_use]
    pub fn pencil(&self) -> &TrussModalPencil {
        &self.pencil
    }

    /// Differentiate the single eigenvalue isolated by one native window.
    ///
    /// The requested window must contain exactly one inertia-counted mode, and
    /// that mode's original-pencil interval must lie strictly inside it.
    /// Empty windows, repeated/clustered modes and unresolved endpoint overlap
    /// refuse: an individual eigenvalue derivative is not unique at a cluster.
    ///
    /// The one adjoint solves K_e,ff z_f = [B diag(EA/L) c]_f with homogeneous
    /// prescribed motion. For e_i = total strain - free strain, a_i = B_i^T phi,
    /// c_i = (T_i^T phi)^2/L_i, and d_i = B_i^T z, the derivatives include
    /// E_i a_i²/L_i + E_i e_i(c_i-d_i) - lambda rho_i M_i/A_i.
    /// Thus direct stiffness, preload redistribution and consistent mass all
    /// contribute. No global derivative matrix or mode-shape derivative forms.
    ///
    /// Node-coordinate derivatives also differentiate the original member
    /// lengths/directions, compatible preload and full consistent mass. They
    /// hold global load/support-motion components and support flags fixed.
    /// The same adjoint accounts for geometric changes of equilibrium; no
    /// second solve or assembled geometry derivative is needed. Moving a node
    /// changes its reference position, not its prescribed displacement value.
    ///
    /// Loads and support motion may depend on the design: the caller composes
    /// the returned covectors with their own derivatives, including self-weight.
    /// The result is a local derivative of the frozen-reference initial-stress
    /// approximation, not of finite deformation or a code-based stability check.
    /// Spectral/adjoint residuals are numerical diagnostics, not a certified
    /// enclosure on the derivative or a global stability assertion.
    ///
    /// # Errors
    /// Original spectral refusal, non-isolation, exhausted payload/iterations,
    /// unresolved adjoint mechanisms, nonfinite arithmetic and cancellation
    /// publish no partial derivative result.
    #[allow(clippy::too_many_lines)]
    pub fn solve_window(
        &self,
        plan: &WindowPlan,
        options: &SliceOptions,
        controls: ModalSensitivityControls,
        cx: &Cx<'_>,
    ) -> Result<ModalSensitivity, TrussModalError> {
        checkpoint(cx)?;
        if plan.windows.len() != 1 {
            return Err(invalid(
                "truss modal sensitivity window",
                "requires exactly one requested isolated window",
            )
            .into());
        }
        if !(controls.elastic.relative_tolerance > 0.0 && controls.elastic.relative_tolerance < 1.0)
            || controls.elastic.maximum_iterations > 1_000_000
            || controls.maximum_bytes == 0
        {
            return Err(invalid(
                "truss modal sensitivity controls",
                "require tolerance in (0,1), at most 1000000 iterations and positive bytes",
            )
            .into());
        }
        let nodes = self.ground.nodes().len();
        let count = self.ground.members().len();
        // Includes native matrix-free member storage, CG/defect-correction
        // vectors and history, recovery arrays, copied case and final gradients.
        let bytes = 528_u128 * nodes as u128
            + 192_u128 * count as u128
            + 16_u128 * controls.elastic.maximum_iterations as u128
            + 4096;
        if bytes > controls.maximum_bytes as u128 {
            return Err(TrussConstructionError::WorkBudget {
                resource: "truss modal sensitivity additional payload",
                limit: controls.maximum_bytes,
                observed: usize::try_from(bytes).unwrap_or(usize::MAX),
            }
            .into());
        }
        let spectrum = self.pencil.slice_windows(plan, options, cx)?;
        let report = &spectrum.reports[0];
        if report.expected != 1 || report.modes.len() != 1 {
            return Err(invalid(
                "truss modal sensitivity isolation",
                "requires exactly one inertia-counted eigenvalue, not an empty window or cluster",
            )
            .into());
        }
        let mode = &report.modes[0];
        if !(mode.interval.0 > report.window.0 && mode.interval.1 < report.window.1) {
            return Err(invalid(
                "truss modal sensitivity isolation",
                "the original-pencil mode interval must lie strictly inside its window",
            )
            .into());
        }
        // Re-evaluate normalization on the retained original CSR; no assumption
        // that the native eigenvector's rounded modal mass is exactly one.
        let mut modal_mass = 0.0;
        for (row, &phi) in mode.phi.iter().enumerate() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            let (columns, values) = self.pencil.mass.row(row);
            let mut value = 0.0;
            for (index, (&column, &entry)) in columns.iter().zip(values).enumerate() {
                if index % 256 == 0 {
                    checkpoint(cx)?;
                }
                value += entry * mode.phi[column];
            }
            modal_mass += phi * value;
        }
        if !modal_mass.is_finite() || modal_mass <= 0.0 {
            return Err(invalid(
                "truss modal sensitivity normalization",
                "requires finite positive original modal mass",
            )
            .into());
        }
        let scale = modal_mass.sqrt().recip();
        let value = |slot: usize| self.pencil.dof_map[slot].map_or(0.0, |i| mode.phi[i] * scale);
        let mut loads = vec![[0.0; 2]; nodes];
        for index in 0..count {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            if self.areas_m2[index] == 0.0 {
                continue;
            }
            let member = MemberTangent::new(self.ground, index);
            let (_, transverse, _) = mode_work(&member, &value);
            let length = self.ground.lengths()[index];
            let c = transverse * transverse / length;
            let k = self.young_pa[index] * (self.areas_m2[index] / length);
            let force = k * c;
            for local in 0..4 {
                let slot = member.ends[local];
                loads[slot / 2][slot % 2] += member.axial[local] * force;
            }
        }
        checkpoint(cx)?;
        let adjoint_case = LayoutCase::try_new(self.case.supported().to_vec(), loads, nodes)?;
        let adjoint = crate::elastic::solve(
            self.ground,
            &adjoint_case,
            self.areas_m2,
            self.young_pa,
            controls.elastic,
            cx,
        )?;
        let lambda = mode.lambda;
        let mut members = Vec::with_capacity(count);
        let mut node_coordinate_gradient = vec![[0.0; 2]; nodes];
        for index in 0..count {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            let area = self.areas_m2[index];
            if area == 0.0 {
                members.push(None);
                continue;
            }
            let member = MemberTangent::new(self.ground, index);
            let (axial, transverse, inertia_work) = mode_work(&member, &value);
            let length = self.ground.lengths()[index];
            let young = self.young_pa[index];
            let density = self.inertia.member_density_kg_m3[index];
            let c = transverse * transverse / length;
            let d = length * adjoint.member_strain[index];
            let redistribution = c - d;
            let mass_per_area_density = length * inertia_work / 3.0;
            let e = self.pencil.equilibrium.member_elastic_strain[index];
            let derivative = MemberModalSensitivity {
                area_derivative: young * (axial * axial / length)
                    + self.pencil.equilibrium.member_stress_pa[index] * redistribution
                    - lambda * density * mass_per_area_density,
                young_modulus_derivative: area * (axial * axial / length)
                    + area * e * redistribution,
                stress_free_strain_derivative: -(young * area) * redistribution,
                density_derivative: -lambda * area * mass_per_area_density,
            };
            if [
                derivative.area_derivative,
                derivative.young_modulus_derivative,
                derivative.stress_free_strain_derivative,
                derivative.density_derivative,
            ]
            .iter()
            .any(|v| !v.is_finite())
            {
                return Err(invalid("truss modal design derivatives", "must remain finite").into());
            }
            members.push(Some(derivative));

            // With r = x_b-x_a = L n and t perpendicular to n, write
            // dr = n dL + L t dtheta. At fixed u, dN = k(u_t dtheta
            // - epsilon_total dL). The adjoint subtracts
            // z^T d(BN) = d*dN + N*z_t*dtheta, so total strain, not
            // elastic strain, enters the length derivative. The final term
            // is the length derivative of the original consistent mass.
            let (start, end) = self.ground.members()[index];
            let transverse_extension = |field: &[[f64; 2]]| {
                member.transverse[2] * (field[end][0] - field[start][0])
                    + member.transverse[3] * (field[end][1] - field[start][1])
            };
            let u_t = transverse_extension(&self.pencil.equilibrium.displacements_m);
            let z_t = transverse_extension(&adjoint.displacements_m);
            let k = young * (area / length);
            let force = self.pencil.equilibrium.member_force_n[index];
            let radial = -k * (axial * axial / length)
                - redistribution * k * self.pencil.equilibrium.member_strain[index]
                - force * (c / length)
                - lambda * density * area * (inertia_work / 3.0);
            let angular = (2.0 * (k - force / length) * axial * transverse
                + redistribution * k * u_t
                - force * z_t)
                / length;
            for axis in 0..2 {
                let contribution =
                    radial * member.axial[axis + 2] + angular * member.transverse[axis + 2];
                node_coordinate_gradient[start][axis] -= contribution;
                node_coordinate_gradient[end][axis] += contribution;
                if !contribution.is_finite()
                    || !node_coordinate_gradient[start][axis].is_finite()
                    || !node_coordinate_gradient[end][axis].is_finite()
                {
                    return Err(invalid(
                        "truss modal coordinate derivatives",
                        "must remain finite",
                    )
                    .into());
                }
            }
        }
        let mut support_gradient = adjoint.support_reactions_n;
        let mut joint_mass_derivative = Vec::with_capacity(nodes);
        for (node, gradient) in support_gradient.iter_mut().enumerate() {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            for entry in gradient {
                *entry = -*entry;
            }
            let x = value(2 * node);
            let y = value(2 * node + 1);
            let derivative = -lambda * (x * x + y * y);
            if !derivative.is_finite() {
                return Err(
                    invalid("truss modal joint mass derivative", "must remain finite").into(),
                );
            }
            joint_mass_derivative.push(derivative);
        }
        checkpoint(cx)?;
        Ok(ModalSensitivity {
            spectrum,
            eigenvalue_s2: lambda,
            members,
            load_gradient: adjoint.displacements_m,
            prescribed_displacement_gradient: support_gradient,
            node_coordinate_gradient,
            joint_mass_derivative,
            adjoint_relative_equilibrium_residual: adjoint.relative_equilibrium_residual,
            adjoint_iterations: adjoint.iterations,
        })
    }
}

// Consistent bar mass has phi^T M phi = rho A L/3 *
// sum_axis(phi_a² + phi_a phi_b + phi_b²), including transverse motion.
fn mode_work(member: &MemberTangent, value: &impl Fn(usize) -> f64) -> (f64, f64, f64) {
    let a = [value(member.ends[0]), value(member.ends[1])];
    let b = [value(member.ends[2]), value(member.ends[3])];
    let delta = [b[0] - a[0], b[1] - a[1]];
    let axial = member.axial[2] * delta[0] + member.axial[3] * delta[1];
    let transverse = member.transverse[2] * delta[0] + member.transverse[3] * delta[1];
    let inertia =
        (a[0] * a[0] + a[0] * b[0] + b[0] * b[0]) + (a[1] * a[1] + a[1] * b[1] + b[1] * b[1]);
    (axial, transverse, inertia)
}
