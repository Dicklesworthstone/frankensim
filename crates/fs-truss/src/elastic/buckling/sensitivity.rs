//! Total design derivatives of the native fixed/reference buckling family.
//!
//! At one positive root, N*=N_fixed+lambda*N_reference and likewise for u*.
//! For the original mode phi, s=-phi^T K_g(N_reference)phi must be positive.
//! The numerator is d[phi^T(K_e+K_g(N*))phi] at fixed lambda, including the
//! compatible equilibrium derivative through one ordinary elastic adjoint.
//! Divide by s; the fixed/reference load covectors are z/s and lambda*z/s.
//! Reference-density mass has no role in this static buckling pencil.

use super::{BucklingControls, BucklingError, BucklingLoadFamily, TrussBuckling, solve};
use crate::elastic::{ElasticControls, MemberTangent, checkpoint, invalid};
use crate::{GroundStructure, LayoutCase, TrussConstructionError};
use fs_exec::Cx;

/// Additional work admitted for one critical-load derivative query.
#[derive(Debug, Clone, Copy)]
pub struct BucklingSensitivityControls {
    /// Native compatible-adjoint tolerance and total iteration allowance.
    pub elastic: ElasticControls,
    /// Conservative numerical payload beyond the prepared native analysis.
    pub maximum_bytes: usize,
}
impl Default for BucklingSensitivityControls {
    fn default() -> Self {
        Self {
            elastic: ElasticControls::default(),
            maximum_bytes: 512 * 1024 * 1024,
        }
    }
}

/// Fixed-geometry partials for one retained original member.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MemberBucklingSensitivity {
    /// d(lambda)/dA [1/m²].
    pub area_derivative: f64,
    /// d(lambda)/dE [1/Pa], with other material properties frozen.
    pub young_modulus_derivative: f64,
    /// d(lambda)/d(fixed free engineering strain).
    pub stress_free_strain_derivative: f64,
}

/// One numerically separated positive load factor and its design covectors.
#[derive(Debug, Clone)]
pub struct BucklingSensitivity {
    /// Index in the prepared native positive-load list, excluding its guard.
    pub mode_index: usize,
    pub load_factor: f64,
    /// Recomputed -phi^T K_g(N_reference)phi on the original member work.
    pub geometric_modal_work: f64,
    /// Native dimensionless relative residual of this original pencil mode.
    pub relative_eigen_residual: f64,
    /// Minimum relative gap to adjacent RETURNED positive roots, when present.
    /// Numerical separation is not a certificate excluding missed modes.
    pub returned_root_relative_gap: Option<f64>,
    /// Original member order; None denotes removed topology, not activation.
    pub members: Vec<Option<MemberBucklingSensitivity>>,
    /// d(lambda)/d(fixed nodal force) [1/N]; zero at supports.
    pub fixed_load_gradient: Vec<[f64; 2]>,
    /// d(lambda)/d(reference nodal force) [1/N]; zero at supports.
    pub reference_load_gradient: Vec<[f64; 2]>,
    /// d(lambda)/d(fixed prescribed displacement) [1/m]; zero on free slots.
    pub prescribed_displacement_gradient: Vec<[f64; 2]>,
    /// d(lambda)/d(reference node coordinate) [1/m], including supported nodes.
    /// Global force/motion components, A, E, free strain and supports stay fixed.
    pub node_coordinate_gradient: Vec<[f64; 2]>,
    pub adjoint_relative_equilibrium_residual: f64,
    pub adjoint_iterations: usize,
}

/// Immutable source design with the native fixed/reference states and modes.
///
/// Repeated derivative queries reuse both physical equilibria and the native
/// buckling analysis. Each query adds exactly one matrix-free elastic adjoint.
pub struct PreparedBucklingDesign<'a> {
    ground: &'a GroundStructure,
    reference_case: &'a LayoutCase,
    areas_m2: &'a [f64],
    young_pa: &'a [f64],
    analysis: TrussBuckling,
    requested_modes: usize,
}
impl<'a> PreparedBucklingDesign<'a> {
    /// Prepare the native physical family, requesting one extra guard mode.
    ///
    /// controls.modes is the desired derivative range. The additional root,
    /// if returned, is used only to check the upper boundary. Native spectral
    /// dimension/work caps still apply and are never silently expanded.
    /// Fewer returned positive roots do not prove spectrum completeness.
    ///
    /// # Errors
    /// Inherits the native physical, positive-base, work and spectral refusals.
    /// A zero or unrepresentable requested mode count also refuses.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        ground: &'a GroundStructure,
        reference_case: &'a LayoutCase,
        areas_m2: &'a [f64],
        young_pa: &'a [f64],
        family: BucklingLoadFamily<'_>,
        controls: BucklingControls,
        cx: &Cx<'_>,
    ) -> Result<Self, BucklingError> {
        checkpoint(cx)?;
        let requested_modes = controls.modes;
        let guarded = requested_modes
            .checked_add(1)
            .filter(|_| requested_modes > 0)
            .ok_or_else(|| {
                invalid(
                    "truss buckling sensitivity modes",
                    "require positive bounded count",
                )
            })?;
        let analysis = solve(
            ground,
            reference_case,
            areas_m2,
            young_pa,
            family,
            BucklingControls {
                modes: guarded,
                ..controls
            },
            cx,
        )?;
        Ok(Self {
            ground,
            reference_case,
            areas_m2,
            young_pa,
            analysis,
            requested_modes,
        })
    }

    /// Native modes and both compatible states, including any upper guard.
    #[must_use]
    pub fn analysis(&self) -> &TrussBuckling {
        &self.analysis
    }

    /// Available derivative roots; a returned upper guard is never included.
    #[must_use]
    pub fn retained_modes(&self) -> usize {
        self.requested_modes.min(self.analysis.buckling.loads.len())
    }

    /// Differentiate an admitted positive root of the original load family.
    ///
    /// As in the native solid sensitivity owner, adjacent positive returned
    /// roots must have relative separation greater than 1e-6. Both neighbors
    /// are checked; the threshold is also at least ten times the sum of their
    /// dimensionless original-pencil relative residuals. This is a numerical
    /// branch policy, not a certified eigenvalue enclosure or missed-mode test.
    ///
    /// One adjoint solves K_e,ff z_f = [B diag(EA/L)c]_f, c=(T^T phi)^2/L.
    /// At fixed lambda, combined displacement, total/elastic strain and force
    /// are fixed-state plus lambda times reference-state values. The existing
    /// bar geometry contraction differentiates length, direction and preload.
    /// All partials divide by original positive reference geometric modal work.
    ///
    /// Loads, prescribed motion, source state, geometry parameterizations and
    /// self-weight may need further caller-owned chain rules. Only variations
    /// preserving topology, support flags, axial-equilibrium rank and solvability
    /// are admissible. These are local initial-stress derivatives, not finite
    /// deformation, member Euler capacity, postbuckling or nonlinear collapse.
    ///
    /// # Errors
    /// Missing/guard/clustered roots, invalid controls, nonpositive reference
    /// work, failed adjoint, nonfinite derivatives or cancellation publish none.
    #[allow(clippy::too_many_lines)]
    pub fn solve_mode(
        &self,
        mode_index: usize,
        controls: BucklingSensitivityControls,
        cx: &Cx<'_>,
    ) -> Result<BucklingSensitivity, BucklingError> {
        checkpoint(cx)?;
        if !(controls.elastic.relative_tolerance > 0.0 && controls.elastic.relative_tolerance < 1.0)
            || controls.elastic.maximum_iterations > 1_000_000
            || controls.maximum_bytes == 0
        {
            return Err(invalid(
                "truss buckling sensitivity controls",
                "require tolerance in (0,1), at most 1000000 iterations and positive bytes",
            )
            .into());
        }
        if mode_index >= self.retained_modes() {
            return Err(invalid(
                "truss buckling sensitivity root",
                "must name an available derivative root, not the upper guard",
            )
            .into());
        }
        let nodes = self.ground.nodes().len();
        let count = self.ground.members().len();
        // Adjoint member/CG/recovery arrays, case, covectors and bounded history.
        let bytes = 544_u128 * nodes as u128
            + 192_u128 * count as u128
            + 16_u128 * controls.elastic.maximum_iterations as u128
            + 4096;
        if bytes > controls.maximum_bytes as u128 {
            return Err(TrussConstructionError::WorkBudget {
                resource: "truss buckling sensitivity additional payload",
                limit: controls.maximum_bytes,
                observed: usize::try_from(bytes).unwrap_or(usize::MAX),
            }
            .into());
        }
        let buckling = &self.analysis.buckling;
        let lambda = buckling.loads[mode_index];
        let residual = buckling.relative_residuals[mode_index];
        let mut separation: Option<f64> = None;
        for neighbor in [mode_index.checked_sub(1), mode_index.checked_add(1)]
            .into_iter()
            .flatten()
        {
            if let Some(&other) = buckling.loads.get(neighbor) {
                let gap = (other - lambda).abs() / other.max(lambda);
                let required =
                    (10.0 * (residual + buckling.relative_residuals[neighbor])).max(1e-6);
                if !gap.is_finite() || gap <= required {
                    return Err(invalid(
                        "truss buckling sensitivity isolation",
                        "requires separated adjacent returned positive roots including the guard",
                    )
                    .into());
                }
                separation = Some(separation.map_or(gap, |previous| previous.min(gap)));
            }
        }
        let mode = &buckling.modes[mode_index];
        let value = |slot: usize| buckling.dof_map[slot].map_or(0.0, |i| mode[i]);
        let mut loads = vec![[0.0; 2]; nodes];
        let mut denominator = 0.0;
        for index in 0..count {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            if self.areas_m2[index] == 0.0 {
                continue;
            }
            let member = MemberTangent::new(self.ground, index);
            let (_, transverse) = mode_work(&member, &value);
            let length = self.ground.lengths()[index];
            let c = transverse * transverse / length;
            denominator -= c * self.analysis.reference_load.member_force_n[index];
            let force = self.young_pa[index] * (self.areas_m2[index] / length) * c;
            for local in 0..4 {
                let slot = member.ends[local];
                loads[slot / 2][slot % 2] += member.axial[local] * force;
            }
        }
        if !denominator.is_finite() || denominator <= 0.0 {
            return Err(invalid(
                "truss buckling sensitivity normalization",
                "requires positive finite original reference geometric modal work",
            )
            .into());
        }
        checkpoint(cx)?;
        let adjoint_case =
            LayoutCase::try_new(self.reference_case.supported().to_vec(), loads, nodes)?;
        let adjoint = crate::elastic::solve(
            self.ground,
            &adjoint_case,
            self.areas_m2,
            self.young_pa,
            controls.elastic,
            cx,
        )?;
        let baseline = &self.analysis.baseline;
        let reference = &self.analysis.reference_load;
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
            let (axial, transverse) = mode_work(&member, &value);
            let length = self.ground.lengths()[index];
            let young = self.young_pa[index];
            let c = transverse * transverse / length;
            let d = length * adjoint.member_strain[index];
            let redistribution = c - d;
            let elastic_strain = baseline.member_elastic_strain[index]
                + lambda * reference.member_elastic_strain[index];
            let total_strain =
                baseline.member_strain[index] + lambda * reference.member_strain[index];
            let force = baseline.member_force_n[index] + lambda * reference.member_force_n[index];
            let derivative = MemberBucklingSensitivity {
                area_derivative: (young * (axial * axial / length)
                    + young * elastic_strain * redistribution)
                    / denominator,
                young_modulus_derivative: (area * (axial * axial / length)
                    + area * elastic_strain * redistribution)
                    / denominator,
                stress_free_strain_derivative: -(young * area) * redistribution / denominator,
            };
            if [
                derivative.area_derivative,
                derivative.young_modulus_derivative,
                derivative.stress_free_strain_derivative,
            ]
            .iter()
            .any(|v| !v.is_finite())
            {
                return Err(
                    invalid("truss buckling member derivatives", "must remain finite").into(),
                );
            }
            members.push(Some(derivative));
            let (start, end) = self.ground.members()[index];
            let transverse_extension = |field: &[[f64; 2]]| {
                member.transverse[2] * (field[end][0] - field[start][0])
                    + member.transverse[3] * (field[end][1] - field[start][1])
            };
            let u_t = transverse_extension(&baseline.displacements_m)
                + lambda * transverse_extension(&reference.displacements_m);
            let z_t = transverse_extension(&adjoint.displacements_m);
            let [radial, angular] = MemberTangent::coordinate_work(
                length,
                young * (area / length),
                force,
                [axial, transverse],
                redistribution,
                total_strain,
                u_t,
                z_t,
            );
            for axis in 0..2 {
                let derivative = (radial * member.axial[axis + 2]
                    + angular * member.transverse[axis + 2])
                    / denominator;
                node_coordinate_gradient[start][axis] -= derivative;
                node_coordinate_gradient[end][axis] += derivative;
                if !derivative.is_finite()
                    || !node_coordinate_gradient[start][axis].is_finite()
                    || !node_coordinate_gradient[end][axis].is_finite()
                {
                    return Err(invalid(
                        "truss buckling coordinate derivatives",
                        "must remain finite",
                    )
                    .into());
                }
            }
        }
        let mut fixed_load_gradient = adjoint.displacements_m;
        let mut support_gradient = adjoint.support_reactions_n;
        let mut reference_load_gradient = vec![[0.0; 2]; nodes];
        for node in 0..nodes {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            for axis in 0..2 {
                fixed_load_gradient[node][axis] /= denominator;
                reference_load_gradient[node][axis] = lambda * fixed_load_gradient[node][axis];
                support_gradient[node][axis] = -support_gradient[node][axis] / denominator;
                if !fixed_load_gradient[node][axis].is_finite()
                    || !reference_load_gradient[node][axis].is_finite()
                    || !support_gradient[node][axis].is_finite()
                {
                    return Err(invalid(
                        "truss buckling load/support derivatives",
                        "must remain finite",
                    )
                    .into());
                }
            }
        }
        checkpoint(cx)?;
        Ok(BucklingSensitivity {
            mode_index,
            load_factor: lambda,
            geometric_modal_work: denominator,
            relative_eigen_residual: residual,
            returned_root_relative_gap: separation,
            members,
            fixed_load_gradient,
            reference_load_gradient,
            prescribed_displacement_gradient: support_gradient,
            node_coordinate_gradient,
            adjoint_relative_equilibrium_residual: adjoint.relative_equilibrium_residual,
            adjoint_iterations: adjoint.iterations,
        })
    }
}

fn mode_work(member: &MemberTangent, value: &impl Fn(usize) -> f64) -> (f64, f64) {
    let dx = value(member.ends[2]) - value(member.ends[0]);
    let dy = value(member.ends[3]) - value(member.ends[1]);
    (
        member.axial[2] * dx + member.axial[3] * dy,
        member.transverse[2] * dx + member.transverse[3] * dy,
    )
}
