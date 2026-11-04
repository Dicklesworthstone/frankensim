//! Held uniform pressure on the entire boundary of a Tet10 material domain.
//!
//! Positive q is compressive: the external force is -q grad(V), with loading
//! potential q V and contribution +q Hessian(V) to the equilibrium Hessian.
//! Volume, force and tangent use the SAME original P2 map and quadrature. No
//! boundary flattening, frozen normal or separately added prestress matrix is
//! used. Internal material pressure p retains its original tensile convention;
//! q is never substituted into the material's bulk dual or source conditions.
//!
//! The load acts on every exterior face of every original component, including
//! cavity surfaces with the material-domain orientation. It is not pressure on
//! an arbitrarily selected patch or a gas-cavity constitutive law. Conforming
//! shared midsides and a closed oriented face/edge incidence are admitted;
//! sampled positive Jacobians do not prove geometric injectivity. Native
//! geometry/quadrature, load stepping, cancellation and solver budgets remain
//! authoritative. This adapter covers equilibrium and held-load vibration,
//! not a new time integrator.

use std::collections::BTreeMap;

use super::vibration::incompressible::IncompressibleMixedVibration;
use super::vibration::{MixedUpVibration, MixedVibrationError, MixedVibrationSettings};
use super::{
    MixedUpError, MixedUpEvaluation, MixedUpProblem, MixedUpSolution, PreparedMixedUpProblem,
    checkpoint, invalid, point_deformation,
};
use fs_exec::Cx;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_material::hyper::mixed::MixedFiniteStrainPotential;

/// One pressure on the complete material-domain boundary. Positive means
/// compression, negative means outward traction; the value is held in modes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniformBoundaryPressure {
    pressure_pa: f64,
}

impl UniformBoundaryPressure {
    pub fn new(pressure_pa: f64) -> Result<Self, MixedUpError> {
        if !pressure_pa.is_finite() {
            return Err(invalid("uniform boundary pressure must be finite"));
        }
        Ok(Self { pressure_pa })
    }

    pub const fn pressure_pa(self) -> f64 {
        self.pressure_pa
    }
}

/// The original prepared mixed problem with a separately identified external
/// pressure. All native material, support, pressure-gauge and mass semantics
/// are retained; the prepared owner is not exposed as an unloaded problem.
pub struct UniformPressureMixedProblem<'a, L: MixedFiniteStrainPotential = Hyperelastic> {
    prepared: PreparedMixedUpProblem<'a, L>,
}

/// Immutable pressure-loaded equilibrium with separate material and load
/// potentials. The original source/problem remains available for inspection.
#[derive(Debug)]
pub struct MixedPressureSolution<'a, L: MixedFiniteStrainPotential = Hyperelastic> {
    problem: MixedUpProblem<'a, L>,
    boundary_pressure: UniformBoundaryPressure,
    equilibrium: MixedUpSolution,
    volume_m3: f64,
    loading_potential_j: f64,
}

impl<'a, L: MixedFiniteStrainPotential> MixedUpProblem<'a, L> {
    /// Prepare a uniform pressure on the WHOLE material-domain boundary. Its
    /// magnitude ramps with the original support/dead-load continuation factor.
    /// This does not change the constitutive law or prescribe internal p.
    pub fn prepare_with_uniform_boundary_pressure(
        &self,
        pressure_pa: f64,
        cx: &Cx<'_>,
    ) -> Result<UniformPressureMixedProblem<'a, L>, MixedUpError> {
        checkpoint(cx)?;
        let pressure = UniformBoundaryPressure::new(pressure_pa)?;
        let mut prepared = self.prepare(cx)?;
        admit_closed_domain(self, cx)?;
        prepared.boundary_pressure = Some(pressure);
        Ok(UniformPressureMixedProblem { prepared })
    }
}

impl<'a, L: MixedFiniteStrainPotential> UniformPressureMixedProblem<'a, L> {
    pub const fn problem(&self) -> &MixedUpProblem<'a, L> {
        &self.prepared.problem
    }

    pub fn boundary_pressure(&self) -> UniformBoundaryPressure {
        self.prepared
            .boundary_pressure
            .expect("bound pressure load")
    }

    pub fn pressure_nodes(&self) -> &[usize] {
        self.prepared.pressure_nodes()
    }

    pub fn reference_volume_m3(&self) -> f64 {
        self.prepared.reference_volume_m3()
    }

    /// Volume of the deformed material domain using its original isoparametric
    /// volume rule. Material/source admission remains with evaluate/solve.
    pub fn volume_m3(&self, displacement_m: &[[f64; 3]], cx: &Cx<'_>) -> Result<f64, MixedUpError> {
        checkpoint(cx)?;
        if displacement_m.len() != self.prepared.mesh.nodes_count
            || displacement_m.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid("invalid mixed pressure volume displacement"));
        }
        let mut volume = 0.0;
        for element in &self.prepared.mesh.elements {
            for point in &element.points {
                checkpoint(cx)?;
                let f = point_deformation(element, point, displacement_m);
                let cof = cofactor(&f);
                let j = (0..3).map(|i| f[i] * cof[i]).sum::<f64>();
                if !j.is_finite() || j <= 0.0 {
                    return Err(invalid(
                        "mixed pressure volume requires positive sampled det(F)",
                    ));
                }
                volume += point.weight_m3 * j;
            }
        }
        if !volume.is_finite() || volume <= 0.0 {
            return Err(invalid("mixed pressure volume is unrepresentable"));
        }
        Ok(volume)
    }

    /// Full-load residual and Hessian. `internal_forces_n` is material force
    /// MINUS external pressure force; `potential_j` remains material-only.
    /// Pressure residual/coupling/compliance retain their original meaning.
    pub fn evaluate(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: &[f64],
        with_tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<MixedUpEvaluation, MixedUpError> {
        self.prepared
            .evaluate(displacement_m, pressure_pa, with_tangent, cx)
    }

    pub fn solve(&self, cx: &Cx<'_>) -> Result<MixedPressureSolution<'a, L>, MixedUpError> {
        self.retain(self.prepared.solve(cx)?, cx)
    }

    /// Continue from an iterate, enforcing the original absolute support and
    /// pressure-gauge targets at every native load step.
    pub fn solve_from(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: &[f64],
        cx: &Cx<'_>,
    ) -> Result<MixedPressureSolution<'a, L>, MixedUpError> {
        self.retain(
            self.prepared.solve_from(displacement_m, pressure_pa, cx)?,
            cx,
        )
    }

    fn retain(
        &self,
        equilibrium: MixedUpSolution,
        cx: &Cx<'_>,
    ) -> Result<MixedPressureSolution<'a, L>, MixedUpError> {
        let volume_m3 = self.volume_m3(&equilibrium.displacement_m, cx)?;
        let boundary_pressure = self.boundary_pressure();
        let loading_potential_j = boundary_pressure.pressure_pa() * volume_m3;
        if !loading_potential_j.is_finite() {
            return Err(invalid(
                "mixed pressure loading potential is unrepresentable",
            ));
        }
        checkpoint(cx)?;
        Ok(MixedPressureSolution {
            problem: self.prepared.problem,
            boundary_pressure,
            equilibrium,
            volume_m3,
            loading_potential_j,
        })
    }
}

impl<L: MixedFiniteStrainPotential> MixedPressureSolution<'_, L> {
    pub const fn problem(&self) -> &MixedUpProblem<'_, L> {
        &self.problem
    }
    pub const fn boundary_pressure(&self) -> UniformBoundaryPressure {
        self.boundary_pressure
    }
    /// The native accepted equilibrium. Its potential_j is the original
    /// material saddle potential, with neither q V nor dead-force work added.
    pub const fn equilibrium(&self) -> &MixedUpSolution {
        &self.equilibrium
    }
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    /// External loading potential q V, separate from material stored energy.
    pub const fn loading_potential_j(&self) -> f64 {
        self.loading_potential_j
    }
}

impl<'a, L: HyperelasticMaterial + MixedFiniteStrainPotential> UniformPressureMixedProblem<'a, L> {
    /// Original finite-compliance modes with the pressure held at its full
    /// value. The SAME prepared owner equilibrates and assembles the tangent.
    pub fn solve_vibration(
        &self,
        density: &[f64],
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<MixedUpVibration<'a, L>, MixedVibrationError> {
        self.prepared
            .solve_vibration_impl(None, density, settings, cx)
    }

    pub fn solve_vibration_from(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: &[f64],
        density: &[f64],
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<MixedUpVibration<'a, L>, MixedVibrationError> {
        self.prepared.solve_vibration_impl(
            Some((displacement_m, pressure_pa)),
            density,
            settings,
            cx,
        )
    }

    /// Original exact-constraint modes, including its declared pressure gauge
    /// and full original weak-volume residual checks after mode lifting.
    pub fn solve_incompressible_vibration(
        &self,
        density: &[f64],
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<IncompressibleMixedVibration<'a, L>, MixedVibrationError> {
        self.prepared
            .solve_incompressible_vibration_impl(None, density, settings, cx)
    }

    pub fn solve_incompressible_vibration_from(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: &[f64],
        density: &[f64],
        settings: MixedVibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<IncompressibleMixedVibration<'a, L>, MixedVibrationError> {
        self.prepared.solve_incompressible_vibration_impl(
            Some((displacement_m, pressure_pa)),
            density,
            settings,
            cx,
        )
    }
}

fn admit_closed_domain<L: MixedFiniteStrainPotential>(
    problem: &MixedUpProblem<'_, L>,
    cx: &Cx<'_>,
) -> Result<(), MixedUpError> {
    // Four face entries and at most twelve exterior edge entries per Tet10.
    // Both temporary maps fit under the existing assembly payload budget.
    let bytes = problem
        .elements
        .len()
        .checked_mul(16 * 128)
        .ok_or_else(|| invalid("pressure boundary incidence size overflow"))?;
    if bytes > problem.settings.maximum_assembly_bytes {
        return Err(invalid(
            "pressure boundary incidence exceeds mixed assembly budget",
        ));
    }
    let mut faces = BTreeMap::<[usize; 3], (u8, i8, [usize; 3])>::new();
    for tet in problem.elements {
        checkpoint(cx)?;
        // Outward local orientation of the admitted positive Tet10 map.
        for local in [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]] {
            let face = local.map(|a| tet[a]);
            let mut key = face;
            key.sort_unstable();
            let inversions = usize::from(face[0] > face[1])
                + usize::from(face[0] > face[2])
                + usize::from(face[1] > face[2]);
            let orientation = if inversions % 2 == 0 { 1 } else { -1 };
            let entry = faces.entry(key).or_insert((0, 0, face));
            entry.0 += 1;
            entry.1 += orientation;
            if entry.0 > 2 || (entry.0 == 2 && entry.1 != 0) {
                return Err(invalid(
                    "pressure domain has nonmanifold or inconsistently oriented faces",
                ));
            }
        }
    }
    let mut edges = BTreeMap::<[usize; 2], (u8, i8)>::new();
    for &(count, _, face) in faces.values() {
        checkpoint(cx)?;
        if count != 1 {
            continue;
        }
        for (a, b) in [(face[0], face[1]), (face[1], face[2]), (face[2], face[0])] {
            let (key, orientation) = if a < b { ([a, b], 1) } else { ([b, a], -1) };
            let entry = edges.entry(key).or_insert((0, 0));
            entry.0 += 1;
            entry.1 += orientation;
            if entry.0 > 2 {
                return Err(invalid("pressure exterior has a nonmanifold edge"));
            }
        }
    }
    if edges.is_empty() {
        return Err(invalid("pressure domain requires a nonempty exterior"));
    }
    for &(count, orientation) in edges.values() {
        checkpoint(cx)?;
        if count != 2 || orientation != 0 {
            return Err(invalid(
                "pressure exterior is not closed and consistently oriented",
            ));
        }
    }
    Ok(())
}

fn cofactor(f: &[f64; 9]) -> [f64; 9] {
    [
        f[4] * f[8] - f[5] * f[7],
        f[5] * f[6] - f[3] * f[8],
        f[3] * f[7] - f[4] * f[6],
        f[2] * f[7] - f[1] * f[8],
        f[0] * f[8] - f[2] * f[6],
        f[1] * f[6] - f[0] * f[7],
        f[1] * f[5] - f[2] * f[4],
        f[2] * f[3] - f[0] * f[5],
        f[0] * f[4] - f[1] * f[3],
    ]
}

pub(super) fn add_volume_gradient(f: &[f64; 9], q: f64, stress: &mut [f64; 9]) {
    if q != 0.0 {
        for (value, cof) in stress.iter_mut().zip(cofactor(f)) {
            *value += q * cof;
        }
    }
}

pub(super) fn add_volume_hessian(f: &[f64; 9], q: f64, tangent: &mut [[f64; 9]; 9]) {
    if q == 0.0 {
        return;
    }
    // d²J/dF_iJ dF_kL = epsilon_ikm epsilon_JLn F_mn. This is
    // geometric loading only: original up/pu/pp blocks remain untouched.
    for (ij, row) in tangent.iter_mut().enumerate() {
        let (i, j) = (ij / 3, ij % 3);
        for (kl, value) in row.iter_mut().enumerate() {
            let (k, l) = (kl / 3, kl % 3);
            if i != k && j != l {
                let sign_i = if (i + 1) % 3 == k { 1.0 } else { -1.0 };
                let sign_j = if (j + 1) % 3 == l { 1.0 } else { -1.0 };
                *value += q * sign_i * sign_j * f[3 * (3 - i - k) + (3 - j - l)];
            }
        }
    }
}
