//! Geometrically exact Cosserat rods (bead tfz.14): Lie-group nodal
//! state — positions in R³ plus canonical `fs-ga::So3` nodal frames updated
//! MULTIPLICATIVELY on the group (never additive quaternion arithmetic) —
//! with the full strain set: axial/shear
//! `Γ = Rᵀ r′ − e₁` and bending/torsion `κ = Log(Rᵢ⁻¹Rᵢ₊₁)/L₀`,
//! both built from RELATIVE quantities so rigid motions produce
//! exactly zero strain (the objectivity battery checks the energy is
//! invariant, not just small).
//!
//! Statics: finite-difference internal residual with exact dead tip loads
//! (body/right rotational perturbations — `R <- R Exp(delta)`, never
//! component nudging) and FD
//! tangent, dense-LU Newton with load stepping. Fixture-scale by
//! design (≤ a few hundred DOFs); analytic tangents and SE(3)
//! DYNAMICS under fs-time's symplectic integrators are the recorded
//! successor scope. [`MaterialRod`] binds the shared source-resolved concentric
//! section to this same energy and equilibrium solver, retaining all layers.

use crate::SolidError;
use fs_ga::{GaError, So3, Vec3};

/// Spatial material assignments using this same rod kinematics/Newton owner.
pub mod material_field;

/// Diagonal section stiffness of a Cosserat rod.
#[derive(Debug, Clone, Copy)]
pub struct RodSection {
    /// Axial stiffness EA.
    pub ea: f64,
    /// Shear stiffness GA (both transverse directions).
    pub ga: f64,
    /// Torsional stiffness GJ.
    pub gj: f64,
    /// Bending stiffness EI (both bending directions).
    pub ei: f64,
}

/// A discrete rod: reference = straight along +x with uniform segment
/// length; state = nodal positions plus canonical `SO(3)` frames (body→world).
#[derive(Debug, Clone)]
pub struct Rod {
    /// Nodal positions.
    pub positions: Vec<[f64; 3]>,
    /// Validated, canonical nodal frames mapping body vectors into the world.
    pub frames: Vec<So3>,
    /// Reference segment length.
    pub l0: f64,
    /// Section stiffness.
    pub section: RodSection,
}

/// Dead end-loading at the rod tip.
#[derive(Debug, Clone, Copy, Default)]
pub struct TipLoad {
    /// Dead force on the last node (world frame).
    pub force: [f64; 3],
    /// Dead moment on the last node (body frame of the tip).
    pub moment: [f64; 3],
}

fn group_error(error: GaError) -> SolidError {
    SolidError::InternalInvariant {
        what: format!("canonical rod-frame operation refused: {error}"),
    }
}

// Translation probes have length units; rotation probes are angles. Material
// rods supply the segment length, avoiding an artificial 3% strain probe on
// a 1 cm segment. Legacy raw rods retain their original 1 metre probe scale.
fn fd_step(coordinate: usize, translation_scale_m: f64, relative_step: f64) -> Result<f64, SolidError> {
    let step = if coordinate % 6 < 3 { relative_step * translation_scale_m } else { relative_step };
    if !step.is_finite() || step <= 0.0 {
        return Err(SolidError::InvalidInput { what: "rod derivative probe is unrepresentable".into() });
    }
    Ok(step)
}

impl Rod {
    /// A straight reference rod along +x with `segments` segments of
    /// total length `length`.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn straight(length: f64, segments: usize, section: RodSection) -> Rod {
        let l0 = length / segments as f64;
        Rod {
            positions: (0..=segments).map(|i| [i as f64 * l0, 0.0, 0.0]).collect(),
            frames: vec![So3::identity(); segments + 1],
            l0,
            section,
        }
    }

    /// Segment strains: (Γ − e₁-relative, κ), both in the segment
    /// frame.
    /// # Errors
    /// Returns an internal-invariant refusal if a canonical group operation
    /// detects accumulated non-unit drift.
    pub fn strains(&self, seg: usize) -> Result<([f64; 3], [f64; 3]), SolidError> {
        let (a, b) = (seg, seg + 1);
        let relative = self.frames[b]
            .body_minus(self.frames[a])
            .map_err(group_error)?;
        let midpoint = self.frames[a]
            .body_plus(relative.scale(0.5))
            .map_err(group_error)?;
        let dr = Vec3::new(
            (self.positions[b][0] - self.positions[a][0]) / self.l0,
            (self.positions[b][1] - self.positions[a][1]) / self.l0,
            (self.positions[b][2] - self.positions[a][2]) / self.l0,
        );
        // Γ = Rᵀ r′ − e₁ (rotate world tangent into the mid frame).
        let local = midpoint.inverse().rotate(dr).map_err(group_error)?;
        let gamma = [local.x - 1.0, local.y, local.z];
        let angular = relative.angular;
        let kappa = [
            angular.x / self.l0,
            angular.y / self.l0,
            angular.z / self.l0,
        ];
        Ok((gamma, kappa))
    }

    /// Total internal strain energy.
    /// # Errors
    /// Propagates canonical frame-operation refusals from [`Rod::strains`].
    pub fn energy(&self) -> Result<f64, SolidError> {
        let s = &self.section;
        let mut e = 0.0;
        for seg in 0..self.positions.len() - 1 {
            let (g, k) = self.strains(seg)?;
            e += 0.5
                * self.l0
                * s.twice_energy_per_length(g, k);
        }
        Ok(e)
    }

    /// Free DOFs: everything except node 0 (clamped position + frame);
    /// layout per free node: [dx, dy, dz, θx, θy, θz].
    fn ndof(&self) -> usize {
        6 * (self.positions.len() - 1)
    }

    fn apply_increment(&mut self, delta: &[f64], scale: f64) -> Result<(), SolidError> {
        for node in 1..self.positions.len() {
            let k = 6 * (node - 1);
            for c in 0..3 {
                self.positions[node][c] += scale * delta[k + c];
            }
            let increment = Vec3::new(
                scale * delta[k + 3],
                scale * delta[k + 4],
                scale * delta[k + 5],
            );
            self.frames[node] = self.frames[node]
                .body_plus(fs_ga::So3Tangent::new(increment))
                .map_err(group_error)?;
        }
        Ok(())
    }

    /// FD internal residual in body/right perturbation coordinates, with
    /// exact subtraction of the applied tip force and body-frame tip moment.
    /// The energy callback may also contain distributed-load potentials.
    fn residual<E: From<SolidError>>(
        &self, load: &TipLoad, factor: f64,
        energy: &impl Fn(&Self) -> Result<f64, E>,
        checkpoint: &impl Fn() -> Result<(), E>,
        translation_scale_m: f64,
    ) -> Result<Vec<f64>, E> {
        let n = self.ndof();
        let mut r = vec![0.0f64; n];
        // FD scales: the tangent is FD-of-FD — the residual step must
        // sit well above energy roundoff or the nested difference is
        // noise (measured: a junk Newton direction that creeps).
        let mut probe = self.clone();
        for k in 0..n {
            checkpoint()?;
            let eps = fd_step(k, translation_scale_m, 3e-6)?;
            let mut d = vec![0.0f64; n];
            d[k] = eps;
            probe.clone_from(self);
            probe.apply_increment(&d, 1.0)?;
            let ep = energy(&probe)?;
            probe.clone_from(self);
            probe.apply_increment(&d, -1.0)?;
            let em = energy(&probe)?;
            r[k] = (ep - em) / (2.0 * eps);
        }
        // Differencing the external potential against a much larger stored
        // energy can erase a real load. These derivatives are known exactly.
        let tipk = n - 6;
        for c in 0..3 {
            r[tipk + c] -= factor * load.force[c];
            r[tipk + 3 + c] -= factor * load.moment[c];
        }
        Ok(r)
    }

    /// Newton statics under load stepping; returns residual norms per
    /// step (evidence).
    ///
    /// # Errors
    /// [`SolidError::NewtonStalled`] with the history on failure.
    pub fn solve_static(
        &mut self,
        load: &TipLoad,
        steps: usize,
        tol: f64,
    ) -> Result<Vec<f64>, SolidError> {
        self.solve_static_with_energy(load, steps, tol, &Self::energy, &|| Ok(()), None, 1.0)
    }

    // Both uniform and spatial material paths use the SAME nodal update,
    // residual, tangent, load stepping and line search. Additional supports
    // project the same residual onto their free coordinates. No solver fork.
    #[allow(clippy::too_many_lines, clippy::too_many_arguments)] // one Newton/line-search owner
    fn solve_static_with_energy<E: From<SolidError>>(
        &mut self, load: &TipLoad, steps: usize, tol: f64,
        energy: &impl Fn(&Self) -> Result<f64, E>,
        checkpoint: &impl Fn() -> Result<(), E>,
        free_dofs: Option<&[usize]>,
        translation_scale_m: f64,
    ) -> Result<Vec<f64>, E> {
        checkpoint()?;
        let full_n = self.ndof();
        // These indices are admitted by the spatial support owner. The legacy
        // path uses the full identity projection and keeps its operation order.
        let free = match free_dofs {
            Some(dofs) => dofs.to_vec(),
            None => (0..full_n).collect(),
        };
        let n = free.len();
        let mut fixed = vec![true; full_n];
        for &k in &free { fixed[k] = false; }
        let anchored_positions = self.positions.clone();
        let anchored_frames = self.frames.clone();
        // Even a zero group increment may canonicalize a frame's bits. Restore
        // fixed coordinates explicitly rather than letting a clamp drift.
        let constrain = |rod: &mut Self| {
            for node in 1..rod.positions.len() {
                let k = 6 * (node - 1);
                for c in 0..3 {
                    if fixed[k + c] { rod.positions[node][c] = anchored_positions[node][c]; }
                }
                if fixed[k + 3] && fixed[k + 4] && fixed[k + 5] {
                    rod.frames[node] = anchored_frames[node];
                }
            }
        };
        if n == 0 {
            energy(self)?;
            checkpoint()?;
            return Ok(vec![0.0; steps]);
        }
        let residual_norm = |r: &[f64]| free.iter().map(|&k| r[k] * r[k]).sum::<f64>().sqrt();
        let mut finals = Vec::new();
        for step in 1..=steps {
            #[allow(clippy::cast_precision_loss)]
            let factor = step as f64 / steps as f64;
            let mut history = Vec::new();
            let mut converged = false;
            for _ in 0..40 {
                checkpoint()?;
                let r = self.residual(load, factor, energy, checkpoint, translation_scale_m)?;
                let rn = residual_norm(&r);
                history.push(rn);
                if rn < tol {
                    converged = true;
                    break;
                }
                // FD tangent (fixture-scale dense).
                let mut kmat = vec![0.0f64; n * n];
                let mut probe = self.clone();
                for col in 0..n {
                    checkpoint()?;
                    let eps = fd_step(free[col], translation_scale_m, 3e-4)?;
                    let mut d = vec![0.0f64; full_n];
                    d[free[col]] = eps;
                    probe.clone_from(self);
                    probe.apply_increment(&d, 1.0)?;
                    constrain(&mut probe);
                    let rp = probe.residual(load, factor, energy, checkpoint, translation_scale_m)?;
                    probe.clone_from(self);
                    probe.apply_increment(&d, -1.0)?;
                    constrain(&mut probe);
                    let rm = probe.residual(load, factor, energy, checkpoint, translation_scale_m)?;
                    for row in 0..n {
                        kmat[row * n + col] = (rp[free[row]] - rm[free[row]]) / (2.0 * eps);
                    }
                }
                checkpoint()?;
                let f = fs_la::factor::lu(&kmat, n).map_err(|_| SolidError::SolveFailed {
                    iters: 0,
                    rel_residual: f64::INFINITY,
                })?;
                let mut reduced: Vec<f64> = free.iter().map(|&k| -r[k]).collect();
                f.solve(&mut reduced);
                let mut d = vec![0.0; full_n];
                for (&k, &value) in free.iter().zip(&reduced) { d[k] = value; }
                checkpoint()?;
                // Backtracking on the residual norm (dead moments make
                // Π alone an incomplete merit function).
                let mut alpha = 1.0f64;
                let mut accepted = false;
                for _ in 0..20 {
                    checkpoint()?;
                    let mut trial = self.clone();
                    trial.apply_increment(&d, alpha)?;
                    constrain(&mut trial);
                    let rt = trial.residual(load, factor, energy, checkpoint, translation_scale_m)?;
                    let rtn = residual_norm(&rt);
                    if rtn < rn {
                        *self = trial;
                        accepted = true;
                        break;
                    }
                    alpha *= 0.5;
                }
                if !accepted {
                    return Err(SolidError::NewtonStalled { history }.into());
                }
            }
            if !converged {
                return Err(SolidError::NewtonStalled { history }.into());
            }
            finals.push(*history.last().expect("nonempty"));
        }
        Ok(finals)
    }
}

impl RodSection {
    // Keep the original operation tree for the uniform path. The spatial
    // path changes coefficients, not the constitutive equation or kinematics.
    fn twice_energy_per_length(&self, g: [f64; 3], k: [f64; 3]) -> f64 {
        self.ea * g[0] * g[0]
            + self.ga * (g[1] * g[1] + g[2] * g[2])
            + self.gj * k[0] * k[0]
            + self.ei * (k[1] * k[1] + k[2] * k[2])
    }

    /// Project the canonical layered material section into this rod's existing
    /// diagonal energy law. EA, EI, corrected GA and circular GJ are all taken
    /// from the same section; no isotropic preset or effective radius is inferred.
    /// The returned numerical coefficients alone do not retain source authority;
    /// use [`MaterialRod`] to keep the material bundles with the solved geometry.
    ///
    /// # Errors
    /// Refuses a sampled frequency/other-state curve as a static elastic law.
    /// Coefficients must be constants or T-only curves at fixed positive absolute
    /// T, and any explicit frequency support must include zero (static loading).
    pub fn try_from_concentric_material(
        section: &fs_material::section::ConcentricShearSection,
    ) -> Result<Self, SolidError> {
        use fs_material::state_point::{DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY};
        for (index, layer) in section.elastic().layers().iter().enumerate() {
            for name in [DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY, section.shear_property()] {
                let refuse = || SolidError::MaterialRefused {
                    what: format!("rod layer {index} property {name} must support a frozen static elastic law"),
                };
                let property = layer.material.property(name).ok_or_else(refuse)?;
                if !property.is_constant_at_fixed_temperature() {
                    return Err(refuse());
                }
                if let Some((lo, hi)) = property.angular_frequency_band_rad_s().map_err(|_| refuse())?
                    && (lo > 0.0 || hi < 0.0)
                {
                    return Err(refuse());
                }
            }
        }
        Ok(Self {
            ea: section.elastic().axial_stiffness_n(),
            ga: section.shear_stiffness_n(),
            gj: section.torsional_stiffness_n_m2(),
            ei: section.elastic().bending_stiffness_n_m2(),
        })
    }
}

/// The retained dense/finite-difference rod solver is a fixture-scale solver.
/// This material entry point admits at most 384 free degrees of freedom.
pub const MAX_MATERIAL_ROD_SEGMENTS: usize = 64;

/// The existing Cosserat rod and its immutable, source-bound layered section.
///
/// Each layer describes the stress-free reference configuration at its retained
/// query state. Density supplies reference mass; it does not silently add gravity
/// or dynamics. This is the existing diagonal linear elastic rod law with exact
/// rotational kinematics, not helical winding, interlayer slip, plasticity,
/// thermal evolution or a validated 3D composite material. Replacing a material
/// constructs a new specimen; no old stress/history is transplanted.
#[derive(Clone, Debug)]
pub struct MaterialRod {
    rod: Rod,
    section: fs_material::section::ConcentricShearSection,
    reference_length_m: f64,
    reference_mass_kg: f64,
}

impl MaterialRod {
    /// Build the actual rod, clamped at node zero, without hand-copied stiffness.
    /// Work and allocation remain within the existing fixture-scale solver.
    ///
    /// # Errors
    /// Refuses nonpositive/nonfinite length, zero or oversized segment counts,
    /// unrepresentable segment length/reference mass and unsupported source laws.
    #[allow(clippy::cast_precision_loss)]
    pub fn try_straight(
        length_m: f64,
        segments: usize,
        section: fs_material::section::ConcentricShearSection,
    ) -> Result<Self, SolidError> {
        if !length_m.is_finite() || length_m <= 0.0
            || segments == 0 || segments > MAX_MATERIAL_ROD_SEGMENTS
        {
            return Err(SolidError::InvalidInput {
                what: "material rod needs finite positive length and 1..=64 segments".into(),
            });
        }
        let segment_length = length_m / segments as f64;
        let mass = length_m * section.elastic().line_density_kg_m();
        if !segment_length.is_finite() || segment_length <= 0.0 || !mass.is_finite() || mass <= 0.0 {
            return Err(SolidError::InvalidInput {
                what: "material rod segment length or reference mass is unrepresentable".into(),
            });
        }
        let coefficients = RodSection::try_from_concentric_material(&section)?;
        Ok(Self {
            rod: Rod::straight(length_m, segments, coefficients), section,
            reference_length_m: length_m, reference_mass_kg: mass,
        })
    }

    /// Solved geometry, nodal frames, energy and strains from the existing rod.
    /// Read-only access prevents silently replacing its bound section coefficients.
    #[must_use]
    pub const fn rod(&self) -> &Rod { &self.rod }
    /// Original ordered material bundles, reference geometry and shear correction.
    #[must_use]
    pub const fn material_section(&self) -> &fs_material::section::ConcentricShearSection { &self.section }
    /// Physical stress-free reference length [m].
    #[must_use]
    pub const fn reference_length_m(&self) -> f64 { self.reference_length_m }
    /// Reference mass [kg], unchanged by elastic extension in this static model.
    #[must_use]
    pub const fn reference_mass_kg(&self) -> f64 { self.reference_mass_kg }

    /// Solve with the existing energy/Newton/Lie-group implementation.
    /// Only a fully successful, finite result replaces the previous geometry.
    /// This preserves the legacy solver's numerical and no-dynamics boundaries;
    /// a small residual is not a model-accuracy or stability certificate.
    ///
    /// # Errors
    /// Refuses nonfinite loads, tolerance outside (0,infinity), load-step count
    /// outside 1..=256, or any solver failure. Material and prior pose are unchanged.
    pub fn solve_static(
        &mut self,
        load: &TipLoad,
        steps: usize,
        tolerance: f64,
    ) -> Result<Vec<f64>, SolidError> {
        if steps == 0 || steps > 256 || !tolerance.is_finite() || tolerance <= 0.0
            || !load.force.iter().chain(&load.moment).all(|value| value.is_finite())
        {
            return Err(SolidError::InvalidInput {
                what: "material rod solve needs finite loads, positive tolerance and 1..=256 load steps".into(),
            });
        }
        let mut candidate = self.rod.clone();
        let residuals = candidate.solve_static(load, steps, tolerance)?;
        if !candidate.positions.iter().flatten().all(|value| value.is_finite())
            || !candidate.energy()?.is_finite()
        {
            return Err(SolidError::InvalidInput {
                what: "material rod solve produced an unrepresentable configuration".into(),
            });
        }
        self.rod = candidate;
        Ok(residuals)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
    use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection};
    use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point};
    use fs_qty::{Density, Pressure};

    fn section(sleeve_scale: f64, band: Option<(f64, f64)>) -> ConcentricShearSection {
        let layers = [(0.05, 1.0), (0.1, sleeve_scale)].into_iter().map(|(radius, scale)| {
            let mut claims = ClaimSet::new();
            let mut requirements = Vec::new();
            for (name, value, dims) in [("density", 1000.0*scale, Density::DIMS),
                ("young_modulus", 1e6*scale, Pressure::DIMS), ("shear_modulus", 4e5*scale, Pressure::DIMS)] {
                claims.insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
                    validity: band.map_or_else(ValidityDomain::unconstrained,
                        |(lo, hi)| ValidityDomain::unconstrained().with("omega", lo, hi)),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
                    provenance: Provenance { source: "synthetic layered rod".into(), license: "test-only".into(), artifact: None },
                }).unwrap();
                requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive).unwrap());
            }
            let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(),
                phase: "solid".into(), process: "test".into(), revision: 0 }, claims, vec![]).unwrap();
            let point = band.map_or_else(QueryPoint::new,
                |(lo, hi)| QueryPoint::new().with("omega", (lo + hi) / 2.0).unwrap());
            ConcentricMaterialLayer { outer_radius_m: radius, material: resolve_material_state_point(
                &card, &point, &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap() }
        }).collect();
        ConcentricElasticSection::try_new(0.0, layers).unwrap().with_shear_modulus("shear_modulus", 0.8).unwrap()
    }
    fn near(a: f64, b: f64) { assert!((a-b).abs() < 2e-11 * a.abs().max(b.abs()), "{a} != {b}"); }

    #[test]
    fn shared_layer_integrals_drive_actual_axial_shear_torsion_and_bending_energy() {
        let section = section(2.0, None);
        let coefficients = RodSection::try_from_concentric_material(&section).unwrap();
        assert_eq!(coefficients.ea, section.elastic().axial_stiffness_n());
        assert_eq!(coefficients.ei, section.elastic().bending_stiffness_n_m2());
        let mut rod = Rod::straight(1.0, 1, coefficients);
        rod.positions[1] = [1.001, 0.002, 0.003];
        near(rod.energy().unwrap(), 0.5 * (coefficients.ea*0.001_f64.powi(2)
            + coefficients.ga*(0.002_f64.powi(2) + 0.003_f64.powi(2))));
        rod.positions[1] = [1.0, 0.0, 0.0];
        rod.frames[1] = So3::identity().body_plus(fs_ga::So3Tangent::new(Vec3::new(0.01, 0.0, 0.0))).unwrap();
        near(rod.energy().unwrap(), 0.5 * coefficients.gj * 0.01_f64.powi(2));
        rod.frames[1] = So3::identity().body_plus(fs_ga::So3Tangent::new(Vec3::new(0.0, 0.0, 0.01))).unwrap();
        let midpoint = So3::identity().body_plus(fs_ga::So3Tangent::new(Vec3::new(0.0, 0.0, 0.005))).unwrap();
        let tip = midpoint.rotate(Vec3::new(1.0, 0.0, 0.0)).unwrap();
        rod.positions[1] = [tip.x, tip.y, tip.z];
        near(rod.energy().unwrap(), 0.5 * coefficients.ei * 0.01_f64.powi(2));
    }

    #[test]
    fn replacing_sleeve_data_changes_solved_extension_and_twist_not_the_solver() {
        let mut a = MaterialRod::try_straight(1.0, 1, section(1.0, None)).unwrap();
        let mut b = MaterialRod::try_straight(1.0, 1, section(2.0, None)).unwrap();
        let original_section = b.material_section().clone();
        let mass = b.reference_mass_kg();
        let load = TipLoad { force: [1.0, 0.0, 0.0], moment: [0.1, 0.0, 0.0] };
        for rod in [&mut a, &mut b] {
            rod.solve_static(&load, 1, 1e-8).unwrap();
            let (gamma, kappa) = rod.rod().strains(0).unwrap();
            let expected_extension = load.force[0] / rod.rod().section.ea;
            let expected_twist = load.moment[0] / rod.rod().section.gj;
            assert!((gamma[0] / expected_extension - 1.0).abs() < 1e-5);
            assert!((kappa[0] / expected_twist - 1.0).abs() < 1e-5);
        }
        assert!(b.rod().positions[1][0] < a.rod().positions[1][0]);
        assert!(b.rod().strains(0).unwrap().1[0] < a.rod().strains(0).unwrap().1[0]);
        assert_eq!(b.material_section(), &original_section);
        assert_eq!(b.reference_mass_kg(), mass);
        assert!(mass > a.reference_mass_kg());
    }

    #[test]
    fn static_material_binding_does_not_extrapolate_dynamic_source_support() {
        assert!(MaterialRod::try_straight(1.0, 1, section(1.0, Some((1.0, 100.0)))).is_err());
        assert!(MaterialRod::try_straight(1.0, 1, section(1.0, Some((0.0, 100.0)))).is_ok());
    }

    #[test]
    fn invalid_configuration_or_load_preserves_material_and_pose() {
        let source = section(1.0, None);
        for length in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(MaterialRod::try_straight(length, 1, source.clone()).is_err());
        }
        for segments in [0, MAX_MATERIAL_ROD_SEGMENTS + 1, usize::MAX] {
            assert!(MaterialRod::try_straight(1.0, segments, source.clone()).is_err());
        }
        let mut rod = MaterialRod::try_straight(1.0, 2, source.clone()).unwrap();
        let positions = rod.rod().positions.clone();
        assert!(rod.solve_static(&TipLoad::default(), 0, 1e-8).is_err());
        assert!(rod.solve_static(&TipLoad::default(), 1, f64::NAN).is_err());
        assert!(rod.solve_static(&TipLoad { force: [f64::NAN, 0.0, 0.0], moment: [0.0; 3] }, 1, 1e-8).is_err());
        assert_eq!(rod.rod().positions, positions);
        assert_eq!(rod.material_section(), &source);
        near(rod.reference_mass_kg(), source.elastic().line_density_kg_m());
    }

    #[test]
    fn linear_tip_deflection_matches_beam_theory() {
        let section = RodSection {
            ea: 1e4,
            ga: 1e4,
            gj: 1.0,
            ei: 1.0,
        };
        let mut rod = Rod::straight(1.0, 8, section);
        let p = 0.01; // PL²/EI = 0.01: linear regime
        let load = TipLoad {
            force: [0.0, p, 0.0],
            moment: [0.0; 3],
        };
        let hist = rod.solve_static(&load, 1, 1e-9).expect("linear solve");
        let want = p / (3.0 * section.ei); // PL³/3EI
        let got = rod.positions[8][1];
        assert!(
            (got - want).abs() < 0.02 * want,
            "tip {got} vs beam theory {want}; residual history {hist:?}"
        );
    }

    #[test]
    fn small_dead_tip_load_is_not_lost_against_large_stored_energy() {
        let mut rod = Rod::straight(
            1.0,
            2,
            RodSection {
                ea: 1e18,
                ga: 1e18,
                gj: 1e15,
                ei: 1e15,
            },
        );
        for position in &mut rod.positions {
            position[0] *= 1.01;
        }
        assert!(rod.energy().unwrap() > 1e12);

        let load = TipLoad {
            force: [0.0, 0.0, 1e-3],
            moment: [0.0; 3],
        };
        let factor = 0.25;
        let unloaded = rod
            .residual(
                &TipLoad::default(),
                factor,
                &Rod::energy,
                &|| Ok(()),
                1.0,
            )
            .unwrap();
        let loaded = rod
            .residual(&load, factor, &Rod::energy, &|| Ok(()), 1.0)
            .unwrap();
        let tip_z = rod.ndof() - 4;

        // Reflection symmetry makes this internal force exactly zero. The old
        // total-potential stencil rounds away F * eps beside the axial energy.
        assert_eq!(unloaded[tip_z], 0.0);
        assert_eq!(loaded[tip_z], -factor * load.force[2]);
        for k in 0..rod.ndof() {
            if k != tip_z {
                assert_eq!(loaded[k], unloaded[k]);
            }
        }
    }

}
