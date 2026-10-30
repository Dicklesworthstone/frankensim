//! Finite bonded patches on the original native mechanical network.
//!
//! Each quadrature site retains its own material, reference area and irreversible
//! damage. Complete work columns include the original translational, rotational
//! and elastic modes. All sites solve together with every cross-compliance term,
//! so peeling at one edge redistributes load to the rest of the bond. Columns
//! describe fixed local frames and linearized geometry; this is not a finite-
//! rotation contact search or an automatic choice of quadrature resolution.

use super::{
    BilinearCohesive, CancelGate, CohesiveError, CohesivePoint, CoupledModalSystem,
    MaterialStateId, MaterialStatePointError, ModalAcousticTimeModel, ModalCouplingError,
    ResolvedBilinearCohesive, dot, finite, invalid, limit, network_response, poll,
};
use crate::cohesive_contact::patch::{
    CohesivePatchConfig, CohesivePatchError, CohesivePatchLoadLine, MAX_COHESIVE_PATCH_POINTS,
};

/// Joint progressive-damage friction with actual-body work and atomic history.
pub mod frictional_damage;

/// One original bonded area and its work-conjugate local motion on all modes.
/// Geometry, quadrature weight and bonding are explicit caller declarations.
#[derive(Debug)]
pub struct CohesiveModalSite<'a> {
    law: &'a BilinearCohesive,
    source: Option<&'a ResolvedBilinearCohesive>,
    area_m2: f64,
    columns: [Vec<f64>; 3],
}
impl<'a> CohesiveModalSite<'a> {
    /// Declare a numerical material. The area is a reference quadrature weight,
    /// not an additional total patch area or a damage-scaled contact fraction.
    /// Column order is positive opening, first shear, second shear.
    pub fn new(
        law: &'a BilinearCohesive,
        area_m2: f64,
        columns: [Vec<f64>; 3],
    ) -> Result<Self, CohesivePatchModalError> {
        if !area_m2.is_finite()
            || area_m2 <= 0.0
            || columns.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "cohesive site needs positive finite area and finite work columns",
            )
            .into());
        }
        Ok(Self {
            law,
            source: None,
            area_m2,
            columns,
        })
    }

    /// Retain the actual ordered interface and fixed source query at this site.
    /// The two supplied identities must be those assigned to its receiving faces.
    pub fn from_source(
        source: &'a ResolvedBilinearCohesive,
        surface_materials: [&MaterialStateId; 2],
        area_m2: f64,
        columns: [Vec<f64>; 3],
    ) -> Result<Self, CohesivePatchModalError> {
        source.require_surface_materials(surface_materials[0], surface_materials[1])?;
        let mut site = Self::new(source.law(), area_m2, columns)?;
        site.source = Some(source);
        Ok(site)
    }
    #[must_use]
    pub const fn law(&self) -> &'a BilinearCohesive {
        self.law
    }
    #[must_use]
    pub const fn source(&self) -> Option<&'a ResolvedBilinearCohesive> {
        self.source
    }
    #[must_use]
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    #[must_use]
    pub const fn work_columns(&self) -> &[Vec<f64>; 3] {
        &self.columns
    }
}

/// Physical solve limits and an explicit dense setup-work bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CohesivePatchModalConfig {
    /// The total material-call ceiling includes one actual-motion audit per site.
    pub patch: CohesivePatchConfig,
    /// Bound on `(3 * sites)^2 * modes` projection terms before allocation.
    pub max_setup_terms: usize,
    /// Absolute held-force mismatch at each local component [N].
    pub force_absolute_tolerance_n: f64,
    /// Relative held-force mismatch, in (0,1), scaled by actual forces.
    pub force_relative_tolerance: f64,
}

/// Refusals preserve every accepted point, mechanical state and observation.
#[derive(Debug)]
pub enum CohesivePatchModalError {
    Network(ModalCouplingError),
    Contact(CohesivePatchError),
    Friction(crate::cohesive_contact::patch::frictional_damage::DamageFrictionPatchError),
    Source(MaterialStatePointError),
    Material {
        point: usize,
        source: CohesiveError,
    },
    JumpMismatch {
        point: usize,
        component: usize,
        residual_m: f64,
        tolerance_m: f64,
    },
    ForceMismatch {
        point: usize,
        component: usize,
        residual_n: f64,
        tolerance_n: f64,
    },
}
impl core::fmt::Display for CohesivePatchModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native cohesive patch refused: {self:?}")
    }
}
impl std::error::Error for CohesivePatchModalError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Network(e) => Some(e),
            Self::Contact(e) => Some(e),
            Self::Friction(e) => Some(e),
            Self::Source(e) => Some(e),
            Self::Material { source, .. } => Some(source),
            _ => None,
        }
    }
}
impl From<ModalCouplingError> for CohesivePatchModalError {
    fn from(value: ModalCouplingError) -> Self {
        Self::Network(value)
    }
}
impl From<CohesivePatchError> for CohesivePatchModalError {
    fn from(value: CohesivePatchError) -> Self {
        Self::Contact(value)
    }
}
impl From<MaterialStatePointError> for CohesivePatchModalError {
    fn from(value: MaterialStatePointError) -> Self {
        Self::Source(value)
    }
}
impl From<crate::cohesive_contact::patch::frictional_damage::DamageFrictionPatchError> for CohesivePatchModalError {
    fn from(value: crate::cohesive_contact::patch::frictional_damage::DamageFrictionPatchError) -> Self {
        Self::Friction(value)
    }
}

/// One accepted site's physical response; all energy values include its area.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CohesivePatchPointFrame {
    pub area_m2: f64,
    pub jump_m: [f64; 3],
    pub damage: f64,
    /// Actual held force used to drive the bodies [N].
    pub mean_force_n: [f64; 3],
    /// Timed endpoint traction, not substituted for the work-conjugate mean [Pa].
    pub endpoint_traction_pa: [f64; 3],
    pub cohesive_energy_j: f64,
    pub fracture_dissipation_j: f64,
    pub viscous_dissipation_j: f64,
    pub damage_friction: Option<super::frictional_damage::DamageFrictionFrame>,
    pub cohesive_work_j: f64,
    pub compatibility_residual_m: [f64; 3],
    pub compatibility_tolerance_m: f64,
    pub constitutive_residual_n: [f64; 3],
    pub force_tolerance_n: [f64; 3],
}

/// Whole-system work balance, published only with all sites and bodies.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CohesivePatchModalFrame {
    pub sample: u64,
    pub observer_pressure_pa: f64,
    pub network_energy_j: f64,
    pub cohesive_energy_j: f64,
    pub external_work_j: f64,
    pub prescribed_work_j: f64,
    pub network_dissipation_j: f64,
    /// Intrinsic fracture is not silently assigned to heat.
    pub fracture_dissipation_j: f64,
    pub viscous_dissipation_j: f64,
    /// Sum of the damaged-area sliding losses, without an inferred heat sink.
    pub friction_dissipation_j: f64,
    pub cohesive_work_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub failed_points: usize,
    pub iterations: usize,
    pub kinetic_switches: usize,
    /// Includes all final actual-motion constitutive checks.
    pub material_evaluations: usize,
}

/// Heterogeneous finite bonded interface with one coupled native transaction.
///
/// All contact scratch and observations are allocated at construction. The
/// original network retains its own propagation allocation and capacity rules.
/// Source coefficients and reference frames stay fixed; failure does not heal
/// on closure. Post-failure compression remains with the original material law.
pub struct CohesivePatchModalSystem<'a> {
    network: CoupledModalSystem,
    sites: Vec<CohesiveModalSite<'a>>,
    points: Vec<CohesivePoint<'a>>,
    staged_points: Vec<CohesivePoint<'a>>,
    load_line: CohesivePatchLoadLine,
    config: CohesivePatchModalConfig,
    forces: Vec<f64>,
    free_jump: Vec<[f64; 3]>,
    prescribed_delta: Vec<[f64; 3]>,
    observations: Vec<CohesivePatchPointFrame>,
    staged_observations: Vec<CohesivePatchPointFrame>,
    last: Option<CohesivePatchModalFrame>,
    damage_friction: Option<frictional_damage::DamageFrictionPatchBinding<'a>>,
}
impl<'a> CohesivePatchModalSystem<'a> {
    /// Bond the current sample-zero configuration at zero local jump everywhere.
    /// Each site's complete column supplies both motion and reciprocal force;
    /// moments enter through its actual rotational/elastic mode shapes.
    pub fn new(
        network: CoupledModalSystem,
        sites: Vec<CohesiveModalSite<'a>>,
        config: CohesivePatchModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, CohesivePatchModalError> {
        poll(Some(gate))?;
        let n = sites.len();
        if network.samples_rendered() != 0
            || n == 0
            || !(1..=MAX_COHESIVE_PATCH_POINTS).contains(&config.patch.max_points)
            || n > config.patch.max_points
            || sites
                .iter()
                .any(|site| site.columns.iter().any(|c| c.len() != network.mode_count()))
            || !config.force_absolute_tolerance_n.is_finite()
            || config.force_absolute_tolerance_n <= 0.0
            || !config.force_relative_tolerance.is_finite()
            || !(0.0 < config.force_relative_tolerance && config.force_relative_tolerance < 1.0)
            || config.patch.contact.max_material_evaluations < 2 * n
        {
            return Err(invalid("cohesive patch requires a sample-zero network, complete columns and positive bounded work/force admission").into());
        }
        let dimension = 3 * n;
        let square = dimension
            .checked_mul(dimension)
            .ok_or_else(|| invalid("cohesive patch compliance extent overflow"))?;
        let terms = square
            .checked_mul(network.mode_count())
            .ok_or_else(|| invalid("cohesive patch setup work overflow"))?;
        if terms > config.max_setup_terms {
            return Err(invalid("cohesive patch projection exceeds max_setup_terms").into());
        }
        let mut compliance = Vec::new();
        compliance
            .try_reserve_exact(square)
            .map_err(|_| invalid("cohesive patch compliance allocation refused"))?;
        compliance.resize(square, 0.0);
        for j in 0..dimension {
            let response = network_response(&network, &sites[j / 3].columns[j % 3], gate)?;
            for i in 0..dimension {
                if i % 3 == 0 {
                    poll(Some(gate))?;
                }
                compliance[i * dimension + j] = dot(&sites[i / 3].columns[i % 3], &response)?;
            }
        }
        let mut solve_config = config.patch;
        solve_config.contact.max_material_evaluations -= n;
        let areas = sites.iter().map(|s| s.area_m2).collect();
        let load_line = CohesivePatchLoadLine::new(areas, compliance, solve_config)?;
        let points: Vec<_> = sites.iter().map(|s| s.law.bonded_point()).collect();
        let owner = Self {
            forces: vec![0.0; network.mode_count()],
            staged_points: points.clone(),
            points,
            sites,
            network,
            load_line,
            config,
            free_jump: vec![[0.0; 3]; n],
            prescribed_delta: vec![[0.0; 3]; n],
            observations: vec![CohesivePatchPointFrame::default(); n],
            staged_observations: vec![CohesivePatchPointFrame::default(); n],
            last: None,
            damage_friction: None,
        };
        limit(
            "initial energy including cohesive patch",
            owner.total_energy_j()?,
            owner.network.config.maximum_total_energy_j,
        )?;
        poll(Some(gate))?;
        Ok(owner)
    }
    #[must_use]
    pub fn sites(&self) -> &[CohesiveModalSite<'a>] {
        &self.sites
    }
    #[must_use]
    pub fn points(&self) -> &[CohesivePoint<'a>] {
        &self.points
    }
    #[must_use]
    pub fn components(&self) -> &[ModalAcousticTimeModel] {
        self.network.components()
    }
    #[must_use]
    pub fn mode_count(&self) -> usize {
        self.network.mode_count()
    }
    #[must_use]
    pub fn sample_period_s(&self) -> f64 {
        self.network.sample_period_s()
    }
    #[must_use]
    pub fn samples_rendered(&self) -> u64 {
        self.network.samples_rendered()
    }
    #[must_use]
    pub fn compliance_m_per_n(&self) -> &[f64] {
        self.load_line.compliance_m_per_n()
    }
    #[must_use]
    pub const fn config(&self) -> CohesivePatchModalConfig {
        self.config
    }
    #[must_use]
    pub const fn last_frame(&self) -> Option<&CohesivePatchModalFrame> {
        self.last.as_ref()
    }
    #[must_use]
    pub fn last_point_frames(&self) -> Option<&[CohesivePatchPointFrame]> {
        self.last.as_ref().map(|_| self.observations.as_slice())
    }

    /// Original native storage plus the sum over reference bonded quadrature.
    pub fn total_energy_j(&self) -> Result<f64, CohesivePatchModalError> {
        let mut energy = self.network.total_energy_j()?;
        for (site, point) in self.sites.iter().zip(&self.points) {
            energy = finite(energy + site.area_m2 * point.stored_energy_j_m2())?;
        }
        Ok(energy)
    }

    /// Hold external modal forces and each site's declared jump velocity [m/s]
    /// for one physical sample. Supply zero velocities for internal bonds.
    pub fn step(
        &mut self,
        external: &[f64],
        prescribed_jump_rate_m_s: &[[f64; 3]],
    ) -> Result<&CohesivePatchModalFrame, CohesivePatchModalError> {
        self.step_inner(external, prescribed_jump_rate_m_s, None)
    }
    pub fn step_under_gate(
        &mut self,
        external: &[f64],
        prescribed_jump_rate_m_s: &[[f64; 3]],
        gate: &CancelGate,
    ) -> Result<&CohesivePatchModalFrame, CohesivePatchModalError> {
        self.step_inner(external, prescribed_jump_rate_m_s, Some(gate))
    }

    fn step_inner(
        &mut self,
        external: &[f64],
        prescribed_rate: &[[f64; 3]],
        gate: Option<&CancelGate>,
    ) -> Result<&CohesivePatchModalFrame, CohesivePatchModalError> {
        if self.damage_friction.is_some() {
            return self.step_damage_friction(external, prescribed_rate, gate);
        }
        poll(gate)?;
        if prescribed_rate.len() != self.points.len() {
            return Err(invalid("prescribed cohesive velocities must match all sites").into());
        }
        let before = self.total_energy_j()?;
        let dt = self.network.dt;
        for (delta, rate) in self.prescribed_delta.iter_mut().zip(prescribed_rate) {
            for i in 0..3 {
                delta[i] = finite(rate[i] * dt)?;
            }
        }
        self.network.prepare_prediction(external, gate)?;
        for p in 0..self.points.len() {
            poll(gate)?;
            let old = self.points[p].jump_m();
            for i in 0..3 {
                self.free_jump[p][i] = finite(
                    old[i]
                        + self.prescribed_delta[p][i]
                        + dot(&self.sites[p].columns[i], &self.network.free_delta)?,
                )?;
            }
        }
        let solved = self
            .load_line
            .trial(&self.points, &self.free_jump, dt, || {
                !gate.is_some_and(CancelGate::is_requested)
            })?;
        self.forces.copy_from_slice(external);
        for p in 0..self.points.len() {
            poll(gate)?;
            for i in 0..3 {
                for (force, &weight) in self.forces.iter_mut().zip(&self.sites[p].columns[i]) {
                    *force = finite(*force - weight * solved.mean_force_n[p][i])?;
                }
            }
        }
        self.network.stage_inner(&self.forces, gate)?;
        let mut cohesive_energy = 0.0;
        let mut fracture = 0.0;
        let mut viscous = 0.0;
        let mut cohesive_work = 0.0;
        let mut prescribed_work = 0.0;
        let mut failed_points = 0;
        for p in 0..self.points.len() {
            poll(gate)?;
            let old = self.points[p].jump_m();
            let mut actual_jump = [0.0; 3];
            let mut compatibility_residual = [0.0; 3];
            for i in 0..3 {
                actual_jump[i] = finite(
                    old[i]
                        + self.prescribed_delta[p][i]
                        + dot(&self.sites[p].columns[i], &self.network.free_delta)?,
                )?;
                compatibility_residual[i] = finite(actual_jump[i] - solved.jump_m[p][i])?;
                if compatibility_residual[i].abs() > solved.compatibility_tolerance_m[p] {
                    return Err(CohesivePatchModalError::JumpMismatch {
                        point: p,
                        component: i,
                        residual_m: compatibility_residual[i],
                        tolerance_m: solved.compatibility_tolerance_m[p],
                    });
                }
                limit(
                    "actual cohesive patch jump increment",
                    (actual_jump[i] - old[i]).abs(),
                    self.config.patch.contact.maximum_jump_increment_m,
                )?;
            }
            let (actual, mean) = self.points[p]
                .trial_over_time_with_mean_traction(
                    actual_jump,
                    dt,
                    self.config.patch.contact.material_energy_tolerance_j_m2,
                )
                .map_err(|source| CohesivePatchModalError::Material { point: p, source })?;
            let area = self.sites[p].area_m2;
            let mut force_residual = [0.0; 3];
            let mut force_tolerance = [0.0; 3];
            for i in 0..3 {
                let expected = finite(area * mean[i])?;
                limit(
                    "actual cohesive patch mean force",
                    expected.abs(),
                    self.config.patch.contact.maximum_force_n,
                )?;
                force_residual[i] = finite(solved.mean_force_n[p][i] - expected)?;
                force_tolerance[i] = finite(
                    self.config.force_absolute_tolerance_n
                        + self.config.force_relative_tolerance
                            * expected.abs().max(solved.mean_force_n[p][i].abs()),
                )?;
                if force_residual[i].abs() > force_tolerance[i] {
                    return Err(CohesivePatchModalError::ForceMismatch {
                        point: p,
                        component: i,
                        residual_n: force_residual[i],
                        tolerance_n: force_tolerance[i],
                    });
                }
            }
            let point_energy = finite(area * actual.equilibrium.point.stored_energy_j_m2())?;
            let point_fracture = finite(area * actual.equilibrium.energy.dissipated)?;
            let point_viscous = finite(area * actual.viscous_dissipated_j_m2)?;
            let point_work = finite(area * actual.total_work_j_m2())?;
            if point_fracture < 0.0 || point_viscous < 0.0 {
                return Err(invalid(
                    "cohesive patch fracture and viscosity must not create energy",
                )
                .into());
            }
            cohesive_energy = finite(cohesive_energy + point_energy)?;
            fracture = finite(fracture + point_fracture)?;
            viscous = finite(viscous + point_viscous)?;
            cohesive_work = finite(cohesive_work + point_work)?;
            prescribed_work =
                finite(prescribed_work + dot(&solved.mean_force_n[p], &self.prescribed_delta[p])?)?;
            failed_points += usize::from(actual.equilibrium.point.failed());
            self.staged_points[p] = actual.equilibrium.point;
            self.staged_observations[p] = CohesivePatchPointFrame {
                area_m2: area,
                jump_m: actual_jump,
                damage: actual.equilibrium.point.damage(),
                mean_force_n: solved.mean_force_n[p],
                endpoint_traction_pa: actual.traction_pa,
                cohesive_energy_j: point_energy,
                fracture_dissipation_j: point_fracture,
                viscous_dissipation_j: point_viscous,
                damage_friction: None,
                cohesive_work_j: point_work,
                compatibility_residual_m: compatibility_residual,
                compatibility_tolerance_m: solved.compatibility_tolerance_m[p],
                constitutive_residual_n: force_residual,
                force_tolerance_n: force_tolerance,
            };
        }
        let frame = &self.network.staged_frame;
        let network_energy = finite(frame.total_energy_j())?;
        let network_loss = finite(frame.total_dissipation_j())?;
        let external_work = dot(external, &self.network.free_delta)?;
        let after = finite(network_energy + cohesive_energy)?;
        limit(
            "total energy including finite cohesive patch",
            after,
            self.network.config.maximum_total_energy_j,
        )?;
        let loss = finite(network_loss + fracture + viscous)?;
        let residual = finite((after - before) + loss - external_work - prescribed_work)?;
        let scale = before
            .max(after)
            .max(loss.abs())
            .max(external_work.abs())
            .max(prescribed_work.abs());
        let tolerance = finite(
            self.network.config.energy_absolute_tolerance_j
                + self.network.config.energy_relative_tolerance * scale,
        )?;
        if residual.abs() > tolerance {
            return Err(ModalCouplingError::EnergyBalance {
                residual_j: residual,
                tolerance_j: tolerance,
            }
            .into());
        }
        let accepted = CohesivePatchModalFrame {
            sample: frame.sample,
            observer_pressure_pa: frame.observer_pressure_pa,
            network_energy_j: network_energy,
            cohesive_energy_j: cohesive_energy,
            external_work_j: external_work,
            prescribed_work_j: prescribed_work,
            network_dissipation_j: network_loss,
            fracture_dissipation_j: fracture,
            viscous_dissipation_j: viscous,
            friction_dissipation_j: 0.0,
            cohesive_work_j: cohesive_work,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            failed_points,
            iterations: solved.iterations,
            kinetic_switches: 0,
            material_evaluations: solved.material_evaluations + self.points.len(),
        };
        self.network.validate_staged_material_points(gate)?;
        poll(gate)?;
        self.network.publish_staged();
        std::mem::swap(&mut self.points, &mut self.staged_points);
        std::mem::swap(&mut self.observations, &mut self.staged_observations);
        Ok(self.last.insert(accepted))
    }
}
