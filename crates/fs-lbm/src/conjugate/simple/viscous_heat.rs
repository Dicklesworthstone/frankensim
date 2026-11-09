//! Viscous power from the SAME source-resolved staggered stress operator.
//!
//! This entry solves steady Stokes flow in a closed Cartesian box, including
//! fixed no-slip solid voxels. Molecular viscosity may vary by fluid cell;
//! exterior walls may move tangentially. A declared stationary linear Darcy
//! frame and explicitly calibrated separable Forchheimer resistance may absorb
//! work as local heat. The separate pressure-work opt-in admits one opposed
//! reservoir normal-traction pair with zero tangential velocity trace. Inertia,
//! momentum convection, moving solids, turbulence and fans remain excluded.
//! Exterior walls, reservoirs and explicit body-force actuators supply work; neither
//! pressure-constraint work nor momentum residual is heat.
//! It has no kinetic-energy store and makes no Navier--Stokes energy claim.
//!
//! Normal stress uses `2 mu (du_a/dx_a)^2 V` at each cell. Each reciprocal MAC
//! shear edge contributes `mu_edge (du_a/dx_d + du_d/dx_a)^2 V_edge` once. The
//! original harmonic edge viscosity supplies the common stress; its half-cell
//! resistance fractions partition the heat. Only four-fluid-cell edges carry
//! that reciprocal shear. Blocked transverse faces use the original stationary
//! half-cell no-slip walls, including staircase corners. The same parallel wall
//! conductances supply traction work and local fluid heat; fixed internal walls
//! supply zero work. Solid voxels receive no fluid dissipation. No residual,
//! upwind diffusion or pressure-constraint work is converted into heat.
//! The material-momentum bulk opt-in replaces the cell normal contribution with
//! `2 mu sum_a(dev(D)_aa^2) V + zeta div(u)^2 V`; reciprocal shear-edge and
//! tangential-wall terms retain their original coefficients and work account.

use super::*;

pub(in crate::conjugate) mod material_momentum;
mod material_sgs;

pub(in crate::conjugate) mod pressure;
pub use pressure::creeping_flow_with_pressure_work_heat;

#[cfg(feature = "enthalpy-coupling")]
mod unsteady;
#[cfg(feature = "enthalpy-coupling")]
pub use unsteady::{UnsteadyHeatingAudit, UnsteadyMechanicalStep};

struct MolecularPower {
    heat: Vec<f64>,
    quadrature_w: f64,
    wall_power_w: [f64; 6],
    wall_tangential_force_n: [[f64; 3]; 6],
    molecular_heat_w: f64,
    /// Empty for the original molecular sampler; modeled transfer otherwise.
    sgs_transfer: Vec<f64>,
    sgs_transfer_w: f64,
}

#[cfg(test)]
mod strain_tests {
    use super::*;
    use fs_material::fluid_transport::viscous::{DeformationRate, NewtonianViscosity};

    #[test]
    fn mac_bulk_power_matches_the_original_material_response() {
        let domain = VoxelDomain::new(5, 5, 5, 1.0).unwrap();
        let config = SimpleConfig::new([FvBoundary::Symmetry; 6]);
        let fluid = FluidProperties {
            density_kg_m3: 1.0,
            specific_heat_j_kg_k: 1.0,
            conductivity_w_m_k: 1.0,
            kinematic_viscosity_m2_s: 1.0,
        };
        let gate = CancelGate::new_clock_free();
        let mut sampler =
            material_solver(&domain, &fluid, &[1.0; 125], None, None, &config, &gate).unwrap();
        let center = 2 + 5 * (2 + 5 * 2);
        for gradient in [
            [[0.003, 0.0, 0.0], [0.0, 0.003, 0.0], [0.0, 0.0, 0.003]],
            [[0.0, -0.004, 0.0], [0.004, 0.0, 0.0], [0.0, 0.0, 0.0]],
            [[-0.002, 0.006, 0.0], [0.0, -0.002, 0.0], [0.0, 0.0, -0.002]],
            [[0.004, 0.0, 0.0], [0.0, -0.002, 0.0], [0.0, 0.0, -0.002]],
        ] {
            for a in 0..3 {
                for index in 0..sampler.vel[a].len() {
                    let at = sampler.comps[a].coords(index);
                    sampler.vel[a][index] = (0..3)
                        .map(|d| gradient[a][d] * (at[d] as f64 + if a == d { 0.0 } else { 0.5 }))
                        .sum();
                }
            }
            let rate = DeformationRate::from_velocity_gradient(gradient).unwrap();
            for bulk in [0.0, 0.4, 7.0] {
                sampler.molecular_bulk_viscosity = vec![bulk; 125];
                let expected = NewtonianViscosity::new(1.0, bulk)
                    .unwrap()
                    .evaluate(rate)
                    .unwrap()
                    .dissipation_w_m3();
                let actual = sampler.sample_molecular_power::<false>(&gate).unwrap().heat[center];
                assert!((actual - expected).abs() < 2e-14);
                let observed = sampler.molecular_shear_rates(&gate).unwrap()[center];
                assert!((observed - rate.shear_rate_per_s()).abs() < 2e-14);
            }
        }
    }

    /// G0/G3: source observation removes spherical strain before squaring,
    /// while the existing full-stress operator retains its volumetric power.
    #[test]
    fn native_mac_shear_query_excludes_volume_without_losing_stress_power() {
        let domain = VoxelDomain::new(5, 5, 5, 1.0).unwrap();
        let config = SimpleConfig::new([FvBoundary::Symmetry; 6]);
        let fluid = FluidProperties {
            density_kg_m3: 1.0,
            specific_heat_j_kg_k: 1.0,
            conductivity_w_m_k: 1.0,
            kinematic_viscosity_m2_s: 1.0,
        };
        let gate = CancelGate::new_clock_free();
        let mut sampler =
            material_solver(&domain, &fluid, &[1.0; 125], None, None, &config, &gate).unwrap();
        for (gradient, expected) in [
            (
                [[0.003, 0.0, 0.0], [0.0, 0.003, 0.0], [0.0, 0.0, 0.003]],
                0.0,
            ),
            (
                [[0.0, -0.004, 0.0], [0.004, 0.0, 0.0], [0.0, 0.0, 0.0]],
                0.0,
            ),
            (
                [[-0.002, 0.006, 0.0], [0.0, -0.002, 0.0], [0.0, 0.0, -0.002]],
                0.006,
            ),
            (
                [[0.004, 0.0, 0.0], [0.0, -0.002, 0.0], [0.0, 0.0, -0.002]],
                0.004 * 3.0_f64.sqrt(),
            ),
        ] {
            for a in 0..3 {
                for index in 0..sampler.vel[a].len() {
                    let at = sampler.comps[a].coords(index);
                    sampler.vel[a][index] = (0..3)
                        .map(|d| gradient[a][d] * (at[d] as f64 + if a == d { 0.0 } else { 0.5 }))
                        .sum();
                }
            }
            // Interior quadrature is independent of the manufactured field's
            // boundary trace; no pressure projection changes the observation.
            let center = 2 + 5 * (2 + 5 * 2);
            let rate = sampler.molecular_shear_rates(&gate).unwrap()[center];
            assert!((rate - expected).abs() < 2e-14);
            let full = sampler.sample_molecular_power::<false>(&gate).unwrap().heat[center];
            let divergence = gradient[0][0] + gradient[1][1] + gradient[2][2];
            assert!(
                (full - expected * expected - (2.0 / 3.0) * divergence * divergence).abs() < 2e-14
            );
        }
    }
}

/// The checked value itself: the shared `finite` gate returns only the verdict.
fn finite_value(field: &'static str, value: f64) -> Result<f64, ChtError> {
    finite(field, value)?;
    Ok(value)
}

/// Explicit creeping-flow and algebraic work gates, not accuracy certificates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscousHeatControls {
    /// Based on the full box's longest dimension, maximum wall/face speed, and
    /// minimum molecular kinematic viscosity. Must be in (0, 0.1].
    pub maximum_reynolds: f64,
    /// Absolute W allowance for each independent power/partition constraint.
    pub power_tolerance_w: f64,
}

/// Actual molecular stress and stationary porous-loss power at one epoch.
#[derive(Debug, Clone, PartialEq)]
pub struct ViscousPower {
    /// Nonnegative physical heat in native voxel order [W].
    pub cell_heat_w: Vec<f64>,
    /// Total molecular, Darcy and Forchheimer heat [W], deposited exactly once.
    pub heat_w: f64,
    pub molecular_heat_w: f64,
    pub darcy_heat_w: f64,
    /// Physical cubic loss from the native rho*beta*|u|u force, not time loss.
    pub forchheimer_heat_w: f64,
    /// Darcy contribution in native voxel order. None without a declared drag
    /// frame; zero entries at zero drag. Already included in cell_heat_w.
    pub cell_darcy_heat_w: Option<Vec<f64>>,
    /// Cubic-loss subaccount, already included in cell_heat_w.
    /// None when all native quadratic coefficients are zero.
    pub cell_forchheimer_heat_w: Option<Vec<f64>>,
    /// Signed actuator supply at each Face3::ALL boundary [W].
    pub wall_power_w: [f64; 6],
    /// Tangential force exerted by each actuator on the fluid [N]. Normal
    /// reactions of impermeable walls are not included and do no work.
    pub wall_tangential_force_n: [[f64; 3]; 6],
    pub supplied_wall_power_w: f64,
    /// Signed reservoir normal-traction supply at each Face3::ALL mouth [W].
    /// Uses -(p_res-p_ref)*Q_out; zero outside the explicit pressure-work mode.
    pub pressure_power_w: [f64; 6],
    /// Mechanical supply only; physical stress/drag heat is deposited once.
    pub supplied_pressure_power_w: f64,
    /// Signed total body-force work, partitioned with the original MAC half
    /// cells. This is supply, not an additional heat source. None without forcing.
    /// A material unsteady receipt separately reports its buoyancy subset.
    pub cell_body_power_w: Option<Vec<f64>>,
    /// Total signed body-force power, including any material buoyancy [W].
    pub supplied_body_power_w: f64,
    /// Discrete sum p div(u) V [W]; independently required to be negligible.
    pub pressure_constraint_power_w: f64,
    /// Wall, reservoir and body supply plus pressure constraint minus heat [W].
    /// Steady flow requires this to vanish. In a transient mechanical receipt
    /// it is inertial power, checked against the native kinetic/time account.
    pub residual_w: f64,
    pub reynolds: f64,
}

#[cfg(feature = "enthalpy-coupling")]
impl ViscousPower {
    /// The exactly constrained, zero-drive system has no mechanical power.
    pub(in crate::conjugate) fn stationary(domain: &VoxelDomain, darcy: bool, body: bool) -> Self {
        Self {
            cell_heat_w: vec![0.0; domain.cell_count()],
            heat_w: 0.0,
            molecular_heat_w: 0.0,
            darcy_heat_w: 0.0,
            forchheimer_heat_w: 0.0,
            cell_darcy_heat_w: darcy.then(|| vec![0.0; domain.cell_count()]),
            cell_forchheimer_heat_w: None,
            wall_power_w: [0.0; 6],
            wall_tangential_force_n: [[0.0; 3]; 6],
            supplied_wall_power_w: 0.0,
            pressure_power_w: [0.0; 6],
            supplied_pressure_power_w: 0.0,
            cell_body_power_w: body.then(|| vec![0.0; domain.cell_count()]),
            supplied_body_power_w: 0.0,
            pressure_constraint_power_w: 0.0,
            residual_w: 0.0,
            reynolds: 0.0,
        }
    }
}

fn invalid_power(reason: &str) -> ChtError {
    ChtError::InvalidInput {
        field: "simple.viscous_heat",
        reason: reason.into(),
    }
}

pub(crate) fn admit_heating(
    domain: &VoxelDomain,
    config: &SimpleConfig,
    controls: ViscousHeatControls,
) -> Result<(), ChtError> {
    if domain.fluid_count() == 0
        || config
            .faces
            .iter()
            .any(|f| !matches!(f, FvBoundary::Wall { .. } | FvBoundary::Symmetry))
        || config.turbulence != Turbulence::Laminar
        || config.fan.is_some()
        || !config.internal_fans.is_empty()
        || !config.resistances.is_empty()
    {
        return Err(invalid_power(
            "viscous heat requires a closed box containing fluid with fixed no-slip solids and walls/symmetry, without fans, configured resistances or turbulence",
        ));
    }
    if !controls.maximum_reynolds.is_finite()
        || controls.maximum_reynolds <= 0.0
        || controls.maximum_reynolds > 0.1
        || !controls.power_tolerance_w.is_finite()
        || controls.power_tolerance_w <= 0.0
    {
        return Err(invalid_power(
            "positive finite power gate and Reynolds ceiling at most 0.1 required",
        ));
    }
    admit_molecular_boundaries(domain, config)
}

// Geometry/actuator admission is shared with the molecular-viscosity split,
// whose convection and pressure are owned by other existing substeps.
fn admit_molecular_boundaries(domain: &VoxelDomain, config: &SimpleConfig) -> Result<(), ChtError> {
    for (side, face) in config.faces.iter().enumerate() {
        if let FvBoundary::Wall { velocity } = face {
            if velocity.iter().any(|v| !v.is_finite())
                || velocity[side / 2] != 0.0
                || (0..3).any(|a| velocity[a] != 0.0 && domain.dims()[a] < 2)
            {
                return Err(invalid_power(
                    "walls must move tangentially with at least two cells along every driven component",
                ));
            }
            // A solid obstacle can remove every free tangential MAC face of a
            // driven wall patch even when the full box has several columns.
            // Refuse an unrepresented actuator instead of returning zero work.
            for c in 0..domain.cell_count() {
                if !domain.is_fluid(c) || domain.neighbor(c, side).is_some() {
                    continue;
                }
                for a in 0..3 {
                    if velocity[a] != 0.0
                        && !(0..2).any(|s| {
                            domain
                                .neighbor(c, 2 * a + s)
                                .is_some_and(|n| domain.is_fluid(n))
                        })
                    {
                        return Err(invalid_power(
                            "every moving wall patch needs adjacent fluid cells along its driven component",
                        ));
                    }
                }
            }
        }
    }
    Ok(())
}

/// Original MAC/SIMPLE solve with momentum convection explicitly omitted.
/// The power sampler consumes the RETURNED projected face field, so the thermal
/// owner never uses velocities different from those whose heat was audited.
/// Density, material claims and caloric state remain the caller's original owners.
pub fn creeping_flow_with_viscous_heat(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    viscosity_pa_s: &[f64],
    config: &SimpleConfig,
    controls: ViscousHeatControls,
    gate: &CancelGate,
) -> Result<(FvFlow, ViscousPower), ChtError> {
    creeping_material_heat(
        domain,
        fluid,
        viscosity_pa_s,
        None,
        None,
        config,
        controls,
        gate,
    )
}

/// The same native Stokes/stress owner with declared stationary Darcy drag.
///
/// `drag_pa_s_m2` is the coefficient beta of `-beta*u` [N/m3], one finite
/// nonnegative value per voxel. Native half-cell volumes supply both momentum
/// resistance and local dissipation. A stationary frame performs no reaction
/// work: every absorbed watt becomes heat in the mapped bulk/mixture cell.
/// The caloric model must enclose that frame at local thermal equilibrium;
/// separate solid/fluid temperatures or a moving porous skeleton need another
/// accounting boundary. No nonlinear inertial/Forchheimer loss is admitted.
pub fn creeping_flow_with_darcy_heat(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    viscosity_pa_s: &[f64],
    drag_pa_s_m2: &[f64],
    config: &SimpleConfig,
    controls: ViscousHeatControls,
    gate: &CancelGate,
) -> Result<(FvFlow, ViscousPower), ChtError> {
    creeping_material_heat(
        domain,
        fluid,
        viscosity_pa_s,
        Some(drag_pa_s_m2),
        None,
        config,
        controls,
        gate,
    )
}

/// Externally driven Stokes flow and its actual molecular/Darcy heat.
///
/// Force density [N/m3] is held over this source epoch. Its work uses exactly
/// the half-cell volumes and projected face velocities in native momentum.
/// Negative body work is allowed: an actuator can absorb wall-supplied power.
/// Only stress/Darcy dissipation heats cells; force work is never deposited a
/// second time. Optional Darcy resistance retains the stationary-frame model.
/// This does not debit a finite electrical source or gravitational potential;
/// the caller must include the reported work in any larger accounting boundary.
#[allow(clippy::too_many_arguments)]
pub fn creeping_flow_with_body_force_heat(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    viscosity_pa_s: &[f64],
    force_density_n_m3: &[[f64; 3]],
    drag_pa_s_m2: Option<&[f64]>,
    config: &SimpleConfig,
    controls: ViscousHeatControls,
    gate: &CancelGate,
) -> Result<(FvFlow, ViscousPower), ChtError> {
    creeping_material_heat(
        domain,
        fluid,
        viscosity_pa_s,
        drag_pa_s_m2,
        Some(force_density_n_m3),
        config,
        controls,
        gate,
    )
}

#[allow(clippy::too_many_arguments)]
fn creeping_material_heat(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    viscosity_pa_s: &[f64],
    drag_pa_s_m2: Option<&[f64]>,
    force_density_n_m3: Option<&[[f64; 3]]>,
    config: &SimpleConfig,
    controls: ViscousHeatControls,
    gate: &CancelGate,
) -> Result<(FvFlow, ViscousPower), ChtError> {
    creeping_material_heat_with_darcy_tensor(
        domain,
        fluid,
        viscosity_pa_s,
        drag_pa_s_m2,
        None,
        None,
        force_density_n_m3,
        config,
        controls,
        false,
        gate,
    )
}

/// Symmetric material drag shares the native force and heat quadrature.
#[allow(clippy::too_many_arguments)]
pub(in crate::conjugate) fn creeping_material_heat_with_darcy_tensor(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    viscosity_pa_s: &[f64],
    drag_pa_s_m2: Option<&[f64]>,
    darcy_resistance_pa_s_m2: Option<&[[[f64; 3]; 3]]>,
    forchheimer_resistance: Option<&[Option<ForchheimerResistance>]>,
    force_density_n_m3: Option<&[[f64; 3]]>,
    config: &SimpleConfig,
    controls: ViscousHeatControls,
    pressure_work: bool,
    gate: &CancelGate,
) -> Result<(FvFlow, ViscousPower), ChtError> {
    poll(gate)?;
    let pressure_work_axis = if pressure_work {
        Some(pressure::admit_pressure_work_heating(
            domain, config, controls,
        )?)
    } else {
        admit_heating(domain, config, controls)?;
        None
    };
    let mut solver = material_solver_with_darcy_tensor(
        domain,
        fluid,
        viscosity_pa_s,
        force_density_n_m3,
        drag_pa_s_m2,
        darcy_resistance_pa_s_m2,
        forchheimer_resistance,
        config,
        gate,
    )?;
    solver.creeping = true;
    solver.pressure_work_axis = pressure_work_axis;
    let flow = solve_steady(solver, fluid, config, gate)?;
    // Reconstruct only the original coefficient/geometry owner from its public
    // projected field. This is an audit, not another flow solve or projection.
    let mut sampler = material_solver_with_darcy_tensor(
        domain,
        fluid,
        viscosity_pa_s,
        force_density_n_m3,
        drag_pa_s_m2,
        darcy_resistance_pa_s_m2,
        forchheimer_resistance,
        config,
        gate,
    )?;
    sampler.pressure_work_axis = pressure_work_axis;
    let arrays = flow.field.face_arrays();
    for (a, fluxes) in arrays.into_iter().enumerate() {
        for (i, flux) in fluxes.iter().enumerate() {
            if i % 256 == 0 {
                poll(gate)?;
            }
            sampler.vel[a][i] = finite_value(
                "viscous_heat.velocity",
                *flux / sampler.comp_face_area(a, i),
            )?;
        }
    }
    let power = sampler.viscous_power(&flow, controls, gate)?;
    Ok((flow, power))
}

/// Kinematic observation of an already solved physical field. The unit
/// molecular coefficient removes viscosity weighting while preserving every
/// original cell, shear edge, wall and opt-in reservoir mouth quadrature.
/// This neither solves momentum nor creates a thermal source.
#[cfg(feature = "enthalpy-coupling")]
pub(in crate::conjugate) fn physical_flow_shear_rates(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    field: &FlowField,
    config: &SimpleConfig,
    controls: ViscousHeatControls,
    pressure_work: bool,
    gate: &CancelGate,
) -> Result<Vec<f64>, ChtError> {
    let pressure_work_axis = if pressure_work {
        Some(pressure::admit_pressure_work_heating(
            domain, config, controls,
        )?)
    } else {
        admit_heating(domain, config, controls)?;
        None
    };
    let mut unit = *fluid;
    unit.kinematic_viscosity_m2_s = 1.0 / unit.density_kg_m3;
    let mut sampler = material_solver(
        domain,
        &unit,
        &vec![1.0; domain.cell_count()],
        None,
        None,
        config,
        gate,
    )?;
    sampler.pressure_work_axis = pressure_work_axis;
    for (a, fluxes) in field.face_arrays().into_iter().enumerate() {
        if fluxes.len() != sampler.vel[a].len() {
            return Err(invalid_power(
                "physical shear field differs from the original native grid",
            ));
        }
        for (i, flux) in fluxes.iter().enumerate() {
            if i % 256 == 0 {
                poll(gate)?;
            }
            sampler.vel[a][i] = finite_value(
                "rheology.physical_velocity",
                *flux / sampler.comp_face_area(a, i),
            )?;
        }
    }
    sampler.molecular_shear_rates(gate)
}

impl Solver<'_> {
    // One source-based Reynolds bound for old and new native face fields.
    fn source_reynolds(&self, velocities: &[Vec<f64>; 3]) -> Result<f64, ChtError> {
        let domain = self.domain;
        let mut speed = velocities
            .iter()
            .flatten()
            .map(|u| u.abs())
            .fold(0.0, f64::max)
            * 3.0_f64.sqrt();
        for face in self.faces {
            if let FvBoundary::Wall { velocity } = face {
                speed = speed.max(velocity.iter().map(|u| u * u).sum::<f64>().sqrt());
            }
        }
        let length = (0..3)
            .map(|a| domain.face_coord(a, self.n[a]) - domain.face_coord(a, 0))
            .fold(0.0, f64::max);
        let minimum_mu = self
            .molecular_viscosity
            .iter()
            .enumerate()
            .filter(|(c, _)| domain.is_fluid(*c))
            .map(|(_, mu)| *mu)
            .fold(f64::INFINITY, f64::min);
        finite_value(
            "viscous_heat.reynolds",
            self.rho * speed * length / minimum_mu,
        )
    }

    fn viscous_power(
        &self,
        flow: &FvFlow,
        controls: ViscousHeatControls,
        gate: &CancelGate,
    ) -> Result<ViscousPower, ChtError> {
        let power =
            self.sample_viscous_power::<true>(&flow.field, &flow.pressure_pa, controls, gate)?;
        if power.residual_w.abs() > controls.power_tolerance_w
            || (power.supplied_wall_power_w
                + power.supplied_pressure_power_w
                + power.supplied_body_power_w
                - power.heat_w)
                .abs()
                > controls.power_tolerance_w
        {
            return Err(invalid_power(
                "independent steady wall/reservoir/body work and physical heat do not close",
            ));
        }
        Ok(power)
    }

    // Work of a total force or an explicitly retained force component on the
    // SAME endpoint faces and half-cell volumes used by momentum(). Keeping
    // this operator shared makes signed buoyancy attribution independent of
    // the physical stress/drag heat partition and cell-centred diagnostics.
    fn body_force_power(
        &self,
        force_n_m3: &[[f64; 3]],
        controls: ViscousHeatControls,
        gate: &CancelGate,
    ) -> Result<(Vec<f64>, f64), ChtError> {
        let domain = self.domain;
        if force_n_m3.len() != domain.cell_count() {
            return Err(invalid_power("body-force work needs one vector per voxel"));
        }
        let add = |value: &mut f64, amount: f64| -> Result<(), ChtError> {
            *value = finite_value("viscous_heat.power", *value + amount)?;
            Ok(())
        };
        let mut body_work = vec![0.0; domain.cell_count()];
        for (c, work) in body_work.iter_mut().enumerate() {
            poll(gate)?;
            if force_n_m3[c]
                .iter()
                .any(|f| !f.is_finite() || (!domain.is_fluid(c) && *f != 0.0))
            {
                return Err(invalid_power(
                    "body-force work needs finite, zero-solid forces",
                ));
            }
            if !domain.is_fluid(c) {
                continue;
            }
            let at = domain.coords(c);
            for a in 0..3 {
                for plus in [false, true] {
                    add(
                        work,
                        0.5 * domain.volume(c)
                            * force_n_m3[c][a]
                            * self.vel[a][self.cell_face(a, at, plus)],
                    )?;
                }
            }
        }
        // Independently traverse the exact RHS assembled in momentum(), on
        // its unknown faces. Fixed normal walls have zero velocity.
        let mut staggered_work = 0.0;
        for a in 0..3 {
            for &index in &self.comps[a].unknowns {
                poll(gate)?;
                let (minus, plus, weight) = self.power_face_cells(a, index);
                let fm = force_n_m3[cell_of(self.n, minus)][a];
                let fp = force_n_m3[cell_of(self.n, plus)][a];
                let force = self.comp_face_area(a, index)
                    * (0.5 * domain.width(a, minus[a]))
                        .mul_add(fm, 0.5 * domain.width(a, plus[a]) * fp);
                add(&mut staggered_work, weight * force * self.vel[a][index])?;
            }
        }
        let total = finite_value("viscous_heat.body_supply", body_work.iter().sum())?;
        if (total - staggered_work).abs() > controls.power_tolerance_w {
            return Err(invalid_power(
                "body-force work and original-cell partition disagree",
            ));
        }
        Ok((body_work, staggered_work))
    }

    // The stress/drag quadratures are shared by steady and inertial Stokes.
    // Only the caller owns the storage/time balance; no steady power equality
    // is imposed on an accelerating field.
    fn sample_molecular_power<const DEVIATORIC: bool>(
        &self,
        gate: &CancelGate,
    ) -> Result<MolecularPower, ChtError> {
        self.sample_stress_power::<DEVIATORIC, false>(gate)
    }

    /// The SAME native stress quadrature, with physical molecular dissipation
    /// and modeled SGS transfer kept in distinct original-cell accounts.
    fn sample_material_sgs_power(&self, gate: &CancelGate) -> Result<MolecularPower, ChtError> {
        if !self.has_sgs()
            || self.molecular_viscosity.len() != self.domain.cell_count()
            || self.molecular_bulk_viscosity.len() != self.domain.cell_count()
            || self.eddy_viscosity.len() != self.domain.cell_count()
        {
            return Err(invalid_power(
                "material SGS power requires original molecular, bulk and SGS fields",
            ));
        }
        self.sample_stress_power::<false, true>(gate)
    }

    fn sample_stress_power<const DEVIATORIC: bool, const SGS: bool>(
        &self,
        gate: &CancelGate,
    ) -> Result<MolecularPower, ChtError> {
        let domain = self.domain;
        let mut heat = vec![0.0; domain.cell_count()];
        let mut sgs_transfer = if SGS {
            vec![0.0; domain.cell_count()]
        } else {
            Vec::new()
        };
        let mut quadrature_w = 0.0;
        let mut wall_power_w = [0.0; 6];
        let mut wall_tangential_force_n = [[0.0; 3]; 6];
        let add = |value: &mut f64, amount: f64| -> Result<(), ChtError> {
            *value = finite_value("viscous_heat.power", *value + amount)?;
            Ok(())
        };
        for (c, cell_heat) in heat.iter_mut().enumerate() {
            poll(gate)?;
            if !domain.is_fluid(c) {
                continue;
            }
            let at = domain.coords(c);
            let normal: [f64; 3] = core::array::from_fn(|a| {
                (self.vel[a][self.cell_face(a, at, true)]
                    - self.vel[a][self.cell_face(a, at, false)])
                    / domain.width(a, at[a])
            });
            // Observe shear without spherical strain. An explicit bulk law
            // likewise separates its nonnegative heat from deviatoric heat;
            // the legacy full-stress path retains its original arithmetic.
            if !DEVIATORIC && let Some(&bulk) = self.molecular_bulk_viscosity.get(c) {
                let divergence: f64 = normal.iter().sum();
                let power = (bulk * divergence) * divergence * domain.volume(c);
                add(cell_heat, power)?;
                add(&mut quadrature_w, power)?;
            }
            let normal = if DEVIATORIC || !self.molecular_bulk_viscosity.is_empty() {
                core::array::from_fn(|a| {
                    (normal[a] - normal[(a + 1) % 3]) / 3.0
                        + (normal[a] - normal[(a + 2) % 3]) / 3.0
                })
            } else {
                normal
            };
            for gradient in normal {
                let power = 2.0 * self.mu_eff(at) * domain.volume(c) * gradient * gradient;
                if SGS {
                    let [molecular, modeled] =
                        material_sgs::split_power(self, c, self.eddy_viscosity[c], power)?;
                    add(cell_heat, molecular)?;
                    add(&mut sgs_transfer[c], modeled)?;
                } else {
                    add(cell_heat, power)?;
                }
                add(&mut quadrature_w, power)?;
            }
        }
        // One edge for each unordered pair of velocity components, including
        // graded cells. Normal fixed boundary DOFs carry no normal wall work.
        for a in 0..3 {
            for d in a + 1..3 {
                let third = 3 - a - d;
                for c in 0..domain.cell_count() {
                    poll(gate)?;
                    let plus = domain.coords(c);
                    if plus[a] == 0 || plus[d] == 0 {
                        continue;
                    }
                    let mut minus_a = plus;
                    minus_a[a] -= 1;
                    let mut minus_d = plus;
                    minus_d[d] -= 1;
                    let mut minus_both = minus_a;
                    minus_both[d] -= 1;
                    let cells = [minus_both, minus_d, minus_a, plus];
                    // Both reciprocal component pairs are active precisely
                    // when all four cells are fluid. A missing cell instead
                    // blocks the neighbouring MAC face: momentum uses its
                    // no-slip wall coefficient, accounted below.
                    if cells.iter().any(|&p| !domain.is_fluid(cell_of(self.n, p))) {
                        continue;
                    }
                    let la = 0.5 * (domain.width(a, minus_a[a]) + domain.width(a, plus[a]));
                    let ld = 0.5 * (domain.width(d, minus_d[d]) + domain.width(d, plus[d]));
                    let gradient = (self.vel[a][self.comps[a].index(plus)]
                        - self.vel[a][self.comps[a].index(minus_d)])
                        / ld
                        + (self.vel[d][self.comps[d].index(plus)]
                            - self.vel[d][self.comps[d].index(minus_a)])
                            / la;
                    let mu = self.edge_mu(minus_a, plus, d, false);
                    let power =
                        mu * la * ld * domain.width(third, plus[third]) * gradient * gradient;
                    add(&mut quadrature_w, power)?;
                    if SGS {
                        let eddy =
                            self.sgs_edge_eddy(minus_a, plus, d, false).ok_or_else(|| {
                                invalid_power("material SGS shear edge lacks its native closure")
                            })?;
                        let partition = material_sgs::edge_power(self, cells, [a, d], eddy, power)?;
                        for (p, [molecular, modeled]) in cells.into_iter().zip(partition) {
                            let cell = cell_of(self.n, p);
                            add(&mut heat[cell], molecular)?;
                            add(&mut sgs_transfer[cell], modeled)?;
                        }
                        continue;
                    }
                    let minimum = cells
                        .iter()
                        .map(|&p| self.mu_eff(p))
                        .fold(f64::INFINITY, f64::min);
                    let area: f64 = cells
                        .iter()
                        .map(|p| domain.width(a, p[a]) * domain.width(d, p[d]))
                        .sum();
                    let weights = cells.map(|p| {
                        domain.width(a, p[a]) * domain.width(d, p[d]) / area
                            * (minimum / self.mu_eff(p))
                    });
                    let total: f64 = weights.iter().sum();
                    for (p, weight) in cells.into_iter().zip(weights) {
                        add(&mut heat[cell_of(self.n, p)], power * (weight / total))?;
                    }
                }
            }
        }
        // Use exactly the parallel half-cell wall coefficient in momentum(),
        // on both exterior walls and stationary blocked transverse faces.
        // No fluid loss is deposited into the immobile obstacle.
        for a in 0..3 {
            for &index in &self.comps[a].unknowns {
                poll(gate)?;
                let f = self.comps[a].coords(index);
                let (minus, plus, weight) = self.power_face_cells(a, index);
                let wm = domain.width(a, minus[a]);
                let wp = domain.width(a, plus[a]);
                for d in 0..3 {
                    if d == a {
                        continue;
                    }
                    for side in [0, 1] {
                        let outside = f[d] == if side == 0 { 0 } else { self.n[d] - 1 };
                        let (boundary, wall_velocity) = if outside {
                            let face = 2 * d + side;
                            let FvBoundary::Wall { velocity } = self.faces[face] else {
                                continue;
                            };
                            (Some(face), velocity[a])
                        } else {
                            let mut next = f;
                            if side == 0 {
                                next[d] -= 1;
                            } else {
                                next[d] += 1;
                            }
                            if !matches!(
                                self.comps[a].kind[self.comps[a].index(next)],
                                Kind::Fixed(_)
                            ) {
                                continue;
                            }
                            // A blocked face has exactly zero velocity. The
                            // native staggered wall closure also applies at
                            // three-fluid-cell staircase corners.
                            (None, 0.0)
                        };
                        let wall_mu = self.wall_mu(minus, plus);
                        let coefficient =
                            weight * 2.0 * self.wall_shear_base_conductance(minus, plus, a, d);
                        let difference = wall_velocity - self.vel[a][index];
                        let power = coefficient * difference * difference;
                        add(&mut quadrature_w, power)?;
                        if let Some(face) = boundary {
                            add(
                                &mut wall_tangential_force_n[face][a],
                                coefficient * difference,
                            )?;
                            add(
                                &mut wall_power_w[face],
                                coefficient * difference * wall_velocity,
                            )?;
                        }
                        for (p, width) in [(minus, wm), (plus, wp)] {
                            // Resolved SGS has no wall function: the native
                            // wall coefficient and ALL its heat stay molecular.
                            let molecular = if SGS {
                                self.molecular_viscosity[cell_of(self.n, p)]
                            } else {
                                self.mu_eff(p)
                            };
                            add(
                                &mut heat[cell_of(self.n, p)],
                                power * (width / (wm + wp)) * (molecular / wall_mu),
                            )?;
                        }
                    }
                }
            }
        }
        self.pressure_mouth_power::<SGS>(&mut heat, &mut sgs_transfer, &mut quadrature_w, gate)?;
        let molecular_heat_w = finite_value("viscous_heat.molecular", heat.iter().sum())?;
        let sgs_transfer_w = finite_value("viscous_heat.sgs_transfer", sgs_transfer.iter().sum())?;
        Ok(MolecularPower {
            heat,
            quadrature_w,
            wall_power_w,
            wall_tangential_force_n,
            molecular_heat_w,
            sgs_transfer,
            sgs_transfer_w,
        })
    }

    /// Root-mean-square sqrt(2 dev(D):dev(D)) on the original symmetric-stress
    /// quadrature, including native half-cell walls and graded-cell edges.
    /// The observation must use a unit molecular coefficient, so its strain
    /// measure contains no viscosity, Darcy loss, pressure, forcing or LES.
    pub(in crate::conjugate) fn molecular_shear_rates(
        &self,
        gate: &CancelGate,
    ) -> Result<Vec<f64>, ChtError> {
        if !self.eddy_viscosity.is_empty()
            || self.molecular_viscosity.len() != self.domain.cell_count()
            || self.molecular_viscosity.iter().any(|mu| *mu != 1.0)
        {
            return Err(invalid_power(
                "kinematic strain observation requires the unit molecular operator",
            ));
        }
        let power = self.sample_molecular_power::<true>(gate)?;
        power
            .heat
            .into_iter()
            .enumerate()
            .map(|(cell, integral)| {
                poll(gate)?;
                let square = finite_value(
                    "rheology.strain_square",
                    integral / self.domain.volume(cell),
                )?;
                if square < 0.0 {
                    return Err(invalid_power("negative molecular strain invariant"));
                }
                finite_value("rheology.shear_rate", fs_math::det::sqrt(square))
            })
            .collect()
    }

    fn sample_viscous_power<const PHYSICAL_PRESSURE: bool>(
        &self,
        field: &FlowField,
        pressure_pa: &[f64],
        controls: ViscousHeatControls,
        gate: &CancelGate,
    ) -> Result<ViscousPower, ChtError> {
        let domain = self.domain;
        let MolecularPower {
            mut heat,
            mut quadrature_w,
            wall_power_w,
            wall_tangential_force_n,
            molecular_heat_w,
            ..
        } = self.sample_molecular_power::<false>(gate)?;
        let add = |value: &mut f64, amount: f64| -> Result<(), ChtError> {
            *value = finite_value("viscous_heat.power", *value + amount)?;
            Ok(())
        };
        let mut cell_forchheimer_heat_w = self
            .porous_nonlinear
            .then(|| vec![0.0; domain.cell_count()]);
        let mut staggered_forchheimer_w = 0.0;
        let (cell_darcy_heat_w, darcy_heat_w) = if self.porous.is_empty() {
            (None, 0.0)
        } else {
            let mut drag_heat = vec![0.0; domain.cell_count()];
            for (c, cell_heat) in drag_heat.iter_mut().enumerate() {
                poll(gate)?;
                if !domain.is_fluid(c) {
                    continue;
                }
                let at = domain.coords(c);
                for a in 0..3 {
                    let (beta, inertial) = self.porous[c][a];
                    // The two MAC half cells, not the square of an averaged
                    // cell-centre velocity. Fixed normal boundary faces are zero.
                    for plus in [false, true] {
                        let u = self.vel[a][self.cell_face(a, at, plus)];
                        add(cell_heat, (0.5 * domain.volume(c)) * beta * u * u)?;
                        if let Some(cubic) = &mut cell_forchheimer_heat_w {
                            add(
                                &mut cubic[c],
                                (0.5 * domain.volume(c)) * inertial * u.abs() * u * u,
                            )?;
                        }
                    }
                }
                if !self.porous_cross.is_empty() {
                    add(cell_heat, self.darcy_cross_power(c))?;
                    if *cell_heat < 0.0 {
                        return Err(invalid_power(
                            "the full symmetric Darcy loss is not representably nonnegative",
                        ));
                    }
                }
                if *cell_heat != 0.0 {
                    add(&mut heat[c], *cell_heat)?;
                }
                if let Some(cubic) = &mut cell_forchheimer_heat_w {
                    if !self.porous_framed.is_empty() {
                        add(&mut cubic[c], self.framed_forchheimer_power(c)?)?;
                    }
                    add(&mut heat[c], cubic[c])?;
                }
            }
            // Independent traversal of the exact implicit momentum coefficient
            // confirms its staggered work equals the original-cell heat split.
            let mut staggered_drag_w = 0.0;
            for a in 0..3 {
                for &index in &self.comps[a].unknowns {
                    poll(gate)?;
                    let (minus, plus, weight) = self.power_face_cells(a, index);
                    let wm = domain.width(a, minus[a]);
                    let wp = domain.width(a, plus[a]);
                    let bm = self.porous[cell_of(self.n, minus)][a].0;
                    let bp = self.porous[cell_of(self.n, plus)][a].0;
                    let coefficient =
                        self.comp_face_area(a, index) * (0.5 * wm).mul_add(bm, 0.5 * wp * bp);
                    let u = self.vel[a][index];
                    add(&mut staggered_drag_w, weight * coefficient * u * u)?;
                    if self.porous_nonlinear {
                        let cm = self.porous[cell_of(self.n, minus)][a].1;
                        let cp = self.porous[cell_of(self.n, plus)][a].1;
                        // This is the coefficient used by momentum at THIS
                        // frozen velocity, after the nonlinear residual gate.
                        let coefficient = self.comp_face_area(a, index)
                            * (0.5 * wm).mul_add(cm * u.abs(), 0.5 * wp * cp * u.abs());
                        add(&mut staggered_forchheimer_w, weight * coefficient * u * u)?;
                        if !self.porous_framed.is_empty() {
                            add(
                                &mut staggered_forchheimer_w,
                                weight * self.framed_forchheimer_force(a, minus, plus)? * u,
                            )?;
                        }
                    }
                    if !self.porous_cross.is_empty() {
                        add(
                            &mut staggered_drag_w,
                            weight * self.darcy_cross_force(a, minus, plus) * u,
                        )?;
                    }
                }
            }
            let total = finite_value("viscous_heat.darcy", drag_heat.iter().sum())?;
            if (total - staggered_drag_w).abs() > controls.power_tolerance_w {
                return Err(invalid_power(
                    "Darcy momentum work and original-cell heat partition disagree",
                ));
            }
            add(&mut quadrature_w, staggered_drag_w)?;
            (Some(drag_heat), total)
        };
        let forchheimer_heat_w = finite_value(
            "viscous_heat.forchheimer",
            cell_forchheimer_heat_w
                .as_ref()
                .map_or(0.0, |cells| cells.iter().sum()),
        )?;
        if (forchheimer_heat_w - staggered_forchheimer_w).abs() > controls.power_tolerance_w {
            return Err(invalid_power(
                "Forchheimer momentum work and original-cell heat partition disagree",
            ));
        }
        add(&mut quadrature_w, staggered_forchheimer_w)?;
        let heat_w = finite_value("viscous_heat.total", heat.iter().sum())?;
        let supplied_wall_power_w = finite_value("viscous_heat.supply", wall_power_w.iter().sum())?;
        let (pressure_power_w, supplied_pressure_power_w) = self.pressure_work(field, gate)?;
        let (cell_body_power_w, supplied_body_power_w) = if self.force.is_empty() {
            (None, 0.0)
        } else {
            let (body_work, staggered_work) = self.body_force_power(&self.force, controls, gate)?;
            (Some(body_work), staggered_work)
        };
        let mut pressure_constraint_power_w = 0.0;
        for c in 0..domain.cell_count() {
            poll(gate)?;
            if !domain.is_fluid(c) {
                continue;
            }
            let outward: f64 = (0..6).map(|face| field.outward(domain, c, face)).sum();
            add(
                &mut pressure_constraint_power_w,
                (if PHYSICAL_PRESSURE && self.pressure_work_axis.is_some() {
                    pressure_pa[c] - self.physical_pressure_reference_pa(c)
                } else {
                    pressure_pa[c]
                }) * outward,
            )?;
        }
        let residual_w = finite_value(
            "viscous_heat.residual",
            supplied_wall_power_w
                + supplied_pressure_power_w
                + supplied_body_power_w
                + pressure_constraint_power_w
                - heat_w,
        )?;
        let reynolds = self.source_reynolds(&self.vel)?;
        if reynolds > controls.maximum_reynolds {
            return Err(invalid_power(
                "source viscosity and actual wall/flow speed exceed the declared creeping-flow bound",
            ));
        }
        if (heat_w - quadrature_w).abs() > controls.power_tolerance_w
            || (heat_w - molecular_heat_w - darcy_heat_w - forchheimer_heat_w).abs()
                > controls.power_tolerance_w
            || pressure_constraint_power_w.abs() > controls.power_tolerance_w
        {
            return Err(invalid_power(
                "independent constraint, molecular, Darcy and Forchheimer power partitions do not close",
            ));
        }
        poll(gate)?;
        Ok(ViscousPower {
            cell_heat_w: heat,
            heat_w,
            molecular_heat_w,
            darcy_heat_w,
            cell_darcy_heat_w,
            forchheimer_heat_w,
            cell_forchheimer_heat_w,
            wall_power_w,
            wall_tangential_force_n,
            supplied_wall_power_w,
            pressure_power_w,
            supplied_pressure_power_w,
            cell_body_power_w,
            supplied_body_power_w,
            pressure_constraint_power_w,
            residual_w,
            reynolds,
        })
    }
}

#[cfg(test)]
mod opening_pressure_tests {
    use super::*;

    /// G1/G3: the actual native traction-mouth solve and its original work
    /// quadrature see the same independently driven reservoir pressures.
    #[test]
    fn independent_opening_loads_preserve_native_stokes_work_and_heat() {
        let domain =
            VoxelDomain::graded([vec![0.5, 1.0, 1.5], vec![0.75, 1.25], vec![1.0]]).unwrap();
        let fluid = FluidProperties {
            density_kg_m3: 1.0,
            specific_heat_j_kg_k: 2.0,
            conductivity_w_m_k: 0.02,
            kinematic_viscosity_m2_s: 100.0,
        };
        let mut faces = [FvBoundary::Symmetry; 6];
        faces[0] = FvBoundary::Outlet;
        faces[1] = FvBoundary::Outlet;
        faces[2] = FvBoundary::wall();
        faces[3] = FvBoundary::wall();
        let mut config = SimpleConfig::new(faces);
        config.tolerance = 1e-12;
        config.momentum_tolerance = 1e-10;
        config.pressure_tolerance = 1e-12;
        config.pressure_profiles_pa[0] = Some(vec![2.0, 6.0]);
        config.pressure_profiles_pa[1] = Some(vec![0.0; 2]);
        let viscosity = vec![100.0; domain.cell_count()];
        let controls = ViscousHeatControls {
            maximum_reynolds: 0.1,
            power_tolerance_w: 1e-7,
        };
        let gate = CancelGate::new_clock_free();
        let solve = |common: f64| {
            let offsets = [common + 1.0, common - 2.0, 0.0, 0.0, 0.0, 0.0];
            let mut solver =
                material_solver(&domain, &fluid, &viscosity, None, None, &config, &gate).unwrap();
            solver.creeping = true;
            solver.pressure_work_axis = Some(0);
            solver.set_pressure_loads(0.5, offsets, &gate).unwrap();
            let flow = solve_steady(solver, &fluid, &config, &gate).unwrap();
            // Reuse the original audit owner on the returned physical field.
            let mut sampler =
                material_solver(&domain, &fluid, &viscosity, None, None, &config, &gate).unwrap();
            sampler.pressure_work_axis = Some(0);
            sampler.set_pressure_loads(0.5, offsets, &gate).unwrap();
            for (axis, fluxes) in flow.field.face_arrays().into_iter().enumerate() {
                for (face, flux) in fluxes.iter().enumerate() {
                    sampler.vel[axis][face] = flux / sampler.comp_face_area(axis, face);
                }
            }
            let power = sampler.viscous_power(&flow, controls, &gate).unwrap();
            (flow, power)
        };
        let (flow, heat) = solve(0.0);
        // Actual mouths are [3,5] Pa on XMin and -1 Pa on XMax.
        // Their current reference is 3 Pa, so inlet heads are [0,2],
        // and both outlet heads are -4 Pa.
        let supply = 2.0 * flow.field.flux_x(0, 1, 0)
            + 4.0 * (flow.field.flux_x(3, 0, 0) + flow.field.flux_x(3, 1, 0));
        assert!(supply > 0.0);
        assert!((heat.supplied_pressure_power_w - supply).abs() < 1e-10);
        assert!((heat.heat_w - supply).abs() < controls.power_tolerance_w);
        let (shifted, shifted_heat) = solve(1_099_511_627_776.0);
        assert_eq!(flow.velocity_m_s, shifted.velocity_m_s);
        assert_eq!(flow.field, shifted.field);
        assert_eq!(heat.cell_heat_w, shifted_heat.cell_heat_w);
        assert_eq!(heat.pressure_power_w, shifted_heat.pressure_power_w);
        assert!(shifted_heat.pressure_constraint_power_w.abs() < controls.power_tolerance_w);
        assert!(shifted_heat.residual_w.abs() < controls.power_tolerance_w);
    }
}
