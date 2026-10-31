//! Coupled discrete vapor pockets at native DG endpoint and interior ports.
//!
//! An interior pocket occupies an exact predeclared element face. Its shared
//! numerical pressure is the original DG pressure trace, and its volume rate
//! is Q = q_right - q_left, the original side-port flow jump. It approximates
//! a short pure-vapor separation between liquid columns, with no added liquid
//! inertia or finite-length moving mesh. The optional gas inventory adds its
//! original nonlinear volume response; surface tension and phase heat remain
//! outside this negligible-vapor-mass model.
//!
//! All active pressures are constrained TOGETHER by the Gram system
//! (B_A^T B_A + diag(K_g)) Q_A = -B_A^T L x, where K_g = -dp_g/dV.
//! Pure vapor has K_g=0 and retains its cached constant solve. One native
//! dense trajectory locates every rupture/collapse. Overlapping event
//! brackets are sampled at one common
//! time and resolved by the coupled contact complementarity, not port order.
//! No event resets pressure, liquid coordinates or integrated cavity volume.

use super::*;
use fs_duct::dg::continuous::state::DgTrajectoryPort;
use fs_la::factor::{Cholesky, cholesky};
use fs_time::adaptive::events::multiple::{
    EventOrder, EventSetAdvance, EventSetOptions, EventSpec, rk45_until_any_event,
};
use fs_time::adaptive::events::{EventDirection, InitialEvent};
use fs_time::adaptive::{AdaptiveState, AdaptiveStatus, PiController};
use std::cell::RefCell;

/// Fixed noncondensable gas on the SAME coupled pocket trajectory.
pub mod gas;
use gas::{GasFlowWorkspace, TrappedGasPocket};

/// Bounds the coupled Gram solve and the 2^k contact-branch search.
pub const MAX_DUCT_POCKETS: usize = 8;

/// One cavity on an ALREADY declared native endpoint or interior side port.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctPocket {
    pub port: PortId,
    pub initial_volume_m3: f64,
    /// Local ceiling, also limited to 10% of smaller bore area times distance
    /// to the nearest neighboring pocket or physical column endpoint.
    pub maximum_volume_m3: f64,
}

/// Accepted pocket coordinate, with signed dense-root residuals left visible.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctPocketState {
    pub port: PortId,
    pub volume_m3: f64,
    pub mode: DuctCavityMode,
}

/// Local observation and displacement on the accepted common trajectory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctPocketStep {
    pub port: PortId,
    pub mode: DuctCavityMode,
    pub volume_m3: f64,
    pub displaced_liquid_volume_m3: f64,
    pub end_absolute_pressure_pa: f64,
    pub end_inflow_m3_s: f64,
    pub rupture_events: usize,
    pub collapse_events: usize,
    pub last_rupture_time_s: Option<f64>,
    pub last_collapse_time_s: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct DuctPocketsStep {
    pub duct: DgStepOutput,
    pub pockets: Vec<DuctPocketStep>,
    /// Vapor displacement potential plus original gas barotropic work potential.
    /// This is not a gas thermal or latent-energy report.
    pub cavity_potential_j: f64,
    pub total_stored_energy_j: f64,
    pub energy_balance_residual_j: f64,
    /// Groups whose coupled solve changed more than one contact branch.
    pub simultaneous_groups: usize,
    pub accepted_substeps: usize,
    pub rejected_substeps: usize,
}

// Exact sparse rows sampled from the ORIGINAL native mass maps once at bind.
// Keeping exact zeros out changes no arithmetic coefficients or physical basis.
#[derive(Clone)]
struct ObservationRow(Vec<(usize, f64)>);
impl ObservationRow {
    fn value(&self, x: &[f64]) -> f64 {
        self.0.iter().map(|&(i, value)| value * x[i]).sum()
    }
}
struct Site {
    specification: DuctPocket,
    station_m: f64,
    pressure_node: usize,
    upstream_node: Option<usize>,
    upstream_area_m2: f64,
    downstream_area_m2: f64,
    upstream_flow: ObservationRow,
    gas: Option<TrappedGasPocket>,
}
struct ConstraintGain {
    active: u16,
    coefficients: Vec<f64>,
}
impl ConstraintGain {
    fn flow(&self, index: usize, x: &[f64]) -> f64 {
        -dot(
            &self.coefficients[index * x.len()..(index + 1) * x.len()],
            x,
        )
    }
}

/// Same-material liquid columns with simultaneous endpoint/interior pockets.
/// Every original declared port must have exactly one pocket; fixed native
/// closed, matched and pressure-release boundaries are retained unchanged.
pub struct CavitatingDuctPockets<F: CavitationProperties + FluidAcoustics + Clone> {
    duct: DgDuct<DgFluidMaterial<F>>,
    material: CavitationState,
    pub(super) limits: DuctCavityLimits,
    sites: Vec<Site>,
    has_gas: bool,
    gram: Vec<f64>,
    column_norms: Vec<f64>,
    pressure_rates: Vec<f64>,
    pressure_rows: Vec<ObservationRow>,
    state: DuctCavityState,
}

fn active(mask: u16, index: usize) -> bool {
    mask & (1 << index) != 0
}
fn mode(mask: u16, index: usize) -> DuctCavityMode {
    if active(mask, index) {
        DuctCavityMode::VaporCavity
    } else {
        DuctCavityMode::LiquidFilled
    }
}
fn mask(states: &[DuctPocketState]) -> u16 {
    states.iter().enumerate().fold(0, |mask, (i, state)| {
        if state.mode == DuctCavityMode::VaporCavity {
            mask | (1 << i)
        } else {
            mask
        }
    })
}

impl<F: CavitationProperties + FluidAcoustics + Clone> CavitatingDuctPockets<F> {
    pub fn new(
        duct: DgDuct<DgFluidMaterial<F>>,
        initial: DgState,
        time_s: f64,
        specifications: &[DuctPocket],
        limits: DuctCavityLimits,
    ) -> Result<Self, DuctCavityError> {
        Self::new_with_trapped_gas(duct, initial, time_s, specifications, &[], limits)
    }

    /// Bind optional fixed gas inventories without changing the original liquid
    /// owner. Exactly zero gas inventory retains the pure-vapor event path.
    /// Positive gas never closes and must stay above its declared minimum V.
    pub fn new_with_trapped_gas(
        duct: DgDuct<DgFluidMaterial<F>>,
        initial: DgState,
        time_s: f64,
        specifications: &[DuctPocket],
        gas_pockets: &[TrappedGasPocket],
        limits: DuctCavityLimits,
    ) -> Result<Self, DuctCavityError> {
        let material = duct.material().cavitation_state()?;
        if material.acoustic_state()? != duct.acoustic_material() {
            return Err(DuctCavityError::InvalidInput(
                "native acoustics and cavitation must retain the same material state",
            ));
        }
        let view = duct.continuous_operator();
        let x = view.energy_coordinates(&initial)?;
        let count = view.ports().len();
        if !(1..=MAX_DUCT_POCKETS).contains(&count)
            || specifications.len() != count
            || view.dimension() > MAX_CAVITATING_DUCT_DOFS
            || !time_s.is_finite()
            || time_s < 0.0
            || !limits.reference_pressure_pa.is_finite()
            || limits.reference_pressure_pa < material.saturation_pressure_pa()
            || !positive(limits.maximum_cavity_volume_m3)
            || !fraction(limits.maximum_liquid_strain)
            || !fraction(limits.maximum_interface_mach)
            || !fraction(limits.maximum_vapor_density_ratio)
            || x.iter().any(|value| !value.is_finite())
        {
            return Err(DuctCavityError::InvalidInput(
                "invalid native pocket count, clock, complete dimension or physical limits",
            ));
        }
        if material.saturated_vapor().is_some_and(|v| {
            v.density_kg_m3() / material.liquid_density_kg_m3() > limits.maximum_vapor_density_ratio
        }) {
            return Err(DuctCavityError::OutsideDomain(
                "vapor mass is not negligible relative to liquid",
            ));
        }
        if gas_pockets.len() > count {
            return Err(DuctCavityError::InvalidInput(
                "too many trapped-gas specifications",
            ));
        }
        for (i, gas) in gas_pockets.iter().enumerate() {
            if !positive(gas.minimum_volume_m3)
                || !view.ports().iter().any(|(id, _)| *id == gas.port)
                || gas_pockets[..i].iter().any(|other| other.port == gas.port)
            {
                return Err(DuctCavityError::InvalidInput(
                    "unknown/duplicate gas port or invalid minimum gas volume",
                ));
            }
        }
        let mut sites = Vec::with_capacity(count);
        let n_local = duct.mesh().order + 1;
        for (id, _) in view.ports() {
            let mut matching = specifications.iter().filter(|s| s.port == *id);
            let specification = *matching.next().ok_or(DuctCavityError::InvalidInput(
                "an original port has no vapor-pocket specification",
            ))?;
            if matching.next().is_some()
                || !positive(specification.maximum_volume_m3)
                || !specification.initial_volume_m3.is_finite()
                || specification.initial_volume_m3 < 0.0
                || specification.initial_volume_m3 > specification.maximum_volume_m3
            {
                return Err(DuctCavityError::InvalidInput(
                    "duplicate pocket or invalid local cavity volume",
                ));
            }
            let gas = gas_pockets.iter().find(|gas| gas.port == *id).copied();
            if gas.is_some_and(|g| {
                g.has_inventory()
                    && (specification.initial_volume_m3 < g.minimum_volume_m3
                        || g.minimum_volume_m3 >= specification.maximum_volume_m3)
            }) {
                return Err(DuctCavityError::InvalidInput(
                    "positive gas inventory requires an admitted nonzero initial volume",
                ));
            }
            let (station_m, pressure_node, upstream_node, upstream_area_m2, downstream_area_m2) =
                if duct.inlet() == DgBoundary::Port(*id) {
                    (0.0, 0, None, duct.inlet_area_m2(), duct.inlet_area_m2())
                } else if duct.outlet() == DgBoundary::Port(*id) {
                    (
                        duct.total_length(),
                        duct.field_dofs() - 1,
                        None,
                        duct.outlet_area_m2(),
                        duct.outlet_area_m2(),
                    )
                } else {
                    let geometry = duct.side_port_geometry(*id)?;
                    let node = geometry.port.downstream_element * n_local;
                    (
                        geometry.station_m,
                        node,
                        Some(node - 1),
                        geometry.upstream_area_m2,
                        geometry.downstream_area_m2,
                    )
                };
            sites.push(Site {
                specification,
                station_m,
                pressure_node,
                upstream_node,
                upstream_area_m2,
                downstream_area_m2,
                upstream_flow: ObservationRow(Vec::new()),
                gas,
            });
        }
        let total_volume: f64 = sites
            .iter()
            .map(|s| s.specification.initial_volume_m3)
            .sum();
        if !total_volume.is_finite() || total_volume > limits.maximum_cavity_volume_m3 {
            return Err(DuctCavityError::InvalidInput(
                "initial total cavity volume exceeds the declared limit",
            ));
        }
        for (i, site) in sites.iter().enumerate() {
            let mut span = duct.total_length();
            for boundary in [0.0, duct.total_length()] {
                let distance = (site.station_m - boundary).abs();
                if distance > 0.0 {
                    span = span.min(distance);
                }
            }
            for (j, other) in sites.iter().enumerate() {
                if i != j {
                    span = span.min((site.station_m - other.station_m).abs());
                }
            }
            if !positive(span)
                || site.specification.maximum_volume_m3
                    > 0.1 * site.upstream_area_m2.min(site.downstream_area_m2) * span
            {
                return Err(DuctCavityError::OutsideDomain(
                    "pocket volume is not small relative to its adjacent liquid columns",
                ));
            }
        }
        let n = view.dimension();
        let mut gram = vec![0.0; count * count];
        let mut column_norms = vec![0.0; count];
        let mut pressure_rates = vec![0.0; count * n];
        for i in 0..count {
            let column = &view.ports()[i].1;
            for j in 0..count {
                gram[i * count + j] = dot(column, &view.ports()[j].1);
            }
            column_norms[i] = gram[i * count + i].sqrt();
            if !positive(column_norms[i]) {
                return Err(DuctCavityError::ActiveConstraints(
                    "unrepresentable native pressure column",
                ));
            }
            for j in 0..n {
                pressure_rates[i * n + j] = (0..n)
                    .map(|r| column[r] * view.generator_per_s()[r * n + j])
                    .sum();
            }
        }
        if gram.iter().chain(&pressure_rates).any(|x| !x.is_finite()) {
            return Err(DuctCavityError::ActiveConstraints(
                "native coupled pressure coefficients overflow",
            ));
        }
        let mut pressure_rows: Vec<_> = (0..duct.field_dofs())
            .map(|_| ObservationRow(Vec::new()))
            .collect();
        let mut basis = vec![0.0; n];
        for column in 0..n {
            basis[column] = 1.0;
            let (pressures, flows) = duct.from_energy_coords(&basis);
            for (row, value) in pressure_rows.iter_mut().zip(pressures) {
                if value != 0.0 {
                    row.0.push((column, value));
                }
            }
            for site in &mut sites {
                if let Some(node) = site.upstream_node {
                    if flows[node] != 0.0 {
                        site.upstream_flow.0.push((column, flows[node]));
                    }
                }
            }
            basis[column] = 0.0;
        }
        let pocket_states: Vec<_> = sites
            .iter()
            .map(|s| DuctPocketState {
                port: s.specification.port,
                volume_m3: s.specification.initial_volume_m3,
                mode: DuctCavityMode::LiquidFilled,
            })
            .collect();
        let first = pocket_states[0];
        let has_gas = sites
            .iter()
            .any(|s| s.gas.is_some_and(TrappedGasPocket::has_inventory));
        let mut model = Self {
            duct,
            material,
            limits,
            sites,
            has_gas,
            gram,
            column_norms,
            pressure_rates,
            pressure_rows,
            state: DuctCavityState {
                duct: initial,
                time_s,
                volume_m3: first.volume_m3,
                mode: first.mode,
                pockets: pocket_states,
            },
        };
        let initial_x = model
            .duct
            .continuous_operator()
            .energy_coordinates(&model.state.duct)?;
        let volumes: Vec<_> = model.state.pockets.iter().map(|p| p.volume_m3).collect();
        let roundoff = 128.0 * f64::EPSILON * limits.reference_pressure_pa.max(1.0);
        let gain = model.resolve_contacts(initial_x, &volumes, 0.0, roundoff, &mut || false)?;
        let flows = if model.has_gas {
            model.gas_flows(
                gain.active,
                initial_x,
                1.0,
                &volumes,
                &mut GasFlowWorkspace::new(count)?,
            )?
        } else {
            let mut flows = [0.0; MAX_DUCT_POCKETS];
            for (i, flow) in flows[..count].iter_mut().enumerate() {
                *flow = gain.flow(i, initial_x);
            }
            flows
        };
        model.check_field(initial_x, 1.0, &flows[..count], roundoff)?;
        for (i, pocket) in model.state.pockets.iter_mut().enumerate() {
            pocket.mode = mode(gain.active, i);
        }
        model.state.mode = model.state.pockets[0].mode;
        Ok(model)
    }

    pub const fn duct(&self) -> &DgDuct<DgFluidMaterial<F>> {
        &self.duct
    }
    pub fn source(&self) -> &F {
        self.duct.material()
    }
    pub const fn material(&self) -> CavitationState {
        self.material
    }
    pub const fn state(&self) -> &DuctCavityState {
        &self.state
    }
    pub const fn limits(&self) -> DuctCavityLimits {
        self.limits
    }

    fn factor(
        &self,
        active_mask: u16,
        stiffness: &[f64],
    ) -> Result<(Vec<usize>, Option<Cholesky>), DuctCavityError> {
        let count = self.sites.len();
        let indices: Vec<_> = (0..count).filter(|&i| active(active_mask, i)).collect();
        if indices.is_empty() {
            return Ok((indices, None));
        }
        let mut normalized = vec![0.0; indices.len() * indices.len()];
        for (i, &row) in indices.iter().enumerate() {
            for (j, &column) in indices.iter().enumerate() {
                normalized[i * indices.len() + j] = self.gram[row * count + column]
                    / self.column_norms[row]
                    / self.column_norms[column];
            }
            normalized[i * indices.len() + i] +=
                stiffness[row] / self.column_norms[row] / self.column_norms[row];
        }
        let factor = cholesky(&normalized, indices.len()).map_err(|_| {
            DuctCavityError::ActiveConstraints(
                "dependent or indefinite native pressure constraints",
            )
        })?;
        Ok((indices, Some(factor)))
    }

    fn gain(&self, active_mask: u16) -> Result<ConstraintGain, DuctCavityError> {
        if self.has_gas {
            // Its volume-dependent coefficients are evaluated in the reused LU
            // workspace; pure-vapor models retain their cached constant gain.
            return Ok(ConstraintGain {
                active: active_mask,
                coefficients: Vec::new(),
            });
        }
        let n = self.duct.total_dofs();
        let count = self.sites.len();
        let (indices, factor) = self.factor(active_mask, &[0.0; MAX_DUCT_POCKETS])?;
        let mut coefficients = vec![0.0; count * n];
        if let Some(factor) = factor {
            let mut rhs = vec![0.0; indices.len()];
            for column in 0..n {
                for (j, &i) in indices.iter().enumerate() {
                    rhs[j] = self.pressure_rates[i * n + column] / self.column_norms[i];
                }
                factor.solve(&mut rhs);
                for (j, &i) in indices.iter().enumerate() {
                    coefficients[i * n + column] = rhs[j] / self.column_norms[i];
                }
            }
        }
        if coefficients.iter().any(|x| !x.is_finite()) {
            return Err(DuctCavityError::ActiveConstraints(
                "coupled constraint gain overflows",
            ));
        }
        Ok(ConstraintGain {
            active: active_mask,
            coefficients,
        })
    }

    // At a common root, contact branches solve q>=0, p'>=0, q*p'=0.
    // Positive-volume pockets remain mandatory constraints even when shrinking.
    // A strictly positive-definite Gram matrix makes the physical flow unique.
    // Multiple zero-flow branch labels are harmless; prefer the filled label.
    fn resolve_contacts(
        &self,
        x: &[f64],
        volumes: &[f64],
        volume_tolerance: f64,
        pressure_tolerance: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<ConstraintGain, DuctCavityError> {
        let count = self.sites.len();
        let n = x.len();
        let view = self.duct.continuous_operator();
        let mut mandatory = 0u16;
        let mut contacts = Vec::new();
        let mut rates = [0.0; MAX_DUCT_POCKETS];
        let stiffness = self.gas_stiffness(volumes)?;
        for i in 0..count {
            let gap = self.limits.reference_pressure_pa + dot(&view.ports()[i].1, x)
                - self.material.saturation_pressure_pa()
                - self.gas_response(i, volumes[i])?.pressure_pa;
            if !volumes[i].is_finite()
                || volumes[i] < -volume_tolerance
                || gap < -pressure_tolerance
            {
                return Err(DuctCavityError::RootResolution(
                    "joint pressure/volume contact is outside root tolerance",
                ));
            }
            if self.has_inventory(i) || volumes[i] > volume_tolerance {
                if gap.abs() > pressure_tolerance {
                    return Err(DuctCavityError::InvalidInput(
                        "nonempty pocket pressure disagrees with original vapor/gas pressure",
                    ));
                }
                mandatory |= 1 << i;
            } else if gap.abs() <= pressure_tolerance {
                contacts.push(i);
            }
            rates[i] = dot(&self.pressure_rates[i * n..(i + 1) * n], x);
        }
        let mut chosen: Option<(u16, [f64; MAX_DUCT_POCKETS], f64)> = None;
        for subset in 0..(1usize << contacts.len()) {
            if cancelled() {
                return Err(DuctCavityError::Cancelled);
            }
            let mut active_mask = mandatory;
            for (bit, &i) in contacts.iter().enumerate() {
                if subset & (1 << bit) != 0 {
                    active_mask |= 1 << i;
                }
            }
            let (indices, factor) = self.factor(active_mask, &stiffness)?;
            let mut flows = [0.0; MAX_DUCT_POCKETS];
            if let Some(factor) = factor {
                let mut rhs: Vec<_> = indices
                    .iter()
                    .map(|&i| -rates[i] / self.column_norms[i])
                    .collect();
                factor.solve(&mut rhs);
                for (j, &i) in indices.iter().enumerate() {
                    flows[i] = rhs[j] / self.column_norms[i];
                }
            }
            if flows[..count].iter().any(|q| !q.is_finite()) {
                return Err(DuctCavityError::ActiveConstraints(
                    "joint contact flow overflows",
                ));
            }
            let scale = (0..count)
                .map(|i| {
                    flows[i]
                        .abs()
                        .max((rates[i] / (self.gram[i * count + i] + stiffness[i])).abs())
                })
                .fold(1e-300, f64::max);
            let flow_roundoff = 512.0 * f64::EPSILON * scale;
            let admissible = contacts.iter().all(|&i| {
                if active(active_mask, i) {
                    flows[i] >= -flow_roundoff
                } else {
                    let mut rate = rates[i];
                    let mut magnitude = rates[i].abs();
                    for j in 0..count {
                        let term = self.gram[i * count + j] * flows[j];
                        rate += term;
                        magnitude += term.abs();
                    }
                    rate >= -512.0 * f64::EPSILON * magnitude.max(1e-300)
                }
            });
            if !admissible {
                continue;
            }
            if let Some((_, previous, tolerance)) = chosen {
                if (0..count)
                    .any(|i| (flows[i] - previous[i]).abs() > 4.0 * flow_roundoff.max(tolerance))
                {
                    return Err(DuctCavityError::ActiveConstraints(
                        "numerically ambiguous coupled contact flow",
                    ));
                }
            } else {
                chosen = Some((active_mask, flows, flow_roundoff));
            }
        }
        let (active_mask, _, _) = chosen.ok_or(DuctCavityError::ActiveConstraints(
            "no consistent coupled contact branch",
        ))?;
        self.gain(active_mask)
    }

    fn check_field(
        &self,
        x: &[f64],
        scale: f64,
        flows: &[f64],
        pressure_tolerance: f64,
    ) -> Result<(), DuctCavityError> {
        let bulk = self.duct.acoustic_material().adiabatic_bulk_modulus_pa();
        for (node, row) in self.pressure_rows.iter().enumerate() {
            let p = scale * row.value(x);
            if !p.is_finite() || (p / bulk).abs() > self.limits.maximum_liquid_strain {
                return Err(DuctCavityError::OutsideDomain(
                    "native liquid pressure exceeds declared small-strain limit",
                ));
            }
            // A designated face uses its original SHARED numerical trace.
            // The two discontinuous polynomial endpoint values are not two
            // independently imposed cavity pressures.
            let at_site = self
                .sites
                .iter()
                .any(|s| s.pressure_node == node || s.upstream_node == Some(node));
            if !at_site
                && self.limits.reference_pressure_pa + p
                    < self.material.saturation_pressure_pa() - pressure_tolerance
            {
                return Err(DuctCavityError::OutsideDomain(
                    "an undeclared liquid station reaches vapor pressure",
                ));
            }
        }
        let c = self.material.liquid_sound_speed_m_s();
        for (i, site) in self.sites.iter().enumerate() {
            let (left, right) = if site.upstream_node.is_some() {
                let left = scale * site.upstream_flow.value(x);
                (left, left + flows[i])
            } else {
                (flows[i], flows[i])
            };
            if !left.is_finite()
                || !right.is_finite()
                || (left / (site.upstream_area_m2 * c)).abs() > self.limits.maximum_interface_mach
                || (right / (site.downstream_area_m2 * c)).abs()
                    > self.limits.maximum_interface_mach
            {
                return Err(DuctCavityError::OutsideDomain(
                    "a liquid-interface trace exceeds the declared Mach limit",
                ));
            }
        }
        Ok(())
    }

    pub fn advance<Cancel: FnMut() -> bool>(
        &mut self,
        end_time_s: f64,
        numerics: &DuctCavityNumerics,
        cancelled: &mut Cancel,
    ) -> Result<DuctPocketsStep, DuctCavityError> {
        self.advance_kernel(end_time_s, numerics, cancelled)
    }

    fn advance_kernel(
        &mut self,
        end_time_s: f64,
        numerics: &DuctCavityNumerics,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<DuctPocketsStep, DuctCavityError> {
        numerics.admit()?;
        if cancelled() {
            return Err(DuctCavityError::Cancelled);
        }
        let dt = end_time_s - self.state.time_s;
        if !end_time_s.is_finite() || !positive(dt) {
            return Err(DuctCavityError::InvalidInput(
                "the physical endpoint must be finite and strictly later",
            ));
        }
        let view = self.duct.continuous_operator();
        let original = view.energy_coordinates(&self.state.duct)?;
        let n = original.len();
        let count = self.sites.len();
        let vapor = self.material.saturation_pressure_pa();
        let pressure_scale = self.limits.reference_pressure_pa.max(1.0);
        let energy_scale = (2.0 * self.state.duct.stored_energy_j())
            .max(pressure_scale * self.limits.maximum_cavity_volume_m3)
            .max(numerics.absolute_energy_tolerance_j);
        let coordinate_scale = energy_scale.sqrt();
        let impulse_scale = pressure_scale * dt;
        let mut atol = (numerics.absolute_energy_tolerance_j / energy_scale)
            .min(numerics.absolute_pressure_tolerance_pa / pressure_scale);
        for site in &self.sites {
            if numerics.absolute_volume_tolerance_m3 >= site.specification.maximum_volume_m3 {
                return Err(DuctCavityError::InvalidInput(
                    "volume tolerance cannot resolve a pocket",
                ));
            }
            atol = atol
                .min(numerics.absolute_volume_tolerance_m3 / site.specification.maximum_volume_m3);
        }
        if !positive(energy_scale) || !positive(impulse_scale) || !positive(atol) {
            return Err(DuctCavityError::InvalidInput(
                "pocket integration scales overflow",
            ));
        }
        let mut y: Vec<_> = original.iter().map(|x| x / coordinate_scale).collect();
        for (state, site) in self.state.pockets.iter().zip(&self.sites) {
            y.push(state.volume_m3 / site.specification.maximum_volume_m3);
        }
        y.resize(n + 3 * count + 1, 0.0);
        let mut integration =
            AdaptiveState::new(self.state.time_s, &y, numerics.maximum_step_s.min(dt));
        let mut gain = self.gain(mask(&self.state.pockets))?;
        let flow_workspace = if self.has_gas {
            Some(RefCell::new(GasFlowWorkspace::new(count)?))
        } else {
            None
        };
        let mut ruptures = [0usize; MAX_DUCT_POCKETS];
        let mut collapses = [0usize; MAX_DUCT_POCKETS];
        let mut last_ruptures = [None; MAX_DUCT_POCKETS];
        let mut last_collapses = [None; MAX_DUCT_POCKETS];
        let mut event_count = 0;
        let mut simultaneous_groups = 0;
        let specs: Vec<_> = self
            .sites
            .iter()
            .filter(|s| !s.gas.is_some_and(TrappedGasPocket::has_inventory))
            .map(|s| EventSpec {
                id: u64::from(s.specification.port.0),
                direction: EventDirection::Falling,
                initial: InitialEvent::Ignore,
            })
            .collect();
        let options = EventSetOptions {
            max_step: numerics.maximum_step_s,
            scan_substeps: numerics.event_scan_substeps,
            time_tolerance: numerics.event_time_tolerance_s,
            max_iterations: 128,
            // Only chooses the common sample within overlapping brackets.
            // Every contender is checked and the contact solve is simultaneous.
            order: EventOrder::LowestId,
        };
        let pi = PiController::default();
        while integration.t < end_time_s {
            if cancelled() {
                return Err(DuctCavityError::Cancelled);
            }
            let attempts = integration
                .accepted
                .checked_add(integration.rejected)
                .ok_or(DuctCavityError::AttemptLimit)?;
            if attempts >= numerics.maximum_attempts {
                return Err(DuctCavityError::AttemptLimit);
            }
            let failure = RefCell::new(None);
            let rhs = |_: f64, u: &[f64], out: &mut [f64]| {
                multiply(view.generator_per_s(), &u[..n], &mut out[..n]);
                let mut volumes = [0.0; MAX_DUCT_POCKETS];
                for (i, volume) in volumes[..count].iter_mut().enumerate() {
                    *volume = u[n + i] * self.sites[i].specification.maximum_volume_m3;
                }
                let flows = if let Some(workspace) = &flow_workspace {
                    match self.gas_flows(
                        gain.active,
                        &u[..n],
                        coordinate_scale,
                        &volumes[..count],
                        &mut workspace.borrow_mut(),
                    ) {
                        Ok(flows) => flows,
                        Err(error) => {
                            if failure.borrow().is_none() {
                                *failure.borrow_mut() = Some(error);
                            }
                            out.fill(f64::NAN);
                            return;
                        }
                    }
                } else {
                    let mut flows = [0.0; MAX_DUCT_POCKETS];
                    for (i, flow) in flows[..count].iter_mut().enumerate() {
                        *flow = coordinate_scale * gain.flow(i, &u[..n]);
                    }
                    flows
                };
                let domain = self
                    .check_field(
                        &u[..n],
                        coordinate_scale,
                        &flows[..count],
                        numerics.absolute_pressure_tolerance_pa,
                    )
                    .and_then(|()| {
                        let mut total = 0.0;
                        for (i, site) in self.sites.iter().enumerate() {
                            let volume = u[n + i] * site.specification.maximum_volume_m3;
                            total += volume.max(0.0);
                            if volume > site.specification.maximum_volume_m3 {
                                return Err(DuctCavityError::OutsideDomain(
                                    "a cavity exceeds its local fixed-column volume limit",
                                ));
                            }
                            let pressure = self.limits.reference_pressure_pa
                                + coordinate_scale * dot(&view.ports()[i].1, &u[..n]);
                            // Nonlinear RK stage points need not lie exactly
                            // on the gas manifold. Its residual is checked on
                            // accepted endpoints/events, never projected.
                            if active(gain.active, i)
                                && !self.has_inventory(i)
                                && (pressure - vapor).abs()
                                    > numerics.absolute_pressure_tolerance_pa
                            {
                                return Err(DuctCavityError::RootResolution(
                                    "an active shared pressure is outside vapor-pressure tolerance",
                                ));
                            }
                        }
                        if total > self.limits.maximum_cavity_volume_m3 {
                            return Err(DuctCavityError::OutsideDomain(
                                "total cavity volume exceeds the declared limit",
                            ));
                        }
                        Ok(())
                    });
                if let Err(error) = domain {
                    if failure.borrow().is_none() {
                        *failure.borrow_mut() = Some(error);
                    }
                    out.fill(f64::NAN);
                    return;
                }
                for i in 0..count {
                    for (value, b) in out[..n].iter_mut().zip(&view.ports()[i].1) {
                        *value += b * (flows[i] / coordinate_scale);
                    }
                    let pressure = coordinate_scale * dot(&view.ports()[i].1, &u[..n]);
                    out[n + i] = flows[i] / self.sites[i].specification.maximum_volume_m3;
                    out[n + count + i] = pressure / impulse_scale;
                    out[n + 2 * count + i] = pressure * flows[i] / energy_scale;
                }
                out[n + 3 * count] = quadratic_loss(view.generator_per_s(), &u[..n]);
            };
            let guard = |id: u64, _: f64, u: &[f64]| {
                let i = self
                    .sites
                    .iter()
                    .position(|s| u64::from(s.specification.port.0) == id)
                    .expect("native event ID");
                if active(gain.active, i) {
                    u[n + i]
                } else {
                    (self.limits.reference_pressure_pa
                        + coordinate_scale * dot(&view.ports()[i].1, &u[..n])
                        - vapor)
                        / pressure_scale
                }
            };
            let start = integration.t;
            let result = rk45_until_any_event(
                &mut integration,
                &rhs,
                &guard,
                &specs,
                end_time_s,
                numerics.relative_tolerance,
                atol,
                &pi,
                numerics.maximum_attempts - attempts,
                &options,
                &mut || cancelled(),
            );
            if let Some(error) = failure.into_inner() {
                return Err(error);
            }
            match result? {
                EventSetAdvance::Stopped(report) => match report.status {
                    AdaptiveStatus::ReachedEnd => break,
                    AdaptiveStatus::StepLimit => return Err(DuctCavityError::AttemptLimit),
                    AdaptiveStatus::Cancelled => return Err(DuctCavityError::Cancelled),
                },
                EventSetAdvance::Event(hit) => {
                    if integration.t <= start {
                        return Err(DuctCavityError::RootResolution(
                            "distinct simultaneous-contact time is unrepresentable",
                        ));
                    }
                    for contender in &hit.contenders {
                        let i = self
                            .sites
                            .iter()
                            .position(|s| u64::from(s.specification.port.0) == contender.id)
                            .expect("native event contender");
                        let residual = if active(gain.active, i) {
                            integration.u[n + i] * self.sites[i].specification.maximum_volume_m3
                                / numerics.absolute_volume_tolerance_m3
                        } else {
                            (self.limits.reference_pressure_pa
                                + coordinate_scale * dot(&view.ports()[i].1, &integration.u[..n])
                                - vapor)
                                / numerics.absolute_pressure_tolerance_pa
                        };
                        if residual.abs() > 1.0 {
                            return Err(DuctCavityError::RootResolution(
                                "overlapping roots do not share an admissible common contact state",
                            ));
                        }
                    }
                    let x: Vec<_> = integration.u[..n]
                        .iter()
                        .map(|x| x * coordinate_scale)
                        .collect();
                    let volumes: Vec<_> = self
                        .sites
                        .iter()
                        .enumerate()
                        .map(|(i, s)| integration.u[n + i] * s.specification.maximum_volume_m3)
                        .collect();
                    let next = self.resolve_contacts(
                        &x,
                        &volumes,
                        numerics.absolute_volume_tolerance_m3,
                        numerics.absolute_pressure_tolerance_pa,
                        cancelled,
                    )?;
                    let changed = (gain.active ^ next.active).count_ones() as usize;
                    if changed == 0 {
                        return Err(DuctCavityError::RootResolution(
                            "grazing contact has no admissible branch change",
                        ));
                    }
                    event_count += changed;
                    if event_count > numerics.maximum_events {
                        return Err(DuctCavityError::EventLimit);
                    }
                    if changed > 1 {
                        simultaneous_groups += 1;
                    }
                    for i in 0..count {
                        if active(gain.active, i) != active(next.active, i) {
                            if active(next.active, i) {
                                ruptures[i] += 1;
                                last_ruptures[i] = Some(integration.t);
                            } else {
                                collapses[i] += 1;
                                last_collapses[i] = Some(integration.t);
                            }
                        }
                    }
                    gain = next;
                }
            }
        }
        let endpoint: Vec<_> = integration.u[..n]
            .iter()
            .map(|x| x * coordinate_scale)
            .collect();
        let mut end_volumes = [0.0; MAX_DUCT_POCKETS];
        for (i, volume) in end_volumes[..count].iter_mut().enumerate() {
            *volume = integration.u[n + i] * self.sites[i].specification.maximum_volume_m3;
        }
        let end_flows = if let Some(workspace) = &flow_workspace {
            self.gas_flows(
                gain.active,
                &endpoint,
                1.0,
                &end_volumes[..count],
                &mut workspace.borrow_mut(),
            )?
        } else {
            let mut flows = [0.0; MAX_DUCT_POCKETS];
            for (i, flow) in flows[..count].iter_mut().enumerate() {
                *flow = gain.flow(i, &endpoint);
            }
            flows
        };
        let mut receipts = Vec::with_capacity(count);
        let mut pockets = Vec::with_capacity(count);
        let mut states = Vec::with_capacity(count);
        let mut potential = 0.0;
        let mut before = self.state.duct.stored_energy_j();
        for (i, site) in self.sites.iter().enumerate() {
            let volume = integration.u[n + i] * site.specification.maximum_volume_m3;
            let pressure = self.limits.reference_pressure_pa + dot(&view.ports()[i].1, &endpoint);
            let target_pressure = vapor + self.gas_response(i, volume)?.pressure_pa;
            if volume < -numerics.absolute_volume_tolerance_m3
                || volume > site.specification.maximum_volume_m3
                || (active(gain.active, i)
                    && (pressure - target_pressure).abs() > numerics.absolute_pressure_tolerance_pa)
                || (!active(gain.active, i)
                    && (volume.abs() > numerics.absolute_volume_tolerance_m3
                        || pressure < vapor - numerics.absolute_pressure_tolerance_pa))
            {
                return Err(DuctCavityError::RootResolution(
                    "a final pocket violates contact tolerances",
                ));
            }
            receipts.push(DgTrajectoryPort {
                id: site.specification.port,
                pressure_impulse_pa_s: integration.u[n + count + i] * impulse_scale,
                work_into_duct_j: integration.u[n + 2 * count + i] * energy_scale,
            });
            pockets.push(DuctPocketStep {
                port: site.specification.port,
                mode: mode(gain.active, i),
                volume_m3: volume,
                displaced_liquid_volume_m3: volume - self.state.pockets[i].volume_m3,
                end_absolute_pressure_pa: pressure,
                end_inflow_m3_s: end_flows[i],
                rupture_events: ruptures[i],
                collapse_events: collapses[i],
                last_rupture_time_s: last_ruptures[i],
                last_collapse_time_s: last_collapses[i],
            });
            states.push(DuctPocketState {
                port: site.specification.port,
                volume_m3: volume,
                mode: mode(gain.active, i),
            });
            potential += self.pocket_potential(i, volume)?;
            before += self.pocket_potential(i, self.state.pockets[i].volume_m3)?;
        }
        let loss = integration.u[n + 3 * count] * energy_scale;
        let (candidate, duct_report) = view.candidate_from_trajectory(
            &self.state.duct,
            endpoint,
            dt,
            &receipts,
            loss,
            numerics.absolute_energy_tolerance_j,
            numerics.relative_tolerance,
        )?;
        let total = candidate.stored_energy_j() + potential;
        let residual = total - before + loss;
        let tolerance = numerics.absolute_energy_tolerance_j
            + numerics.relative_tolerance * (total.abs() + before.abs() + loss.abs());
        if !residual.is_finite() || residual.abs() > tolerance {
            return Err(DuctCavityError::EnergyBalance);
        }
        let report = DuctPocketsStep {
            duct: duct_report,
            pockets,
            cavity_potential_j: potential,
            total_stored_energy_j: total,
            energy_balance_residual_j: residual,
            simultaneous_groups,
            accepted_substeps: integration.accepted,
            rejected_substeps: integration.rejected,
        };
        if cancelled() {
            return Err(DuctCavityError::Cancelled);
        }
        self.state = DuctCavityState {
            duct: candidate,
            time_s: end_time_s,
            volume_m3: states[0].volume_m3,
            mode: states[0].mode,
            pockets: states,
        };
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
