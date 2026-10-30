//! Finite third-body heating in the original enthalpy network.
//!
//! A physical layer is an EXISTING thermal cell between two existing body
//! contacts. Friction work enters that cell once. Its actual mass and caloric
//! law retain heat; conduction into the bodies follows the two original
//! geometries, material conductivities and source-resolved contact resistances.
//! The same trial layer temperature is supplied to the mechanical/material
//! callback. There is no second heat store, authored heat fraction or artificial
//! interface capacitance.
//!
//! This is an orthogonal, cell-centred finite-volume rung for a homogeneous
//! layer cell with volumetric heat generation. The cell centreline joins its
//! opposite face centres. Its two neighbours are distinct body cells; further
//! links incident on the layer refuse. Body networks may contain other cells.
//! Time and spatial uniformity require independent refinement. Loaded elastic
//! material changes still require their constitutive owner's thermoelastic
//! accounting; this route never relabels reversible coefficient work as heat.

use super::{FrictionHeatEntropy, entropy::InterfaceEntropy, finite, invalid};
use fs_math::det;
use fs_thermal::enthalpy::backward_euler::{
    BackwardEulerStep, BackwardEulerTrial, StepCoefficients,
};
use fs_thermal::enthalpy::{
    CellGeometry, CellId, ContactPatch, DepositionError, DepositionSimulation, EnergyAudit,
    MaterialId, series_contact_conductance_w_k,
};
use std::cell::RefCell;

/// Work-conjugate finite-thickness shear on the original rigid bodies.
pub mod rigid;

/// Original two-body/layer topology and finite material storage.
#[derive(Clone, Debug)]
pub struct ThirdBodyHeatRoute {
    /// Body A, body B, layer.
    cells: [CellId; 3],
    patches: [ContactPatch; 2],
    geometry: [CellGeometry; 3],
    material: [MaterialId; 3],
    mass_kg: [f64; 3],
    /// [body half-distance, layer half-distance] for each contact.
    distances_m: [[f64; 2]; 2],
}

impl ThirdBodyHeatRoute {
    /// Bind existing contacts in body A/B order and the original layer cell.
    /// Reference conductivities are [body A, body B, layer]. The original
    /// conductances must equal the derived bulk-plus-contact series values;
    /// a supplied contact resistance is therefore preserved, never discarded.
    pub fn from_links(
        thermal: &DepositionSimulation<'_>,
        layer: CellId,
        links: [u64; 2],
        conductivity_w_m_k: [f64; 3],
        area_resistance_m2_k_w: [f64; 2],
    ) -> Result<Self, DepositionError> {
        if links[0] == links[1] {
            return Err(invalid("third body needs two distinct physical contacts"));
        }
        let patches = [
            thermal
                .contact_patch(links[0])
                .ok_or(invalid("inactive third-body contact"))?,
            thermal
                .contact_patch(links[1])
                .ok_or(invalid("inactive third-body contact"))?,
        ];
        let mut cells = [layer; 3];
        for side in 0..2 {
            let patch = patches[side];
            if patch.covered_area_m2 != [patch.link.area_m2; 2] {
                return Err(invalid("third-body faces need matching physical areas"));
            }
            cells[side] = if patch.link.a == layer {
                patch.link.b
            } else if patch.link.b == layer {
                patch.link.a
            } else {
                return Err(invalid("both contacts must meet the original layer cell"));
            };
        }
        if cells[0] == cells[1] || cells[0] == layer || cells[1] == layer {
            return Err(invalid("third-body route needs three distinct cells"));
        }
        let original = cells.map(|id| &thermal.cells()[&id]);
        let geometry = original.map(|cell| cell.geometry);
        let (distance_a, face_a) =
            half_distances(patches[0].link.area_m2, geometry[0], geometry[2])?;
        let (distance_b, face_b) =
            half_distances(patches[1].link.area_m2, geometry[1], geometry[2])?;
        if face_a == face_b {
            return Err(invalid("body contacts must meet opposite layer faces"));
        }
        let route = Self {
            cells,
            patches,
            geometry,
            material: original.map(|cell| cell.material_id),
            mass_kg: original.map(|cell| cell.mass_kg),
            distances_m: [distance_a, distance_b],
        };
        route.validate(thermal)?;
        let conductances = route.conductances(conductivity_w_m_k, area_resistance_m2_k_w)?;
        for (side, expected) in conductances.into_iter().enumerate() {
            let actual = patches[side].link.conductance_w_k;
            if (actual - expected).abs() > 64.0 * f64::EPSILON * actual.abs().max(expected.abs()) {
                return Err(invalid(
                    "third-body link differs from its source material/geometry resistance",
                ));
            }
        }
        Ok(route)
    }

    /// Ordered original body destinations.
    pub const fn body_cells(&self) -> [CellId; 2] {
        [self.cells[0], self.cells[1]]
    }
    /// Original finite layer storage; no mass is minted by this route.
    pub const fn layer_cell(&self) -> CellId {
        self.cells[2]
    }
    /// [body half-length, layer half-length] in body A/B order.
    pub const fn distances_m(&self) -> [[f64; 2]; 2] {
        self.distances_m
    }

    fn validate(&self, thermal: &DepositionSimulation<'_>) -> Result<(), DepositionError> {
        for patch in self.patches {
            if thermal.contact_patch(patch.link.id) != Some(patch) {
                return Err(invalid("third-body route has stale contacts or coverage"));
            }
        }
        for (i, id) in self.cells.into_iter().enumerate() {
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(invalid("inactive third-body cell"))?;
            if cell.geometry != self.geometry[i]
                || cell.material_id != self.material[i]
                || cell.mass_kg != self.mass_kg[i]
            {
                return Err(invalid(
                    "third-body route changed original geometry, material or mass",
                ));
            }
        }
        for link in thermal.links().values() {
            if (link.a == self.layer_cell() || link.b == self.layer_cell())
                && !self.patches.iter().any(|patch| patch.link.id == link.id)
            {
                return Err(invalid(
                    "third-body cell has an unaccounted conductive face",
                ));
            }
        }
        Ok(())
    }

    fn conductances(
        &self,
        conductivity: [f64; 3],
        resistance: [f64; 2],
    ) -> Result<[f64; 2], DepositionError> {
        let mut result = [0.0; 2];
        for side in 0..2 {
            result[side] = series_contact_conductance_w_k(
                self.patches[side].link.area_m2,
                self.distances_m[side],
                [conductivity[side], conductivity[2]],
                resistance[side],
            )?;
        }
        // Unlike a massless interface, a finite insulated layer can store heat.
        Ok(result)
    }
}

fn half_distances(
    area: f64,
    body: CellGeometry,
    layer: CellGeometry,
) -> Result<([f64; 2], [f64; 3]), DepositionError> {
    let mut shared = None;
    for a in [body.start_m, body.end_m] {
        for b in [layer.start_m, layer.end_m] {
            if a == b && shared.replace(a).is_some() {
                return Err(invalid(
                    "third-body cells overlap instead of sharing one face",
                ));
            }
        }
    }
    let face = shared.ok_or(invalid("third-body cells need a shared face centre"))?;
    let mut distances = [0.0; 2];
    let mut normal = [[0.0; 3]; 2];
    for (i, geometry) in [body, layer].into_iter().enumerate() {
        let delta: [f64; 3] = std::array::from_fn(|j| geometry.end_m[j] - geometry.start_m[j]);
        let length = finite(det::hypot(det::hypot(delta[0], delta[1]), delta[2]))?;
        if length <= 0.0 || area > geometry.volume_m3 / length {
            return Err(invalid("third-body contact exceeds its prismatic face"));
        }
        distances[i] = 0.5 * length;
        for j in 0..3 {
            normal[i][j] = if face == geometry.end_m {
                delta[j] / length
            } else {
                -delta[j] / length
            };
        }
    }
    let cosine: f64 = normal[0]
        .into_iter()
        .zip(normal[1])
        .map(|(a, b)| a * b)
        .sum();
    if (cosine + 1.0).abs() > 256.0 * f64::EPSILON {
        return Err(invalid(
            "third-body route requires opposing orthogonal cell normals",
        ));
    }
    Ok((distances, face))
}

/// One original physical loss and material evaluation at the trial layer T.
pub struct ThirdBodyHeatInput<P> {
    /// Physical irreversible loss [J], not elastic storage or coefficient work.
    pub dissipated_j: f64,
    /// Original endpoint conductivity [body A, body B, layer], W/(m K).
    pub conductivity_w_m_k: [f64; 3],
    /// Original source-resolved area-specific contact resistances, m² K/W.
    pub area_resistance_m2_k_w: [f64; 2],
    /// OTHER sources and links. Neither owned contact may be overridden here.
    /// The original staged mechanical state belongs in the payload.
    pub transport: StepCoefficients<P>,
}

/// Reversible mechanical/caloric exchange around the same physical friction
/// input. The constitutive owner supplies both terms from its admitted potential;
/// neither is an independently authored heating fraction or material damping.
pub struct ThirdBodyThermoelasticInput<P> {
    pub heat: ThirdBodyHeatInput<P>,
    /// Signed transfer into the existing caloric cell [J].
    pub reversible_heat_j: f64,
    /// Actual change of the mechanical potential's entropy [J/K].
    pub mechanical_entropy_change_j_k: f64,
}

/// Independent thermal cell, whole-network, layer and combined joule gates.
#[derive(Clone, Copy, Debug)]
pub struct ThirdBodyHeatTolerances {
    pub cell_j: f64,
    pub thermal_j: f64,
    pub layer_j: f64,
    pub total_j: f64,
}

#[derive(Clone, Debug)]
pub struct ThirdBodyHeatReport {
    pub layer_temperature_k: f64,
    pub conductances_w_k: [f64; 2],
    /// Signed conduction into each body, including ordinary unequal-T exchange.
    pub body_heat_j: [f64; 2],
    pub dissipated_j: f64,
    /// Signed internal exchange from the mechanical potential, separate from D.
    pub reversible_heat_j: f64,
    /// Actual change in the EXISTING layer's mass times specific enthalpy.
    pub layer_stored_change_j: f64,
    /// Other heat reaching the layer from its original boundary and source.
    pub layer_external_heat_j: f64,
    /// Layer storage + body conduction - friction/reversible/other layer heat.
    pub layer_energy_residual_j: f64,
    pub external_heat_j: f64,
    pub thermal: EnergyAudit,
    /// Original caloric storage - outside heat - friction - reversible transfer.
    pub total_energy_residual_j: f64,
    /// Optional exact caloric storage and endpoint-quadrature entropy account.
    /// Includes the finite layer once; its contacts remain network conduction.
    pub entropy: Option<FrictionHeatEntropy>,
}

struct Endpoint<P> {
    payload: P,
    temperature: f64,
    conductances: [f64; 2],
    body_heat: [f64; 2],
    dissipated: f64,
    reversible: f64,
    mechanical_entropy_change: f64,
    layer_storage: f64,
    layer_external: f64,
    external: f64,
}

/// Residual over existing cell energies only; the finite layer owns its T.
pub struct ThirdBodyHeatStep<'s, 'm> {
    thermal: BackwardEulerStep<'s, 'm>,
    route: &'s ThirdBodyHeatRoute,
    initial_layer_h_j_kg: f64,
}

impl<'s, 'm> ThirdBodyHeatStep<'s, 'm> {
    pub fn new(
        thermal: &'s DepositionSimulation<'m>,
        target_s: f64,
        route: &'s ThirdBodyHeatRoute,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        route.validate(thermal)?;
        Ok(Self {
            thermal: BackwardEulerStep::new(thermal, target_s, checkpoint)?,
            route,
            initial_layer_h_j_kg: thermal.cells()[&route.layer_cell()].specific_h_j_kg,
        })
    }
    pub fn cell_ids(&self) -> &[CellId] {
        self.thermal.cell_ids()
    }
    pub fn duration_s(&self) -> f64 {
        self.thermal.duration_s()
    }

    /// Resolve mechanics and materials at the SAME tentative layer temperature.
    /// A caller can use `tangential_heat_candidate` with any original
    /// `TangentialLaw`. Every invocation must start from accepted history;
    /// return the candidate as payload and publish it only after `accept`.
    pub fn evaluate<P, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(&DepositionSimulation<'m>, f64) -> Result<ThirdBodyHeatInput<P>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ThirdBodyHeatTrial<'m, P>, E> {
        self.evaluate_entropy_option(energy_increments_j, |state, temperature| {
            Ok(ThirdBodyThermoelasticInput { heat: resolve(state, temperature)?,
                reversible_heat_j: 0.0, mechanical_entropy_change_j_k: 0.0 })
        }, false, checkpoint)
    }

    /// Audit the SAME finite-layer candidate without another material callback.
    /// Original caloric laws supply all storage, including layer latent heat.
    /// The original thermal links own bulk/contact conductive entropy; adding
    /// their body heat again as massless-interface production would double count.
    /// Physical friction produces D/T_layer on the endpoint quadrature. Its
    /// stored fraction need not leave the layer within this step. Independent
    /// sources remain unclassified, and temporal/algebraic error stays numerical.
    /// Every active law must provide its entropy integral. Use at a final trial
    /// to avoid entropy work inside the caller's nonlinear iterations.
    pub fn evaluate_with_entropy<P, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(&DepositionSimulation<'m>, f64) -> Result<ThirdBodyHeatInput<P>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ThirdBodyHeatTrial<'m, P>, E> {
        self.evaluate_entropy_option(energy_increments_j, |state, temperature| {
            Ok(ThirdBodyThermoelasticInput { heat: resolve(state, temperature)?,
                reversible_heat_j: 0.0, mechanical_entropy_change_j_k: 0.0 })
        }, true, checkpoint)
    }

    /// The same thermal residual with signed reciprocal heat from an admitted
    /// mechanical free energy. Physical loss remains nonnegative and is counted
    /// once. Reversible exchange is removed from external-source classification;
    /// its mechanical entropy is included separately when the audit is requested.
    pub fn evaluate_thermoelastic<P, E: From<DepositionError>>(
        &self, energy_increments_j: &[f64],
        resolve: impl FnOnce(&DepositionSimulation<'m>, f64)
            -> Result<ThirdBodyThermoelasticInput<P>, E>,
        audit_entropy: bool,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ThirdBodyHeatTrial<'m, P>, E> {
        self.evaluate_entropy_option(energy_increments_j, resolve, audit_entropy, checkpoint)
    }

    fn evaluate_entropy_option<P, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(&DepositionSimulation<'m>, f64) -> Result<ThirdBodyThermoelasticInput<P>, E>,
        with_entropy: bool,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ThirdBodyHeatTrial<'m, P>, E> {
        let checkpoint = RefCell::new(checkpoint);
        let check = || (*checkpoint.borrow_mut())();
        let dt = self.duration_s();
        let mut entropy = None;
        let coefficients =
            |state: &DepositionSimulation<'m>| -> Result<StepCoefficients<Endpoint<P>>, E> {
                let layer_id = self.route.layer_cell();
                let temperature = state.temperature_k(layer_id)?;
                let ThirdBodyThermoelasticInput { heat: input, reversible_heat_j,
                    mechanical_entropy_change_j_k } = resolve(state, temperature)?;
                finite(reversible_heat_j)?;
                finite(mechanical_entropy_change_j_k)?;
                if !input.dissipated_j.is_finite() || input.dissipated_j < 0.0 {
                    return Err(invalid(
                        "third-body heat needs finite nonnegative physical friction work",
                    )
                    .into());
                }
                let StepCoefficients {
                    mut sources_w,
                    mut conductances_w_k,
                    payload,
                } = input.transport;
                if with_entropy {
                    entropy = Some(InterfaceEntropy::new(state, &sources_w, dt, &check)?);
                }
                let conductances = self
                    .route
                    .conductances(input.conductivity_w_m_k, input.area_resistance_m2_k_w)?;
                let mut body_heat = [0.0; 2];
                for side in 0..2 {
                    if conductances_w_k
                        .insert(self.route.patches[side].link.id, conductances[side])
                        .is_some()
                    {
                        return Err(
                            invalid("third-body conduction has two competing owners").into()
                        );
                    }
                    body_heat[side] = finite(
                        dt * finite(
                            conductances[side]
                                * (temperature - state.temperature_k(self.route.cells[side])?),
                        )?,
                    )?;
                }
                let external = finite(
                    dt * state.external_step_power_w(&sources_w, &conductances_w_k, &check)?,
                )?;
                let layer = &state.cells()[&layer_id];
                let layer_storage = finite(
                    layer.mass_kg * finite(layer.specific_h_j_kg - self.initial_layer_h_j_kg)?,
                )?;
                // Independent local source account, using the original boundary's
                // exposed area (including both original contact coverages).
                let boundary = layer.boundary;
                let mut layer_power = finite(
                    finite(boundary.source_w + sources_w.get(&layer_id).copied().unwrap_or(0.0))?
                        + finite(state.exposed_area_m2(layer_id)? * boundary.htc_w_m2_k)?
                            * (boundary.ambient_k - temperature),
                )?;
                if let Some(reservoir) = boundary.reservoir {
                    layer_power = finite(
                        layer_power
                            + finite(
                                reservoir.conductance_w_k * (reservoir.temperature_k - temperature),
                            )?,
                    )?;
                }
                let layer_external = finite(dt * layer_power)?;
                let caloric_heat_j = finite(input.dissipated_j + reversible_heat_j)?;
                let power = finite(caloric_heat_j / dt)?;
                if caloric_heat_j != 0.0 && power == 0.0 {
                    return Err(DepositionError::Arithmetic.into());
                }
                let source = sources_w.entry(layer_id).or_insert(0.0);
                *source = finite(*source + power)?;
                Ok(StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload: Endpoint {
                        payload,
                        temperature,
                        conductances,
                        body_heat,
                        dissipated: input.dissipated_j,
                        reversible: reversible_heat_j,
                        mechanical_entropy_change: mechanical_entropy_change_j_k,
                        layer_storage,
                        layer_external,
                        external,
                    },
                })
            };
        let trial = if with_entropy {
            self.thermal.evaluate_with_entropy(energy_increments_j, coefficients, &check)?
        } else {
            self.thermal.evaluate(energy_increments_j, coefficients, &check)?
        };
        let endpoint = trial.payload();
        let entropy = entropy.map(|other_sources| {
            // This is a finite thermal CELL, not a massless interface. Its
            // energy residual already belongs to the original cell audit.
            FrictionHeatEntropy::from_trial(
                trial.entropy_audit().ok_or(DepositionError::Arithmetic)?,
                other_sources, endpoint.dissipated, endpoint.temperature, 0.0,
            )?.with_mechanical_storage(endpoint.mechanical_entropy_change,
                endpoint.reversible, endpoint.temperature)
        }).transpose()?;
        check()?;
        let report = ThirdBodyHeatReport {
            layer_temperature_k: endpoint.temperature,
            conductances_w_k: endpoint.conductances,
            body_heat_j: endpoint.body_heat,
            dissipated_j: endpoint.dissipated,
            reversible_heat_j: endpoint.reversible,
            layer_stored_change_j: endpoint.layer_storage,
            layer_external_heat_j: endpoint.layer_external,
            layer_energy_residual_j: finite(
                endpoint.layer_storage + finite(endpoint.body_heat[0] + endpoint.body_heat[1])?
                    - endpoint.dissipated
                    - endpoint.reversible
                    - endpoint.layer_external,
            )?,
            external_heat_j: endpoint.external,
            thermal: *trial.audit(),
            total_energy_residual_j: finite(
                trial.audit().stored_energy_change_j - endpoint.external - endpoint.dissipated
                    - endpoint.reversible,
            )?,
            entropy,
        };
        Ok(ThirdBodyHeatTrial { trial, report })
    }
}

pub struct ThirdBodyHeatTrial<'m, P> {
    trial: BackwardEulerTrial<'m, Endpoint<P>>,
    report: ThirdBodyHeatReport,
}
impl<'m, P> ThirdBodyHeatTrial<'m, P> {
    pub fn residual_j(&self) -> &[f64] {
        self.trial.residual_j()
    }
    pub fn payload(&self) -> &P {
        &self.trial.payload().payload
    }
    pub const fn report(&self) -> &ThirdBodyHeatReport {
        &self.report
    }
    pub fn accept(
        self,
        tolerance: ThirdBodyHeatTolerances,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, P, ThirdBodyHeatReport), DepositionError> {
        checkpoint()?;
        if [
            tolerance.cell_j,
            tolerance.thermal_j,
            tolerance.layer_j,
            tolerance.total_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid(
                "third-body tolerances must be finite nonnegative joules",
            ));
        }
        if self.report.layer_energy_residual_j.abs() > tolerance.layer_j {
            return Err(DepositionError::EnergyBalance(
                self.report.layer_energy_residual_j,
            ));
        }
        if self.report.total_energy_residual_j.abs() > tolerance.total_j {
            return Err(DepositionError::EnergyBalance(
                self.report.total_energy_residual_j,
            ));
        }
        let (thermal, endpoint, _) =
            self.trial
                .accept(tolerance.cell_j, tolerance.thermal_j, &mut checkpoint)?;
        checkpoint()?;
        Ok((thermal, endpoint.payload, self.report))
    }
}
