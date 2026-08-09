//! The lumped-network reduced transient rung: the cheap tier of the fidelity
//! graph, with an explicit validity gate.
//!
//! Bead `frankensim-extreal-program-f85xj.5.13`, the third item staged by
//! `f85xj.5.9`. The method-of-lines march in [`crate::transient`] is the
//! expensive answer; most early design questions ("does this die exceed its
//! limit during a thirty-second burst?") are answered by a small RC network in
//! microseconds.
//!
//! # Cost and authority are separate axes
//!
//! The point of a fidelity graph is that a cheaper model is not merely a
//! faster version of an expensive one — it is a WEAKER one, and both facts
//! have to be explicit. This rung is cheaper AND it is only valid in a stated
//! regime. It therefore refuses outside that regime rather than returning a
//! cheap wrong number, because a cheap wrong number is worse than no number:
//! it looks like an answer.
//!
//! # The validity gate is the Biot number
//!
//! Lumping asserts that a body is isothermal — that internal conduction is
//! fast compared with surface transfer. The Biot number `Bi = h Lc / k`
//! measures exactly that, and the `fs-vvreg` Level-A lumped row states the
//! admitted context as `Bi <= 0.1`. [`LumpedNode`] therefore carries the
//! quantities that determine Biot and [`BiotGate`] adjudicates them, so a
//! caller cannot construct a network that silently lies about its own
//! applicability.
//!
//! # One power declaration, two rungs
//!
//! Node power comes from the same [`crate::power::PowerMap`] the full rung
//! consumes — the caller supplies the per-node watts that map already
//! validated and audited. Two rungs describing the same hardware differently
//! is the failure mode a fidelity graph exists to prevent, so the reduced
//! rung deliberately has no power vocabulary of its own.

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_matdb::{MaterialCard, SelectionPolicy};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};

use crate::ConductionError;
use crate::material::ConductivityTable;
use crate::radiation::{STEFAN_BOLTZMANN_W_M2_K4, SurfaceEmissivity};

/// The admitted Biot ceiling for lumped treatment.
///
/// Matches the `fs-vvreg` Level-A `thermal-a-lumped-transient` row's declared
/// context. Duplicated as a constant rather than read at runtime because this
/// crate must not depend on the corpus registry; `tests/lumped.rs` asserts the
/// two agree, so the duplication cannot drift silently.
pub const LUMPED_BIOT_CEILING: f64 = 0.1;

/// Maximum nodes admitted in one reduced network.
pub const MAX_LUMPED_NODES: usize = 1_024;
/// Maximum accepted implicit enthalpy steps in one reduced march.
pub const MAX_LUMPED_ENTHALPY_STEPS: usize = 10_000_000;
/// Fixed bisection ceiling for each monotone backward-Euler enthalpy step.
pub const LUMPED_ENTHALPY_BISECTION_ITERATIONS: usize = 96;

#[derive(Debug, Clone, PartialEq)]
enum LumpedEmissivityModel {
    Declared(f64),
    MaterialCard(Vec<SurfaceEmissivity>),
}

/// Temperature-dependent transport data used by a lumped enthalpy body.
///
/// A declared transport says explicitly that its two constants have no
/// material-database provenance. A card-backed transport retains every
/// conductivity and emissivity query receipt over one shared temperature
/// grid, and refuses outside that grid rather than extrapolating toward a
/// phase boundary.
#[derive(Debug, Clone, PartialEq)]
pub struct LumpedThermalTransport {
    conductivity: ConductivityTable,
    emissivity: LumpedEmissivityModel,
    material_card_identity: Option<ContentHash>,
    minimum_conductivity_w_per_m_k: f64,
    maximum_emissivity: f64,
    identity: ContentHash,
}

impl LumpedThermalTransport {
    /// Construct an explicitly unproven constant transport model.
    pub fn try_declared(
        conductivity_w_per_m_k: f64,
        emissivity: f64,
    ) -> Result<Self, ConductionError> {
        if !(emissivity.is_finite() && (0.0..=1.0).contains(&emissivity)) {
            return Err(lumped_error(
                "declared thermal transport",
                format!("emissivity {emissivity} is outside [0,1]"),
            ));
        }
        let conductivity = ConductivityTable::declared(conductivity_w_per_m_k)?;
        let mut identity =
            DomainHasher::new("org.frankensim.fs-conduction.lumped-thermal-transport.v1");
        identity.update(&[0]);
        identity.update(&conductivity_w_per_m_k.to_bits().to_le_bytes());
        identity.update(&emissivity.to_bits().to_le_bytes());
        Ok(Self {
            conductivity,
            emissivity: LumpedEmissivityModel::Declared(emissivity),
            material_card_identity: None,
            minimum_conductivity_w_per_m_k: conductivity_w_per_m_k,
            maximum_emissivity: emissivity,
            identity: identity.finalize(),
        })
    }

    /// Resolve conductivity and hemispherical-total emissivity from one card.
    ///
    /// Both properties are queried at every strictly increasing temperature
    /// in `grid_k`. All receipts remain load-bearing. The grid must cover the
    /// entire phase curve later attached to the body.
    pub fn from_material_card(
        card: &MaterialCard,
        conductivity_property: &str,
        grid_k: &[f64],
        policy: SelectionPolicy,
    ) -> Result<Self, ConductionError> {
        let conductivity =
            ConductivityTable::from_claims(card.claims(), conductivity_property, grid_k, policy)?;
        let mut emissivity = Vec::new();
        emissivity.try_reserve_exact(grid_k.len()).map_err(|_| {
            lumped_error(
                "card-backed thermal transport",
                format!(
                    "retaining {} card-backed emissivity samples exceeded capacity",
                    grid_k.len()
                ),
            )
        })?;
        for &temperature_k in grid_k {
            emissivity.push(SurfaceEmissivity::from_card(
                "lumped-body whole boundary",
                card,
                temperature_k,
                policy,
            )?);
        }
        let minimum_conductivity_w_per_m_k = conductivity
            .knots()
            .iter()
            .map(|(_, value)| *value)
            .fold(f64::INFINITY, f64::min);
        let maximum_emissivity = emissivity
            .iter()
            .map(SurfaceEmissivity::value)
            .fold(0.0_f64, f64::max);
        let card_identity = card.content_hash();
        let mut identity =
            DomainHasher::new("org.frankensim.fs-conduction.lumped-thermal-transport.v1");
        identity.update(&[1]);
        identity.update(card_identity.as_bytes());
        identity.update(&(conductivity_property.len() as u64).to_le_bytes());
        identity.update(conductivity_property.as_bytes());
        identity.update(&(conductivity.knots().len() as u64).to_le_bytes());
        for ((temperature_k, conductivity_value), emissivity_sample) in
            conductivity.knots().iter().zip(&emissivity)
        {
            identity.update(&temperature_k.to_bits().to_le_bytes());
            identity.update(&conductivity_value.to_bits().to_le_bytes());
            identity.update(emissivity_sample.receipt().content_hash().as_bytes());
            identity.update(&emissivity_sample.value().to_bits().to_le_bytes());
        }
        for receipt in conductivity.receipts() {
            identity.update(receipt.content_hash().as_bytes());
        }
        Ok(Self {
            conductivity,
            emissivity: LumpedEmissivityModel::MaterialCard(emissivity),
            material_card_identity: Some(card_identity),
            minimum_conductivity_w_per_m_k,
            maximum_emissivity,
            identity: identity.finalize(),
        })
    }

    /// Card that supplied both properties, or `None` for declared constants.
    #[must_use]
    pub const fn material_card_identity(&self) -> Option<ContentHash> {
        self.material_card_identity
    }

    /// Complete transport-property and receipt identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    fn conductivity_at(&self, temperature_k: f64) -> Result<f64, ConductionError> {
        self.conductivity.eval(temperature_k)
    }

    fn emissivity_at(&self, temperature_k: f64) -> Result<f64, ConductionError> {
        match &self.emissivity {
            LumpedEmissivityModel::Declared(value) => Ok(*value),
            LumpedEmissivityModel::MaterialCard(samples) => {
                let low = samples
                    .first()
                    .expect("card-backed transport has a nonempty admitted grid")
                    .temperature_k();
                let high = samples
                    .last()
                    .expect("card-backed transport has a nonempty admitted grid")
                    .temperature_k();
                if !(low..=high).contains(&temperature_k) {
                    return Err(ConductionError::Radiation {
                        surface: "lumped-body whole boundary".to_owned(),
                        what: format!(
                            "temperature {temperature_k} K is outside sampled emissivity span [{low}, {high}] K"
                        ),
                        fix: "extend the material-card temperature grid; emissivity is never extrapolated"
                            .to_owned(),
                    });
                }
                match samples
                    .binary_search_by(|sample| sample.temperature_k().total_cmp(&temperature_k))
                {
                    Ok(index) => Ok(samples[index].value()),
                    Err(upper) => {
                        let lower = upper - 1;
                        let left = &samples[lower];
                        let right = &samples[upper];
                        let alpha = (temperature_k - left.temperature_k())
                            / (right.temperature_k() - left.temperature_k());
                        Ok(alpha.mul_add(right.value() - left.value(), left.value()))
                    }
                }
            }
        }
    }
}

/// One isothermal body whose constitutive coordinate is specific enthalpy.
///
/// Unlike [`LumpedNode`], this body does not assume a constant heat capacity.
/// Its temperature, density, and phase fraction come from an admitted
/// [`EquilibriumEnthalpyPhaseCurve`], allowing an isothermal latent-heat
/// interval without a fabricated effective heat capacity.
#[derive(Debug, Clone, PartialEq)]
pub struct LumpedEnthalpyBody<'a> {
    name: String,
    mass_kg: f64,
    surface_area_m2: f64,
    convection_w_per_m2_k: f64,
    characteristic_length_m: f64,
    transport: LumpedThermalTransport,
    phase_curve: &'a EquilibriumEnthalpyPhaseCurve,
    identity: ContentHash,
}

impl<'a> LumpedEnthalpyBody<'a> {
    /// Admit a uniform-body enthalpy model and the quantities needed to gate
    /// its isothermal assumption.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        name: impl Into<String>,
        mass_kg: f64,
        surface_area_m2: f64,
        convection_w_per_m2_k: f64,
        emissivity: f64,
        characteristic_length_m: f64,
        conductivity_w_per_m_k: f64,
        phase_curve: &'a EquilibriumEnthalpyPhaseCurve,
    ) -> Result<Self, ConductionError> {
        let transport = LumpedThermalTransport::try_declared(conductivity_w_per_m_k, emissivity)?;
        Self::try_new_with_transport(
            name,
            mass_kg,
            surface_area_m2,
            convection_w_per_m2_k,
            characteristic_length_m,
            transport,
            phase_curve,
        )
    }

    /// Admit a body with receipt-bearing temperature-dependent transport.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_transport(
        name: impl Into<String>,
        mass_kg: f64,
        surface_area_m2: f64,
        convection_w_per_m2_k: f64,
        characteristic_length_m: f64,
        transport: LumpedThermalTransport,
        phase_curve: &'a EquilibriumEnthalpyPhaseCurve,
    ) -> Result<Self, ConductionError> {
        let name = name.into();
        if name.trim().is_empty() {
            return Err(lumped_error(
                "<unnamed>",
                "enthalpy body has a blank name".to_owned(),
            ));
        }
        for (value, field) in [
            (mass_kg, "mass"),
            (surface_area_m2, "surface area"),
            (characteristic_length_m, "characteristic length"),
        ] {
            if !(value.is_finite() && value > 0.0) {
                return Err(lumped_error(
                    &name,
                    format!("{field} {value} is not finite and positive"),
                ));
            }
        }
        if !(convection_w_per_m2_k.is_finite() && convection_w_per_m2_k >= 0.0) {
            return Err(lumped_error(
                &name,
                format!(
                    "convection coefficient {convection_w_per_m2_k} W/(m2 K) is not finite and non-negative"
                ),
            ));
        }
        if let Some(card_identity) = transport.material_card_identity()
            && card_identity != phase_curve.material_card_identity()
        {
            return Err(lumped_error(
                &name,
                "thermal transport and equilibrium phase curve come from different material cards"
                    .to_owned(),
            ));
        }
        for knot in phase_curve.knots() {
            transport.conductivity_at(knot.temperature_k)?;
            transport.emissivity_at(knot.temperature_k)?;
        }
        let mut identity =
            DomainHasher::new("org.frankensim.fs-conduction.lumped-enthalpy-body.v2");
        identity.update(&(name.len() as u64).to_le_bytes());
        identity.update(name.as_bytes());
        for value in [
            mass_kg,
            surface_area_m2,
            convection_w_per_m2_k,
            characteristic_length_m,
        ] {
            identity.update(&value.to_bits().to_le_bytes());
        }
        identity.update(transport.identity().as_bytes());
        identity.update(phase_curve.identity().as_bytes());
        Ok(Self {
            name,
            mass_kg,
            surface_area_m2,
            convection_w_per_m2_k,
            characteristic_length_m,
            transport,
            phase_curve,
            identity: identity.finalize(),
        })
    }

    /// Complete body, transfer, and phase-curve identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Invariant material mass advanced by the specific-enthalpy balance [kg].
    #[must_use]
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }

    /// Complete external area participating in convection and radiation [m2].
    #[must_use]
    pub const fn surface_area_m2(&self) -> f64 {
        self.surface_area_m2
    }

    /// Volume-to-area length used by the Biot applicability gate [m].
    #[must_use]
    pub const fn characteristic_length_m(&self) -> f64 {
        self.characteristic_length_m
    }

    /// Temperature-dependent conductivity/emissivity authority.
    #[must_use]
    pub const fn transport(&self) -> &LumpedThermalTransport {
        &self.transport
    }

    /// Equilibrium phase curve that closes temperature from enthalpy.
    #[must_use]
    pub const fn phase_curve(&self) -> &EquilibriumEnthalpyPhaseCurve {
        self.phase_curve
    }

    /// Conservative Biot number including convection and the maximum
    /// linearized radiative transfer coefficient over the admitted curve and
    /// declared ambient.
    #[must_use]
    pub fn maximum_biot(&self, ambient_temperature_k: f64) -> f64 {
        let maximum_curve_temperature = self
            .phase_curve
            .knots()
            .iter()
            .map(|knot| knot.temperature_k)
            .fold(0.0_f64, f64::max);
        let maximum_temperature = maximum_curve_temperature.max(ambient_temperature_k);
        let radiation = 4.0
            * self.transport.maximum_emissivity
            * STEFAN_BOLTZMANN_W_M2_K4
            * maximum_temperature
            * maximum_temperature
            * maximum_temperature;
        (self.convection_w_per_m2_k + radiation) * self.characteristic_length_m
            / self.transport.minimum_conductivity_w_per_m_k
    }
}

/// Controls for one constant-environment uniform enthalpy march.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LumpedEnthalpyMarchConfig {
    /// Initial specific enthalpy [J/kg].
    pub initial_specific_enthalpy_j_kg: f64,
    /// Constant environmental radiation/convection temperature [K].
    pub ambient_temperature_k: f64,
    /// Constant internally deposited power [W], positive into the body.
    pub internal_power_w: f64,
    /// Requested physical horizon [s].
    pub duration_s: f64,
    /// Maximum backward-Euler step [s]; the final step is shortened exactly.
    pub maximum_step_s: f64,
    /// Maximum accepted step count.
    pub maximum_steps: usize,
    /// Absolute bisection interval tolerance in specific enthalpy [J/kg].
    pub enthalpy_tolerance_j_kg: f64,
}

/// One accepted state and the discrete energy ledger for the preceding step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LumpedEnthalpySample {
    /// Physical time [s].
    pub time_s: f64,
    /// Equilibrium thermodynamic state at this boundary.
    pub phase_state: EquilibriumPhaseState,
    /// Convective power into the body [W].
    pub convection_into_body_w: f64,
    /// Net radiative power into the body [W].
    pub radiation_into_body_w: f64,
    /// Declared internal power [W].
    pub internal_power_w: f64,
    /// Total power into the body [W].
    pub net_power_into_body_w: f64,
    /// Backward-Euler energy closure `m delta h - dt P_end` [J].
    pub step_energy_residual_j: f64,
}

/// Complete admitted reduced enthalpy trajectory.
#[derive(Debug, Clone, PartialEq)]
pub struct LumpedEnthalpyMarch {
    /// Exact admitted body that supplied mass, area, transport, and phase data.
    body_identity: ContentHash,
    /// Conservative Biot number that admitted the isothermal rung.
    maximum_biot: f64,
    /// Ordered initial and accepted endpoint samples.
    samples: Vec<LumpedEnthalpySample>,
    /// Sum of absolute per-step discrete energy residuals [J].
    cumulative_absolute_energy_residual_j: f64,
    /// Identity binding body, controls, and every accepted phase state.
    identity: ContentHash,
}

impl LumpedEnthalpyMarch {
    /// Exact admitted body that supplied mass, area, transport, and phase data.
    #[must_use]
    pub const fn body_identity(&self) -> ContentHash {
        self.body_identity
    }

    /// Conservative Biot number that admitted the isothermal rung.
    #[must_use]
    pub const fn maximum_biot(&self) -> f64 {
        self.maximum_biot
    }

    /// Ordered initial and accepted endpoint samples.
    #[must_use]
    pub fn samples(&self) -> &[LumpedEnthalpySample] {
        &self.samples
    }

    /// Sum of absolute per-step discrete energy residuals [J].
    #[must_use]
    pub const fn cumulative_absolute_energy_residual_j(&self) -> f64 {
        self.cumulative_absolute_energy_residual_j
    }

    /// Identity binding body, controls, and every accepted phase state.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// March a Biot-admitted isothermal body in specific enthalpy.
///
/// Each step solves the monotone backward-Euler balance by deterministic
/// bisection over the phase curve's admitted enthalpy interval. Exhausting
/// that interval refuses rather than extrapolating a solid or liquid law.
pub fn solve_lumped_enthalpy(
    cx: &Cx<'_>,
    body: &LumpedEnthalpyBody<'_>,
    gate: BiotGate,
    config: LumpedEnthalpyMarchConfig,
) -> Result<LumpedEnthalpyMarch, ConductionError> {
    validate_enthalpy_config(body, config)?;
    let maximum_biot = body.maximum_biot(config.ambient_temperature_k);
    if maximum_biot > gate.ceiling() {
        return Err(lumped_error(
            &body.name,
            format!(
                "convection-plus-radiation Biot {maximum_biot} exceeds ceiling {}; escalate to spatial transient enthalpy transport",
                gate.ceiling()
            ),
        ));
    }
    let requested_steps_f64 = (config.duration_s / config.maximum_step_s).ceil();
    if !requested_steps_f64.is_finite()
        || requested_steps_f64 > config.maximum_steps as f64
        || requested_steps_f64 > MAX_LUMPED_ENTHALPY_STEPS as f64
    {
        return Err(lumped_error(
            &body.name,
            format!(
                "enthalpy march requests {requested_steps_f64} steps, above caller maximum {} or hard maximum {MAX_LUMPED_ENTHALPY_STEPS}",
                config.maximum_steps,
            ),
        ));
    }
    let requested_steps = requested_steps_f64 as usize;
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "lumped-enthalpy",
        at: 0,
    })?;
    let initial = phase_state(body, config.initial_specific_enthalpy_j_kg)?;
    let (convection, radiation, net) = enthalpy_power(body, config, initial.temperature_k())?;
    let mut samples = Vec::with_capacity(requested_steps.saturating_add(1));
    samples.push(LumpedEnthalpySample {
        time_s: 0.0,
        phase_state: initial,
        convection_into_body_w: convection,
        radiation_into_body_w: radiation,
        internal_power_w: config.internal_power_w,
        net_power_into_body_w: net,
        step_energy_residual_j: 0.0,
    });
    let mut time_s = 0.0;
    let mut specific_enthalpy = config.initial_specific_enthalpy_j_kg;
    let mut cumulative_absolute_energy_residual_j = 0.0;
    for step in 0..requested_steps {
        cx.checkpoint().map_err(|_| ConductionError::Cancelled {
            stage: "lumped-enthalpy",
            at: step,
        })?;
        let dt_s = config.maximum_step_s.min(config.duration_s - time_s);
        let next = solve_enthalpy_step(body, config, specific_enthalpy, dt_s)?;
        let state = phase_state(body, next)?;
        let (convection, radiation, net) = enthalpy_power(body, config, state.temperature_k())?;
        let residual = body.mass_kg * (next - specific_enthalpy) - dt_s * net;
        cumulative_absolute_energy_residual_j += residual.abs();
        time_s += dt_s;
        specific_enthalpy = next;
        samples.push(LumpedEnthalpySample {
            time_s,
            phase_state: state,
            convection_into_body_w: convection,
            radiation_into_body_w: radiation,
            internal_power_w: config.internal_power_w,
            net_power_into_body_w: net,
            step_energy_residual_j: residual,
        });
    }
    let mut identity = DomainHasher::new("org.frankensim.fs-conduction.lumped-enthalpy-march.v1");
    identity.update(body.identity().as_bytes());
    for value in [
        config.initial_specific_enthalpy_j_kg,
        config.ambient_temperature_k,
        config.internal_power_w,
        config.duration_s,
        config.maximum_step_s,
        config.enthalpy_tolerance_j_kg,
        maximum_biot,
    ] {
        identity.update(&value.to_bits().to_le_bytes());
    }
    identity.update(&(config.maximum_steps as u64).to_le_bytes());
    for sample in &samples {
        identity.update(&sample.time_s.to_bits().to_le_bytes());
        identity.update(sample.phase_state.identity().as_bytes());
        identity.update(&sample.step_energy_residual_j.to_bits().to_le_bytes());
    }
    Ok(LumpedEnthalpyMarch {
        body_identity: body.identity(),
        maximum_biot,
        samples,
        cumulative_absolute_energy_residual_j,
        identity: identity.finalize(),
    })
}

fn validate_enthalpy_config(
    body: &LumpedEnthalpyBody<'_>,
    config: LumpedEnthalpyMarchConfig,
) -> Result<(), ConductionError> {
    for (value, field, strictly_positive) in [
        (config.ambient_temperature_k, "ambient temperature", true),
        (config.duration_s, "duration", false),
        (config.maximum_step_s, "maximum step", true),
        (config.enthalpy_tolerance_j_kg, "enthalpy tolerance", true),
    ] {
        if !value.is_finite()
            || (strictly_positive && value <= 0.0)
            || (!strictly_positive && value < 0.0)
        {
            return Err(lumped_error(
                &body.name,
                format!("{field} {value} is outside its finite admitted domain"),
            ));
        }
    }
    if !config.initial_specific_enthalpy_j_kg.is_finite() || !config.internal_power_w.is_finite() {
        return Err(lumped_error(
            &body.name,
            "initial enthalpy and internal power must be finite".to_owned(),
        ));
    }
    if config.maximum_steps == 0 {
        return Err(lumped_error(
            &body.name,
            "maximum_steps must be positive".to_owned(),
        ));
    }
    Ok(())
}

fn solve_enthalpy_step(
    body: &LumpedEnthalpyBody<'_>,
    config: LumpedEnthalpyMarchConfig,
    previous_h: f64,
    dt_s: f64,
) -> Result<f64, ConductionError> {
    let knots = body.phase_curve.knots();
    let mut lower = knots[0].specific_enthalpy_j_kg;
    let mut upper = knots[knots.len() - 1].specific_enthalpy_j_kg;
    let residual = |specific_enthalpy: f64| -> Result<f64, ConductionError> {
        let state = phase_state(body, specific_enthalpy)?;
        let (_, _, power) = enthalpy_power(body, config, state.temperature_k())?;
        Ok(body.mass_kg * (specific_enthalpy - previous_h) - dt_s * power)
    };
    let lower_residual = residual(lower)?;
    let upper_residual = residual(upper)?;
    if lower_residual > 0.0 || upper_residual < 0.0 {
        return Err(lumped_error(
            &body.name,
            format!(
                "one implicit step leaves the admitted phase-curve enthalpy interval [{lower}, {upper}] J/kg; extend evidence-backed phase data or reduce the step/environment"
            ),
        ));
    }
    for _ in 0..LUMPED_ENTHALPY_BISECTION_ITERATIONS {
        let midpoint = lower + 0.5 * (upper - lower);
        if midpoint.to_bits() == lower.to_bits() || midpoint.to_bits() == upper.to_bits() {
            break;
        }
        let midpoint_residual = residual(midpoint)?;
        if midpoint_residual <= 0.0 {
            lower = midpoint;
        } else {
            upper = midpoint;
        }
        if upper - lower <= config.enthalpy_tolerance_j_kg {
            break;
        }
    }
    Ok(lower + 0.5 * (upper - lower))
}

fn phase_state(
    body: &LumpedEnthalpyBody<'_>,
    specific_enthalpy_j_kg: f64,
) -> Result<EquilibriumPhaseState, ConductionError> {
    body.phase_curve
        .state_at_specific_enthalpy(specific_enthalpy_j_kg)
        .map_err(|error| lumped_error(&body.name, error.to_string()))
}

fn enthalpy_power(
    body: &LumpedEnthalpyBody<'_>,
    config: LumpedEnthalpyMarchConfig,
    body_temperature_k: f64,
) -> Result<(f64, f64, f64), ConductionError> {
    let convection = body.convection_w_per_m2_k
        * body.surface_area_m2
        * (config.ambient_temperature_k - body_temperature_k);
    let radiation = body.transport.emissivity_at(body_temperature_k)?
        * STEFAN_BOLTZMANN_W_M2_K4
        * body.surface_area_m2
        * (fourth_power(config.ambient_temperature_k) - fourth_power(body_temperature_k));
    let net = config.internal_power_w + convection + radiation;
    Ok((convection, radiation, net))
}

fn fourth_power(value: f64) -> f64 {
    let square = value * value;
    square * square
}

/// One lumped node: an isothermal body with a surface path to ambient.
#[derive(Debug, Clone, PartialEq)]
pub struct LumpedNode {
    name: String,
    capacitance_j_per_k: f64,
    conductance_w_per_k: f64,
    characteristic_length_m: f64,
    conductivity_w_per_m_k: f64,
    surface_area_m2: f64,
}

impl LumpedNode {
    /// Declare a node from the quantities that determine both its dynamics
    /// and its validity.
    ///
    /// `conductance_w_per_k` is the surface path to ambient, `h·A`. The
    /// characteristic length, conductivity and area are not redundant with
    /// it: they are what make the Biot number computable, and a node that
    /// could not state its own Biot could not be gated.
    ///
    /// # Errors
    /// [`ConductionError::ScenarioRow`] for a blank name or any non-positive
    /// or non-finite quantity.
    pub fn new(
        name: impl Into<String>,
        capacitance_j_per_k: f64,
        conductance_w_per_k: f64,
        characteristic_length_m: f64,
        conductivity_w_per_m_k: f64,
        surface_area_m2: f64,
    ) -> Result<Self, ConductionError> {
        let name = name.into();
        if name.trim().is_empty() {
            return Err(lumped_error(
                "<unnamed>",
                "lumped node has a blank name".to_string(),
            ));
        }
        for (value, field) in [
            (capacitance_j_per_k, "capacitance"),
            (conductance_w_per_k, "surface conductance"),
            (characteristic_length_m, "characteristic length"),
            (conductivity_w_per_m_k, "conductivity"),
            (surface_area_m2, "surface area"),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(lumped_error(
                    &name,
                    format!("{field} {value} is not finite and positive"),
                ));
            }
        }
        Ok(Self {
            name,
            capacitance_j_per_k,
            conductance_w_per_k,
            characteristic_length_m,
            conductivity_w_per_m_k,
            surface_area_m2,
        })
    }

    /// Node name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Thermal capacitance, J/K.
    #[must_use]
    pub const fn capacitance_j_per_k(&self) -> f64 {
        self.capacitance_j_per_k
    }

    /// Surface conductance to ambient, W/K.
    #[must_use]
    pub const fn conductance_w_per_k(&self) -> f64 {
        self.conductance_w_per_k
    }

    /// Effective surface transfer coefficient `h = (hA) / A`, W/(m²·K).
    #[must_use]
    pub fn transfer_coefficient_w_per_m2_k(&self) -> f64 {
        self.conductance_w_per_k / self.surface_area_m2
    }

    /// The node's Biot number, `h Lc / k`.
    ///
    /// This is the quantity that decides whether lumping is admissible at
    /// all, so it is derived from declared inputs rather than accepted as one.
    #[must_use]
    pub fn biot(&self) -> f64 {
        self.transfer_coefficient_w_per_m2_k() * self.characteristic_length_m
            / self.conductivity_w_per_m_k
    }

    /// The node's time constant `C / (hA)`, s.
    #[must_use]
    pub fn time_constant_s(&self) -> f64 {
        self.capacitance_j_per_k / self.conductance_w_per_k
    }
}

/// The verdict of the lumped validity gate.
#[derive(Debug, Clone, PartialEq)]
pub enum ValidityVerdict {
    /// Every node is inside the declared Biot ceiling.
    Admitted {
        /// Largest Biot number in the network.
        worst_biot: f64,
    },
    /// At least one node is outside it. The rung REFUSES rather than
    /// returning a cheap wrong number.
    Refused {
        /// The offending node.
        node: String,
        /// Its Biot number.
        biot: f64,
        /// The ceiling it exceeded.
        ceiling: f64,
    },
}

impl ValidityVerdict {
    /// Whether the network may be solved on this rung.
    #[must_use]
    pub const fn admitted(&self) -> bool {
        matches!(self, Self::Admitted { .. })
    }
}

/// The Biot gate, with its ceiling declared rather than assumed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BiotGate {
    ceiling: f64,
}

impl BiotGate {
    /// The gate at the corpus-declared ceiling.
    #[must_use]
    pub const fn corpus_default() -> Self {
        Self {
            ceiling: LUMPED_BIOT_CEILING,
        }
    }

    /// A gate at a caller-declared ceiling.
    ///
    /// Loosening this is a modelling decision the caller owns; the returned
    /// verdict always reports the ceiling actually applied, so a loosened
    /// gate cannot be mistaken for the corpus one.
    ///
    /// # Errors
    /// A non-positive or non-finite ceiling.
    pub fn at(ceiling: f64) -> Result<Self, ConductionError> {
        if !ceiling.is_finite() || ceiling <= 0.0 {
            return Err(lumped_error(
                "biot gate",
                format!("ceiling {ceiling} is not finite and positive"),
            ));
        }
        Ok(Self { ceiling })
    }

    /// The applied ceiling.
    #[must_use]
    pub const fn ceiling(self) -> f64 {
        self.ceiling
    }

    /// Adjudicate a network.
    #[must_use]
    pub fn adjudicate(self, network: &LumpedNetwork) -> ValidityVerdict {
        let mut worst = 0.0f64;
        for node in &network.nodes {
            let biot = node.biot();
            if biot > self.ceiling {
                return ValidityVerdict::Refused {
                    node: node.name.clone(),
                    biot,
                    ceiling: self.ceiling,
                };
            }
            worst = worst.max(biot);
        }
        ValidityVerdict::Admitted { worst_biot: worst }
    }
}

/// A reduced thermal network: isothermal nodes coupled to one ambient.
#[derive(Debug, Clone, PartialEq)]
pub struct LumpedNetwork {
    nodes: Vec<LumpedNode>,
    ambient_k: f64,
}

impl LumpedNetwork {
    /// Admit a network.
    ///
    /// # Errors
    /// [`ConductionError::ScenarioRow`] for an empty or oversized node set, a
    /// duplicated node name, or a non-finite ambient.
    pub fn new(nodes: Vec<LumpedNode>, ambient_k: f64) -> Result<Self, ConductionError> {
        if nodes.is_empty() {
            return Err(lumped_error(
                "lumped network",
                "declares no nodes".to_string(),
            ));
        }
        if nodes.len() > MAX_LUMPED_NODES {
            return Err(lumped_error(
                "lumped network",
                format!(
                    "declares {} nodes, above the admitted maximum {MAX_LUMPED_NODES}",
                    nodes.len()
                ),
            ));
        }
        if !ambient_k.is_finite() {
            return Err(lumped_error(
                "lumped network",
                format!("ambient temperature {ambient_k} K is not finite"),
            ));
        }
        let mut nodes = nodes;
        nodes.sort_by(|a, b| a.name.cmp(&b.name));
        for pair in nodes.windows(2) {
            if pair[0].name == pair[1].name {
                return Err(lumped_error(
                    &pair[0].name,
                    "lumped node name is duplicated".to_string(),
                ));
            }
        }
        Ok(Self { nodes, ambient_k })
    }

    /// Nodes in deterministic name order.
    #[must_use]
    pub fn nodes(&self) -> &[LumpedNode] {
        &self.nodes
    }

    /// Declared ambient temperature, K.
    #[must_use]
    pub const fn ambient_k(&self) -> f64 {
        self.ambient_k
    }

    /// Steady temperature of every node under a constant power vector, K.
    ///
    /// Nodes are uncoupled here (each has its own path to one ambient), so
    /// the steady state is `T_amb + P/(hA)` node by node.
    ///
    /// # Errors
    /// A power vector of the wrong length, or a non-finite entry.
    pub fn steady(&self, power_w: &[f64]) -> Result<Vec<f64>, ConductionError> {
        self.admit_power(power_w)?;
        Ok(self
            .nodes
            .iter()
            .zip(power_w.iter())
            .map(|(node, power)| self.ambient_k + power / node.conductance_w_per_k)
            .collect())
    }

    /// The closed-form response at time `t` from an initial temperature.
    ///
    /// Each node is a first-order system, so the reduced rung has an ANALYTIC
    /// solution — no time stepping and no step-size choice. That is exactly
    /// what makes it the cheap tier, and it also means the rung introduces no
    /// integration error of its own to confuse with model error.
    ///
    /// # Errors
    /// A mismatched power or initial vector, a non-finite entry, or a
    /// negative time.
    pub fn response_at(
        &self,
        power_w: &[f64],
        initial_k: &[f64],
        time_s: f64,
    ) -> Result<Vec<f64>, ConductionError> {
        self.admit_power(power_w)?;
        if initial_k.len() != self.nodes.len() {
            return Err(lumped_error(
                "lumped network",
                format!(
                    "initial vector has {} entries for {} nodes",
                    initial_k.len(),
                    self.nodes.len()
                ),
            ));
        }
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(lumped_error(
                "lumped network",
                format!("time {time_s} s is not finite and non-negative"),
            ));
        }
        let steady = self.steady(power_w)?;
        let mut out = Vec::with_capacity(self.nodes.len());
        for ((node, start), settle) in self.nodes.iter().zip(initial_k.iter()).zip(steady.iter()) {
            let decay = fs_math::det::exp(-time_s / node.time_constant_s());
            let value = settle + (start - settle) * decay;
            if !value.is_finite() {
                return Err(lumped_error(
                    &node.name,
                    format!("response at {time_s} s left the finite range"),
                ));
            }
            out.push(value);
        }
        Ok(out)
    }

    fn admit_power(&self, power_w: &[f64]) -> Result<(), ConductionError> {
        if power_w.len() != self.nodes.len() {
            return Err(lumped_error(
                "lumped network",
                format!(
                    "power vector has {} entries for {} nodes",
                    power_w.len(),
                    self.nodes.len()
                ),
            ));
        }
        for (node, power) in self.nodes.iter().zip(power_w.iter()) {
            if !power.is_finite() || *power < 0.0 {
                return Err(lumped_error(
                    &node.name,
                    format!("node power {power} W is not finite and non-negative"),
                ));
            }
        }
        Ok(())
    }
}

/// A solved reduced-rung response, with the verdict that admitted it.
#[derive(Debug, Clone, PartialEq)]
pub struct LumpedSolution {
    /// Nodal temperatures at the requested time, K.
    pub temperature_k: Vec<f64>,
    /// Steady temperatures the response is approaching, K.
    pub steady_k: Vec<f64>,
    /// The validity verdict under which this was produced.
    pub verdict: ValidityVerdict,
}

/// Solve the reduced rung, GATED.
///
/// The gate runs first and a refusal short-circuits: outside the declared
/// Biot regime this returns an error rather than a number, because the whole
/// value of a cheap rung is destroyed if it answers questions it cannot
/// answer.
///
/// # Errors
/// [`ConductionError::ScenarioRow`] when the gate refuses or any input is
/// inadmissible.
pub fn solve_gated(
    network: &LumpedNetwork,
    gate: BiotGate,
    power_w: &[f64],
    initial_k: &[f64],
    time_s: f64,
) -> Result<LumpedSolution, ConductionError> {
    let verdict = gate.adjudicate(network);
    match &verdict {
        ValidityVerdict::Refused {
            node,
            biot,
            ceiling,
        } => Err(lumped_error(
            node,
            format!(
                "Biot {biot} exceeds the admitted ceiling {ceiling}: the body is not isothermal, so the lumped rung cannot answer this and escalating to the full transient is the correct move rather than loosening the gate"
            ),
        )),
        ValidityVerdict::Admitted { .. } => Ok(LumpedSolution {
            temperature_k: network.response_at(power_w, initial_k, time_s)?,
            steady_k: network.steady(power_w)?,
            verdict,
        }),
    }
}

/// Build a single-node reduced model from a measured steady response.
///
/// Given the steady temperature rise a body reaches under a known power, the
/// surface conductance is `hA = P / dT`. That is an EXTRACTION, not a
/// derivation: it assumes the steady rise is dominated by the surface path
/// the reduced model represents, and it inherits every modelling assumption
/// of the run it was extracted from. Callers should record it as model-form
/// evidence, not as a measured property.
///
/// # Errors
/// [`ConductionError::ScenarioRow`] for a non-positive rise or power, or any
/// non-finite input.
#[allow(clippy::too_many_arguments)]
pub fn extract_node_from_steady_rise(
    name: impl Into<String>,
    power_w: f64,
    steady_rise_k: f64,
    capacitance_j_per_k: f64,
    characteristic_length_m: f64,
    conductivity_w_per_m_k: f64,
    surface_area_m2: f64,
) -> Result<LumpedNode, ConductionError> {
    let name = name.into();
    for (value, field) in [
        (power_w, "extraction power"),
        (steady_rise_k, "steady rise"),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(lumped_error(
                &name,
                format!("{field} {value} is not finite and positive"),
            ));
        }
    }
    let conductance = power_w / steady_rise_k;
    LumpedNode::new(
        name,
        capacitance_j_per_k,
        conductance,
        characteristic_length_m,
        conductivity_w_per_m_k,
        surface_area_m2,
    )
}

fn lumped_error(node: &str, what: String) -> ConductionError {
    ConductionError::ScenarioRow {
        region: node.to_string(),
        what,
        fix: "correct the declared lumped network, or escalate to the full transient rung"
            .to_string(),
    }
}
