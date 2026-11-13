//! Finite material-body heat feeding simultaneous moving-gas molecular flow.
//!
//! All body groups first exchange heat with their own ORIGINAL finite gas,
//! using the existing simultaneous implicit multi-body owner. The resulting
//! gas states then feed one original molecular network transaction. Both
//! stages use the SAME physical duration; this is first-order Lie splitting.
//!
//! Thermal contacts change internal storage, not gas mass, volume or momentum.
//! Their mechanical pressure/drag tractions remain separate. Fixed-volume,
//! constant-density material charts give wall delta-h its internal-energy
//! meaning here; expansion requires its own work owner. A convection law must
//! qualify the actual rarefied-flow/thermal contact: no coefficient is inferred
//! from the molecular aperture. All gas and material endpoints publish together.

/// Adaptive accepted thermal/molecular histories with fixed physical budgets.
pub mod trajectory;

use super::{
    MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS, MolecularExchangeError, MolecularNetworkControls,
    MolecularNetworkError, MolecularNetworkLink, MolecularNetworkStep, MovingGasBoundaryError,
    MovingGasBoundarySolve, MovingGasParcel,
};
use crate::{
    cavity_phs::piston::thermal::HeatExchangeSolve,
    finite_gas_thermal::{
        ConvectionCoefficient, ConvectionError,
        bodies::{
            FiniteGasBodiesCandidate, FiniteGasBodiesEnvironment, FiniteGasBodiesError,
            FiniteGasBody, MAX_BODIES,
        },
    },
    liquid_gas::evaporation::{EvaporatingLiquidGasError, momentum::boundary::kinetic_change},
};
use fs_exec::Cx;
use fs_material::{
    phase::EquilibriumPhaseState,
    thermochemical_gas::{
        ResolvedThermochemicalGasState,
        energy::{IsochoricGasParcel, parcel::ClosedGasParcel},
    },
};
use fs_math::dd::Dd;
use std::collections::{BTreeMap, BTreeSet};

/// The native per-gas and per-body caps also bound this complete transaction.
pub const MAX_NETWORK_THERMAL_BODIES: usize = MAX_SNAPSHOT_PARCELS * MAX_BODIES;

/// Distinct physical bodies sharing ONE original gas endpoint during heat flow.
/// Gas IDs are unique across groups, body IDs across the whole network. Bodies
/// may share immutable material/environment objects, but not physical inventory.
#[derive(Clone, Copy, Debug)]
pub struct MolecularGasThermalBoundary<'a> {
    pub gas: u64,
    pub bodies: &'a [FiniteGasBody<'a>],
    /// Authored gas heat [J], excluding all reciprocal body transfers.
    pub gas_heat_j: f64,
    pub exchange: HeatExchangeSolve,
    pub maximum_gas_temperature_increment_k: f64,
}

/// Native molecular controls plus one fixed budget for the complete proposal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularNetworkThermalControls {
    pub exchange: MolecularNetworkControls,
    /// Total distinct body count, at most MAX_NETWORK_THERMAL_BODIES.
    /// Zero permits only an empty thermal boundary list.
    pub maximum_bodies: usize,
    /// Noncancelling thermal/molecular defects in both energy frames [J].
    pub energy_tolerance_j: f64,
    /// Noncancelling actual/requested momentum defects, per component [N s].
    pub momentum_tolerance_n_s: f64,
}

/// Original owner failures remain inspectable and never publish partial state.
#[derive(Debug)]
pub enum MolecularNetworkThermalError<E = core::convert::Infallible> {
    Invalid(&'static str),
    Thermal {
        gas: u64,
        error: FiniteGasBodiesError<E>,
    },
    Motion {
        gas: u64,
        error: EvaporatingLiquidGasError,
    },
    Molecular(MolecularNetworkError),
    Accounting(MolecularExchangeError),
    EnergyBalance {
        residual_j: f64,
        laboratory_residual_j: f64,
        absolute_defect_j: f64,
    },
    MomentumBalance {
        residual_n_s: [f64; 3],
        absolute_defect_n_s: [f64; 3],
    },
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for MolecularNetworkThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-body molecular network refused: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for MolecularNetworkThermalError<E> {}
type Error<E> = MolecularNetworkThermalError<E>;

#[derive(Clone, Debug, PartialEq)]
pub struct MolecularNetworkThermalStep<S = ()> {
    pub duration_s: f64,
    /// Canonical gas order. Native candidates retain the heated gas and all
    /// actual body endpoints, constitutive evaluations and reciprocal accounts.
    pub thermal: BTreeMap<u64, FiniteGasBodiesCandidate<S, IsochoricGasParcel>>,
    pub molecular: MolecularNetworkStep,
    /// Final combined-endpoint source admission, keyed by physical body ID.
    /// These receipts apply no additional heat; heat-stage samples remain above.
    pub final_convection: BTreeMap<u64, ConvectionCoefficient<S>>,
    /// Gas/body authored sources plus enclosing radiation, never internal Q.
    pub external_heat_j: f64,
    /// Original gas caloric integrals, at unchanged composition and mass.
    pub gas_heat_storage_change_j: f64,
    /// Actual returned finite material storage, counted once per body [J].
    pub wall_internal_energy_change_j: f64,
    pub thermal_momentum_change_n_s: [f64; 3],
    pub thermal_kinetic_change_j: f64,
    pub thermal_laboratory_kinetic_change_j: f64,
    /// Native heat defects plus actual storage/kinetic reconstruction corrections.
    /// Retained separately so adaptive histories never subtract rounded totals.
    pub thermal_absolute_energy_defect_j: f64,
    pub thermal_absolute_momentum_defect_n_s: [f64; 3],
    /// Reconstructed thermal defect plus native composition-aware molecular
    /// closure. Original species formation references remain canceled.
    pub combined_energy_residual_j: f64,
    pub combined_laboratory_energy_residual_j: f64,
    pub absolute_energy_defect_j: f64,
    /// Actual original-to-final gas momentum in the native inertial axes.
    pub momentum_residual_n_s: [f64; 3],
    pub absolute_momentum_defect_n_s: [f64; 3],
}

#[derive(Clone, Debug, PartialEq)]
pub struct MolecularNetworkThermalCandidate<S = ()> {
    gases: BTreeMap<u64, MovingGasParcel>,
    bodies: BTreeMap<u64, EquilibriumPhaseState>,
    molecular_initial_gases: BTreeMap<u64, MovingGasParcel>,
    step: MolecularNetworkThermalStep<S>,
}
impl<S> MolecularNetworkThermalCandidate<S> {
    pub fn gases(&self) -> &BTreeMap<u64, MovingGasParcel> {
        &self.gases
    }
    pub fn bodies(&self) -> &BTreeMap<u64, EquilibriumPhaseState> {
        &self.bodies
    }
    /// Exact heated snapshot sampled by the native molecular stage. Adaptive
    /// consumers use its temperatures and drifts, without reconstructing heat.
    pub fn molecular_initial_gases(&self) -> &BTreeMap<u64, MovingGasParcel> {
        &self.molecular_initial_gases
    }
    pub const fn step(&self) -> &MolecularNetworkThermalStep<S> {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, MovingGasParcel>,
        BTreeMap<u64, EquilibriumPhaseState>,
        MolecularNetworkThermalStep<S>,
    ) {
        (self.gases, self.bodies, self.step)
    }
}

fn poll<E>(cx: &Cx<'_>) -> Result<(), Error<E>> {
    cx.checkpoint().map_err(|_| Error::Cancelled)
}
fn network_error<E>(error: MolecularNetworkError) -> Error<E> {
    match error {
        MolecularNetworkError::Cancelled => Error::Cancelled,
        error => Error::Molecular(error),
    }
}
fn number<E>(x: Dd) -> Result<f64, Error<E>> {
    let value = super::finite(x).map_err(Error::Accounting)?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::Invalid(
            "complete thermal/molecular account is unrepresentable",
        ))
    }
}
fn vector<E>(x: [Dd; 3]) -> Result<[f64; 3], Error<E>> {
    Ok([number(x[0])?, number(x[1])?, number(x[2])?])
}
fn maximum(a: Dd, b: Dd) -> Dd {
    if a.lt(b) { b } else { a }
}

fn admit<E>(
    cx: &Cx<'_>,
    gases: &BTreeMap<u64, MovingGasParcel>,
    links: &[MolecularNetworkLink<'_>],
    duration: f64,
    endpoints: &BTreeMap<u64, MovingGasBoundarySolve>,
    boundaries: &[MolecularGasThermalBoundary<'_>],
    controls: MolecularNetworkThermalControls,
) -> Result<(), Error<E>> {
    poll(cx)?;
    if gases.is_empty()
        || gases.len() > MAX_SNAPSHOT_PARCELS
        || links.len() > MAX_SNAPSHOT_TRANSFERS
        || boundaries.len() > gases.len()
        || !gases.keys().eq(endpoints.keys())
        || controls.maximum_bodies > MAX_NETWORK_THERMAL_BODIES
        || ![
            duration,
            controls.exchange.maximum_duration_s,
            controls.energy_tolerance_j,
            controls.momentum_tolerance_n_s,
        ]
        .into_iter()
        .all(|v| v.is_finite() && v > 0.0)
    {
        return Err(Error::Invalid(
            "bounded original inventories, exact solve IDs, positive clocks and fixed physical budgets required",
        ));
    }
    if duration > controls.exchange.maximum_duration_s {
        return Err(Error::Molecular(MolecularNetworkError::Accounting(
            MolecularExchangeError::Duration {
                duration_s: duration,
                maximum_duration_s: controls.exchange.maximum_duration_s,
            },
        )));
    }
    let mut group_ids = BTreeSet::new();
    let mut body_ids = BTreeSet::new();
    for group in boundaries {
        poll(cx)?;
        if !gases.contains_key(&group.gas)
            || !group_ids.insert(group.gas)
            || group.bodies.is_empty()
            || group.bodies.len() > MAX_BODIES
            || !group.gas_heat_j.is_finite()
            || !group.maximum_gas_temperature_increment_k.is_finite()
            || group.maximum_gas_temperature_increment_k <= 0.0
        {
            return Err(Error::Invalid(
                "unique existing thermal gas IDs, nonempty bounded bodies and finite heat/increment required",
            ));
        }
        for body in group.bodies {
            poll(cx)?;
            if !body_ids.insert(body.id) || body_ids.len() > controls.maximum_bodies {
                return Err(Error::Invalid(
                    "physical thermal body IDs must be globally unique and fit the shared budget",
                ));
            }
            let input = body.input;
            let density = input.initial.bulk_density_kg_m3();
            let volume = input.mass_kg / density;
            if input.duration_s != duration
                || !input.mass_kg.is_finite()
                || input.mass_kg <= 0.0
                || !volume.is_finite()
                || volume <= 0.0
                || !input.volume_m3.is_finite()
                || input.volume_m3 <= 0.0
                || (input.volume_m3 - volume).abs()
                    > 64.0 * f64::EPSILON * input.volume_m3.max(volume)
            {
                return Err(Error::Invalid(
                    "thermal bodies require one physical clock and fixed mass/density/volume",
                ));
            }
            for knot in input.curve.knots() {
                poll(cx)?;
                if knot.bulk_density_kg_m3 != density {
                    return Err(Error::Invalid(
                        "body expansion needs its work owner; this network requires constant-density calorics",
                    ));
                }
            }
        }
    }
    // Heating must not conceal an initially unsupported active source or drift.
    // Native exchange still owns normal, reciprocal and endpoint admission.
    let mut ordered = BTreeMap::new();
    for link in links {
        poll(cx)?;
        if link.gases[0] == link.gases[1]
            || link.gases.iter().any(|id| !gases.contains_key(id))
            || !link.surface.area_m2.is_finite()
            || link.surface.area_m2 < 0.0
            || ordered.insert(link.id, link).is_some()
        {
            return Err(Error::Invalid(
                "unique aperture IDs, existing distinct gases and nonnegative areas required",
            ));
        }
    }
    for link in ordered.values().filter(|link| link.surface.area_m2 > 0.0) {
        for id in link.gases {
            poll(cx)?;
            let gas = &gases[&id];
            link.interface
                .admit_relative_velocity(core::array::from_fn(|d| {
                    gas.velocity_m_s()[d] - link.surface.velocity_world_m_s[d]
                }))
                .map_err(|error| network_error(super::source(link.id, error)))?;
            link.interface
                .momentum_interface()
                .interface()
                .admit_parcel(gas.parcel(), || cx.checkpoint().is_err())
                .map_err(|error| network_error(super::source(link.id, error)))?;
        }
    }
    Ok(())
}

struct HeatAccount {
    gas: Dd,
    wall: Dd,
    external: Dd,
    momentum: [Dd; 3],
    momentum_debit: [Dd; 3],
    kinetic: Dd,
    laboratory_kinetic: Dd,
    residual: Dd,
    laboratory_residual: Dd,
    debit: Dd,
}
impl Default for HeatAccount {
    fn default() -> Self {
        Self {
            gas: Dd::ZERO,
            wall: Dd::ZERO,
            external: Dd::ZERO,
            momentum: [Dd::ZERO; 3],
            momentum_debit: [Dd::ZERO; 3],
            kinetic: Dd::ZERO,
            laboratory_kinetic: Dd::ZERO,
            residual: Dd::ZERO,
            laboratory_residual: Dd::ZERO,
            debit: Dd::ZERO,
        }
    }
}
impl HeatAccount {
    #[allow(clippy::too_many_arguments)]
    fn record<S, E>(
        &mut self,
        cx: &Cx<'_>,
        original: &MovingGasParcel,
        heated: &MovingGasParcel,
        group: MolecularGasThermalBoundary<'_>,
        candidate: &FiniteGasBodiesCandidate<S, IsochoricGasParcel>,
        frame: [f64; 3],
        bodies: &mut BTreeMap<u64, EquilibriumPhaseState>,
    ) -> Result<(), Error<E>> {
        poll(cx)?;
        let gas_error = |error| {
            network_error(super::gas_error(
                group.gas,
                MovingGasBoundaryError::Gas(error),
            ))
        };
        // The original native state-difference owner admits identical mass,
        // composition, sources and auxiliary coordinates before integrating Cv.
        let gas = Dd::from_f64(
            original
                .parcel()
                .stored_energy_change_to(heated.parcel(), || cx.checkpoint().is_err())
                .map_err(gas_error)?,
        );
        let kinetic = Dd::from_f64(
            kinetic_change(
                original.parcel().mass_kg(),
                heated.parcel().mass_kg(),
                original.momentum_kg_m_s(),
                heated.momentum_kg_m_s(),
                frame,
            )
            .map_err(gas_error)?,
        );
        let laboratory = Dd::from_f64(
            kinetic_change(
                original.parcel().mass_kg(),
                heated.parcel().mass_kg(),
                original.momentum_kg_m_s(),
                heated.momentum_kg_m_s(),
                [0.0; 3],
            )
            .map_err(gas_error)?,
        );
        let inputs: BTreeMap<_, _> = group
            .bodies
            .iter()
            .map(|body| (body.id, body.input))
            .collect();
        let (mut wall, mut external) = (Dd::ZERO, Dd::from_f64(group.gas_heat_j));
        if candidate.bodies().len() != inputs.len() {
            return Err(Error::Invalid(
                "native thermal owner changed physical body inventory",
            ));
        }
        for body in candidate.bodies() {
            poll(cx)?;
            let input = inputs.get(&body.report.id).ok_or(Error::Invalid(
                "native thermal owner changed physical body identity",
            ))?;
            if bodies.insert(body.report.id, body.state).is_some() {
                return Err(Error::Invalid("native thermal body was published twice"));
            }
            wall = wall
                + Dd::from_f64(input.mass_kg)
                    * (Dd::from_f64(body.state.specific_enthalpy_j_kg())
                        - Dd::from_f64(input.initial.specific_enthalpy_j_kg()));
            external = external
                + Dd::from_f64(input.internal_heat_j)
                + Dd::from_f64(body.report.radiation_heat_j);
        }
        let residual = gas + wall + kinetic - external;
        let laboratory_residual = gas + wall + laboratory - external;
        let native = Dd::from_f64(candidate.report().combined_energy_residual_j);
        // A disagreement with the native thermal account is an explicit debit,
        // including actual bulk motion; it never becomes invented physical heat.
        let correction = maximum(
            (residual - native).abs(),
            (laboratory_residual - native).abs(),
        );
        self.debit =
            self.debit + Dd::from_f64(candidate.report().absolute_energy_defect_j) + correction;
        self.residual = self.residual + residual;
        self.laboratory_residual = self.laboratory_residual + laboratory_residual;
        self.gas = self.gas + gas;
        self.wall = self.wall + wall;
        self.external = self.external + external;
        self.kinetic = self.kinetic + kinetic;
        self.laboratory_kinetic = self.laboratory_kinetic + laboratory;
        for d in 0..3 {
            let dp = Dd::from_f64(heated.momentum_kg_m_s()[d])
                - Dd::from_f64(original.momentum_kg_m_s()[d]);
            self.momentum[d] = self.momentum[d] + dp;
            self.momentum_debit[d] = self.momentum_debit[d] + dp.abs();
        }
        Ok(())
    }
}

/// Use each original body's declared whole-surface heat-transfer coefficient.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate(
    cx: &Cx<'_>,
    gases: &BTreeMap<u64, MovingGasParcel>,
    links: &[MolecularNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, MovingGasBoundarySolve>,
    boundaries: &[MolecularGasThermalBoundary<'_>],
    controls: MolecularNetworkThermalControls,
) -> Result<MolecularNetworkThermalCandidate, MolecularNetworkThermalError> {
    step_candidate_with_convection(
        cx,
        gases,
        links,
        duration_s,
        endpoints,
        boundaries,
        controls,
        |_, body, _, _| {
            Ok::<_, core::convert::Infallible>(ConvectionCoefficient {
                coefficient_w_per_m2_k: body.environment.convection_w_per_m2_k,
                source: (),
            })
        },
    )
}

/// Source-qualified state-dependent h at the original implicit gas endpoint.
/// During heating the callback receives gas ID, original body descriptor,
/// actual trial gas state and current gas velocity in the original inertial
/// axes. Final source admission replaces only the descriptor's initial material
/// state with the accepted wall endpoint, retaining its curve, geometry, sources
/// and duration. The callback owns the wall/relative-flow convention and must
/// be repeatable; callback side effects are not rolled back. The same law must
/// also admit the actual combined gas
/// endpoint and post-exchange velocity; those final receipts apply no extra heat.
/// No viscosity, conductivity or rarefied-flow law is copied.
///
/// All thermal groups are independent proposals from the same original gas
/// snapshot. Multiple bodies sharing one gas are solved simultaneously by the
/// existing thermal owner, not heated sequentially. One subsequent molecular
/// transaction sees ALL heated states. Empty thermal boundaries reproduce the
/// native network, subject to the additional fixed combined defect budgets.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate_with_convection<S, E>(
    cx: &Cx<'_>,
    gases: &BTreeMap<u64, MovingGasParcel>,
    links: &[MolecularNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, MovingGasBoundarySolve>,
    boundaries: &[MolecularGasThermalBoundary<'_>],
    controls: MolecularNetworkThermalControls,
    mut resolve: impl FnMut(
        u64,
        &FiniteGasBody<'_>,
        &ResolvedThermochemicalGasState,
        [f64; 3],
    ) -> Result<ConvectionCoefficient<S>, E>,
) -> Result<MolecularNetworkThermalCandidate<S>, Error<E>> {
    admit(
        cx, gases, links, duration_s, endpoints, boundaries, controls,
    )?;
    let ordered: BTreeMap<_, _> = boundaries.iter().map(|group| (group.gas, *group)).collect();
    let mut heated = gases.clone();
    let mut thermal = BTreeMap::new();
    let mut bodies = BTreeMap::new();
    let mut heat = HeatAccount::default();
    for (&id, group) in &ordered {
        poll(cx)?;
        let original = &gases[&id];
        let parcel = original.parcel().clone().into_isochoric();
        let candidate = FiniteGasBodiesEnvironment {
            gas: &parcel,
            gas_heat_j: group.gas_heat_j,
            gas_solve: endpoints[&id].gas.gas,
            exchange: group.exchange,
            maximum_gas_temperature_increment_k: group.maximum_gas_temperature_increment_k,
            maximum_bodies: MAX_BODIES,
        }
        .advance_with_convection(cx, group.bodies, |body, state| {
            resolve(id, body, state, original.velocity_m_s())
        })
        .map_err(|error| match error {
            FiniteGasBodiesError::Cancelled => Error::Cancelled,
            error => Error::Thermal { gas: id, error },
        })?;
        let next = ClosedGasParcel::from_isochoric(candidate.gas().clone());
        if next.mass_kg().to_bits() != original.parcel().mass_kg().to_bits()
            || next.volume_m3().to_bits() != original.parcel().volume_m3().to_bits()
        {
            return Err(Error::Invalid(
                "thermal owner changed original fixed gas inventory or volume",
            ));
        }
        let next = MovingGasParcel::from_momentum(next, original.momentum_kg_m_s())
            .map_err(|error| Error::Motion { gas: id, error })?;
        heat.record(
            cx,
            original,
            &next,
            *group,
            &candidate,
            controls.exchange.energy_frame_velocity_m_s,
            &mut bodies,
        )?;
        heated.insert(id, next);
        thermal.insert(id, candidate);
    }
    poll(cx)?;
    let molecular = super::step_candidate(
        &heated,
        links,
        duration_s,
        endpoints,
        controls.exchange,
        || cx.checkpoint().is_err(),
    )
    .map_err(network_error)?;
    // Molecular transfer changes temperature, composition and mean velocity.
    // Qualify that ACTUAL endpoint and accepted wall state too, without another
    // heat or body solve. Only the body descriptor's initial state is replaced.
    let mut final_convection = BTreeMap::new();
    for (&id, group) in &ordered {
        let gas = &molecular.gases()[&id];
        let mut ordered_bodies: Vec<_> = group.bodies.iter().collect();
        ordered_bodies.sort_unstable_by_key(|body| body.id);
        for body in ordered_bodies {
            poll(cx)?;
            let mut endpoint_body = *body;
            endpoint_body.input.initial = bodies[&body.id];
            let sample = resolve(id, &endpoint_body, gas.parcel().state(), gas.velocity_m_s())
                .map_err(|error| Error::Thermal {
                    gas: id,
                    error: FiniteGasBodiesError::Body {
                        id: body.id,
                        error: ConvectionError::Convection(error),
                    },
                })?;
            poll(cx)?;
            if !sample.coefficient_w_per_m2_k.is_finite() || sample.coefficient_w_per_m2_k < 0.0 {
                return Err(Error::Invalid(
                    "final convection coefficient must be finite and nonnegative",
                ));
            }
            final_convection.insert(body.id, sample);
        }
    }
    let step = molecular.step();
    let mut molecular_debit = Dd::from_f64(step.pooled_stream_energy_residual_j).abs();
    let mut requested = [Dd::ZERO; 3];
    let mut momentum_debit = heat.momentum_debit;
    for gas in step.gas.values() {
        poll(cx)?;
        molecular_debit = molecular_debit + Dd::from_f64(gas.total_energy_residual_j).abs();
        for d in 0..3 {
            requested[d] = requested[d] + Dd::from_f64(gas.requested_impulse_n_s[d]);
            momentum_debit[d] =
                momentum_debit[d] + Dd::from_f64(gas.momentum_residual_n_s[d]).abs();
        }
    }
    let mut momentum = [Dd::ZERO; 3];
    for (&id, original) in gases {
        poll(cx)?;
        for d in 0..3 {
            momentum[d] = momentum[d] + Dd::from_f64(molecular.gases()[&id].momentum_kg_m_s()[d])
                - Dd::from_f64(original.momentum_kg_m_s()[d]);
        }
    }
    let mut conversion = Dd::ZERO;
    for d in 0..3 {
        let represented = heat.momentum[d] + Dd::from_f64(step.momentum_residual_n_s[d]);
        momentum_debit[d] =
            momentum_debit[d] + requested[d].abs() + (momentum[d] - represented).abs();
        let velocity = Dd::from_f64(controls.exchange.energy_frame_velocity_m_s[d]);
        conversion = conversion + velocity * Dd::from_f64(step.momentum_residual_n_s[d])
            - (Dd::from_f64(step.mass_residual_kg) * Dd::from_f64(0.5) * velocity) * velocity;
    }
    let momentum_residual_n_s = vector(momentum)?;
    let absolute_momentum_defect_n_s = vector(momentum_debit)?;
    if momentum_residual_n_s
        .iter()
        .any(|v| v.abs() > controls.momentum_tolerance_n_s)
        || absolute_momentum_defect_n_s
            .iter()
            .any(|v| *v > controls.momentum_tolerance_n_s)
    {
        return Err(Error::MomentumBalance {
            residual_n_s: momentum_residual_n_s,
            absolute_defect_n_s: absolute_momentum_defect_n_s,
        });
    }
    let combined_energy_residual_j = number(heat.residual + Dd::from_f64(step.energy_residual_j))?;
    let combined_laboratory_energy_residual_j =
        number(heat.laboratory_residual + Dd::from_f64(step.laboratory_energy_residual_j))?;
    let absolute_energy_defect_j = number(heat.debit + molecular_debit + conversion.abs())?;
    if combined_energy_residual_j.abs() > controls.energy_tolerance_j
        || combined_laboratory_energy_residual_j.abs() > controls.energy_tolerance_j
        || absolute_energy_defect_j > controls.energy_tolerance_j
    {
        return Err(Error::EnergyBalance {
            residual_j: combined_energy_residual_j,
            laboratory_residual_j: combined_laboratory_energy_residual_j,
            absolute_defect_j: absolute_energy_defect_j,
        });
    }
    let (gases, molecular) = molecular.into_parts();
    let step = MolecularNetworkThermalStep {
        duration_s,
        thermal,
        molecular,
        final_convection,
        external_heat_j: number(heat.external)?,
        gas_heat_storage_change_j: number(heat.gas)?,
        wall_internal_energy_change_j: number(heat.wall)?,
        thermal_momentum_change_n_s: vector(heat.momentum)?,
        thermal_kinetic_change_j: number(heat.kinetic)?,
        thermal_laboratory_kinetic_change_j: number(heat.laboratory_kinetic)?,
        thermal_absolute_energy_defect_j: number(heat.debit)?,
        thermal_absolute_momentum_defect_n_s: vector(heat.momentum_debit)?,
        combined_energy_residual_j,
        combined_laboratory_energy_residual_j,
        absolute_energy_defect_j,
        momentum_residual_n_s,
        absolute_momentum_defect_n_s,
    };
    poll(cx)?;
    Ok(MolecularNetworkThermalCandidate {
        gases,
        bodies,
        molecular_initial_gases: heated,
        step,
    })
}
