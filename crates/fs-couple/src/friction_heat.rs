//! General friction work coupled to the existing finite enthalpy network.
//!
//! An existing contact patch becomes a zero-capacity, common-temperature
//! interface. Its two conductances use the ORIGINAL cell geometry, endpoint
//! material conductivities and explicitly supplied contact resistances. The
//! original direct link is suppressed only in the trial, retaining its surface
//! coverage; otherwise conduction would
//! be counted twice. Signed body heat includes ordinary hot-to-cold exchange.
//! Only the sum equals exported frictional work. No heat fraction or store is
//! introduced, and interface temperature is a solved coordinate, not an average.
//!
//! This is a cell-centred, orthogonal finite-volume discretization: cell
//! centrelines run between opposite face centres of uniform prismatic cells.
//! The contact is their shared face. The default has zero interfacial resistance;
//! a finite resistance split locates the massless heat-generation plane.
//! Spatial and temporal refinement remain necessary. Moving footprints, thermal
//! expansion, finite third-body storage and a flash-temperature approximation
//! are not supplied by this module. The source owns any such constitutive work.

use fs_math::det;
use fs_thermal::enthalpy::backward_euler::{
    BackwardEulerStep, BackwardEulerTrial, StepCoefficients,
};
use fs_thermal::enthalpy::{
    CellGeometry, CellId, ContactPatch, DepositionError, DepositionSimulation, EnergyAudit,
    MaterialId, series_contact_conductance_w_k,
};
use fs_tribo::{TangentialLaw, TangentialStep, TangentialStepInput, TriboError};
use std::cell::RefCell;

mod entropy;
pub use entropy::FrictionHeatEntropy;

mod area;
pub use area::FrictionHeatAreaTrial;

/// Physical finite-layer storage and source-derived conduction into both bodies.
pub mod third_body;

/// Conservative moving contact patches with locally resolved temperatures.
pub mod moving;

/// Bounded pulse-history flash-temperature comparison against resolved conduction.
pub mod flash;

/// Native rigid momenta, retained partial slip and finite heat accept together.
pub mod rigid;

/// Original directional material friction, rigid recoil and resolved heat partition.
pub mod directional;

/// Simultaneous friction sources on one original thermal network.
pub mod network;

/// Prepared native articulated contacts and original finite-cell heat together.
pub mod articulated;

/// Original adhesive junctions, simultaneous native recoil and finite-cell heat.
pub mod adhesive;
/// Native deformable cohesive reclosure, physical interface loss and finite heat.
pub mod cohesive;
/// Native SDF solid contact, source-temperature feedback and finite-cell heat.
#[cfg(feature = "experimental-solid-contact-heat")]
pub mod solid;
/// Native sustained finite-solid contact, original source histories and finite heat.
pub mod finite_solid;
/// Native rigid adhesive joints and their original finite thermal receivers.
pub mod rigid_cohesive;
/// Source-backed Hertz sliding/spin recoil and finite-cell temperature feedback.
pub mod hertz;
/// Event-resolved spherical spin braking with original finite-body heat feedback.
pub mod sphere_plane;
pub mod hertz_rolling;
pub mod rolling;

/// Continuous native joint friction, finite-cell heat and endpoint source feedback.
#[path = "friction_heat/joints.rs"]
pub mod joints;

/// Retained source Mindlin histories, shared native bodies and finite heat.
pub mod mindlin;

/// Original hardness-limited sliding and wear with finite-body thermal feedback.
pub mod plastic_wear;

fn invalid(message: &'static str) -> DepositionError {
    DepositionError::InvalidInput(message)
}
fn finite(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic)
    }
}

/// A fixed physical contact in the original network, with no added heat store.
#[derive(Clone, Debug)]
pub struct FrictionHeatRoute {
    patch: ContactPatch,
    geometry: [CellGeometry; 2],
    material: [MaterialId; 2],
    mass_kg: [f64; 2],
    distances_m: [f64; 2],
    area_resistance_m2_k_w: [f64; 2],
}
impl FrictionHeatRoute {
    /// Bind an existing perfect-contact link between adjoining uniform cells.
    ///
    /// Each centreline endpoint is a face centre. The common endpoint is the
    /// interface; each half distance comes from THAT body's own cell length.
    /// Conductivities must be resolved from the two original materials at the
    /// current thermal state. The original link must have their series
    /// conductance, so an authored contact resistance is never silently erased.
    /// Nonmatching reference-area patches and nonorthogonal cells refuse.
    pub fn from_link(
        thermal: &DepositionSimulation<'_>,
        link_id: u64,
        reference_k_w_m_k: [f64; 2],
    ) -> Result<Self, DepositionError> {
        Self::from_link_with_resistance(thermal, link_id, reference_k_w_m_k, [0.0; 2])
    }

    /// Bind finite area-specific resistances [m² K/W] on either side of the
    /// massless friction heat-generation plane. Their sum must be the original
    /// direct link's contact resistance. The split is an explicit physical
    /// assumption about that plane, not an arbitrary prescribed heat fraction.
    /// Source identity and pressure/temperature evolution belong to the caller.
    pub fn from_link_with_resistance(
        thermal: &DepositionSimulation<'_>,
        link_id: u64,
        reference_k_w_m_k: [f64; 2],
        area_resistance_m2_k_w: [f64; 2],
    ) -> Result<Self, DepositionError> {
        let resistance = Self::total_resistance(area_resistance_m2_k_w)?;
        let patch = thermal
            .contact_patch(link_id)
            .ok_or(invalid("inactive friction heat contact"))?;
        if patch.covered_area_m2 != [patch.link.area_m2; 2] {
            return Err(invalid("friction heat needs matching physical face areas"));
        }
        let cells = [
            &thermal.cells()[&patch.link.a],
            &thermal.cells()[&patch.link.b],
        ];
        let geometry = [cells[0].geometry, cells[1].geometry];
        let mut shared = None;
        for a in [geometry[0].start_m, geometry[0].end_m] {
            for b in [geometry[1].start_m, geometry[1].end_m] {
                if a == b {
                    if shared.replace(a).is_some() {
                        return Err(invalid(
                            "friction heat cells overlap instead of sharing one face",
                        ));
                    }
                }
            }
        }
        let interface = shared.ok_or(invalid(
            "friction heat cells need a shared physical face centre",
        ))?;
        let mut distances_m = [0.0; 2];
        let mut normals = [[0.0; 3]; 2];
        for side in 0..2 {
            let g = geometry[side];
            let delta = std::array::from_fn::<_, 3, _>(|i| g.end_m[i] - g.start_m[i]);
            let length = finite(det::hypot(det::hypot(delta[0], delta[1]), delta[2]))?;
            if length <= 0.0
                || {
                    // Subtracting decimal endpoints and reconstructing A=V/L
                    // can differ from the original physical face by roundoff.
                    // This admits arithmetic noise only, not geometric overlap.
                    let face_area = finite(g.volume_m3 / length)?;
                    patch.link.area_m2 - face_area
                        > 64.0 * f64::EPSILON * patch.link.area_m2.max(face_area)
                }
            {
                return Err(invalid(
                    "friction heat patch exceeds its prismatic cell face",
                ));
            }
            distances_m[side] = 0.5 * length;
            for i in 0..3 {
                normals[side][i] = if interface == g.end_m {
                    delta[i] / length
                } else {
                    -delta[i] / length
                };
            }
        }
        let cosine: f64 = normals[0].iter().zip(normals[1]).map(|(&a, b)| a * b).sum();
        if (cosine + 1.0).abs() > 256.0 * f64::EPSILON {
            return Err(invalid(
                "friction heat requires opposing orthogonal cell normals",
            ));
        }
        let expected = series_contact_conductance_w_k(
            patch.link.area_m2,
            distances_m,
            reference_k_w_m_k,
            resistance,
        )?;
        if (patch.link.conductance_w_k - expected).abs()
            > 64.0 * f64::EPSILON * expected.abs().max(patch.link.conductance_w_k.abs())
        {
            return Err(invalid(
                "friction heat link differs from its original material/geometry conductance",
            ));
        }
        let route = Self {
            patch,
            geometry,
            material: [cells[0].material_id, cells[1].material_id],
            mass_kg: [cells[0].mass_kg, cells[1].mass_kg],
            distances_m,
            area_resistance_m2_k_w,
        };
        route.conductances(reference_k_w_m_k)?;
        Ok(route)
    }
    /// Ordered existing destinations, following the original link's A/B order.
    pub const fn cells(&self) -> [CellId; 2] {
        [self.patch.link.a, self.patch.link.b]
    }
    /// Identity of the original admitted thermal contact link.
    pub const fn link_id(&self) -> u64 {
        self.patch.link.id
    }
    /// Physical area of the original admitted contact face [m²].
    pub const fn contact_area_m2(&self) -> f64 {
        self.patch.link.area_m2
    }
    /// Actual per-side centre-to-interface distances [m].
    pub const fn distances_m(&self) -> [f64; 2] {
        self.distances_m
    }
    /// Accepted contact resistance from the heat-generation plane to each body.
    pub const fn area_resistance_m2_k_w(&self) -> [f64; 2] {
        self.area_resistance_m2_k_w
    }
    fn total_resistance(resistance: [f64; 2]) -> Result<f64, DepositionError> {
        if resistance.iter().any(|r| !r.is_finite() || *r < 0.0) {
            return Err(invalid(
                "contact resistances must be finite and nonnegative",
            ));
        }
        finite(resistance[0] + resistance[1])
    }
    fn validate(&self, thermal: &DepositionSimulation<'_>) -> Result<(), DepositionError> {
        if thermal.contact_patch(self.patch.link.id) != Some(self.patch) {
            return Err(invalid(
                "friction heat route has stale topology or contact coverage",
            ));
        }
        for (i, id) in self.cells().into_iter().enumerate() {
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(invalid("inactive friction heat cell"))?;
            if cell.geometry != self.geometry[i]
                || cell.material_id != self.material[i]
                || cell.mass_kg != self.mass_kg[i]
            {
                return Err(invalid(
                    "friction heat route changed its original cell geometry, material or mass",
                ));
            }
        }
        Ok(())
    }
    fn conductances(&self, conductivity: [f64; 2]) -> Result<[f64; 2], DepositionError> {
        let mut result = [0.0; 2];
        for i in 0..2 {
            if !conductivity[i].is_finite() || conductivity[i] < 0.0 {
                return Err(invalid(
                    "friction heat conductivity must be finite nonnegative material data",
                ));
            }
            result[i] = if self.area_resistance_m2_k_w[i] == 0.0 {
                // Preserve the original perfect-contact arithmetic exactly.
                finite(self.patch.link.area_m2 * (conductivity[i] / self.distances_m[i]))?
            } else if conductivity[i] == 0.0 {
                0.0
            } else {
                let bulk = finite(self.distances_m[i] / conductivity[i])?;
                let total = finite(bulk + self.area_resistance_m2_k_w[i])?;
                if bulk <= 0.0 || total <= 0.0 {
                    return Err(DepositionError::Arithmetic);
                }
                finite(self.patch.link.area_m2 / total)?
            };
            if conductivity[i] > 0.0 && result[i] == 0.0 {
                return Err(DepositionError::Arithmetic);
            }
        }
        if finite(result[0] + result[1])? <= 0.0 {
            return Err(invalid(
                "friction heat interface needs at least one conducting body",
            ));
        }
        Ok(result)
    }
}

/// One endpoint's original mechanical heat and material/thermal source result.
/// The callback must stage physical history, never publish it during a trial.
pub struct FrictionHeatInput<P> {
    /// Exported physical friction loss [J], excluding stored layer energy and
    /// numerical damping. It must not also appear in `transport.sources_w`.
    pub dissipated_j: f64,
    /// Original material conductivities at each body's own endpoint state.
    pub conductivity_w_m_k: [f64; 2],
    /// OTHER heat/conduction and the exact staged mechanical/source payload.
    /// The bound interface link must not be overridden here.
    pub transport: StepCoefficients<P>,
}

/// Independent energy gates, all in joules; none relaxes the other gates.
#[derive(Clone, Copy, Debug)]
pub struct FrictionHeatTolerances {
    pub cell_j: f64,
    pub thermal_j: f64,
    pub interface_j: f64,
    pub total_j: f64,
}

/// Original thermal audit plus independently reconstructed interface accounts.
#[derive(Clone, Debug)]
pub struct FrictionHeatReport {
    pub interface_temperature_k: f64,
    pub conductances_w_k: [f64; 2],
    /// Signed net heat into each body, including hot-to-cold redistribution.
    pub body_heat_j: [f64; 2],
    pub dissipated_j: f64,
    pub interface_residual_j: f64,
    /// Authored/callback/environment heat, measured without friction transfers.
    pub external_heat_j: f64,
    pub thermal: EnergyAudit,
    /// Actual thermal storage minus outside supply minus exported friction work.
    pub total_energy_residual_j: f64,
    /// Optional physical/source/numerical entropy from this same trial. Present
    /// only for `evaluate_with_entropy`; accepted only with the complete trial.
    pub entropy: Option<FrictionHeatEntropy>,
}

/// Coupled residual; the caller's existing nonlinear solver owns iteration.
/// Unknowns are the original cell energy increments and one interface Kelvin
/// coordinate. No interface capacitance or artificial material is introduced.
pub struct FrictionHeatStep<'s, 'm> {
    thermal: BackwardEulerStep<'s, 'm>,
    route: &'s FrictionHeatRoute,
}
impl<'s, 'm> FrictionHeatStep<'s, 'm> {
    pub fn new(
        thermal: &'s DepositionSimulation<'m>,
        target_s: f64,
        route: &'s FrictionHeatRoute,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        route.validate(thermal)?;
        Ok(Self {
            thermal: BackwardEulerStep::new(thermal, target_s, checkpoint)?,
            route,
        })
    }
    pub fn cell_ids(&self) -> &[CellId] {
        self.thermal.cell_ids()
    }
    pub fn duration_s(&self) -> f64 {
        self.thermal.duration_s()
    }

    /// Evaluate mechanics, material transport and heat at the SAME trial.
    /// All original accepted states remain untouched, including on refusal.
    pub fn evaluate<P, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        interface_temperature_k: f64,
        resolve: impl FnOnce(&DepositionSimulation<'m>, f64) -> Result<FrictionHeatInput<P>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<FrictionHeatTrial<'m, P>, E> {
        self.evaluate_entropy_option(energy_increments_j, interface_temperature_k,
            resolve, false, checkpoint)
    }

    /// Evaluate the SAME mechanical/thermal candidate with exact caloric entropy.
    /// The actual interface temperature classifies exported friction work as
    /// D/T_interface and each signed interface-to-body heat as conductive entropy
    /// production. Other sources retain their original unclassified status.
    /// Every active thermal law must provide its exact entropy integral.
    /// Use this on a final candidate to keep entropy work outside root iterations.
    pub fn evaluate_with_entropy<P, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        interface_temperature_k: f64,
        resolve: impl FnOnce(&DepositionSimulation<'m>, f64) -> Result<FrictionHeatInput<P>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<FrictionHeatTrial<'m, P>, E> {
        self.evaluate_entropy_option(energy_increments_j, interface_temperature_k,
            resolve, true, checkpoint)
    }

    fn evaluate_entropy_option<P, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        interface_temperature_k: f64,
        resolve: impl FnOnce(&DepositionSimulation<'m>, f64) -> Result<FrictionHeatInput<P>, E>,
        with_entropy: bool,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<FrictionHeatTrial<'m, P>, E> {
        if !interface_temperature_k.is_finite() || interface_temperature_k <= 0.0 {
            return Err(
                invalid("friction interface temperature must be finite positive kelvin").into(),
            );
        }
        let checkpoint = RefCell::new(checkpoint);
        let check = || (*checkpoint.borrow_mut())();
        let coefficients =
            |state: &DepositionSimulation<'m>| -> Result<StepCoefficients<TrialPayload<P>>, E> {
                let input = resolve(state, interface_temperature_k)?;
                self.coefficients_for_input(self.route, state, interface_temperature_k,
                    input, with_entropy, &check)
            };
        let trial = if with_entropy {
            self.thermal.evaluate_with_entropy(energy_increments_j, coefficients, &check)?
        } else {
            self.thermal.evaluate(energy_increments_j, coefficients, &check)?
        };
        self.finish_trial(trial, interface_temperature_k).map_err(E::from)
    }

    // Fixed and changing-area routes use this same physical transfer, source
    // classification and entropy assembly. Only endpoint geometry differs.
    #[allow(clippy::too_many_arguments)]
    fn coefficients_for_input<P, E: From<DepositionError>>(
        &self,
        route: &FrictionHeatRoute,
        state: &DepositionSimulation<'m>,
        interface_temperature_k: f64,
        input: FrictionHeatInput<P>,
        with_entropy: bool,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<StepCoefficients<TrialPayload<P>>, E> {
        let dt = self.duration_s();
        let FrictionHeatInput { dissipated_j, conductivity_w_m_k, transport } = input;
        if !dissipated_j.is_finite() || dissipated_j < 0.0 {
            return Err(invalid("only finite nonnegative exported friction work may heat bodies").into());
        }
        let StepCoefficients { mut sources_w, mut conductances_w_k, payload } = transport;
        if conductances_w_k.insert(route.patch.link.id, 0.0).is_some() {
            return Err(invalid("friction interface conduction has two competing owners").into());
        }
        let external_heat_j = finite(
            dt * state.external_step_power_w(&sources_w, &conductances_w_k, &mut checkpoint)?,
        )?;
        let conductances = route.conductances(conductivity_w_m_k)?;
        let mut body_heat_j = [0.0; 2];
        let mut interface_entropy = if with_entropy {
            Some(entropy::InterfaceEntropy::new(state, &sources_w, dt, &mut checkpoint)?)
        } else { None };
        for (i, id) in route.cells().into_iter().enumerate() {
            let body_temperature_k = state.temperature_k(id)?;
            let power = finite(conductances[i] * (interface_temperature_k - body_temperature_k))?;
            body_heat_j[i] = finite(dt * power)?;
            if power != 0.0 && body_heat_j[i] == 0.0 {
                return Err(DepositionError::Arithmetic.into());
            }
            if let Some(entropy) = &mut interface_entropy {
                entropy.add_transfer(body_heat_j[i], body_temperature_k, interface_temperature_k)?;
            }
            let source = sources_w.entry(id).or_insert(0.0);
            *source = finite(*source + power)?;
        }
        let interface_residual_j = finite(finite(body_heat_j[0] + body_heat_j[1])? - dissipated_j)?;
        Ok(StepCoefficients {
            sources_w,
            conductances_w_k,
            payload: (payload, conductances, body_heat_j, dissipated_j,
                interface_residual_j, external_heat_j, interface_entropy),
        })
    }

    fn finish_trial<P>(
        &self,
        trial: BackwardEulerTrial<'m, TrialPayload<P>>,
        interface_temperature_k: f64,
    ) -> Result<FrictionHeatTrial<'m, P>, DepositionError> {
        let (_, conductances_w_k, body_heat_j, dissipated_j, interface_residual_j, external_heat_j,
            interface_entropy) =
            *trial.payload();
        let entropy = interface_entropy.map(|interface| {
            let thermal = trial.entropy_audit().ok_or(DepositionError::Arithmetic)?;
            FrictionHeatEntropy::from_trial(thermal, interface, dissipated_j,
                interface_temperature_k, interface_residual_j)
        }).transpose()?;
        let report = FrictionHeatReport {
            interface_temperature_k,
            conductances_w_k,
            body_heat_j,
            dissipated_j,
            interface_residual_j,
            external_heat_j,
            thermal: *trial.audit(),
            total_energy_residual_j: finite(
                trial.audit().stored_energy_change_j - external_heat_j - dissipated_j,
            )?,
            entropy,
        };
        let mut residual = trial.residual_j().to_vec();
        residual.push(interface_residual_j);
        Ok(FrictionHeatTrial {
            trial,
            residual,
            report,
        })
    }
}

type TrialPayload<P> = (P, [f64; 2], [f64; 2], f64, f64, f64, Option<entropy::InterfaceEntropy>);
/// A tentative pair of original owner states; acceptance is the only publication.
pub struct FrictionHeatTrial<'m, P> {
    trial: BackwardEulerTrial<'m, TrialPayload<P>>,
    residual: Vec<f64>,
    report: FrictionHeatReport,
}
impl<'m, P> FrictionHeatTrial<'m, P> {
    /// N cell energy residuals followed by the interface energy residual [J].
    pub fn residual_j(&self) -> &[f64] {
        &self.residual
    }
    pub fn payload(&self) -> &P {
        &self.trial.payload().0
    }
    pub const fn report(&self) -> &FrictionHeatReport {
        &self.report
    }
    pub fn accept(
        self,
        tolerance: FrictionHeatTolerances,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, P, FrictionHeatReport), DepositionError> {
        checkpoint()?;
        if [
            tolerance.cell_j,
            tolerance.thermal_j,
            tolerance.interface_j,
            tolerance.total_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid(
                "friction heat tolerances must be finite nonnegative joules",
            ));
        }
        if self.report.interface_residual_j.abs() > tolerance.interface_j {
            return Err(DepositionError::EnergyBalance(
                self.report.interface_residual_j,
            ));
        }
        if self.report.total_energy_residual_j.abs() > tolerance.total_j {
            return Err(DepositionError::EnergyBalance(
                self.report.total_energy_residual_j,
            ));
        }
        let (thermal, (payload, ..), _) =
            self.trial
                .accept(tolerance.cell_j, tolerance.thermal_j, &mut checkpoint)?;
        checkpoint()?;
        Ok((thermal, payload, self.report))
    }
}

/// Advance the original tangential owner at a coupled interface temperature.
/// The returned candidate can be the low-level heat callback's staged payload.
/// The law owns temperature dependence and any thermodynamic storage; this
/// helper neither refreshes a loaded elastic law nor invents layer heat.
pub fn tangential_heat_candidate<L: TangentialLaw>(
    law: &L,
    accepted: &L::State,
    mut input: TangentialStepInput,
    interface_temperature_k: f64,
    work_tolerance_j: f64,
) -> Result<TangentialStep<L::State>, TriboError> {
    if !work_tolerance_j.is_finite()
        || work_tolerance_j < 0.0
        || !interface_temperature_k.is_finite()
        || interface_temperature_k <= 0.0
    {
        return Err(TriboError::InvalidInput {
            field: "friction_heat_work_tolerance_j",
        });
    }
    input.contact_temperature_k = Some(interface_temperature_k);
    let step = law.advance(accepted, &input)?;
    let work = step.traction_n[0] * step.slip_increment_m[0]
        + step.traction_n[1] * step.slip_increment_m[1];
    let storage = law.stored_energy_j(&step.state) - law.stored_energy_j(accepted);
    let residual = storage + step.dissipated_j - work;
    if !work.is_finite()
        || !storage.is_finite()
        || !residual.is_finite()
        || !step.dissipated_j.is_finite()
        || step.dissipated_j < 0.0
        || !step.stored_change_j.is_finite()
        || !step.closure_residual_j.is_finite()
        || residual.abs() > work_tolerance_j
        || (storage - step.stored_change_j).abs() > work_tolerance_j
    {
        return Err(TriboError::InvalidInput {
            field: "friction_heat_constitutive_work_balance",
        });
    }
    Ok(step)
}

