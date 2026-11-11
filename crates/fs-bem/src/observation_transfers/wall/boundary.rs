//! Local surface laws evaluated by their original mechanical/material owners.
use super::*;
use crate::helmholtz::bc::{PorousLayersBoundary, PorousRigidBoundary, PorousTransmittingBoundary};
use fs_material::fluid_acoustics::{FluidAcousticState, thermoviscous::ThermoviscousFluidState};
use fs_material::porous::{
    PorousFluidLaw,
    layers::{
        Layer, MAX_LAYERS, RigidBackedLayers, TransmittingLayers,
        absorption::{LayerAbsorptionBudget, LayerAbsorptionError},
    },
};
use fs_material::relaxation::RelaxationState;
use fs_phs::WallPin;

/// An owned layer in a locally reacting treatment, ordered from
/// the exposed face toward its backing. The original material owner retains
/// pressure/normal-velocity continuity and coherent propagation between layers.
#[derive(Clone, Debug)]
pub enum LocalWallLayer {
    Porous {
        law: PorousFluidLaw,
        thickness_m: f64,
    },
    /// Lossless gap containing the wall's coherent saturating/exterior fluid.
    FluidGap { thickness_m: f64 },
    /// Lossless fluid with its own original density and bulk modulus.
    FluidLayer {
        fluid: FluidAcousticState,
        thickness_m: f64,
    },
    /// Free gas with its original molecular dispersion and compression loss.
    /// No pore drag or viscothermal propagation is inferred. The layer owner
    /// retains the exact lossless limit when no reservoirs are active.
    RelaxingFluid {
        relaxation: RelaxationState,
        thickness_m: f64,
    },
    /// Original porous skeleton with molecular calorics coupled to pore heat
    /// diffusion. Its equilibrium gas must match the wall's incident/pore
    /// snapshot, and the skeleton's original frequency band remains binding.
    RelaxingPorous {
        law: PorousFluidLaw,
        relaxation: RelaxationState,
        thickness_m: f64,
    },
}

impl LocalWallLayer {
    fn source(&self) -> Layer<'_> {
        match self {
            Self::Porous { law, thickness_m } => Layer::Porous {
                law,
                thickness_m: *thickness_m,
            },
            Self::FluidGap { thickness_m } => Layer::FluidGap {
                thickness_m: *thickness_m,
            },
            Self::FluidLayer { fluid, thickness_m } => Layer::FluidLayer {
                fluid,
                thickness_m: *thickness_m,
            },
            Self::RelaxingFluid {
                relaxation,
                thickness_m,
            } => Layer::RelaxingFluid {
                relaxation,
                thickness_m: *thickness_m,
            },
            Self::RelaxingPorous {
                law,
                relaxation,
                thickness_m,
            } => Layer::RelaxingPorous {
                law,
                relaxation,
                thickness_m: *thickness_m,
            },
        }
    }
}

/// A passive locally reacting surface, in specific-impedance units [Pa s/m].
/// These are boundary models, not a new shell or moving porous-frame solver.
#[derive(Clone, Debug)]
pub enum LocalWall {
    Rigid,
    /// Original surface mass, stiffness and resistance, with exp(-i omega t).
    Mechanical(WallPin),
    /// Original JCA(L) source and one coherent saturating-fluid snapshot.
    /// The BEM adapter performs its existing time-convention conversion once.
    PorousRigid {
        law: PorousFluidLaw,
        fluid: ThermoviscousFluidState,
        thickness_m: f64,
    },
    /// Original coherent multilayer response, used at NORMAL incidence by
    /// the existing local Robin adapter. No lateral propagation, curvature
    /// correction or moving frame is inferred. At most MAX_LAYERS are admitted;
    /// zero thickness is an identity, and an empty stack is rigid.
    LayeredRigid {
        fluid: ThermoviscousFluidState,
        layers: Vec<LocalWallLayer>,
    },
    /// Ordered layers terminated by an original lossless receiving fluid.
    /// The normal-incidence Robin law represents a local outgoing planar load,
    /// not a second 3D receiving field. Accepted work includes transmission;
    /// only the explicit power-partition query evaluates irreversible loss.
    LayeredTransmitting {
        fluid: ThermoviscousFluidState,
        receiving_fluid: FluidAcousticState,
        layers: Vec<LocalWallLayer>,
    },
}

/// Coordinates for evaluating the original local wall work. Near a pressure
/// node, x obeys p = (1+r)*x/2 and rho*c*v_in = (1-r)*x/2. Its bounded inverse
/// uses both retained traces, without dividing by the vanishing pressure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum WallWorkCoordinates {
    Pressure(C64),
    Reflection { reflection: C64, conductance: f64 },
}

impl WallWorkCoordinates {
    pub(super) fn accepted_conductance(self) -> f64 {
        match self {
            Self::Pressure(beta) => beta.re,
            Self::Reflection { conductance, .. } => conductance,
        }
    }
}

#[derive(Clone, Copy)]
pub(super) struct WallPowerCoefficients {
    pub(super) coordinates: WallWorkCoordinates,
    pub(super) irreversible_conductance: f64,
    pub(super) transmitted_conductance: f64,
}

/// An immutable material table and one table index per physical boundary panel.
/// Each distinct law is evaluated once per frequency, irrespective of panel or
/// source count. Every declared law must be used; at most 256 are admitted.
#[derive(Clone, Debug)]
pub struct LocalBoundary {
    medium: Medium,
    laws: Vec<LocalWall>,
    panel_laws: Vec<usize>,
    supported_hz: (f64, f64),
    hash: ContentHash,
}

impl LocalBoundary {
    /// Bind physical panel assignments without defaults, interpolation or
    /// extrapolation. Porous exterior density and speed must match the actual
    /// saturating fluid; the material's positive-frequency support is retained.
    pub fn try_new(
        medium: Medium,
        laws: Vec<LocalWall>,
        panel_laws: Vec<usize>,
    ) -> Result<Self, TransferError> {
        Self::with_panel_cap(medium, laws, panel_laws, helmholtz::MAX_DENSE_PANELS)
    }

    /// The same physical material binding for an explicit iterative observer.
    /// Only the dense allocation ceiling changes; constitutive and geometric
    /// source support remain the original owners' responsibilities.
    pub fn try_new_iterative(
        medium: Medium,
        laws: Vec<LocalWall>,
        panel_laws: Vec<usize>,
    ) -> Result<Self, TransferError> {
        Self::with_panel_cap(medium, laws, panel_laws, crate::panel3d::MAX_SURFACE_PANELS)
    }

    fn with_panel_cap(
        medium: Medium,
        laws: Vec<LocalWall>,
        panel_laws: Vec<usize>,
        maximum_panels: usize,
    ) -> Result<Self, TransferError> {
        if laws.is_empty()
            || laws.len() > 256
            || panel_laws.is_empty()
            || panel_laws.len() > maximum_panels
            || !medium.density.is_finite()
            || medium.density <= 0.0
            || !medium.sound_speed.is_finite()
            || medium.sound_speed <= 0.0
            || !(medium.density * medium.sound_speed).is_finite()
        {
            return Err(TransferError::Invalid(
                "invalid local boundary dimensions or medium",
            ));
        }
        let mut used = vec![false; laws.len()];
        for &index in &panel_laws {
            let Some(flag) = used.get_mut(index) else {
                return Err(TransferError::Invalid(
                    "local boundary panel names an absent law",
                ));
            };
            *flag = true;
        }
        if used.iter().any(|used| !used) {
            return Err(TransferError::Invalid(
                "local boundary contains an unused law",
            ));
        }
        let mut supported_hz = (0.0_f64, f64::INFINITY);
        for law in &laws {
            match law {
                LocalWall::Rigid => {}
                LocalWall::Mechanical(pin) => {
                    // This original owner admits the physical coefficients.
                    fs_phs::wall_specific_impedance(pin, 1.0)?;
                }
                LocalWall::PorousRigid {
                    law,
                    fluid,
                    thickness_m,
                } => {
                    let band = law.reference_card().valid_hz;
                    let prepared = PorousRigidBoundary::new(
                        law,
                        fluid,
                        core::f64::consts::TAU * band.0,
                        *thickness_m,
                    )?;
                    if prepared.medium() != medium {
                        return Err(TransferError::Invalid(
                            "porous and exterior fluid coefficients differ",
                        ));
                    }
                    supported_hz.0 = supported_hz.0.max(band.0);
                    supported_hz.1 = supported_hz.1.min(band.1);
                }
                LocalWall::LayeredRigid { fluid, layers }
                | LocalWall::LayeredTransmitting { fluid, layers, .. } => {
                    if layers.len() > MAX_LAYERS {
                        return Err(TransferError::Invalid(
                            "local wall exceeds the material layer ceiling",
                        ));
                    }
                    let mut band = (0.0_f64, f64::INFINITY);
                    for layer in layers {
                        if let LocalWallLayer::Porous { law, .. }
                        | LocalWallLayer::RelaxingPorous { law, .. } = layer
                        {
                            let source_band = law.reference_card().valid_hz;
                            band.0 = band.0.max(source_band.0);
                            band.1 = band.1.min(source_band.1);
                        }
                    }
                    if band.0 >= band.1 {
                        return Err(TransferError::Invalid(
                            "local wall layers have no common frequency band",
                        ));
                    }
                    // Validate thicknesses and original sources together at
                    // one supported frequency, including zero-thickness laws.
                    // Robin pressure-node refusals belong to the actual query,
                    // not an arbitrary validation frequency for lossless gaps.
                    let probe_hz = if band.0 > 0.0 { band.0 } else { 1.0 };
                    let omega = core::f64::consts::TAU * probe_hz;
                    let acoustic = if let LocalWall::LayeredTransmitting { receiving_fluid, .. } = law {
                        prepare_transmitting(fluid, receiving_fluid, layers, omega)?
                            .prepared_layers().fluid().acoustic()
                    } else {
                        prepare_layers(fluid, layers, omega)?.fluid().acoustic()
                    };
                    if acoustic.density_kg_m3() != medium.density
                        || acoustic.sound_speed_m_s() != medium.sound_speed
                    {
                        return Err(TransferError::Invalid(
                            "layered wall and exterior fluid coefficients differ",
                        ));
                    }
                    supported_hz.0 = supported_hz.0.max(band.0);
                    supported_hz.1 = supported_hz.1.min(band.1);
                }
            }
        }
        if supported_hz.0 >= supported_hz.1 {
            return Err(TransferError::Invalid(
                "local wall materials have no common frequency band",
            ));
        }
        let mut boundary = Self {
            medium,
            laws,
            panel_laws,
            supported_hz,
            hash: ContentHash([0; 32]),
        };
        boundary.hash = boundary.compute_hash();
        Ok(boundary)
    }

    pub const fn medium(&self) -> Medium {
        self.medium
    }
    pub fn panel_count(&self) -> usize {
        self.panel_laws.len()
    }
    /// Common source support [Hz]; infinity means no finite constitutive upper
    /// cutoff, not that mesh resolution or floating-point range is unlimited.
    pub const fn supported_band_hz(&self) -> (f64, f64) {
        self.supported_hz
    }
    pub const fn content_hash(&self) -> ContentHash {
        self.hash
    }

    pub(in crate::observation_transfers) fn rigid_at(&self, panel: usize) -> bool {
        matches!(self.laws[self.panel_laws[panel]], LocalWall::Rigid)
    }

    pub(super) fn admit_band(
        &self,
        highest_omega: f64,
        lower: Option<LowerBandWindow>,
    ) -> Result<(), TransferError> {
        let tau = core::f64::consts::TAU;
        let lowest = lower.map_or(0.0, |window| window.stop_hz);
        if lowest < self.supported_hz.0
            || highest_omega > tau * self.supported_hz.1
            || lower.is_some_and(|window| tau * window.pass_hz >= highest_omega)
        {
            return Err(TransferError::Invalid(
                "observation band exceeds local-wall source support",
            ));
        }
        Ok(())
    }

    pub(super) fn admit_admittance_band(
        &self,
        highest_omega: f64,
        lower: Option<LowerBandWindow>,
    ) -> Result<(), TransferError> {
        self.admit_band(highest_omega, lower)?;
        self.admit_resonances(
            core::f64::consts::TAU * lower.map_or(0.0, |window| window.stop_hz),
            highest_omega,
            false,
        )
    }

    /// Constitutive support for owners that retain finite reflection coordinates
    /// through a pressure-release node. This does not relax source bandwidth.
    pub(in crate::observation_transfers) fn admit_physical_reflection_band_hz(
        &self,
        band: [f64; 2],
    ) -> Result<(), TransferError> {
        if !band.iter().all(|f| f.is_finite() && *f > 0.0)
            || band[0] > band[1]
            || band[0] < self.supported_hz.0
            || band[1] > self.supported_hz.1
        {
            return Err(TransferError::Invalid(
                "physical load band exceeds local-wall source support",
            ));
        }
        Ok(())
    }

    pub(in crate::observation_transfers) fn admit_physical_band_hz(
        &self,
        band: [f64; 2],
    ) -> Result<(), TransferError> {
        self.admit_physical_reflection_band_hz(band)?;
        let tau = core::f64::consts::TAU;
        self.admit_resonances(tau * band[0], tau * band[1], true)
    }

    fn admit_resonances(
        &self,
        lowest: f64,
        highest: f64,
        include_lowest: bool,
    ) -> Result<(), TransferError> {
        for law in &self.laws {
            if let LocalWall::Mechanical(pin) = law
                && pin.resistance == 0.0
            {
                let resonance = det::sqrt(pin.stiffness_per_area) / det::sqrt(pin.surface_density);
                if (resonance > lowest || (include_lowest && resonance == lowest))
                    && resonance <= highest
                {
                    return Err(TransferError::Invalid(
                        "lossless wall resonance requires a pressure-release boundary owner",
                    ));
                }
            }
        }
        Ok(())
    }

    /// Dimensionless beta=rho*c/Z in exp(-i omega t), in original panel order.
    /// Native constitutive refusals propagate. A zero-impedance lossless
    /// resonance refuses; admittance is never clipped or replaced with rigid.
    pub fn normalized_admittance(
        &self,
        omega: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, TransferError> {
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        if !omega.is_finite() || omega <= 0.0 {
            return Err(TransferError::Invalid(
                "local wall frequency must be positive and finite",
            ));
        }
        let rho_c = self.medium.density * self.medium.sound_speed;
        let mut coefficients = Vec::with_capacity(self.laws.len());
        for law in &self.laws {
            if cancelled() {
                return Err(TransferError::Cancelled);
            }
            let beta = law.admittance(rho_c, omega)?;
            coefficients.push(beta);
        }
        Ok(self
            .panel_laws
            .iter()
            .map(|&index| coefficients[index])
            .collect())
    }

    /// Normal-incidence pressure reflection in exp(-i omega t), in original
    /// panel order. The bounded coefficients retain zero-impedance mechanical
    /// resonances and layered pressure nodes without an admittance division.
    /// Every original source frequency domain and constitutive refusal remains.
    pub fn pressure_reflection(
        &self,
        omega: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, TransferError> {
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        if !omega.is_finite() || omega <= 0.0 {
            return Err(TransferError::Invalid("local wall frequency must be positive and finite"));
        }
        let rho_c = self.medium.density * self.medium.sound_speed;
        let mut coefficients = Vec::with_capacity(self.laws.len());
        for law in &self.laws {
            if cancelled() {
                return Err(TransferError::Cancelled);
            }
            let (reflection, _) = law.reflected_work(rho_c, omega)?;
            coefficients.push(reflection);
        }
        let mut panels = Vec::with_capacity(self.panel_laws.len());
        for (panel, &index) in self.panel_laws.iter().enumerate() {
            if panel % 256 == 0 && cancelled() {
                return Err(TransferError::Cancelled);
            }
            panels.push(coefficients[index]);
        }
        Ok(panels)
    }

    /// Whether any physical law is a mechanical pin. Only then does the near-node reflection route
    /// use native pressure/flow coordinates; every other boundary keeps the reflection pair it always had.
    pub(super) fn has_mechanical_law(&self) -> bool {
        self.laws.iter().any(|law| matches!(law, LocalWall::Mechanical(_)))
    }

    /// Native `[p, rho c v_inward]` coordinates per panel, in original panel order, for the reflection
    /// constructors' `from_pressure_flow`. A mechanical law gives `[z, rho c]` divided by their common scale
    /// and never forms the reflection `r`; every other law gives `[(1 + r)/2, (1 - r)/2]`, the pair the
    /// reflection constructors form. Both are the same boundary equation up to a common complex amplitude.
    /// Every constitutive, frequency and passivity refusal of `pressure_reflection` remains.
    ///
    /// Why not through `r`: for a lossless pin at detuning X from its resonance, r = -1 + i theta with
    /// theta = 2|X|/(rho c). The true real part is -1 + theta^2/2, which for theta below about 1e-8 is under
    /// half an ulp of 1, so `r.re` rounds to exactly -1 and |r|^2 = 1 + theta^2: a hair active. The pair
    /// `[(1 + r)/2, (1 - r)/2]` then has conductance -theta^2/4, that is -theta/2 of |p||w|, and the passive-work
    /// check (-1e-10) refused the wall for detunings from about 2e-11 (theta/2 above 1e-10) up to about 3e-9,
    /// where |beta| falls below 1/sqrt(eps) and the admittance chart takes over. Here the pair has
    /// `Re(p conj w) = z.re rho c / scale^2`, exactly zero for a lossless pin.
    pub(super) fn pressure_flow_coordinates(
        &self,
        omega: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<[C64; 2]>, TransferError> {
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        if !omega.is_finite() || omega <= 0.0 {
            return Err(TransferError::Invalid("local wall frequency must be positive and finite"));
        }
        let rho_c = self.medium.density * self.medium.sound_speed;
        let mut coordinates = Vec::with_capacity(self.laws.len());
        for law in &self.laws {
            if cancelled() {
                return Err(TransferError::Cancelled);
            }
            coordinates.push(law.pressure_flow_coordinates(rho_c, omega)?);
        }
        let mut panels = Vec::with_capacity(self.panel_laws.len());
        for (panel, &index) in self.panel_laws.iter().enumerate() {
            if panel % 256 == 0 && cancelled() {
                return Err(TransferError::Cancelled);
            }
            panels.push(coordinates[index]);
        }
        Ok(panels)
    }

    /// Retain ordinary admittance arithmetic where well scaled, and the
    /// original bounded wall law through mechanical/layered pressure nodes.
    /// Distinct material laws are queried once before scattering to panels.
    pub(super) fn work_coordinates(
        &self,
        omega: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<WallWorkCoordinates>, TransferError> {
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        if !omega.is_finite() || omega <= 0.0 {
            return Err(TransferError::Invalid("local wall frequency must be positive and finite"));
        }
        let rho_c = self.medium.density * self.medium.sound_speed;
        let mut coefficients = Vec::with_capacity(self.laws.len());
        for law in &self.laws {
            if cancelled() {
                return Err(TransferError::Cancelled);
            }
            coefficients.push(law.work_coordinates(rho_c, omega)?);
        }
        let mut panels = Vec::with_capacity(self.panel_laws.len());
        for (panel, &index) in self.panel_laws.iter().enumerate() {
            if panel % 256 == 0 && cancelled() {
                return Err(TransferError::Cancelled);
            }
            panels.push(coefficients[index]);
        }
        Ok(panels)
    }

    /// Query expensive positive layer-loss integrals only for a requested
    /// physical power partition. The allowance counts declared layers across
    /// distinct transmitting wall laws, including zero-thickness identities.
    /// Each native query admits that count before its propagation/integration
    /// passes and polls cancellation per layer. No panel repeats a wall query.
    pub(super) fn power_coefficients(
        &self,
        omega: f64,
        maximum_absorption_layers: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Vec<WallPowerCoefficients>, usize), TransferError> {
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        if !omega.is_finite() || omega <= 0.0 {
            return Err(TransferError::Invalid("local wall frequency must be positive and finite"));
        }
        let layer_count = self.laws.iter().try_fold(0_usize, |total, law| {
            let count = match law {
                LocalWall::LayeredTransmitting { layers, .. } => layers.len(),
                _ => 0,
            };
            total.checked_add(count).ok_or(TransferError::Invalid("wall loss extent overflow"))
        })?;
        if layer_count > maximum_absorption_layers {
            return Err(TransferError::Invalid("wall loss exceeds declared material-layer budget"));
        }
        let rho_c = self.medium.density * self.medium.sound_speed;
        let mut coefficients = Vec::with_capacity(self.laws.len());
        for law in &self.laws {
            if cancelled() {
                return Err(TransferError::Cancelled);
            }
            let coefficient = if let LocalWall::LayeredTransmitting { fluid, receiving_fluid, layers } = law {
                let prepared = prepare_transmitting(fluid, receiving_fluid, layers, omega)?;
                let power = prepared.normal_surface_power(
                    LayerAbsorptionBudget { max_layers: layers.len() },
                    &mut cancelled,
                );
                match power {
                    Ok(power) if well_scaled_admittance(power.normalized_admittance) => {
                        WallPowerCoefficients {
                            coordinates: WallWorkCoordinates::Pressure(power.normalized_admittance.conj()),
                            irreversible_conductance: power.absorbed_conductance,
                            transmitted_conductance: power.transmitted_conductance,
                        }
                    }
                    Err(LayerAbsorptionError::Cancelled) => return Err(TransferError::Cancelled),
                    _ => {
                        // Native positive layer integrals in the incident-wave
                        // basis stay finite at a total-pressure node. x = 2*A,
                        // hence each unit-incident-power fraction contributes
                        // fraction/4 to |x|^2/(2*rho*c). Never infer weak heat
                        // from 1-|r|^2-T or from the BEM work residual.
                        let (reflection, conductance) = law.reflected_work(rho_c, omega)?;
                        let power = prepared.layer_absorption(
                            1.0, 1.0,
                            LayerAbsorptionBudget { max_layers: layers.len() },
                            &mut cancelled,
                        ).map_err(absorption_error)?;
                        WallPowerCoefficients {
                            coordinates: WallWorkCoordinates::Reflection { reflection, conductance },
                            irreversible_conductance: 0.25 * power.absorbed_fraction,
                            transmitted_conductance: 0.25 * power.transmitted_fraction,
                        }
                    }
                }
            } else {
                let coordinates = law.work_coordinates(rho_c, omega)?;
                WallPowerCoefficients {
                    coordinates,
                    irreversible_conductance: coordinates.accepted_conductance(),
                    transmitted_conductance: 0.0,
                }
            };
            coefficients.push(coefficient);
        }
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        Ok((self.panel_laws.iter().map(|&i| coefficients[i]).collect(), layer_count))
    }

    fn compute_hash(&self) -> ContentHash {
        let mut hash = DomainHasher::new("fs-bem:local-boundary:v1");
        hash.update(&(self.laws.len() as u64).to_le_bytes());
        hash.update(&(self.panel_laws.len() as u64).to_le_bytes());
        hash_scalars(&mut hash, &[self.medium.density, self.medium.sound_speed]);
        for law in &self.laws {
            match law {
                LocalWall::Rigid => hash.update(&[0]),
                LocalWall::Mechanical(pin) => {
                    hash.update(&[1]);
                    hash_scalars(
                        &mut hash,
                        &[pin.surface_density, pin.stiffness_per_area, pin.resistance],
                    );
                }
                LocalWall::PorousRigid {
                    law,
                    fluid,
                    thickness_m,
                } => {
                    hash.update(&[2]);
                    let card = law.reference_card();
                    for text in [card.name, card.provenance] {
                        hash.update(&(text.len() as u64).to_le_bytes());
                        hash.update(text.as_bytes());
                    }
                    hash_scalars(
                        &mut hash,
                        &[
                            card.porosity,
                            card.flow_resistivity,
                            card.tortuosity,
                            card.viscous_length_m,
                            card.thermal_length_m,
                            card.valid_hz.0,
                            card.valid_hz.1,
                            law.reference_dynamic_viscosity_pa_s(),
                            *thickness_m,
                            fluid.acoustic().density_kg_m3(),
                            fluid.acoustic().sound_speed_m_s(),
                            fluid.acoustic().adiabatic_bulk_modulus_pa(),
                            fluid.transport().density_kg_m3(),
                            fluid.transport().specific_heat_j_kg_k(),
                            fluid.transport().thermal_conductivity_w_m_k(),
                            fluid.transport().dynamic_viscosity_pa_s(),
                            fluid.heat_capacity_ratio(),
                            fluid.specific_heat_cv_j_kg_k(),
                        ],
                    );
                    match card.thermal_permeability_m2 {
                        None => hash.update(&[0]),
                        Some(value) => {
                            hash.update(&[1]);
                            hash_scalars(&mut hash, &[value]);
                        }
                    }
                }
                LocalWall::LayeredRigid { fluid, layers }
                | LocalWall::LayeredTransmitting { fluid, layers, .. } => {
                    if let LocalWall::LayeredTransmitting { receiving_fluid, .. } = law {
                        hash.update(&[4]);
                        hash_scalars(&mut hash, &[
                            receiving_fluid.density_kg_m3(),
                            receiving_fluid.sound_speed_m_s(),
                            receiving_fluid.adiabatic_bulk_modulus_pa(),
                        ]);
                    } else {
                        hash.update(&[3]);
                    }
                    hash_scalars(
                        &mut hash,
                        &[
                            fluid.acoustic().density_kg_m3(),
                            fluid.acoustic().sound_speed_m_s(),
                            fluid.acoustic().adiabatic_bulk_modulus_pa(),
                            fluid.transport().density_kg_m3(),
                            fluid.transport().specific_heat_j_kg_k(),
                            fluid.transport().thermal_conductivity_w_m_k(),
                            fluid.transport().dynamic_viscosity_pa_s(),
                            fluid.heat_capacity_ratio(),
                            fluid.specific_heat_cv_j_kg_k(),
                        ],
                    );
                    hash.update(&(layers.len() as u64).to_le_bytes());
                    for layer in layers {
                        match layer {
                            LocalWallLayer::Porous { law, thickness_m } => {
                                hash.update(&[0]);
                                hash_scalars(&mut hash, &[*thickness_m]);
                                hash_porous_source(&mut hash, law);
                            }
                            LocalWallLayer::FluidGap { thickness_m } => {
                                hash.update(&[1]);
                                hash_scalars(&mut hash, &[*thickness_m]);
                            }
                            LocalWallLayer::FluidLayer { fluid, thickness_m } => {
                                hash.update(&[2]);
                                hash_scalars(
                                    &mut hash,
                                    &[
                                        *thickness_m,
                                        fluid.density_kg_m3(),
                                        fluid.sound_speed_m_s(),
                                        fluid.adiabatic_bulk_modulus_pa(),
                                    ],
                                );
                            }
                            LocalWallLayer::RelaxingFluid {
                                relaxation,
                                thickness_m,
                            } => {
                                hash.update(&[3]);
                                hash_scalars(&mut hash, &[*thickness_m]);
                                hash_relaxation_source(&mut hash, relaxation);
                            }
                            LocalWallLayer::RelaxingPorous {
                                law,
                                relaxation,
                                thickness_m,
                            } => {
                                hash.update(&[4]);
                                hash_scalars(&mut hash, &[*thickness_m]);
                                hash_porous_source(&mut hash, law);
                                hash_relaxation_source(&mut hash, relaxation);
                            }
                        }
                    }
                }
            }
        }
        for &index in &self.panel_laws {
            hash.update(&(index as u64).to_le_bytes());
        }
        hash.finalize()
    }
}

fn prepare_layers<'a>(
    fluid: &'a ThermoviscousFluidState,
    layers: &'a [LocalWallLayer],
    omega: f64,
) -> Result<RigidBackedLayers<'a>, TransferError> {
    let sources = layer_sources(layers)?;
    Ok(RigidBackedLayers::new(
        fluid,
        omega,
        &sources[..layers.len()],
    )?)
}

impl LocalWall {
    fn work_coordinates(&self, rho_c: f64, omega: f64) -> Result<WallWorkCoordinates, TransferError> {
        if let Ok(beta) = self.admittance(rho_c, omega)
            && well_scaled_admittance(beta)
        {
            return Ok(WallWorkCoordinates::Pressure(beta));
        }
        // The original source still admits the frequency and computes the
        // reflection. A failed admittance is not replaced by a rigid boundary.
        let (reflection, conductance) = self.reflected_work(rho_c, omega)?;
        Ok(WallWorkCoordinates::Reflection { reflection, conductance })
    }

    fn pressure_flow_coordinates(&self, rho_c: f64, omega: f64) -> Result<[C64; 2], TransferError> {
        // The reflection is formed first, so every refusal (nonfinite, active, constitutive) is the old one.
        let (reflection, _) = self.reflected_work(rho_c, omega)?;
        if let Self::Mechanical(pin) = self {
            let z = fs_phs::wall_specific_impedance(pin, omega)?;
            let scale = z.re.abs().max(z.im.abs()).max(rho_c);
            return Ok([C64::new(z.re / scale, z.im / scale), C64::from_re(rho_c / scale)]);
        }
        Ok([(C64::ONE + reflection).scale(0.5), (C64::ONE - reflection).scale(0.5)])
    }

    fn reflected_work(&self, rho_c: f64, omega: f64) -> Result<(C64, f64), TransferError> {
        let (reflection, conductance) = match self {
            Self::Rigid => (C64::ONE, 0.0),
            Self::Mechanical(pin) => {
                let z = fs_phs::wall_specific_impedance(pin, omega)?;
                if !finite(z) {
                    return Err(TransferError::Invalid("wall impedance is unrepresentable"));
                }
                let scale = z.re.abs().max(z.im.abs()).max(rho_c);
                let z = C64::new(z.re / scale, z.im / scale);
                let ambient = C64::from_re(rho_c / scale);
                let denominator = z + ambient;
                // Original resistance retains weak heat even when |r| rounds
                // to one. A lossless mechanical resonance is exactly zero heat.
                ((z - ambient) / denominator, z.re * ambient.re / denominator.norm_sq())
            }
            Self::PorousRigid { law, fluid, thickness_m } => {
                let response = RigidBackedLayers::new(
                    fluid, omega, &[Layer::Porous { law, thickness_m: *thickness_m }],
                )?.reflection(1.0)?;
                (response.pressure_reflection.conj(), 0.25 * response.absorbed_fraction)
            }
            Self::LayeredRigid { fluid, layers } => {
                let response = prepare_layers(fluid, layers, omega)?.reflection(1.0)?;
                (response.pressure_reflection.conj(), 0.25 * response.absorbed_fraction)
            }
            Self::LayeredTransmitting { fluid, receiving_fluid, layers } => {
                let response = prepare_transmitting(fluid, receiving_fluid, layers, omega)?.scattering(1.0)?;
                (response.pressure_reflection.conj(),
                    0.25 * (response.absorbed_fraction + response.transmitted_fraction))
            }
        };
        if !finite(reflection) || reflection.abs() > 1.0 + 64.0 * f64::EPSILON
            || !conductance.is_finite()
        {
            return Err(TransferError::Invalid("local wall reflection is nonfinite or active"));
        }
        Ok((reflection, conductance))
    }

    fn admittance(&self, rho_c: f64, omega: f64) -> Result<C64, TransferError> {
        let beta = match self {
            LocalWall::Rigid => C64::ZERO,
            LocalWall::Mechanical(pin) => {
                let impedance = fs_phs::wall_specific_impedance(pin, omega)?;
                if !finite(impedance) || impedance == C64::ZERO {
                    return Err(TransferError::Invalid(
                        "wall impedance is singular or unrepresentable",
                    ));
                }
                C64::from_re(rho_c) / impedance
            }
            LocalWall::PorousRigid {
                law,
                fluid,
                thickness_m,
            } => PorousRigidBoundary::new(law, fluid, omega, *thickness_m)?
                .normalized_admittance(),
            LocalWall::LayeredRigid { fluid, layers } => {
                let prepared = prepare_layers(fluid, layers, omega)?;
                PorousLayersBoundary::new(&prepared)?.normalized_admittance()?
            }
            LocalWall::LayeredTransmitting { fluid, receiving_fluid, layers } => {
                let prepared = prepare_transmitting(fluid, receiving_fluid, layers, omega)?;
                PorousTransmittingBoundary::new(&prepared)?.normalized_admittance()
            }
        };
        if !finite(beta) || beta.re < 0.0 {
            return Err(TransferError::Invalid(
                "local wall admittance is nonfinite or active",
            ));
        }
        Ok(beta)
    }
}

fn well_scaled_admittance(beta: C64) -> bool {
    finite(beta) && beta.abs() < 1.0 / det::sqrt(f64::EPSILON)
}

fn layer_sources(layers: &[LocalWallLayer]) -> Result<[Layer<'_>; MAX_LAYERS], TransferError> {
    if layers.len() > MAX_LAYERS {
        return Err(TransferError::Invalid(
            "local wall exceeds the material layer ceiling",
        ));
    }
    // Borrow original sources without a second heap allocation per frequency.
    let mut sources = [Layer::FluidGap { thickness_m: 0.0 }; MAX_LAYERS];
    for (source, layer) in sources.iter_mut().zip(layers) {
        *source = layer.source();
    }
    Ok(sources)
}

fn prepare_transmitting<'a>(
    fluid: &'a ThermoviscousFluidState,
    receiving_fluid: &'a FluidAcousticState,
    layers: &'a [LocalWallLayer],
    omega: f64,
) -> Result<TransmittingLayers<'a>, TransferError> {
    let sources = layer_sources(layers)?;
    Ok(TransmittingLayers::new(fluid, receiving_fluid, omega, &sources[..layers.len()])?)
}

fn absorption_error(error: LayerAbsorptionError) -> TransferError {
    match error {
        LayerAbsorptionError::Material(error) => TransferError::Material(error),
        LayerAbsorptionError::Invalid(message) => TransferError::Invalid(message),
        LayerAbsorptionError::Cancelled => TransferError::Cancelled,
    }
}

fn hash_porous_source(hash: &mut DomainHasher, law: &PorousFluidLaw) {
    let card = law.reference_card();
    for text in [card.name, card.provenance] {
        hash.update(&(text.len() as u64).to_le_bytes());
        hash.update(text.as_bytes());
    }
    hash_scalars(
        hash,
        &[
            card.porosity,
            card.flow_resistivity,
            card.tortuosity,
            card.viscous_length_m,
            card.thermal_length_m,
            card.valid_hz.0,
            card.valid_hz.1,
            law.reference_dynamic_viscosity_pa_s(),
        ],
    );
    match card.thermal_permeability_m2 {
        None => hash.update(&[0]),
        Some(value) => {
            hash.update(&[1]);
            hash_scalars(hash, &[value]);
        }
    }
}

fn hash_scalars(hash: &mut DomainHasher, values: &[f64]) {
    for value in values {
        hash.update(&value.to_bits().to_le_bytes());
    }
}

// Original equilibrium state and every evaluated molecular reservoir determine
// the wall response. Equal rho/c alone cannot identify dispersive treatments.
fn hash_relaxation_source(hash: &mut DomainHasher, source: &RelaxationState) {
    let gas = source.gas();
    hash_scalars(
        hash,
        &[
            gas.temperature,
            gas.pressure,
            gas.density,
            gas.sound_speed,
            gas.dynamic_viscosity,
            gas.thermal_conductivity,
            gas.gamma,
            gas.specific_gas_constant,
            gas.specific_heat_cp,
            gas.prandtl,
            gas.characteristic_impedance,
            gas.water_mole_fraction,
            source.frozen_cv_j_kg_k(),
        ],
    );
    hash.update(&(source.processes().len() as u64).to_le_bytes());
    for process in source.processes() {
        let name = process.species.as_str();
        hash.update(&(name.len() as u64).to_le_bytes());
        hash.update(name.as_bytes());
        hash.update(&(process.index as u64).to_le_bytes());
        hash_scalars(
            hash,
            &[process.heat_capacity_j_kg_k, process.relaxation_time_s],
        );
    }
}
