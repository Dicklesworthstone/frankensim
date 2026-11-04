//! Causal coupled reflection from an original smooth lossless duct.
//!
//! The spatial source is fs-duct's full finite-mode Galerkin image of its
//! continuous Riccati operator. Every state and every inlet mode enters ONE
//! shared exact multiport scattering flow. For each pressure wave a held
//! throughout an interval, U=(2a-P)/Z0 and b=P-a. No fitted damping, scalar
//! substitution, discarded higher modes, extra latency or midpoint map.
//!
//! Warm push uses only prepared buffers. It returns interval MEAN b;
//! independently integrated b(t)^2/Z0 remains in the port accounts. The
//! square of the mean is not the reflected energy. This linear ideal-end,
//! lossless spatial rung does not qualify a nonlinear junction, real exterior
//! radiation or a spatial/clock accuracy band without separate refinement.
//! A port mouth exposes both faces; an attached physical admittance load is
//! interconnected before preparation, retaining its full state and losses.

use fs_duct::RadiusProfile;
use fs_duct::modal::ModeCount;
use fs_duct::modal::time_domain::{
    ModalEnd, ModalGalerkinSource, ModalTimeError, assemble_lossless, assemble_lossless_field,
    assemble_lossless_regional,
};
use fs_material::gas::GasState;
use fs_material::gas_field::AxialGasField;
use fs_phs::PortHamiltonian;
use fs_phs::linear_exact::{
    LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, MultiportScatteringStepRecord,
    PreparedMultiportScatteringZoh, ScatteringPortRecord,
};

#[path = "continuous/load.rs"]
mod load;
pub use load::{ModalAdmittanceLoad, ModalLoadSource};

#[path = "continuous/molecular_field.rs"]
mod molecular_field;
mod viscothermal;

/// Original continuous duct and full rational impedance in one exact flow.
pub mod rational;

/// Prepared full-state modal region and its accepted physical history.
///
/// Output pressures and waves use the original mode normalization. The
/// positive real reference at each face is rho*c/S from its original gas
/// and area. This is a wave-coordinate choice, not a matched evanescent
/// characteristic.
pub struct ContinuousModalLine {
    source: ModalGalerkinSource,
    load_source: Option<ModalLoadSource>,
    flow: PreparedMultiportScatteringZoh,
    state: Vec<f64>,
    next: Vec<f64>,
    mean_pressures: Vec<f64>,
    reflected_means: Vec<f64>,
    ports: Vec<ScatteringPortRecord>,
    held_inputs: Vec<f64>,
}

impl ContinuousModalLine {
    /// Prepare original geometry, every retained coordinate and all inlet
    /// modes, plus every mouth mode for `ModalEnd::Ports`. The shared exact
    /// owner admits the complete dense work; no
    /// resolution is reduced to fit. Cancellation is polled during spatial
    /// assembly and after bounded exact-flow preparation.
    ///
    /// # Errors
    /// Original spatial, state-capacity, cancellation and exact-flow refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        profile: &RadiusProfile,
        gas: &GasState,
        modes: ModeCount,
        elements: usize,
        end: ModalEnd,
        duration: f64,
        budget: LinearSplitBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ModalTimeError> {
        let image = assemble_lossless(
            profile,
            gas,
            modes,
            elements,
            end,
            budget.max_modes,
            &mut cancelled,
        )?;
        let (source, system) = image.into_parts();
        Self::prepare(source, system, None, duration, budget, &mut cancelled)
    }

    /// Prepare the same exact coupled flow with one original gas state per
    /// uniform axial element. Density, compressibility and local cutoffs enter
    /// the original Galerkin operator; neither the state nor the response is
    /// averaged or retuned. Both external faces use their own gas reference.
    ///
    /// All regions share equilibrium pressure and remain frozen. Interfaces
    /// lie on element boundaries; thermal evolution, mixing and mean flow are
    /// not supplied by this acoustic image. Every retained modal coordinate,
    /// source state and independent wave-energy account remains inspectable.
    ///
    /// # Errors
    /// Original regional-source, capacity, cancellation and exact-flow refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn new_regional(
        profile: &RadiusProfile,
        gases: &[GasState],
        modes: ModeCount,
        end: ModalEnd,
        duration: f64,
        budget: LinearSplitBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ModalTimeError> {
        let image = assemble_lossless_regional(
            profile,
            gases,
            modes,
            end,
            budget.max_modes,
            &mut cancelled,
        )?;
        let (source, system) = image.into_parts();
        Self::prepare(source, system, None, duration, budget, &mut cancelled)
    }

    /// Prepare the original full modal time evolution in a source-backed
    /// temperature/composition field. The shared assembler evaluates the actual
    /// gas at each physical Gauss point and retains its source receipts; this
    /// runtime uses the resulting consistent mass and stiffness without a fit.
    /// Both external faces retain their own gas/area wave references.
    ///
    /// The background is prescribed, common-pressure and frozen. Material
    /// stations must align with axial element boundaries. The complete state,
    /// all retained modes and original exact-flow resource limits remain active.
    ///
    /// # Errors
    /// Original material, field/grid, capacity, cancellation and exact-flow refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn new_field(
        profile: &RadiusProfile,
        field: &AxialGasField,
        modes: ModeCount,
        elements: usize,
        end: ModalEnd,
        duration: f64,
        budget: LinearSplitBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ModalTimeError> {
        let image = assemble_lossless_field(
            profile,
            field,
            modes,
            elements,
            end,
            budget.max_modes,
            &mut cancelled,
        )?;
        let (source, system) = image.into_parts();
        Self::prepare(source, system, None, duration, budget, &mut cancelled)
    }

    /// Attach a quiescent complete original modal load to every mouth mode, then
    /// prepare ONE exact coupled evolution. The load sees actual mouth pressure
    /// and returns outward duct flow: U_duct=-U_load, P_load=P_duct. Every
    /// mutual term, material memory coordinate and physical loss is retained.
    ///
    /// The source mouth must have the same original radial basis and radius;
    /// both gases must share static pressure. Their density/composition may
    /// differ. Only the inlet modes remain externally driven wave ports.
    /// Accepted state is [duct q; duct M*P; complete original load state].
    ///
    /// # Errors
    /// Mismatched source coordinates/equilibrium, complete-state capacity,
    /// nonquadratic or nonpassive source structure, cancellation or preparation.
    #[allow(clippy::too_many_arguments)]
    pub fn new_loaded(
        profile: &RadiusProfile,
        gas: &GasState,
        modes: ModeCount,
        elements: usize,
        load: ModalAdmittanceLoad,
        duration: f64,
        budget: LinearSplitBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ModalTimeError> {
        Self::validate_load_binding(profile, gas, modes, &load)?;
        let image = assemble_lossless(
            profile,
            gas,
            modes,
            elements,
            ModalEnd::Ports,
            budget.max_modes,
            &mut cancelled,
        )?;
        let (source, duct) = image.into_parts();
        Self::prepare_loaded(source, duct, load, duration, budget, &mut cancelled)
    }

    /// Attach a complete original modal mouth load to a regional gas duct.
    /// Every local density and compressibility, mutual load term and memory
    /// coordinate enters the same exact evolution as [`Self::new_loaded`].
    /// Only the inlet modes remain external incident/reflected wave ports.
    ///
    /// Regional gases remain frozen and share equilibrium pressure. The
    /// load binds to the actual mouth gas and radial coordinates; its own gas
    /// may have a different density/composition at that same static pressure.
    /// No state reduction, independent feedback clock or fitted load is used.
    ///
    /// # Errors
    /// Original regional-source and load-coordinate refusals, complete-state
    /// capacity, nonpassive structure, cancellation or exact-flow preparation.
    #[allow(clippy::too_many_arguments)]
    pub fn new_regional_loaded(
        profile: &RadiusProfile,
        gases: &[GasState],
        modes: ModeCount,
        load: ModalAdmittanceLoad,
        duration: f64,
        budget: LinearSplitBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ModalTimeError> {
        let image = assemble_lossless_regional(
            profile,
            gases,
            modes,
            ModalEnd::Ports,
            budget.max_modes,
            &mut cancelled,
        )?;
        let (source, duct) = image.into_parts();
        let mouth_gas =
            source
                .gas_at_element(source.elements() - 1)
                .ok_or(ModalTimeError::Invalid(
                    "regional modal source requires a mouth gas",
                ))?;
        Self::validate_load_binding(profile, mouth_gas, modes, &load)?;
        Self::prepare_loaded(source, duct, load, duration, budget, &mut cancelled)
    }

    /// Interconnect a complete original mouth load with the source-backed
    /// primitive-field duct before preparing one exact coupled evolution.
    /// The load binds to the actual endpoint gas, radial basis and aperture;
    /// every load state, mutual term and loss enters the same physical clock.
    ///
    /// # Errors
    /// The refusals of [`Self::new_field`] and the original mouth/load binding
    /// and complete coupled-state capacity checks.
    #[allow(clippy::too_many_arguments)]
    pub fn new_field_loaded(
        profile: &RadiusProfile,
        field: &AxialGasField,
        modes: ModeCount,
        elements: usize,
        load: ModalAdmittanceLoad,
        duration: f64,
        budget: LinearSplitBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ModalTimeError> {
        let image = assemble_lossless_field(
            profile,
            field,
            modes,
            elements,
            ModalEnd::Ports,
            budget.max_modes,
            &mut cancelled,
        )?;
        let (source, duct) = image.into_parts();
        Self::validate_load_binding(profile, source.outlet_gas(), modes, &load)?;
        Self::prepare_loaded(source, duct, load, duration, budget, &mut cancelled)
    }

    fn validate_load_binding(
        profile: &RadiusProfile,
        gas: &GasState,
        modes: ModeCount,
        load: &ModalAdmittanceLoad,
    ) -> Result<(), ModalTimeError> {
        let declaration = load.source();
        let radius = profile.outlet_radius();
        if declaration.modes() != modes
            || (declaration.radius_m() - radius).abs() > 128.0 * f64::EPSILON * radius
            || declaration.gas().pressure.to_bits() != gas.pressure.to_bits()
        {
            return Err(ModalTimeError::Invalid(
                "modal mouth/load require identical mode coordinates and radius and equal static pressure",
            ));
        }
        Ok(())
    }

    fn prepare_loaded(
        source: ModalGalerkinSource,
        duct: PortHamiltonian,
        load: ModalAdmittanceLoad,
        duration: f64,
        budget: LinearSplitBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ModalTimeError> {
        let modes = source.modes();
        let (load_source, load, direct) = load.into_parts();
        let states = duct
            .state_dim()
            .checked_add(load.as_ref().map_or(0, PortHamiltonian::state_dim))
            .ok_or(ModalTimeError::Invalid(
                "complete loaded state count overflow",
            ))?;
        let inputs = if load.is_some() { modes.0 } else { 2 * modes.0 };
        let block = states
            .checked_add(inputs)
            .and_then(|n| n.checked_add(1))
            .and_then(|n| n.checked_mul(2))
            .ok_or(ModalTimeError::Invalid(
                "complete loaded exact-flow block overflow",
            ))?;
        if states > budget.max_modes || block > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM {
            return Err(ModalTimeError::Invalid(
                "complete duct and original load exceed authored/exact-flow capacity",
            ));
        }
        let mouth: Vec<_> = (modes.0..2 * modes.0).collect();
        let duct = duct.with_port_conductance(&mouth, &direct)?;
        let system = if let Some(load) = load {
            let pairs: Vec<_> = (0..modes.0).map(|mode| (modes.0 + mode, mode)).collect();
            fs_phs::interconnect(duct, load, &pairs)?
        } else {
            // These original mouth source inputs remain present and held at
            // zero. Their physical direct load already acts through R.
            duct
        };
        Self::prepare(
            source,
            system,
            Some(load_source),
            duration,
            budget,
            cancelled,
        )
    }

    fn prepare(
        source: ModalGalerkinSource,
        system: PortHamiltonian,
        load_source: Option<ModalLoadSource>,
        duration: f64,
        budget: LinearSplitBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ModalTimeError> {
        let states = system.state_dim();
        let inputs = system.port_dim();
        let count = if load_source.is_some() {
            source.modes().0
        } else {
            inputs
        };
        let selected: Vec<_> = (0..count)
            .map(|port| {
                (
                    port,
                    if port < source.modes().0 {
                        source.inlet_reference_impedance()
                    } else {
                        source.mouth_reference_impedance()
                    },
                )
            })
            .collect();
        let flow = PreparedMultiportScatteringZoh::new_owned(system, &selected, duration, budget)?;
        if cancelled() {
            return Err(fs_duct::DuctError::Cancelled {
                completed_steps: source.elements(),
            }
            .into());
        }
        Ok(Self {
            source,
            load_source,
            flow,
            state: vec![0.0; states],
            next: vec![0.0; states],
            mean_pressures: vec![0.0; inputs],
            reflected_means: vec![0.0; count],
            ports: vec![ScatteringPortRecord::default(); count],
            held_inputs: vec![0.0; inputs],
        })
    }

    /// Original continuous source and the actual full spatial operators.
    pub const fn source(&self) -> &ModalGalerkinSource {
        &self.source
    }
    /// Original attached load's physical mouth and gas, if present.
    pub const fn load_source(&self) -> Option<&ModalLoadSource> {
        self.load_source.as_ref()
    }
    /// Complete interconnected physical system used by the exact flow.
    /// Frequency-domain consumers can use the same source through fs-phs's
    /// harmonic owner without reconstructing or refitting the load.
    pub fn system(&self) -> &PortHamiltonian {
        self.flow.system()
    }

    /// Prepare the CONTINUOUS harmonic wave response of this complete original
    /// duct/load system. All inlet modes and every unbound mouth mode use the
    /// SAME reference feedback as the exact runtime, closed before factorization.
    /// No time-step transfer, mean-wave substitution or bare impedance inverse.
    ///
    /// Returned source slots are incident pressure waves at `wave_ports()` and
    /// ordinary volume-flow inputs elsewhere. With a static mouth conductance,
    /// its retained mouth input slots must be zero to reproduce an unforced load.
    /// Original duct/load state, physical heat and boundary wave powers stay
    /// separate. This query is independent of the accepted physical history and
    /// runtime clock; it does not establish a spatial or material accuracy band.
    pub fn prepare_harmonic_scattering(
        &self,
        harmonic_budget: fs_phs::harmonic::HarmonicBudget,
        scattering_budget: fs_phs::harmonic::scattering::ScatteringBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<fs_phs::harmonic::scattering::ScatteringModel, fs_phs::harmonic::HarmonicError>
    {
        let model = fs_phs::harmonic::HarmonicModel::from_ode(
            self.system(),
            harmonic_budget,
            &mut cancelled,
        )?;
        fs_phs::harmonic::scattering::ScatteringModel::new(
            model,
            self.flow.ports(),
            scattering_budget,
            cancelled,
        )
    }
    /// Accepted complete state: duct [potential; M*pressure; physical memory],
    /// followed by every original load coordinate when attached.
    pub fn state(&self) -> &[f64] {
        &self.state
    }
    /// Complete original duct coordinates, including its physical loss memory.
    pub fn duct_state(&self) -> &[f64] {
        &self.state[..self.source.state_dimension()]
    }
    /// Complete accepted original load coordinates; empty for an unbound mouth.
    pub fn load_state(&self) -> &[f64] {
        &self.state[self.source.state_dimension()..]
    }
    /// Actual exact-flow interval [s].
    pub fn duration(&self) -> f64 {
        self.flow.duration()
    }
    /// Total stored physical energy [J], with original duct and load storage.
    pub fn stored_energy_j(&self) -> f64 {
        self.flow.system().hamiltonian(&self.state)
    }
    /// Accepted interval-mean reflected pressure waves [Pa]. External order is
    /// inlet modes, then mouth modes only when that face remains a wave port.
    pub fn incoming(&self) -> &[f64] {
        &self.reflected_means
    }
    /// Accepted interval-mean pressure at every remaining original system port
    /// [Pa]. A purely conductive load retains the unforced mouth source slots
    /// after the inlet slots; a dynamic interconnection consumes those slots.
    pub fn mean_pressures(&self) -> &[f64] {
        &self.mean_pressures
    }
    /// Number of external incident/reflected wave channels.
    pub fn port_count(&self) -> usize {
        self.ports.len()
    }
    /// Independently integrated incident/outgoing energy and boundary work.
    /// Zero entries before the first accepted push represent no interval.
    pub fn port_accounts(&self) -> &[ScatteringPortRecord] {
        &self.ports
    }

    /// Same-interval affine relation b_bar=free+direct*a, using accepted
    /// physical history. The caller supplies P and P*P output buffers, where
    /// P is N for a loaded/ideal mouth or 2N for an external port mouth.
    /// There is no extra sample of delay or history commit during this query.
    ///
    /// # Errors
    /// Shared exact owner's shape/nonfinite/refusal; buffers stay unchanged.
    pub fn reflection_relation_into(
        &self,
        free: &mut [f64],
        direct: &mut [f64],
    ) -> Result<(), ModalTimeError> {
        // Selected entries are ignored by the relation. Any remaining mouth
        // forcing inputs are zero; their actual direct load already acts in R.
        self.flow
            .average_relation_into(&self.state, &self.held_inputs, free, direct)?;
        let count = self.ports.len();
        for mode in 0..count {
            direct[mode * count + mode] -= 1.0;
        }
        Ok(())
    }

    /// Consume one held incident pressure wave per ORIGINAL external port.
    /// Commit physical history only after the complete coupled candidate
    /// and every independent energy account succeed. No warm allocation.
    ///
    /// # Errors
    /// Shared shape, finite-arithmetic, exact-flow or physical-work refusal.
    /// On any refusal accepted state, pressure, reflection and accounts remain
    /// unchanged; prepared numerical scratch is not physical history.
    pub fn push(
        &mut self,
        incident: &[f64],
    ) -> Result<MultiportScatteringStepRecord, ModalTimeError> {
        if incident.len() != self.ports.len() || incident.iter().any(|v| !v.is_finite()) {
            return Err(ModalTimeError::Invalid(
                "one finite incident wave is required per external modal port",
            ));
        }
        self.held_inputs[..incident.len()].copy_from_slice(incident);
        let record = self.flow.step_into(
            &self.state,
            &self.held_inputs,
            &mut self.next,
            &mut self.mean_pressures,
            &mut self.ports,
        )?;
        core::mem::swap(&mut self.state, &mut self.next);
        for (wave, port) in self.reflected_means.iter_mut().zip(&self.ports) {
            *wave = port.outgoing_wave;
        }
        Ok(record)
    }
}
