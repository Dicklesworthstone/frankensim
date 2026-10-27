//! Native DC current on an explicitly bounded, fully occupied liquid core.
//!
//! The existing voxel current owner supplies half-cell/contact I^2 R. The
//! existing thermal timeline deposits it once over each accepted interval;
//! its actual T and original liquid shear state refresh the complete source
//! bundle before the next solve. No second current or thermal solver exists.

use super::*;
use crate::conjugate::{ElectricSetup, ElectricSolution, ElectricalContact, VoxelDomain};
use crate::d3q19::freesurface3::Cell3;
use fs_material::conductor::ResolvedConductorState;

/// Imposed-field motional current and native Lorentz forcing.
pub mod magnetic;
pub use magnetic::{FreeSurfaceMagneticInput, FreeSurfaceMagneticState, FreeSurfaceMagneticMotion};

/// A fixed Cartesian electrical subdomain inside the moving liquid.
///
/// Every selected cell, including ideal electrode cells, must remain a native
/// `Fluid` cell. The native current solver uses complete dx^3 cells and dx^2
/// faces; a partial VOF cell needs a cut-cell electrical closure and refuses.
/// Origin is in global liquid coordinates; electrodes and contact cell indices
/// are LOCAL to this core. Its exterior faces explicitly impose zero normal
/// current, as insulating or symmetry boundaries. They do not infer electrical
/// insulation for the surrounding liquid or create mechanical walls.
///
/// Electrode regions keep the original zero-bulk-resistance constraint. Their
/// liquid thermal mass remains real; explicitly partitioned contact heat can
/// enter it. No finite electrode body, electrochemistry or contact-area law is
/// inferred. Every liquid bundle must contain exactly one scalar electrical
/// conductivity or resistivity source. Changing that declared source basis is
/// a data change; missing or ambiguous transport never becomes insulation.
#[derive(Clone, Debug)]
pub struct FreeSurfaceElectricalInput {
    pub origin: [usize; 3],
    pub dimensions: [usize; 3],
    pub setup: ElectricSetup,
    pub contacts: Vec<ElectricalContact>,
    /// Checked before allocating the native core or querying its material cells.
    pub maximum_cells: usize,
    /// All applied and endpoint current solves in ONE native flow interval.
    pub maximum_solves: usize,
    /// Absolute electrode-power minus Joule and sampled mechanical power allowance [W].
    /// The conduction energy tolerance additionally bounds its integrated error.
    pub power_tolerance_w: f64,
}

/// One native frozen current solve and the exact material epoch it consumed.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceElectricalState {
    /// Time within this thermal/flow interval, beginning at zero [s].
    pub time_s: f64,
    /// Core-local voxel order, each entry naming the actual liquid cell.
    pub coordinates: Vec<[usize; 3]>,
    pub specific_enthalpy_j_kg: Vec<f64>,
    pub temperature_k: Vec<f64>,
    /// Complete original T/shear/thermal/electrical source bundles and receipts.
    pub materials: Vec<ResolvedConductorState>,
    pub solution: ElectricSolution,
    /// Electrode power minus Joule power minus sampled Lorentz power [W].
    pub power_defect_w: f64,
    /// Present for an imposed-B solve; Joule heat remains conductor-frame I^2 R.
    pub magnetic: Option<FreeSurfaceMagneticState>,
}

/// Electrical supply and heat on exactly the accepted thermal intervals.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceElectricalExchange {
    /// Sum dt * native I^2 R; included ONCE in liquid external heat [J].
    pub joule_heat_j: f64,
    /// Sum dt * sum_e(V_e I_e), independently measured by the current owner [J].
    pub electrode_energy_j: f64,
    /// Signed sum dt * F dot sampled initial velocity; negative for braking.
    /// This is not an assertion about the subsequent native kinetic-energy change.
    pub sampled_lorentz_work_j: f64,
    /// Core-local integrated native Lorentz forces [N s], empty without B.
    pub lorentz_impulse_n_s: Vec<[f64; 3]>,
    /// Electrode energy minus Joule heat minus sampled Lorentz work [J].
    pub energy_defect_j: f64,
    pub substeps: usize,
    pub current_solves: usize,
    /// Last epoch actually used to deposit heat, not an endpoint diagnostic.
    pub last_applied: FreeSurfaceElectricalState,
    /// Current at the accepted thermal endpoint BEFORE subsequent VOF advection.
    /// This diagnostic contributes no extra heat or electrode energy.
    pub thermal_endpoint: FreeSurfaceElectricalState,
}

#[derive(Clone, Debug)]
pub(crate) struct ElectricalCore {
    input: FreeSurfaceElectricalInput,
    domain: VoxelDomain,
    coordinates: Vec<[usize; 3]>,
    magnetic: Option<FreeSurfaceMagneticInput>,
}

impl ElectricalCore {
    fn new(input: FreeSurfaceElectricalInput, dx: f64, limit: usize) -> Result<Self, Error> {
        let count = input.dimensions.iter().try_fold(1usize, |n, &width| n.checked_mul(width))
            .ok_or(Error::Invalid("electrical core cell count overflows"))?;
        if count == 0 || count > input.maximum_cells || input.maximum_cells > limit
            || input.maximum_solves == 0 || input.maximum_solves > 1_000_000
            || input.contacts.len() > count.saturating_mul(3)
            || !input.power_tolerance_w.is_finite() || input.power_tolerance_w < 0.0
            || (0..3).any(|axis| input.origin[axis].checked_add(input.dimensions[axis]).is_none())
        {
            return Err(Error::Invalid("invalid electrical core geometry, work or power bounds"));
        }
        let [nx, ny, nz] = input.dimensions;
        let domain = VoxelDomain::new(nx, ny, nz, dx)?;
        let coordinates = (0..count).map(|cell| {
            let p = domain.coords(cell);
            core::array::from_fn(|axis| input.origin[axis] + p[axis])
        }).collect();
        Ok(Self { input, domain, coordinates, magnetic: None })
    }

    pub(crate) fn admit(
        &self,
        lattice: &FreeSurface3,
        materials: &BTreeMap<MassLocation3, ResolvedLiquidState>,
        gate: &CancelGate,
    ) -> Result<(), Error> {
        for &p in &self.coordinates {
            poll(gate)?;
            let sample = lattice.sample_at(p[0], p[1], p[2])
                .ok_or(Error::Invalid("electrical core has no native liquid cell"))?;
            if sample.cell != Cell3::Fluid || sample.fill_fraction != 1.0 {
                return Err(Error::Invalid("electrical core requires fully occupied native Fluid cells"));
            }
            let location = FreeSurface3::transport_location(p)
                .ok_or(Error::Invalid("electrical core coordinates cannot be represented"))?;
            let material = materials.get(&location)
                .ok_or(Error::Invalid("electrical core has no accepted liquid source bundle"))?;
            ResolvedConductorState::try_from_material(material.material().clone())
                .map_err(Error::Conductor)?;
        }
        Ok(())
    }
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Bind source-resolved DC current to every subsequent conduction timeline.
    /// Existing finite-wall, environmental, finite-gas, rheology and wetting
    /// entry points all consume this same binding automatically.
    ///
    /// Each native thermal substep solves current using its actual full liquid
    /// T/shear source bundle, deposits I^2 R once, and checks the electrical
    /// source again at its endpoint. This is explicit thermal feedback followed
    /// by the original thermal/flow split; refine the interval for nonlinear
    /// resistive feedback. Source support bounds are not stability guarantees.
    ///
    /// Pure flow and instantaneous-heat entry points refuse an attached core:
    /// use `step_with_conduction` or its coupled variants so the drive cannot
    /// silently disappear. A core that dewets, loses its source support or
    /// exhausts its budget aborts the whole liquid/wall/gas transaction.
    ///
    /// Without [`Self::with_imposed_magnetic_field`], this zero-B path supplies
    /// no motional EMF or Lorentz force. That opt-in reuses the original moving
    /// current and native forcing owners. Neither path adds induction, Maxwell
    /// traction, phase conversion or partial-cell current. At zero B,
    /// electrical work equals Joule heat to the native solve tolerance; no
    /// pressure work, capillary energy or numerical dissipation is thermalized.
    pub fn with_electrical_core(
        mut self,
        input: FreeSurfaceElectricalInput,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        poll(gate)?;
        let core = ElectricalCore::new(input, self.flow.scaling.cell_size_m(), self.controls.max_cells)?;
        core.admit(&self.flow.lattice, &self.states, gate)?;
        self.electrical = Some(core);
        Ok(self)
    }

    /// Retained geometry, drive, contacts and work bounds of this thermal binding.
    #[must_use]
    pub fn electrical_core(&self) -> Option<&FreeSurfaceElectricalInput> {
        self.electrical.as_ref().map(|core| &core.input)
    }
}

#[derive(Clone)]
pub(super) struct ElectricalTimeline {
    core: ElectricalCore,
    pub(super) cells: Vec<CellId>,
    pub(super) last: Option<FreeSurfaceElectricalState>,
    endpoint: Option<FreeSurfaceElectricalState>,
    solves: usize,
    substeps: usize,
    joule_heat_j: f64,
    electrode_energy_j: f64,
    velocity_m_s: Vec<[f64; 3]>,
    sampled_lorentz_work_j: f64,
    lorentz_impulse_n_s: Vec<[f64; 3]>,
}

impl ElectricalTimeline {
    pub(super) fn new(
        core: &ElectricalCore,
        cell_id: &impl Fn([usize; 3]) -> Result<CellId, Error>,
        flow: &MaterialFreeSurface3<'_>,
    ) -> Result<Self, Error> {
        let cells = core.coordinates.iter().map(|&p| cell_id(p)).collect::<Result<_, _>>()?;
        let velocity_m_s = core.velocities(&flow.lattice, &flow.scaling)?;
        let lorentz_impulse_n_s = if core.magnetic.is_some() { vec![[0.0; 3]; core.coordinates.len()] }
            else { Vec::new() };
        Ok(Self {
            core: core.clone(), cells, last: None, endpoint: None,
            solves: 0, substeps: 0, joule_heat_j: 0.0, electrode_energy_j: 0.0,
            velocity_m_s, sampled_lorentz_work_j: 0.0, lorentz_impulse_n_s,
        })
    }

    pub(super) fn accept(
        &mut self,
        start_s: f64,
        duration_s: f64,
        endpoint: FreeSurfaceElectricalState,
    ) -> Result<(), Error> {
        let applied = self.last.as_ref()
            .ok_or(Error::Invalid("electrical interval has no applied current epoch"))?;
        if applied.time_s.to_bits() != start_s.to_bits()
            || !duration_s.is_finite() || duration_s <= 0.0
        {
            return Err(Error::Invalid("electrical heat does not match its accepted thermal interval"));
        }
        self.joule_heat_j = finite(self.joule_heat_j + duration_s * applied.solution.total_joule_w)?;
        self.electrode_energy_j = finite(self.electrode_energy_j + duration_s * applied.solution.delivered_w)?;
        if let Some(magnetic) = &applied.magnetic {
            self.sampled_lorentz_work_j = finite(self.sampled_lorentz_work_j
                + duration_s * magnetic.mechanical_power_w)?;
            for (cell, (impulse, force)) in self.lorentz_impulse_n_s.iter_mut()
                .zip(&magnetic.lorentz_force_density_n_m3).enumerate()
            {
                for axis in 0..3 {
                    impulse[axis] = finite(impulse[axis]
                        + duration_s * self.core.domain.volume(cell) * force[axis])?;
                }
            }
        }
        self.substeps = self.substeps.checked_add(1)
            .ok_or(Error::Invalid("electrical substep count overflow"))?;
        self.endpoint = Some(endpoint);
        Ok(())
    }

    pub(super) fn report(&self, tolerance_j: f64) -> Result<FreeSurfaceElectricalExchange, Error> {
        let defect = finite(self.electrode_energy_j - self.joule_heat_j - self.sampled_lorentz_work_j)?;
        if defect.abs() > tolerance_j {
            return Err(Error::Invalid("electrical supply, Joule heat and sampled Lorentz work do not close"));
        }
        Ok(FreeSurfaceElectricalExchange {
            joule_heat_j: self.joule_heat_j,
            electrode_energy_j: self.electrode_energy_j,
            sampled_lorentz_work_j: self.sampled_lorentz_work_j,
            lorentz_impulse_n_s: self.lorentz_impulse_n_s.clone(),
            energy_defect_j: defect,
            substeps: self.substeps,
            current_solves: self.solves,
            last_applied: self.last.clone()
                .ok_or(Error::Invalid("electrical timeline has no applied current"))?,
            thermal_endpoint: self.endpoint.clone()
                .ok_or(Error::Invalid("electrical timeline has no endpoint current"))?,
        })
    }
}

impl Conduction<'_, '_> {
    pub(super) fn sample_current(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<FreeSurfaceElectricalState>, Error> {
        check()?;
        let Some(electrical) = self.electrical.as_mut() else { return Ok(None); };
        if electrical.solves >= electrical.core.input.maximum_solves {
            return Err(Error::Invalid("electrical current-solve budget exhausted"));
        }
        electrical.solves += 1;
        let electrical = self.electrical.as_ref().expect("bound electrical core");
        let mut materials = Vec::with_capacity(electrical.cells.len());
        let mut enthalpy = Vec::with_capacity(electrical.cells.len());
        let mut temperature = Vec::with_capacity(electrical.cells.len());
        for &id in &electrical.cells {
            check()?;
            let h = thermal.cells()[&id].specific_h_j_kg;
            let liquid = self.sample_liquid(id, h, check)?;
            materials.push(ResolvedConductorState::try_from_material(liquid.material().clone())
                .map_err(Error::Conductor)?);
            enthalpy.push(h);
            temperature.push(self.caloric.state_at_specific_enthalpy(h)?.temperature_k());
        }
        let core = &electrical.core;
        let (solution, magnetic) = electrical.solve_current(&temperature, &materials, self.gate)?;
        let work = magnetic.as_ref().map_or(0.0, |state| state.mechanical_power_w);
        let defect = finite(solution.delivered_w - solution.total_joule_w - work)?;
        if defect.abs() > core.input.power_tolerance_w
            || solution.joule_w.iter().any(|watts| !watts.is_finite() || *watts < 0.0)
        {
            return Err(Error::Invalid("native electrical power balance or Joule heat is inadmissible"));
        }
        check()?;
        Ok(Some(FreeSurfaceElectricalState {
            time_s: thermal.time_s(),
            coordinates: core.coordinates.clone(),
            specific_enthalpy_j_kg: enthalpy,
            temperature_k: temperature,
            materials,
            solution,
            power_defect_w: defect,
            magnetic,
        }))
    }
}

pub(super) fn finite(value: f64) -> Result<f64, Error> {
    if value.is_finite() { Ok(value) }
    else { Err(Error::Invalid("electrical thermal account is not finite")) }
}
