//! Reciprocal resolved mouths on the complete FE and finite-duct assembly.
//! Geometry, exterior memory, nonlinear solve and conservative NSF predictors
//! remain owned by their original implementations. Only the port basis changes.
use super::*;
use crate::nsf_dg::{NsfDgRadiationSolve, NsfDgReservoirRadiation};
use crate::nsf_dg::reservoir::{
    exchange, radiation_surface, multiple::solve_joint_interval,
};
use fs_lbm::nsf3::MAX_NSF_RESERVOIR_APERTURES;

impl<'a, 'c, L: Lattice3<27>> NsfThermoviscousDg<'a, 'c, L> {
    /// Attach independent resolved escape apertures to the primary finite
    /// duct's original computed self/mutual exterior. That one exterior now
    /// feeds pressure back through both the complete FE region and each mouth.
    /// Additional finite ducts retain their own authored terminations.
    ///
    /// The authored order fixes the exterior coordinates. Each rectangle has
    /// one uniform acoustic-pressure mode and an explicit admissible ghost
    /// pressure bracket. This preserves the existing zero-reference-flow
    /// projection; it does not transport mean gas or enthalpy into the field.
    pub fn with_reservoir_radiations(
        mut self, inputs: Vec<(usize, NsfDgReservoirRadiation)>,
        solve: NsfDgRadiationSolve, gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        poll(gate)?;
        solve.validate()?;
        if !(1..=MAX_NSF_RESERVOIR_APERTURES).contains(&inputs.len())
            || self.steps != 0 || self.reservoir_radiations.is_some()
            || inputs.iter().enumerate().any(|(i, (aperture, _))|
                self.finite_reservoir_aperture() == Some(*aperture)
                    || inputs[..i].iter().any(|(prior, _)| prior == aperture))
        {
            return Err(NsfDgError::Input(
                "complete-field radiation needs distinct cold passive apertures and one exclusive exterior owner"));
        }
        let mut surfaces = Vec::with_capacity(inputs.len());
        let mut temperatures_k = Vec::with_capacity(inputs.len());
        let mut rate = None;
        for &(aperture, input) in &inputs {
            let (surface, temperature, source_rate) = radiation_surface(
                &self.zone, self.reference, self.load.system().load(), aperture, input, gate,
            )?;
            if rate.is_some_and(|previous| previous != source_rate) {
                return Err(NsfDgError::Input("complete-field mouths changed the original exterior bake clock"));
            }
            rate = Some(source_rate);
            surfaces.push(surface);
            temperatures_k.push(temperature);
        }
        self.load = self.load.with_independent_mouths(
            &surfaces, rate.expect("nonempty admitted mouths"), gate,
        ).map_err(NsfDgError::Thermoviscous)?;
        poll(gate)?;
        self.reservoir_radiations = Some(ReservoirRadiations {
            inputs, temperatures_k, solve, previous_coordinates: Vec::new(),
        });
        Ok(self)
    }

    /// Select one passive mouth explicitly, independently of its breath supply.
    pub fn with_reservoir_radiation_at(
        self, aperture: usize, input: NsfDgReservoirRadiation,
        solve: NsfDgRadiationSolve, gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        self.with_reservoir_radiations(vec![(aperture, input)], solve, gate)
    }

    /// Original passive declarations, in shared exterior coordinate order.
    pub fn reservoir_radiations(&self) -> &[(usize, NsfDgReservoirRadiation)] {
        self.reservoir_radiations.as_ref().map_or(&[], |owner| owner.inputs.as_slice())
    }

    pub fn reservoir_radiation_at(&self, aperture: usize) -> Option<&NsfDgReservoirRadiation> {
        self.reservoir_radiations().iter().find(|(index, _)| *index == aperture)
            .map(|(_, input)| input)
    }

    pub fn radiation_solve(&self) -> Option<&NsfDgRadiationSolve> {
        self.reservoir_radiations.as_ref().map(|owner| &owner.solve)
    }

    pub(super) fn step_with_radiations(
        &mut self, drive: NsfDgDrive, external: NsfBoundaryDrive,
        breath: Option<(usize, NsfOpenState)>, gate: &CancelGate,
    ) -> Result<NsfThermoviscousDgFrame, NsfDgError> {
        let owner = self.reservoir_radiations.as_ref()
            .ok_or(NsfDgError::Input("complete-field mouths are not attached"))?;
        let inputs = owner.inputs.clone();
        let temperatures = owner.temperatures_k.clone();
        let solve = owner.solve;
        let previous_coordinates = owner.previous_coordinates.clone();
        let relation = self.load.independent_mouths_response(gate)
            .map_err(NsfDgError::Thermoviscous)?;
        let width = inputs.len() + 1;
        let limits = self.limits;
        let side = usize::from(self.high);
        let reference = self.reference;
        let step = self.steps + 1;
        let dt = self.load.duration_s();
        let load = &mut self.load;
        let finite = &mut self.finite_reservoir;
        let (frame, normalized) = solve_joint_interval(
            &mut self.zone, &inputs, &temperatures, external, breath, side,
            self.patch_background, reference.pressure, limits, solve,
            &previous_coordinates, &relation, None, gate,
            |fluid, physical, evaluations| {
                let field = coupled_exchange(&fluid, side)?;
                let mouths = inputs.iter().map(|(aperture, _)| exchange(&fluid, *aperture))
                    .collect::<Result<Vec<_>, _>>()?;
                let mouth_flows: Vec<_> = mouths.iter()
                    .map(|mouth| mouth.normal_volume_flow_m3_s).collect();
                let mut flows = Vec::with_capacity(width);
                flows.push(field.normal_volume_flow_m3_s);
                flows.extend_from_slice(&mouth_flows);
                let scalar_relation = |row: usize| {
                    let mut history = relation.history_pressure_pa[row];
                    for (column, flow) in flows.iter().enumerate() {
                        if column != row {
                            history += relation.impedance_pa_s_m3[row * width + column] * flow;
                        }
                    }
                    DgPortResponse { history_pressure_pa: history,
                        impedance_pa_s_m3: relation.impedance_pa_s_m3[row * width + row] }
                };
                let acoustic = load.preview_step_with_mouth_flows(
                    field.normal_volume_flow_m3_s, &mouth_flows, gate,
                ).map_err(NsfDgError::Thermoviscous)?;
                let radiation = acoustic.independent_mouths_frame()
                    .ok_or(NsfDgError::Input("complete field lost its original mouth receipts"))?;
                if radiation.mouths.len() != mouths.len() {
                    return Err(NsfDgError::Input("complete-field exterior receipt count changed"));
                }
                let mut reservoir_radiations = Vec::with_capacity(inputs.len());
                for (index, (&mouth, receipt)) in mouths.iter().zip(&radiation.mouths).enumerate() {
                    let audit = upstream::audit(mouth, receipt.mean_pressure_pa, receipt.load_work_j,
                        0.0, scalar_relation(index + 1), reference, dt, limits)?;
                    let resolved_exchange = *fluid.reservoir_exchange_at(inputs[index].0)
                        .ok_or(NsfDgError::Input("complete-field mouth lost its conservative NSF exchange"))?;
                    reservoir_radiations.push(NsfDgReservoirFrame {
                        resolved_exchange, ghost_pressure_offset_pa: physical[index + 1],
                        mean_acoustic_pressure_pa: mouth.mean_acoustic_pressure_pa,
                        outward_flow_m3_s: mouth.normal_volume_flow_m3_s,
                        radiation_pressure_pa: receipt.mean_pressure_pa,
                        radiation_work_j: receipt.load_work_j,
                        pressure_residual_pa: audit.pressure, projection_covariance_j: audit.covariance,
                        acoustic_work_residual_j: audit.work, mass_projection_difference_kg: audit.mass,
                    });
                }
                accept_trial(fluid, acoustic, finite, side, scalar_relation(0), reference, limits,
                    dt, physical[0], drive, evaluations, step, reservoir_radiations, gate)
            },
        )?;
        self.reservoir_radiations.as_mut().expect("original complete-field mouth owner")
            .previous_coordinates = normalized;
        self.steps = step;
        Ok(frame)
    }
}
