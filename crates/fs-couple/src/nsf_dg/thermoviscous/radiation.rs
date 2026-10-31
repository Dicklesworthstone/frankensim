//! Reciprocal resolved mouths on the complete FE and finite-duct assembly.
//! Geometry, exterior memory, nonlinear solve and conservative NSF predictors
//! remain owned by their original implementations. Only the port basis changes.
use super::*;
use crate::nsf_dg::{NsfDgRadiationSolve, NsfDgReservoirRadiation};
use crate::nsf_dg::reservoir::radiation_surface;
use fs_lbm::nsf3::MAX_NSF_RESERVOIR_APERTURES;

/// Physical mouth bindings only. The complete assembly owns one nonlinear
/// policy and accepted warm start, including when an upstream field is added.
pub(super) struct ReservoirRadiations {
    pub(super) inputs: Vec<(usize, NsfDgReservoirRadiation)>,
    pub(super) temperatures_k: Vec<f64>,
}

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
    /// The explicit nonlinear policy also governs an optional upstream field;
    /// both cold attachments must agree with the same supplied policy.
    pub fn with_reservoir_radiations(
        mut self, inputs: Vec<(usize, NsfDgReservoirRadiation)>,
        solve: NsfDgRadiationSolve, gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        poll(gate)?;
        self.admit_joint_solve(solve)?;
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
        self.reservoir_radiations = Some(ReservoirRadiations { inputs, temperatures_k });
        self.bind_joint_solve(solve);
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
        self.reservoir_radiations.as_ref().and_then(|_| self.joint_solve())
    }
}
