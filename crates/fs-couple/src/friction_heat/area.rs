//! Actual endpoint contact area on the original finite caloric receivers.

use super::*;
use fs_thermal::enthalpy::ContactPatchUpdate;
use fs_thermal::enthalpy::backward_euler::ContactPatchEndpoint;

impl FrictionHeatRoute {
    fn patch_at_area(&self, area_m2: f64) -> Result<ContactPatch, DepositionError> {
        if !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(invalid(
                "a changing friction interface needs positive finite contact area",
            ));
        }
        for side in 0..2 {
            let face = finite(self.geometry[side].volume_m3 / (2.0 * self.distances_m[side]))?;
            if area_m2 - face > 64.0 * f64::EPSILON * area_m2.max(face) {
                return Err(invalid(
                    "changing friction contact exceeds its original prismatic face",
                ));
            }
        }
        let mut patch = self.patch;
        patch.link.area_m2 = area_m2;
        // The direct path is absent while a massless common-temperature
        // interface supplies both half-cell fluxes. This is not zero material k.
        patch.link.conductance_w_k = 0.0;
        patch.covered_area_m2 = [area_m2; 2];
        Ok(patch)
    }

    fn at_area_with_resistance(
        &self,
        area_m2: f64,
        conductivity_w_m_k: [f64; 2],
        area_resistance_m2_k_w: [f64; 2],
    ) -> Result<Self, DepositionError> {
        let resistance = Self::total_resistance(area_resistance_m2_k_w)?;
        let mut next = self.clone();
        next.patch = self.patch_at_area(area_m2)?;
        next.area_resistance_m2_k_w = area_resistance_m2_k_w;
        next.patch.link.conductance_w_k = series_contact_conductance_w_k(
            area_m2,
            self.distances_m,
            conductivity_w_m_k,
            resistance,
        )?;
        next.conductances(conductivity_w_m_k)?;
        Ok(next)
    }
}

impl<'s, 'm> FrictionHeatStep<'s, 'm> {
    /// Solve heat on an actual endpoint area within the original uniform faces.
    ///
    /// `prepare` stages original mechanical history and returns its positive
    /// contact area. The native thermal owner then updates actual coverage.
    /// `resolve` sees that covered endpoint, with the direct exchange suppressed,
    /// and supplies original material k, exported friction work and other heat.
    /// No rescaled conductivity or compensating ambient source is introduced.
    ///
    /// Accepted contact coverage AND the source-derived series conductance are
    /// retained in the returned thermal state and successor route. Cell geometry,
    /// inventory and caloric objects remain unchanged. Adopt both together.
    ///
    /// This is an endpoint, uniform-face approximation, not swept contact
    /// residence or lateral heat spreading. Discontinuous finite-cell contact
    /// changes may prevent nonlinear convergence and require time/grid refinement.
    /// Complete detachment refuses rather than creating an arbitrary temperature
    /// for a vanished massless interface.
    pub fn evaluate_with_area<A, P, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        interface_temperature_k: f64,
        prepare: impl FnOnce(&DepositionSimulation<'m>, f64) -> Result<(f64, A), E>,
        resolve: impl FnOnce(&DepositionSimulation<'m>, f64, A) -> Result<FrictionHeatInput<P>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<FrictionHeatAreaTrial<'m, P>, E> {
        self.evaluate_area_option(
            energy_increments_j,
            interface_temperature_k,
            prepare,
            |state, temperature, staged| {
                resolve(state, temperature, staged)
                    .map(|input| (input, self.route.area_resistance_m2_k_w))
            },
            false,
            checkpoint,
        )
    }

    /// The same changing-area endpoint with exact native caloric entropy.
    /// Actual uncovered area supplies ambient heat and entropy exactly once.
    pub fn evaluate_with_area_and_entropy<A, P, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        interface_temperature_k: f64,
        prepare: impl FnOnce(&DepositionSimulation<'m>, f64) -> Result<(f64, A), E>,
        resolve: impl FnOnce(&DepositionSimulation<'m>, f64, A) -> Result<FrictionHeatInput<P>, E>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<FrictionHeatAreaTrial<'m, P>, E> {
        self.evaluate_area_option(
            energy_increments_j,
            interface_temperature_k,
            prepare,
            |state, temperature, staged| {
                resolve(state, temperature, staged)
                    .map(|input| (input, self.route.area_resistance_m2_k_w))
            },
            true,
            checkpoint,
        )
    }

    #[allow(clippy::too_many_arguments)]
    /// Evaluate actual area and two source-resolved contact resistances together.
    /// The callback returns resistance from the massless heat-generation plane
    /// to each original body [m² K/W]. The same native cell/interface residuals
    /// and acceptance gates own all heat. Adopt the returned route and thermal
    /// state together; finite resistance is retained in the accepted direct link.
    pub fn evaluate_with_area_and_resistance<A, P, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        interface_temperature_k: f64,
        prepare: impl FnOnce(&DepositionSimulation<'m>, f64) -> Result<(f64, A), E>,
        resolve: impl FnOnce(
            &DepositionSimulation<'m>,
            f64,
            A,
        ) -> Result<(FrictionHeatInput<P>, [f64; 2]), E>,
        with_entropy: bool,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<FrictionHeatAreaTrial<'m, P>, E> {
        self.evaluate_area_option(
            energy_increments_j,
            interface_temperature_k,
            prepare,
            resolve,
            with_entropy,
            checkpoint,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate_area_option<A, P, E: From<DepositionError>>(
        &self,
        energy_increments_j: &[f64],
        interface_temperature_k: f64,
        prepare: impl FnOnce(&DepositionSimulation<'m>, f64) -> Result<(f64, A), E>,
        resolve: impl FnOnce(
            &DepositionSimulation<'m>,
            f64,
            A,
        ) -> Result<(FrictionHeatInput<P>, [f64; 2]), E>,
        with_entropy: bool,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<FrictionHeatAreaTrial<'m, P>, E> {
        if !interface_temperature_k.is_finite() || interface_temperature_k <= 0.0 {
            return Err(
                invalid("friction interface temperature must be finite positive kelvin").into(),
            );
        }
        let checkpoint = RefCell::new(checkpoint);
        let check = || (*checkpoint.borrow_mut())();
        let mut successor_route = None;
        let coefficients = |endpoint: &mut ContactPatchEndpoint<'_, 'm>|
            -> Result<StepCoefficients<TrialPayload<P>>, E> {
            let (area_m2, staged) = prepare(endpoint.thermal(), interface_temperature_k)?;
            let provisional = self.route.patch_at_area(area_m2)?;
            endpoint.replace_contact_patches(&[ContactPatchUpdate {
                id: self.route.link_id(),
                expected: Some(self.route.patch),
                replacement: Some(provisional),
            }])?;
            let (input, resistance) = resolve(endpoint.thermal(), interface_temperature_k, staged)?;
            let route = self.route.at_area_with_resistance(area_m2, input.conductivity_w_m_k, resistance)?;
            endpoint.replace_contact_patches(&[ContactPatchUpdate {
                id: self.route.link_id(),
                expected: Some(provisional),
                replacement: Some(route.patch),
            }])?;
            let coefficients = self.coefficients_for_input(
                &route, endpoint.thermal(), interface_temperature_k,
                input, with_entropy, &check,
            )?;
            successor_route = Some(route);
            Ok(coefficients)
        };
        let native = if with_entropy {
            self.thermal.evaluate_with_contact_patches_and_entropy(
                energy_increments_j,
                coefficients,
                &check,
            )?
        } else {
            self.thermal
                .evaluate_with_contact_patches(energy_increments_j, coefficients, &check)?
        };
        let route = successor_route.ok_or(DepositionError::Arithmetic)?;
        let trial = self.finish_trial(native, interface_temperature_k)?;
        Ok(FrictionHeatAreaTrial { trial, route })
    }
}

/// One native heat trial and the same endpoint's original-material contact route.
pub struct FrictionHeatAreaTrial<'m, P> {
    trial: FrictionHeatTrial<'m, P>,
    route: FrictionHeatRoute,
}

impl<'m, P> FrictionHeatAreaTrial<'m, P> {
    pub fn residual_j(&self) -> &[f64] {
        self.trial.residual_j()
    }
    pub fn payload(&self) -> &P {
        self.trial.payload()
    }
    pub const fn report(&self) -> &FrictionHeatReport {
        self.trial.report()
    }
    pub const fn route(&self) -> &FrictionHeatRoute {
        &self.route
    }

    pub fn accept(
        self,
        tolerance: FrictionHeatTolerances,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (
            DepositionSimulation<'m>,
            P,
            FrictionHeatReport,
            FrictionHeatRoute,
        ),
        DepositionError,
    > {
        let (thermal, payload, report) = self.trial.accept(tolerance, checkpoint)?;
        Ok((thermal, payload, report, self.route))
    }

    pub(super) fn into_trial_and_route(self) -> (FrictionHeatTrial<'m, P>, FrictionHeatRoute) {
        (self.trial, self.route)
    }
}
