//! Explicit magnetic temperature laws on the retained tetrahedral A–phi solve.

use super::*;

impl SpatialMagneticThermalNetwork<'_> {
    /// Evolve sourced nonlinear magnetism, induced currents and finite heat.
    ///
    /// Every thermal region declares its original linear mu(T) or separable
    /// H(B,T)=s(T)H_ref(B) law and admitted B envelope. Source refresh and the
    /// magnetic entropy calculation are shared with planar induction. The
    /// original tetrahedral field retains A, impressed-current history, gauge,
    /// conductor topology and mechanical-work convention across every trial.
    ///
    /// Parameter work P=psi_new(B_old)-psi_old(B_old) and delta(T*S) jointly
    /// supply reversible heat Q_rev=-P-delta(T*S). Local sigma*E_conductor² is
    /// deposited once; backward-Euler convexity loss remains numerical export.
    /// The physical balance uses U=psi+T*S, including fixed exterior regions.
    ///
    /// The reference nonlinear curve is temperature independent and its
    /// positive scale affine over the declared interval. Density, geometry and
    /// solid occupancy stay fixed. Prescribed velocity exchanges Lorentz work
    /// with an external mechanical receiver; it does not advance mechanics.
    /// General H(B,T), recoil entropy, hysteresis and phase change require
    /// additional constitutive owners. Source or convergence failure publishes
    /// neither accepted history. Coupled time accuracy remains first order.
    #[allow(clippy::too_many_arguments)]
    pub fn step_nonlinear_thermal_candidate<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialNonlinearInduction<'m>,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        bindings: &[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>],
        controls: CircuitEnthalpySolve,
        field_controls: SpatialNonlinearInductionControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearInductionThermalCandidate<'t, 'm>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_constitutive_resolved_candidate(
            thermal,
            target_s,
            field,
            endpoint_current_a_m2,
            velocity_m_s,
            Some(bindings),
            controls,
            field_controls,
            initial_energy_increments_j,
            |trial| {
                prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            |sample| sample,
            &cancel,
        )
    }
}
