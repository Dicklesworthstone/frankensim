//! Private design-map calculus shared by dead-load and thermal compliance.
//! The original filter, SIMP chains and inventory operators retain their order.

use super::finite;
use crate::elasticity::second_order::{ComplianceError, ComplianceSolveOptions};
use crate::filter::{DensityFilter, heaviside, heaviside_derivative, heaviside_second_derivative};
use crate::pipeline::{LoadCase, SimpParams};
use crate::{DensityElasticity, SolveControl};

pub(super) struct DesignMap<'a> {
    filter: &'a DensityFilter,
    reference_cell_masses_kg: Option<&'a [f64]>,
    pub(super) projected: Vec<f64>,
    pub(super) volume_projected: Vec<f64>,
    first: Vec<f64>,
    second: Vec<f64>,
    // The mechanical realization's occupancy slope; robust inventory uses a
    // different projection. Retain the original arithmetic instead of dividing
    // a tiny SIMP slope to reconstruct it in frequency derivatives.
    stiffness_projection_first: Vec<f64>,
    projection_first: Vec<f64>,
    projection_second: Vec<f64>,
}

impl<'a> DesignMap<'a> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn prepare(
        filter: &'a DensityFilter,
        p: SimpParams,
        volume_eta: f64,
        elasticity: &'a DensityElasticity,
        rho: &[f64],
        loads: &[LoadCase<'_>],
        options: ComplianceSolveOptions,
        control: &mut SolveControl<'_>,
    ) -> Result<(Self, Vec<f64>), ComplianceError> {
        control.checkpoint("simp-hessian-admission")?;
        let n = elasticity.cells();
        if n == 0 || rho.len() != n || filter.cells() != n {
            return Err(ComplianceError::InvalidInput(
                "raw design, filter and elasticity shapes must agree",
            ));
        }
        if !(p.e_min.is_finite()
            && p.e_min > 0.0
            && p.e_min < 1.0
            && p.penal.is_finite()
            && p.penal >= 1.0
            && p.beta.is_finite()
            && p.beta >= 0.0
            && p.eta.is_finite()
            && (0.0..=1.0).contains(&p.eta)
            && volume_eta.is_finite()
            && (0.0..=1.0).contains(&volume_eta))
        {
            return Err(ComplianceError::InvalidInput(
                "invalid SIMP/projection parameters",
            ));
        }
        if !(options.tolerance.is_finite() && options.tolerance > 0.0 && options.tolerance < 1.0) {
            return Err(ComplianceError::InvalidInput(
                "relative tolerance must lie in (0, 1)",
            ));
        }
        for (e, &r) in rho.iter().enumerate() {
            if e % 256 == 0 {
                control.checkpoint("simp-hessian-admission")?;
            }
            if !r.is_finite() || !(0.0..=1.0).contains(&r) {
                return Err(ComplianceError::InvalidInput(
                    "raw densities must lie in [0, 1]",
                ));
            }
        }
        if loads.is_empty() {
            return Err(ComplianceError::InvalidInput(
                "at least one structural load case required",
            ));
        }
        for load in loads {
            if load.force.len() != elasticity.n()
                || !(load.weight.is_finite() && load.weight >= 0.0)
            {
                return Err(ComplianceError::InvalidInput(
                    "invalid structural load shape or weight",
                ));
            }
            for (d, &f) in load.force.iter().enumerate() {
                if d % 256 == 0 {
                    control.checkpoint("simp-hessian-admission")?;
                }
                finite(f, "structural force")?;
            }
        }
        let filtered = filter.try_apply(rho, control)?;
        let mut projected = Vec::with_capacity(n);
        let mut volume_projected = Vec::with_capacity(n);
        let mut moduli = Vec::with_capacity(n);
        let mut first = Vec::with_capacity(n);
        let mut second = Vec::with_capacity(n);
        let mut stiffness_projection_first = Vec::with_capacity(n);
        let mut projection_first = Vec::with_capacity(n);
        let mut projection_second = Vec::with_capacity(n);
        for (e, &rt) in filtered.iter().enumerate() {
            if e % 256 == 0 {
                control.checkpoint("simp-curvature")?;
            }
            let rb = finite(heaviside(rt, p.beta, p.eta), "stiffness projection")?;
            if !(rb > 1e-12 && rb < 1.0) {
                return Err(ComplianceError::InvalidInput(
                    "SIMP Hessian requires projected density strictly in (1e-12, 1)",
                ));
            }
            projected.push(rb);
            // Same native projection and SIMP arithmetic as try_forward. Since
            // rb is strictly inside the smooth interval, clamp/floor are idle.
            moduli.push(finite(
                p.e_min + (1.0 - p.e_min) * fs_math::det::pow(rb, p.penal),
                "SIMP modulus",
            )?);
            let dh = finite(heaviside_derivative(rt, p.beta, p.eta), "projection slope")?;
            stiffness_projection_first.push(dh);
            let ddh = finite(
                heaviside_second_derivative(rt, p.beta, p.eta),
                "projection curvature",
            )?;
            let dm = (1.0 - p.e_min) * p.penal * fs_math::det::pow(rb, p.penal - 1.0);
            let ddm = if p.penal == 1.0 {
                0.0
            } else {
                (1.0 - p.e_min) * p.penal * (p.penal - 1.0) * fs_math::det::pow(rb, p.penal - 2.0)
            };
            first.push(finite(dm * dh, "SIMP first chain")?);
            second.push(finite(ddm * dh * dh + dm * ddh, "SIMP second chain")?);
            let (vb, dv, ddv) = if volume_eta == p.eta {
                (rb, dh, ddh)
            } else {
                (
                    finite(heaviside(rt, p.beta, volume_eta), "volume projection")?,
                    finite(
                        heaviside_derivative(rt, p.beta, volume_eta),
                        "volume projection slope",
                    )?,
                    finite(
                        heaviside_second_derivative(rt, p.beta, volume_eta),
                        "volume projection curvature",
                    )?,
                )
            };
            if !(0.0..=1.0).contains(&vb) {
                return Err(ComplianceError::InvalidInput(
                    "volume projection must lie in [0,1]",
                ));
            }
            volume_projected.push(vb);
            projection_first.push(dv);
            projection_second.push(ddv);
        }
        Ok((
            Self {
                filter,
                reference_cell_masses_kg: elasticity.reference_cell_masses_kg(),
                projected,
                volume_projected,
                first,
                second,
                stiffness_projection_first,
                projection_first,
                projection_second,
            },
            moduli,
        ))
    }

    pub(super) fn pullback_gradient(
        &self,
        gradient: &[f64],
        stage: &'static str,
        quantity: &'static str,
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        let mut local = Vec::with_capacity(self.first.len());
        for e in 0..self.first.len() {
            if e % 256 == 0 {
                control.checkpoint(stage)?;
            }
            local.push(finite(self.first[e] * gradient[e], quantity)?);
        }
        Ok(self.filter.try_apply_transpose(&local, control)?)
    }

    /// Pull both stiffness and actual mechanical occupancy through ONE
    /// transpose filter. Robust inventory dilation is not an inertial field.
    pub(super) fn pullback_stiffness_and_mass_gradient(
        &self,
        modulus: &[f64],
        mass_fraction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        if modulus.len() != self.first.len() || mass_fraction.len() != self.first.len() {
            return Err(ComplianceError::InvalidInput(
                "modal design gradient shape mismatch",
            ));
        }
        let mut local = Vec::with_capacity(self.first.len());
        for e in 0..self.first.len() {
            if e % 256 == 0 {
                control.checkpoint("simp-thermal-frequency-pullback")?;
            }
            local.push(finite(
                self.first[e] * modulus[e] + self.stiffness_projection_first[e] * mass_fraction[e],
                "SIMP thermal frequency gradient",
            )?);
        }
        Ok(self.filter.try_apply_transpose(&local, control)?)
    }

    pub(super) fn tangent(
        &self,
        direction: &[f64],
        stage: &'static str,
        quantity: &'static str,
        control: &mut SolveControl<'_>,
    ) -> Result<(Vec<f64>, Vec<f64>), ComplianceError> {
        let filtered = self.filter.try_apply(direction, control)?;
        let mut tangent = Vec::with_capacity(direction.len());
        for (e, &v) in filtered.iter().enumerate() {
            if e % 256 == 0 {
                control.checkpoint(stage)?;
            }
            tangent.push(finite(self.first[e] * v, quantity)?);
        }
        Ok((filtered, tangent))
    }

    pub(super) fn pullback_hessian(
        &self,
        filtered: &[f64],
        product: &[f64],
        gradient: &[f64],
        stage: &'static str,
        quantity: &'static str,
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        let mut local = Vec::with_capacity(filtered.len());
        for e in 0..filtered.len() {
            if e % 256 == 0 {
                control.checkpoint(stage)?;
            }
            local.push(finite(
                self.first[e] * product[e] + self.second[e] * filtered[e] * gradient[e],
                quantity,
            )?);
        }
        Ok(self.filter.try_apply_transpose(&local, control)?)
    }

    pub(super) fn direction(
        &self,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<(), ComplianceError> {
        control.checkpoint("simp-direction")?;
        if direction.len() != self.first.len() {
            return Err(ComplianceError::InvalidInput(
                "raw-density Hessian direction shape mismatch",
            ));
        }
        for (e, &v) in direction.iter().enumerate() {
            if e % 256 == 0 {
                control.checkpoint("simp-direction")?;
            }
            finite(v, "raw-density direction")?;
        }
        Ok(())
    }

    fn weights(
        &self,
        volumes: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        control.checkpoint("simp-volume-admission")?;
        if volumes.len() != self.first.len() {
            return Err(ComplianceError::InvalidInput(
                "one physical volume per design cell required",
            ));
        }
        let mut total = 0.0;
        for (e, &volume) in volumes.iter().enumerate() {
            if e % 256 == 0 {
                control.checkpoint("simp-volume-admission")?;
            }
            if !(volume.is_finite() && volume > 0.0) {
                return Err(ComplianceError::InvalidInput(
                    "cell volumes must be finite and positive",
                ));
            }
            total = finite(total + volume, "total volume")?;
        }
        Ok(volumes.iter().map(|volume| volume / total).collect())
    }

    /// Physical projected volume. For robust goals this is DILATED volume,
    /// not eroded stiffness, nominal material, or an unweighted cell mean.
    pub(super) fn volume_value_gradient(
        &self,
        volumes: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<(f64, Vec<f64>), ComplianceError> {
        let weights = self.weights(volumes, control)?;
        self.projected_quantity_value_gradient(&weights, "simp-volume-gradient", control)
    }

    /// Material inventory in kilograms, with its exact raw-occupancy gradient.
    /// Uses the stiffness owner's reference cell masses, without normalization
    /// or modal void-mass interpolation. Robust goals use DILATED occupancy.
    pub(super) fn mass_value_gradient(
        &self,
        control: &mut SolveControl<'_>,
    ) -> Result<(f64, Vec<f64>), ComplianceError> {
        let masses = self
            .reference_cell_masses_kg
            .ok_or(ComplianceError::InvalidInput(
                "physical mass requires an assigned material field",
            ))?;
        self.projected_quantity_value_gradient(masses, "simp-mass-gradient", control)
    }

    fn projected_quantity_value_gradient(
        &self,
        weights: &[f64],
        stage: &'static str,
        control: &mut SolveControl<'_>,
    ) -> Result<(f64, Vec<f64>), ComplianceError> {
        let mut value = 0.0;
        let mut local = Vec::with_capacity(weights.len());
        for (e, &w) in weights.iter().enumerate() {
            if e % 256 == 0 {
                control.checkpoint(stage)?;
            }
            value += w * self.volume_projected[e];
            local.push(finite(w * self.projection_first[e], stage)?);
        }
        let value = finite(value, stage)?;
        Ok((value, self.filter.try_apply_transpose(&local, control)?))
    }

    /// F^T diag(V/sum(V) * H_volume'') F v; no SIMP term belongs in volume.
    pub(super) fn volume_hessian_vector_product(
        &self,
        volumes: &[f64],
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        self.direction(direction, control)?;
        let weights = self.weights(volumes, control)?;
        self.projected_quantity_hessian_vector_product(
            &weights,
            direction,
            "simp-volume-hessian",
            control,
        )
    }

    /// F^T diag(m_ref * H_inventory'') F v, in kilograms. Includes the same
    /// robust dilation as the mass value; no stiffness or modal mass derivative.
    pub(super) fn mass_hessian_vector_product(
        &self,
        direction: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        self.direction(direction, control)?;
        let masses = self
            .reference_cell_masses_kg
            .ok_or(ComplianceError::InvalidInput(
                "physical mass requires an assigned material field",
            ))?;
        self.projected_quantity_hessian_vector_product(
            masses,
            direction,
            "simp-mass-hessian",
            control,
        )
    }

    fn projected_quantity_hessian_vector_product(
        &self,
        weights: &[f64],
        direction: &[f64],
        stage: &'static str,
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, ComplianceError> {
        if direction.iter().all(|&v| v == 0.0) || self.projection_second.iter().all(|&h| h == 0.0) {
            return Ok(vec![0.0; direction.len()]);
        }
        let filtered = self.filter.try_apply(direction, control)?;
        let mut local = Vec::with_capacity(direction.len());
        for e in 0..direction.len() {
            if e % 256 == 0 {
                control.checkpoint(stage)?;
            }
            local.push(finite(
                weights[e] * self.projection_second[e] * filtered[e],
                stage,
            )?);
        }
        Ok(self.filter.try_apply_transpose(&local, control)?)
    }
}
