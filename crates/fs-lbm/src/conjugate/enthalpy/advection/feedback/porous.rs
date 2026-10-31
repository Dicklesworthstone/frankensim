//! Source-bound, stationary, grid-aligned porous drag and its original heat.

use super::*;
use fs_material::porous::permeability::temperature::ResolvedPermeability;
use fs_qty::Temperature;

/// One accepted source result; the original material bundle retains all keys,
/// selected claims, fixed conditions and actual local temperature.
#[derive(Clone, Debug, PartialEq)]
pub struct PorousFlowSample {
    pub voxel: usize,
    pub cell: CellId,
    pub material: ResolvedPermeability,
}

impl<'a, F: Clone, R, O: DepositionObserver, T: FlowTransportSnapshot>
    MaterialFlowFeedback<'a, F, R, O, T>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    /// Bind one original permeability law to each declared porous fluid voxel.
    /// None is an explicitly clear-fluid voxel; missing material data must not
    /// be converted to None. Principal axes must be exact signed grid axes or
    /// their permutation. Oblique tensors require a different momentum stencil
    /// and refuse, even if their sampled eigenvalues happen to coincide.
    ///
    /// Actual local mu(T) divided by K_i(T) enters the SAME implicit native
    /// momentum coefficient and optional Darcy heating quadrature. Scalar mushy
    /// resistance remains additive. Both source-start and accepted endpoint T
    /// must fit the original claims; no clipping or replacement source is used.
    ///
    /// The rigid skeleton has zero velocity, so drag work is physical loss.
    /// With a heating mode it enters the original cell's registered caloric
    /// storage once under local thermal equilibrium. That chart must represent
    /// the intended participating heat capacity; this binding invents no extra
    /// skeleton store or mass. Fixed geometry/density and linear Darcy flow
    /// remain assumptions. Sealed axes refuse without a velocity constraint.
    /// Attach once before the first advance; no live skeleton replacement or
    /// momentum reset is supplied.
    pub fn with_porous_materials(
        mut self,
        bindings: Vec<Option<PermeabilityTemperatureLaw<'a>>>,
    ) -> Result<Self, MaterialFlowError<O::Error>> {
        check_gate(self.gate, &mut || Ok(()))?;
        if self.porous.is_some()
            || self.applied.is_some()
            || bindings.len() != self.domain.cell_count()
        {
            return Err(invalid(
                "one bounded original porous map is required before the first advance",
            )
            .into());
        }
        let mut has_source = false;
        for (voxel, law) in bindings.iter().enumerate() {
            if voxel % 256 == 0 {
                check_gate(self.gate, &mut || Ok(()))?;
            }
            if law.is_some() && !self.domain.is_fluid(voxel) {
                return Err(
                    invalid("porous skeleton bindings belong to mapped fluid voxels").into(),
                );
            }
            has_source |= law.is_some();
        }
        if !has_source {
            return Err(
                invalid("a porous binding needs at least one original skeleton source").into(),
            );
        }
        self.porous = Some(bindings);
        self.check_porous_endpoint(&mut || Ok(()))?;
        Ok(self)
    }

    pub fn endpoint_porous_samples(&self) -> &[PorousFlowSample] {
        &self.porous_endpoint
    }

    pub(super) fn check_porous_endpoint(
        &mut self,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), MaterialFlowError<O::Error>> {
        self.porous_endpoint = self.porous_drag(&self.endpoint, None, check)?.0;
        Ok(())
    }

    pub(super) fn porous_drag(
        &self,
        samples: &[FlowMaterialSample<F, T>],
        mushy: Option<&[f64]>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Vec<PorousFlowSample>, Option<Vec<[f64; 3]>>), MaterialFlowError<O::Error>> {
        let Some(bindings) = &self.porous else {
            return Ok((Vec::new(), None));
        };
        let mut drag = vec![[0.0; 3]; self.domain.cell_count()];
        if let Some(mushy) = mushy {
            for (voxel, (row, value)) in drag.iter_mut().zip(mushy).enumerate() {
                if voxel % 256 == 0 {
                    check_gate(self.gate, check)?;
                }
                *row = [*value; 3];
            }
        }
        let mut sources = Vec::new();
        for sample in samples {
            check_gate(self.gate, check)?;
            let Some(law) = &bindings[sample.voxel] else {
                continue;
            };
            let material = law
                .at_temperature(Temperature::new(sample.temperature_k))
                .map_err(|source| MaterialFlowError::Material {
                    cell: sample.cell,
                    source,
                })?;
            let permeability = material.permeability();
            let mut grid_k = [0.0; 3];
            for (axis, &k) in permeability
                .material_axes_in_world()
                .iter()
                .zip(permeability.principal_m2())
            {
                let direction = axis
                    .iter()
                    .position(|x| x.abs() == 1.0)
                    .filter(|&a| (0..3).all(|i| i == a || axis[i] == 0.0))
                    .ok_or(invalid(
                        "native porous momentum requires exact signed grid principal axes",
                    ))?;
                if !k.is_finite() || k <= 0.0 {
                    return Err(invalid("native Darcy drag requires positive intrinsic permeability; sealed axes need explicit constraints").into());
                }
                grid_k[direction] = k;
            }
            if let Some(transport) = sample.transport.mobile_transport() {
                for (coefficient, k) in drag[sample.voxel].iter_mut().zip(grid_k) {
                    let value = transport.dynamic_viscosity_pa_s() / k;
                    *coefficient += value;
                    if !value.is_finite() || value <= 0.0 || !coefficient.is_finite() {
                        return Err(invalid(
                            "source-bound directional Darcy drag is not representable",
                        )
                        .into());
                    }
                }
            }
            sources.push(PorousFlowSample {
                voxel: sample.voxel,
                cell: sample.cell,
                material,
            });
        }
        Ok((sources, Some(drag)))
    }
}
