//! Source-resolved in-plane Kelvin heat on each original thermal sublamina.

use super::*;
use fs_material::plane_stress_viscosity::{MaterialPlane, PlaneStressViscosity};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_plate::shell::mitc::viscosity::{ShearShellViscosity, projection_terms};
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{Frequency, QuantitySpec};

/// Applicability and cold assembly allowance for explicit material loss.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminateViscosityControls {
    /// Declared motion bandwidth, in rad/s. Every selected source must cover
    /// this interval at all prepared temperatures. This caller declaration is
    /// not a spectral certificate for the full native transient trajectory.
    pub angular_frequency_band_rad_s: [f64; 2],
    pub maximum_assembly_terms: usize,
}

pub(super) struct LayerViscosity {
    reference: Vec<PlaneStressViscosity>,
    controls: LaminateViscosityControls,
    pub(super) accepted: PreparedLayerViscosity,
}

#[derive(Clone)]
pub(super) struct PreparedLayerViscosity {
    pub(super) sources: Arc<[PlaneStressViscosity]>,
    pub(super) facets: Arc<[ShearShellViscosity]>,
}

impl<'s, 'b, 'a> DynamicLayeredThermalShell<'s, 'b, 'a> {
    /// Bind one source-resolved in-plane viscosity to each ORIGINAL thermal
    /// cell / mechanical sublamina. Original cards, physical reference query,
    /// layer order, orientation and reference faces remain authoritative.
    /// Frequency applicability coordinates may accompany that same physical
    /// query; they do not change the temperature or other material conditions.
    ///
    /// Each interval evaluates the pinned coefficients at the arithmetic mean
    /// of its accepted and trial temperatures, then freezes that one passive
    /// operator throughout the native midpoint step. The source law is
    /// sigma_v = V * total engineering strain rate. No thermal eigenstrain-rate
    /// dashpot or hereditary state is introduced. Reduced sources characterize
    /// in-plane loss only: transverse viscous resultants are explicitly zero,
    /// not inferred from elastic shear moduli or a layerwise correction factor.
    ///
    /// Physical heat uses each source's original positive work contraction on
    /// the actual midpoint strain-rate field, and returns to that same thermal
    /// cell. It is separately checked against the native mechanical heat.
    /// Rebinding adds no stiffness, inertia or stored energy. Failure changes
    /// no accepted thermomechanical state.
    pub fn with_layer_viscosity(
        mut self,
        sources: Vec<PlaneStressViscosity>,
        controls: LaminateViscosityControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LaminateThermoelasticError> {
        check(&mut cancelled)?;
        let [lo, hi] = controls.angular_frequency_band_rad_s;
        if sources.len() != self.laws.len()
            || sources.is_empty()
            || !lo.is_finite()
            || !hi.is_finite()
            || lo < 0.0
            || hi < lo
            || controls.maximum_assembly_terms == 0
            || projection_terms(
                self.source.discretization.mesh.nodes.len(),
                self.source.discretization.mesh.tris.len(),
                0,
            )
            .is_none_or(|work| work > controls.maximum_assembly_terms)
        {
            return Err(invalid(
                "bounded source viscosity required for every original thermal cell",
            ));
        }
        for ((law, layer), thermoelastic) in sources
            .iter()
            .zip(self.source.source.layers())
            .zip(&self.laws)
        {
            check(&mut cancelled)?;
            if law.plane() != MaterialPlane::Axes12
                || law.source().card_identity() != layer.input().card.content_hash()
                || !same_physical_query(law.source(), layer.reference_material())
            {
                return Err(invalid(
                    "layer viscosity must retain its original card, reference state and material axes 1,2",
                ));
            }
            for temperature in thermoelastic.temperature_range_k() {
                let _ = refresh(law, layer, temperature, controls)?;
            }
        }
        let accepted = prepare(
            self.source,
            &sources,
            self.state.temperatures_k(),
            controls,
            &mut cancelled,
        )?;
        let native = self
            .native
            .with_viscosity(
                &accepted.facets,
                controls.maximum_assembly_terms,
                &mut cancelled,
            )
            .map_err(LaminateThermoelasticError::Native)?;
        check(&mut cancelled)?;
        self.native = native;
        self.midpoint = None;
        self.viscosity = Some(LayerViscosity {
            reference: sources,
            controls,
            accepted,
        });
        Ok(self)
    }

    /// Last accepted interval's actual frozen sources, in original cell order.
    /// Before the first interval these are evaluated at the launch temperature.
    pub fn layer_viscosity_sources(&self) -> Option<&[PlaneStressViscosity]> {
        self.viscosity
            .as_ref()
            .map(|loss| loss.accepted.sources.as_ref())
    }
}

impl LayerViscosity {
    pub(super) fn prepare_midpoint(
        &self,
        source: &LayeredThermalShell<'_, '_>,
        before: &[f64],
        after: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedLayerViscosity, LaminateThermoelasticError> {
        if before.len() != self.reference.len() || after.len() != before.len() {
            return Err(invalid(
                "viscosity temperatures must cover every original cell",
            ));
        }
        let temperatures: Vec<_> = before
            .iter()
            .zip(after)
            .map(|(a, b)| 0.5 * a + 0.5 * b)
            .collect();
        if temperatures
            .iter()
            .zip(self.accepted.sources.iter())
            .all(|(t, law)| {
                law.source()
                    .state_coordinate("T")
                    .is_some_and(|old| old.to_bits() == t.to_bits())
            })
        {
            check(cancelled)?;
            return Ok(self.accepted.clone());
        }
        prepare(
            source,
            &self.reference,
            &temperatures,
            self.controls,
            cancelled,
        )
    }

    pub(super) fn interval_native(
        &self,
        before: &EnrichedMitcShell,
        prepared: &PreparedLayerViscosity,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EnrichedMitcShell, LaminateThermoelasticError> {
        if prepared.same_operator(&self.accepted) {
            return Ok(before.clone());
        }
        before
            .with_viscosity(
                &prepared.facets,
                self.controls.maximum_assembly_terms,
                cancelled,
            )
            .map_err(LaminateThermoelasticError::Native)
    }
}

impl PreparedLayerViscosity {
    pub(super) fn same_operator(&self, other: &Self) -> bool {
        self.facets.len() == other.facets.len()
            && self
                .facets
                .iter()
                .zip(other.facets.iter())
                .all(|(a, b)| a.coupled_matrices() == b.coupled_matrices())
    }

    pub(super) fn cell_heat(
        &self,
        source: &LayeredThermalShell<'_, '_>,
        rates: &[RealElementStrain],
        quadrature: &[([f64; 3], f64)],
        dt_s: f64,
        native_heat_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, LaminateThermoelasticError> {
        use fs_math::dd::Dd;
        if rates.len() != self.facets.len() {
            return Err(invalid(
                "native rates must cover every original viscous facet",
            ));
        }
        let section = source.source.reference_section();
        let mut heat = vec![Dd::ZERO; self.sources.len()];
        for (element, field) in rates.iter().enumerate() {
            check(cancelled)?;
            let area = source
                .discretization
                .mesh
                .facet(element)
                .map_err(LaminateThermoelasticError::Native)?
                .area_m2;
            let section_angle = source.discretization.section_angles_rad[element];
            for &(lambda, weight) in quadrature {
                let point = field
                    .at_barycentric(lambda)
                    .map_err(LaminateThermoelasticError::Native)?;
                for (layer, ((law, geometry), total)) in self
                    .sources
                    .iter()
                    .zip(section.layer_integrals())
                    .zip(&mut heat)
                    .enumerate()
                {
                    check(cancelled)?;
                    let angle = section_angle + section.layers()[layer].orientation_rad;
                    let rotation = [det::cos(angle), det::sin(angle)];
                    let half = 0.5 * (geometry.z_top_m - geometry.z_bottom_m);
                    let middle = geometry.z_bottom_m + half;
                    let offset = half / det::sqrt(3.0);
                    for z in [middle - offset, middle + offset] {
                        let (_, power) = law
                            .stress_and_power(principal(point, z, rotation))
                            .map_err(LaminateThermoelasticError::Viscosity)?;
                        let joules = finite(dt_s * area * weight * half * power)?;
                        if joules < 0.0 {
                            return Err(invalid("negative source viscous heat"));
                        }
                        *total = *total + Dd::from_f64(joules);
                    }
                }
            }
        }
        let heat: Vec<_> = heat.into_iter().map(Dd::to_f64).collect();
        let total = heat
            .iter()
            .fold(Dd::ZERO, |sum, &v| sum + Dd::from_f64(v))
            .to_f64();
        let scale = total.abs().max(native_heat_j.abs());
        // Independent original volume/source power versus integrated native
        // A/B/D power. This allowance scales with physical heat, never with
        // stored energy or a user-selected thermomechanical error tolerance.
        let roundoff = 16384.0 * f64::EPSILON * scale + 64.0 * f64::from_bits(1);
        if !total.is_finite()
            || !native_heat_j.is_finite()
            || native_heat_j < 0.0
            || (total - native_heat_j).abs() > roundoff
        {
            return Err(invalid(
                "original layer viscous heat and native mechanical loss disagree",
            ));
        }
        check(cancelled)?;
        Ok(heat)
    }
}

fn prepare(
    source: &LayeredThermalShell<'_, '_>,
    reference: &[PlaneStressViscosity],
    temperatures: &[f64],
    controls: LaminateViscosityControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedLayerViscosity, LaminateThermoelasticError> {
    let mut sources = Vec::with_capacity(reference.len());
    for ((law, layer), &temperature) in reference
        .iter()
        .zip(source.source.layers())
        .zip(temperatures)
    {
        check(cancelled)?;
        sources.push(refresh(law, layer, temperature, controls)?);
    }
    let mut facets = Vec::with_capacity(source.discretization.section_angles_rad.len());
    for &angle in &source.discretization.section_angles_rad {
        check(cancelled)?;
        let matrices = sources
            .iter()
            .map(|law| {
                law.rotated_matrix_pa_s(angle)
                    .map_err(LaminateThermoelasticError::Viscosity)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let (a, b, d) = source
            .source
            .reference_section()
            .integrate_layer_matrices(&matrices)
            .map_err(|_| invalid("unrepresentable original layer viscosity resultants"))?;
        facets.push(
            ShearShellViscosity::new_coupled(a, b, d, [0.0; 3])
                .map_err(LaminateThermoelasticError::Native)?,
        );
    }
    check(cancelled)?;
    Ok(PreparedLayerViscosity {
        sources: sources.into(),
        facets: facets.into(),
    })
}

fn refresh(
    original: &PlaneStressViscosity,
    layer: &fs_material::laminate::thermal::source::BoundLaminaSource,
    temperature: f64,
    controls: LaminateViscosityControls,
) -> Result<PlaneStressViscosity, LaminateThermoelasticError> {
    let point = original
        .source()
        .query()
        .with("T", temperature)
        .map_err(|_| invalid("invalid original viscosity temperature coordinate"))?;
    let source = original
        .source()
        .at_query_point(&layer.input().card, &point)
        .map_err(LaminateThermoelasticError::ViscosityQuery)?;
    let law = PlaneStressViscosity::from_material(
        source,
        MaterialPlane::Axes12,
        core::array::from_fn(|i| original.properties()[i].as_str()),
    )
    .map_err(LaminateThermoelasticError::Viscosity)?;
    for omega in controls.angular_frequency_band_rad_s {
        law.validate_frequency(omega)
            .map_err(LaminateThermoelasticError::Viscosity)?;
    }
    Ok(law)
}

fn same_physical_query(a: &ResolvedMaterialStatePoint, b: &ResolvedMaterialStatePoint) -> bool {
    let physical = |state: &ResolvedMaterialStatePoint, axis: &str| {
        let quantity = state.axis_quantities().get(axis).copied();
        let frequency = quantity
            .and_then(QuantitySpec::semantic_type)
            .is_some_and(|ty| {
                ty.form() == ValueForm::Static && matches!(ty.kind(), QuantityKind::Frequency(_))
            })
            || axis == "omega"
                && quantity.is_none_or(|q| q == QuantitySpec::dimensional(Frequency::DIMS));
        !frequency
    };
    let mut left = a.query_point().iter().filter(|(axis, _)| physical(a, axis));
    let mut right = b.query_point().iter().filter(|(axis, _)| physical(b, axis));
    loop {
        match (left.next(), right.next()) {
            (None, None) => return true,
            (Some((ka, va)), Some((kb, vb)))
                if ka == kb
                    && va.to_bits() == vb.to_bits()
                    && a.axis_quantities().get(ka) == b.axis_quantities().get(kb) => {}
            _ => return false,
        }
    }
}
