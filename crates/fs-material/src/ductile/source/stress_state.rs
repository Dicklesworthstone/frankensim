//! Lode-dependent initiation on the original two-coordinate source cells.

use fs_matdb::{InterpolationPolicy, PreparedSurfaceResponse, SurfaceDifferential};

use super::*;

/// Prepared local initiation strain on original triaxiality/Lode source cells.
/// Cell lookup and analytical gradients reuse the common material-data owner.
/// At shared edges and vertices its first incident source cell selects a
/// semismooth branch; no classical derivative there is asserted.
#[derive(Debug, Clone)]
pub struct DuctileStressSurface {
    response: PreparedSurfaceResponse,
    domain: DuctileDomain,
    triaxiality_index: usize,
    lode_bounds: [f64; 2],
}

impl PartialEq for DuctileStressSurface {
    fn eq(&self, other: &Self) -> bool {
        self.domain == other.domain
            && self.triaxiality_index == other.triaxiality_index
            && self.lode_bounds == other.lode_bounds
            && self.response.axes() == other.response.axes()
            && self.response.samples() == other.response.samples()
            && self.response.triangles() == other.response.triangles()
    }
}

impl DuctileStressSurface {
    /// Original source-axis order, sample coordinates, and explicit cells.
    pub const fn response(&self) -> &PreparedSurfaceResponse {
        &self.response
    }

    pub const fn lode_bounds(&self) -> [f64; 2] {
        self.lode_bounds
    }

    fn evaluate(
        &self,
        stress: DuctileStressState,
        rate: f64,
        temperature: f64,
    ) -> Result<SurfaceDifferential, DuctileError> {
        self.domain.check(stress.triaxiality, rate, temperature)?;
        super::super::in_domain("lode_parameter", stress.lode_parameter, self.lode_bounds)?;
        let mut point = [0.0; 2];
        point[self.triaxiality_index] = stress.triaxiality;
        point[1 - self.triaxiality_index] = stress.lode_parameter;
        let value = self.response.evaluate(point).map_err(|_| {
            DuctileError::Invalid(
                "actual stress state has no representable response on the original ductile source cells",
            )
        })?;
        if !super::super::positive(value.value) || value.gradient.iter().any(|v| !v.is_finite()) {
            return Err(DuctileError::Numerical(
                "source-cell initiation strain or derivative is not representable",
            ));
        }
        Ok(value)
    }
}

impl DuctileInitiationLaw for DuctileStressSurface {
    fn domain(&self) -> &DuctileDomain {
        &self.domain
    }

    fn fracture_strain(&self, _: f64, _: f64, _: f64) -> Result<f64, DuctileError> {
        Err(DuctileError::Invalid(
            "a Lode-dependent initiation locus requires the complete stress state",
        ))
    }

    fn fracture_strain_at(
        &self,
        stress: DuctileStressState,
        rate: f64,
        temperature: f64,
    ) -> Result<f64, DuctileError> {
        Ok(self.evaluate(stress, rate, temperature)?.value)
    }
}

impl DifferentiableDuctileInitiationLaw for DuctileStressSurface {
    fn fracture_strain_gradient(
        &self,
        _: f64,
        _: f64,
        _: f64,
    ) -> Result<DuctileLocusGradient, DuctileError> {
        Err(DuctileError::Invalid(
            "a Lode-dependent initiation gradient requires the complete stress state",
        ))
    }

    fn fracture_strain_gradient_at(
        &self,
        stress: DuctileStressState,
        rate: f64,
        temperature: f64,
    ) -> Result<DuctileLocusGradient, DuctileError> {
        let evaluated = self.evaluate(stress, rate, temperature)?;
        Ok(DuctileLocusGradient {
            triaxiality: evaluated.gradient[self.triaxiality_index],
            lode_parameter: evaluated.gradient[1 - self.triaxiality_index],
            strain_rate_s_inv: 0.0,
        })
    }
}

/// Prepare an original two-coordinate ductile initiation locus.
///
/// The exact scalar ordinate is local equivalent plastic strain at initiation;
/// axes are tensile-positive stress triaxiality and the named Lode parameter
/// `xi = cos(3 theta)`. Either source-axis order is accepted without changing
/// sample order, triangles, or evidence. Only declared `LinearInside` cells
/// supply a runtime response; holes and isolated samples gain no support.
/// The original source must bound the positive equivalent plastic strain rate
/// and a semantic absolute temperature. Temperature is fixed at the selected
/// observation; the rate interval is applicability, not a fitted rate law.
///
/// The request bounds retained source samples. The shared surface owner also
/// enforces its existing sample/cell caps. Preparation queries each original
/// sample once and verifies its receipt; runtime trials never query the card.
pub fn resolve_ductile_stress_locus(
    card: &MaterialCard,
    point: &QueryPoint,
    request: &DuctileSourceRequest,
    selection: ClaimSelection,
) -> Result<ResolvedDuctileLocus, DuctileSourceError> {
    require(
        request.maximum_points >= 3
            && !request.temperature_axis.trim().is_empty()
            && ![
                STRESS_TRIAXIALITY_AXIS,
                LODE_PARAMETER_AXIS,
                PLASTIC_STRAIN_RATE_AXIS,
            ]
            .contains(&request.temperature_axis.as_str()),
        "declare a distinct temperature axis and at least three stress-state samples",
    )?;
    let key = PropertyKey::new(FRACTURE_PLASTIC_STRAIN_PROPERTY, Dims::NONE);
    let selected = query(card, &key, point, selection)?;
    let claim =
        card.claims()
            .claim(selected.receipt.selected)
            .ok_or(DuctileSourceError::Invalid(
                "selected ductile source is absent",
            ))?;
    let PropertyValue::Surface2D {
        axes,
        axis_dims,
        samples,
        triangles,
        dims,
    } = &claim.value
    else {
        return Err(DuctileSourceError::Invalid(
            "original equivalent-plastic-strain surface on triaxiality and Lode parameter required",
        ));
    };
    let eta_index = axes
        .iter()
        .position(|axis| axis == STRESS_TRIAXIALITY_AXIS)
        .ok_or(DuctileSourceError::Invalid(
            "stress-state source lacks triaxiality",
        ))?;
    require(
        axes[1 - eta_index] == LODE_PARAMETER_AXIS
            && *axis_dims == [Dims::NONE; 2]
            && *dims == Dims::NONE
            && claim.interpolation == InterpolationPolicy::LinearInside,
        "dimensionless triaxiality/cos(3 theta) axes and original affine-cell interpolation required",
    )?;
    require(
        samples.len() >= 3 && samples.len() <= request.maximum_points && !triangles.is_empty(),
        "stress-state sample count or original cell support violates the declared preparation budget",
    )?;
    for (axis, quantity) in [
        (
            STRESS_TRIAXIALITY_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
        ),
        (LODE_PARAMETER_AXIS, QuantitySpec::dimensional(Dims::NONE)),
        (
            PLASTIC_STRAIN_RATE_AXIS,
            QuantitySpec::dimensional(RATE_DIMS),
        ),
    ] {
        require(
            point.axis_quantities().get(axis) == Some(&quantity)
                && claim.validity.axis_quantities().get(axis) == Some(&quantity)
                && claim.validity.bounds().contains_key(axis),
            "source stress invariants and equivalent-plastic-strain rate need their exact typed bounds",
        )?;
    }
    let temperature = point
        .axis_quantities()
        .get(&request.temperature_axis)
        .ok_or(DuctileSourceError::Invalid("source temperature is untyped"))?;
    require(
        temperature.semantic_type().is_some_and(|q| {
            q.kind() == QuantityKind::AbsoluteTemperature && q.form() == ValueForm::Static
        }) && claim
            .validity
            .axis_quantities()
            .get(&request.temperature_axis)
            == Some(temperature)
            && claim
                .validity
                .bounds()
                .contains_key(&request.temperature_axis),
        "source temperature must be declared static absolute kelvin",
    )?;
    let temperature_k = *point
        .axes()
        .get(&request.temperature_axis)
        .ok_or(DuctileSourceError::Invalid("source temperature is absent"))?;
    let anchor = [
        *point
            .axes()
            .get(&axes[0])
            .ok_or(DuctileSourceError::Invalid(
                "first stress coordinate is absent",
            ))?,
        *point
            .axes()
            .get(&axes[1])
            .ok_or(DuctileSourceError::Invalid(
                "second stress coordinate is absent",
            ))?,
    ];
    require(
        samples.iter().any(|(at, _)| *at == anchor),
        "select an original stress-state sample, not an interpolated source anchor",
    )?;
    for axis in claim.validity.bounds().keys() {
        if [
            STRESS_TRIAXIALITY_AXIS,
            LODE_PARAMETER_AXIS,
            PLASTIC_STRAIN_RATE_AXIS,
        ]
        .contains(&axis.as_str())
            || axis == &request.temperature_axis
        {
            continue;
        }
        require(
            is_fixed_external_coordinate(
                axis,
                point.axis_quantities().get(axis).copied(),
                point.axes().get(axis).copied(),
            ) && point.axis_quantities().get(axis).is_none_or(|quantity| {
                quantity
                    .semantic_type()
                    .is_none_or(|q| q.kind() != QuantityKind::AbsoluteTemperature)
            }) && axis != "T",
            "a stress-state ductile locus cannot freeze another evolving or temperature axis",
        )?;
    }
    let &(minimum_rate, maximum_rate) = claim
        .validity
        .bounds()
        .get(PLASTIC_STRAIN_RATE_AXIS)
        .ok_or(DuctileSourceError::Invalid(
            "plastic strain rate bounds are absent",
        ))?;
    let mut eta_bounds = [f64::INFINITY, f64::NEG_INFINITY];
    let mut lode_bounds = [f64::INFINITY, f64::NEG_INFINITY];
    let mut observations = Vec::with_capacity(samples.len());
    for &(coordinates, strain) in samples {
        let eta = coordinates[eta_index];
        let xi = coordinates[1 - eta_index];
        require(
            super::super::positive(strain) && (-1.0..=1.0).contains(&xi),
            "every original sample needs positive initiation strain and cos(3 theta) in [-1,1]",
        )?;
        eta_bounds = [eta_bounds[0].min(eta), eta_bounds[1].max(eta)];
        lode_bounds = [lode_bounds[0].min(xi), lode_bounds[1].max(xi)];
        let at = point
            .clone()
            .with(&axes[0], coordinates[0])?
            .with(&axes[1], coordinates[1])?;
        let answer = query(card, &key, &at, selection)?;
        require(
            answer.receipt.selected == selected.receipt.selected
                && answer.evidence.value.value == strain,
            "every original stress-state sample must retain its selected source and exact plastic strain",
        )?;
        observations.push((at, answer));
    }
    let domain = DuctileDomain {
        triaxiality: eta_bounds,
        strain_rate_s_inv: [minimum_rate, maximum_rate],
        temperature_k: [temperature_k; 2],
    };
    domain.validate()?;
    let response = claim.value.prepare_surface_response()?;
    let law = DuctileStressSurface {
        response,
        domain,
        triaxiality_index: eta_index,
        lode_bounds,
    };
    Ok(ResolvedDuctileLocus {
        law: PreparedDuctileLocus::StressState(law),
        source: DuctileCalibrationSource {
            card_identity: card.content_hash(),
            material_state: card.id().clone(),
            claim: claim.clone(),
            request: request.clone(),
            query: point.clone(),
            selection,
            observations,
        },
    })
}
