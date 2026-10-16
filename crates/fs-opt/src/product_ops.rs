//! Differential operations on ordered products of the existing manifolds.
//! Factor-local geometry remains owned by `Manifold`; this module only slices,
//! delegates and concatenates. No aggregate output is returned on a refusal.

use crate::{
    OptError, ProductCoordinate, ProductFactorLayout, ProductManifold,
    ProductManifoldError, RetractionCurve,
};
use core::ops::Range;

/// A factor/shape refusal or failure to reserve a complete differential output.
#[derive(Debug, Clone, PartialEq)]
pub enum ProductDifferentialError {
    /// Original product refusal, including stable factor identity and index.
    Geometry(ProductManifoldError),
    /// Aggregate output allocation failed before factor evaluation.
    Allocation {
        /// Requested scalar elements.
        elements: usize,
    },
}

impl From<ProductManifoldError> for ProductDifferentialError {
    fn from(error: ProductManifoldError) -> Self { Self::Geometry(error) }
}

impl core::fmt::Display for ProductDifferentialError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Geometry(error) => write!(f, "{error}"),
            Self::Allocation { elements } => write!(f, "product differential allocation refused for {elements} scalars"),
        }
    }
}

impl std::error::Error for ProductDifferentialError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Geometry(error) => Some(error), _ => None }
    }
}

fn storage(len: usize) -> Result<Vec<f64>, ProductDifferentialError> {
    let mut result = Vec::new();
    result.try_reserve_exact(len)
        .map_err(|_| ProductDifferentialError::Allocation { elements: len })?;
    Ok(result)
}

fn point_range(block: &ProductFactorLayout) -> Range<usize> {
    let start = block.point_offset().get() as usize;
    start..start + block.manifold_layout().point_dim().get() as usize
}

fn parameter_range(block: &ProductFactorLayout) -> Range<usize> {
    let start = block.param_offset().get() as usize;
    start..start + block.manifold_layout().param_dim().get() as usize
}

fn factor_error(block: &ProductFactorLayout, source: OptError) -> ProductDifferentialError {
    ProductManifoldError::FactorOperation {
        id: block.factor().id(), index: block.index(), source,
    }.into()
}

impl ProductManifold {
    /// Pull back an ambient gradient to the product's retraction parameters.
    ///
    /// Point and ambient payloads must both match the total point dimension.
    /// SO(3) contributes three body coordinates for four quaternion coordinates;
    /// Sphere/Stiefel use the existing embedded-metric tangent projections.
    /// Reductions, domain checks and scaling are exactly the factor operations.
    /// This is a product metric, not a user-selected physical preconditioner.
    pub fn parameter_gradient(
        &self, point: &[f64], ambient: &[f64],
    ) -> Result<Vec<f64>, ProductDifferentialError> {
        self.layout().validate_payload_len(ProductCoordinate::Point, point)?;
        self.layout().validate_payload_len(ProductCoordinate::Point, ambient)?;
        let mut output = storage(self.layout().param_dim().get() as usize)?;
        for block in self.layout().factors() {
            let range = point_range(block);
            let gradient = block.factor().manifold()
                .parameter_gradient(&point[range.clone()], &ambient[range])
                .map_err(|e| factor_error(block, e))?;
            output.extend_from_slice(&gradient);
        }
        Ok(output)
    }

    /// Check each parameter block against its factor's actual tangent space.
    /// Shapes are checked for the whole product before slicing any factor.
    pub fn validate_parameter_tangent(
        &self, point: &[f64], parameter: &[f64],
    ) -> Result<(), ProductDifferentialError> {
        self.layout().validate_payload_len(ProductCoordinate::Point, point)?;
        self.layout().validate_payload_len(ProductCoordinate::Parameter, parameter)?;
        for block in self.layout().factors() {
            block.factor().manifold().validate_parameter_tangent(
                &point[point_range(block)], &parameter[parameter_range(block)],
            ).map_err(|e| factor_error(block, e))?;
        }
        Ok(())
    }

    /// Retract every factor at the SAME curve parameter and concatenate its
    /// landing and parameter-coordinate curve velocity. Pairing the velocity
    /// with a parameter gradient gives the derivative along this product curve.
    /// The point has total point length; velocity has total parameter length.
    /// This preserves the exact public factor landings, including SO(3)'s
    /// canonical representative and the production Stiefel QR operation order.
    pub fn retract_curve(
        &self, point: &[f64], direction: &[f64], alpha: f64,
    ) -> Result<RetractionCurve, ProductDifferentialError> {
        self.layout().validate_payload_len(ProductCoordinate::Point, point)?;
        self.layout().validate_payload_len(ProductCoordinate::Parameter, direction)?;
        let mut output = storage(self.layout().point_dim().get() as usize)?;
        let mut velocity = storage(self.layout().param_dim().get() as usize)?;
        for block in self.layout().factors() {
            let curve = block.factor().manifold().retract_curve(
                &point[point_range(block)], &direction[parameter_range(block)], alpha,
            ).map_err(|e| factor_error(block, e))?;
            output.extend_from_slice(&curve.point);
            velocity.extend_from_slice(&curve.velocity);
        }
        Ok(RetractionCurve { point: output, velocity })
    }

    /// Move a parameter vector along the declared blockwise retraction.
    /// Every destination block must equal its authoritative retraction bit for
    /// bit. Near-antipodal Sphere and other factor-local refusals propagate;
    /// SO(3)'s body-coordinate and Stiefel's differentiated-QR transport keep
    /// their original scope. A product of these transports is not generally
    /// an isometry or a Levi-Civita parallel-transport claim.
    pub fn transport_parameter(
        &self, from: &[f64], step: &[f64], to: &[f64], vector: &[f64],
    ) -> Result<Vec<f64>, ProductDifferentialError> {
        self.layout().validate_payload_len(ProductCoordinate::Point, from)?;
        self.layout().validate_payload_len(ProductCoordinate::Parameter, step)?;
        self.layout().validate_payload_len(ProductCoordinate::Point, to)?;
        self.layout().validate_payload_len(ProductCoordinate::Parameter, vector)?;
        let mut output = storage(self.layout().param_dim().get() as usize)?;
        for block in self.layout().factors() {
            let pr = point_range(block);
            let tr = parameter_range(block);
            let transported = block.factor().manifold().transport_parameter(
                &from[pr.clone()], &step[tr.clone()], &to[pr], &vector[tr],
            ).map_err(|e| factor_error(block, e))?;
            output.extend_from_slice(&transported);
        }
        Ok(output)
    }
}
