//! Exact face integration of the original temperature-dependent film sources.

use super::*;
use crate::material::{ConductivityTable, TemperatureSpan, integrate_scalar_path};

#[derive(Debug, Clone, PartialEq)]
pub(super) struct TemperatureSheet {
    pub(super) axes: Option<[[f64; 3]; 2]>,
    tables: Vec<ConductivityTable>,
    thickness_m: f64,
}

impl TemperatureSheet {
    fn action(&self, values: [f64; 2], gradient: [f64; 3]) -> [f64; 3] {
        match self.axes {
            None => gradient.map(|g| values[0] * g),
            Some(axes) => {
                let q = std::array::from_fn::<_, 2, _>(|i| values[i] * dot(axes[i], gradient));
                std::array::from_fn(|i| axes[0][i].mul_add(q[0], axes[1][i] * q[1]))
            }
        }
    }
}

impl SurfaceHeatConduction {
    /// Bind an isotropic in-plane conductivity table [W/(m K)] and fixed
    /// positive film thickness [m]. Original knots, support and receipts are
    /// retained. A declared table gains no material-source authority.
    ///
    /// This is the tangential response of a thin perfectly bonded film; the
    /// caller justifies negligible through-thickness temperature variation.
    /// Storage and normal/contact resistance remain separate physical owners.
    pub fn from_isotropic_conductivity(
        boundary_face: usize,
        conductivity: ConductivityTable,
        thickness_m: f64,
    ) -> Result<Self, ConductionError> {
        Self::from_temperature_tables(boundary_face, None, vec![conductivity], thickness_m)
    }

    /// Bind two original principal in-plane conductivity tables [W/(m K)],
    /// fixed orthonormal world-space tangent axes and positive thickness [m].
    /// Each principal law retains its own temperature dependence and source.
    /// These are intrinsic tangential components: no reduction of a bulk
    /// tensor with normal/tangential coupling or evolving frame is inferred.
    pub fn from_orthotropic_conductivity(
        boundary_face: usize,
        tangent_axes: [[f64; 3]; 2],
        conductivity: [ConductivityTable; 2],
        thickness_m: f64,
    ) -> Result<Self, ConductionError> {
        check_axes(tangent_axes)?;
        Self::from_temperature_tables(
            boundary_face,
            Some(tangent_axes),
            Vec::from(conductivity),
            thickness_m,
        )
    }

    fn from_temperature_tables(
        boundary_face: usize,
        axes: Option<[[f64; 3]; 2]>,
        tables: Vec<ConductivityTable>,
        thickness_m: f64,
    ) -> Result<Self, ConductionError> {
        crate::require_finite("film thickness", thickness_m)?;
        if thickness_m <= 0.0 {
            return Err(invalid("table-backed film thickness must be positive"));
        }
        let mut span = TemperatureSpan::Unbounded;
        for table in &tables {
            span = span.intersect(table.span());
            for &(_, k) in table.knots() {
                sheet_product(k, thickness_m, 1.0)?;
            }
        }
        if matches!(span, TemperatureSpan::Sampled { low, high } if low > high) {
            return Err(invalid(
                "principal film sources have disjoint temperature support",
            ));
        }
        Ok(Self {
            boundary_face,
            law: SheetLaw::Temperature(Arc::new(TemperatureSheet {
                axes,
                tables,
                thickness_m,
            })),
        })
    }

    /// Retained original conductivity sources, in principal-axis order (one
    /// table for isotropy). Existing nodal coefficient laws return `None`.
    pub fn conductivity_tables(&self) -> Option<&[ConductivityTable]> {
        match &self.law {
            SheetLaw::Temperature(law) => Some(&law.tables),
            _ => None,
        }
    }

    /// Whether the sheet actually varies with its trace temperature.
    #[must_use]
    pub fn is_temperature_dependent(&self) -> bool {
        self.conductivity_tables().is_some_and(|tables| {
            tables
                .iter()
                .any(ConductivityTable::is_temperature_dependent)
        })
    }

    /// Common source support. Even a flat sampled table retains its bounds.
    #[must_use]
    pub fn temperature_span(&self) -> TemperatureSpan {
        self.conductivity_tables()
            .map_or(TemperatureSpan::Unbounded, |tables| {
                tables
                    .iter()
                    .fold(TemperatureSpan::Unbounded, |span, table| {
                        span.intersect(table.span())
                    })
            })
    }

    // The constant carrier is only the exact CURRENT face moment for assembly;
    // the original table-backed boundary is never replaced or frozen in time.
    pub(super) fn at_temperature(
        &self,
        cx: Option<&Cx<'_>>,
        nodal: [f64; 3],
        with_derivatives: bool,
    ) -> Result<(Self, [[f64; 2]; 3]), ConductionError> {
        let SheetLaw::Temperature(law) = &self.law else {
            return Ok((self.clone(), [[0.0; 2]; 3]));
        };
        for t in nodal {
            crate::require_finite("film temperature", t)?;
            self.temperature_span().check(t)?;
        }
        let mut values = [0.0; 2];
        let mut derivatives = [[0.0; 2]; 3];
        for (axis, table) in law.tables.iter().enumerate() {
            let (mean, slopes) = integrate_scalar_path(cx, table.knots(), nodal, with_derivatives)?;
            values[axis] = sheet_product(mean, law.thickness_m, 1.0)?;
            for (column, slope) in derivatives.iter_mut().zip(slopes) {
                column[axis] = sheet_product(slope, law.thickness_m, 1.0)?;
            }
        }
        let frozen = match law.axes {
            None => SheetLaw::Isotropic([values[0]; 3]),
            Some(axes) => SheetLaw::Tensor {
                axes,
                nodal: [[[values[0], 0.0], [0.0, values[1]]]; 3],
            },
        };
        Ok((
            Self {
                boundary_face: self.boundary_face,
                law: frozen,
            },
            derivatives,
        ))
    }

    /// Same face graph for the Picard operator or its full temperature tangent.
    /// `J_ac = A g_a.<K_s>g_c + A g_a.<lambda_c K_s'>grad_s(T)`.
    pub(crate) fn matrix_at(
        &self,
        cx: &Cx<'_>,
        mesh: &ConductionMesh,
        temperature: &[f64],
        jacobian: bool,
    ) -> Result<([usize; 3], [[f64; 3]; 3]), ConductionError> {
        let SheetLaw::Temperature(law) = &self.law else {
            return self.matrix(mesh);
        };
        let (vertices, gradients, area) = self.geometry(mesh)?;
        let nodal = vertices.map(|v| temperature[v]);
        let (frozen, derivatives) = self.at_temperature(Some(cx), nodal, jacobian)?;
        let (_, mut block) = frozen.matrix_on(vertices, gradients, area)?;
        if jacobian {
            let gradient = std::array::from_fn(|i| {
                (nodal[1] - nodal[0])
                    .mul_add(gradients[1][i], (nodal[2] - nodal[0]) * gradients[2][i])
            });
            for (c, derivative) in derivatives.into_iter().enumerate() {
                let action = law.action(derivative, gradient);
                let first = area * dot(gradients[1], action);
                let second = area * dot(gradients[2], action);
                // Internal transfer has a zero column sum, including the
                // coefficient derivative; it cannot create a heat source.
                for (row, correction) in block.iter_mut().zip([-first - second, first, second]) {
                    row[c] += correction;
                    crate::require_finite("film conductivity tangent", row[c])?;
                }
            }
        }
        Ok((vertices, block))
    }

    pub(crate) fn is_smooth(
        &self,
        mesh: &ConductionMesh,
        temperature: &[f64],
    ) -> Result<bool, ConductionError> {
        let Some(tables) = self.conductivity_tables() else {
            return Ok(true);
        };
        let (vertices, _, _) = self.geometry(mesh)?;
        let nodal = vertices.map(|v| temperature[v]);
        let span = self.temperature_span();
        for t in nodal {
            crate::require_finite("film temperature", t)?;
            span.check(t)?;
        }
        if let TemperatureSpan::Sampled { low, high } = span {
            if nodal.iter().any(|&t| t == low || t == high) {
                return Ok(false);
            }
        }
        if nodal.iter().any(|&t| t != nodal[0]) {
            return Ok(true);
        }
        for table in tables {
            let slope = table.derivative(nodal[0])?;
            crate::require_finite("film conductivity derivative", slope)?;
            for neighbour in [nodal[0].next_down(), nodal[0].next_up()] {
                match table.derivative(neighbour) {
                    Ok(other) if other == slope => {}
                    Ok(_) | Err(ConductionError::OutsideTemperatureSpan { .. }) => {
                        return Ok(false);
                    }
                    Err(error) => return Err(error),
                }
            }
        }
        Ok(true)
    }
}

impl ThermalBoundary {
    /// Linear-only owners must not silently freeze a temperature-dependent film.
    #[must_use]
    pub fn has_temperature_dependent_surface_conduction(&self) -> bool {
        self.surface_conduction
            .iter()
            .any(SurfaceHeatConduction::is_temperature_dependent)
    }
}
