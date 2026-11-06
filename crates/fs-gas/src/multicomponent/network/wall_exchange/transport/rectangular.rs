//! No-slip rectangular Poiseuille resistance on the original gas-wall owner.
//!
//! Mortensen, Okkels and Bruus, Phys. Rev. E 71, 057301 (2005), equations
//! (3), (10), (11), <https://doi.org/10.1103/PhysRevE.71.057301>.

use super::{Length, Mesh1d, WallExchangeError, positive, same_geometry};
use fs_math::det;

/// Explicit quasi-steady, developed-laminar model of a rectangular gas passage.
///
/// The complete no-slip cross-sectional solution supplies resistance. A
/// circular Poiseuille formula with the hydraulic diameter substituted does
/// not describe this shape. Geometry is prepared once; source intervals still
/// refresh viscosity and density from the original gas state.
///
/// All four walls share the original contact's temperature and axial velocity.
/// The authored Nusselt number uses `Dh = 2 width height / (width + height)`
/// and must apply to this rectangle's aspect ratio and thermal wall conditions.
/// No circular-pipe Nusselt value, turbulent law, entrance development,
/// rarefaction/slip or unsteady transverse velocity profile is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RectangularLaminarWall {
    width: Length,
    height: Length,
    area_m2: f64,
    hydraulic_diameter: Length,
    conductance_factor: f64,
    nusselt_number: f64,
    reynolds_number: [f64; 2],
    maximum_mach_number: f64,
}

impl RectangularLaminarWall {
    /// Prepare the exact rectangular shape factor and explicit applicability.
    ///
    /// Width and height are positive finite lengths. `nusselt_number = h Dh/k`
    /// is finite and nonnegative; zero disables wall heat exchange. The
    /// inclusive Reynolds interval must lie within `[0, 2300]`, with positive
    /// upper limit. This is a screening ceiling, not proof of developed flow.
    /// The wall-relative Mach ceiling is positive and at most 0.3. Geometric
    /// area and hydraulic diameter must remain representable.
    pub fn try_new(
        width: Length,
        height: Length,
        nusselt_number: f64,
        reynolds_number: [f64; 2],
        maximum_mach_number: f64,
    ) -> Result<Self, WallExchangeError> {
        let [minimum, maximum] = reynolds_number;
        if !positive(width.value())
            || !positive(height.value())
            || !nusselt_number.is_finite()
            || nusselt_number < 0.0
            || !minimum.is_finite()
            || minimum < 0.0
            || !positive(maximum)
            || minimum > maximum
            || maximum > 2300.0
            || !positive(maximum_mach_number)
            || maximum_mach_number > 0.3
        {
            return Err(WallExchangeError::Invalid(
                "rectangular wall geometry or applicability",
            ));
        }
        let h = width.value().min(height.value());
        let w = width.value().max(height.value());
        let aspect = h / w;
        let area_m2 = h * w;
        let hydraulic_diameter = Length::new(h * (2.0 / (1.0 + aspect)));
        if !positive(area_m2) || !positive(hydraulic_diameter.value()) {
            return Err(WallExchangeError::Invalid(
                "unrepresentable rectangular wall geometry",
            ));
        }
        Ok(Self {
            width,
            height,
            area_m2,
            hydraulic_diameter,
            conductance_factor: conductance_factor(aspect),
            nusselt_number,
            reynolds_number,
            maximum_mach_number,
        })
    }

    /// Original width [m]. Interchanging width and height preserves the law.
    #[must_use]
    pub const fn width(self) -> Length {
        self.width
    }

    /// Original height [m].
    #[must_use]
    pub const fn height(self) -> Length {
        self.height
    }

    /// Hydraulic diameter [m], used for Reynolds, Nusselt and wetted area.
    #[must_use]
    pub const fn hydraulic_diameter(self) -> Length {
        self.hydraulic_diameter
    }

    /// Authored rectangular-wall `h Dh/k`, including its boundary assumptions.
    #[must_use]
    pub const fn nusselt_number(self) -> f64 {
        self.nusselt_number
    }

    /// Inclusive admitted Reynolds interval based on the hydraulic diameter.
    #[must_use]
    pub const fn reynolds_number(self) -> [f64; 2] {
        self.reynolds_number
    }

    /// Original wall-relative Mach ceiling.
    #[must_use]
    pub const fn maximum_mach_number(self) -> f64 {
        self.maximum_mach_number
    }

    pub(super) fn admit_mesh(self, mesh: &Mesh1d) -> Result<(), WallExchangeError> {
        let Some(area) = mesh.uniform_cross_section() else {
            return Err(WallExchangeError::Invalid(
                "developed rectangular wall requires constant area",
            ));
        };
        if !same_geometry(self.area_m2, area) {
            return Err(WallExchangeError::Invalid(
                "wall width and height differ from original branch area",
            ));
        }
        for (faces, &volume) in mesh.positions().windows(2).zip(mesh.volumes()) {
            if !same_geometry(volume, area * (faces[1] - faces[0])) {
                return Err(WallExchangeError::Invalid(
                    "rectangular wall volume differs from original geometry",
                ));
            }
        }
        Ok(())
    }

    // Q = w h^3 C dp/(12 mu L), hence dU/dt = -12 mu U/(rho h^2 C).
    // The geometry series is never reevaluated in this hot source path.
    pub(super) fn linear_rate(
        self,
        viscosity: f64,
        density: f64,
    ) -> Result<f64, WallExchangeError> {
        let h = self.width.value().min(self.height.value());
        let rate = (12.0 / self.conductance_factor) * (viscosity / density) / h / h;
        if !positive(rate) {
            return Err(WallExchangeError::Invalid(
                "unrepresentable rectangular laminar wall rate",
            ));
        }
        Ok(rate)
    }
}

fn conductance_factor(aspect: f64) -> f64 {
    // C = 1 - 192 a/pi^5 sum_odd(tanh(n pi/(2a))/n^5).
    // Sum the exponential differences from (31/32) zeta(5), avoiding the slow
    // algebraic tail. For 0 < a <= 1, the omitted correction from n=13 is
    // <1e-23. a=0 is the representable parallel-plate limit; exp(-infinity)=0.
    const ODD_ZETA_FIVE: f64 = 1.004_523_762_795_139_6;
    let mut difference = 0.0;
    for n in [1.0_f64, 3.0, 5.0, 7.0, 9.0, 11.0] {
        let e = det::exp(-core::f64::consts::PI * n / aspect);
        difference += (2.0 * e / (1.0 + e)) / n.powi(5);
    }
    1.0 - (192.0 / core::f64::consts::PI.powi(5)) * aspect * (ODD_ZETA_FIVE - difference)
}
