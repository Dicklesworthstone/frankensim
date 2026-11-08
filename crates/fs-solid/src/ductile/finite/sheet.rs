//! Ductile initiation at every original finite-sheet material point.
//!
//! The native sheet solver still owns transverse traction balance, thickness,
//! plastic/creep history, constitutive kinetics and Newton iteration. This gate
//! reads the complete returned Piola/deformation pair once, so triaxiality and
//! Lode dependence see the actual spatial stress, including the solved thickness.
//! Neither element-average stress nor an assumed isochoric total thickness is
//! substituted. Rejected endpoints leave both mechanics and initiation intact.
//!
//! This supplies a local calibrated initiation limit, not a necking/forming-limit
//! diagram, progressive softening, crack topology or post-failure continuation.
//! Endpoint-held integration retains the original criterion's refinement limits.

use fs_material::plastic::finite::FiniteJ2;
use fs_material::plastic::finite::combined::FiniteNortonJ2;
use fs_material::plastic::finite::source::ResolvedFiniteJ2;
use fs_material::plastic::finite::tabulated::source::ResolvedTabulatedFiniteJ2;
use fs_material::plastic::finite::plane_stress::FinitePlaneStressMaterial;
use fs_material::plastic::finite::rate::FiniteOverstressJ2;

use super::*;
use crate::hyper2d::plastic::thermal::{
    FiniteThermoplasticSheet, FiniteThermoplasticSheetError, FiniteThermoplasticSheetStep,
};
use crate::hyper2d::plastic::{
    FinitePlasticSheet, FinitePlasticSheetError, FinitePlasticSheetStep,
};
use crate::hyper3::plastic::FinitePlasticLawField;

/// Sheet laws whose native `alpha` is accumulated equivalent PLASTIC strain.
///
/// Implementing this admission promises that independent creep is excluded
/// from `alpha`. Pure Norton creep deliberately has no implementation: its
/// scalar accumulation is not a calibrated ductile-plastic history. Thermal
/// families use the separately checked thermoplastic entry point below.
pub trait DuctileSheetMaterial: FinitePlaneStressMaterial {}

impl DuctileSheetMaterial for FiniteJ2 {}
impl DuctileSheetMaterial for ResolvedFiniteJ2 {}
impl DuctileSheetMaterial for ResolvedTabulatedFiniteJ2 {}
impl<M: DuctileSheetMaterial + FiniteJ2Material + ?Sized> DuctileSheetMaterial for &M {}
impl DuctileSheetMaterial for FiniteOverstressJ2 {}
impl DuctileSheetMaterial for FiniteNortonJ2 {}

fn check_duration<M: DuctileSheetMaterial>(
    laws: FinitePlasticLawField<'_, M>,
    duration_s: f64,
    cx: &Cx<'_>,
) -> Result<(), DuctileSolidError> {
    let laws = match laws {
        FinitePlasticLawField::Uniform(law) => core::slice::from_ref(law),
        FinitePlasticLawField::PerElement(laws) => laws,
    };
    for (element, law) in laws.iter().enumerate() {
        if element % 256 == 0 {
            checkpoint(cx)?;
        }
        if law.duration_s().is_some_and(|native| native != duration_s) {
            return Err(DuctileSolidError::Input(
                "ductility duration must equal the native sheet constitutive interval",
            ));
        }
    }
    Ok(())
}

impl<'a, L: DuctileInitiationLaw> DuctileHistory<'a, L> {
    /// Advance one original P1/Q1 sheet endpoint, retaining every quadrature
    /// history and its full three-dimensional stress. For kinetic materials,
    /// the supplied duration must be the actual native return-map interval.
    /// Combined Norton/J2 laws consume only plastic `alpha`, never `creep_alpha`.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_finite_sheet<'m, M: DuctileSheetMaterial>(
        &self,
        sheet: &FinitePlasticSheet<'m, M>,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        conditions: DuctileStepConditions<'_>,
        settings: IncrementSettings,
        controls: FinitePlasticStepControls,
        cx: &Cx<'_>,
    ) -> Result<
        DuctileAccepted<'a, FinitePlasticSheetStep<'m, M>, L>,
        DuctileAdvanceError<FinitePlasticSheetError>,
    > {
        if matches!(conditions.temperatures_k, DuctileTemperatures::PerHexElement(_)) {
            return Err(DuctileSolidError::Input(
                "sheet temperatures must be uniform or cover original quadrature points",
            ).into());
        }
        self.prepare(sheet.quadrature_states(), conditions, cx)?;
        check_duration(sheet.laws(), conditions.duration_s, cx)?;
        let mechanics = sheet
            .advance(prescribed_m, nodal_forces_n, settings, controls, cx)
            .map_err(DuctileAdvanceError::Mechanics)?;
        let next = &mechanics.sheet;
        let ductility = self.trial_finite(
            next.quadrature_states(),
            next.quadrature_first_piola_pa(),
            next.quadrature_deformation(),
            conditions,
            cx,
        )?;
        Ok(DuctileAccepted { mechanics, ductility })
    }

    /// Advance one thermoplastic sheet interval with actual endpoint
    /// temperatures in original element/quadrature order. P1 contributes one
    /// point per triangle and Q1 four per quadrilateral; an element mean is not
    /// a substitute for these temperatures. The same physical duration drives
    /// native rate/creep evolution and the ductile initiation law.
    ///
    /// Pure-creep families refuse before any solve; combined creep/plastic
    /// families retain both histories and spend only their plastic accumulation.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_thermoplastic_sheet<'m, M: FiniteThermalJ2Material>(
        &self,
        sheet: &FiniteThermoplasticSheet<'m, M>,
        quadrature_temperatures_k: &[f64],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        duration_s: f64,
        settings: IncrementSettings,
        controls: FinitePlasticStepControls,
        cx: &Cx<'_>,
    ) -> Result<
        DuctileAccepted<'a, FiniteThermoplasticSheetStep<'m, M>, L>,
        DuctileAdvanceError<FiniteThermoplasticSheetError>,
    > {
        let families = match sheet.families() {
            FinitePlasticLawField::Uniform(family) => core::slice::from_ref(family),
            FinitePlasticLawField::PerElement(families) => families,
        };
        check_thermal_families(families, cx)?;
        let conditions = DuctileStepConditions {
            duration_s,
            temperatures_k: DuctileTemperatures::PerPoint(quadrature_temperatures_k),
        };
        self.prepare(sheet.quadrature_states(), conditions, cx)?;
        let mechanics = sheet
            .advance_with_quadrature_temperatures(
                quadrature_temperatures_k,
                prescribed_m,
                nodal_forces_n,
                duration_s,
                settings,
                controls,
                cx,
            )
            .map_err(DuctileAdvanceError::Mechanics)?;
        let next = &mechanics.sheet;
        let ductility = self.trial_finite(
            next.quadrature_states(),
            next.quadrature_first_piola_pa(),
            next.quadrature_deformation(),
            DuctileStepConditions {
                duration_s,
                temperatures_k: DuctileTemperatures::PerPoint(next.quadrature_temperatures_k()),
            },
            cx,
        )?;
        Ok(DuctileAccepted { mechanics, ductility })
    }
}
