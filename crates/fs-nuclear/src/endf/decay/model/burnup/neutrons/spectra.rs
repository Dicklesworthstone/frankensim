//! Source-law dispatch into the original spectral moment owners.
//!
//! LF5 uses the probability-coordinate change E = theta*x. Its density is
//! g(x)/theta, so counts are unchanged and energy moments scale by theta.
//! Original g, theta, U and source grids stay in the parsed evaluation; the
//! transformed intervals are streamed through the same LF1/MF8 integrator.

use crate::endf::decay::model::emissions::{
    DecayEmissionLimits, integrate_continuum, integrate_continuum_segments,
};
use crate::endf::decay::spectra::{DecayRadiation, DecaySpectrumPoint};
use crate::endf::spectrum::PositiveIntegral;
use crate::endf::spontaneous_neutrons::analytic::{self, AnalyticFissionSpectrum};
use crate::endf::spontaneous_neutrons::{
    SpontaneousFissionShapePoint, SpontaneousFissionSpectrumComponent,
    SpontaneousFissionSpectrumLaw,
};
use crate::endf::{EndfError, convert_si};
use crate::{FissionError, positive};

pub(super) fn integrate_component(
    component: &SpontaneousFissionSpectrumComponent,
    edges: &[f64],
    counts: &mut [PositiveIntegral],
    energies: &mut [PositiveIntegral],
    work: &mut usize,
    limits: DecayEmissionLimits,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(), EndfError> {
    match &component.law {
        SpontaneousFissionSpectrumLaw::Tabulated {
            points,
            interpolation,
        } => integrate_continuum(
            DecayRadiation::Neutron,
            points,
            interpolation,
            1.0,
            edges,
            counts,
            energies,
            work,
            limits,
            checkpoint,
        ),
        SpontaneousFissionSpectrumLaw::GeneralEvaporation {
            temperature_j,
            shape,
            interpolation,
        } => {
            let maximum = outgoing_limit(component)?;
            let segments = shape.windows(2).map(|pair| {
                let points = [
                    scaled_point(pair[0], *temperature_j)?,
                    scaled_point(pair[1], *temperature_j)?,
                ];
                if pair[0].coordinate < pair[1].coordinate
                    && points[0].energy_j >= points[1].energy_j
                {
                    return Err(FissionError::Arithmetic.into());
                }
                Ok(points)
            });
            integrate_continuum_segments(
                DecayRadiation::Neutron,
                segments,
                interpolation,
                1.0,
                maximum,
                edges,
                counts,
                energies,
                work,
                limits,
                checkpoint,
            )
        }
        SpontaneousFissionSpectrumLaw::Maxwellian { temperature_j }
        | SpontaneousFissionSpectrumLaw::Evaporation { temperature_j } => {
            let shape = if matches!(
                &component.law,
                SpontaneousFissionSpectrumLaw::Maxwellian { .. }
            ) {
                AnalyticFissionSpectrum::Maxwellian
            } else {
                AnalyticFissionSpectrum::Evaporation
            };
            analytic::integrate(
                shape,
                *temperature_j,
                outgoing_limit(component)?,
                edges,
                counts,
                energies,
                work,
                limits,
                checkpoint,
            )
        }
    }
}

fn outgoing_limit(component: &SpontaneousFissionSpectrumComponent) -> Result<f64, EndfError> {
    let restriction = component
        .restriction_energy_j
        .ok_or(FissionError::InvalidInput(
            "evaporation spectra require their original restriction energy",
        ))?;
    let [reference] = component.incident_energies_j.as_ref() else {
        return Err(EndfError::Unsupported(
            "spontaneous evaporation requires one explicitly supplied reference-energy row",
        ));
    };
    let maximum = reference - restriction;
    if !positive(maximum) {
        return Err(FissionError::InvalidInput(
            "finite positive spontaneous-neutron outgoing cutoff required",
        )
        .into());
    }
    Ok(maximum)
}

fn scaled_point(
    point: SpontaneousFissionShapePoint,
    temperature_j: f64,
) -> Result<DecaySpectrumPoint, EndfError> {
    let energy_j = convert_si(point.coordinate, temperature_j)?;
    let density_per_j = point.density / temperature_j;
    if !density_per_j.is_finite() || (point.density > 0.0 && density_per_j == 0.0) {
        return Err(FissionError::Arithmetic.into());
    }
    Ok(DecaySpectrumPoint {
        energy_j,
        density_per_j,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::endf::neutron_yield::FissionNeutronKind;
    use crate::endf::spontaneous_neutrons::{
        EndfSpontaneousNeutronSelection, SpontaneousFissionNeutronSpectrum,
        SpontaneousFissionNeutronYield, SpontaneousNeutronReadLimits,
    };
    use crate::irradiation::NuclideId;

    #[test]
    fn g2_cf252_evaluated_family_spectra_retain_original_counts_and_energy() {
        // Unchanged MF1/MF5 records extracted from the primary ENDF/B-VIII.0
        // decay_3644_98-Cf-252.zip archive. Other files are outside this fixture.
        // https://www-nds.iaea.org/public/download-endf/ENDF-B-VIII.0/decay/decay_3644_98-Cf-252.zip
        // Expected moments below were independently integrated from original
        // decimal histogram records at 50-digit precision, before SI scaling.
        let text =
            include_str!("../../../../spontaneous_neutrons/tests/cf252_endfb_viii_0_neutrons.endf");
        let selection = EndfSpontaneousNeutronSelection {
            material_number: 3644,
            nuclide: NuclideId::new(98, 252, 0).unwrap(),
        };
        let limits = SpontaneousNeutronReadLimits::default();
        let prompt = SpontaneousFissionNeutronSpectrum::parse(
            "ENDF-B-VIII.0-Cf252-prompt".into(),
            selection,
            FissionNeutronKind::Prompt,
            text,
            limits,
            || Ok(()),
        )
        .unwrap();
        assert_eq!(prompt.components().len(), 1);
        assert!(matches!(
            prompt.components()[0].law,
            SpontaneousFissionSpectrumLaw::Tabulated { .. }
        ));
        let delayed = SpontaneousFissionNeutronSpectrum::parse(
            "ENDF-B-VIII.0-Cf252-delayed".into(),
            selection,
            FissionNeutronKind::Delayed,
            text,
            limits,
            || Ok(()),
        )
        .unwrap();
        let yield_ = SpontaneousFissionNeutronYield::parse(
            "ENDF-B-VIII.0-Cf252-delayed-yield".into(),
            selection,
            FissionNeutronKind::Delayed,
            text,
            limits,
            || Ok(()),
        )
        .unwrap();
        assert_eq!(delayed.components().len(), 6);
        assert_eq!(yield_.decay_constants_s_inv().len(), 6);
        let expected = [
            (1.0000000120373, 525881.5257291035),
            (1.000000035537044, 587883.4445469742),
            (1.000000006764, 380927.25689426),
            (1.00000000646, 449114.1779536),
            (1.00000000055, 421187.90633315),
            (1.00000001174, 493750.5845413),
        ];
        const EV_J: f64 = 1.602_176_634e-19;
        let edges = [0.0, 0.5e6 * EV_J, 1.5e6 * EV_J, 3.0e6 * EV_J, 30.0e6 * EV_J];
        let mut work = 0;
        for (component, (number, energy_ev)) in delayed.components().iter().zip(expected) {
            let mut counts = [PositiveIntegral::new(); 4];
            let mut energies = counts;
            integrate_component(
                component,
                &edges,
                &mut counts,
                &mut energies,
                &mut work,
                DecayEmissionLimits::default(),
                &mut || Ok(()),
            )
            .unwrap();
            let actual_count: f64 = counts.into_iter().map(|x| x.value().unwrap()).sum();
            let actual_energy: f64 = energies
                .into_iter()
                .map(|x| x.value().unwrap() / EV_J)
                .sum();
            assert!(
                (actual_count - number).abs() < 3e-14,
                "{actual_count} != {number}"
            );
            assert!(
                (actual_energy - energy_ev).abs() < 2e-8,
                "{actual_energy} != {energy_ev}"
            );
            assert!(
                (actual_count - 1.0).abs() > 1e-10,
                "source integral was repaired"
            );
        }
        let probability: f64 = delayed.components().iter().map(|c| c.probability).sum();
        assert!((probability - 1.00000014).abs() < 1e-15);
    }
}
