//! Borrowed original one- and two-sided Maxwell source lines.
//!
//! The temporal observer uses one grouping and contraction engine for both
//! source kinds. Native modal fields are summed before the six spatial FFTs.

use super::*;

#[derive(Clone, Copy)]
pub(super) enum ExposureLines<'b, 's, 'm> {
    OneSided(&'b [PlanarAngularSpectrum<'s, 'm>]),
    TwoSided(&'b [PreparedAngularSpectrumIllumination<'b, 's, 'm>]),
}

#[derive(Clone, Copy)]
pub(super) enum ExposureLine<'b, 's, 'm> {
    OneSided(&'b PlanarAngularSpectrum<'s, 'm>),
    TwoSided(&'b PreparedAngularSpectrumIllumination<'b, 's, 'm>),
}

#[derive(Clone, Copy)]
pub(super) enum ExposureSide {
    OneSided(AngularSpectrumSide),
    TwoSided(AngularSpectrumIlluminationSide),
}

impl ExposureSide {
    fn illumination(self) -> AngularSpectrumIlluminationSide {
        use AngularSpectrumIlluminationSide as I;
        match self {
            Self::TwoSided(side) => side,
            Self::OneSided(side) => match side {
                AngularSpectrumSide::Incident => I::IncomingLeft,
                AngularSpectrumSide::Reflected => I::OutgoingLeft,
                AngularSpectrumSide::TotalIncidentSide => I::TotalLeft,
                AngularSpectrumSide::Transmitted => I::OutgoingRight,
            },
        }
    }
}

impl<'b, 's, 'm> ExposureLines<'b, 's, 'm> {
    pub(super) fn len(self) -> usize {
        match self {
            Self::OneSided(lines) => lines.len(),
            Self::TwoSided(lines) => lines.len(),
        }
    }

    pub(super) fn get(self, index: usize) -> ExposureLine<'b, 's, 'm> {
        match self {
            Self::OneSided(lines) => ExposureLine::OneSided(&lines[index]),
            Self::TwoSided(lines) => ExposureLine::TwoSided(&lines[index]),
        }
    }
}

impl<'b, 's, 'm> ExposureLine<'b, 's, 'm> {
    pub(super) fn beam(self) -> &'b PlanarAngularSpectrum<'s, 'm> {
        match self {
            Self::OneSided(beam) => beam,
            Self::TwoSided(illumination) => illumination.left_beam(),
        }
    }

    pub(super) fn is_two_sided(self) -> bool {
        matches!(self, Self::TwoSided(_))
    }

    pub(super) fn workspace_bytes(self) -> usize {
        match self {
            Self::OneSided(beam) => beam.workspace_bytes(),
            Self::TwoSided(illumination) => illumination.workspace_bytes(),
        }
    }

    pub(super) fn modal_field(
        self,
        side: ExposureSide,
        distance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<[C64; 3]>, Vec<[C64; 3]>), AngularSpectrumError> {
        let side = side.illumination();
        match self {
            Self::TwoSided(illumination) => illumination.modal_field(side, distance, cancelled),
            Self::OneSided(beam) => {
                use AngularSpectrumIlluminationSide as I;
                let side = match side {
                    I::IncomingLeft => AngularSpectrumSide::Incident,
                    I::OutgoingLeft => AngularSpectrumSide::Reflected,
                    I::TotalLeft => AngularSpectrumSide::TotalIncidentSide,
                    I::OutgoingRight | I::TotalRight => AngularSpectrumSide::Transmitted,
                    I::IncomingRight => {
                        // A one-sided line declares no source at the right port.
                        poll(cancelled)?;
                        let n = beam.modes.len();
                        let mut electric = Vec::new();
                        let mut magnetic = Vec::new();
                        electric
                            .try_reserve_exact(n)
                            .map_err(|_| bad("dark incident field allocation failed"))?;
                        magnetic
                            .try_reserve_exact(n)
                            .map_err(|_| bad("dark incident field allocation failed"))?;
                        for index in 0..n {
                            if index % 256 == 0 {
                                poll(cancelled)?;
                            }
                            electric.push([C64::ZERO; 3]);
                            magnetic.push([C64::ZERO; 3]);
                        }
                        return Ok((electric, magnetic));
                    }
                };
                beam.modal_field(side, distance, cancelled)
            }
        }
    }

    pub(super) fn modal_layer_field(
        self,
        cached_one_sided: Option<&PreparedAngularSpectrumLayer<'b, 's, 'm>>,
        layer_index: usize,
        depth: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<[C64; 3]>, Vec<[C64; 3]>), AngularSpectrumError> {
        match self {
            Self::OneSided(beam) => cached_one_sided
                .ok_or_else(|| bad("original one-sided film frequency cache is missing"))?
                .modal_field_for_beam(beam, depth, cancelled),
            Self::TwoSided(illumination) => {
                illumination.modal_layer_field(layer_index, depth, cancelled)
            }
        }
    }
}
