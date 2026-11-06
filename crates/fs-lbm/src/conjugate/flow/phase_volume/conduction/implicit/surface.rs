//! Original exposed-face laws in the same conservative fluid/wall heat solve.

use super::*;
use crate::d3q19::Face3;
use fs_conduction::enthalpy::backward_euler::PrescribedBathHeat;
use fs_conduction::enthalpy::conduction::implicit::boundary::{
    ImplicitBoundaryExchangeLaw, ImplicitBoundaryGeometry, ImplicitBoundaryLaw,
};
use fs_conduction::enthalpy::conduction::implicit::source::ImplicitHeatSourceState;

/// Distinct fixed baths and original source laws on exposed thermal faces.
/// Each physical face has one owner, including across these two slices.
///
/// Source laws include their actual conductive half-cell and all surface
/// mechanisms; the bridge supplies the original graded area and half-width.
/// No second film, radiation term, surface storage or time integrator is added.
/// They join the existing native local-source residual, whose entropy report
/// remains explicitly unclassified. Separate returned bath heats retain their
/// actual temperatures and signs without asserting a complete entropy audit.
#[derive(Clone, Copy, Default)]
pub struct ThermalVolumeBoundaries<'a> {
    pub baths: &'a [ThermalVolumeBathFace],
    pub sources: &'a [ThermalVolumeBoundaryFace<'a>],
}

impl<'a> From<&'a [ThermalVolumeBathFace]> for ThermalVolumeBoundaries<'a> {
    fn from(baths: &'a [ThermalVolumeBathFace]) -> Self {
        Self {
            baths,
            sources: &[],
        }
    }
}

impl<'a, const N: usize> From<&'a [ThermalVolumeBathFace; N]> for ThermalVolumeBoundaries<'a> {
    fn from(baths: &'a [ThermalVolumeBathFace; N]) -> Self {
        Self {
            baths,
            sources: &[],
        }
    }
}

/// A complete original-law boundary on one exposed fluid or finite-wall face.
///
/// The descriptor supplies the selected bath's temperature and fixed areal
/// resistance. Its resistance has the ORIGINAL LAW's meaning: an integrated
/// Robin law may own the external film there; a gray surface law instead owns
/// a coating/contact before the surface, with its convective film supplied
/// separately. Do not put that film into both fields. A two-bath owner binds
/// the second actual temperature itself; neither is an effective bath.
/// Material source selection and external pressure/flow context remain fixed
/// in the borrowed law. Its caloric object must be the original cell's object.
#[derive(Clone, Copy)]
pub enum ThermalVolumeBoundaryFace<'a> {
    Law {
        face: ThermalVolumeBathFace,
        law: &'a dyn ImplicitBoundaryLaw,
    },
    Exchange {
        face: ThermalVolumeBathFace,
        law: &'a dyn ImplicitBoundaryExchangeLaw,
    },
}

/// One accepted reservoir exchange, even when another channel cancels it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalVolumeBathExchange {
    pub cell: usize,
    pub face: Face3,
    /// Original source order 0/1; a single-bath law uses channel 0.
    pub channel: usize,
    pub temperature_k: f64,
    /// Signed physical heat INTO the original fluid or wall inventory [J].
    pub heat_j: f64,
}

#[derive(Clone, Copy)]
enum FaceLaw<'a> {
    Resistance,
    Law(&'a dyn ImplicitBoundaryLaw),
    Exchange(&'a dyn ImplicitBoundaryExchangeLaw),
}

pub(super) struct BoundBath<'a> {
    pub(super) side: usize,
    pub(super) face: Face3,
    temperature: f64,
    area: f64,
    distance: f64,
    resistance: f64,
    law: FaceLaw<'a>,
}

pub(super) struct FaceResponse {
    pub(super) heat: ImplicitHeatSourceState,
    /// Original signed OUTWARD powers. Only `bath_count` entries participate.
    pub(super) baths: [PrescribedBathHeat; 2],
    pub(super) bath_count: usize,
}

impl<'a> BoundBath<'a> {
    fn new(domain: &VoxelDomain, face: ThermalVolumeBathFace, law: FaceLaw<'a>) -> Self {
        let side = face.face as usize;
        Self {
            side,
            face: face.face,
            temperature: face.temperature_k,
            area: domain.face_area(face.cell, side / 2),
            distance: 0.5 * domain.widths(face.cell)[side / 2],
            resistance: face.resistance_m2_k_w,
            law,
        }
    }

    pub(super) fn response(
        &self,
        sample: PhaseConductivitySample,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<FaceResponse, DepositionError> {
        checkpoint()?;
        let geometry = ImplicitBoundaryGeometry {
            area_m2: self.area,
            half_distance_m: self.distance,
            surface_resistance_m2_k_w: self.resistance,
        };
        let (heat, baths, bath_count) = match self.law {
            FaceLaw::Resistance => {
                let heat = self.fixed_response(sample)?;
                let baths = self.single_bath(-heat.power_w);
                (heat, baths, 1)
            }
            FaceLaw::Law(law) => {
                // Query before any insulating-bulk shortcut. The complete
                // source owns its domain and the actual half-cell physics.
                let state = law.state_at_specific_enthalpy(
                    sample.specific_enthalpy_j_kg(),
                    self.temperature,
                    geometry,
                    checkpoint,
                )?;
                let g = state.conductance_w_k;
                if !g.is_finite() || g < 0.0 || !state.power_enthalpy_derivative_w_kg_j.is_finite()
                {
                    return Err(DepositionError::MaterialDomain);
                }
                let delta = sample.temperature_k() - self.temperature;
                let outward = g * delta;
                if !outward.is_finite() || (g > 0.0 && delta != 0.0 && outward == 0.0) {
                    return Err(DepositionError::Arithmetic);
                }
                (
                    ImplicitHeatSourceState {
                        power_w: -outward,
                        power_enthalpy_derivative_w_kg_j: -state.power_enthalpy_derivative_w_kg_j,
                    },
                    self.single_bath(outward),
                    1,
                )
            }
            FaceLaw::Exchange(law) => {
                let state = law.state_at_specific_enthalpy(
                    sample.specific_enthalpy_j_kg(),
                    self.temperature,
                    geometry,
                    checkpoint,
                )?;
                if !state.power_enthalpy_derivative_w_kg_j.is_finite()
                    || state.baths.iter().any(|bath| {
                        !bath.temperature_k.is_finite()
                            || bath.temperature_k <= 0.0
                            || !bath.outward_power_w.is_finite()
                    })
                {
                    return Err(DepositionError::MaterialDomain);
                }
                let outward = state.baths[0].outward_power_w + state.baths[1].outward_power_w;
                if !outward.is_finite() {
                    return Err(DepositionError::Arithmetic);
                }
                (
                    ImplicitHeatSourceState {
                        power_w: -outward,
                        power_enthalpy_derivative_w_kg_j: -state.power_enthalpy_derivative_w_kg_j,
                    },
                    state.baths,
                    2,
                )
            }
        };
        checkpoint()?;
        Ok(FaceResponse {
            heat,
            baths,
            bath_count,
        })
    }

    fn single_bath(&self, outward_power_w: f64) -> [PrescribedBathHeat; 2] {
        [
            PrescribedBathHeat {
                temperature_k: self.temperature,
                outward_power_w,
            },
            PrescribedBathHeat {
                temperature_k: self.temperature,
                outward_power_w: 0.0,
            },
        ]
    }

    // Preserve the original fixed-bath scalar arithmetic and selected tangent.
    fn fixed_response(
        &self,
        sample: PhaseConductivitySample,
    ) -> Result<ImplicitHeatSourceState, DepositionError> {
        let k = sample.conductivity_w_m_k();
        let dk = sample.conductivity_enthalpy_derivative_kg_w_m_k_j();
        if k == 0.0 {
            // A nonzero selected tangent at zero k needs its finite one-sided
            // limit; the original prepared material source refuses that case.
            if dk != 0.0 {
                return Err(DepositionError::MaterialDomain);
            }
            return Ok(ImplicitHeatSourceState {
                power_w: 0.0,
                power_enthalpy_derivative_w_kg_j: 0.0,
            });
        }
        let bulk = self.distance / k;
        if !bulk.is_finite() || bulk <= 0.0 {
            return Err(DepositionError::Arithmetic);
        }
        let largest = bulk.max(self.resistance);
        let scaled_total = 1.0 + bulk.min(self.resistance) / largest;
        let g = (self.area / largest) / scaled_total;
        let delta = self.temperature - sample.temperature_k();
        let dg_term = if dk == 0.0 || delta == 0.0 {
            0.0
        } else {
            g * ((bulk / largest) / scaled_total) * (dk / k) * delta
        };
        let power = g * delta;
        let derivative = dg_term - g * sample.temperature_enthalpy_derivative_kg_k_j();
        if !g.is_finite()
            || g <= 0.0
            || !power.is_finite()
            || !derivative.is_finite()
            || (delta != 0.0 && power == 0.0)
        {
            return Err(DepositionError::Arithmetic);
        }
        Ok(ImplicitHeatSourceState {
            power_w: power,
            power_enthalpy_derivative_w_kg_j: derivative,
        })
    }
}

pub(super) struct BathSource<'a, 'm> {
    pub(super) conductivity: &'a PreparedPhaseConductivity<'m>,
    pub(super) faces: Vec<BoundBath<'a>>,
}

impl ImplicitHeatSourceLaw for BathSource<'_, '_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.conductivity.caloric_curve()
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitHeatSourceState, DepositionError> {
        checkpoint()?;
        let sample = self
            .conductivity
            .state_at_specific_enthalpy(h)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let mut power = Sum::default();
        let mut derivative = Sum::default();
        for face in &self.faces {
            let state = face.response(sample, checkpoint)?.heat;
            power.add(state.power_w);
            derivative.add(state.power_enthalpy_derivative_w_kg_j);
        }
        if !power.value.is_finite() || !derivative.value.is_finite() {
            return Err(DepositionError::Arithmetic);
        }
        checkpoint()?;
        Ok(ImplicitHeatSourceState {
            power_w: power.value,
            power_enthalpy_derivative_w_kg_j: derivative.value,
        })
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn prepare<'a, 'm>(
    domain: &VoxelDomain,
    plain: BTreeMap<(usize, usize), &ThermalVolumeBathFace>,
    selections: &[ThermalVolumeBoundaryFace<'a>],
    samples: &[Option<PhaseConductivitySample>],
    source: impl Fn(usize) -> &'a PreparedPhaseConductivity<'m>,
    max_faces: usize,
    gate: &CancelGate,
) -> Result<BTreeMap<usize, BathSource<'a, 'm>>, ChtError> {
    poll(gate)?;
    if plain
        .len()
        .checked_add(selections.len())
        .is_none_or(|count| count > max_faces)
    {
        return Err(refuse(
            "conductive face-work budget exhausted by thermal boundaries",
        ));
    }
    let mut grouped = BTreeMap::<usize, BathSource<'a, 'm>>::new();
    for (index, ((cell, _), face)) in plain.into_iter().enumerate() {
        poll_every(gate, index)?;
        grouped
            .entry(cell)
            .or_insert_with(|| BathSource {
                conductivity: source(cell),
                faces: Vec::new(),
            })
            .faces
            .push(BoundBath::new(domain, *face, FaceLaw::Resistance));
    }
    for (index, selection) in selections.iter().enumerate() {
        poll_every(gate, index)?;
        let (face, law, caloric) = match *selection {
            ThermalVolumeBoundaryFace::Law { face, law } => {
                (face, FaceLaw::Law(law), law.caloric_law())
            }
            ThermalVolumeBoundaryFace::Exchange { face, law } => {
                (face, FaceLaw::Exchange(law), law.caloric_law())
            }
        };
        let (cell, side) = super::super::boundary::admit_face(domain, &face, |cell| {
            samples.get(cell).is_some_and(Option::is_some)
        })?;
        let group = grouped.entry(cell).or_insert_with(|| BathSource {
            conductivity: source(cell),
            faces: Vec::new(),
        });
        if !std::ptr::addr_eq(
            group.conductivity.caloric_curve() as &dyn ThermalLaw,
            caloric,
        ) {
            return Err(refuse(
                "thermal boundary source requires the original cell caloric object",
            ));
        }
        if group.faces.iter().any(|old| old.side == side) || group.faces.len() >= 6 {
            return Err(refuse(
                "one exposed thermal face cannot own duplicate bath or source laws",
            ));
        }
        group.faces.push(BoundBath::new(domain, face, law));
    }
    // Plain faces already follow the exact original BTreeMap traversal order.
    if !selections.is_empty() {
        for group in grouped.values_mut() {
            poll(gate)?;
            group.faces.sort_by_key(|face| face.side);
        }
    }
    poll(gate)?;
    Ok(grouped)
}
