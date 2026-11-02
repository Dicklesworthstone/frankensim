//! One physical wind clock shared by the wall and homogeneous spectral models.

use fs_atmo::evolving::{EvolvingTurbulenceCheckpoint, EvolvingTurbulenceField};
use fs_atmo::mann_field::{HomogeneousMannField, MannFieldCheckpoint};
use fs_atmo::ou::OuCheckpoint;
use fs_atmo::{FlowSample, Refusal};

#[derive(Clone, Debug, PartialEq)]
pub(in crate::simloop) enum WindField {
    Wall(EvolvingTurbulenceField),
    Homogeneous(HomogeneousMannField),
}

#[derive(Clone, Debug, PartialEq)]
pub(in crate::simloop) enum WindCheckpoint {
    Wall(EvolvingTurbulenceCheckpoint),
    Homogeneous(MannFieldCheckpoint),
}

impl WindField {
    pub fn advance_to(&mut self, tick: i64) -> Result<(), Refusal> {
        match self {
            Self::Wall(field) => field.advance_to(tick),
            Self::Homogeneous(field) => field.advance_to(tick),
        }
    }

    #[cfg(test)]
    pub fn tick(&self) -> i64 {
        match self {
            Self::Wall(field) => field.tick(),
            Self::Homogeneous(field) => field.tick(),
        }
    }

    pub fn sample(&self, x: f64, y: f64, h: f64) -> Result<FlowSample, Refusal> {
        match self {
            Self::Wall(field) => field.sample(x, y, h),
            Self::Homogeneous(field) => field.sample([x, y, h]),
        }
    }

    pub fn checkpoint(&self) -> WindCheckpoint {
        match self {
            Self::Wall(field) => WindCheckpoint::Wall(field.checkpoint()),
            Self::Homogeneous(field) => WindCheckpoint::Homogeneous(field.checkpoint()),
        }
    }

    pub fn resume(checkpoint: WindCheckpoint) -> Result<Self, Refusal> {
        match checkpoint {
            WindCheckpoint::Wall(state) => EvolvingTurbulenceField::resume(state).map(Self::Wall),
            WindCheckpoint::Homogeneous(state) => {
                HomogeneousMannField::resume(state).map(Self::Homogeneous)
            }
        }
    }
}

impl WindCheckpoint {
    pub fn seed(&self) -> u64 {
        match self {
            Self::Wall(state) => state.cosine.seed,
            Self::Homogeneous(state) => state.evolution_seed,
        }
    }

    pub fn tick(&self) -> i64 {
        match self {
            Self::Wall(state) => state.cosine.tick,
            Self::Homogeneous(state) => state.paths[0].tick,
        }
    }

    pub fn paths(&self) -> impl Iterator<Item = &OuCheckpoint> {
        let paths = match self {
            Self::Wall(state) => [Some(&state.cosine), Some(&state.sine), None, None],
            Self::Homogeneous(state) => state.paths.each_ref().map(Some),
        };
        paths.into_iter().flatten()
    }

    pub fn paths_mut(&mut self) -> impl Iterator<Item = &mut OuCheckpoint> {
        let paths = match self {
            Self::Wall(state) => [Some(&mut state.cosine), Some(&mut state.sine), None, None],
            Self::Homogeneous(state) => state.paths.each_mut().map(Some),
        };
        paths.into_iter().flatten()
    }
}
