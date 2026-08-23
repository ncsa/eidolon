//! In DNA sequencing, a fragment is a bit of DNA, roughly uniform in length that is sequenced
//! by the machine. Sometimes these fragments have special molecules attached to the end for ID
//! purposes. How this is done is a process called Chemistry Magic. For our purposes, we expect
//! the data to be uniform enough that a mean and standard deviation will describe the set.

use crate::models::lib::{model_reader, model_writer};
use crate::rng::NeatRngError;
use crate::structs::distributions::{DiscreteDistribution, DistributionErrors, NormalDistribution};
use flate2::read::GzDecoder;
use serde::{Deserialize, Serialize};
use serde_json;
use std::io;
use std::path::PathBuf;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FragmentModelError {
    #[error("Fragment model error: {0}")]
    FragModelError(&'static str),
    #[error("Fragment model returned an RNG error: {0}")]
    RngError(#[from] NeatRngError),
    #[error("Fragment model returned an IO error: {0}")]
    IoError(#[from] io::Error),
    #[error("Fragment Model attempted to load a file that it could not find: {0}")]
    FileNotFound(String),
    #[error("Fragment model reported a distribution initiation error: {0}")]
    DistributionInitError(#[from] DistributionErrors),
    #[error("Error building default model!")]
    SerdeError(#[from] serde_json::Error),
    #[error("Error creating fragments of sufficient size. Check fragment length model.")]
    FragGenerationError,
}

#[derive(Debug, Serialize, Deserialize)]
pub enum FragmentLengthModel {
    Discrete {
        distribution: DiscreteDistribution<usize>,
    },
    Normal {
        mean: f64,
        st_dev: f64,
    },
}

static DATA_FILE: &[u8] = include_bytes!("model_data/default_fragment_length_model.json.gz");

impl FragmentLengthModel {
    pub fn new_discrete(
        lengths: Vec<usize>,
        weights: Vec<f64>,
    ) -> Result<Self, FragmentModelError> {
        // These were numbers routinely used for testing in NEAT genReads
        let fragment_dist = DiscreteDistribution::new(&weights, &lengths)?;
        Ok(FragmentLengthModel::Discrete {
            distribution: fragment_dist,
        })
    }

    pub fn new_normal(fragment_mean: f64, fragment_std: f64) -> Result<Self, FragmentModelError> {
        Ok(FragmentLengthModel::Normal {
            mean: fragment_mean,
            st_dev: fragment_std,
        })
    }

    // Returns Result because it deserializes an embedded model file; std::Default
    // requires infallible `fn default() -> Self`, which doesn't fit.
    #[allow(clippy::should_implement_trait)]
    pub fn default() -> Result<Self, FragmentModelError> {
        // The parameters of the default model from the original neat
        // The lengths range from 1 to 799, though it skips from 1 to 32 before counting up.
        // The weights are just numbers between 0 and 1.
        // These are data gathered from publicly availble human data, but should reflect
        // whatever chemistry was used at the time
        let reader = GzDecoder::new(DATA_FILE);
        let data: FragmentLengthModel =
            serde_json::from_reader(reader).map_err(FragmentModelError::SerdeError)?;
        Ok(data)
    }

    pub fn default_normal() -> Result<Self, FragmentModelError> {
        // These were numbers routinely used for testing in NEAT genReads
        Ok(FragmentLengthModel::Normal {
            mean: 300.0,
            st_dev: 30.0,
        })
    }

    pub fn discrete_from_file(filename: &PathBuf) -> Result<Self, FragmentModelError> {
        // The baseline model is really just a mathematical equation and can be reconstructed by the two input parameters
        // But reading from data, it may be better to use the discrete distribution to maintain outliers. This will load such
        // a model from a json file. The file must be of the format:
        // {
        //   "fragment_lengths": [ ... ],
        //   "fragment_weights": [ ... ]
        // }
        // Where fragment lengths are of type usize and fragment weights are of (or can be cast as) type f64.
        if !filename.exists() {
            return Err(FragmentModelError::FileNotFound(
                filename.display().to_string(),
            ));
        }

        let data: FragmentLengthModel = model_reader(filename).unwrap();

        Ok(data)
    }

    /// Mean fragment length, exact and deterministic (no RNG draw) for either
    /// variant -- `Normal`'s own parameter, or `Discrete`'s weighted average
    /// over its value/weight table. Used to scale a boundary-extension
    /// dilution correction to the model actually in use, including a trained
    /// (real-BAM-derived) discrete model, not just the synthetic Normal case.
    pub fn mean(&self) -> Result<f64, FragmentModelError> {
        match self {
            FragmentLengthModel::Normal { mean, .. } => Ok(*mean),
            FragmentLengthModel::Discrete { distribution } => {
                let values = distribution.values()?;
                let weights = distribution.weights()?;
                let total_weight: f64 = weights.iter().sum();
                if total_weight <= 0.0 {
                    return Err(FragmentModelError::FragModelError(
                        "Discrete fragment model has zero total weight",
                    ));
                }
                let weighted_sum: f64 = values
                    .iter()
                    .zip(weights.iter())
                    .map(|(&v, &w)| v as f64 * w)
                    .sum();
                Ok(weighted_sum / total_weight)
            }
        }
    }

    pub fn normal_params(&self) -> Result<(f64, f64), FragmentModelError> {
        // This returns the parameters used to initiate the model
        match self {
            FragmentLengthModel::Discrete { distribution: _ } => {
                Err(FragmentModelError::FragModelError(
                    "Called normal_params on a discrete fragment model",
                ))
            }
            FragmentLengthModel::Normal { mean, st_dev } => Ok((*mean, *st_dev)),
        }
    }

    pub fn generate_fragment(&self, rand: f64) -> Result<usize, FragmentModelError> {
        // This function generates a fragment length based on mean and standard deviation,
        // or based on a discrete distribution.
        match self {
            // The discrete one is pretty easy
            Self::Discrete { distribution } => Ok(distribution.sample(rand)?),
            // for normal we have to build the distribution, then sample. Not sure if this will be a pinch point.
            Self::Normal { mean, st_dev } => {
                let distribution = NormalDistribution::new(*mean, *st_dev)?;
                Ok(distribution.sample(rand)?.trunc() as usize)
            }
        }
    }

    pub fn write_file(&self, filename: &PathBuf) -> Result<(), FragmentModelError> {
        // serialize a model with serde and write it to file
        model_writer(self, filename)?;
        Ok(())
    }

    #[allow(dead_code)]
    fn is_discrete(&self) -> bool {
        // This is mainly for testing purposes
        if let FragmentLengthModel::Discrete { distribution } = self {
            !distribution.weights.is_empty()
        } else {
            false
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normal_default() {
        let model = FragmentLengthModel::default_normal().unwrap();
        match model {
            FragmentLengthModel::Normal { mean, st_dev } => {
                assert_eq!(mean, 300.0);
                assert_eq!(st_dev, 30.0);
            }
            _ => panic!("Wrong type!!"),
        }
    }

    #[test]
    fn test_discrete_default() {
        let model = FragmentLengthModel::default().unwrap();
        match model {
            FragmentLengthModel::Discrete { distribution } => {
                assert_eq!(
                    distribution.values().unwrap(),
                    [
                        1, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
                        50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68,
                        69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87,
                        88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104,
                        105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119,
                        120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134,
                        135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149,
                        150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164,
                        165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179,
                        180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194,
                        195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209,
                        210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224,
                        225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
                        240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254,
                        255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269,
                        270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284,
                        285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299,
                        300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314,
                        315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329,
                        330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344,
                        345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359,
                        360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374,
                        375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389,
                        390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404,
                        405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419,
                        420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434,
                        435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449,
                        450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464,
                        465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479,
                        480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494,
                        495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509,
                        510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524,
                        525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539,
                        540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554,
                        555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569,
                        570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584,
                        585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599,
                        600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614,
                        615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629,
                        630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644,
                        645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659,
                        660, 661, 662, 663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 673, 674,
                        675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689,
                        690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704,
                        705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719,
                        720, 721, 722, 723, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 734,
                        735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749,
                        750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 761, 762, 763, 764,
                        765, 766, 767, 768, 769, 770, 771, 772, 773, 774, 775, 776, 777, 778, 779,
                        780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 790, 791, 793, 795, 796,
                        797, 799
                    ]
                )
            }
            _ => panic!("Wrong type!!"),
        }
    }

    #[test]
    fn test_new_from_mean() {
        let mean = 34.33;
        let std_dev = 1.232;
        let model = FragmentLengthModel::new_normal(mean, std_dev).unwrap();
        match model {
            FragmentLengthModel::Normal { mean, st_dev } => {
                assert_eq!(mean, 34.33);
                assert_eq!(st_dev, 1.232);
            }
            _ => panic!("Wrong type!!"),
        }
    }

    #[test]
    fn test_new_discrete() {
        let l_vec = vec![1, 8, 9, 10];
        let w_vec = vec![1.0, 3.0, 2.0, 1.2];
        let model = FragmentLengthModel::new_discrete(l_vec.clone(), w_vec.clone()).unwrap();
        match model {
            FragmentLengthModel::Discrete { distribution } => {
                assert_eq!(distribution.values().unwrap(), l_vec);
                assert_eq!(
                    distribution.weights().unwrap(),
                    [
                        0.1388888888888889,
                        0.5555555555555556,
                        0.8333333333333334,
                        1.0
                    ]
                );
            }
            _ => panic!("Wrong type!!"),
        }
    }

    #[test]
    fn test_new_normal() {
        let mean = 34.33;
        let std_dev = 1.232;
        let model = FragmentLengthModel::new_normal(mean, std_dev).unwrap();
        match model {
            FragmentLengthModel::Normal { mean, st_dev } => {
                assert_eq!((mean, st_dev), (34.33, 1.232));
            }
            _ => panic!("Wrong type!!"),
        }
    }

    #[test]
    fn test_generate_fragment() {
        let model = FragmentLengthModel::default().unwrap();
        let frag = model.generate_fragment(0.1).unwrap();
        assert_eq!(frag, 295)
    }

    #[test]
    fn test_from_file() {
        let model = FragmentLengthModel::default().unwrap();
        let temp_dir = tempfile::tempdir().unwrap();
        let mut temp_file = PathBuf::from(temp_dir.path());
        let filename = "model_test.json";
        temp_file.push(filename);
        model
            .write_file(&temp_file)
            .expect("write_file should succeed");
        let model2 = FragmentLengthModel::discrete_from_file(&temp_file).unwrap();
        assert_eq!(model.is_discrete(), model2.is_discrete());
        temp_dir.close().unwrap();
    }
}
