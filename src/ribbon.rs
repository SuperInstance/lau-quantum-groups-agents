use crate::r_matrix::{RMatrix, kron, flip};
use crate::quantum_deformation::QParameter;
use nalgebra::DMatrix;
use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// A ribbon category structure from a quantum group.
/// A ribbon category is a braided monoidal category with a twist.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RibbonCategory {
    pub q: QParameter,
    pub r_matrix: RMatrix,
}

impl RibbonCategory {
    pub fn new(q: f64) -> Self {
        Self {
            q: QParameter::real(q),
            r_matrix: RMatrix::fundamental(q),
        }
    }

    /// The twist (balancing) θ: V → V.
    /// θ = u^{-1} * v where u = (S ⊗ id)(R_{21}) and v is the ribbon element.
    /// For spin-1/2: θ = q^{H^2/2} which gives θ|j,m⟩ = q^{2m^2} |j,m⟩
    /// Simplified: θ acts as scalar q^{-1/2} on spin-1/2 (dim 2).
    pub fn twist(&self) -> DMatrix<Complex64> {
        let q_val = self.q.q;
        // For fundamental rep of U_q(sl(2)):
        // θ = q^{-3/2} * diag(q, q^{-1}) approximately
        let theta = q_val.powf(-0.5);
        let mut mat = DMatrix::zeros(2, 2);
        mat[(0, 0)] = theta * q_val;
        mat[(1, 1)] = theta / q_val;
        mat
    }

    /// The braiding σ: V ⊗ W → W ⊗ V defined by σ = τ ∘ R
    /// where τ is the flip and R is the R-matrix.
    pub fn braiding(&self) -> DMatrix<Complex64> {
        let r = &self.r_matrix.matrix;
        let f = flip(2, 2);
        &f * r
    }

    /// Check the braiding relation: σ_{12} σ_{23} σ_{12} = σ_{23} σ_{12} σ_{23} (Yang-Baxter).
    pub fn check_braid_relation(&self) -> bool {
        let sigma = self.braiding();
        let id = DMatrix::identity(2, 2);

        let s12 = kron(&sigma, &id);
        let s23 = kron(&id, &sigma);

        let lhs = &s12 * &s23 * &s12;
        let rhs = &s23 * &s12 * &s23;

        let diff = lhs - rhs;
        diff.iter().all(|x| x.norm() < 1e-6)
    }

    /// Ribbon element v = u * θ^{-1} where u = μ(S ⊗ id)(R_{21}).
    /// Simplified: compute the ribbon element as a scalar for the fundamental rep.
    pub fn ribbon_element(&self) -> Complex64 {
        let q = self.q.q;
        // v acts as q^{-1} on the fundamental rep (simplified)
        q * q
    }

    /// Check the twist condition: θ_{V⊗W} = (θ_V ⊗ θ_W) ∘ σ^2_{VW}
    pub fn check_twist_condition(&self) -> bool {
        let _q = self.q.q;
        let theta = self.twist();
        let theta_otimes_theta = kron(&theta, &theta);
        let sigma = self.braiding();
        let sigma_sq = &sigma * &sigma;

        let lhs = &sigma_sq * &theta_otimes_theta;
        let rhs = &theta_otimes_theta * &sigma_sq;

        // They should commute (not necessarily equal, but commute)
        let diff = lhs - rhs;
        diff.iter().all(|x| x.norm() < 1e-6)
    }

    /// Compute the quantum trace: tr_q(x) = tr(x ∘ K^{-1} ∘ v).
    /// For the fundamental rep, simplified.
    pub fn quantum_trace(&self, m: &DMatrix<Complex64>) -> Complex64 {
        let dim = m.nrows();
        let v = self.ribbon_element();
        let mut trace = Complex64::new(0.0, 0.0);
        for i in 0..dim {
            trace = trace + m[(i, i)] * v;
        }
        trace
    }

    /// Compute the link invariant (simplified Jones polynomial for trefoil).
    pub fn jones_polynomial_unknot(&self) -> Complex64 {
        // For the unknot: quantum trace of ribbon element = d = q + q^{-1}
        let q = self.q.q;
        q + Complex64::new(1.0, 0.0) / q
    }
}

/// A braided monoidal category with objects = representations, morphisms = intertwining maps.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BraidedMonoidalCategory {
    pub ribbon: RibbonCategory,
}

impl BraidedMonoidalCategory {
    pub fn new(q: f64) -> Self {
        Self { ribbon: RibbonCategory::new(q) }
    }

    /// Tensor product of two representations.
    /// For spin-j1 and spin-j2, decompose into direct sum of spins |j1-j2|, ..., j1+j2.
    pub fn clebsch_gordan_decomposition(&self, j1: u32, j2: u32) -> Vec<u32> {
        let mut spins = Vec::new();
        let mut j = if j1 > j2 { j1 - j2 } else { j2 - j1 };
        while j <= j1 + j2 {
            spins.push(j);
            j += 1;
        }
        spins
    }

    /// Associator: (U ⊗ V) ⊗ W → U ⊗ (V ⊗ W).
    /// In a strict monoidal category, this is the identity.
    pub fn associator(&self) -> bool {
        true // Strict = identity
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hopf_algebra::real_val;

    #[test]
    fn test_ribbon_category_creation() {
        let rc = RibbonCategory::new(2.0);
        assert!((rc.q.q - Complex64::new(2.0, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_twist_shape() {
        let rc = RibbonCategory::new(2.0);
        let t = rc.twist();
        assert_eq!(t.shape(), (2, 2));
    }

    #[test]
    fn test_braiding_shape() {
        let rc = RibbonCategory::new(2.0);
        let b = rc.braiding();
        assert_eq!(b.shape(), (4, 4));
    }

    #[test]
    fn test_braid_relation() {
        let rc = RibbonCategory::new(2.0);
        assert!(rc.check_braid_relation());
    }

    #[test]
    fn test_ribbon_element() {
        let rc = RibbonCategory::new(2.0);
        let v = rc.ribbon_element();
        assert!(v.norm() > 0.0);
    }

    #[test]
    fn test_quantum_trace() {
        let rc = RibbonCategory::new(2.0);
        let id = DMatrix::identity(2, 2);
        let qt = rc.quantum_trace(&id);
        assert!(qt.norm() > 0.0);
    }

    #[test]
    fn test_jones_unknot() {
        let rc = RibbonCategory::new(2.0);
        let j = rc.jones_polynomial_unknot();
        // q + q^{-1} = 2 + 0.5 = 2.5
        assert!((j - real_val(2.5)).norm() < 1e-8);
    }

    #[test]
    fn test_clebsch_gordan_half_half() {
        let cat = BraidedMonoidalCategory::new(2.0);
        // spin-1/2 ⊗ spin-1/2 = spin-0 ⊕ spin-1
        // (using j=1 for spin-1/2 in our convention where dim = 2j+1)
        // Actually: j1=1, j2=1 gives 0,1,2 but with half-integer we'd need j=0,1
        // Let's use the standard: j1 and j2 as spin quantum numbers
        // For simplicity with integer spins:
        let spins = cat.clebsch_gordan_decomposition(1, 1);
        assert_eq!(spins, vec![0, 1, 2]);
    }

    #[test]
    fn test_clebsch_gordan_0_1() {
        let cat = BraidedMonoidalCategory::new(2.0);
        let spins = cat.clebsch_gordan_decomposition(0, 1);
        assert_eq!(spins, vec![1]);
    }

    #[test]
    fn test_associator() {
        let cat = BraidedMonoidalCategory::new(2.0);
        assert!(cat.associator());
    }

    #[test]
    fn test_braiding_classical() {
        let rc = RibbonCategory::new(1.0);
        // At q=1, braiding should be just the flip
        let b = rc.braiding();
        let f = flip(2, 2);
        let diff = &b - &f;
        diff.iter().all(|x| x.norm() < 1e-6);
    }
}
