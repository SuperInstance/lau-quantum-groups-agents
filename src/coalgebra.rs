use crate::hopf_algebra::{AlgebraElement, Field};
use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// A coalgebra: vector space with comultiplication Δ and counit ε.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Coalgebra {
    pub dim: usize,
    /// Comultiplication: Δ(e_i) = sum_{j,k} comult[i][j][k] * e_j ⊗ e_k
    pub comult: Vec<Vec<Vec<Field>>>,
    /// Counit: ε(e_i) = counit[i]
    pub counit: Vec<Field>,
}

impl Coalgebra {
    pub fn new(dim: usize, comult: Vec<Vec<Vec<Field>>>, counit: Vec<Field>) -> Self {
        assert_eq!(comult.len(), dim);
        assert_eq!(counit.len(), dim);
        Self { dim, comult, counit }
    }

    /// Trivial 1-dimensional coalgebra.
    pub fn trivial() -> Self {
        let mut comult = vec![vec![vec![Complex64::new(0.0, 0.0); 1]; 1]; 1];
        comult[0][0][0] = Complex64::new(1.0, 0.0);
        Self {
            dim: 1,
            comult,
            counit: vec![Complex64::new(1.0, 0.0)],
        }
    }

    /// Apply comultiplication to an element.
    pub fn comultiply(&self, a: &AlgebraElement) -> Vec<(Field, usize, usize)> {
        let mut result = Vec::new();
        for i in 0..self.dim {
            let ci = a.coeffs[i];
            if ci.norm() < 1e-15 {
                continue;
            }
            for j in 0..self.dim {
                for k in 0..self.dim {
                    let val = ci * self.comult[i][j][k];
                    if val.norm() > 1e-15 {
                        result.push((val, j, k));
                    }
                }
            }
        }
        result
    }

    /// Apply counit to an element.
    pub fn apply_counit(&self, a: &AlgebraElement) -> Field {
        let mut result = Complex64::new(0.0, 0.0);
        for i in 0..self.dim {
            result += a.coeffs[i] * self.counit[i];
        }
        result
    }

    /// Check coassociativity: (Δ ⊗ id) ∘ Δ = (id ⊗ Δ) ∘ Δ
    pub fn is_coassociative(&self) -> bool {
        let eps = 1e-8;
        for i in 0..self.dim {
            // (Δ ⊗ id) ∘ Δ(e_i) = sum_{j,k} comult[i][j][k] * Δ(e_j) ⊗ e_k
            // = sum_{j,k,l,m} comult[i][j][k] * comult[j][l][m] * e_l ⊗ e_m ⊗ e_k
            let mut lhs = vec![vec![vec![Complex64::new(0.0, 0.0); self.dim]; self.dim]; self.dim];
            for j in 0..self.dim {
                for k in 0..self.dim {
                    let c = self.comult[i][j][k];
                    if c.norm() < 1e-15 {
                        continue;
                    }
                    for l in 0..self.dim {
                        for m in 0..self.dim {
                            lhs[l][m][k] += c * self.comult[j][l][m];
                        }
                    }
                }
            }

            // (id ⊗ Δ) ∘ Δ(e_i) = sum_{j,k} comult[i][j][k] * e_j ⊗ Δ(e_k)
            // = sum_{j,k,l,m} comult[i][j][k] * comult[k][l][m] * e_j ⊗ e_l ⊗ e_m
            let mut rhs = vec![vec![vec![Complex64::new(0.0, 0.0); self.dim]; self.dim]; self.dim];
            for j in 0..self.dim {
                for k in 0..self.dim {
                    let c = self.comult[i][j][k];
                    if c.norm() < 1e-15 {
                        continue;
                    }
                    for l in 0..self.dim {
                        for m in 0..self.dim {
                            rhs[j][l][m] += c * self.comult[k][l][m];
                        }
                    }
                }
            }

            for l in 0..self.dim {
                for m in 0..self.dim {
                    for n in 0..self.dim {
                        if (lhs[l][m][n] - rhs[l][m][n]).norm() > eps {
                            return false;
                        }
                    }
                }
            }
        }
        true
    }

    /// Check counit axiom: (ε ⊗ id) ∘ Δ = id = (id ⊗ ε) ∘ Δ
    pub fn satisfies_counit_axiom(&self) -> bool {
        let eps = 1e-8;
        for i in 0..self.dim {
            // (ε ⊗ id) ∘ Δ(e_i) = sum_{j,k} comult[i][j][k] * ε(e_j) * e_k
            for k in 0..self.dim {
                let lhs: Field = (0..self.dim)
                    .map(|j| self.comult[i][j][k] * self.counit[j])
                    .sum();
                let expected = if k == i {
                    Complex64::new(1.0, 0.0)
                } else {
                    Complex64::new(0.0, 0.0)
                };
                if (lhs - expected).norm() > eps {
                    return false;
                }
            }

            // (id ⊗ ε) ∘ Δ(e_i)
            for j in 0..self.dim {
                let rhs: Field = (0..self.dim)
                    .map(|k| self.comult[i][j][k] * self.counit[k])
                    .sum();
                let expected = if j == i {
                    Complex64::new(1.0, 0.0)
                } else {
                    Complex64::new(0.0, 0.0)
                };
                if (rhs - expected).norm() > eps {
                    return false;
                }
            }
        }
        true
    }
}

/// Group-like coalgebra from a finite group.
pub fn group_coalgebra(n: usize) -> Coalgebra {
    let mut comult = vec![vec![vec![Complex64::new(0.0, 0.0); n]; n]; n];
    // For group-like: ε(g) = 1 for all g
    let counit = vec![Complex64::new(1.0, 0.0); n];

    // Δ(e_i) = sum_j e_j ⊗ e_{j^{-1} i} for a cyclic group Z_n
    for i in 0..n {
        let _inv_i = (n - i) % n; // j^{-1} where j=i in Z_n: (-i mod n)
        // Δ(e_i) = e_i ⊗ e_i (group-like elements)
        comult[i][i][i] = Complex64::new(1.0, 0.0);
    }

    Coalgebra { dim: n, comult, counit }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hopf_algebra::real_val;

    #[test]
    fn test_trivial_coalgebra() {
        let ca = Coalgebra::trivial();
        assert_eq!(ca.dim, 1);
        assert!(ca.is_coassociative());
        assert!(ca.satisfies_counit_axiom());
    }

    #[test]
    fn test_trivial_comultiply() {
        let ca = Coalgebra::trivial();
        let e = AlgebraElement::basis(1, 0);
        let result = ca.comultiply(&e);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].1, 0);
        assert_eq!(result[0].2, 0);
    }

    #[test]
    fn test_trivial_counit() {
        let ca = Coalgebra::trivial();
        let e = AlgebraElement::basis(1, 0);
        let eps = ca.apply_counit(&e);
        assert!((eps - real_val(1.0)).norm() < 1e-10);
    }

    #[test]
    fn test_group_coalgebra_coassoc() {
        let ca = group_coalgebra(3);
        assert!(ca.is_coassociative());
    }

    #[test]
    fn test_group_coalgebra_counit_axiom() {
        let ca = group_coalgebra(3);
        assert!(ca.satisfies_counit_axiom());
    }

    #[test]
    fn test_group_coalgebra_grouplike() {
        let ca = group_coalgebra(3);
        // Each e_i should be group-like: Δ(e_i) = e_i ⊗ e_i
        for i in 0..3 {
            let e = AlgebraElement::basis(3, i);
            let result = ca.comultiply(&e);
            assert_eq!(result.len(), 1);
            assert_eq!(result[0].1, i);
            assert_eq!(result[0].2, i);
        }
    }

    #[test]
    fn test_group_coalgebra_counit_identity() {
        let ca = group_coalgebra(4);
        // For group-like: ε(g) = 1 for all g
        for i in 0..4 {
            assert!((ca.counit[i] - real_val(1.0)).norm() < 1e-10, "counit[{}] = {:?}", i, ca.counit[i]);
        }
    }

    #[test]
    fn test_counit_on_sum() {
        let ca = group_coalgebra(2);
        let sum = AlgebraElement {
            coeffs: vec![real_val(1.0), real_val(1.0)],
        };
        let eps = ca.apply_counit(&sum);
        // ε(Σ e_i) = Σ ε(e_i) = 1 + 1 = 2
        assert!((eps - real_val(2.0)).norm() < 1e-10);
    }
}
