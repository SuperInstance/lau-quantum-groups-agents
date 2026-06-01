use crate::hopf_algebra::{Algebra, AlgebraElement, Field};
use num_complex::Complex64;
use crate::coalgebra::Coalgebra;
use serde::{Deserialize, Serialize};

/// A Hopf algebra: algebra + coalgebra + antipode, satisfying compatibility conditions.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HopfAlgebra {
    pub algebra: Algebra,
    pub coalgebra: Coalgebra,
    /// Antipode: S(e_i) = sum_j antipode[i][j] * e_j
    pub antipode: Vec<Vec<Field>>,
}

impl HopfAlgebra {
    pub fn new(algebra: Algebra, coalgebra: Coalgebra, antipode: Vec<Vec<Field>>) -> Self {
        assert_eq!(algebra.dim, coalgebra.dim);
        assert_eq!(antipode.len(), algebra.dim);
        Self { algebra, coalgebra, antipode }
    }

    /// Trivial 1-dimensional Hopf algebra (the base field k).
    pub fn trivial() -> Self {
        let alg = Algebra::trivial();
        let ca = Coalgebra::trivial();
        let antipode = vec![vec![Complex64::new(1.0, 0.0)]];
        Self { algebra: alg, coalgebra: ca, antipode }
    }

    /// Apply antipode to an element.
    pub fn apply_antipode(&self, a: &AlgebraElement) -> AlgebraElement {
        let mut result = vec![Complex64::new(0.0, 0.0); self.algebra.dim];
        for i in 0..self.algebra.dim {
            for j in 0..self.algebra.dim {
                result[j] = result[j] + a.coeffs[i] * self.antipode[i][j];
            }
        }
        AlgebraElement { coeffs: result }
    }

    /// Check antipode axiom: S * id = u ∘ ε = id * S
    /// where * is convolution product.
    pub fn satisfies_antipode_axiom(&self) -> bool {
        let eps = 1e-8;
        let dim = self.algebra.dim;

        for i in 0..dim {
            // S * id: sum_j S(e_j)_right * id(e_j)_left via convolution
            // Convolution: (f * g)(x) = mu ∘ (f ⊗ g) ∘ Δ(x)
            // For f=S, g=id: sum_j,k comult[i][j][k] * S(e_j) * e_k
            let mut lhs = vec![Complex64::new(0.0, 0.0); dim];
            for j in 0..dim {
                for k in 0..dim {
                    let c = self.coalgebra.comult[i][j][k];
                    if c.norm() < 1e-15 {
                        continue;
                    }
                    // S(e_j) * e_k
                    let s_j = &self.antipode[j];
                    let e_k = AlgebraElement::basis(dim, k);
                    let s_elem = AlgebraElement { coeffs: s_j.to_vec() };
                    let prod = self.algebra.multiply(&s_elem, &e_k);
                    for m in 0..dim {
                        lhs[m] = lhs[m] + c * prod.coeffs[m];
                    }
                }
            }
            // Should equal ε(e_i) * 1 = counit[i] * unit
            let expected = self.coalgebra.counit[i];
            let unit = self.algebra.unit();
            for m in 0..dim {
                let exp_m = expected * unit.coeffs[m];
                if (lhs[m] - exp_m).norm() > eps {
                    return false;
                }
            }

            // id * S: sum_j,k comult[i][j][k] * e_j * S(e_k)
            let mut rhs = vec![Complex64::new(0.0, 0.0); dim];
            for j in 0..dim {
                for k in 0..dim {
                    let c = self.coalgebra.comult[i][j][k];
                    if c.norm() < 1e-15 {
                        continue;
                    }
                    let e_j = AlgebraElement::basis(dim, j);
                    let s_k = AlgebraElement { coeffs: self.antipode[k].to_vec() };
                    let prod = self.algebra.multiply(&e_j, &s_k);
                    for m in 0..dim {
                        rhs[m] = rhs[m] + c * prod.coeffs[m];
                    }
                }
            }
            for m in 0..dim {
                let exp_m = expected * unit.coeffs[m];
                if (rhs[m] - exp_m).norm() > eps {
                    return false;
                }
            }
        }
        true
    }

    /// Check that multiplication is a coalgebra morphism (bialgebra condition).
    pub fn is_bialgebra(&self) -> bool {
        let eps = 1e-10;
        let dim = self.algebra.dim;

        // Check: Δ(xy) = Δ(x) * Δ(y) (using tensor product multiplication)
        for i in 0..dim {
            for j in 0..dim {
                // x = e_i, y = e_j
                // xy = sum_m mult[i][j][m] * e_m
                // Δ(xy) = sum_m mult[i][j][m] * Δ(e_m)
                let mut lhs = vec![vec![Complex64::new(0.0, 0.0); dim]; dim];
                for m in 0..dim {
                    let cij = self.algebra.mult[i][j][m];
                    if cij.norm() < 1e-15 {
                        continue;
                    }
                    for a in 0..dim {
                        for b in 0..dim {
                            lhs[a][b] = lhs[a][b] + cij * self.coalgebra.comult[m][a][b];
                        }
                    }
                }

                // Δ(x) * Δ(y) = (sum_{a,b} comult[i][a][b] e_a ⊗ e_b) *
                //               (sum_{c,d} comult[j][c][d] e_c ⊗ e_d)
                // = sum_{a,b,c,d} comult[i][a][b] * comult[j][c][d] * (e_a*e_c) ⊗ (e_b*e_d)
                let mut rhs = vec![vec![Complex64::new(0.0, 0.0); dim]; dim];
                for a in 0..dim {
                    for b in 0..dim {
                        let ci = self.coalgebra.comult[i][a][b];
                        if ci.norm() < 1e-15 {
                            continue;
                        }
                        for c in 0..dim {
                            for d in 0..dim {
                                let cj = self.coalgebra.comult[j][c][d];
                                if cj.norm() < 1e-15 {
                                    continue;
                                }
                                // e_a * e_c = sum_p mult[a][c][p] e_p
                                // e_b * e_d = sum_q mult[b][d][q] e_q
                                for p in 0..dim {
                                    for q in 0..dim {
                                        rhs[p][q] = rhs[p][q]
                                            + ci * cj * self.algebra.mult[a][c][p] * self.algebra.mult[b][d][q];
                                    }
                                }
                            }
                        }
                    }
                }

                for a in 0..dim {
                    for b in 0..dim {
                        if (lhs[a][b] - rhs[a][b]).norm() > eps {
                            return false;
                        }
                    }
                }
            }
        }
        true
    }

    /// Check that this is a valid Hopf algebra (bialgebra + antipode axiom).
    pub fn is_valid_hopf(&self) -> bool {
        self.algebra.is_associative()
            && self.coalgebra.is_coassociative()
            && self.coalgebra.satisfies_counit_axiom()
            && self.is_bialgebra()
            && self.satisfies_antipode_axiom()
    }

    /// Group Hopf algebra k[G] for a cyclic group of order n.
    pub fn cyclic_group(n: usize) -> Self {
        let _dim = n;
        let alg = cyclic_algebra(n);
        let ca = cyclic_coalgebra(n);
        let antipode = cyclic_antipode(n);
        Self { algebra: alg, coalgebra: ca, antipode }
    }
}

fn cyclic_algebra(n: usize) -> Algebra {
    let mut mult = vec![vec![vec![Complex64::new(0.0, 0.0); n]; n]; n];
    // Multiplication: e_i * e_j = e_{(i+j) mod n} (convolution algebra of Z_n)
    for i in 0..n {
        for j in 0..n {
            let k = (i + j) % n;
            mult[i][j][k] = Complex64::new(1.0, 0.0);
        }
    }
    // Unit is e_0
    Algebra { dim: n, mult }
}

fn cyclic_coalgebra(n: usize) -> Coalgebra {
    let mut comult = vec![vec![vec![Complex64::new(0.0, 0.0); n]; n]; n];
    let mut counit = vec![Complex64::new(0.0, 0.0); n];
    // Group-like: Δ(e_i) = e_i ⊗ e_i
    for i in 0..n {
        comult[i][i][i] = Complex64::new(1.0, 0.0);
    }
    // Counit: ε(e_0) = 1, ε(e_i) = 1 for all (since it's the group algebra)
    // Actually for group algebra k[Z_n], counit is ε(g) = 1 for all g
    for i in 0..n {
        counit[i] = Complex64::new(1.0, 0.0);
    }
    Coalgebra { dim: n, comult, counit }
}

fn cyclic_antipode(n: usize) -> Vec<Vec<Field>> {
    let mut antipode = vec![vec![Complex64::new(0.0, 0.0); n]; n];
    // S(e_i) = e_{n-i mod n} (inverse in Z_n)
    for i in 0..n {
        let inv = (n - i) % n;
        antipode[i][inv] = Complex64::new(1.0, 0.0);
    }
    antipode
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hopf_algebra::real_val;

    #[test]
    fn test_trivial_hopf() {
        let h = HopfAlgebra::trivial();
        assert!(h.is_valid_hopf());
    }

    #[test]
    fn test_trivial_antipode() {
        let h = HopfAlgebra::trivial();
        let e = AlgebraElement::basis(1, 0);
        let s = h.apply_antipode(&e);
        assert!((s.coeffs[0] - real_val(1.0)).norm() < 1e-10);
    }

    #[test]
    fn test_cyclic_z2() {
        let h = HopfAlgebra::cyclic_group(2);
        assert!(h.is_valid_hopf());
    }

    #[test]
    fn test_cyclic_z3() {
        let h = HopfAlgebra::cyclic_group(3);
        assert!(h.is_valid_hopf());
    }

    #[test]
    fn test_cyclic_z4() {
        let h = HopfAlgebra::cyclic_group(4);
        assert!(h.is_valid_hopf());
    }

    #[test]
    fn test_cyclic_z5() {
        let h = HopfAlgebra::cyclic_group(5);
        assert!(h.is_valid_hopf());
    }

    #[test]
    fn test_antipode_z3_involution() {
        let h = HopfAlgebra::cyclic_group(3);
        let e = AlgebraElement::basis(3, 1);
        let s = h.apply_antipode(&e);
        // S(e_1) = e_2 (since 3-1 = 2)
        assert!((s.coeffs[2] - real_val(1.0)).norm() < 1e-10);
        assert!(s.coeffs[0].norm() < 1e-10);
        assert!(s.coeffs[1].norm() < 1e-10);
    }

    #[test]
    fn test_antipode_squared_z3() {
        let h = HopfAlgebra::cyclic_group(3);
        let e = AlgebraElement::basis(3, 1);
        let s1 = h.apply_antipode(&e);
        let s2 = h.apply_antipode(&s1);
        // S(S(e_1)) = S(e_2) = e_1
        assert!((s2.coeffs[1] - real_val(1.0)).norm() < 1e-10);
    }

    #[test]
    fn test_bialgebra_condition_z2() {
        let h = HopfAlgebra::cyclic_group(2);
        assert!(h.is_bialgebra());
    }

    #[test]
    fn test_multiply_in_cyclic() {
        let h = HopfAlgebra::cyclic_group(3);
        let e1 = AlgebraElement::basis(3, 1);
        let e2 = AlgebraElement::basis(3, 2);
        let prod = h.algebra.multiply(&e1, &e2);
        // e_1 * e_2 = e_0 (since (1+2)%3 = 0)
        assert!((prod.coeffs[0] - real_val(1.0)).norm() < 1e-10);
    }
}
