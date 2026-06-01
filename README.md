# lau-quantum-groups-agents

**Quantum groups, Hopf algebras, and q-deformed representation theory in Rust.**

This crate implements the algebraic machinery of quantum groups: associative algebras with structure constants, coalgebras with comultiplication and counit, Hopf algebras with antipodes, universal enveloping algebras, q-deformations of U(sl(2)), universal R-matrices and the Yang-Baxter equation, ribbon categories with quantum trace and link invariants, and representation theory with Clebsch-Gordan decomposition — all backed by **88 tests**.

---

## What This Does

A **quantum group** is not a group — it's a Hopf algebra (or more precisely, a deformation of a universal enveloping algebra). The "quantum" refers to the deformation parameter *q*: when q = 1 you recover the classical (Lie algebra) theory, and when q ≠ 1 you get a rich non-commutative geometry with braided categories, knot invariants, and topological quantum field theory.

This crate gives you:

- **Finite-dimensional algebras** with structure constants and associativity checking
- **Coalgebras** with comultiplication Δ, counit ε, and coassociativity verification
- **Hopf algebras** combining algebra + coalgebra + antipode S with full axiom checking
- **Lie algebras** with bracket, antisymmetry, and Jacobi identity (includes sl(2))
- **Universal enveloping algebras** U(g) with PBW basis
- **q-deformed U_q(sl(2))** with q-numbers, q-factorials, q-binomials, and representation matrices
- **Universal R-matrices** satisfying the Yang-Baxter equation
- **Ribbon categories** with twist, braiding, quantum trace, and Jones polynomial computations
- **Representation theory** with irreps, characters, tensor product decomposition, and Wigner 3j symbols

---

## Key Idea

Deform the commutation relations of a Lie algebra by introducing a parameter q. For U_q(sl(2)):

```
KEK⁻¹ = q²E
KFK⁻¹ = q⁻²F
[E, F] = (K − K⁻¹) / (q − q⁻¹)
```

When q → 1, these reduce to the ordinary sl(2) relations. The R-matrix provides a solution to the Yang-Baxter equation R₁₂ R₁₃ R₂₃ = R₂₃ R₁₃ R₁₂, which is the algebraic backbone of braided categories and knot invariants.

---

## Install

```toml
[dependencies]
lau-quantum-groups-agents = "0.1.0"
```

Requires **Rust 2021 edition**. Dependencies: `nalgebra` (with serde), `num-complex` (with serde), `serde`, `serde_json`.

---

## Quick Start

### Associative algebra

```rust
use lau_quantum_groups_agents::hopf_algebra::*;

// Matrix algebra M_2(ℂ) — automatically constructed
let alg = Algebra::matrix_algebra(2);
assert!(alg.is_associative());
```

### Hopf algebra from a cyclic group

```rust
use lau_quantum_groups_agents::antipode::*;

// Hopf algebra k[ℤ₃] — the group algebra of the cyclic group
let h = HopfAlgebra::cyclic_group(3);
assert!(h.is_valid_hopf());

// Antipode is the group inverse
let e1 = AlgebraElement::basis(3, 1);
let s = h.apply_antipode(&e1);
// S(e₁) = e₂ in ℤ₃ (since 3 - 1 = 2)
```

### Lie algebra sl(2)

```rust
use lau_quantum_groups_agents::universal_enveloping::*;

let sl2 = LieAlgebra::sl2();
assert!(sl2.is_antisymmetric());
assert!(sl2.satisfies_jacobi());

// [e, f] = h
let e = LieAlgebraElement::basis(3, 0);
let f = LieAlgebraElement::basis(3, 1);
let comm = sl2.bracket(&e, &f);
```

### q-deformed quantum group U_q(sl(2))

```rust
use lau_quantum_groups_agents::quantum_deformation::*;

let uq = QuantumUqSL2::new(2.0, 2);

// q-numbers
let q = QParameter::real(2.0);
println!("[2]_q = {:?}", q.q_number(2));       // 2.5
println!("[2]!_q = {:?}", q.q_factorial(2));   // 2.5
println!("[3 choose 1]_q = {:?}", q.q_binomial(3, 1));

// Representation matrices for spin-1
let rep = uq.representation(1);
println!("K = {:?}", rep.k_matrix(1));
println!("dim = {}", rep.dim()); // 3

// Verify quantum relations
assert!(uq.check_kek_relation(1)); // KEK⁻¹ = q²E
assert!(uq.check_kfk_relation(1)); // KFK⁻¹ = q⁻²F
```

### R-matrix and Yang-Baxter equation

```rust
use lau_quantum_groups_agents::r_matrix::*;

let r = RMatrix::fundamental(2.0);
assert!(r.check_yang_baxter()); // R₁₂R₁₃R₂₃ = R₂₃R₁₃R₁₂

let rinv = r.inverse();
let det = r.quantum_determinant();
```

### Ribbon category and link invariants

```rust
use lau_quantum_groups_agents::ribbon::*;

let rc = RibbonCategory::new(2.0);

// Braid relation
assert!(rc.check_braid_relation());

// Quantum trace and Jones polynomial
let unknot = rc.jones_polynomial_unknot(); // q + q⁻¹

// Clebsch-Gordan decomposition
let cat = BraidedMonoidalCategory::new(2.0);
let decomp = cat.clebsch_gordan_decomposition(1, 1); // [0, 1, 2]
```

### Representation theory

```rust
use lau_quantum_groups_agents::representation::*;

// All irreps of U_q(sl(2)) up to spin-1
let irreps = RepTheory::irreps_uqsl2(2.0, 1);
assert_eq!(irreps[0].dim, 1); // spin-0
assert_eq!(irreps[1].dim, 3); // spin-1

// Tensor product decomposition: V₁ ⊗ V₁ = V₀ ⊕ V₁ ⊕ V₂
let decomp = RepTheory::tensor_product_decomposition(1, 1);
assert_eq!(decomp, vec![0, 1, 2]);

// Casimir eigenvalue
let q = QParameter::real(2.0);
let c = RepTheory::casimir_eigenvalue(&q, 1); // [1]_q [2]_q = 2.5
```

---

## API Reference

### `hopf_algebra` — Algebras and Algebra Elements

| Type / Function | Description |
|---|---|
| `AlgebraElement` | Element of a finite-dimensional algebra (vector of complex coefficients) |
| `Algebra` | Associative algebra with structure constants; supports multiply, unit, associativity check |
| `Algebra::matrix_algebra(n)` | M_n(ℂ) — the full matrix algebra |
| `tensor_multiply` | Multiplication in tensor product algebras |

### `coalgebra` — Coalgebras

| Type / Function | Description |
|---|---|
| `Coalgebra` | Coalgebra with comultiplication Δ and counit ε; coassociativity and counit axiom checks |
| `group_coalgebra(n)` | Group-like coalgebra for ℤ_n with Δ(g) = g ⊗ g |

### `antipode` — Hopf Algebras

| Type / Function | Description |
|---|---|
| `HopfAlgebra` | Full Hopf algebra with algebra + coalgebra + antipode; validity check (bialgebra + antipode axiom) |
| `HopfAlgebra::cyclic_group(n)` | Group Hopf algebra k[ℤ_n] with multiplication, group-like coproduct, and inverse antipode |

### `universal_enveloping` — Lie Algebras & Enveloping Algebras

| Type / Function | Description |
|---|---|
| `LieAlgebra` | Lie algebra with structure constants; antisymmetry and Jacobi checks |
| `LieAlgebra::sl2()` | sl(2,ℂ) with basis {e, f, h}: [h,e]=2e, [h,f]=-2f, [e,f]=h |
| `UniversalEnvelopingAlgebra` | U(g) with PBW monomials up to a truncation degree |

### `quantum_deformation` — q-Deformed Quantum Groups

| Type / Function | Description |
|---|---|
| `QParameter` | q with q^n, q-number [n]_q, q-factorial, q-binomial, root-of-unity check |
| `QuantumUqSL2` | U_q(sl(2)) with K, E, F matrices for any spin-j representation |
| `QuantumRepresentation` | Spin-j rep matrices with Casimir element |

### `r_matrix` — Universal R-Matrix

| Type / Function | Description |
|---|---|
| `RMatrix` | R-matrix for U_q(sl(2)) in fundamental rep; Yang-Baxter check, inverse, quantum determinant |
| `kron` | Kronecker product of complex matrices |
| `flip` | Flip (swap) operator σ: V⊗W → W⊗V |

### `ribbon` — Ribbon Categories

| Type / Function | Description |
|---|---|
| `RibbonCategory` | Ribbon category with twist θ, braiding σ, quantum trace, Jones polynomial |
| `BraidedMonoidalCategory` | Monoidal category with Clebsch-Gordan decomposition |

### `representation` — Representation Theory

| Type / Function | Description |
|---|---|
| `Irrep` | Irreducible representation with matrices, character, dimension |
| `RepTheory` | Irrep enumeration, tensor product decomposition, Casimir eigenvalue, Wigner 3j, direct sum |

---

## How It Works

1. **Start with a Lie algebra** g (e.g., sl(2)) and its universal enveloping algebra U(g).
2. **Deform** the commutation relations by introducing a parameter q, getting U_q(g).
3. **Build the R-matrix** — an element of U_q(g) ⊗ U_q(g) that satisfies the Yang-Baxter equation.
4. **Construct representations** using q-numbers instead of ordinary integers.
5. **Form a ribbon category** from the representation theory, enabling quantum traces and knot invariants.

---

## The Math

### Hopf Algebras

A **Hopf algebra** H is simultaneously an algebra (A, μ, η) and a coalgebra (A, Δ, ε) with an antipode S, satisfying:

- **Multiplicativity of Δ**: Δ(xy) = Δ(x) · Δ(y) (using tensor product multiplication)
- **Antipode axiom**: S ∗ id = id ∗ S = η ∘ ε (convolution product)

The group algebra k[G] of a finite group G is a Hopf algebra with Δ(g) = g ⊗ g, ε(g) = 1, S(g) = g⁻¹.

### q-Numbers

The **q-number** [n]_q is defined as:

```
[n]_q = (q^n − q^{−n}) / (q − q^{−1})
```

When q → 1, [n]_q → n. The q-factorial and q-binomial are built from q-numbers exactly as their classical counterparts.

### U_q(sl(2))

The quantum group U_q(sl(2)) has generators E, F, K, K⁻¹ with relations:

```
KK⁻¹ = K⁻¹K = 1
KEK⁻¹ = q²E
KFK⁻¹ = q⁻²F
[E, F] = (K − K⁻¹)/(q − q⁻¹)
```

The spin-j representation has dimension 2j+1 with:

```
K|j,m⟩ = q^{2m}|j,m⟩
E|j,m⟩ = √([j−m]_q [j+m+1]_q) |j,m+1⟩
F|j,m⟩ = √([j+m]_q [j−m+1]_q) |j,m−1⟩
```

### R-Matrix and Yang-Baxter Equation

The universal R-matrix for U_q(sl(2)) in the fundamental representation:

```
R = q^{1/2}|00⟩⟨00| + q^{−1/2}|11⟩⟨11|
  + q^{−1/2}|01⟩⟨01| + q^{−1/2}|10⟩⟨10|
  + (q − q⁻¹)q^{−1/2}|10⟩⟨01|
```

This satisfies the **Yang-Baxter equation**:

```
R₁₂ R₁₃ R₂₃ = R₂₃ R₁₃ R₁₂
```

which is the algebraic foundation of braiding in quantum groups.

### Ribbon Categories

A **ribbon category** is a braided monoidal category with a twist θ satisfying:

```
θ_{V⊗W} = σ²_{VW} ∘ (θ_V ⊗ θ_W)
```

The **quantum trace** is:

```
tr_q(x) = tr(x ∘ K⁻¹ ∘ v)
```

where v is the ribbon element. This gives rise to **knot invariants** like the Jones polynomial.

---

## Test Coverage

**88 tests** across all modules:

| Module | Tests | What's covered |
|---|---|---|
| `hopf_algebra` | 10 | Element operations, trivial algebra, M_2, M_3 associativity |
| `coalgebra` | 8 | Trivial coalgebra, comultiplication, counit, group coalgebra, coassociativity |
| `antipode` | 12 | Trivial/cyclic Hopf algebras (ℤ₂–ℤ₅), antipode axiom, bialgebra condition, multiplication |
| `universal_enveloping` | 11 | sl(2) brackets, antisymmetry, Jacobi, UEA construction, PBW monomials |
| `quantum_deformation` | 15 | q-numbers, q-factorials, q-binomials, classical limit, root-of-unity, U_q(sl(2)) relations |
| `r_matrix` | 12 | R-matrix entries, inverse, Yang-Baxter equation, Kronecker product, flip, determinant |
| `ribbon` | 12 | Twist, braiding, braid relation, quantum trace, Jones polynomial, Clebsch-Gordan |
| `representation` | 14 | Irreps, characters, tensor decomposition, Casimir, Wigner 3j, direct sum |

Run them with:

```bash
cargo test
```

---

## License

MIT
