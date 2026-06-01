# lau-contact-agents

**Contact geometry for autonomous agents** — contact manifolds, Reeb vector fields, contactomorphisms, Legendrian submanifolds, contact Hamiltonian dynamics, Darboux charts, Gray stability, the Weinstein conjecture, contact homology, and thermodynamic agent dynamics on contact manifolds.

## What This Does

This crate implements the core structures of **contact geometry** — the odd-dimensional sibling of symplectic geometry — and applies them to model autonomous agents whose dynamics are intrinsically dissipative. Every agent lives on a (2n+1)-dimensional contact manifold where position variables track extensive quantities, momentum variables track intensive quantities, and an extra coordinate z accumulates entropy-like displacement.

You get:

- **Contact forms** α and their **Reeb vector fields** on R^{2n+1}
- **Legendrian submanifolds** (maximal submanifolds tangent to the contact structure)
- **Contactomorphisms** (structure-preserving diffeomorphisms)
- **Contact Hamiltonian dynamics** with built-in dissipation
- **Darboux charts** proving local normal forms
- **Gray stability** for isotopies of contact structures
- **Weinstein conjecture** verification (closed Reeb orbits)
- **Contact homology** (a Floer-theoretic invariant built from Reeb orbits)
- **Thermodynamic agents** that dissipate, produce entropy, and interact through coupling

## Key Idea

A **contact form** α on a (2n+1)-dimensional manifold satisfies α ∧ (dα)ⁿ ≠ 0 everywhere. This single condition guarantees that the "contact distribution" ξ = ker(α) is maximally non-integrable — you *cannot* foliate the manifold with n+1-dimensional slices tangent to ξ. This rigidity is what makes contact geometry so powerful: it constrains dynamics in ways that ordinary differential geometry cannot.

For agents, this is *natural*. Thermodynamic systems are irreversible. Entropy always increases. The contact framework captures this through the Reeb vector field and contact Hamiltonian equations, where energy is generically *not* conserved — it dissipates at a rate proportional to the Hamiltonian itself.

## Install

```toml
[dependencies]
lau-contact-agents = "0.1.0"
```

Requires `nalgebra` (with `serde-serialize`) and `serde` + `serde_json`.

## Quick Start

```rust
use lau_contact_agents::contact_form::{ContactForm, ReebVectorField};
use lau_contact_agents::legendrian::LegendrianSubmanifold;
use lau_contact_agents::thermo_dynamics::ThermoAgent;
use nalgebra::DVector;

// Create the standard contact form on R^3: α = dz - y dx
let form = ContactForm::standard(1);
let p = DVector::from_vec(vec![1.0, 2.0, 3.0]);
assert!(form.is_contact_at(&p)); // α ∧ (dα) ≠ 0

// The Reeb vector field is ∂/∂z
let reeb = ReebVectorField::new(1);
assert!(reeb.verify_alpha_condition(&p, &form));  // α(R) = 1
assert!(reeb.verify_dalpha_condition(&p, &form)); // dα(R, ·) = 0

// Standard Legendrian submanifold: y = 0, z = 0
let legendrian = LegendrianSubmanifold::standard(1);
let alpha_at_origin = DVector::from_vec(vec![0.0, 0.0, 1.0]);
assert!(legendrian.is_legendrian(&alpha_at_origin)); // α|_L = 0

// Thermodynamic agent that dissipates over time
let state = DVector::from_vec(vec![0.0, 1.0, 0.0]); // (x=0, y=1, z=0)
let mut agent = ThermoAgent::new(1, state);
let gamma = DVector::from_vec(vec![1.0]);
agent.evolve(0.1, &gamma);
// x moves: ẋ = γy = 1, so x → 0.1
```

## API Reference

### `contact_form` — Contact Forms and Reeb Fields

| Type | Description |
|------|-------------|
| `ContactForm` | A contact 1-form α on a (2n+1)-manifold. Stores dimension parameter `n` and local coefficients. |
| `ReebVectorField` | The unique vector field R satisfying α(R) = 1 and dα(R, ·) = 0. |

**`ContactForm` methods:**

- `standard(n)` — Create the canonical contact form α = dz − Σ yᵢ dxᵢ on R^{2n+1}.
- `from_coefficients(n, coefficients)` — Build from a covector of length 2n+1.
- `evaluate_at(&p)` — Evaluate α at point p, returning the covector α_p.
- `exterior_derivative(&p)` — Compute dα as a skew-symmetric (2n+1)×(2n+1) matrix.
- `volume_form(&p) → f64` — Compute α ∧ (dα)ⁿ (nonzero iff contact).
- `is_contact_at(&p) → bool` — Check the contact condition at a point.
- `dimension() → usize` — Returns 2n+1.

**`ReebVectorField` methods:**

- `new(n)` — Create for an n-dimensional contact structure.
- `evaluate(&p) → DVector<f64>` — Compute R at point p (for standard form: R = ∂/∂z).
- `verify_alpha_condition(&p, &form) → bool` — Check α(R) = 1.
- `verify_dalpha_condition(&p, &form) → bool` — Check dα(R, ·) = 0.

Both types derive `Serialize`/`Deserialize`.

---

### `contactomorphism` — Structure-Preserving Diffeomorphisms

| Type | Description |
|------|-------------|
| `Contactomorphism` | A diffeomorphism φ with φ*α = f·α (strict if f = 1). Stores the Jacobian and strictness flag. |

**Methods:**

- `identity(n)` — The identity map.
- `from_jacobian(n, jacobian, strict)` — Build from a Jacobian matrix.
- `compose(&other)` — Compose two contactomorphisms (matrix multiplication of Jacobians).
- `inverse() → Option<Self>` — Invert (fails if Jacobian is singular).
- `preserves_contact_structure(&alpha_coeffs) → bool` — Check if the pullback of α is proportional to α.
- `z_translation(n, dz)` — Translation in the z-direction (strictly contact).
- `xy_rotation(n, i, theta)` — Rotation in the (xᵢ, yᵢ) plane (strictly contact).
- `dimension() → usize` — Returns 2n+1.

---

### `legendrian` — Legendrian Submanifolds and Knots

| Type | Description |
|------|-------------|
| `LegendrianSubmanifold` | An n-dimensional submanifold of (M^{2n+1}, α) with α\|_L = 0. |
| `LegendrianKnot` | A Legendrian knot in a 3D contact manifold. |

**`LegendrianSubmanifold` methods:**

- `new(n, tangent_basis)` — Create from tangent vectors stored as columns of a (2n+1)×n matrix.
- `standard(n)` — The "flat" Legendrian: y = 0, z = 0, tangent space = span{∂/∂x₁, …, ∂/∂xₙ}.
- `from_generating_function(n, &grad_q, q_val)` — Build from a generating function q: Rⁿ → R.
- `is_legendrian(&alpha_coeffs) → bool` — Verify α vanishes on all tangent vectors.
- `submanifold_dim() → usize` — Returns n.
- `ambient_dim() → usize` — Returns 2n+1.
- `projection_matrix() → DMatrix<f64>` — Orthogonal projection onto the tangent space.
- `is_tangent(&v) → bool` — Check if a vector lies in TL.

**`LegendrianKnot` methods:**

- `from_points(points)` — Build from sampled 3D points.
- `unknot(num_points)` — The standard Legendrian unknot.
- `verify_legendrian() → bool` — Approximate Legendrian check.
- `num_points() → usize`
- `thurston_bennequin() → i32` — Thurston-Bennequin invariant (simplified).

---

### `contact_hamiltonian` — Contact Hamiltonian Dynamics

| Type | Description |
|------|-------------|
| `ContactHamiltonian` | A Hamiltonian function H: M → R with its gradient, producing a contact Hamiltonian vector field X_H. |
| `DissipativeContactSystem` | A contact system with quadratic dissipation H = ½ Σ γᵢ yᵢ². |

**`ContactHamiltonian` methods:**

- `new(n, h, grad_h)` — Create at a point with value h and gradient.
- `vector_field() → DVector<f64>` — Compute X_H (ẋᵢ = ∂H/∂yᵢ, ẏᵢ = −∂H/∂xᵢ).
- `vector_field_full(reeb_h) → DVector<f64>` — Full equations including Reeb coupling.
- `is_conserved() → bool` — Always false for generic contact Hamiltonians (they dissipate).
- `dissipation_rate(reeb_h) → f64` — Returns −R(H)·H.
- `dimension() → usize`

**`DissipativeContactSystem` methods:**

- `new(n, gamma)` — Create with positive dissipation coefficients γ.
- `hamiltonian_value(&y) → f64` — Evaluate H = ½ Σ γᵢ yᵢ².
- `equations_of_motion(&state) → DVector<f64>` — Full (ẋ, ẏ, ż) system.
- `step(&state, dt) → DVector<f64>` — Forward Euler step.
- `total_gamma() → f64` — Sum of dissipation coefficients.
- `is_equilibrium(&state) → bool` — Check if ẋ = ẏ = ż = 0.

---

### `darboux` — Darboux Theorem (Local Normal Form)

| Type | Description |
|------|-------------|
| `DarbouxChart` | A local coordinate system where the contact form becomes α₀ = dz − Σ yᵢ dxᵢ. |

**Methods:**

- `new(n, center)` — Create an identity chart at a center point.
- `with_transform(n, center, transform)` — Create with a specific transformation matrix.
- `to_darboux(&point) → DVector<f64>` — Map from original to Darboux coordinates.
- `from_darboux(&darboux_point) → DVector<f64>` — Inverse map.
- `verify_standard_form(&alpha_coeffs) → bool` — Check that α is standard in chart coordinates.
- `dimension() → usize`

**Free function:**

- `construct_darboux_chart(n, point, &alpha_coeffs, &dalpha_matrix) → DarbouxChart` — Algorithmically construct a Darboux chart from a contact form at a point by finding a symplectic basis for dα on ker(α).

---

### `gray_stability` — Gray Stability Theorem

| Type | Description |
|------|-------------|
| `GrayStability` | Result of checking whether two contact forms are isotopic through contact forms. |

**Methods:**

- `check_path(n, &alpha0, &alpha1, &dalpha0, &dalpha1, steps) → Self` — Interpolate αₜ = (1−t)α₀ + tα₁ and verify the contact condition at each step.
- `moser_vector_field(&alpha_t, &alpha_dot, &reeb) → DVector<f64>` — Compute the Moser-type vector field for the Gray isotopy.
- `verify(n) → bool` — Quick verification for the standard contact form.

Fields: `n`, `is_isotopic: bool`, `interpolation_steps: usize`.

---

### `weinstein` — Weinstein Conjecture

| Type | Description |
|------|-------------|
| `ReebOrbit` | A closed orbit of the Reeb vector field, with period, sample points, and action. |
| `WeinsteinVerifier` | Collects Reeb orbits and checks if the Weinstein conjecture is satisfied. |

**`ReebOrbit` methods:**

- `new(n, period, points)` — Create from sampled points.
- `is_closed() → bool` — Check start ≈ end.
- `action() → f64` — The period (= ∫ α along the orbit).
- `conley_zehnder_index() → i32` — CZ index (simplified).
- `num_points() → usize`

**`WeinsteinVerifier` methods:**

- `new(n)` — Create empty verifier.
- `standard_reeb_orbits(n, num_points) → Vec<ReebOrbit>` — Generate standard Hopf-type orbits.
- `weinstein_satisfied() → bool` — At least one closed Reeb orbit found.
- `verify_s3() → Self` — Verify on S³ with the standard contact form.
- `verify_s2n1(n) → Self` — Verify on S^{2n+1}.

---

### `contact_homology` — Contact Homology (Floer-Type Invariant)

| Type | Description |
|------|-------------|
| `HomologyGenerator` | A closed Reeb orbit as a chain generator, with grading (CZ index) and action. |
| `Differential` | The contact homology differential ∂, counting J-holomorphic curves. Stored as a directed graph. |
| `ContactHomologyComplex` | The full chain complex (generators + differential). |
| `ContactHomology` | The computed homology H_* (ker ∂ / im ∂). |

**`HomologyGenerator`:**

- `new(id, grading, action)` — Create a generator.
- `degree() → i32` — The grading.
- `parity() → usize` — Grading mod 2 (for Z/2 coefficients).

**`Differential`:**

- `zero()` — The zero differential.
- `add_edge(source, target)` — Add a differential coefficient.
- `apply(&[generators]) → Vec<usize>` — Apply ∂ with Z/2 coefficients (toggle presence).
- `check_d_squared_zero() → bool` — Verify ∂² = 0.

**`ContactHomologyComplex`:**

- `new(generators, differential)` — Build the complex.
- `compute_homology() → ContactHomology` — Compute H_* = ker ∂ / im ∂.
- `betti_numbers() → HashMap<Grading, usize>` — Betti numbers by grading.
- `euler_characteristic() → i32` — Σ (−1)^k rank H_k.

**`ContactHomology`:**

- `zero()` — Trivial homology.
- `rank() → usize`
- `is_trivial() → bool`

---

### `thermo_dynamics` — Thermodynamic Agent Dynamics

| Type | Description |
|------|-------------|
| `ThermoAgent` | A single agent on a contact manifold with extensive (x), intensive (y), and entropy (z) variables. |
| `ThermoMultiAgent` | A system of coupled thermodynamic agents. |

**`ThermoAgent` methods:**

- `new(n, state)` — Create with default energy coefficients (all 1).
- `with_energy(n, state, energy_coeffs)` — Custom energy function E(x) = Σ cᵢ xᵢ².
- `extensive() → DVector<f64>` — Get x variables.
- `intensive() → DVector<f64>` — Get y variables.
- `entropy() → f64` — Get z variable.
- `internal_energy() → f64` — Compute E(x) = Σ cᵢ xᵢ².
- `free_energy(temperature) → f64` — Compute F = E − T·S.
- `entropy_production(&gamma) → f64` — Rate of entropy production: Σ γᵢ yᵢ².
- `evolve(dt, &gamma)` — Step forward using contact Hamiltonian dynamics with dissipation γ.
- `is_equilibrium() → bool` — True when all y = 0.
- `thermodynamic_length(&other) → f64` — Euclidean distance between agent states.
- `dimension() → usize` — Returns 2n+1.

**`ThermoMultiAgent` methods:**

- `new(agents, coupling_matrix)` — Create with a coupling matrix between agents.
- `decoupled(agents)` — No inter-agent coupling.
- `num_agents() → usize`
- `total_entropy_production(&gamma) → f64` — Sum across all agents.
- `evolve(dt, &gamma)` — Evolve all agents (coupling + individual dynamics).
- `all_equilibrium() → bool` — Check if every agent is at equilibrium.
- `total_energy() → f64` — Sum of internal energies.

## How It Works

The crate models contact geometry in local coordinates on R^{2n+1} using `nalgebra` dense vectors and matrices. Every structure is concrete and computable:

1. **Contact forms** are represented by their coefficients. The standard form α = dz − Σ yᵢ dxᵢ has point-dependent coefficients (the dxᵢ coefficient at point p is −yᵢ).

2. **Exterior derivatives** are stored as skew-symmetric matrices. For the standard form, dα = Σ dxᵢ ∧ dyᵢ.

3. **The contact condition** α ∧ (dα)ⁿ ≠ 0 is checked by computing the determinant of a (2n+1) × (2n+1) matrix whose first row is α and remaining rows encode dα.

4. **Reeb vector fields** satisfy a pair of linear conditions: α(R) = 1 and dα(R, ·) = 0. For the standard form, this gives R = ∂/∂z.

5. **Contact Hamiltonian dynamics** use the equations ẋᵢ = ∂H/∂yᵢ, ẏᵢ = −∂H/∂xᵢ + yᵢ ∂H/∂z, ż = −H + Σ yᵢ ∂H/∂yᵢ. Unlike symplectic dynamics, H is generically not conserved.

6. **Darboux charts** are constructed by Gram-Schmidt orthogonalization against α to find a basis for the contact distribution ξ = ker(α), then finding a symplectic basis for dα|_ξ.

7. **Gray stability** linearly interpolates between two contact forms and checks the contact condition at each step.

8. **Thermodynamic agents** evolve under contact Hamiltonian flow with quadratic dissipation. The dissipation coefficients γᵢ control how fast each intensive variable decays, and the z variable accumulates the thermodynamic displacement (analogous to entropy production).

## The Math

### Contact Manifolds

A **contact manifold** (M, α) is a (2n+1)-dimensional manifold equipped with a 1-form α such that α ∧ (dα)ⁿ is a volume form (nowhere vanishing). The **contact distribution** ξ = ker(α) is a maximally non-integrable hyperplane distribution.

### Reeb Vector Field

The **Reeb vector field** R is uniquely determined by α(R) = 1 and dα(R, ·) = 0. Its flow preserves the contact structure and is the "trivial" dynamics on a contact manifold.

### Legendrian Submanifolds

A **Legendrian submanifold** L ⊂ (M, α) satisfies dim L = n and α|_L = 0. These are the maximal submanifolds tangent to ξ. In R³ with α = dz − y dx, the x-axis (y = 0, z = 0) is Legendrian.

### Contact Hamiltonian Equations

For a Hamiltonian H: M → R, the contact Hamiltonian vector field X_H satisfies:
- α(X_H) = −H
- dα(X_H, ·) = dH − R(H)α

This gives: dH/dt = −R(H)·H, so energy dissipates generically.

### Darboux Theorem

Every contact form is locally equivalent to α₀ = dz − Σ yᵢ dxᵢ. There are no local invariants — contact geometry is rigid.

### Gray Stability

If α₀ and α₁ are contact forms on a closed manifold with ker(α₀) = ker(α₁), they are isotopic through contact forms. Contact structures on closed manifolds are "stable."

### Weinstein Conjecture

On any closed contact manifold, the Reeb vector field has at least one closed orbit. Proved by Taubes (2007) for dimension 3, known in many higher-dimensional cases.

### Contact Homology

A Floer-type invariant built from closed Reeb orbits. Generators are Reeb orbits, graded by the Conley-Zehnder index. The differential counts rigid J-holomorphic curves in the symplectization. The homology H_* (ker ∂ / im ∂) is an invariant of the contact structure.

### Thermodynamic Contact Geometry

Thermodynamics naturally lives on contact manifolds. The first law of thermodynamics dU = TdS − pdV is a contact form, and the Gibbs contact manifold (with coordinates U, T, S, p, V, …) carries a natural contact structure. Dissipative dynamics on this manifold automatically satisfy the second law (entropy production ≥ 0).

## License

MIT
