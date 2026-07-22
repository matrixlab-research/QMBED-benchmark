# quspin — Reference unit tests (visible examples)

> Generated from frozen SpecBundle `sha256:a9093788a8b2687689113102daaef231e58d33f6188b74927e9b60c097ee1721` — do not hand-edit.

Worked examples for orientation and local self-testing — feel free to turn
them into your own test suite. They are NOT the gate: the private verification
suite additionally holds edge-stratified held-out values, randomized property
checks, and real integration usage (see CONTRACT.md).

| symbol | input | expected | comparison |
|---|---|---|---|
| `error::QuSpinError` | operation=Basis.state;index_out_of_range=true | category=StateNotInBasis;panic=false | `custom:error-category` |
| `basis::Basis` | SpinBasis1D(sites=4, up=2); state=0b0101 | dimension=6;round_trip=true | `custom:basis-observation` |
| `basis::SpinBasis1D` | sites=4;spin_twice=1;up=2;no_symmetry;pauli=false | dimension=6;states_have_weight=2 | `custom:basis-observation` |
| `basis::BosonBasis1D` | sites=3;particles=2;states_per_site=3 | dimension=6;occupations_sum_to=2 | `custom:basis-observation` |
| `basis::SpinlessFermionBasis1D` | sites=4;particles=2;no_symmetry | dimension=6;states_have_weight=2 | `custom:basis-observation` |
| `basis::SpinfulFermionBasis1D` | sites=3;particles_up=1;particles_down=1 | dimension=9;sector_counts=[1,1] | `custom:basis-observation` |
| `basis::UserBasis` | sites=4;filter=periodic_no_adjacent_ones;operators=[x,z] | dimension=7;round_trip=true;universal_builder=true | `custom:basis-observation` |
| `operator::OperatorTerm` | operator=zz;couplings=[(1.0,[0,1]),(-0.5,[1,2])] | arity=2;coupling_count=2;parsed_once=true | `custom:parsed-term` |
| `operator::Coupling` | coefficient=1-0.5i;sites=[2,0] | coefficient=1-0.5i;sites=[2,0] | `isequal` |
| `operator::MatrixFormat` | requested=[Dense,Csc,Csr,Dia,MatrixFree] | variants=[Dense,Csc,Csr,Dia,MatrixFree] | `isequal` |
| `operator::LinearOperator` | shape=[2,3];matrix=[[1,0,2],[0,-1,1]];input=[1,2,3] | output=[7,1];shape=[2,3] | `num:1e-12:1e-10` |
| `operator::OperatorBuilder` | spin_half_dimer;terms=[zz:1,+-:0.5,-+:0.5];pauli=false;format=Csc | shape=[4,4];eigenvalues=[-0.75,0.25,0.25,0.25];hermitian=true | `num:1e-12:1e-10` |
| `solve::EigshOptions` | eigenpairs=2;target=SmallestAlgebraic;krylov_dimension=8;tolerance=1e-10;max_iterations=500;seed=7 | valid=true | `isequal` |
| `solve::eigsh` | spin_half_dimer_heisenberg;k=2;target=SmallestAlgebraic;tolerance=1e-12 | eigenvalues=[-0.75,0.25];residuals_below=1e-12 | `num:1e-09:1e-08` |
| `solve::EvolutionOptions` | times=[0,0.5,1];krylov_dimension=12;tolerance=1e-9;max_substeps=100;hamiltonian=true | valid=true;state_count=3 | `isequal` |
| `solve::evolve` | H=diag[0,2];initial=[1,1]/sqrt(2);times=[0,pi/2];hamiltonian=true | states=[[1,1]/sqrt(2),[1,-1]/sqrt(2)];norms=[1,1] | `num:1e-09:1e-08` |
| `dynamics::DriveStep` | hamiltonian=diag[0,2];duration=0.5 | valid=true;shape=[2,2];duration=0.5 | `custom:drive-step` |
| `dynamics::Floquet` | steps=[(H=diag[0,2],duration=pi/2),(H=zero2,duration=3)];input=[1,1]/sqrt(2) | output=[1,-1]/sqrt(2);norm=1 | `num:1e-09:1e-08` |
| `dynamics::SpectrumOptions` | frequencies=[-1,0,1];reference_energy=0.2;broadening=0.1;krylov_dimension=16;tolerance=1e-8 | valid=true;output_length=3 | `isequal` |
| `dynamics::spectral_function` | H_target=diag[1];source=[1];probe=[[1]];E0=0;omega=[1];eta=0.5 | spectrum=[0.6366197723675814] | `num:1e-08:1e-07` |
| `measure::Subspace` | ambient=3;rank=2;columns=[e1,e2] | rank=2;projector=diag[1,1,0] | `num:1e-12:1e-10` |
| `measure::subspace_fidelity` | left=span(e1,e2);right=span((e1+e2)/sqrt(2),(e1-e2)/sqrt(2)) | 1.0 | `num:1e-12:1e-10` |
| `workflow::LindbladGenerator` | H=zero2;jump=sqrt(1)*lowering;rho0=excited_state;t=ln(2) | excited_population=0.5;trace=1;hermitian=true | `num:1e-09:1e-08` |
